/* -*- c++ -*- */
/* 
 * Copyright 2016 Bastille Networks.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include "demod_impl.h"

#define DEBUG_OFF     0
#define DEBUG_INFO    1
#define DEBUG_VERBOSE 2
#define DEBUG         DEBUG_OFF

#define DUMP_IQ       0

#define OVERLAP_DEFAULT 1
#define OVERLAP_FACTOR  16

namespace gr {
  namespace lora {

    demod::sptr
    demod::make( uint8_t   spreading_factor,
                 bool      header,
                 uint8_t   payload_len,
                 uint8_t   cr,
                 bool      crc,
                 bool      low_data_rate,
                 float     beta,
                 uint16_t  fft_factor,
                 uint8_t   peak_search_algorithm,
                 uint16_t  peak_search_phase_k,
                 float     fs_bw_ratio)
    {
      return gnuradio::get_initial_sptr
        (new demod_impl(spreading_factor, header, payload_len, cr, crc, low_data_rate, beta, fft_factor, peak_search_algorithm, peak_search_phase_k, fs_bw_ratio));
    }

    /*
     * The private constructor
     */
    demod_impl::demod_impl( uint8_t   spreading_factor,
                            bool      header,
                            uint8_t   payload_len,
                            uint8_t   cr,
                            bool      crc,
                            bool      low_data_rate,
                            float     beta,
                            uint16_t  fft_factor,
                            uint8_t   peak_search_algorithm,
                            uint16_t  peak_search_phase_k,
                            float     fs_bw_ratio)
      : gr::block("demod",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(0, 0, 0)),
        f_raw("raw.out", std::ios::out),
        f_fft("fft.out", std::ios::out),
        f_up_windowless("up_windowless.out", std::ios::out),
        f_up("up.out", std::ios::out),
        f_down("down.out", std::ios::out),
        d_sf(spreading_factor),
        d_header(header),
        d_payload_len(payload_len),
        d_cr(cr),
        d_crc(crc),
        d_ldr(low_data_rate),
        d_beta(beta),
        d_fft_size_factor(fft_factor),
        d_peak_search_algorithm(peak_search_algorithm),
        d_peak_search_phase_k(peak_search_phase_k)
    {
      assert((d_sf > 5) && (d_sf < 13));
      if (d_sf == 6) assert(!header);
      assert(d_fft_size_factor > 0);
      assert(((int)fs_bw_ratio) == fs_bw_ratio);
      d_p = (int) fs_bw_ratio;

      if (!header) // implicit header mode
      {
        // calculate the total number of symbols in a packet
        d_packet_symbol_len = 8 + std::max((4+d_cr)*(int)std::ceil((2.0*d_payload_len-d_sf+7+4*d_crc-5*!d_header)/(d_sf-2*d_ldr)), 0);
      }

      d_header_port = pmt::mp("header");
      message_port_register_in(d_header_port);
      d_out_port = pmt::mp("out");
      message_port_register_out(d_out_port);

      set_msg_handler(d_header_port, boost::bind(&demod_impl::parse_header, this, _1));

      d_state = S_RESET;

      d_num_symbols = (1 << d_sf);
      d_num_samples = d_p*d_num_symbols;
      d_bin_size = d_fft_size_factor*d_num_symbols;
      d_fft_size = d_fft_size_factor*d_num_samples;
      d_fft = new fft::fft_complex(d_fft_size, true, 1);
      d_overlaps = OVERLAP_DEFAULT;
      d_offset = 0;
      d_preamble_drift_max = d_fft_size_factor * (d_ldr ? 2 : 1);

      d_window = fft::window::build(fft::window::WIN_KAISER, d_num_samples, d_beta);

      // Create local chirp tables.  Each table is 2 chirps long to allow memcpying from arbitrary offsets.
      for (int i = 0; i < d_num_samples; i++) {
        double phase = M_PI/d_p*(i-i*i/(float)d_num_samples);
        d_downchirp.push_back(gr_complex(std::polar(1.0, phase)));
        d_upchirp.push_back(gr_complex(std::polar(1.0, -phase)));
      }

      set_history(DEMOD_HISTORY_DEPTH*d_num_samples);  // Sync is 2.25 chirp periods long
    }

    /*
     * Our virtual destructor.
     */
    demod_impl::~demod_impl()
    {
      delete d_fft;
    }

    uint32_t
    demod_impl::argmax_32f(float *fft_result, float *max_val_p)
    {
      float mag   = abs(fft_result[0]);
      float max_val = mag;
      uint32_t   max_idx = 0;

      for (uint32_t i = 0; i < d_bin_size; i++)
      {
        mag = abs(fft_result[i]);
        if (mag > max_val)
        {
          max_idx = i;
          max_val = mag;
        }
      }

      *max_val_p = max_val;
      return max_idx;
    }

    uint32_t
    demod_impl::search_fft_peak(const lv_32fc_t *fft_result,
                                float *buffer1, float *buffer2,
                                gr_complex *buffer_c, float *max_val_p)
    {
      // size of buffer1:   d_fft_size (float)
      // size of buffer2:   d_bin_size  (float)
      // size of buffer_c:  d_bin_size  (complex)
      uint32_t max_idx = 0;
      *max_val_p = 0;
      if (d_peak_search_algorithm == FFT_PEAK_SEARCH_ABS)
      {
        // fft result magnitude summation
        volk_32fc_magnitude_32f(buffer1, fft_result, d_fft_size);
        volk_32f_x2_add_32f(buffer2, buffer1, &buffer1[d_fft_size-d_bin_size], d_bin_size);

        // Take argmax of returned FFT (similar to MFSK demod)
        max_idx = argmax_32f(buffer2, max_val_p);
      }
      else if (d_peak_search_algorithm == FFT_PEAK_SEARCH_PHASE)
      {
        uint32_t tmp_max_idx;
        float tmp_max_val;
        for (int i = 0; i < d_peak_search_phase_k; i++)
        {
          float phase_offset = 2*M_PI/d_peak_search_phase_k*i;
          tmp_max_idx = fft_add(fft_result, buffer2, buffer_c, &tmp_max_val, phase_offset);
          if (tmp_max_val > *max_val_p)
          {
            *max_val_p = tmp_max_val;
            max_idx = tmp_max_idx;
          }
        }
      }
      else
      {
        max_idx = fft_add(fft_result, buffer2, buffer_c, max_val_p, 0);
      }
      
      return max_idx;
    }

    uint32_t
    demod_impl::fft_add(const lv_32fc_t *fft_result, float *buffer, gr_complex *buffer_c,
                        float *max_val_p, float phase_offset)
    {
      lv_32fc_t s = lv_cmake((float)std::cos(phase_offset), (float)std::sin(phase_offset));
      volk_32fc_s32fc_multiply_32fc(buffer_c, fft_result, s, d_bin_size);
      volk_32fc_x2_add_32fc(buffer_c, buffer_c, &fft_result[d_fft_size-d_bin_size], d_bin_size);
      volk_32fc_magnitude_32f(buffer, buffer_c, d_bin_size);
      return argmax_32f(buffer, max_val_p); 
    }

    uint16_t
    demod_impl::argmax(gr_complex *fft_result)
    {
      float magsq   = pow(real(fft_result[0]), 2) + pow(imag(fft_result[0]), 2);
      float max_val = magsq;
      uint16_t   max_idx = 0;


      for (uint16_t i = 0; i < d_fft_size; i++)
      {
        magsq = pow(real(fft_result[i]), 2) + pow(imag(fft_result[i]), 2);
        if (magsq > max_val)
        {
          max_idx = i;
          max_val = magsq;
        }
      }

      return max_idx;
    }

    void
    demod_impl::parse_header(pmt::pmt_t dict)
    {
      pmt::pmt_t not_found  = pmt::from_bool(false);

      std::string symbol_id = pmt::symbol_to_string(pmt::dict_ref(dict, pmt::intern("id"), not_found));
      d_header_valid        = pmt::to_bool(pmt::dict_ref(dict, pmt::intern("is_valid"), not_found));
      d_header_received     = true;

      if (d_header_valid)
      {
        d_payload_len       = pmt::to_long(pmt::dict_ref(dict, pmt::intern("payload_len"), not_found));
        d_cr                = pmt::to_long(pmt::dict_ref(dict, pmt::intern("cr"), not_found));
        d_crc               = pmt::to_bool(pmt::dict_ref(dict, pmt::intern("crc"), not_found));
        d_packet_symbol_len = 8 + std::max((4+d_cr)*(int)std::ceil((2.0*d_payload_len-d_sf+7+4*d_crc-5*!d_header)/(d_sf-2*d_ldr)), 0); 

        #if DEBUG >= DEBUG_INFO
          std::cout << "PARSE HEADER" << std::endl;
          std::cout << "id: " << symbol_id << std::endl;
          std::cout << "payload_len: " << int(d_payload_len) << std::endl;
          std::cout << "cr: " << int(d_cr) << std::endl;
          std::cout << "crc: " << int(d_crc) << std::endl;
          std::cout << "packet_symbol_len: " << int(d_packet_symbol_len) << std::endl;
        #endif
      }
    }

    void
    demod_impl::dynamic_compensation(std::vector<uint16_t>& compensated_symbols)
    {
      float modulus   = 4.0;
      float bin_drift = 0;
      float bin_comp  = 0;
      float v         = 0;
      float v_last    = 1;
      for (int i = 0; i < d_symbols.size(); i++)
      {
        v = d_symbols[i];

        bin_drift = gr::lora::fpmod(v - v_last, modulus);

        // compensate bin drift
        if (bin_drift < modulus / 2) bin_comp -= bin_drift;
        else bin_comp -= (bin_drift - modulus);
        bin_comp = d_ldr ? bin_comp : 0;
        v_last = v;
        compensated_symbols.push_back(gr::lora::pmod(round(gr::lora::fpmod(v + bin_comp, d_num_symbols)), d_num_symbols));
      }
    }

    void
    demod_impl::forecast (int noutput_items,
                          gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items * (1 << d_sf) * 2;
    }

    int
    demod_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      if (ninput_items[0] < DEMOD_HISTORY_DEPTH*d_num_samples) return 0;
      const gr_complex *in0 = (const gr_complex *) input_items[0];
      const gr_complex *in  = &in0[(DEMOD_HISTORY_DEPTH-1)*d_num_samples];
      uint32_t  *out    = (uint32_t   *) output_items[0];


      uint32_t num_consumed   = d_num_samples;
      uint32_t max_idx        = 0;
      uint32_t max_idx_sfd    = 0;
      bool         preamble_found = false;
      bool         sfd_found      = false;
      float        max_val        = 0;
      float        max_val_sfd    = 0;

      // Nomenclature:
      //  up_block   == de-chirping buffer to contain upchirp features: the preamble, sync word, and data chirps
      //  down_block == de-chirping buffer to contain downchirp features: the SFD
      gr_complex *up_block   = (gr_complex *)volk_malloc(d_fft_size*sizeof(gr_complex), volk_get_alignment());
      gr_complex *down_block = (gr_complex *)volk_malloc(d_fft_size*sizeof(gr_complex), volk_get_alignment());
      float *fft_res_mag = (float*)volk_malloc(d_fft_size*sizeof(float), volk_get_alignment());
      float *fft_res_add = (float*)volk_malloc(d_bin_size*sizeof(float), volk_get_alignment());
      gr_complex *fft_res_add_c = (gr_complex*)volk_malloc(d_bin_size*sizeof(gr_complex), volk_get_alignment());

      if (up_block == NULL || down_block == NULL ||
          fft_res_mag == NULL || fft_res_add == NULL || fft_res_add_c == NULL)
      {
        std::cerr << "Unable to allocate processing buffer!" << std::endl;
      }

      // Dechirp the incoming signal
      volk_32fc_x2_multiply_32fc(up_block, in, &d_downchirp[0], d_num_samples);

      if (d_state == S_SFD_SYNC)
      {
        volk_32fc_x2_multiply_32fc(down_block, in, &d_upchirp[0], d_num_samples);
      }

      // Enable to write IQ to disk for debugging
      #if DUMP_IQ
        f_up_windowless.write((const char*)&up_block[0], d_num_samples*sizeof(gr_complex));
      #endif

      // Windowing
      // volk_32fc_32f_multiply_32fc(up_block, up_block, &d_window[0], d_num_samples);

      #if DUMP_IQ
        if (d_state != S_SFD_SYNC) f_down.write((const char*)&down_block[0], d_num_samples*sizeof(gr_complex));
        f_up.write((const char*)&up_block[0], d_num_samples*sizeof(gr_complex));
      #endif

      // Preamble and Data FFT
      // If d_fft_size_factor is greater than 1, the rest of the sample buffer will be zeroed out and blend into the window
      memset(d_fft->get_inbuf(),            0, d_fft_size*sizeof(gr_complex));
      memcpy(d_fft->get_inbuf(), &up_block[0], d_num_samples*sizeof(gr_complex));
      d_fft->execute();
      #if DUMP_IQ
        f_fft.write((const char*)d_fft->get_outbuf(), d_fft_size*sizeof(gr_complex));
      #endif

      // Take argmax of returned FFT (similar to MFSK demod)
      max_idx = search_fft_peak(d_fft->get_outbuf(), fft_res_mag, fft_res_add, fft_res_add_c, &max_val);

      d_argmax_history.insert(d_argmax_history.begin(), max_idx);

      if (d_argmax_history.size() > REQUIRED_PREAMBLE_CHIRPS)
      {
        d_argmax_history.pop_back();
      }

      switch (d_state) {
      case S_RESET:
      {
        d_overlaps = OVERLAP_DEFAULT;
        d_offset = 0;
        d_symbols.clear();
        d_argmax_history.clear();
        d_sfd_history.clear();
        d_sync_recovery_counter = 0;
        d_header_received = false;

        d_state = S_PREFILL;

        #if DEBUG >= DEBUG_INFO
          std::cout << "Next state: S_PREFILL" << std::endl;
        #endif

        break;
      }



      case S_PREFILL:
      {
        if (d_argmax_history.size() >= REQUIRED_PREAMBLE_CHIRPS)
        {
          d_state = S_DETECT_PREAMBLE;

          #if DEBUG >= DEBUG_INFO
            std::cout << "Next state: S_DETECT_PREAMBLE" << std::endl;
          #endif
        }
        break;
      }



      // Looks for the same symbol appearing consecutively, signifying the LoRa preamble
      case S_DETECT_PREAMBLE:
      {
        d_preamble_idx = d_argmax_history[0];

        #if DEBUG >= DEBUG_VERBOSE
          std::cout << "PREAMBLE " << d_argmax_history[0] << std::endl;
        #endif

        // Check for discontinuities that exceed some tolerance
        preamble_found = true;
        for (int i = 1; i < REQUIRED_PREAMBLE_CHIRPS; i++)
        {
          uint32_t dis = gr::lora::pmod(int(d_preamble_idx) - int(d_argmax_history[i]), d_bin_size);
          if (dis > d_preamble_drift_max && dis < d_bin_size-d_preamble_drift_max)
          {
            preamble_found = false;
          }
        }

        // Advance to SFD/sync discovery if a contiguous preamble is found
        if (preamble_found)
        {
          d_state = S_SFD_SYNC;

          // move preamble peak to bin zero
          num_consumed = d_num_samples - d_p*d_preamble_idx/d_fft_size_factor;

          #if DEBUG >= DEBUG_INFO
            std::cout << "Next state: S_SFD_SYNC" << std::endl;
          #endif
        }
        break;
      }



      // Accurately synchronize to the SFD by computing overlapping FFTs of the downchirp/SFD IQ stream
      // Effectively increases FFT's time-based resolution, allowing for a better sync
      case S_SFD_SYNC:
      {
        d_overlaps = OVERLAP_FACTOR;

        // Recover if the SFD is missed, or if we wind up in this state erroneously (false positive on preamble)
        if (d_sync_recovery_counter++ > DEMOD_SYNC_RECOVERY_COUNT)
        {
          d_state = S_RESET;
          d_overlaps = OVERLAP_DEFAULT;

          #if DEBUG >= DEBUG_INFO
            std::cout << "Bailing out of sync loop"   << std::endl;
            std::cout << "Next state: S_RESET" << std::endl;
          #endif
        }

        // Dechirp
        volk_32fc_x2_multiply_32fc(down_block, in, &d_upchirp[0], d_num_samples);

        // Enable to write out overlapped chirps to disk for debugging
        #if DUMP_IQ
          f_down.write((const char*)&down_block[0], d_num_samples*sizeof(gr_complex));
        #endif

        memset(d_fft->get_inbuf(),          0, d_fft_size*sizeof(gr_complex));
        memcpy(d_fft->get_inbuf(), down_block, d_num_samples*sizeof(gr_complex));
        d_fft->execute();

        // Take argmax of downchirp FFT
        max_idx_sfd = search_fft_peak(d_fft->get_outbuf(), fft_res_mag, fft_res_add, fft_res_add_c, &max_val_sfd);

        // If SFD is detected
        if (max_val_sfd > max_val)
        {
          int idx = max_idx_sfd;
          if (max_idx_sfd > d_bin_size / 2) {
            idx = max_idx_sfd - d_bin_size;
          }
          num_consumed = (int)round(2.25*d_num_samples + d_p*idx/2.0/d_fft_size_factor);

          // refine CFO
          volk_32fc_x2_multiply_32fc(up_block, 
            &in0[(int)round((DEMOD_HISTORY_DEPTH-1-5.25)*d_num_samples) + num_consumed],
            &d_downchirp[0], d_num_samples);
          memset(d_fft->get_inbuf(),        0, d_fft_size*sizeof(gr_complex));
          memcpy(d_fft->get_inbuf(), up_block, d_num_samples*sizeof(gr_complex));
          d_fft->execute();
          d_cfo = (float)search_fft_peak(d_fft->get_outbuf(), fft_res_mag, fft_res_add, fft_res_add_c, &max_val);

          d_state = S_READ_HEADER;

          #if DEBUG >= DEBUG_INFO
            std::cout << "Next state: S_READ_HEADER" << std::endl;
            std::cout << "CFO: " << d_cfo << std::endl;
          #endif

          break;
        }

        break;
      }



      case S_READ_HEADER:
      {
        /* Preamble + modulo operation normalizes the symbols about the preamble; preamble symbol == value 0
         * Dividing by d_fft_size_factor reduces symbols to [0:(2**sf)-1] range
         * Dividing by 4 to further reduce symbol set to [0:(2**(sf-2)-1)], since header is sent at SF-2
         */
        float bin_idx = gr::lora::fpmod((max_idx - d_cfo)/(float)d_fft_size_factor, d_num_symbols);
        #if DEBUG >= DEBUG_INFO
          std::cout << "MIDX: " << bin_idx << ", MV: " << max_val << std::endl;
        #endif
        d_symbols.push_back( bin_idx );

        if (d_symbols.size() == 8)   // Symbols [0:7] contain 2**(SF-2) bits/symbol, symbols [8:] have the full 2**(SF) bits
        {
          std::vector<uint16_t> compensated_symbols;
          dynamic_compensation(compensated_symbols);
          pmt::pmt_t header   = pmt::init_u16vector(compensated_symbols.size(), compensated_symbols);
          pmt::pmt_t dict     = pmt::make_dict();
          dict = pmt::dict_add(dict, pmt::intern("id"), pmt::intern("header"));
          pmt::pmt_t msg_pair = pmt::cons(dict, header);
          message_port_pub(d_out_port, msg_pair); 
        }
        else if (d_symbols.size() > 8)
        {
          if (!d_header || (d_header && d_header_received))
          {
            if (d_header_received && !d_header_valid)
            {
              d_state = S_RESET;

              #if DEBUG >= DEBUG_INFO
                std::cout << "Invalid header received" << std::endl;
              #endif
            }
            else
            {
              d_state = S_READ_PAYLOAD;

              #if DEBUG >= DEBUG_INFO
                std::cout << "Next state: S_READ_PAYLOAD" << std::endl;
              #endif
            }
          }
          // wait for header parsing
        }

        break;
      }


      case S_READ_PAYLOAD:
      {
        if (d_symbols.size() >= d_packet_symbol_len)
        {
          d_state = S_OUT;

          #if DEBUG >= DEBUG_INFO
            std::cout << "Next state: S_OUT" << std::endl;
          #endif
        }
        else {
          /* Preamble + modulo operation normalizes the symbols about the preamble; preamble symbol == value 0
          * Dividing by d_fft_size_factor reduces symbols to [0:(2**sf)-1] range
          */
          float bin_idx = gr::lora::fpmod((max_idx - d_cfo)/(float)d_fft_size_factor, d_num_symbols);
          #if DEBUG >= DEBUG_INFO
            std::cout << "MIDX: " << bin_idx << ", MV: " << max_val << std::endl;
          #endif
          d_symbols.push_back( bin_idx );
        }

        break;
      }



      // Emit a PDU to the decoder
      case S_OUT:
      {
        std::vector<uint16_t> compensated_symbols;
        dynamic_compensation(compensated_symbols);
        pmt::pmt_t output = pmt::init_u16vector(compensated_symbols.size(), compensated_symbols);
        pmt::pmt_t dict = pmt::make_dict();
        dict = pmt::dict_add(dict, pmt::intern("id"), pmt::intern("packet"));
        pmt::pmt_t msg_pair = pmt::cons(dict, output);
        message_port_pub(d_out_port, msg_pair);

        d_state = S_RESET;
        #if DEBUG >= DEBUG_INFO
          std::cout << "Next state: S_RESET" << std::endl;
          std::cout << "d_symbols size: " << d_symbols.size() << std::endl;
          std::cout << "compensated_symbols: ";
          for (auto i: compensated_symbols) {
            std::cout << i << " ";
          }
          std::cout << std::endl;
        #endif

        break;
      }



      default:
        break;
      }

      #if DUMP_IQ
        f_raw.write((const char*)&in[0], num_consumed*sizeof(gr_complex));
      #endif

      consume_each (num_consumed);

      volk_free(down_block);
      volk_free(up_block);
      volk_free(fft_res_mag);
      volk_free(fft_res_add);
      volk_free(fft_res_add_c);

      return noutput_items;
    }

  } /* namespace lora */
} /* namespace gr */

