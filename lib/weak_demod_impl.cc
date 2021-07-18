/* -*- c++ -*- */
/*
 * Copyright 2021 jkadbear.
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
#include "weak_demod_impl.h"

#define DEBUG_OFF     0
#define DEBUG_INFO    1
#define DEBUG_VERBOSE 2
#define DEBUG         DEBUG_OFF

#define OVERLAP_DEFAULT 1
#define OVERLAP_FACTOR  16

namespace gr {
  namespace lora {

    weak_demod::sptr
    weak_demod::make( uint8_t   spreading_factor,
                 bool      header,
                 uint8_t   payload_len,
                 uint8_t   cr,
                 bool      crc,
                 bool      low_data_rate,
                 uint32_t  sym_num,
                 float     beta,
                 uint16_t  fft_factor,
                 uint8_t   peak_search_algorithm,
                 uint16_t  peak_search_phase_k,
                 float     fs_bw_ratio)
    {
      return gnuradio::get_initial_sptr
        (new weak_demod_impl(spreading_factor, header, payload_len, cr, crc, low_data_rate, sym_num, beta, fft_factor, peak_search_algorithm, peak_search_phase_k, fs_bw_ratio));
    }


    /*
     * The private constructor
     */
    weak_demod_impl::weak_demod_impl(uint8_t   spreading_factor,
                            bool      header,
                            uint8_t   payload_len,
                            uint8_t   cr,
                            bool      crc,
                            bool      low_data_rate,
                            uint32_t  sym_num,
                            float     beta,
                            uint16_t  fft_factor,
                            uint8_t   peak_search_algorithm,
                            uint16_t  peak_search_phase_k,
                            float     fs_bw_ratio)
      : gr::block("weak_demod",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(0, 0, 0)),
        d_sf(spreading_factor),
        d_header(header),
        d_payload_len(payload_len),
        d_cr(cr),
        d_crc(crc),
        d_ldr(low_data_rate),
        d_sym_num(sym_num),
        d_beta(beta),
        d_fft_size_factor(fft_factor),
        d_p((int)fs_bw_ratio),
        f_raw("raw.out", std::ios::out),
        f_fft("fft.out", std::ios::out),
        f_down("down.out", std::ios::out),
        d_peak_search_algorithm(peak_search_algorithm),
        d_peak_search_phase_k(peak_search_phase_k)
    {
      assert((d_sf > 5) && (d_sf < 13));
      if (d_sf == 6) assert(!header);
      assert(d_fft_size_factor > 0);
      assert(((int)fs_bw_ratio) == fs_bw_ratio);

      // if (!header) // implicit header mode
      // {
      //   // calculate the total number of symbols in a packet
      //   d_packet_symbol_len = 8 + std::max((4+d_cr)*(int)std::ceil((2.0*d_payload_len-d_sf+7+4*d_crc-5*!d_header)/(d_sf-2*d_ldr)), 0);
      // }

      // d_header_port = pmt::mp("header");
      // message_port_register_in(d_header_port);
      d_out_port = pmt::mp("out");
      message_port_register_out(d_out_port);

      // set_msg_handler(d_header_port, [this](pmt::pmt_t msg) { this->parse_header(msg); });
      // set_msg_handler(d_header_port, boost::bind(&weak_demod_impl::parse_header, this, _1));

      d_state = WS_RESET;

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

      set_history(WEAK_DEMOD_BUFFER_SIZE*d_num_samples);  // Sync is 2.25 chirp periods long
    }

    /*
     * Our virtual destructor.
     */
    weak_demod_impl::~weak_demod_impl()
    {
    }

    void
    weak_demod_impl::parse_header(pmt::pmt_t dict)
    {
    }

    void
    weak_demod_impl::dechirp(bool is_up,
                        const gr_complex *in,
                        gr_complex *block,
                        float *fft_mag,
                        float *fft_add)
    {
      if (is_up) {
        volk_32fc_x2_multiply_32fc(block, in, &d_downchirp[0], d_num_samples);
      }
      else {
        volk_32fc_x2_multiply_32fc(block, in, &d_upchirp[0], d_num_samples);
      }

      memset(d_fft->get_inbuf(),            0, d_fft_size*sizeof(gr_complex));
      memcpy(d_fft->get_inbuf(), &block[0], d_num_samples*sizeof(gr_complex));
      d_fft->execute();
      volk_32fc_magnitude_32f(fft_mag, d_fft->get_outbuf(), d_fft_size);
      volk_32f_x2_add_32f(fft_add, fft_mag, &fft_mag[d_fft_size-d_bin_size], d_bin_size);
      #if DEBUG >= DEBUG_VERBOSE
        float max_val = 0;
        uint32_t max_idx = gr::lora::argmax_32f(fft_add, &max_val, d_bin_size);
        std::cout << "[dechirp] max_idx: " << max_idx << ", max_val: " << max_val << std::endl;
      #endif
    }

    uint32_t
    weak_demod_impl::search_fft_peak(bool is_up,
                                const gr_complex *in,
                                gr_complex *block1,
                                gr_complex *block2,
                                float *fft_mag1,
                                float *fft_mag2,
                                float *fft_add1,
                                float *fft_add2,
                                float *p_max_val)
    {
      // Dechirp the incoming signal
      dechirp(is_up, in, block1, fft_mag1, fft_add1);
      dechirp(is_up, &in[d_num_samples], block2, fft_mag2, fft_add2);

      if (*p_max_val == 999) {
        f_fft.write((const char*)&fft_add1[0], d_bin_size*sizeof(float));
        f_down.write((const char*)&fft_add2[0], d_bin_size*sizeof(float));
      }

      volk_32f_x2_add_32f(fft_add1, fft_add1, fft_add2, d_bin_size);
      return gr::lora::argmax_32f(fft_add1, p_max_val, d_bin_size);
    }

    void
    weak_demod_impl::dynamic_compensation(std::vector<uint16_t>& compensated_symbols)
    {
      float modulus   = d_ldr ? 4.0 : 1.0;
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
        v_last = v;
        compensated_symbols.push_back(gr::lora::pmod(round(gr::lora::fpmod(v + bin_comp, d_num_symbols)), d_num_symbols));
      }
    }

    void
    weak_demod_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items * 2 * d_num_samples;
    }

    int
    weak_demod_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      if (ninput_items[0] < WEAK_DEMOD_BUFFER_SIZE*d_num_samples) return 0;
      const gr_complex *in0 = (const gr_complex *) input_items[0];
      const gr_complex *in  = &in0[WEAK_DEMOD_HISTORY*d_num_samples];
      uint32_t *out = (uint32_t *) output_items[0];
      bool preamble_found = false;

      uint32_t max_idx = 0;
      float max_val = 0;

      gr_complex *block1 = (gr_complex *)volk_malloc(d_fft_size*sizeof(gr_complex), volk_get_alignment());
      gr_complex *block2 = (gr_complex *)volk_malloc(d_fft_size*sizeof(gr_complex), volk_get_alignment());
      float *fft_mag1 = (float*)volk_malloc(d_fft_size*sizeof(float), volk_get_alignment());
      float *fft_mag2 = (float*)volk_malloc(d_fft_size*sizeof(float), volk_get_alignment());
      float *fft_add1 = (float*)volk_malloc(d_bin_size*sizeof(float), volk_get_alignment());
      float *fft_add2 = (float*)volk_malloc(d_bin_size*sizeof(float), volk_get_alignment());

      if (block1 == NULL ||
          block2 == NULL ||
          fft_mag1 == NULL ||
          fft_mag2 == NULL ||
          fft_add1 == NULL ||
          fft_add2 == NULL
      )
      {
        std::cerr << "Unable to allocate processing buffer!" << std::endl;
      }

      max_idx = search_fft_peak(true, in, block1, block2, fft_mag1, fft_mag2, fft_add1, fft_add2, &max_val);

      uint32_t num_consumed = d_num_samples;
      static uint32_t sym_cnt = 0;
      static uint32_t idx_cnt = 0;

      #if DEBUG >= DEBUG_VERBOSE
        std::cout << "idx: " << max_idx << ", val: " << max_val << std::endl;
      #endif

      if (max_val > 0) {
        d_argmax_history.insert(d_argmax_history.begin(), max_idx);

        if (d_argmax_history.size() > WEAK_REQUIRED_PREAMBLE_CHIRPS)
        {
          d_argmax_history.pop_back();
        }
      }

      switch (d_state) {
      case WS_RESET:
      {
        d_overlaps = OVERLAP_DEFAULT;
        d_offset = 0;
        d_symbols.clear();
        d_argmax_history.clear();
        d_sfd_history.clear();
        d_sync_recovery_counter = 0;
        d_header_received = false;

        d_state = WS_PREFILL;

        sym_cnt = 0;

        #if DEBUG >= DEBUG_INFO
          std::cout << "Next state: WS_PREFILL" << std::endl;
        #endif

        break;
      }

      case WS_PREFILL:
      {
        if (d_argmax_history.size() == WEAK_REQUIRED_PREAMBLE_CHIRPS)
        {
          d_state = WS_DETECT_PREAMBLE;

          #if DEBUG >= DEBUG_INFO
            std::cout << "Next state: WS_DETECT_PREAMBLE" << std::endl;
          #endif
        }
        break;
      }

      case WS_DETECT_PREAMBLE:
      {
        d_preamble_idx = d_argmax_history[0];


        // Check for discontinuities that exceed some tolerance
        preamble_found = true;
        for (int i = 1; i < WEAK_REQUIRED_PREAMBLE_CHIRPS; i++)
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
          #if DEBUG >= DEBUG_VERBOSE
            std::cout << "PREAMBLE: ";
            for (auto idx: d_argmax_history) {
              std::cout << idx << ",";
            }
            std::cout << std::endl;
          #endif

          d_state = WS_SFD_SYNC;

          // move preamble peak to bin zero
          num_consumed = d_num_samples - d_p*d_preamble_idx/d_fft_size_factor;

          #if DEBUG >= DEBUG_INFO
            std::cout << "Next state: WS_SFD_SYNC" << std::endl;
          #endif
        }
        break;
      }

      case WS_SFD_SYNC:
      {
        if (d_sync_recovery_counter++ > WEAK_DEMOD_SYNC_RECOVERY_COUNT)
        {
          d_state = WS_RESET;

          #if DEBUG >= DEBUG_INFO
            std::cout << "Bailing out of sync loop"   << std::endl;
            std::cout << "Next state: WS_RESET" << std::endl;
          #endif
        }

        float max_val_up[2] = {max_val, 0};
        float max_val_down[2] = {0};
        uint32_t max_idx_up[2] = {max_idx, 0};
        uint32_t max_idx_down[2] = {0};
        max_idx_up[1] = search_fft_peak(true, &in[d_num_samples], block1, block2, fft_mag1, fft_mag2, fft_add1, fft_add2, &max_val_up[1]);
        max_idx_down[0] = search_fft_peak(false, in, block1, block2, fft_mag1, fft_mag2, fft_add1, fft_add2, &max_val_down[0]);
        max_idx_down[1] = search_fft_peak(false, &in[d_num_samples], block1, block2, fft_mag1, fft_mag2, fft_add1, fft_add2, &max_val_down[1]);

        #if DEBUG >= DEBUG_VERBOSE
          std::cout << "[SYNC] max_idx_up[0]: " << max_idx_up[0] << ", max_idx_up[1]: " << max_idx_up[1] << std::endl;
          std::cout << "[SYNC] max_idx_down[0]: " << max_idx_down[0] << ", max_idx_down[1]: " << max_idx_down[1] << std::endl;
        #endif

        float tmp = 0;
        uint32_t i = gr::lora::argmax_32f(max_val_down, &tmp, 2);
        // If SFD is detected
        if (i == 0 && max_val_down[i] > max_val_up[i]) {
            int32_t offset = (max_idx_down[i] > d_bin_size / 2) ? ((int32_t)max_idx_down[i] - d_bin_size) : max_idx_down[i];
            num_consumed = (int)round((2.25+i)*d_num_samples + d_p*offset/2.0/d_fft_size_factor);
            // refine CFO
            float max_val_cfo = 999;
            d_cfo = (float) search_fft_peak(true, &in0[(int)round((WEAK_DEMOD_HISTORY-6.25)*d_num_samples + num_consumed)], block1, block2, fft_mag1, fft_mag2, fft_add1, fft_add2, &max_val_cfo);

            d_state = WS_READ_PAYLOAD;

            #if DEBUG >= DEBUG_INFO
              std::cout << "Next state: WS_READ_PAYLOAD" << std::endl;
              std::cout << "CFO: " << d_cfo << ", CFO peak height: " << max_val_cfo << std::endl;
              std::cout << "num_consumed: " << num_consumed << ", offset: " << offset << std::endl;
            #endif
        }

        break;
      }

      case WS_READ_PAYLOAD:
      {
        if (idx_cnt == 0) idx_cnt = 1;
        // TODO
        if (d_symbols.size() >= d_sym_num) {
          d_state = WS_OUT;
          #if DEBUG >= DEBUG_INFO
            std::cout << "Next state: WS_OUT" << std::endl;
          #endif
        }
        else {
          float bin_idx = gr::lora::fpmod((max_idx - d_cfo)/(float)d_fft_size_factor, d_num_symbols);
          if (sym_cnt < 2) {
            num_consumed = 2*d_num_samples;
            d_symbols.push_back( bin_idx );
            #if DEBUG >= DEBUG_INFO
              std::cout << "MIDX: " << bin_idx << ", MV: " << max_val << std::endl;
            #endif
          }
          else if (sym_cnt == 2) {
            // skip the checksum of header symbols
            num_consumed = 4*d_num_samples;
          }
          else {
            if ((sym_cnt-3) % 3 != 2) {
              num_consumed = 2*d_num_samples;
              d_symbols.push_back( bin_idx );

              #if DEBUG >= DEBUG_INFO
              std::cout << "MIDX: " << bin_idx << ", max_idx: " << max_idx << ", MV: " << max_val << std::endl;
              #endif
            }
            else {
              #if DEBUG >= DEBUG_INFO
                std::cout << "else MIDX: " << bin_idx << ", MV: " << max_val << std::endl;
              #endif
              num_consumed = d_num_samples;
            }
          }
          sym_cnt++;
        }
        break;
      }

      // Emit a PDU to the decoder
      case WS_OUT:
      {
        std::vector<uint16_t> compensated_symbols;
        dynamic_compensation(compensated_symbols);
        pmt::pmt_t output = pmt::init_u16vector(compensated_symbols.size(), compensated_symbols);
        pmt::pmt_t dict = pmt::make_dict();
        dict = pmt::dict_add(dict, pmt::intern("id"), pmt::intern("packet"));
        pmt::pmt_t msg_pair = pmt::cons(dict, output);
        message_port_pub(d_out_port, msg_pair);

        d_state = WS_RESET;
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

      consume_each (num_consumed);
      if (idx_cnt > 0) idx_cnt += num_consumed;

      volk_free(block1);
      volk_free(block2);
      volk_free(fft_mag1);
      volk_free(fft_mag2);
      volk_free(fft_add1);
      volk_free(fft_add2);

      return noutput_items;
    }

  } /* namespace lora */
} /* namespace gr */

