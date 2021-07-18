/* -*- c++ -*- */
/* 
 * Copyright 2020 jkadbear.
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
#include "pyramid_demod_impl.h"

#define DEBUG_OFF                0
#define DEBUG_INFO               1
#define DEBUG_VERBOSE            2
#define DEBUG_VERBOSE_VERBOSE    3
#define DEBUG                    DEBUG_OFF

#define DUMP_IQ       0

#define OVERLAP_FACTOR  8

namespace gr {
  namespace lora {

    pyramid_demod::sptr
    pyramid_demod::make( uint8_t spreading_factor,
                         bool  low_data_rate,
                         float beta,
                         uint16_t fft_factor,
                         float threshold,
                         float fs_bw_ratio)
    {
      return gnuradio::get_initial_sptr
        (new pyramid_demod_impl(spreading_factor, low_data_rate, beta, fft_factor, threshold, fs_bw_ratio));
    }

    /*
     * The private constructor
     */
    pyramid_demod_impl::pyramid_demod_impl( uint16_t spreading_factor,
                            bool  low_data_rate,
                            float beta,
                            uint16_t fft_factor,
                            float threshold,
                            float fs_bw_ratio)
      : gr::block("pyramid_demod",
              gr::io_signature::make(1, 1, sizeof(gr_complex)),
              gr::io_signature::make(0, 0, 0)),
        f_raw("raw.out", std::ios::out),
        f_fft("fft.out", std::ios::out),
        f_up_windowless("up_windowless.out", std::ios::out),
        f_up("up.out", std::ios::out),
        f_down("down.out", std::ios::out),
        d_sf(spreading_factor),
        d_ldr(low_data_rate),
        d_beta(beta),
        d_fft_size_factor(fft_factor),
        d_threshold(threshold)
    {
      assert((d_sf > 5) && (d_sf < 13));
      if (d_sf == 6) assert(!header);
      assert(d_fft_size_factor > 0);
      assert(((int)fs_bw_ratio) == fs_bw_ratio);
      d_p = (int) fs_bw_ratio;

      d_header_port = pmt::mp("header");
      message_port_register_in(d_header_port);
      d_out_port = pmt::mp("out");
      message_port_register_out(d_out_port);

      // set_msg_handler(d_header_port, [this](pmt::pmt_t msg) { this->parse_header(msg); });

      d_num_symbols = (1 << d_sf);
      d_num_samples = d_p*d_num_symbols;
      d_bin_size = d_fft_size_factor*d_num_symbols;
      d_fft_size = d_fft_size_factor*d_num_samples;
      d_fft = new fft::fft_complex(d_fft_size, true, 1);
      d_overlaps = OVERLAP_FACTOR;
      d_ttl = 6*d_overlaps; // MAGIC
      d_offset = 0;

      d_window = fft::window::build(fft::window::WIN_KAISER, d_num_samples, d_beta);

      d_ts_ref  = 0;
      d_bin_ref = 0;
      d_bin_tolerance = (d_ldr ? d_fft_size_factor * 2 : d_fft_size_factor / 2);

      // Create local chirp tables.  Each table is 2 chirps long to allow memcpying from arbitrary offsets.
      for (int i = 0; i < d_num_samples; i++) {
        double phase = M_PI/d_p*(i-i*i/(float)d_num_samples);
        d_downchirp.push_back(gr_complex(std::polar(1.0, phase)));
        d_upchirp.push_back(gr_complex(std::polar(1.0, -phase)));
      }

      uint16_t track_size = 1000; // MAGIC
      d_num_preamble = 6; // MAGIC
      d_track.resize(track_size);
      for (auto & v : d_track)
      {
        v.reserve(d_overlaps * (d_num_preamble + 2));
      }
      d_bin_track_id_list.reserve(track_size);
      for (uint16_t i = 0; i < track_size; i++)
      {
        d_track_id_pool.push_back(i);
      }

      uint16_t packet_id_size = 40; // MAGIC
      d_packet.resize(packet_id_size);
      d_packet_state_list.reserve(packet_id_size);
      for (uint16_t i = 0; i < packet_id_size; i++)
      {
        d_packet_id_pool.push_back(i);
      }

      set_history(PY_DEMOD_HISTORY_DEPTH*d_num_samples);  // Sync is 2.25 symbols long
    }

    /*
     * Our virtual destructor.
     */
    pyramid_demod_impl::~pyramid_demod_impl()
    {
      delete d_fft;
    }

    uint32_t
    pyramid_demod_impl::argmax_32f(float *fft_result, 
                       bool update_squelch, float *max_val_p)
    {
      float mag   = abs(fft_result[0]);
      float max_val = mag;
      uint32_t max_idx = 0;

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

    uint16_t
    pyramid_demod_impl::argmax(gr_complex *fft_result, 
                       bool update_squelch)
    {
      float magsq   = pow(real(fft_result[0]), 2) + pow(imag(fft_result[0]), 2);
      float max_val = magsq;
      uint16_t max_idx = 0;


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

    float
    pyramid_demod_impl::get_dis(uint32_t ts1, float h1, uint32_t ts2, float h2)
    {
      float dis = gr::lora::pmod(ts1 - ts2, d_num_samples) / (float) d_num_samples;
      // In definition above, "dis->0" and "dis->1" both mean their timestamp differences to the actual timestamp is small
      // therefore we transform dis to let "dis->1" represent larger difference
      dis = (dis > 0.5) ? (1 - dis) * 2 : dis * 2;
      dis += std::abs(h1 - h2) / h2;
      return dis;
    }

    // void
    // pyramid_demod_impl::parse_header(pmt::pmt_t dict)
    // {
    //   pmt::pmt_t not_found = pmt::from_bool(false);

    //   std::string symbol_id = pmt::symbol_to_string(pmt::dict_ref(dict, pmt::intern("id"), not_found));
    //   d_header_valid        = pmt::to_bool(pmt::dict_ref(dict, pmt::intern("is_valid"), not_found));
    //   d_header_received     = true;

    //   if (d_header_valid)
    //   {
    //     d_payload_len       = pmt::to_long(pmt::dict_ref(dict, pmt::intern("payload_len"), not_found));
    //     d_cr                = pmt::to_long(pmt::dict_ref(dict, pmt::intern("cr"), not_found));
    //     d_crc               = pmt::to_bool(pmt::dict_ref(dict, pmt::intern("crc"), not_found));
    //     d_packet_symbol_len = 8 + std::max((4+d_cr)*(int)std::ceil((2.0*d_payload_len-d_sf+7+4*d_crc-5*!d_header)/(d_sf-2*d_ldr)), 0); 

    //     #if DEBUG >= DEBUG_INFO
    //       std::cout << "PARSE HEADER" << std::endl;
    //       std::cout << "id: " << symbol_id << std::endl;
    //       std::cout << "payload_len: " << int(d_payload_len) << std::endl;
    //       std::cout << "cr: " << int(d_cr) << std::endl;
    //       std::cout << "crc: " << int(d_crc) << std::endl;
    //       std::cout << "packet_symbol_len: " << int(d_packet_symbol_len) << std::endl;
    //     #endif
    //   }
    // }

    void
    pyramid_demod_impl::find_and_add_peak(float *fft_add, float *fft_add_w, float *fft_mag)
    {
      for (uint32_t i = 0; i < d_bin_size; i++)
      {
        // find peak: search local maximum larger than d_threshold
        uint32_t l_idx = gr::lora::pmod(i-1, d_bin_size);
        uint32_t r_idx = gr::lora::pmod(i+1, d_bin_size);
        if(fft_add_w[i] > d_threshold && fft_add_w[i] > fft_add_w[l_idx] && fft_add_w[i] > fft_add_w[r_idx])
        {
          // this is a peak, insert it into the peak track
          uint32_t cur_bin = gr::lora::pmod(d_bin_size + i - d_bin_ref, d_bin_size);
          bool found = false;
          uint16_t track_id;
          for (auto & bt: d_bin_track_id_list)
          {
            uint32_t dis = gr::lora::pmod(d_bin_size + cur_bin - bt.bin, d_bin_size);
            #if DEBUG >= DEBUG_VERBOSE_VERBOSE
              std::cout << "dis: " << dis << ", bt.bin: " << bt.bin << std::endl;
            #endif
            // Abs(current_bin - target_bin) < d_bin_tolerance
            if (dis <= d_bin_tolerance || dis >= d_bin_size - d_bin_tolerance)
            {
              found = true;
              track_id = bt.track_id;
              bt.updated = true;
              break;
            }
          }
          if (!found)
          {
            if (d_track_id_pool.empty())
            {
              std::cerr << "Current threshold is low! Increase the threshold or track size" << std::endl;
              exit(-1);
            }
            track_id = d_track_id_pool.front();
            d_track_id_pool.pop_front();
            d_bin_track_id_list.push_back(bin_track_id(cur_bin, track_id, true));
          }

          #if DEBUG >= DEBUG_VERBOSE_VERBOSE
            std::cout << "track id: " << track_id << ", track size: " << d_track[track_id].size() << ", track id pool size: " << d_track_id_pool.size() << ", bin: " << i << ", ref bin: " << d_bin_ref << ", peak height: " << fft_add[l_idx] << " " << fft_add[i] << " " << fft_add[r_idx] << " " << fft_add_w[i] << std::endl;
          #endif
          d_track[track_id].push_back(peak(d_ts_ref, i, fft_add[i], std::max(fft_mag[i], fft_mag[d_fft_size-d_bin_size+i])));
        }
      }
    }

    void
    pyramid_demod_impl::get_apex(std::vector<peak> &track, peak & pk, bool is_preamble)
    {
      float h[OVERLAP_FACTOR*2 + 10];
      for (int i = 0; i < track.size(); i++)
      {
        h[i] = is_preamble ? track[i].h_single : track[i].h;
      }
      float max_h = h[0];
      uint16_t idx = 0;
      for (uint16_t i = 1; i < track.size(); i++)
      {
        if (h[i] > max_h)
        {
          max_h = h[i];
          idx = i;
        }
      }
#if APEX_ALGORITHM == APEX_ALGORITHM_SEGMENT
      pk.ts  = track[idx].ts;
      pk.bin = track[idx].bin;
      pk.h   = h[idx];
#elif APEX_ALGORITHM == APEX_ALGORITHM_LINEAR_REGRESSION
      // linear regression method requires at least 4 points
      if (idx < 1 || idx > track.size() - 2 || track.size() < 4)
      {
        pk.ts  = track[idx].ts;
        pk.bin = track[idx].bin;
        pk.h   = h[idx];
      }
      else
      {
        uint16_t l_idx = h[idx-1] > h[idx+1] ? idx-1 : idx;
        float k1, b1, k2, b2;
        linear_regression(h, 0, l_idx, &k1, &b1);
        linear_regression(h, l_idx+1, track.size()-1, &k2, &b2);
        float x = -(b2-b1)/(k2-k1);
        pk.ts = gr::lora::pmod(track[l_idx].ts + round((x-l_idx)*d_num_samples/d_overlaps), TIMESTAMP_MOD);
        pk.h = k1 * x + b1;
        // std::cout << "k1: " << k1 << ", b1: " << b1 << ", k2: " << k2 << ", b2: " << b2 << ", x: " << x << std::endl;
        pk.bin = gr::lora::pmod(track[l_idx].bin + round((x-l_idx)*d_bin_size/d_overlaps), d_bin_size);
      }
#endif
    }

    symbol_type
    pyramid_demod_impl::get_central_peak(uint16_t track_id, peak & pk)
    {
      auto track = d_track[track_id];
      uint16_t len = track.size();

      #if DEBUG >= DEBUG_VERBOSE
        std::cout << "track id: " << track_id << ", peak height: ";
        for (auto & pk : track)
        {
          std::cout << pk.h << ", ";
        }
        std::cout << std::endl;
      #endif
      if (len >= d_overlaps*(d_num_preamble-1) + 2)
      {
        // // preamble
        // uint16_t l_idx = len/2 - d_overlaps*(d_num_preamble-1)/2;
        // uint16_t r_idx = (len-1)/2 + d_overlaps*(d_num_preamble-1)/2;
        // if (track[l_idx].h > track[r_idx].h)
        // {
        //   pk.ts  = track[l_idx].ts + d_num_samples/4 + (d_num_preamble-1)*d_num_samples;
        //   pk.bin = track[l_idx].bin;
        // }
        // else 
        // {
        //   pk.ts  = track[r_idx].ts + d_num_samples/4;
        //   pk.bin = track[r_idx].bin;
        // }

        // extract the last chirp of preamble
        // determine the timestamp of preamble through the single peak trajectory
        uint16_t r_idx = track.size() - d_overlaps;
        float max_h = -1;
        for (int i = r_idx; i < track.size(); i++)
        {
          if (track[i].h > max_h)
          {
            max_h = track[i].h;
            r_idx = i;
          }
        }
        uint16_t start_idx = r_idx;
        for (; start_idx > r_idx - d_overlaps/2; start_idx--)
        {
          if (track[start_idx-1].h_single > track[start_idx].h_single
              || track[start_idx].h_single < d_threshold)
            break;
        }

        std::vector<peak> track_cut(track.begin()+start_idx, track.end());
        get_apex(track_cut, pk, true);
        pk.ts = gr::lora::pmod(pk.ts + d_num_samples/4, TIMESTAMP_MOD);

        float sum = 0;
        for (uint32_t i = d_overlaps*2; i < d_overlaps*(d_num_preamble-2); i++)
        {
          sum += track[i].h;
        }
        pk.h = sum / (d_overlaps*(d_num_preamble-4));
        return SYMBOL_PREAMBLE;
      }
      // TODO filter noise peak
      else if (len >= 2 && len <= 2*d_overlaps)
      {
        // get apex of this peak tracking
        get_apex(track, pk);
        return SYMBOL_DATA;
      }
      // TODO special case: more than two consecutive same data symbols

      return SYMBOL_BROKEN_DATA;
    }

    bool
    pyramid_demod_impl::add_symbol_to_packet(peak & pk, symbol_type st)
    {
      #if DEBUG >= DEBUG_VERBOSE
        std::cout << "symbol type: ";
        if (st == SYMBOL_PREAMBLE) std::cout << "PREAMBLE" << std::endl;
        else if (st == SYMBOL_DATA) std::cout << "DATA" << std::endl;
        else std::cout << "BROKEN" << std::endl;
      #endif

      if (st == SYMBOL_PREAMBLE)
      {
        // preamble detected, create a new packet
        uint16_t pkt_id = d_packet_id_pool.front();
        d_packet_id_pool.pop_front();
        d_packet[pkt_id].push_back(pk);
        d_packet_state_list.push_back(packet_state(pkt_id, d_ttl));
        #if DEBUG >= DEBUG_INFO
          std::cout << "New preamble detected (ts:" << std::fixed << std::setprecision(2) << pk.ts/(float)d_num_samples << ", bin:" << pk.bin << ", h:" << pk.h << ") Packet#" << pkt_id << std::endl;
        #endif
        return true;
      }
      else if (st == SYMBOL_DATA)
      {
        uint16_t pkt_idx  = 0;
        uint16_t pkt_id   = 0;
        float min_dis = std::numeric_limits<float>::infinity();
        bool found = false;

        for (uint32_t i = 0; i < d_packet_state_list.size(); i++)
        {
          // put peak into the best matched packet
          auto const & ps = d_packet_state_list[i];
          uint32_t ts_dis = gr::lora::pmod(pk.ts - d_packet[ps.packet_id][0].ts, TIMESTAMP_MOD);
          // candidate symbols must have valid timestamp
          if (ts_dis > 4 * d_num_samples && ts_dis < TIMESTAMP_MOD / 2)
          {
            float dis = gr::lora::pmod(ts_dis, d_num_samples) / (float) d_num_samples;
            // In definition above, "dis->0" and "dis->1" both mean their timestamp differences to the actual timestamp is small
            // therefore we transform dis to let "dis->1" represent larger difference
            dis = (dis > 0.5) ? (1 - dis) * 2 : dis * 2;
            // dis += std::abs(d_packet[ps.packet_id][0].h - pk.h) / d_packet[ps.packet_id][0].h;
            float h_dis = std::abs(d_packet[ps.packet_id][0].h - pk.h) / d_packet[ps.packet_id][0].h;
            // dis += h_dis;
            // peak height difference should not be large
            if (dis < min_dis && h_dis < 0.5)
            {
              found     = true;
              pkt_idx   = i;
              pkt_id    = ps.packet_id;
              min_dis   = dis;
            }
          }
        }

        if (found)
        {
          // new symbol added, reset TTL of the packet
          d_packet_state_list[pkt_idx].ttl = d_ttl;
          d_packet[pkt_id].push_back(pk);
          #if DEBUG >= DEBUG_INFO
            std::cout << "Add symbol (ts:" << std::fixed << std::setprecision(2) << pk.ts/(float)d_num_samples << ", bin:" << pk.bin << ", h:" << pk.h << ") to Packet#" << pkt_id << std::endl;
          #endif
          return true;
        }
        else
        {
          #if DEBUG >= DEBUG_INFO
            std::cout << "packet_state_list size: " << d_packet_state_list.size() << ", failed to classify symbol (ts:" << std::fixed << std::setprecision(2) << pk.ts/(float)d_num_samples << ", bin:" << pk.bin << ", h:" << pk.h << ")" << std::endl;
          #endif
          return false;
        }
      }
      else
      {
        #if DEBUG >= DEBUG_INFO
          std::cout << "Unrecognized symbol type!" << std::endl;
        #endif
        return false;
      }
    }

    void
    pyramid_demod_impl::check_and_update_track()
    {
      int erase_cnt = 0;
      for (auto const & bt : d_bin_track_id_list)
      {
        if (!bt.updated)
        {
          erase_cnt++ ;
          // this peak tracking is over, extract the apex
          peak pk(0,0,0,0);
          symbol_type st = get_central_peak(bt.track_id, pk);
          // #if DEBUG >= DEBUG_INFO
          //   if (st == SYMBOL_PREAMBLE || st == SYMBOL_DATA)
          //   {
          //     std::cout << "Add new symbol (ts:" << pk.ts << ", bin:" << pk.bin << ", h:" << pk.h << ")" << std::endl;
          //   }
          // #endif
          if (st == SYMBOL_PREAMBLE || st == SYMBOL_DATA)
          {
            bool res = add_symbol_to_packet(pk, st);
            #if DEBUG >= DEBUG_VERBOSE
              if (!res)
              {
                std::cout << "Failed to add symbol to packet." << std::endl;
              }
            #endif
          }
          d_track_id_pool.push_back(bt.track_id); // id recycle
          d_track[bt.track_id].clear();           // track vector recycle
        }
      }
      // #if DEBUG >= DEBUG_VERBOSE
      //   std::cout << "size before erase: " << d_bin_track_id_list.size() << ", erase_cnt: " << erase_cnt << std::endl;
      // #endif
      d_bin_track_id_list.erase(
        std::remove_if(
          d_bin_track_id_list.begin(),
          d_bin_track_id_list.end(),
          [](bin_track_id & bt) { return !bt.updated; }
        ),
        d_bin_track_id_list.end()
      );
      // #if DEBUG >= DEBUG_VERBOSE
      //   std::cout << "size after erase: " << d_bin_track_id_list.size() << std::endl;
      // #endif
      for (auto & bt : d_bin_track_id_list)
      {
        bt.updated = false;
      }
    }

    void
    pyramid_demod_impl::forecast (int noutput_items,
                          gr_vector_int &ninput_items_required)
    {
      ninput_items_required[0] = noutput_items * (1 << d_sf);
    }

    int
    pyramid_demod_impl::general_work (int noutput_items,
                       gr_vector_int &ninput_items,
                       gr_vector_const_void_star &input_items,
                       gr_vector_void_star &output_items)
    {
      if (ninput_items[0] < 4*d_num_samples) return 0;
      const gr_complex *in        = (const gr_complex *)  input_items[0];
      uint32_t  *out          = (uint32_t   *) output_items[0];
      uint32_t num_consumed   = d_num_samples / d_overlaps;
      float max_val               = 0;
      uint32_t max_index_sfd  = 0;
      float max_val_sfd           = 0;
      uint32_t tmp_idx        = 0;
      // #if DEBUG >= DEBUG_VERBOSE
      //   std::cout << "d_num_samples: " << d_num_samples <<  ", d_overlaps: " << d_overlaps << ", num_consumed: " << num_consumed << ", ts: " << d_ts_ref << std::endl;
      // #endif

      // Nomenclature:
      //  up_block   == de-chirping buffer to contain upchirp features: the preamble, sync word, and data chirps
      //  down_block == de-chirping buffer to contain downchirp features: the SFD
      gr_complex *buffer     = (gr_complex *)volk_malloc(d_fft_size*sizeof(gr_complex), volk_get_alignment());
      gr_complex *up_block   = (gr_complex *)volk_malloc(d_fft_size*sizeof(gr_complex), volk_get_alignment());
      gr_complex *up_block_w = (gr_complex *)volk_malloc(d_fft_size*sizeof(gr_complex), volk_get_alignment());
      gr_complex *down_block = (gr_complex *)volk_malloc(d_fft_size*sizeof(gr_complex), volk_get_alignment());
      float *fft_mag     = (float*)volk_malloc(d_fft_size*sizeof(float), volk_get_alignment());
      float *fft_mag_w   = (float*)volk_malloc(d_fft_size*sizeof(float), volk_get_alignment());
      float *fft_add     = (float*)volk_malloc(d_bin_size*sizeof(float), volk_get_alignment());
      float *fft_add_w   = (float*)volk_malloc(d_fft_size*sizeof(float), volk_get_alignment());

      if (buffer == NULL || up_block == NULL || down_block == NULL)
      {
        std::cerr << "Unable to allocate processing buffer!" << std::endl;
      }

      // Dechirp the incoming signal
      volk_32fc_x2_multiply_32fc(down_block, in, &d_upchirp[0], d_num_samples);
      volk_32fc_x2_multiply_32fc(up_block, in, &d_downchirp[0], d_num_samples);

      // Enable to write IQ to disk for debugging
      #if DUMP_IQ
        f_up_windowless.write((const char*)&up_block[0], d_num_samples*sizeof(gr_complex));
      #endif

      // Windowing
      volk_32fc_32f_multiply_32fc(up_block_w, up_block, &d_window[0], d_num_samples);

      #if DUMP_IQ
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

      // fft result magnitude summation
      volk_32fc_magnitude_32f(fft_mag, d_fft->get_outbuf(), d_fft_size);
      volk_32f_x2_add_32f(fft_add, fft_mag, &fft_mag[d_bin_size], d_bin_size);

      // apply FFT on windowed signal
      memset(d_fft->get_inbuf(),              0, d_fft_size*sizeof(gr_complex));
      memcpy(d_fft->get_inbuf(), &up_block_w[0], d_num_samples*sizeof(gr_complex));
      d_fft->execute();
      volk_32fc_magnitude_32f(fft_mag_w, d_fft->get_outbuf(), d_fft_size);
      volk_32f_x2_add_32f(fft_add_w, fft_mag_w, &fft_mag_w[d_bin_size], d_bin_size);

      // 1. peak tracking
      find_and_add_peak(fft_add, fft_add_w, fft_mag);
      // 2. stop tracking those peaks without updates, and push them into packet list
      check_and_update_track();
      // 3. check packets
      for (auto const & ps : d_packet_state_list)
      {
        #if DEBUG >= DEBUG_VERBOSE_VERBOSE
          std::cout << "packet id: " << ps.packet_id << ", ttl: " << ps.ttl << std::endl;
        #endif
        if (ps.ttl <= 0)
        {
          // send the demodulation result to decoder
          std::vector<uint16_t> symbols;

          auto & pkt = d_packet[ps.packet_id];
          uint32_t pre_ts  = pkt[0].ts;     // preamble timestamp
          uint32_t pre_bin = pkt[0].bin;    // preamble bin
          float        pre_h   = pkt[0].h;      // preamble peak height
          #if DEBUG >= DEBUG_VERBOSE_VERBOSE
            std::cout << "preamble ts: " << pre_ts << ", preamble bin: " << pre_bin << std::endl;
            std::cout << "ts: ";
            for (uint32_t i = 1; i < pkt.size(); i++)
            {
              std::cout << pkt[i].ts << ", ";
            }
            std::cout << std::endl << "bin: ";
            for (uint32_t i = 1; i < pkt.size(); i++)
            {
              std::cout << pkt[i].bin << ", ";
            }
            std::cout << std::endl;

            std::cout << " d_packet size: ";
            for (uint32_t i = 0; i < d_packet.size(); i++)
            {
              std::cout << d_packet[i].size() << ",";
            }
            std::cout << std::endl;

            std::cout << "current packet id: " << ps.packet_id << ", d_packet id: ";
            for (uint32_t i = 0; i < d_packet_id_pool.size(); i++)
            {
              std::cout << d_packet_id_pool[i] << ",";
            }
            std::cout << std::endl;
          #endif

          #if DEBUG >= DEBUG_INFO
            std::cout << "Finished Packet#" << ps.packet_id << ": ";
          #endif

          // set relative timestamp to preamble
          // then ts_preamble = 0
          for (auto & pk : pkt)
          {
            pk.ts = gr::lora::pmod(pk.ts-pre_ts, TIMESTAMP_MOD);
          }
          pre_ts = 0;

          // sort peaks by their timestamps
          sort(pkt.begin(), pkt.end(),
            [](const peak & a, const peak & b) -> bool
          {
            return a.ts < b.ts;
          });

          #if DEBUG >= DEBUG_VERBOSE
            std::cout << std::endl;
            for (auto & pk : pkt)
            {
              std::cout << "(ts: " << pk.ts/(float)d_num_samples << ", bin: " << pk.bin << ", h: " << pk.h << ")" << std::endl;
            }
          #endif

          // LoRa PHY: Preamble + NetID(2) + SFD(2.25) + Data Payload
          // There are 4.25 symbols between preamble and data payload
          // ts_data - ts_preamble = 5*d_num_samples (ts_preamble has 0.25 fix in our implementation)
          // The first data symbol timestamp is in ts_preamble+[4.5,5.5]*d_num_samples
          uint32_t ts_interval_l = 4*d_num_samples + d_num_samples/2;
          // "i" starts from 1, ignoring preamble
          for (uint32_t start_idx = 1; start_idx < pkt.size(); )
          {
            bool is_first = true;
            bool found = false;
            uint32_t end_idx = start_idx;
            for (; end_idx < pkt.size(); end_idx++)
            {
              if (is_first)
              {
                if ( pkt[end_idx].ts > ts_interval_l && pkt[end_idx].ts < ts_interval_l + d_num_samples )
                {
                  start_idx = end_idx;
                  is_first = false;
                  found = true;
                  #if DEBUG >= DEBUG_VERBOSE
                    std::cout << std::endl << "pkt[end_idx].ts: " << pkt[end_idx].ts/(float)d_num_samples << ", ts_interval: " << ts_interval_l/(float)d_num_samples << ", r: " << (ts_interval_l + d_num_samples)/(float)d_num_samples << std::endl;
                  #endif
                }
              }
              else
              {
                if (!( pkt[end_idx].ts > ts_interval_l && pkt[end_idx].ts < ts_interval_l + d_num_samples ))
                {
                  break;
                }
              }
            }

            if (found)
            {
              // search the best matched peak for the packet
              float min_dis = std::numeric_limits<float>::infinity();
              uint32_t idx = start_idx;
              for (uint32_t i = start_idx; i < end_idx; i++)
              {
                float dis = get_dis(pkt[i].ts, pkt[i].h, pre_ts, pre_h);
                if (dis < min_dis)
                {
                  min_dis = dis;
                  idx = i;
                }
              }

              // ts could have overflowed
              int bin_shift = gr::lora::pmod(pkt[idx].ts - pre_ts, d_num_samples) * d_bin_size / d_num_samples;
              uint32_t bin = gr::lora::pmod(pkt[idx].bin - pre_bin - bin_shift, d_bin_size);
              symbols.push_back(bin / d_fft_size_factor);

              #if DEBUG >= DEBUG_VERBOSE
                std::cout << "bin: " << bin / d_fft_size_factor << ", packet bin: " << pkt[idx].bin << ", bin_shift: " << bin_shift << std::endl;
              #elif DEBUG >= DEBUG_INFO
                std::cout << bin / d_fft_size_factor << ",";
              #endif
            }
            else
            {
              symbols.push_back(0);
              #if DEBUG >= DEBUG_INFO
                std::cout << "missing,";
              #endif
            }

            start_idx = end_idx;
            ts_interval_l = gr::lora::pmod(ts_interval_l + d_num_samples, TIMESTAMP_MOD);
          }
          #if DEBUG >= DEBUG_INFO
            std::cout << std::endl;
          #endif

          // LoRa data payload has at least 8 symbols
          if (symbols.size() >= 8)
          {
            pmt::pmt_t dict = pmt::make_dict();
            dict = pmt::dict_add(dict, pmt::intern("id"), pmt::intern("packet"));
            pmt::pmt_t output = pmt::init_u16vector(symbols.size(), symbols);
            pmt::pmt_t msg_pair = pmt::cons(dict, output);
            message_port_pub(d_out_port, msg_pair);
          }

          pkt.clear();
          d_packet_id_pool.push_back(ps.packet_id);
        }
      }

      // packets with negative ttl have been sent to decoder by the above for loop
      d_packet_state_list.erase(
        std::remove_if(
          d_packet_state_list.begin(),
          d_packet_state_list.end(),
          [](packet_state const & ps) { return ps.ttl <= 0; }
        ),
        d_packet_state_list.end()
      );
      for (auto & ps : d_packet_state_list)
      {
        ps.ttl -= 1;
      }

      d_ts_ref   = gr::lora::pmod(d_ts_ref + d_num_samples / d_overlaps, TIMESTAMP_MOD);
      d_bin_ref  = gr::lora::pmod(d_bin_ref + d_bin_size / d_overlaps, d_bin_size);

      #if DUMP_IQ
        f_raw.write((const char*)&in[0], num_consumed*sizeof(gr_complex));
      #endif

      consume_each(num_consumed);

      free(down_block);
      free(up_block);
      free(up_block_w);
      free(buffer);
      free(fft_mag);
      free(fft_mag_w);
      free(fft_add);
      free(fft_add_w);

      return noutput_items;
    }

  } /* namespace lora */
} /* namespace gr */
