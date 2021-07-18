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

#ifndef INCLUDED_LORA_PYRAMID_DEMOD_IMPL_H
#define INCLUDED_LORA_PYRAMID_DEMOD_IMPL_H

#include <cmath>
#include <cstdlib>
#include <vector>
#include <deque>
#include <queue>
#include <complex>
#include <fstream>
#include <limits>
#include <gnuradio/fft/fft.h>
#include <gnuradio/fft/window.h>
#include <volk/volk.h>
#include <lora/pyramid_demod.h>
#include "utilities.h"

namespace gr {
  namespace lora {
    enum symbol_type {
      SYMBOL_PREAMBLE,
      SYMBOL_DATA,
      SYMBOL_BROKEN_DATA,
    };

    struct peak {
      uint32_t     ts;     // timestamp
      uint32_t     bin;    // peak position bin
      float            h;      // peak height
      float            h_single;     // single peak height
      peak(uint32_t ts, uint32_t bin, float h, float h_single)
        : ts(ts), bin(bin), h(h), h_single(h_single)
      {}
    };

    struct bin_track_id {
      uint32_t     bin;           // peak bin
      uint16_t   track_id;      // index to reference peak_track
      bool             updated;       // whether a peak tracking process is over
      bin_track_id(uint32_t bin, uint16_t track_id, bool updated)
        : bin(bin), track_id(track_id), updated(updated)
      {}
      // bin_track_id(bin_track_id && other)
      //   : bin(other.bin), track_id(other.track_id), updated(other.updated)
      // {}
      // bin_track_id& operator=(const bin_track_id & other) = default;
    };

    struct packet_state {
      uint16_t   packet_id;     // index to reference packet list
      short   ttl;           // time to live
      packet_state(uint16_t packet_id, uint16_t ttl)
        : packet_id(packet_id), ttl(ttl)
      {}
    };

    class pyramid_demod_impl : public pyramid_demod
    {
     private:
      pmt::pmt_t d_header_port;
      pmt::pmt_t d_out_port;

      uint16_t  d_sf;
      bool d_ldr;

      uint16_t  d_num_symbols;
      uint16_t  d_fft_size_factor;
      uint32_t    d_fft_size;
      uint16_t  d_overlaps;
      short           d_ttl;
      uint16_t  d_bin_tolerance;
      uint16_t  d_offset;
      uint16_t  d_p;
      uint32_t    d_num_samples;
      uint32_t    d_bin_size;
      uint16_t  d_num_preamble;

      float           d_cfo;
      float           d_power;
      float           d_threshold;
      bool            d_squelched;

      uint32_t    d_preamble_idx;
      uint16_t  d_sfd_idx;
      std::vector<uint32_t>  d_argmax_history;
      std::vector<uint16_t>  d_sfd_history;
      uint16_t  d_sync_recovery_counter;

      uint32_t    d_bin_ref;
      uint32_t    d_ts_ref;
      std::vector<std::vector<peak>>        d_track;
      std::vector<bin_track_id>             d_bin_track_id_list;
      std::deque<uint16_t>            d_track_id_pool;
      std::vector<std::vector<peak>>        d_packet;
      std::vector<packet_state>             d_packet_state_list;
      std::deque<uint16_t>            d_packet_id_pool;

      fft::fft_complex   *d_fft;
      std::vector<float> d_window;
      float              d_beta;

      std::vector<gr_complex> d_upchirp;
      std::vector<gr_complex> d_downchirp;

      std::vector<uint16_t> d_symbols;

      std::ofstream f_raw, f_up_windowless, f_up, f_down, f_fft;

     public:
      pyramid_demod_impl(uint16_t spreading_factor,
                         bool low_data_rate,
                         float beta,
                         uint16_t fft_factor,
                         float threshold,
                         float fs_bw_ratio);
      ~pyramid_demod_impl();

      uint16_t argmax(gr_complex *fft_result, bool update_squelch);
      uint32_t argmax_32f(float *fft_result, bool update_squelch, float *max_val_p);
      float get_dis(uint32_t ts1, float h1, uint32_t ts2, float h2);

      // void parse_header(pmt::pmt_t dict);

      void find_and_add_peak(float *fft_mag_add, float *fft_mag_add_w, float *fft_mag);
      void get_apex(std::vector<peak> & track, peak & pk, bool is_preamble=false);
      symbol_type get_central_peak(uint16_t track_id, peak & pk);
      bool add_symbol_to_packet(peak & pk, symbol_type st);
      void check_and_update_track();

      // Where all the action really happens
      void forecast(int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace lora
} // namespace gr

#endif /* INCLUDED_LORA_PYRAMID_DEMOD_IMPL_H */

