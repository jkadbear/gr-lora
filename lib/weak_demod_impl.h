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

#ifndef INCLUDED_LORA_WEAK_DEMOD_IMPL_H
#define INCLUDED_LORA_WEAK_DEMOD_IMPL_H


#include <cmath>
#include <cstdlib>
#include <vector>
#include <queue>
#include <complex>
#include <fstream>
#include <gnuradio/fft/fft.h>
#include <gnuradio/fft/window.h>
#include <volk/volk.h>
#include <lora/weak_demod.h>
#include "utilities.h"

namespace gr {
  namespace lora {

    class weak_demod_impl : public weak_demod
    {
     private:
      pmt::pmt_t d_header_port;
      pmt::pmt_t d_out_port;

      weak_demod_state_t d_state;
      uint8_t d_sf;
      uint8_t d_cr;
      uint8_t d_payload_len;
      bool d_crc;
      bool d_ldr;
      bool d_header;
      bool d_header_received;
      bool d_header_valid;

      uint16_t d_num_symbols;
      uint16_t d_fft_size_factor;
      uint32_t d_fft_size;
      uint16_t d_overlaps;
      uint16_t d_offset;
      uint16_t d_p;
      uint32_t d_num_samples;
      uint32_t d_bin_size;
      uint32_t d_preamble_drift_max;
      uint32_t d_sym_num;

      uint32_t d_packet_symbol_len;

      float d_cfo;

      uint32_t d_preamble_idx;
      uint16_t d_sfd_idx;
      std::vector<uint32_t> d_argmax_history;
      std::vector<uint16_t> d_sfd_history;
      uint16_t d_sync_recovery_counter;


      uint16_t d_peak_search_algorithm;
      uint16_t d_peak_search_phase_k;

      fft::fft_complex   *d_fft;
      std::vector<float> d_window;
      float              d_beta;

      std::vector<gr_complex> d_upchirp;
      std::vector<gr_complex> d_downchirp;

      std::vector<float> d_symbols;

      std::ofstream f_raw, f_up_windowless, f_up, f_down, f_fft;

     public:
      weak_demod_impl(uint8_t   spreading_factor,
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
                  float     fs_bw_ratio);
      ~weak_demod_impl();

      void parse_header(pmt::pmt_t dict);

      void dechirp(bool is_up,
                        const gr_complex *in,
                        gr_complex *up_block,
                        float *fft_mag,
                        float *fft_add);

      uint32_t search_fft_peak(bool is_up,
                      const gr_complex *in,
                      gr_complex *block1,
                      gr_complex *block2,
                      float *fft_mag1,
                      float *fft_mag2,
                      float *fft_add1,
                      float *fft_add2,
                      float *p_max_val);
      
      void dynamic_compensation(std::vector<uint16_t>& compensated_symbols);

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);

    };

  } // namespace lora
} // namespace gr

#endif /* INCLUDED_LORA_WEAK_DEMOD_IMPL_H */

