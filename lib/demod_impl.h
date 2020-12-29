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

#ifndef INCLUDED_LORA_DEMOD_IMPL_H
#define INCLUDED_LORA_DEMOD_IMPL_H

#include <cmath>
#include <cstdlib>
#include <vector>
#include <queue>
#include <complex>
#include <fstream>
#include <gnuradio/fft/fft.h>
#include <gnuradio/fft/window.h>
#include <volk/volk.h>
#include "lora/demod.h"
#include "utilities.h"

namespace gr {
  namespace lora {

    class demod_impl : public demod
    {
     private:
      pmt::pmt_t d_header_port;
      pmt::pmt_t d_out_port;

      demod_state_t   d_state;
      uint8_t d_sf;
      uint8_t d_cr;
      uint8_t d_payload_len;
      bool d_crc;
      bool d_ldr;
      bool d_header;
      bool d_header_received;
      bool d_header_valid;

      unsigned short  d_num_symbols;
      unsigned short  d_fft_size_factor;
      unsigned int    d_fft_size;
      unsigned short  d_overlaps;
      unsigned short  d_offset;
      unsigned short  d_p;
      unsigned int    d_num_samples;
      unsigned int    d_bin_len;
      unsigned int    d_preamble_drift_max;

      uint32_t d_packet_symbol_len;

      float           d_cfo;

      unsigned int    d_preamble_idx;
      unsigned short  d_sfd_idx;
      std::vector<unsigned int>  d_argmax_history;
      std::vector<unsigned short>  d_sfd_history;
      unsigned short  d_sync_recovery_counter;

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
      demod_impl( uint8_t   spreading_factor,
                  bool      header,
                  uint8_t   payload_len,
                  uint8_t   cr,
                  bool      crc,
                  bool      low_data_rate,
                  float     beta,
                  uint16_t  fft_factor,
                  uint8_t   peak_search_algorithm,
                  uint16_t  peak_search_phase_k,
                  float     fs_bw_ratio);
      ~demod_impl();

      unsigned short argmax(gr_complex *fft_result);
      unsigned int argmax_32f(float *fft_result, float *max_val_p);
      unsigned int search_fft_peak(const lv_32fc_t *fft_result,
                                   float *buffer1, float *buffer2,
                                   gr_complex *buffer_c, float *max_val_p);
      unsigned int fft_add(const lv_32fc_t *fft_result, float *buffer, gr_complex *buffer_c,
                           float *max_val_p, float phase_offset);
      void dynamic_compensation(std::vector<uint16_t>& compensated_symbols);
      
      void parse_header(pmt::pmt_t dict);

      // Where all the action really happens
      void forecast (int noutput_items, gr_vector_int &ninput_items_required);

      int general_work(int noutput_items,
           gr_vector_int &ninput_items,
           gr_vector_const_void_star &input_items,
           gr_vector_void_star &output_items);
    };

  } // namespace lora
} // namespace gr

#endif /* INCLUDED_LORA_DEMOD_IMPL_H */

