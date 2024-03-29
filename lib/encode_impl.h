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

#ifndef INCLUDED_LORA_ENCODE_IMPL_H
#define INCLUDED_LORA_ENCODE_IMPL_H

#include <volk/volk.h>
#include <bitset>
#include <cmath>
#include <lora/encode.h>
#include "utilities.h"

namespace gr {
  namespace lora {

    class encode_impl : public encode
    {
     private:
      const uint8_t *d_whitening_sequence;

      pmt::pmt_t d_in_port;
      pmt::pmt_t d_out_port;
      
      unsigned char d_sf;
      unsigned char d_cr;
      bool          d_crc;
      bool          d_ldr;
      bool          d_header;

      uint16_t d_fft_size;
      unsigned char  d_interleaver_size;

     public:
      encode_impl(  short spreading_factor,
                    short code_rate,
                    bool  crc,
                    bool  low_data_rate,
                    bool  header);
      ~encode_impl();

      void gen_header(std::vector<unsigned char> &nibbheader_nibblesles, uint8_t payload_len);
      uint16_t calc_sym_num(uint8_t payload_len);
      void to_gray(std::vector<uint16_t> &symbols);
      void from_gray(std::vector<uint16_t> &symbols);
      void whiten(std::vector<unsigned char> &bytes, uint8_t len);
      void interleave(std::vector<unsigned char> &codewords, std::vector<uint16_t> &symbols);
      void hamming_encode(std::vector<unsigned char> &nibbles, std::vector<unsigned char> &codewords);
      unsigned char parity(unsigned char c, unsigned char bitmask);
      void print_payload(std::vector<unsigned char> &payload);

      void print_bitwise_u8 (std::vector<unsigned char>  &buffer);
      void print_bitwise_u16(std::vector<uint16_t> &buffer);

      void encode(pmt::pmt_t msg);

    };

  } // namespace lora
} // namespace gr

#endif /* INCLUDED_LORA_ENCODE_IMPL_H */

