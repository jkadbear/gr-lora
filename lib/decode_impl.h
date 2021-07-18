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

#ifndef INCLUDED_LORA_DECODE_IMPL_H
#define INCLUDED_LORA_DECODE_IMPL_H

#include <iostream>
#include <bitset>
#include <lora/decode.h>
#include "utilities.h"

namespace gr {
  namespace lora {

    class decode_impl : public decode
    {
     private:
      pmt::pmt_t d_in_port;
      pmt::pmt_t d_out_port;
      pmt::pmt_t d_header_port;

      const unsigned char *d_whitening_sequence;

      unsigned char d_sf;
      unsigned char d_cr;
      unsigned char d_payload_len;
      bool          d_crc;
      bool          d_ldr;
      bool          d_header;

      uint16_t d_fft_size;
      unsigned char  d_interleaver_size;

      std::vector<uint16_t> d_symbols;
      std::vector<unsigned char> d_codewords;
      std::vector<unsigned char> d_bytes;

     public:
      decode_impl( short spreading_factor,
                   bool  header,
                   short payload_len,
                   short code_rate,
                   bool  crc,
                   bool  low_data_rate);
      ~decode_impl();

      void to_gray(std::vector<uint16_t> &symbols);
      void from_gray(std::vector<uint16_t> &symbols);
      void whiten(std::vector<unsigned char> &codewords);
      void deinterleave(std::vector<uint16_t> &symbols, std::vector<unsigned char> &codewords, unsigned char ppm, unsigned char rdd);
      void hamming_decode(std::vector<unsigned char> &codewords, std::vector<unsigned char> &bytes, unsigned char rdd);
      unsigned char parity(unsigned char c, unsigned char bitmask);
      void print_payload(std::vector<unsigned char> &payload);

      void print_bitwise_u8 (std::vector<unsigned char>  &buffer);
      void print_bitwise_u16(std::vector<uint16_t> &buffer);

      void decode(pmt::pmt_t msg);

    };

  } // namespace lora
} // namespace gr

#endif /* INCLUDED_LORA_DECODE_IMPL_H */

