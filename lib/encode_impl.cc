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
#include "encode_impl.h"

#define HAMMING_P1_BITMASK 0x0D  // 0b00001101
#define HAMMING_P2_BITMASK 0x0B  // 0b00001011
#define HAMMING_P3_BITMASK 0x07  // 0b00000111
#define HAMMING_P4_BITMASK 0x0F  // 0b00001111
#define HAMMING_P5_BITMASK 0x0E  // 0b00001110

#define DEBUG_OUTPUT 0  // Controls debug print statements

namespace gr
{
  namespace lora
  {

    encode::sptr
    encode::make(short spreading_factor,
                 short code_rate,
                 bool  crc,
                 bool  low_data_rate,
                 bool  header)
    {
      return gnuradio::get_initial_sptr
        (new encode_impl(spreading_factor, code_rate, crc, low_data_rate, header));
    }

    /*
     * The private constructor
     */
    encode_impl::encode_impl(short spreading_factor,
                             short code_rate,
                             bool  crc,
                             bool  low_data_rate,
                             bool  header)
        : gr::block("encode",
                    gr::io_signature::make(0, 0, 0),
                    gr::io_signature::make(0, 0, 0)),
          d_sf(spreading_factor),
          d_cr(code_rate),
          d_crc(crc),
          d_ldr(low_data_rate),
          d_header(header)
    {
      assert((d_sf > 5) && (d_sf < 13));
      assert((d_cr > 0) && (d_cr < 5));
      if (d_sf == 6) assert(!header);

      d_in_port = pmt::mp("in");
      d_out_port = pmt::mp("out");

      message_port_register_in(d_in_port);
      message_port_register_out(d_out_port);

      set_msg_handler(d_in_port, boost::bind(&encode_impl::encode, this, _1));

      d_whitening_sequence = whitening_sequence;

      d_interleaver_size = d_sf;

      d_fft_size = (1 << spreading_factor);
    }

    /*
     * Our virtual destructor.
     */
    encode_impl::~encode_impl()
    {
    }

    void
    encode_impl::gen_header(std::vector<unsigned char> &nibbles, uint8_t payload_len)
    {
      uint8_t cr_crc = (d_cr << 1) | d_crc;
      uint8_t cks = gr::lora::header_checksum(payload_len, cr_crc);
      nibbles.push_back(payload_len >> 4);
      nibbles.push_back(payload_len & 0xF);
      nibbles.push_back(cr_crc);
      nibbles.push_back(cks >> 4);
      nibbles.push_back(cks & 0xF);
    }

    uint16_t
    encode_impl::calc_sym_num(uint8_t payload_len)
    {
      double tmp = 2 * payload_len - d_sf + 7 + 4 * d_crc - 5 * (1 - d_header);
      return 8 + std::max((4 + d_cr) * (uint16_t)ceil(tmp / (d_sf - 2 * d_ldr)), 0);
    }

    void
    encode_impl::to_gray(std::vector<uint16_t> &symbols)
    {
      for (int i = 0; i < symbols.size(); i++)
      {
        symbols[i] = (symbols[i] >> 1) ^ symbols[i];
      }
    }

    void
    encode_impl::from_gray(std::vector<uint16_t> &symbols)
    {
      for (int i = 0; i < symbols.size(); i++)
      {
        symbols[i] = symbols[i] ^ (symbols[i] >> 16);
        symbols[i] = symbols[i] ^ (symbols[i] >>  8);
        symbols[i] = symbols[i] ^ (symbols[i] >>  4);
        symbols[i] = symbols[i] ^ (symbols[i] >>  2);
        symbols[i] = symbols[i] ^ (symbols[i] >>  1);
        symbols[i] = (i < 8 || d_ldr) ? (symbols[i] * 4 + 1) % (1 << d_sf) : (symbols[i] + 1) % (1 << d_sf);
      }
    }

    void
    encode_impl::whiten(std::vector<unsigned char> &bytes, uint8_t len)
    {
      for (int i = 0; i < len && i < whitening_sequence_length; i++)
      {
        bytes[i] = ((unsigned char)(bytes[i] & 0xFF) ^ d_whitening_sequence[i]) & 0xFF;
      }
    }

    void
    encode_impl::print_bitwise_u8(std::vector<unsigned char> &buffer)
    {
      for (int i = 0; i < buffer.size(); i++)
      {
        std::cout << i << "\t" << std::bitset<8>(buffer[i] & 0xFF) << "\t";
        std::cout << std::hex << (buffer[i] & 0xFF) << std::endl;
      }
    }

    void
    encode_impl::print_bitwise_u16(std::vector<uint16_t> &buffer)
    {
      for (int i = 0; i < buffer.size(); i++)
      {
        std::cout << i << "\t" << std::bitset<16>(buffer[i] & 0xFFFF) << "\t";
        std::cout << std::hex << (buffer[i] & 0xFFFF) << std::endl;
      }
    }

    // Forward interleaver dimensions:
    //  PPM   == number of bits per symbol OUT of interleaver        AND number of codewords IN to interleaver
    //  RDD+4 == number of bits per codeword IN to interleaver       AND number of interleaved codewords OUT of interleaver
    //
    // bit width in:  (4+rdd)   block length: ppm
    // bit width out: ppm       block length: (4+rdd)
    void
    encode_impl::interleave(std::vector<unsigned char> &codewords,
                            std::vector<uint16_t> &symbols)
    {
      uint32_t bits_per_word = 8;
      uint8_t ppm = d_sf - 2;
      for (uint32_t start_idx = 0; start_idx + ppm - 1< codewords.size();) {
        bits_per_word = (start_idx == 0) ? 8 : (d_cr + 4);
        ppm = (start_idx == 0) ? (d_sf - 2) : (d_sf - 2 * d_ldr);
        std::vector<uint16_t> block(bits_per_word, 0u);

        for (uint32_t i = 0; i < ppm; i++) {
          const uint32_t word = codewords[start_idx + i];
          for (uint32_t j = (1u << (bits_per_word - 1)), x = bits_per_word - 1; j; j >>= 1u, x--) {
            block[x] |= !!(word & j) << i;
          }
        }

        for (uint32_t i = 0; i < bits_per_word; i++)
        {
          // rotate each element to the right by i bits
          block[i] = gr::lora::rotl(block[i], 2*ppm - i, ppm);
        }

        symbols.insert(symbols.end(), block.begin(), block.end());

        start_idx = start_idx + ppm;
      }
    }

    void
    encode_impl::hamming_encode(std::vector<unsigned char> &nibbles,
                                std::vector<unsigned char> &codewords)
    {
      unsigned char p1, p2, p3, p4, p5;
      unsigned char mask;

      for (int i = 0; i < nibbles.size(); i++)
      {
        p1 = parity((unsigned char)nibbles[i], mask = (unsigned char)HAMMING_P1_BITMASK);
        p2 = parity((unsigned char)nibbles[i], mask = (unsigned char)HAMMING_P2_BITMASK);
        p3 = parity((unsigned char)nibbles[i], mask = (unsigned char)HAMMING_P3_BITMASK);
        p4 = parity((unsigned char)nibbles[i], mask = (unsigned char)HAMMING_P4_BITMASK);
        p5 = parity((unsigned char)nibbles[i], mask = (unsigned char)HAMMING_P5_BITMASK);

        uint8_t cr_now = (i < d_sf - 2) ? 4 : d_cr;

        switch (cr_now)
        {
        case 1:
          codewords.push_back((p4 << 4) |
                              (nibbles[i] & 0xF));
          break;
        case 2:
          codewords.push_back((p5 << 5) |
                              (p3 << 4) |
                              (nibbles[i] & 0xF));
          break;
        case 3:
          codewords.push_back((p2 << 6) |
                              (p5 << 5) |
                              (p3 << 4) |
                              (nibbles[i] & 0xF));
          break;
        case 4:
          codewords.push_back((p1 << 7) |
                              (p2 << 6) |
                              (p5 << 5) |
                              (p3 << 4) |
                              (nibbles[i] & 0xF));
          break;
        default:
          // THIS CASE SHOULD NOT HAPPEN
          std::cerr << "Invalid Code Rate  -- this state should never occur." << std::endl;
          break;
        }
      }
    }

    unsigned char
    encode_impl::parity(unsigned char c, unsigned char bitmask)
    {
      unsigned char parity = 0;
      unsigned char shiftme = c & bitmask;

      for (int i = 0; i < 8; i++)
      {
        if (shiftme & 0x1) parity++;
        shiftme = shiftme >> 1;
      }

      return parity % 2;
    }

    void
    encode_impl::print_payload(std::vector<unsigned char> &payload)
    {
      std::cout << "Encoded LoRa packet (hex): ";
      for (int i = 0; i < payload.size(); i++)
      {
        std::cout << std::hex << (uint32_t)payload[i] << " ";
      }
      std::cout << std::endl;
    }

    void
    encode_impl::encode(pmt::pmt_t msg)
    {
      pmt::pmt_t bytes(pmt::cdr(msg));

      size_t pkt_len(0);
      const uint8_t *bytes_in_p = pmt::u8vector_elements(bytes, pkt_len);

      std::vector<uint8_t> bytes_in(bytes_in_p, bytes_in_p + pkt_len);
      std::vector<uint8_t> nibbles;
      std::vector<uint8_t> codewords;
      std::vector<uint8_t> payload_nibbles;
      std::vector<uint16_t> symbols;

      if (d_crc)
      {
        uint16_t checksum = gr::lora::data_checksum(&bytes_in[0], pkt_len);
        bytes_in.push_back(checksum & 0xFF);
        bytes_in.push_back((checksum >> 8) & 0xFF);
      }

      uint16_t sym_num = calc_sym_num(pkt_len);
      uint16_t nibble_num = d_sf - 2 + (sym_num - 8) / (d_cr + 4) * (d_sf - 2 * d_ldr);
      uint16_t redundant_num = ceil((nibble_num - 2 * bytes_in.size()) / 2);
      for (int i = 0; i < redundant_num; i++)
      {
        bytes_in.push_back(0);
      }

      whiten(bytes_in, pkt_len);

      // split bytes into separate data nibbles
      for (int i = 0; i < nibble_num; i++)
      {
        if (i % 2 == 0)
        {
          payload_nibbles.push_back(bytes_in[i / 2] & 0xF);
        }
        else
        {
          payload_nibbles.push_back(bytes_in[i / 2] >> 4);
        }
      }

      if (d_header)
      {
        gen_header(nibbles, pkt_len);
      }

#if DEBUG_OUTPUT
      std::cout << "Header nibbles:" << std::endl;
      print_bitwise_u8(nibbles);
      std::cout << "Payload nibbles:" << std::endl;
      print_bitwise_u8(payload_nibbles);
#endif

      nibbles.insert(nibbles.end(), payload_nibbles.begin(), payload_nibbles.end());
      hamming_encode(nibbles, codewords);

#if DEBUG_OUTPUT
      std::cout << "Codewords:" << std::endl;
      print_bitwise_u8(codewords);
#endif

      interleave(codewords, symbols);

#if DEBUG_OUTPUT
      std::cout << "Interleaved Symbols: " << std::endl;
      print_bitwise_u16(symbols);
#endif

      from_gray(symbols);

#if DEBUG_OUTPUT
      std::cout << "Modulated Symbols: " << std::endl;
      print_bitwise_u16(symbols);
#endif

      pmt::pmt_t output = pmt::init_u16vector(symbols.size(), symbols);
      pmt::pmt_t msg_pair = pmt::cons(pmt::make_dict(), output);

      message_port_pub(d_out_port, msg_pair);
    }

  } /* namespace lora */
} /* namespace gr */
