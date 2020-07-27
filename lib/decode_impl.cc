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
#include "decode_impl.h"

#define MAXIMUM_RDD 4

// #define HAMMING_P1_BITMASK 0xAA  // 0b10101010
// #define HAMMING_P2_BITMASK 0x66  // 0b01100110
// #define HAMMING_P4_BITMASK 0x1E  // 0b00011110
// #define HAMMING_P8_BITMASK 0xFE  // 0b11111110

// not standard hamming code
// LoRa codeword: p4 p2 p1 p3 d1 d2 d4 d3
// pi: the ith parity; di: the ith data bits
#define HAMMING_P1_BITMASK 0x2E  // 0b00101110
#define HAMMING_P2_BITMASK 0x4B  // 0b01001011
#define HAMMING_P3_BITMASK 0x17  // 0b00010111
#define HAMMING_P4_BITMASK 0xFF  // 0b11111111
#define HAMMING_D1_BITMASK 0x08  // 0b00001000
#define HAMMING_D2_BITMASK 0x04  // 0b00000100
#define HAMMING_D3_BITMASK 0x01  // 0b00000001
#define HAMMING_D4_BITMASK 0x02  // 0b00000010

#define INTERLEAVER_BLOCK_SIZE 12

#define DEBUG_OUTPUT 0

namespace gr {
  namespace lora {

    decode::sptr
    decode::make(   short spreading_factor,
                    short code_rate,
                    bool  low_data_rate,
                    bool  header)
    {
      return gnuradio::get_initial_sptr
        (new decode_impl(spreading_factor, code_rate, low_data_rate, header));
    }

    /*
     * The private constructor
     */
    decode_impl::decode_impl( short spreading_factor,
                              short code_rate,
                              bool  low_data_rate,
                              bool  header)
      : gr::block("decode",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0)),
        d_sf(spreading_factor),
        d_cr(code_rate),
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

      set_msg_handler(d_in_port, boost::bind(&decode_impl::decode, this, _1));

      d_whitening_sequence = whitening_sequence;
      if (d_sf < 6 | d_sf > 12)
      {
        std::cerr << "Invalid spreading factor -- this state should never occur." << std::endl;
      }

      if (d_header)
      {
        std::cout << "Warning: Explicit header mode is not yet supported." << std::endl;
        std::cout << "         Using an implicit whitening sequence: demodulation will work correctly; decoding will not." << std::endl;
      }

      d_interleaver_size = d_sf;

      d_fft_size = (1 << spreading_factor);
    }

    /*
     * Our virtual destructor.
     */
    decode_impl::~decode_impl()
    {
    }

    void
    decode_impl::to_gray(std::vector<unsigned short> &symbols)
    {
      for (int i = 0; i < symbols.size(); i++)
      {
        symbols[i] = (symbols[i] >> 1) ^ symbols[i];
      }
    }

    void
    decode_impl::from_gray(std::vector<unsigned short> &symbols)
    {
      for (int i = 0; i < symbols.size(); i++)
      {
        symbols[i] = symbols[i] ^ (symbols[i] >> 16);
        symbols[i] = symbols[i] ^ (symbols[i] >>  8);
        symbols[i] = symbols[i] ^ (symbols[i] >>  4);
        symbols[i] = symbols[i] ^ (symbols[i] >>  2);
        symbols[i] = symbols[i] ^ (symbols[i] >>  1);
      }
    }

    void
    decode_impl::whiten(std::vector<unsigned char> &codewords)
    {
      for (int i = 0; (i < codewords.size()) && (i < whitening_sequence_length); i++)
      {
        codewords[i] ^= d_whitening_sequence[i];
      }
    }

    // Forward interleaver dimensions:
    //  PPM   == number of bits per symbol OUT of interleaver        AND number of codewords IN to interleaver
    //  RDD+4 == number of bits per codeword IN to interleaver       AND number of interleaved codewords OUT of interleaver
    //
    // bit width in:  (4+rdd)   block length: ppm
    // bit width out: ppm       block length: (4+rdd)

    // Reverse interleaver (de-interleaver) dimensions:
    //  PPM   == number of bits per symbol IN to deinterleaver       AND number of codewords OUT of deinterleaver
    //  RDD+4 == number of bits per codeword OUT of deinterleaver    AND number of interleaved codewords IN to deinterleaver
    //
    // bit width in:  ppm       block length: (4+rdd)
    // bit width out: (4+rdd)   block length: ppm
    void
    decode_impl::deinterleave(std::vector <unsigned short> &symbols,
                              std::vector <unsigned char> &codewords,
                              unsigned char ppm,
                              unsigned char rdd)
    {
      const uint32_t bits_per_word = rdd + 4;
      const uint32_t offset_start  = ppm - 1u;

      for (uint32_t start_idx = 0; start_idx < symbols.size(); start_idx += bits_per_word) {
        std::vector<uint8_t> block(ppm, 0u);
        for (uint32_t i = 0; i < bits_per_word; i++) {
          const uint32_t word = gr::lora::rotl(symbols[start_idx + i], i, ppm);

          for (uint32_t j = (1u << offset_start), x = offset_start; j; j >>= 1u, x--) {
            block[x] |= !!(word & j) << i;
          }
        }

        codewords.insert(codewords.end(), block.begin(), block.end());
      }
    }

    void
    decode_impl::hamming_decode(std::vector<unsigned char> &codewords,
                                std::vector<unsigned char> &bytes,
                                unsigned char rdd)
    {
      unsigned char p1, p2, p3;
      for (uint32_t i = 0; i < codewords.size(); i++)
      {
        // first (sf-2) nibbles use CR=4/8
        if (rdd > 2 || i < d_sf - 2) // Hamming(8,4) or Hamming(7,4)
        {
          p1 = parity(codewords[i], (unsigned char)HAMMING_P1_BITMASK);
          p2 = parity(codewords[i], (unsigned char)HAMMING_P2_BITMASK);
          p3 = parity(codewords[i], (unsigned char)HAMMING_P3_BITMASK);
          // p1 covers d1 d2 d4
          // p2 covers d1 d3 d4
          // p3 covers d2 d3 d4
          switch ((p3<<2) | (p2<<1) | p1)
          {
          case 3:
            // p3p2p1 = 011, wrong d1
            codewords[i] ^= HAMMING_D1_BITMASK;
            break;

          case 5:
            // p3p2p1 = 101, wrong d2
            codewords[i] ^= HAMMING_D2_BITMASK;
            break;

          case 6:
            // p3p2p1 = 110, wrong d3
            codewords[i] ^= HAMMING_D3_BITMASK;
            break;

          case 7:
            // p3p2p1 = 111, wrong d4
            codewords[i] ^= HAMMING_D4_BITMASK;
            break;
          
          default:
            // no error, parity error or more than one bit error
            break;
          }
        }
        bytes.push_back(codewords[i] & 0x0F);
      }
    }

    unsigned char
    decode_impl::parity(unsigned char c, unsigned char bitmask)
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
    decode_impl::print_payload(std::vector<unsigned char> &payload)
    {
      std::cout << "Received LoRa packet (hex): ";
      for (int i = 0; i < payload.size(); i++)
      {
        std::cout << std::hex << (unsigned int)payload[i] << " ";
      }
      std::cout << std::endl;
    }

    void
    decode_impl::print_bitwise_u8(std::vector<unsigned char> &buffer)
    {
      for (int i = 0; i < buffer.size(); i++)
      {
        std::cout << i << "\t" << std::bitset<8>(buffer[i] & 0xFF) << "\t";
        std::cout << std::hex << (buffer[i] & 0xFF) << std::endl;
      }
    }

    void
    decode_impl::print_bitwise_u16(std::vector<unsigned short> &buffer)
    {
      for (int i = 0; i < buffer.size(); i++)
      {
        std::cout << i << "\t" << std::bitset<16>(buffer[i] & 0xFFFF) << "\t";
        std::cout << std::hex << (buffer[i] & 0xFFFF) << std::endl;
      }
    }

    void
    decode_impl::decode(pmt::pmt_t msg)
    {
      pmt::pmt_t symbols(pmt::cdr(msg));

      size_t pkt_len(0);
      const uint16_t* symbols_v = pmt::u16vector_elements(symbols, pkt_len);

      std::vector<unsigned short> symbols_in;
      std::vector<unsigned short> header_symbols_in;
      std::vector<unsigned short> payload_symbols_in;
      std::vector<unsigned char> codewords;
      std::vector<unsigned char> nibbles;
      std::vector<unsigned char> combined_bytes;

      short bin_offset = 0;
      unsigned short last_rem;
      unsigned short this_rem;
      bool is_first = true;

      symbols_in.clear();

      for (int i = 0; i < pkt_len; i++)
      {
        unsigned short v = symbols_v[i];
        if (i < 8)
        {
          v = v & ((1<<(d_sf-2))-1);
        }
        // if low data rate optimization is on, give entire packet the header treatment of ppm == SF-2
        else if (d_ldr)
        {
          if (is_first)
          {
            last_rem = v % 4;
            is_first = false;
          }
          this_rem = v % 4;
          // compensate bin drift
          if ((4 + this_rem - last_rem) % 4 == 1) bin_offset -= 1;
          else if ((4 + this_rem - last_rem) % 4 == 3) bin_offset += 1;
          last_rem = this_rem;
          v = (v + (1<<d_sf) + bin_offset) % (1<<d_sf);
          v = v / 4;
        }
        else
        {
          v = (v + (1<<d_sf) - 1) % (1<<d_sf);
        }
        symbols_in.push_back( v );
      }

      to_gray(symbols_in);

      for (int i = 0; (i < symbols_in.size()) && (i < 8); i++) header_symbols_in.push_back(symbols_in[i]);
      for (int i = 8;  i < symbols_in.size()            ; i++) payload_symbols_in.push_back(symbols_in[i]);

      #if DEBUG_OUTPUT
        std::cout << "header syms len " << header_symbols_in.size() << std::endl;
        std::cout << "payload syms len " << payload_symbols_in.size() << std::endl;

        std::cout << "header dewhitened symbols" << std::endl;
        print_bitwise_u16(header_symbols_in);
        std::cout << "payload dewhitened symbols" << std::endl;
        print_bitwise_u16(payload_symbols_in);
      #endif

      // Decode header
      // First 8 symbols are always sent at ppm=d_sf-2, rdd=4 (code rate 4/8), regardless of header mode
      deinterleave(header_symbols_in, codewords, d_sf-2, 4);

      // Decode payload
      // Remaining symbols are at ppm=d_sf, unless sent at the low data rate, in which case ppm=d_sf-2
      deinterleave(payload_symbols_in, codewords, d_ldr ? (d_sf-2) : d_sf, d_cr);
      #if DEBUG_OUTPUT
        std::cout << "deinterleaved codewords" << std::endl;
        print_bitwise_u8(codewords);
      #endif

#if 1 // Disable this #if to derive the whitening sequence

      whiten(codewords);
      #if DEBUG_OUTPUT
        std::cout << "dewhitened codewords" << std::endl;
        print_bitwise_u8(codewords);
      #endif

      hamming_decode(codewords, nibbles, d_cr);
      for (uint32_t i = 0; i < codewords.size(); i++)
      {
        if (i%2 == 0)
        {
          combined_bytes.push_back(nibbles[i] & 0x0F);
        }
        else
        {
          combined_bytes[i/2] |= (nibbles[i] << 4) & 0xF0;
        }
      }

      #if DEBUG_OUTPUT
        std::cout << "data" << std::endl;
        print_bitwise_u8(combined_bytes);
      #endif

      pmt::pmt_t output = pmt::init_u8vector(combined_bytes.size(), combined_bytes);

#else // Whitening sequence derivation

      for (int i = 0; i < codewords.size(); i++)
      {
        std::cout << ", " << std::bitset<8>(codewords[i]);
      }
      std::cout << std::endl;
      std::cout << "Length of above: " << codewords.size() << std::endl;

      pmt::pmt_t output = pmt::init_u8vector(codewords.size(), codewords);

#endif

      pmt::pmt_t msg_pair = pmt::cons(pmt::make_dict(), output);
      message_port_pub(d_out_port, msg_pair);
    }

  } /* namespace lora */
} /* namespace gr */

