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

#define DEBUG_OUTPUT 0

namespace gr {
  namespace lora {

    decode::sptr
    decode::make( int8_t  spreading_factor,
                  bool    header,
                  int16_t payload_len,
                  int8_t  code_rate,
                  bool    crc,
                  bool    low_data_rate)
    {
      return gnuradio::get_initial_sptr
        (new decode_impl(spreading_factor, header, payload_len, code_rate, crc, low_data_rate));
    }

    /*
     * The private constructor
     */
    decode_impl::decode_impl( short spreading_factor,
                              bool  header,
                              short payload_len,
                              short code_rate,
                              bool  crc,
                              bool  low_data_rate)
      : gr::block("decode",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(0, 0, 0)),
        d_sf(spreading_factor),
        d_cr(code_rate),
        d_ldr(low_data_rate),
        d_crc(crc),
        d_payload_len(payload_len),
        d_header(header)
    {
      assert((d_sf > 5) && (d_sf < 13));
      assert((d_cr > 0) && (d_cr < 5));
      if (d_sf == 6) assert(!header);

      d_in_port = pmt::mp("in");
      d_out_port = pmt::mp("out");
      d_header_port = pmt::mp("header");

      message_port_register_in(d_in_port);
      message_port_register_out(d_out_port);
      message_port_register_out(d_header_port);

      set_msg_handler(d_in_port, boost::bind(&decode_impl::decode, this, _1));

      d_whitening_sequence = whitening_sequence;
      if ((d_sf < 6) || (d_sf > 12))
      {
        std::cerr << "Invalid spreading factor -- this state should never occur." << std::endl;
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
    decode_impl::to_gray(std::vector<uint16_t> &symbols)
    {
      for (int i = 0; i < symbols.size(); i++)
      {
        symbols[i] = (symbols[i] >> 1) ^ symbols[i];
      }
    }

    void
    decode_impl::from_gray(std::vector<uint16_t> &symbols)
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
      int offset = d_header ? 3 : 0;
      int crc_offset = d_crc ? 2 : 0;
      for (int i = 0; i + offset < codewords.size() - crc_offset && i < whitening_sequence_length; i++)
      {
        codewords[i+offset] ^= d_whitening_sequence[i];
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
    decode_impl::deinterleave(std::vector <uint16_t> &symbols,
                              std::vector <unsigned char> &codewords,
                              unsigned char ppm,
                              unsigned char rdd)
    {
      const uint32_t bits_per_word = rdd + 4;

      for (uint32_t start_idx = 0; start_idx < symbols.size(); start_idx += bits_per_word) {
        std::vector<uint8_t> block(ppm, 0u);
        for (uint32_t i = 0; i < bits_per_word; i++) {
          const uint32_t word = gr::lora::rotl(symbols[start_idx + i], i, ppm);
          for (uint32_t j = (1u << (ppm - 1)), x = ppm - 1; j; j >>= 1u, x--) {
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
        std::cout << std::hex << (uint32_t)payload[i] << " ";
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
    decode_impl::print_bitwise_u16(std::vector<uint16_t> &buffer)
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
      pmt::pmt_t not_found = pmt::from_bool(false);
      const pmt::pmt_t dict = pmt::car(msg);
      std::string symbol_id = pmt::symbol_to_string(pmt::dict_ref(dict, pmt::intern("id"), not_found));

      pmt::pmt_t symbols(pmt::cdr(msg));

      size_t pkt_len(0);
      const uint16_t* symbols_v = pmt::u16vector_elements(symbols, pkt_len);

      std::vector<uint16_t> symbols_in;
      std::vector<uint16_t> header_symbols_in;
      std::vector<uint16_t> payload_symbols_in;
      std::vector<unsigned char> codewords;
      std::vector<unsigned char> nibbles;
      std::vector<unsigned char> combined_bytes;

      short bin_offset = 0;
      uint16_t last_rem;
      uint16_t this_rem;

      symbols_in.clear();

      for (int i = 0; i < pkt_len; i++)
      {
        uint16_t v = symbols_v[i];
        // if low data rate optimization is on, the ppm of the entire packet is SF-2
        if (i < 8 || d_ldr)
        {
          v /= 4;
        }
        else
        {
          v = gr::lora::pmod(v - 1, 1 << d_sf);
        }
        symbols_in.push_back( v );
      }

      to_gray(symbols_in);

      for (int i = 0; (i < symbols_in.size()) && (i < 8); i++) header_symbols_in.push_back(symbols_in[i]);
      for (int i = 8;  i < symbols_in.size()            ; i++) payload_symbols_in.push_back(symbols_in[i]);

      #if DEBUG_OUTPUT
        std::cout << "header syms len " << header_symbols_in.size() << std::endl;
        std::cout << "payload syms len " << payload_symbols_in.size() << std::endl;

        std::cout << "header symbols" << std::endl;
        print_bitwise_u16(header_symbols_in);
        std::cout << "payload symbols" << std::endl;
        print_bitwise_u16(payload_symbols_in);
      #endif

      // Decode header
      // First 8 symbols are always sent at ppm=d_sf-2, rdd=4 (code rate 4/8), regardless of header mode
      deinterleave(header_symbols_in, codewords, d_sf-2, 4);
      if (d_header) // Explicit Header Mode
      {
        hamming_decode(codewords, nibbles, 4);
        d_payload_len = (nibbles[0] << 4) | nibbles[1];
        d_crc = nibbles[2] & 1;
        d_cr = nibbles[2] >> 1;
        uint8_t checksum = (nibbles[3] << 4) | nibbles[4];
        bool is_valid = true;
        if (checksum != gr::lora::header_checksum(d_payload_len, nibbles[2] & 0xF))
        {
          is_valid = false;
        }
        
        if (symbol_id.substr(0,6) == "header")
        {
          pmt::pmt_t dict = pmt::make_dict();
          dict = pmt::dict_add(dict, pmt::intern("id"), pmt::intern(symbol_id));
          dict = pmt::dict_add(dict, pmt::intern("is_valid"), pmt::from_bool(is_valid));
          dict = pmt::dict_add(dict, pmt::intern("payload_len"), pmt::from_long(d_payload_len));
          dict = pmt::dict_add(dict, pmt::intern("cr"), pmt::from_long(d_cr));
          dict = pmt::dict_add(dict, pmt::intern("crc"), pmt::from_bool(d_crc));
          message_port_pub(d_header_port, dict);
          return;
        }
        if (!is_valid)
        {
          return; // TODO report broken packet
        }
      }

      // Decode payload
      // Remaining symbols are at ppm=d_sf, unless sent at the low data rate, in which case ppm=d_sf-2
      deinterleave(payload_symbols_in, codewords, d_ldr ? (d_sf-2) : d_sf, d_cr);
      #if DEBUG_OUTPUT
        std::cout << "deinterleaved codewords" << std::endl;
        print_bitwise_u8(codewords);
      #endif

      // header has 2.5 bytes, zero-padding to 3 bytes
      if (d_header) codewords.insert(codewords.begin()+5, 0);

      nibbles.clear();
      hamming_decode(codewords, nibbles, d_cr);
      size_t min_len = d_payload_len * 2 + d_header * 6 + d_crc * 4;
      if (nibbles.size() < min_len)
      {
        return; // TODO report broken packet
      }
      for (uint32_t i = 0; i < min_len; i+=2)
      {
        if (d_header && i < 6)
        {
          combined_bytes.push_back((nibbles[i] << 4) | nibbles[i+1]);
        }
        else
        {
          combined_bytes.push_back((nibbles[i+1] << 4) | nibbles[i]);
        }
      }

      #if DEBUG_OUTPUT
        std::cout << "bytes before dewhitening" << std::endl;
        print_bitwise_u8(combined_bytes);
      #endif

#if 1 // Disable this #if to derive the whitening sequence

      whiten(combined_bytes);
      #if DEBUG_OUTPUT
        std::cout << "dewhitened codewords" << std::endl;
        print_bitwise_u8(combined_bytes);
      #endif

      // CRC checksum
      if (d_crc)
      {
        int offset = d_header ? 3 : 0;
        uint16_t checksum = combined_bytes[d_payload_len+offset] | (combined_bytes[d_payload_len+offset+1] << 8);
        combined_bytes.push_back(checksum == gr::lora::data_checksum(&combined_bytes[offset], d_payload_len));
      }

      pmt::pmt_t output = pmt::init_u8vector(combined_bytes.size(), combined_bytes);

#else // Whitening sequence derivation

      for (int i = 0; i < combined_bytes.size(); i++)
      {
        std::cout << ", " << std::bitset<8>(combined_bytes[i]);
      }
      std::cout << std::endl;
      std::cout << "Length of above: " << combined_bytes.size() << std::endl;

      pmt::pmt_t output = pmt::init_u8vector(combined_bytes.size(), combined_bytes);

#endif

      pmt::pmt_t msg_pair = pmt::cons(pmt::make_dict(), output);
      message_port_pub(d_out_port, msg_pair);
    }

  } /* namespace lora */
} /* namespace gr */

