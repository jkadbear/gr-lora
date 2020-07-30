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

#ifndef UTILITIES_H
#define UTILITIES_H

namespace gr {
  namespace lora {

    /**
     *  \brief  Python-style modulo, (i % n) >= 0
     *
     *  \param  x
     *          dividend
     *  \param  n
     *          divisor
     */
    inline unsigned int pos_mod(int x, unsigned int n)
    {
      return ((x % n) + n) % n;
    }

    /**
     *  \brief  Rotate the given bits to the left and return the result.
     *
     *  \param  bits
     *          The value to rotate.
     *  \param  count
     *          The amount of bits to rotate (shift to left and add to right).
     *  \param  size
     *          The size in bits used in `bits`.
     *          <BR>e.g. 1 byte given       => size = 8
     *          <BR>e.g. only 6 bits in use => size = 6, and all bits higher than (1 << size-1) will be zeroed.
     */
    inline uint32_t rotl(uint32_t bits, uint32_t count = 1u, const uint32_t size = 8u) {
        const uint32_t len_mask = (1u << size) - 1u;

        count %= size;      // Limit bit rotate count to size
        bits  &= len_mask;  // Limit given bits to size

        return ((bits << count) & len_mask) | (bits >> (size - count));
    }

    inline uint16_t crc16sx(uint16_t crc, const uint16_t poly) {
        for (int i = 0; i < 8; i++) {
            if (crc & 0x8000) {
                crc = (crc << 1) ^ poly;
            }
            else {
                crc <<= 1;
            }
        }
        return crc;
    }

    inline uint8_t xsum8(uint8_t t) {
        t ^= t >> 4;
        t ^= t >> 2;
        t ^= t >> 1;
        return (t & 1);
    }

    /***********************************************************************
     *  CRC reverse engineered from Sx1272 data stream.
     *  Modified CCITT crc with masking of the output with an 8bit lfsr
     **********************************************************************/
    inline uint16_t data_checksum(const uint8_t *data, int length) {
        uint16_t res = 0;
        uint8_t v = 0xff;
        uint16_t crc = 0;
        for (int i = 0; i < length; i++) {
            crc = crc16sx(res, 0x1021);
            v = xsum8(v & 0xB8) | (v << 1);
            res = crc ^ data[i];
        }
        res ^= v; 
        v = xsum8(v & 0xB8) | (v << 1);
        res ^= v << 8;
        return res;
    }

    inline uint8_t header_checksum(const uint16_t len, const uint8_t cr_crc) {
        auto a0 = (len >> 4) & 0x1;
        auto a1 = (len >> 5) & 0x1;
        auto a2 = (len >> 6) & 0x1;
        auto a3 = (len >> 7) & 0x1;

        auto b0 = (len >> 0) & 0x1;
        auto b1 = (len >> 1) & 0x1;
        auto b2 = (len >> 2) & 0x1;
        auto b3 = (len >> 3) & 0x1;

        auto c0 = (cr_crc >> 0) & 0x1;
        auto c1 = (cr_crc >> 1) & 0x1;
        auto c2 = (cr_crc >> 2) & 0x1;
        auto c3 = (cr_crc >> 3) & 0x1;

        uint8_t res;
        res = (a0 ^ a1 ^ a2 ^ a3) << 4;
        res |= (a3 ^ b1 ^ b2 ^ b3 ^ c0) << 3;
        res |= (a2 ^ b0 ^ b3 ^ c1 ^ c3) << 2;
        res |= (a1 ^ b0 ^ b2 ^ c0 ^ c1 ^ c2) << 1;
        res |= a0 ^ b1 ^ c0 ^ c1 ^ c2 ^ c3;

        return res;
    }
  }
}


#endif /* UTILITIES_H */
