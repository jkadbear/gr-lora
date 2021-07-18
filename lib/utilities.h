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
     *  \brief  Python-style modulo, (x % n) >= 0
     *
     *  \param  x
     *          dividend
     *  \param  n
     *          divisor
     */
    inline uint32_t pmod(int32_t x, int32_t n)
    {
      return ((x % n) + n) % n;
    }

    /**
     *  \brief  Python-style float modulo, (x % n) >= 0
     *
     *  \param  x
     *          dividend
     *  \param  n
     *          divisor
     */
    inline float fpmod(float x, float n)
    {
      return std::fmod(std::fmod(x, n) + n, n);
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

    inline uint16_t data_checksum(const uint8_t *data, int length) {
      uint16_t crc = 0;
      for (int j = 0; j < length - 2; j++)
      {
        uint8_t new_byte = data[j];
        for (int i = 0; i < 8; i++) {
          if (((crc & 0x8000) >> 8) ^ (new_byte & 0x80)) {
            crc = (crc << 1) ^ 0x1021;
          } else {
            crc = (crc << 1);
          }
          new_byte <<= 1;
        }
      }

      // XOR the obtained CRC with the last 2 data bytes
      uint16_t x1 = (length >= 1) ? data[length-1]      : 0;
      uint16_t x2 = (length >= 2) ? data[length-2] << 8 : 0;
      crc = crc ^ x1 ^ x2;
      return crc;
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

    inline void linear_regression(float *h, int start_idx, int end_idx, float *k, float *b)
    {
      // std::cout << "idx: ";
      // for (size_t i = start_idx; i <= end_idx; i++)
      // {
      //   std::cout << i << ",";
      // }
      // std::cout << std::endl << "h: ";
      // for (size_t i = start_idx; i <= end_idx; i++)
      // {
      //   std::cout << h[i] << ",";
      // }
      // std::cout << std::endl;
      
      float avg_x = 0, avg_y = 0;
      for (size_t i = start_idx; i <= end_idx; i++)
      {
        avg_x += i;
        avg_y += h[i];
      }
      avg_x /= (end_idx - start_idx + 1);
      avg_y /= (end_idx - start_idx + 1);

      float numerator = 0, denominator = 0;
      for (size_t i = start_idx; i <= end_idx; i++)
      {
        numerator += (i - avg_x) * (h[i] - avg_y);
        denominator += (i - avg_x) * (i - avg_x);
      }
      
      if (denominator == 0)
      {
        std::cerr << "Invalid denominator -- this state should never occur." << std::endl;
      }

      *k = numerator / denominator;
      *b = avg_y - *k * avg_x;
    }

    inline uint32_t argmax_32f(float *res, float *p_max_val, uint32_t len)
    {
      float mag = abs(res[0]);
      uint32_t max_idx = 0;

      *p_max_val = mag;
      for (int i = 1; i < len; i++)
      {
        mag = abs(res[i]);
        if (mag > *p_max_val)
        {
          max_idx = i;
          *p_max_val = mag;
        }
      }

      return max_idx;
    }
  }
}


#endif /* UTILITIES_H */
