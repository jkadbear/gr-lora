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
     *  \param  i
     *          dividend
     *  \param  n
     *          divisor
     */
    inline unsigned int
    pos_mod(int x, unsigned int n)
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
  }
}


#endif /* UTILITIES_H */
