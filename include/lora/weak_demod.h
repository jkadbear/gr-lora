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

#ifndef INCLUDED_LORA_WEAK_DEMOD_H
#define INCLUDED_LORA_WEAK_DEMOD_H

#include <lora/api.h>
#include <gnuradio/block.h>

#define WEAK_REQUIRED_PREAMBLE_CHIRPS   5
#define WEAK_DEMOD_BUFFER_SIZE          15
#define WEAK_DEMOD_HISTORY              7
#define WEAK_DEMOD_SYNC_RECOVERY_COUNT  7

namespace gr {
  namespace lora {

    enum weak_demod_state_t {
      WS_RESET,
      WS_PREFILL,
      WS_DETECT_PREAMBLE,
      WS_SFD_SYNC,
      WS_READ_HEADER,
      WS_READ_PAYLOAD,
      WS_OUT
    };

    /*!
     * \brief <+description of block+>
     * \ingroup lora
     *
     */
    class LORA_API weak_demod : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<weak_demod> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of lora::weak_demod.
       *
       * To avoid accidental use of raw pointers, lora::weak_demod's
       * constructor is in a private implementation
       * class. lora::weak_demod::make is the public interface for
       * creating new instances.
       */
      static sptr make(uint8_t   spreading_factor,
                  bool      header,
                  uint8_t   payload_len,
                  uint8_t   cr,
                  bool      crc,
                  bool      low_data_rate,
                  uint32_t  d_sym_num,
                  float     beta,
                  uint16_t  fft_factor,
                  uint8_t   peak_search_algorithm,
                  uint16_t  peak_search_phase_k,
                  float     fs_bw_ratio);
    };

  } // namespace lora
} // namespace gr

#endif /* INCLUDED_LORA_WEAK_DEMOD_H */

