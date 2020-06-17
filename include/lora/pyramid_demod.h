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


#ifndef INCLUDED_LORA_PYRAMID_DEMOD_H
#define INCLUDED_LORA_PYRAMID_DEMOD_H

#include <lora/api.h>
#include <gnuradio/block.h>

#define DEMOD_HISTORY_DEPTH        3
#define REQUIRED_PREAMBLE_CHIRPS   4
#define REQUIRED_SFD_CHIRPS        2
#define LORA_SFD_TOLERANCE         1
#define LORA_PREAMBLE_TOLERANCE    1
#define DEMOD_SYNC_RECOVERY_COUNT  (8-REQUIRED_PREAMBLE_CHIRPS)+(2-REQUIRED_SFD_CHIRPS)+8

namespace gr {
  namespace lora {

    enum pyramid_demod_state_t {
      PS_RESET,
      PS_PREFILL,
      PS_DETECT_PREAMBLE,
      PS_SFD_SYNC,
      PS_READ_HEADER,
      PS_READ_PAYLOAD,
      PS_OUT
    };

    /*!
     * \brief <+description of block+>
     * \ingroup lora
     *
     */
    class LORA_API pyramid_demod : virtual public gr::block
    {
     public:
      typedef boost::shared_ptr<pyramid_demod> sptr;

      /*!
       * \brief Return a shared_ptr to a new instance of lora::pyramid_demod.
       *
       * To avoid accidental use of raw pointers, lora::pyramid_demod's
       * constructor is in a private implementation
       * class. lora::pyramid_demod::make is the public interface for
       * creating new instances.
       */
      static sptr make( unsigned short spreading_factor,
                        bool  low_data_rate,
                        float beta,
                        unsigned short fft_factor,
                        float threshold,
                        float fs_bw_ratio);
    };

  } // namespace lora
} // namespace gr

#endif /* INCLUDED_LORA_PYRAMID_DEMOD_H */

