/*  This file is part of the Vc library.

    Copyright (C) 2009 Matthias Kretz <kretz@kde.org>

    Vc is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    Vc is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Vc.  If not, see <http://www.gnu.org/licenses/>.

*/

/**
 * \mainpage
 * The Vc library is a collection of vector classes with existing implementations for SSE, LRBni or
 * a scalar fallback.
 *
 * \li Vc::float_v
 * \li Vc::sfloat_v
 * \li Vc::double_v
 * \li Vc::int_v
 * \li Vc::uint_v
 * \li Vc::short_v
 * \li Vc::ushort_v
 *
 * \li Vc::float_m
 * \li Vc::sfloat_m
 * \li Vc::double_m
 * \li Vc::int_m
 * \li Vc::uint_m
 * \li Vc::short_m
 * \li Vc::ushort_m
 *
 * \todo
 *  write/link example code, document mask classes, document remaining vector functions
 */

/**
 * \brief Vector Classes
 *
 * Depending on preprocessing macros, the vector classes inside the Vc namespace will be implemented
 * with either
 * \li SSE vectors
 * \li LRBni vectors
 * \li scalar fallback vectors (i.e. scalars, but with the same API)
 */
namespace Vc
{
    /**
     * Enum to declare platform specific constants
     */
    enum {
        /**
         * Specifies the byte boundary for memory alignments necessary for aligned loads and stores.
         */
        VectorAlignment,

        /**
         * Special initializer for vector constructors to create a fast initialization to zero.
         */
        Zero,

        /**
         * Special initializer for vector constructors to create a fast initialization to 1.
         */
        One,

        /**
         * Special initializer for vector constructors to create a vector with the entries 0, 1, 2,
         * 3, 4, 5, ... (depending on the vectors size, of course).
         */
        IndexesFromZero
    };

#define INDEX_TYPE uint_v
#define VECTOR_TYPE float_v
#define ENTRY_TYPE float
#define MASK_TYPE float_m
    /**
     * \class float_v dox.h <Vc/vector.h>
     *
     * SIMD Vector of single precision floats.
     */
    class VECTOR_TYPE
    {
        public:
#include "dox-common-ops.h"
    };
#undef VECTOR_TYPE
#undef ENTRY_TYPE
#undef MASK_TYPE

#define VECTOR_TYPE double_v
#define ENTRY_TYPE double
#define MASK_TYPE double_m
    /**
     * \class double_v dox.h <Vc/vector.h>
     *
     * SIMD Vector of double precision floats.
     */
    class VECTOR_TYPE
    {
        public:
#include "dox-common-ops.h"
    };
#undef VECTOR_TYPE
#undef ENTRY_TYPE
#undef MASK_TYPE

#define VECTOR_TYPE int_v
#define ENTRY_TYPE int
#define MASK_TYPE int_m
    /**
     * \class int_v dox.h <Vc/vector.h>
     *
     * SIMD Vector of 32 bit signed integers.
     */
    class VECTOR_TYPE
    {
        public:
#include "dox-common-ops.h"
    };
#undef VECTOR_TYPE
#undef ENTRY_TYPE
#undef MASK_TYPE

#define VECTOR_TYPE uint_v
#define ENTRY_TYPE unsigned int
#define MASK_TYPE uint_m
    /**
     * \class uint_v dox.h <Vc/vector.h>
     *
     * SIMD Vector of 32 bit unsigned integers.
     */
    class VECTOR_TYPE
    {
        public:
#include "dox-common-ops.h"
    };
#undef VECTOR_TYPE
#undef ENTRY_TYPE
#undef MASK_TYPE
#undef INDEX_TYPE

#define INDEX_TYPE ushort_v
#define VECTOR_TYPE short_v
#define ENTRY_TYPE short
#define MASK_TYPE short_m
    /**
     * \class short_v dox.h <Vc/vector.h>
     *
     * SIMD Vector of 16 bit signed integers.
     *
     * \warning Vectors of this type are not supported on all platforms. In that case the vector
     * class will silently fall back to a Vc::int_v.
     */
    class VECTOR_TYPE
    {
        public:
#include "dox-common-ops.h"
    };
#undef VECTOR_TYPE
#undef ENTRY_TYPE
#undef MASK_TYPE

#define VECTOR_TYPE ushort_v
#define ENTRY_TYPE unsigned short
#define MASK_TYPE ushort_m
    /**
     * \class ushort_v dox.h <Vc/vector.h>
     *
     * SIMD Vector of 16 bit unsigned integers.
     *
     * \warning Vectors of this type are not supported on all platforms. In that case the vector
     * class will silently fall back to a Vc::uint_v.
     */
    class VECTOR_TYPE
    {
        public:
#include "dox-common-ops.h"
    };
#undef VECTOR_TYPE
#undef ENTRY_TYPE
#undef MASK_TYPE

#define VECTOR_TYPE sfloat_v
#define ENTRY_TYPE float
#define MASK_TYPE sfloat_m
    /**
     * \class sfloat_v dox.h <Vc/vector.h>
     *
     * SIMD Vector of single precision floats that is guaranteed to have as many entries as a
     * Vc::short_v / Vc::ushort_v.
     */
    class VECTOR_TYPE
    {
        public:
#include "dox-common-ops.h"
    };
#undef VECTOR_TYPE
#undef ENTRY_TYPE
#undef MASK_TYPE
#undef INDEX_TYPE

}