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

#ifndef VECTOR_H
#define VECTOR_H

#ifdef ENABLE_LARRABEE
#include "larrabee/vector.h"
#define VECTOR_NAMESPACE Larrabee
#elif defined(USE_SSE)
#include "sse/vector.h"
#define VECTOR_NAMESPACE SSE
#else
#include "simple/vector.h"
#define VECTOR_NAMESPACE Simple
#endif

#ifdef isfinite
#undef isfinite
#endif
#ifdef isnan
#undef isnan
#endif

namespace Vc
{
  using VECTOR_NAMESPACE::VectorAlignment;
  using namespace VECTOR_NAMESPACE::VectorSpecialInitializerZero;
  using namespace VECTOR_NAMESPACE::VectorSpecialInitializerOne;
  using namespace VECTOR_NAMESPACE::VectorSpecialInitializerRandom;
  using namespace VECTOR_NAMESPACE::VectorSpecialInitializerIndexesFromZero;
  using VECTOR_NAMESPACE::min;
  using VECTOR_NAMESPACE::max;
  using VECTOR_NAMESPACE::sqrt;
  using VECTOR_NAMESPACE::rsqrt;
  using VECTOR_NAMESPACE::abs;
  using VECTOR_NAMESPACE::sin;
  using VECTOR_NAMESPACE::cos;
  using VECTOR_NAMESPACE::log;
  using VECTOR_NAMESPACE::log10;
  using VECTOR_NAMESPACE::isfinite;
  using VECTOR_NAMESPACE::isnan;
} // namespace Vc

#undef VECTOR_NAMESPACE

#endif // VECTOR_H