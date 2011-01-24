/*  This file is part of the Vc library.

    Copyright (C) 2010 Matthias Kretz <kretz@kde.org>

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

#ifndef VC_SSE_DEINTERLEAVE_H
#define VC_SSE_DEINTERLEAVE_H

namespace Vc
{
namespace Internal
{

template<> struct HelperImpl<Vc::SSE2Impl>
{
    static inline void prefetchForOneRead(const void *addr) ALWAYS_INLINE;
    static inline void prefetchForModify(const void *addr) ALWAYS_INLINE;
    static inline void prefetchClose(const void *addr) ALWAYS_INLINE;
    static inline void prefetchMid(const void *addr) ALWAYS_INLINE;
    static inline void prefetchFar(const void *addr) ALWAYS_INLINE;
};

template<> struct HelperImpl<SSE3Impl> : public HelperImpl<SSE2Impl> {};
template<> struct HelperImpl<SSSE3Impl> : public HelperImpl<SSE3Impl> {};
template<> struct HelperImpl<SSE41Impl> : public HelperImpl<SSSE3Impl> {};
template<> struct HelperImpl<SSE42Impl> : public HelperImpl<SSE41Impl> {};
template<> struct HelperImpl<SSE4aImpl> : public HelperImpl<SSE3Impl> {};

} // namespace Internal
} // namespace Vc

#include "prefetches.tcc"
#endif // VC_SSE_DEINTERLEAVE_H
