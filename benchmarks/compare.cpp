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

#include <Vc/Vc>
#include "benchmark.h"
#include "random.h"
#include <cstdio>
#include <cstdlib>

using namespace Vc;

bool blackHole = false;
float_m floatResult;
short_m shortResult;
#if VC_IMPL_SSE
sfloat_m sfloatResult;
#endif

#define unrolled_loop4(_it_, _code_) \
{ enum { _it_ = 0 }; _code_ } \
{ enum { _it_ = 1 }; _code_ } \
{ enum { _it_ = 2 }; _code_ } \
{ enum { _it_ = 3 }; _code_ }

template<typename Vector> class DoCompares
{
    enum {
        Factor = 5120000 / Vector::Size
    };

    public:
        DoCompares(const int Repetitions)
            : a(new Vector[Factor]),
            b(new Vector[Factor])
        {
            setResultPointer();
            for (int i = 0; i < Factor; ++i) {
                a[i] = PseudoRandom<Vector>::next();
                b[i] = PseudoRandom<Vector>::next();
            }

            {
                Benchmark timer("operator<", Vector::Size * Factor, "Op");
                doWork1();
                for (int repetitions = 0; repetitions < Repetitions; ++repetitions) {
                    timer.Start();
                    doWork1();
                    timer.Stop();
                }
                timer.Print(Benchmark::PrintAverage);
            }
            {
                Benchmark timer("masked assign with operator==", Vector::Size * Factor, "Op");
                doWork2();
                for (int repetitions = 0; repetitions < Repetitions; ++repetitions) {
                    timer.Start();
                    doWork2();
                    timer.Stop();
                    for (int i = 0; i < Factor; ++i) {
                        *result = a[i] > Vector(One);
                    }
                }
                timer.Print(Benchmark::PrintAverage);
            }
            {
                Benchmark timer("(operator==).isFull()", Vector::Size * Factor, "Op");
                doWork3();
                for (int repetitions = 0; repetitions < Repetitions; ++repetitions) {
                    timer.Start();
                    doWork3();
                    timer.Stop();
                }
                timer.Print(Benchmark::PrintAverage);
            }
                {
                Benchmark timer("!(operator==).isEmpty()", Vector::Size * Factor, "Op");
                doWork4();
                for (int repetitions = 0; repetitions < Repetitions; ++repetitions) {
                    timer.Start();
                    doWork4();
                    timer.Stop();
                }
                timer.Print(Benchmark::PrintAverage);
            }
        }

        ~DoCompares()
        {
            delete[] a;
            delete[] b;
        }

    private:
        void setResultPointer();
        void doWork1();
        void doWork2();
        void doWork3();
        void doWork4();

        Vector *a;
        Vector *b;
        typename Vector::Mask *result;
};

template<> inline void DoCompares<float_v>::setResultPointer() { result = &floatResult; }
template<> inline void DoCompares<short_v>::setResultPointer() { result = &shortResult; }
#if VC_IMPL_SSE
template<> inline void DoCompares<sfloat_v>::setResultPointer() { result = &sfloatResult; }
#endif

template<typename Vector> inline void DoCompares<Vector>::doWork1()
{
    for (int i = 0; i < Factor; ++i) {
        *result = a[i] < b[i];
    }
}
template<typename Vector> inline void DoCompares<Vector>::doWork2()
{
    const Vector one(One);
    for (int i = 0; i < Factor; ++i) {
        a[i](a[i] == b[i]) = one;
    }
}
template<typename Vector> inline void DoCompares<Vector>::doWork3()
{
    const Vector one(One);
    for (int i = 0; i < Factor; ++i) {
        blackHole = (a[i] == b[i]).isFull();
    }
}
template<typename Vector> inline void DoCompares<Vector>::doWork4()
{
    const Vector one(One);
    for (int i = 0; i < Factor; ++i) {
        blackHole = !(a[i] == b[i]).isEmpty();
    }
}

int bmain(Benchmark::OutputMode out)
{
    const int Repetitions = out == Benchmark::Stdout ? 10 : g_Repetitions > 0 ? g_Repetitions : 100;

    Benchmark::addColumn("datatype");

    Benchmark::setColumnData("datatype", "float_v");
    DoCompares<float_v> a(Repetitions);
    Benchmark::setColumnData("datatype", "short_v");
    DoCompares<short_v> b(Repetitions);
#if VC_IMPL_SSE
    Benchmark::setColumnData("datatype", "sfloat_v");
    DoCompares<sfloat_v> c(Repetitions);
#endif

    return 0;
}
