#include <Vc/Vc>/*{{{*/
#include "benchmark.h"
#include <sys/mman.h>
#include <fstream>
#include <thread>
#include <mutex>
#include <atomic>
#include "../cpuid.h"
#include "cpuset.h"

using namespace Vc;/*}}}*/

using std::mutex;
using std::unique_lock;

static void pinToCpu(int id)
{
    cpu_set_t cpumask;
    cpuZero(&cpumask);
    cpuSet(id, &cpumask);
    sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);
}

class Add1Test
{
    std::mutex m_mutex;
    std::condition_variable m_wait;
    int m_state;
    double m_cycles;
    double *__restrict__ mem;
    int m_iterations;
    static int maxThreadCount() { cpu_set_t cpumask; sched_getaffinity(0, sizeof(cpu_set_t), &cpumask); return cpuCount(&cpumask); }

public:
    Add1Test()
        : m_state(0),
        mem(new double[1024 * 1024 * 64])
    {
        mlockall(MCL_CURRENT);

        const size_t counts[] = {
            CpuId::L1Data() / sizeof(double),
            CpuId::L2Data() / sizeof(double),
            CpuId::L3Data() / sizeof(double)
        };
        for (size_t count : counts) {
            m_iterations = 8000 * 1024 / count;
            const int benchCpu = 0;
            for (int poisonCpu = 0; poisonCpu < maxThreadCount(); ++poisonCpu) {
                m_state = 0;
                m_cycles = 0.;
                std::thread poison(poisonThreadCaller, this, poisonCpu, count);
                std::thread bench(benchmarkThreadCaller, this, benchCpu, count);
                poison.join();
                bench.join();
                std::cout << poisonCpu << /*'\t' << benchCpu <<*/ '\t' << count * sizeof(double) / 1024
                    << '\t' << m_cycles / m_iterations << '\n';
            }
        }
    }

    static void poisonThreadCaller(Add1Test *that, int cpuId, size_t count)
    {
        that->poisonThread(cpuId, count);
    }
    static void benchmarkThreadCaller(Add1Test *that, int cpuId, size_t count)
    {
        that->benchmarkThread(cpuId, count);
    }

    void poisonThread(int cpuId, size_t count)
    {
        pinToCpu(cpuId);
        while (m_state < 2 * m_iterations) {
            unique_lock<mutex> lock(m_mutex);
            if ((m_state & 1) == 1) {
                m_wait.wait(lock);
            }
            // mark the cachelines as modified
            for (size_t i = 0; i < count; i += 4) {
                mem[i] += 1.;
            }
            ++m_state;
            m_wait.notify_one();
        }
    }

    void benchmarkThread(int cpuId, size_t count)
    {
        pinToCpu(cpuId);
        TimeStampCounter tsc;
        const __m128d two = _mm_set1_pd(2.);
        while (m_state < 2 * m_iterations) {
            unique_lock<mutex> lock(m_mutex);
            if ((m_state & 1) == 0) {
                m_wait.wait(lock);
            }
            tsc.Start();
            for (size_t i = 0; i < count; i += 16) {
                const __m128d tmp0 = _mm_add_pd(two, _mm_load_pd(&mem[i +  0]));
                const __m128d tmp2 = _mm_add_pd(two, _mm_load_pd(&mem[i +  4]));
                const __m128d tmp4 = _mm_add_pd(two, _mm_load_pd(&mem[i +  8]));
                const __m128d tmp6 = _mm_add_pd(two, _mm_load_pd(&mem[i + 12]));
                _mm_store_pd(&mem[i +  0], tmp0);
                _mm_store_pd(&mem[i +  4], tmp2);
                _mm_store_pd(&mem[i +  8], tmp4);
                _mm_store_pd(&mem[i + 12], tmp6);
            }
            tsc.Stop();
            ++m_state;
            m_wait.notify_one();
            m_cycles += tsc.Cycles();
        }
    }
};

int bmain()/*{{{*/
{
    Add1Test();
    return 0;
}
