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
    static constexpr int Iterations = 1000;

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
            //std::cout << count << std::endl;
            m_state = 0;
            m_cycles = 0.;
            std::thread poison(poisonThreadCaller, this, 1, count);
            std::thread bench(benchmarkThreadCaller, this, 0, count);
            poison.join();
            bench.join();
            std::cout << "required " << m_cycles / Iterations << " cycles.\n";
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
        while (m_state < 2 * Iterations) {
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
        while (m_state < 2 * Iterations) {
            unique_lock<mutex> lock(m_mutex);
            if ((m_state & 1) == 0) {
                m_wait.wait(lock);
            }
            tsc.Start();
            for (size_t i = 0; i < count; i += 8) {
                const __m128d tmp0 = _mm_add_pd(two, _mm_load_pd(&mem[i + 0]));
                //const __m128d tmp1 = _mm_mul_pd(two, _mm_load_pd(&mem[i + 2]));
                const __m128d tmp2 = _mm_add_pd(two, _mm_load_pd(&mem[i + 4]));
                //const __m128d tmp3 = _mm_mul_pd(two, _mm_load_pd(&mem[i + 6]));
                _mm_store_pd(&mem[i + 0], tmp0);
                //_mm_store_pd(&mem[i + 2], tmp1);
                _mm_store_pd(&mem[i + 4], tmp2);
                //_mm_store_pd(&mem[i + 6], tmp3);
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
