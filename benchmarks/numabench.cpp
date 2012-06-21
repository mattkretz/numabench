/*
    Copyright (C) 2010 Matthias Kretz <kretz@kde.org>

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of
    the License, or (at your option) version 3.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
    02110-1301, USA.

*/

#include <Vc/Vc>/*{{{*/
#include "benchmark.h"
#include <sys/mman.h>
#include <fstream>
#include <thread>
#include <mutex>
#include <atomic>
#include "../cpuid.h"
#ifdef NO_LIBNUMA
#include "cpuset.h"
#else
#include "numa.h"
#endif

extern "C" void numa_init();

using namespace Vc;/*}}}*/
enum Constants {/*{{{*/
    PageSize = 4096,
    CacheLineSize = 64,

    DoublesInPage = PageSize / sizeof(double),
    VectorsInPage = PageSize / sizeof(double_v),
    DoublesInCacheLine = CacheLineSize / sizeof(double),
    VectorsInCacheLine = CacheLineSize / sizeof(double_v)
};
static const size_t GB = 1024ull * 1024ull * 1024ull;
static const size_t step = GB / sizeof(double_v);/*}}}*/
typedef Memory<double_v> MemT;
struct TestArguments/*{{{*/
{
    MemT *__restrict__ mem;
    Timer *__restrict__ timer;
    size_t offset;
    size_t size;
    int repetitions;
};/*}}}*/
typedef void (*TestFunction)(const TestArguments &args);
struct CpuRange/*{{{*/
{
    int first;
    int last;
    int step;
};/*}}}*/
class OneWaitsForN/*{{{*/
{
private:
    std::mutex mutex;
    std::condition_variable wait;
    std::atomic<int> busyCount;
public:
    OneWaitsForN(int count)
        : busyCount(count)
    {
    }

    void oneReady()
    {
        if (busyCount-- == 1) {
            mutex.lock();
            wait.notify_one();
            mutex.unlock();
        }
    }

    void waitForAll()
    {
        if (busyCount > 0) {
            std::unique_lock<std::mutex> lock(mutex);
            if (busyCount > 0) {
                wait.wait(lock);
            }
        }
    }
};/*}}}*/
class ThreadData/*{{{*/
{
    cpu_set_t m_cpumask;
    std::mutex m_mutex;
    std::condition_variable_any &m_wait;
    OneWaitsForN &m_waitForEnd;
    std::atomic<bool> m_exit;
    bool m_disabled;
    Timer m_timer;
    std::thread m_thread;
    TestFunction m_testFunction;
    TestArguments m_arguments;

    public:
        ThreadData(OneWaitsForN *waitForEnd, std::condition_variable_any *wait) // called from main thread
            : m_wait(*wait),
            m_waitForEnd(*waitForEnd),
            m_exit(false),
            m_disabled(false),
            m_thread(ThreadData::callMainLoop, this)
        {}

        const Timer &timer() const { return m_timer; }

        void setPinning(int cpuid) // called from main thread
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_disabled = false;
            cpuZero(&m_cpumask);
            cpuSet(cpuid, &m_cpumask);
        }

        void disable()
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_disabled = true;
        }
        void enable()
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_disabled = false;
        }
        bool isEnabled() const
        {
            return !m_disabled;
        }

        void exit() // called from main thread
        {
            m_exit = true;
        }

        void join() // called from main thread
        {
            m_thread.join();
        }

        void setTestFunction(TestFunction f)
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_testFunction = f;
        }

        void setParameters(TestArguments args)
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_arguments = args;
            m_arguments.timer = &m_timer;
        }

        static void callMainLoop(ThreadData *data) // thread
        {
            data->mainLoop();
        }

    private:
        void mainLoop() // thread
        {
            m_mutex.lock();
            do {
                m_waitForEnd.oneReady();

                // wait for the signal to start
                m_wait.wait(m_mutex);

                if (m_exit) {
                    break;
                } else if (m_disabled) {
                    continue;
                }

                // first pin the thread to a single core/cpu
                sched_setaffinity(0, sizeof(cpu_set_t), &m_cpumask);

                // do the work
                m_testFunction(m_arguments);
            } while (!m_exit);
            m_waitForEnd.oneReady();
            m_mutex.unlock();
        }
};/*}}}*/
class ThreadPool/*{{{*/
{
    OneWaitsForN m_waitForEnd;
    std::condition_variable_any m_waitForStart;
    std::vector<std::shared_ptr<ThreadData>> m_workers;
    static int maxThreadCount() { cpu_set_t cpumask; sched_getaffinity(0, sizeof(cpu_set_t), &cpumask); return cpuCount(&cpumask); }

public:
    ThreadPool(int _size = maxThreadCount())
        : m_waitForEnd(_size),
        m_workers(_size)
    {
        for (int i = 0; i < _size; ++i) {
            m_workers[i] = std::make_shared<ThreadData>(&m_waitForEnd, &m_waitForStart);
        }
    }

    void waitReady()
    {
        m_waitForEnd.waitForAll();
    }

    void setPinning(int firstCpu)
    {
        for (size_t i = 0; i < m_workers.size(); ++i) {
            m_workers[i]->setPinning(firstCpu + i);
        }
    }
    void setPinning(std::vector<int> cpus)
    {
        assert(m_workers.size() == cpus.size());
        for (size_t i = 0; i < m_workers.size(); ++i) {
            m_workers[i]->setPinning(cpus[i]);
        }
    }
    void setPinning(CpuRange range)
    {
        unsigned int i = 0;
        for (int id = range.first; id <= range.last; id += range.step) {
            m_workers[i++]->setPinning(id);
        }
        for (; i < m_workers.size(); ++i) {
            m_workers[i]->disable();
        }
    }

    void setTestFunction(TestFunction f)
    {
        for (auto &t : m_workers) {
            t->setTestFunction(f);
        }
    }

    template<typename OffsetFunction>
    void executeWith(MemT *mem, OffsetFunction offset, size_t _size, int repetitions)
    {
        for (auto &t : m_workers) {
            t->setParameters({ mem, nullptr, offset(), _size, repetitions });
        }
        m_waitForStart.notify_all();
    }

    template<typename F> void eachTimer(F f) const
    {
        for (const auto &t : m_workers) {
            if (t->isEnabled()) {
                f(t->timer());
            }
        }
    }

    ~ThreadPool()
    {
        for (auto &t : m_workers) {
            t->exit();
        }
        m_waitForStart.notify_all();
        for (auto &t : m_workers) {
            t->join();
        }
    }
};/*}}}*/
static size_t largestMemorySize()/*{{{*/
{
    using namespace std;
    fstream meminfo("/proc/meminfo", fstream::in);
    string tmp;
    size_t totalMem, freeMem;
    meminfo >> tmp >> totalMem >> tmp >> tmp >> freeMem;
    meminfo.close();
    return freeMem * 1024;
}/*}}}*/
void testBzero(const TestArguments &args)/*{{{*/
{
    args.timer->start();
    for (int rep = 0; rep < args.repetitions; ++rep) {
        bzero((*args.mem) + args.offset * double_v::Size, args.size * sizeof(double_v));
    }
    args.timer->stop();
}/*}}}*/
void testAddOne(const TestArguments &args)/*{{{*/
{
    const double_v one = 1.;

    double *__restrict__ mStart = args.mem->entries() + args.offset * double_v::Size;

    args.timer->start();
    for (int rep = 0; rep < args.repetitions; ++rep) {
        double *__restrict__ m = mStart;
        for (size_t i = 0; i < args.size; i += 4) {
            (double_v(m + 0) + one).store(m + 0);
            (double_v(m + 2) + one).store(m + 2);
            (double_v(m + 4) + one).store(m + 4);
            (double_v(m + 6) + one).store(m + 6);
            m += 8;
        }
    }
    args.timer->stop();
}/*}}}*/
void testAddOnePrefetch(const TestArguments &args)/*{{{*/
{
    const double_v one = 1.;

    double *__restrict__ mStart = args.mem->entries() + args.offset * double_v::Size;

    args.timer->start();
    for (int rep = 0; rep < args.repetitions; ++rep) {
        double *__restrict__ m = mStart;
        for (size_t i = 0; i < args.size; i += 4) {
            Vc::prefetchForModify(m + 1024);
            (double_v(m + 0) + one).store(m + 0);
            (double_v(m + 2) + one).store(m + 2);
            (double_v(m + 4) + one).store(m + 4);
            (double_v(m + 6) + one).store(m + 6);
            m += 8;
        }
    }
    args.timer->stop();
}/*}}}*/
void testRead(const TestArguments &args)/*{{{*/
{
    double *__restrict__ mStart = args.mem->entries() + args.offset * double_v::Size;
    args.timer->start();
    for (int rep = 0; rep < args.repetitions; ++rep) {
        double *__restrict__ m = mStart;
        for (size_t i = 0; i < args.size; i += 4) {
            const double_v v0(m + 0);
            const double_v v1(m + 2);
            const double_v v2(m + 4);
            const double_v v3(m + 6);
            asm("" :: "x"(v0.data()), "x"(v1.data()), "x"(v2.data()), "x"(v3.data()));
            m += 8;
        }
    }
    args.timer->stop();
}/*}}}*/
void testReadPrefetch(const TestArguments &args)/*{{{*/
{
    double *__restrict__ mStart = args.mem->entries() + args.offset * double_v::Size;
    args.timer->start();
    for (int rep = 0; rep < args.repetitions; ++rep) {
        double *__restrict__ m = mStart;
        for (size_t i = 0; i < args.size; i += 4) {
            Vc::prefetchForOneRead(m + 1024);
            const double_v v0(m + 0);
            const double_v v1(m + 2);
            const double_v v2(m + 4);
            const double_v v3(m + 6);
            asm("" :: "x"(v0.data()), "x"(v1.data()), "x"(v2.data()), "x"(v3.data()));
            m += 8;
        }
    }
    args.timer->stop();
}/*}}}*/
/** testReadLatency {{{
 * We want to measure the latency of a read from memory. To achieve this we read with a stride of
 * PageSize bytes. Then the hardware prefetcher will not do any prefetches and every load will hit a
 * cold cache line. To increase the working size the test then starts over but with an offset of one
 * cache line:
 * [x                               x                               ...]
 * [        x                               x                       ...]
 * [                x                               x               ...]
 * [                        x                               x       ...]
 */
void testReadLatency(const TestArguments &args)
{
    typedef double *__restrict__ Ptr;
    Ptr const mStart = args.mem->entries() + args.offset * double_v::Size;
    Ptr const mEnd = mStart + args.size * 2;
    Ptr const mPageEnd = mStart + DoublesInPage;
    args.timer->start();
    if (((mEnd - mStart) / DoublesInCacheLine) & 1) {
        for (int rep = 0; rep < args.repetitions; ++rep) {
            for (Ptr mCacheLine = mStart; mCacheLine < mPageEnd; mCacheLine += DoublesInCacheLine) {
                for (Ptr m = mCacheLine; m < mEnd; m += DoublesInPage) {
                    //asm volatile("lfence");
                    asm volatile("" :: "d"(*m));
                }
            }
        }
    } else if (((mEnd - mStart) / DoublesInCacheLine) & 2) {
        for (int rep = 0; rep < args.repetitions; ++rep) {
            for (Ptr mCacheLine = mStart; mCacheLine < mPageEnd; mCacheLine += DoublesInCacheLine) {
                for (Ptr m = mCacheLine; m < mEnd - DoublesInPage; m += 2 * DoublesInPage) {
                    //asm volatile("lfence");
                    asm volatile("" :: "d"(*m));
                    asm volatile("" :: "d"(*(m + DoublesInPage)));
                }
            }
        }
    } else {
        for (int rep = 0; rep < args.repetitions; ++rep) {
            for (Ptr mCacheLine = mStart; mCacheLine < mPageEnd; mCacheLine += DoublesInCacheLine) {
                for (Ptr m = mCacheLine; m < mEnd - 3 * DoublesInPage; m += 4 * DoublesInPage) {
                    //asm volatile("lfence");
                    asm volatile("" :: "d"(*m));
                    asm volatile("" :: "d"(*(m + DoublesInPage)));
                    asm volatile("" :: "d"(*(m + 2 * DoublesInPage)));
                    asm volatile("" :: "d"(*(m + 3 * DoublesInPage)));
                }
            }
        }
    }
    args.timer->stop();
}/*}}}*/
template<typename T> struct convertStringTo/*{{{*/
{
    explicit convertStringTo(const std::string &s);
    operator T() { return m_data; }
    T m_data;
};
/*}}}*/
template<> convertStringTo<int>::convertStringTo(const std::string &s) : m_data(atoi(s.c_str())) {}/*{{{*/
template<> convertStringTo<unsigned int>::convertStringTo(const std::string &s) : m_data(atoi(s.c_str())) {}
template<> convertStringTo<long>::convertStringTo(const std::string &s) : m_data(atol(s.c_str())) {}
template<> convertStringTo<unsigned long>::convertStringTo(const std::string &s) : m_data(atol(s.c_str())) {}
template<> convertStringTo<long long>::convertStringTo(const std::string &s) : m_data(atoll(s.c_str())) {}
template<> convertStringTo<unsigned long long>::convertStringTo(const std::string &s) : m_data(atoll(s.c_str())) {}
template<> convertStringTo<std::string>::convertStringTo(const std::string &s) : m_data(s) {}/*}}}*/
template<typename T> static T valueForArgument(const char *name, T defaultValue)/*{{{*/
{
    ArgumentVector::iterator it = std::find(g_arguments.begin(), g_arguments.end(), name);
    if (it != g_arguments.end()) {
        ++it;
        if (it != g_arguments.end()) {
            return convertStringTo<T>(*it);
        }
    }
    return defaultValue;
}/*}}}*/
/*SET_HELP_TEXT{{{*/
#ifdef NO_LIBNUMA
SET_HELP_TEXT(
        "  --firstCpu <id>\n"
        "  --cpuStep <id>\n"
        "  --size <GB>\n"
        "  --only <test function>\n"
        "  --threads <N>\n"
        "  --cores <firstId-lastId[:step][,firstId-lastId[:step][...]]>\n"
        );
#else
SET_HELP_TEXT(
        "  --firstNode <id>\n"
        "  --nodeStep <id>\n"
        "  --size <GB>\n"
        "  --only <test function>\n"
        "  --threads <N>\n"
        "  --cores <firstId-lastId[:step][,firstId-lastId[:step][...]]>\n"
        );
#endif/*}}}*/
class BenchmarkRunner/*{{{*/
{
private:
    const std::vector<CpuRange> m_coreIds;
    int m_threadCount;
    ThreadPool m_threadPool;
    const size_t m_maxMemorySize;
    const size_t m_memorySize;
    const std::string m_only;
    MemT *m_memory;

    void executeTest(const char *name, TestFunction testFun, const char *unit, double interpretFactor);
    void executeAllTests();

public:
    BenchmarkRunner();
};/*}}}*/
void BenchmarkRunner::executeTest(const char *name, TestFunction testFun, const char *unit, double interpretFactor)/*{{{*/
{
    m_threadPool.setTestFunction(testFun);
    if (m_only.empty() || m_only == name) {
        for (size_t offset = 0; offset <= m_memory->vectorsCount() - step; offset += step) {
            std::stringstream ss;
            ss << name << ": " << offset * sizeof(double_v) / GB << " - "
                << (offset + step) * sizeof(double_v) / GB;
            {
                m_threadPool.waitReady();
                Benchmark bench(ss.str().c_str(), step * m_threadCount * interpretFactor, unit);
                Timer timer;
                for (int rep = 0; rep < 2; ++rep) {
                    m_threadPool.executeWith(m_memory, [offset]() -> size_t { return offset; }, step / 16, 1);
                    testFun({m_memory, &timer, offset, step / 16, 1});
                    m_threadPool.waitReady();
                    m_threadPool.eachTimer([&bench](const Timer &t) { bench.addTiming(t); });
                    bench.addTiming(timer);
                }
                bench.Print();
            }
        }
    }

    size_t sizes[] = {
        CpuId::L1Data(),
        CpuId::L2Data(),
        CpuId::L3Data()
    };
    for (int i = 0; i < 3; ++i) {
        size_t size = sizes[i];
        if (size > 0) {
            size /= 2;
            std::stringstream ss0;
            ss0 << name << " (" << size / 1024 << "kB)";
            if (m_only.empty() || m_only == ss0.str()) {
                for (size_t offset = 0; offset <= m_memory->vectorsCount() - step; offset += step) {
                    std::stringstream ss;
                    ss << ss0.str();
                    ss << ": " << offset * sizeof(double_v) / GB << " - "
                        << (offset + size / 16) * sizeof(double_v) / GB;
                    const int repetitions = step / size;
                    Benchmark bench(ss.str().c_str(), size * repetitions * m_threadCount * interpretFactor, unit);
                    Timer timer;
                    for (int rep = 0; rep < 2; ++rep) {
                        size_t offset2 = offset;
                        m_threadPool.executeWith(m_memory, [&offset2, size] { return offset2 += size / 16; }, size / 16, repetitions);
                        testFun({m_memory, &timer, offset, size / 16, repetitions});
                        m_threadPool.waitReady();
                        m_threadPool.eachTimer([&bench](const Timer &t) { bench.addTiming(t); });
                        bench.addTiming(timer);
                    }
                    bench.Print();
                }
            }
        }
    }
}/*}}}*/
inline std::ostream &operator<<(std::ostream &out, const CpuRange &range)/*{{{*/
{
    return out << "Range " << range.first << " - " << range.last << ", step " << range.step << '\n';
}
/*}}}*/
std::vector<CpuRange> parseOnlyCpus()/*{{{*/
{
    enum State {
        ReadFirst, ReadLast, ReadStep
    };
    const std::string cpusStrings = valueForArgument("--cores", std::string());
    std::vector<CpuRange> r;
    CpuRange range = { 0, 0, 1 };
    State state = ReadFirst;
    for (const auto &c : cpusStrings) {
        if (c >= '0' && c <= '9') {
            switch (state) {
            case ReadFirst:
                range.first = range.first * 10 + (c - '0');
                break;
            case ReadLast:
                range.last = range.last * 10 + (c - '0');
                break;
            case ReadStep:
                range.step = range.step * 10 + (c - '0');
                break;
            }
        } else if (c == '-') {
            state = ReadLast;
        } else if (c == ':') {
            state = ReadStep;
            range.step = 0;
        } else if (c == ',') {
            state = ReadFirst;
            r.push_back(range);
            range = { 0, 0, 1 };
        }
    }
    r.push_back(range);
    return r;
}
/*}}}*/
BenchmarkRunner::BenchmarkRunner()/*{{{*/
    : m_coreIds(parseOnlyCpus()),
    m_threadCount(1),
    m_maxMemorySize(largestMemorySize() / GB),
    m_memorySize(valueForArgument("--size", m_maxMemorySize)),
    m_only(valueForArgument("--only", std::string())),
    m_memory(nullptr)
{
#ifndef NO_LIBNUMA/*{{{*/
    if (numa_available() == -1) {
        std::cerr << "NUMA interface does not work. Abort." << std::endl;
        return;
    }
#ifdef LINK_STATICALLY
    numa_init();
#endif

    // first make sure we don't get interleaved memory; this would defeat the purpose of this
    // benchmark
    //numa_set_interleave_mask(numa_no_nodes);
#endif/*}}}*/

    if (m_memorySize < 1) {/*{{{*/
        std::cerr << "Need at least 1GB." << std::endl;
        return;
    }
    if (m_memorySize > m_maxMemorySize) {
        std::cerr << "Not enough memory available. Expect crashes/OOM kills." << std::endl;
    }/*}}}*/
    m_memory = new MemT(m_memorySize * GB / sizeof(double));
    mlockall(MCL_CURRENT);

    if (!m_coreIds.empty()) {/*{{{*/
        Benchmark::addColumn("CPU_ID");
        for (auto cpuRange : m_coreIds) {
            cpu_set_t cpumask;
            sched_getaffinity(0, sizeof(cpu_set_t), &cpumask);
            const int cpuid = cpuRange.first;
            cpuRange.first += cpuRange.step;
            std::ostringstream str;
            str << cpuid;
            m_threadCount = 1;
            for (int id = cpuRange.first; id <= cpuRange.last; id += cpuRange.step) {
                str << ',' << id;
                ++m_threadCount;
            }
            Benchmark::setColumnData("CPU_ID", str.str());
            cpuZero(&cpumask);
            cpuSet(cpuid, &cpumask);
            sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);
            m_threadPool.setPinning(cpuRange);
            executeAllTests();
        }
        return;
    }
/*}}}*/
#ifdef NO_LIBNUMA/*{{{*/
    cpu_set_t cpumask;
    sched_getaffinity(0, sizeof(cpu_set_t), &cpumask);
    int cpucount = cpuCount(&cpumask);
    Benchmark::addColumn("CPU_ID");
    for (int cpuid = valueForArgument("--firstCpu", 1); cpuid < cpucount; cpuid += valueForArgument("--cpuStep", 6)) {
        if (m_threadCount == 1) {
            std::ostringstream str;
            str << cpuid;
            Benchmark::setColumnData("CPU_ID", str.str());
            cpuZero(&cpumask);
            cpuSet(cpuid, &cpumask);
            sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);
        } else {
            std::ostringstream str;
            str << cpuid << " - " << cpuid + m_threadCount - 1;
            Benchmark::setColumnData("CPU_ID", str.str());
            cpuZero(&cpumask);
            cpuSet(cpuid, &cpumask);
            sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);
            m_threadPool.setPinning(cpuid);
        }
        executeAllTests();
    }/*}}}*/
#else/*{{{*/
    // libnuma defines:
    // node: an area where all memory as the same speed as seen from a particular CPU. A node
    //       can contain multiple CPUs
    // cpu: a hardware thread?
    int nodeCount = numa_max_node();
    struct bitmask *nodemask = 0;
    if (nodeCount < 0) {
        std::cerr << "libnuma does not report any NUMA nodes\n";
        nodeCount = 0;
    } else {
        nodemask = numa_allocate_nodemask();
    }
    Benchmark::addColumn("NUMA_ID");
    for (int numaId = valueForArgument("--firstNode", 0); numaId <= nodeCount; numaId += valueForArgument("--nodeStep", 1)) {
        std::ostringstream str;
        str << numaId;
        Benchmark::setColumnData("NUMA_ID", str.str());

        if (nodemask) {
            numa_bitmask_clearall(nodemask);
            numa_bitmask_setbit(nodemask, numaId);
            numa_bind(nodemask);
        }
        executeAllTests();
    }
#endif/*}}}*/
#ifndef NO_LIBNUMA/*{{{*/
    if (nodemask) {
        numa_free_nodemask(nodemask);
    }
#endif/*}}}*/
}/*}}}*/
void BenchmarkRunner::executeAllTests()/*{{{*/
{
    executeTest("bzero"            , &testBzero         , "Byte", 1.);
    executeTest("read"             , &testRead          , "Byte", 1.);
    executeTest("read w/ prefetch" , &testReadPrefetch  , "Byte", 1.);
    executeTest("add 1"            , &testAddOne        , "Byte", 1.);
    executeTest("add 1 w/ prefetch", &testAddOnePrefetch, "Byte", 1.);
    executeTest("read latency"     , &testReadLatency   , "read", 1. / VectorsInCacheLine); // this test only reads one vector out of a cacheline
    Benchmark::finalize();
}
/*}}}*/
int bmain()/*{{{*/
{
    BenchmarkRunner runner;
    return 0;
}/*}}}*/

// vim: foldmethod=marker
