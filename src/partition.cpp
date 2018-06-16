#include "cpucounters.h"
#include "sort/crack.hpp"
#include "sort/partition.hpp"
#include "sort/partition_parallel.hpp"
#include "sort/sort.hpp"
#include "sort/tuple_type.hpp"
#include "util/abrt.hpp"
#include "util/assert.hpp"
#include "util/cpuaffinity.hpp"
#include "util/cpuset.hpp"
#include "util/distribution.hpp"
#include "util/int.hpp"
#include "util/memcpy_fast.hpp"
#include "util/monitor.hpp"
#include "util/rewiring.hpp"
#include "util/system.hpp"
#include "util/typename.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <papi.h>
#include <sstream>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>

#ifdef LINUX
#include <numa.h>
#include <signal.h>
#include <unistd.h>
#endif


using namespace std::chrono;
using clk = high_resolution_clock;
constexpr std::size_t VM_CHUNK_SIZE = 2lu * 1024lu * 1024lu; // 2 MiB
#ifndef NDEBUG
constexpr std::size_t SIZE = 32lu * 1024 * 1024; // Bytes
constexpr std::size_t NUM_CRACKS = 50;
constexpr std::size_t NUM_RUNS = 2;
constexpr std::size_t NUM_STEPS = 4;
constexpr std::initializer_list<double> SELECTIVITIES = {.5, .25, .75, .1, .9};
#else
constexpr std::size_t SIZE = 4lu * 1024 * 1024 * 1024; // Bytes
constexpr std::size_t NUM_CRACKS = 1000;
constexpr std::size_t NUM_RUNS = 5;
constexpr std::size_t NUM_STEPS = 64;
constexpr std::initializer_list<double> SELECTIVITIES = {.50, .40, .60, .33, .67, .25, .75, .10, .90, .01, .99, .001, .999};
#endif
long long shared;


/*-- PAPI events & printing ------------------------------------------------------------------------------------------*/
const EventSet EVENTS{
    PAPI_TOT_CYC,   /* total clock cycles */
    PAPI_TOT_INS,   /* instructions retired */
    PAPI_BR_CN,     /* conditional branches */
    PAPI_BR_MSP,    /* conditional branch mispredictions */
#ifdef PAPI_EXTENDED
    PAPI_RES_STL,   /* cycles stalled on any resource */
    PAPI_PRF_DM,    /* data prefetch cache misses */
    //0x40000003,     /* retired memory loads */
    //0x40000004,     /* retired memory stores */
#endif
};

void print_title(FILE *file)
{
    fprintf(file, "%-6s, %-6s, %-10s, %-15s, %-25s, %-8s, %-12s, %-20s, %-12s, %-35s, %s\n",
                  "Num", "Layout", "Size", "Experiment", "Algorithm", "Threads", "Element", "Distribution", "Selectivity", "Attribute", "Value");
}

#define FMT(X) "%-6lu, %-6s, %10lu, %-15s, %-25s, %8lu, %12lu, %-20s, %-12.3lf, %-35s, %12" X "\n"

template<typename T>
void print_measurement(FILE *file, const Monitor &m, const char *layout, const unsigned long size,
                       const char *experiment, const char *algorithm, const char *distribution,
                       const double selectivity,
                       UncoreCounterState &ucState_before, UncoreCounterState &ucState_after,
                       const std::size_t num_threads)
{
    using std::to_string;
    PCM *pcm = PCM::getInstance();
    const uint64_t nominal_frequency = pcm->getNominalFrequency();

    const double time = double(m.getTSC()) / nominal_frequency;

#define PRINT(FS, NAME, VALUE) \
    fprintf(file, FMT(FS), 1lu, layout, size, experiment, algorithm, num_threads, sizeof(T), distribution, selectivity, (NAME), (VALUE))

    PRINT("lf",     "time",         time * 1000);
    PRINT("lf",     "freq",         m[PAPI_TOT_CYC] / time / 1e9);
    PRINT("lf",     "ipc",          double(m[PAPI_TOT_INS]) / m[PAPI_TOT_CYC]);
    PRINT("lld",    "branch",       m[PAPI_BR_CN]);
    PRINT("lld",    "brmiss",       m[PAPI_BR_MSP]);
#ifdef PAPI_EXTENDED
    PRINT("lld",    "stalled",      m[PAPI_RES_STL]);
    PRINT("lld",    "prefetch",     m[PAPI_PRF_DM]);
    PRINT("lld",    "mem_read",     getBytesReadFromMC(ucState_before, ucState_after));
    PRINT("lld",    "mem_write",    getBytesWrittenToMC(ucState_before, ucState_after));
#endif

#undef PRINT
}

template<typename T>
void print_cracking_per_query(FILE *file, const char *layout, const unsigned long size, const char *algorithm,
                              const char *distribution, const double selectivity,
                              std::vector<uint64_t> &tsc,
                              const std::size_t num_threads)
{
    using std::to_string;
    PCM *pcm = PCM::getInstance();
    const uint64_t nominal_frequency = pcm->getNominalFrequency();

    for (std::size_t i = 0; i != tsc.size(); ++i) {
        const double time = double(tsc[i]) / nominal_frequency * 1000;
        fprintf(file, FMT("lf"), (i % NUM_CRACKS) + 1, layout, size, "per_query", algorithm, num_threads, sizeof(T), distribution, selectivity, "time", time);
    }
}

template<typename T>
void print_time_of_merge(FILE *file, const char *layout, const unsigned long size, const char *algorithm,
                         const char *distribution, const double selectivity,
                         const std::size_t num_threads,
                         const uint64_t tsc)
{
    using std::to_string;
    PCM *pcm = PCM::getInstance();
    const uint64_t nominal_frequency = pcm->getNominalFrequency();

    const double time = double(tsc) / nominal_frequency * 1000;
    fprintf(file, FMT("lf"), 1lu, layout, size, "merge", algorithm, num_threads, sizeof(T), distribution, selectivity, "time", time);
}

#undef FMT


/*-- Benchmarks ------------------------------------------------------------------------------------------------------*/
/*-- Macros ----------------------------------------------------------------------------------------------------------*/
#define RUN(LAYOUT, EXP, ALG, NUM_THREADS, CODE) {\
    Monitor M{EVENTS}; \
    uint64_t sum_TSC = 0; \
    printf("  ` %s (%.2lf MiB): %s (%.3lf) / %-25s (%lu threads)", LAYOUT, 2. * num_tuples * sizeof(T) / (1024 * 1024), EXP, selectivity, ALG, std::size_t(NUM_THREADS)); \
    fflush(stdout); \
    const unsigned long long byte_sum_before = compute_sum((uint8_t*) data, 2 * num_tuples * sizeof(T)); \
    PCM *pcm = PCM::getInstance(); \
    for (std::size_t i = 0; i != NUM_RUNS; ++i) { \
        memfile.resize(num_tuples * 2 * sizeof(T)); \
        auto vm_copy = memfile.map<VM_CHUNK_SIZE>(num_tuples * 2 * sizeof(T), 0); \
        copy = reinterpret_cast<decltype(copy)>(vm_copy.addr); \
        parallel_memcpy(memcpy_fast, copy, data, num_tuples * 2 * sizeof(T), 2); \
        pcm->programServerUncoreMemoryMetrics(); \
        auto ucState_before = pcm->getSocketCounterState(0); \
        M.start(); \
        { CODE } \
        M.stop(); \
        auto ucState_after = pcm->getSocketCounterState(0); \
        print_measurement<T>(file, M, LAYOUT, 2 * num_tuples * sizeof(T), EXP, ALG, distribution, selectivity, ucState_before, ucState_after, NUM_THREADS); \
        printf(" %9.3lf", (double(M.getTSC()) / pcm->getNominalFrequency() * 1000));\
        sum_TSC += M.getTSC(); \
        fflush(stdout); \
        const unsigned long long byte_sum_after = compute_sum((uint8_t*) copy, 2 * num_tuples * sizeof(T)); \
        if (byte_sum_after != byte_sum_before) { \
            std::cerr << "expected " << byte_sum_before << ", got " << byte_sum_after << '\n'; \
            abrt("sums differ"); \
        } \
        vm_copy.unmap(); \
    } \
    printf("    Mean: %9.3lf\n", double(sum_TSC) / NUM_RUNS / pcm->getNominalFrequency() * 1000); \
}


void parallel_memcpy(std::function<void*(void*, const void*, std::size_t)> f,
                     void *dst, const void *src, std::size_t size, std::size_t num_threads)
{
    uint8_t *from = (uint8_t*) src;
    uint8_t *to = (uint8_t*) dst;
    const std::size_t chunk_size = size / num_threads;

    auto g = [=](void *dst, void *src, std::size_t size) { f(dst, src, size); };

    std::thread *threads = new std::thread[num_threads];
    for (unsigned t = 0; t != num_threads; ++t)
        new (&threads[t]) std::thread(g, to + t * chunk_size, from + t * chunk_size, chunk_size);

    for (unsigned t = 0; t != num_threads; ++t)
        threads[t].join();

    delete[] threads;
}

template<typename T>
void run_benchmarks_misc(FILE *file, const unsigned max_num_threads, const std::size_t num_tuples, const T *data)
{
    std::cout << "` Miscellaneous benchmarks.\n";

    const char *distribution = "NA";
    const double selectivity = NAN;
    rewiring::memory_file<false> memfile{"smallpages_misc"};
    T *copy;

    for (std::size_t threads = 1; threads <= max_num_threads; ++threads) {
        RUN("NA", "misc", "memcpy", threads, {
            parallel_memcpy(std::memcpy, copy, data, num_tuples * 2 * sizeof(T), threads);
        });
    }

    for (std::size_t threads : {1, 2, 4}) {
        RUN("NA", "misc", "memcpy (fast)", threads, {
            parallel_memcpy(memcpy_fast, copy, data, num_tuples * 2 * sizeof(T), threads);
        });
    }
}

#define PRINT_MERGE(LAYOUT, ALGORITHM, NUM_THREADS, TSC) \
    print_time_of_merge<T>(file, (LAYOUT), 2 * num_tuples * sizeof(T), (ALGORITHM), distribution, selectivity, (NUM_THREADS), (TSC))

template<typename T>
void run_partitioning_AoS(FILE *file, const unsigned max_num_threads, const char *distribution,
                          const std::size_t num_tuples, const T *sorted, const double selectivity, const T *data)
{
    const tuple_type<T> pivot = sorted[std::size_t(selectivity * num_tuples)];
    tuple_type<T> *copy;

    /* Smallpages */
    {
        rewiring::memory_file<false> memfile{"smallpages_partition_AoS"};

        RUN("AoS", "partition", "branching", 1, {
            partition_branching p;
            tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
        });

#if 0
        RUN("AoS", "partition", "naive", 1, {
            partition_predicated_naive p;
            tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
        });
#endif

        RUN("AoS", "partition", "pirk", 1, {
            partition_predicated_pirk p;
            tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
        });

#if 0
        RUN("AoS", "partition", "decoupled", 1, {
            partition_predicated_decoupled p;
            tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
        });
#endif

        RUN("AoS", "partition", "register", 1, {
            partition_predicated_register p;
            tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
        });

        RUN("AoS", "partition", "rewired (smallpages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired<false> p{memfile};
            tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
            vm_buffer.unmap();
        });

#if __AVX2__
        RUN("AoS", "partition", "rewired simd (smallpages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired_simd<false> p{memfile};
            tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
            vm_buffer.unmap();
        });
#endif

        for (std::size_t threads = 2; threads <= max_num_threads; ++threads) {
            uint64_t tsc_merge;
            RUN("AoS", "partition", "pirk (MT)", threads, {
                tuple_type<T> *mid = refined_partition_and_merge_AoS(pivot, copy, copy + num_tuples, threads, tsc_merge);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
            });
            PRINT_MERGE("AoS", "pirk (MT)", threads, tsc_merge);
        }
    }

    /* Hugepages */
    {
        rewiring::memory_file<true> memfile{"hugepages_partition_AoS"};

        RUN("AoS", "partition", "rewired (hugepages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired<true> p{memfile};
            tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
            vm_buffer.unmap();
        });

#if __AVX2__
        RUN("AoS", "partition", "rewired simd (hugepages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired_simd<true> p{memfile};
            tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
            vm_buffer.unmap();
        });
#endif

        for (std::size_t threads = 2; threads <= max_num_threads; ++threads) {
            uint64_t tsc_merge;
            RUN("AoS", "partition", "rewired (hugepages+MT)", threads, {
                tuple_type<T> *mid = rewired_partition_and_merge_AoS(pivot, copy, copy + num_tuples, memfile, vm_copy, threads, tsc_merge);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
            });
            PRINT_MERGE("AoS", "rewired (hugepages+MT)", threads, tsc_merge);
        }

        for (std::size_t threads = 2; threads <= max_num_threads; ++threads) {
            uint64_t tsc_merge;
            RUN("AoS", "partition", "hybrid (hugepages+MT)", threads, {
                tuple_type<T> *mid = hybrid_partition_and_merge_AoS(pivot, copy, copy + num_tuples, vm_copy, threads, tsc_merge);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
            });
            PRINT_MERGE("AoS", "hybrid (hugepages+MT)", threads, tsc_merge);
        }
    }
}

/**
 * Crack a subsection of the data using the given partitioning kernel and cracker index.
 */
template<typename P, typename T>
void crack(P kernel, crackerindex_t<T> *index, tuple_type<T> *begin, std::size_t size, const T *cracks)
{
    Timer timer;
    for (std::size_t i = 0; i != NUM_CRACKS; ++i)
        crack_AoS<P, T>(*index, begin, size, cracks[i], kernel, timer);
}
template<typename P, typename T, bool UseHugepages>
void crack_rewired(P kernel, crackerindex_t<T> *index, const std::size_t num_tuples,
                   rewiring::memory_mapping<UseHugepages, VM_CHUNK_SIZE> *data,
                   rewiring::memory_mapping<UseHugepages, VM_CHUNK_SIZE> *buffer,
                   const T *cracks)
{
#if 0
    std::ostringstream oss;
    oss << "Thread starts at " << data->addr << " and processes " << num_tuples << " tuples";
    std::cerr << oss.str() << std::endl;
#endif

    Timer timer;
    for (std::size_t i = 0; i != NUM_CRACKS; ++i)
        crack_AoS(*index, *data, num_tuples, cracks[i], kernel, *buffer, timer);
}

/**
 * Divide the data into `num_threads` equi-sized chunks and crack them independently and in parallel.
 */
template<typename P, typename T>
void crack_parallel(P kernel, std::size_t num_threads, std::size_t num_tuples, tuple_type<T> *data, const T *cracks)
{
    std::thread *threads = new std::thread[num_threads];
    crackerindex_t<T> *CIs = new crackerindex_t<T>[num_threads];
    const std::size_t step = num_tuples / num_threads;

    auto start = data;
    for (std::size_t i = 0; i != num_threads - 1; ++i, start += step) {
        assert(start + step < data + num_tuples);
        new (&threads[i]) std::thread(crack<P,T>, kernel, &CIs[i], start, step, cracks);
    }
    new (&threads[num_threads - 1]) std::thread(crack<P,T>,
                                                kernel,
                                                &CIs[num_threads - 1],
                                                start,
                                                data + num_tuples - start,
                                                cracks);


    for (unsigned t = 0; t != num_threads; ++t)
        threads[t].join();

    delete[] CIs;
    delete[] threads;
}
template<typename P, typename T, bool UseHugepages>
void crack_parallel(P kernel, std::size_t num_threads, std::size_t num_tuples,
                    rewiring::memory_file<UseHugepages> &memfile,
                    rewiring::memory_mapping<UseHugepages, VM_CHUNK_SIZE> &data,
                    const T *cracks)
{
    std::thread *threads = new std::thread[num_threads];
    crackerindex_t<T> *CIs = new crackerindex_t<T>[num_threads];
    rewiring::memory_mapping<UseHugepages, VM_CHUNK_SIZE> *slices = new rewiring::memory_mapping<UseHugepages, VM_CHUNK_SIZE>[num_threads];
    rewiring::memory_mapping<UseHugepages, VM_CHUNK_SIZE> *buffers = new rewiring::memory_mapping<UseHugepages, VM_CHUNK_SIZE>[num_threads];

#ifdef VERBOSE
    std::cout << std::endl;
    std::cerr << "Crack parallel " << num_tuples << " tuples at " << data.addr << " with " << num_threads << " threads\n";
#endif

    constexpr std::size_t NUM_TUPLES_PER_CHUNK = VM_CHUNK_SIZE / sizeof(tuple_type<T>);

    const std::size_t num_chunks = data.size / VM_CHUNK_SIZE;
    const std::size_t chunks_per_thread = (num_chunks + num_threads - 1) / num_threads; // CEIL(num_chunks / num_threads)
#ifdef VEROSE
    std::cerr << "we have " << num_chunks << " chunks, and use " << chunks_per_thread << " chunks per thread\n";
#endif

    std::size_t chunks = 0;
    std::size_t thread = 0;
    for (; thread != num_threads and chunks + chunks_per_thread <= num_chunks; ++thread, chunks += chunks_per_thread) {
#ifdef VERBOSE
        std::cerr << "Thread " << thread << ":\n";
#endif
        buffers[thread] = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);

        /* Create the data slice for this thread */
        auto &vm = slices[thread];
        vm.fd = data.fd;
        vm.size = sizeof(tuple_type<T>) * chunks_per_thread * NUM_TUPLES_PER_CHUNK;
        vm.addr = (tuple_type<T>*) data.addr + thread * chunks_per_thread * NUM_TUPLES_PER_CHUNK;

#ifdef VERBOSE
        std::cerr << "  ` addr " << vm.addr << '\n';
        std::cerr << "  ` size " << vm.size << " Bytes\n";
        std::cerr << "  ` #chunks: " << chunks_per_thread << '\n';
#endif

        /* assign the chunks to this thread's buffer */
        vm.chunks = new long[chunks_per_thread];
        for (std::size_t i = 0; i != chunks_per_thread; ++i) {
            vm.chunks[i] = data[i + chunks];
#ifdef VERBOSE
            std::cerr << "    ` assign chunk " << (i + chunks) << " to this thread\n";
#endif
        }

        /* start the cracking thread */
        new (&threads[thread]) std::thread(crack_rewired<P,T, UseHugepages>,
                                           kernel,
                                           &CIs[thread],
                                           chunks_per_thread * NUM_TUPLES_PER_CHUNK,
                                           &slices[thread],
                                           &buffers[thread],
                                           cracks);
    }

    assert(chunks <= num_chunks);
    if (chunks < num_chunks) {
#ifdef VERBOSE
        std::cerr << "Thread " << thread << ":\n";
#endif
        const std::size_t n = num_chunks - chunks;
        buffers[thread] = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);

        /* Create the data slice for this thread */
        auto &vm = slices[thread];
        vm.fd = data.fd;
        vm.size = sizeof(tuple_type<T>) * n * NUM_TUPLES_PER_CHUNK;
        vm.addr = (tuple_type<T>*) data.addr + thread * chunks_per_thread * NUM_TUPLES_PER_CHUNK;

#ifdef VERBOSE
        std::cerr << "  ` addr " << vm.addr << '\n';
        std::cerr << "  ` size " << vm.size << " Bytes\n";
        std::cerr << "  ` #chunks: " << n << '\n';
#endif

        /* assign the chunks to this thread's buffer */
        vm.chunks = new long[n];
        for (std::size_t i = 0; i != n; ++i) {
            vm.chunks[i] = data[i + chunks];
#ifdef VERBOSE
            std::cerr << "    ` assign chunk " << (i + chunks) << " to this thread\n";
#endif
        }

        new (&threads[thread]) std::thread(crack_rewired<P, T, UseHugepages>,
                                           kernel,
                                           &CIs[thread],
                                           n * NUM_TUPLES_PER_CHUNK,
                                           &slices[thread],
                                           &buffers[thread],
                                           cracks);
        ++thread;
    }

#ifdef VERBOSE
    std::cerr << "spawned " << thread << " threads.  Witing for threads to finish..." << std::endl;
#endif

    for (std::size_t i = 0; i != thread; ++i)
        threads[i].join();

#ifdef VERBOSE
    std::cerr << std::endl;
#endif

    delete[] CIs;
    delete[] threads;
    delete[] slices;
}

/**
 * Evaluate how well partitioning algorithms scale up with the number of threads, ignoring the fact that data must be
 * merged afterwards.
 */
template<typename T>
void run_parallel_cracking(FILE *file, const unsigned max_num_threads, const char *distribution,
                           const std::size_t num_tuples, const T *cracks, const T *data)
{
    const double selectivity = NAN;
    rewiring::memory_file<true> memfile{"hugepages_crack_parallel_AoS"};
    tuple_type<T> *copy;

    for (std::size_t threads = 1; threads <= max_num_threads; ++threads) {
#if 1
        RUN("AoS", "scaleup", "branching", threads, {
            partition_branching p;
            crack_parallel(p, threads, num_tuples, copy, cracks);
        });
        RUN("AoS", "scaleup", "pirk", threads, {
            partition_predicated_pirk p;
            crack_parallel(p, threads, num_tuples, copy, cracks);
        });
        RUN("AoS", "scaleup", "register", threads, {
            partition_predicated_register p;
            crack_parallel(p, threads, num_tuples, copy, cracks);
        });
#endif

        RUN("AoS", "scaleup", "rewired (hugepages)", threads, {
            partition_rewired p{memfile};
            crack_parallel(p, threads, num_tuples, memfile, vm_copy, cracks);
        });

#if 1
        RUN("AoS", "scaleup", "rewired simd (hugepages)", threads, {
            partition_rewired_simd p{memfile};
            crack_parallel(p, threads, num_tuples, memfile, vm_copy, cracks);
        });
#endif
    }
}

template<typename T>
void run_partitioning_SoA(FILE *file, const unsigned max_num_threads, const char *distribution,
                          const std::size_t num_tuples, const T *sorted, const double selectivity, const T *data)
{
    const T pivot = sorted[std::size_t(selectivity * num_tuples)];
    T *copy;
    (void) max_num_threads; /* TODO */

    /* Smallpages */
    {
        rewiring::memory_file<false> memfile{"smallpages_partition_SoA"};

        RUN("SoA", "partition", "branching", 1, {
            partition_branching p;
            T *mid = p.SoA(pivot, copy, copy + num_tuples, num_tuples);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
        });

#if 0
        RUN("SoA", "partition", "naive", 1, {
            partition_predicated_naive p;
            T *mid = p.SoA(pivot, copy, copy + num_tuples, num_tuples);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
        });
#endif

        RUN("SoA", "partition", "pirk", 1, {
            partition_predicated_pirk p;
            T *mid = p.SoA(pivot, copy, copy + num_tuples, num_tuples);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
        });

#if 0
        RUN("SoA", "partition", "decoupled", 1, {
            partition_predicated_decoupled p;
            T *mid = p.SoA(pivot, copy, copy + num_tuples, num_tuples);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
        });
#endif

        RUN("SoA", "partition", "register", 1, {
            partition_predicated_register p;
            T *mid = p.SoA(pivot, copy, copy + num_tuples, num_tuples);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
        });

        RUN("SoA", "partition", "rewired (smallpages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired<false> p{memfile};
            T *mid = p.SoA(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
            vm_buffer.unmap();
        });

#if __AVX2__
        RUN("SoA", "partition", "rewired simd (smallpages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired_simd<false> p{memfile};
            T *mid = p.SoA(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
            vm_buffer.unmap();
        });
#endif
    }

    /* Hugepages */
    {
        rewiring::memory_file<true> memfile{"hugepages_partition_SoA"};

        RUN("SoA", "partition", "rewired (hugepages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired<true> p{memfile};
            T *mid = p.SoA(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
            vm_buffer.unmap();
        });

#if __AVX2__
        RUN("SoA", "partition", "rewired simd (hugepages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired_simd<true> p{memfile};
            T *mid = p.SoA(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
            (void) mid;
            assert(verify_partition(copy, mid, copy + num_tuples));
            vm_buffer.unmap();
        });
#endif
    }
}

#undef PRINT_MERGE

#define PRINT_PER_QUERY(LAYOUT, ALGORITHM, TSC, NUM_THREADS) \
    print_cracking_per_query<T>(file, LAYOUT, 2 * num_tuples * sizeof(T), ALGORITHM, distribution, .5, TSC, NUM_THREADS)
#define PRINT_TIMER() \
    printf("      Average time for cracking w/o index lookup: %8.3lf ms\n", std::chrono::duration_cast<std::chrono::nanoseconds>(timer.reset()).count() / 1e6 / NUM_RUNS);

template<typename T>
void run_cracking_AoS(FILE *file, const char *distribution, const std::size_t num_tuples, const T *cracks,
                      const T *data)
{
    const double selectivity = NAN;
    tuple_type<T> *copy;
    std::vector<uint64_t> tsc{NUM_CRACKS};
    Timer timer;

    /* Smallpages */
    {
        rewiring::memory_file<false> memfile{"smallpages_cracking_AoS"};

        tsc.clear();
        RUN("AoS", "crack", "branching", 1, {
            partition_branching p;
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtsc();
                crack_AoS(index, copy, num_tuples, cracks[i], p, timer);
                const uint64_t tsc_stop = rdtscp();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
        });
        PRINT_PER_QUERY("AoS", "branching", tsc, 1);
        PRINT_TIMER();

        tsc.clear();
        RUN("AoS", "crack", "pirk", 1, {
            partition_predicated_pirk p;
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_AoS(index, copy, num_tuples, cracks[i], p, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
        });
        PRINT_PER_QUERY("AoS", "pirk", tsc, 1);
        PRINT_TIMER();

        tsc.clear();
        RUN("AoS", "crack", "register", 1, {
            partition_predicated_register p;
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_AoS(index, copy, num_tuples, cracks[i], p, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
        });
        PRINT_PER_QUERY("AoS", "register", tsc, 1);
        PRINT_TIMER();

#if 1
        tsc.clear();
        RUN("AoS", "crack", "rewired (smallpages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired<false> p{memfile};
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_AoS(index, vm_copy, num_tuples, cracks[i], p, vm_buffer, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
            vm_buffer.unmap();
        });
        PRINT_PER_QUERY("AoS", "rewired (smallpages)", tsc, 1);
        PRINT_TIMER();
#endif

#if __AVX2__
        tsc.clear();
        RUN("AoS", "crack", "rewired simd (smallpages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired_simd<false> p{memfile};
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_AoS(index, vm_copy, num_tuples, cracks[i], p, vm_buffer, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
            vm_buffer.unmap();
        });
        PRINT_PER_QUERY("AoS", "rewired simd (smallpages)", tsc, 1);
        PRINT_TIMER();
#endif
    }

    {
        rewiring::memory_file<true> memfile{"hugepages_cracking_AoS"};

        tsc.clear();
        RUN("AoS", "crack", "rewired (hugepages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired<true> p{memfile};
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_AoS(index, vm_copy, num_tuples, cracks[i], p, vm_buffer, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
            vm_buffer.unmap();
        });
        PRINT_PER_QUERY("AoS", "rewired (hugepages)", tsc, 1);
        PRINT_TIMER();

#if __AVX2__
        tsc.clear();
        RUN("AoS", "crack", "rewired simd (hugepages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired_simd<true> p{memfile};
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_AoS(index, vm_copy, num_tuples, cracks[i], p, vm_buffer, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
            vm_buffer.unmap();
        });
        PRINT_PER_QUERY("AoS", "rewired simd (hugepages)", tsc, 1);
        PRINT_TIMER();
#endif
    }
}

template<typename T>
void run_cracking_SoA(FILE *file, const char *distribution, const std::size_t num_tuples, const T *cracks,
                      const T *data)
{
    const double selectivity = NAN;
    T *copy;
    std::vector<uint64_t> tsc{NUM_CRACKS};
    Timer timer;

    /* Smallpages */
    {
        rewiring::memory_file<false> memfile{"smallpages_cracking_SoA"};

        tsc.clear();
        RUN("SoA", "crack", "branching", 1, {
            partition_branching p;
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_SoA(index, copy, num_tuples, cracks[i], p, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
        });
        PRINT_PER_QUERY("SoA", "branching", tsc, 1);
        PRINT_TIMER();

        tsc.clear();
        RUN("SoA", "crack", "pirk", 1, {
            partition_predicated_pirk p;
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_SoA(index, copy, num_tuples, cracks[i], p, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
        });
        PRINT_PER_QUERY("SoA", "pirk", tsc, 1);
        PRINT_TIMER();

        tsc.clear();
        RUN("SoA", "crack", "register", 1, {
            partition_predicated_register p;
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_SoA(index, copy, num_tuples, cracks[i], p, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
        });
        PRINT_PER_QUERY("SoA", "register", tsc, 1);
        PRINT_TIMER();

#if 1
        tsc.clear();
        RUN("SoA", "crack", "rewired (smallpages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired<false> p{memfile};
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_SoA(index, vm_copy, num_tuples, cracks[i], p, vm_buffer, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
            vm_buffer.unmap();
        });
        PRINT_PER_QUERY("SoA", "rewired (smallpages)", tsc, 1);
        PRINT_TIMER();
#endif

#if __AVX2__
        tsc.clear();
        RUN("SoA", "crack", "rewired simd (smallpages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired_simd<false> p{memfile};
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_SoA(index, vm_copy, num_tuples, cracks[i], p, vm_buffer, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
            vm_buffer.unmap();
        });
        PRINT_PER_QUERY("SoA", "rewired simd (smallpages)", tsc, 1);
        PRINT_TIMER();
#endif
    }

    {
        rewiring::memory_file<true> memfile{"hugepages_cracking_SoA"};

        tsc.clear();
        RUN("SoA", "crack", "rewired (hugepages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired<true> p{memfile};
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_SoA(index, vm_copy, num_tuples, cracks[i], p, vm_buffer, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
            vm_buffer.unmap();
        });
        PRINT_PER_QUERY("SoA", "rewired (hugepages)", tsc, 1);
        PRINT_TIMER();

#if __AVX2__
        tsc.clear();
        RUN("SoA", "crack", "rewired simd (hugepages)", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired_simd<true> p{memfile};
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                crack_SoA(index, vm_copy, num_tuples, cracks[i], p, vm_buffer, timer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
            vm_buffer.unmap();
        });
        PRINT_PER_QUERY("SoA", "rewired simd (hugepages)", tsc, 1);
        PRINT_TIMER();

        if constexpr (sizeof(T) == 4) {
        tsc.clear();
        RUN("SoA", "crack", "mixed cracking", 1, {
            auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
            partition_rewired_simd<true> p{memfile};
            crackerindex_t<T> index;
            for (std::size_t i = 0; i != NUM_CRACKS; ++i) {
                const uint64_t tsc_start = rdtscp();
                mixed_crack_SoA(index, vm_copy, num_tuples, cracks[i], p, vm_buffer);
                const uint64_t tsc_stop = rdtsc();
                tsc.push_back(tsc_stop - tsc_start);
            }
            assert(index.size() == NUM_CRACKS);
            vm_buffer.unmap();
        });
        PRINT_PER_QUERY("SoA", "mixed cracking", tsc, 1);
        }
#endif
    }
}

template<typename T>
void run_increasing_AoS(FILE *file, const char *distribution, const double selectivity, const T *data)
{
    tuple_type<T> *copy;

    rewiring::memory_file<false> smallpages{"smallpages_increasing_AoS"};
    rewiring::memory_file<true> hugepages{"hugepages_increasing_AoS"};
    for (auto size = VM_CHUNK_SIZE; size <= NUM_STEPS * VM_CHUNK_SIZE; size += VM_CHUNK_SIZE) {
        const std::size_t num_tuples = size / 2 / sizeof(T);

        /* Get the exact median of the chunk. */
        tuple_type<T> *sorted = static_cast<tuple_type<T>*>(aligned_alloc(PAGESIZE, size));
        parallel_memcpy(memcpy_fast, sorted, data, size, 2);
        using std::sort;
        sort(sorted, sorted + num_tuples);
        const tuple_type<T> pivot = sorted[std::size_t(selectivity * num_tuples)];
        free(sorted);

        /* Smallpages */
        {
            auto &memfile = smallpages;
            RUN("AoS", "increasing", "register", 1, {
                partition_predicated_register p;
                tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
                assert(mid == copy + std::size_t(selectivity * num_tuples));
            });

            RUN("AoS", "increasing", "rewired (smallpages)", 1, {
                auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);
                partition_rewired<false> p{memfile};
                tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
                vm_buffer.unmap();
            });

#if __AVX2__
            RUN("AoS", "increasing", "rewired simd (smallpages)", 1, {
                auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);
                partition_rewired_simd<false> p{memfile};
                tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
                vm_buffer.unmap();
            });
#endif
        }

        /* Hugepages */
        {
            auto &memfile = hugepages;
            RUN("AoS", "increasing", "rewired (hugepages)", 1, {
                auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);
                partition_rewired<true> p{memfile};
                tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
                vm_buffer.unmap();
            });

#if __AVX2__
            RUN("AoS", "increasing", "rewired (hugepages)", 1, {
                auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(4 * VM_CHUNK_SIZE, memfile.size);
                partition_rewired<true> p{memfile};
                tuple_type<T> *mid = p.AoS(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
                vm_buffer.unmap();
            });
#endif
        }
    }
}

template<typename T>
void run_increasing_SoA(FILE *file, const char *distribution, const double selectivity, const T *data)
{
    T *copy;

    rewiring::memory_file<false> smallpages{"smallpages_increasing_SoA"};
    rewiring::memory_file<true> hugepages{"hugepages_increasing_SoA"};
    for (auto size = 2 * VM_CHUNK_SIZE; size <= NUM_STEPS * VM_CHUNK_SIZE; size += 2 * VM_CHUNK_SIZE) {
        const std::size_t num_tuples = size / 2 / sizeof(T);
        T *sorted = static_cast<T*>(aligned_alloc(PAGESIZE, size / 2));
        parallel_memcpy(memcpy_fast, sorted, data, size / 2, 2);
        using std::sort;
        sort(sorted, sorted + num_tuples);
        const T pivot = sorted[std::size_t(selectivity * num_tuples)];
        free(sorted);

        /* Smallpages */
        {
            auto &memfile = smallpages;
            RUN("SoA", "increasing", "register", 1, {
                partition_predicated_register p;
                T *mid = p.SoA(pivot, copy, copy + num_tuples, num_tuples);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
                assert(mid == copy + std::size_t(selectivity * num_tuples));
            });

            RUN("SoA", "increasing", "rewired (smallpages)", 1, {
                auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
                partition_rewired<false> p{memfile};
                T *mid = p.SoA(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
                vm_buffer.unmap();
            });

#if __AVX2__
            RUN("SoA", "increasing", "rewired simd (smallpages)", 1, {
                auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
                partition_rewired_simd<false> p{memfile};
                T *mid = p.SoA(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
                vm_buffer.unmap();
            });
#endif
        }

        /* Hugepages */
        {
            auto &memfile = hugepages;
            RUN("SoA", "increasing", "rewired (hugepages)", 1, {
                auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
                partition_rewired<true> p{memfile};
                T *mid = p.SoA(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
                vm_buffer.unmap();
            });

#if __AVX2__
            RUN("SoA", "increasing", "rewired simd (hugepages)", 1, {
                auto vm_buffer = memfile.template map<VM_CHUNK_SIZE>(8 * VM_CHUNK_SIZE, memfile.size);
                partition_rewired_simd<true> p{memfile};
                T *mid = p.SoA(pivot, copy, copy + num_tuples, vm_copy, vm_buffer);
                (void) mid;
                assert(verify_partition(copy, mid, copy + num_tuples));
                vm_buffer.unmap();
            });
#endif
        }
    }
}

template<typename T>
void run_benchmarks_AoS(FILE *file, const unsigned max_num_threads, const char *distribution,
                        const std::size_t num_tuples, const T *sorted, const T *cracks, const T *data)
{
    std::cout << "` Use Array of Structs (AoS) layout.\n";

#if 1
    for (double sel : SELECTIVITIES) {
        run_partitioning_AoS(file, max_num_threads, distribution, num_tuples, sorted, sel, data);
        std::cout << '\n';
    }
#endif

#if 0
    run_increasing_AoS(file, distribution, .5, data);
    std::cout << '\n';
#endif

#if 1
    run_cracking_AoS(file, distribution, num_tuples, cracks, data);
    std::cout << '\n';
#endif

#if 1
    run_parallel_cracking(file, max_num_threads, distribution, num_tuples, cracks, data);
#endif
}

template<typename T>
void run_benchmarks_SoA(FILE *file, const unsigned max_num_threads, const char *distribution,
                        const std::size_t num_tuples, const T *sorted, const T *cracks, const T *data)
{
    std::cout << "` Use Struct of Arrays (SoA) layout.\n";

#if 1
    for (double sel : SELECTIVITIES) {
        run_partitioning_SoA(file, max_num_threads, distribution, num_tuples, sorted, sel, data);
        std::cout << '\n';
    }
#endif

#if 0
    run_increasing_SoA(file, distribution, .5, data);
    std::cout << '\n';
#endif

#if 1
    run_cracking_SoA(file, distribution, num_tuples, cracks, data);
    std::cout << '\n';
#endif
}

/**
 * Set up the benchmark environment and invoke benchmarks.
 */
template<typename T>
void run_benchmarks(FILE *file, std::size_t bytes, const unsigned max_num_threads, DISTRIBUTION dist)
{
    const char *distribution_name = DISTRIBUTION_TO_STR.at(dist).c_str();
    const std::size_t num_tuples = (bytes / 2) / sizeof(T);
    bytes = num_tuples * 2 * sizeof(T);

    std::cout << "\nRun benchmarks on data set of " << bytes / (1024 * 1024) << " MiB with " << sizeof(T)
              << " Byte elements and " << distribution_name << " distribution.\n";

    /* Create input data in SoA layout. */
    T *data = static_cast<T*>(aligned_alloc(PAGESIZE, bytes));
    std::cout << "  ` Generating " << num_tuples << " keys: ";
    create_distribution(dist, data, num_tuples, T(-1), /* bits= */ 6, /* max= */ T(T(-1) >> 1));
    std::cout << "  ` Generating " << num_tuples << " payloads: ";
    create_uniform_dense(data + num_tuples, num_tuples, T(-1));

#ifndef NDEBUG
    const long long total_sum = compute_sum(data, 2 * num_tuples);
#endif

    /* Create cracks. */
    T *cracks = static_cast<T*>(aligned_alloc(PAGESIZE, NUM_CRACKS * sizeof(T)));
    std::cout << "  ` Generating " << NUM_CRACKS << " cracks: ";
    switch (dist) {
        case DIST_Ascending:
        case DIST_Uniform_Dense:
            if (std::log2(num_tuples) >= 8 * sizeof(T))
                create_uniform_sparse(cracks, NUM_CRACKS, T(-1));
            else
                create_uniform_sparse(cracks, NUM_CRACKS, T(-1), T(num_tuples));
            break;

        default:
            create_uniform_sparse(cracks, NUM_CRACKS, T(-1));
            break;
    }

    /* Sort keys. */
    std::cout << "  ` Sorting keys.\n";
    T *sorted = static_cast<T*>(aligned_alloc(PAGESIZE, num_tuples * sizeof(T)));
    std::copy_n(data, num_tuples, sorted);
    qsort(sorted, sorted + num_tuples, partition_predicated_register{});
    assert(is_sorted(sorted, sorted + num_tuples));

#if 0
    run_benchmarks_misc(file, max_num_threads, num_tuples, data);
#endif

#if 0
    run_benchmarks_SoA(file, max_num_threads, distribution_name, num_tuples, sorted, cracks, data);
#endif

    /* Convert SoA to AoS layout. */
#ifndef NDEBUG
    auto sum_keys_SoA = compute_sum(data, num_tuples);
#endif
    for (auto runner_keys = data + 1, runner_payloads = data + num_tuples, end = data + num_tuples;
         runner_keys < end; runner_keys += 2, runner_payloads += 2) {
        using std::swap;
        swap(*runner_keys, *runner_payloads);
    }
#ifndef NDEBUG
    decltype(sum_keys_SoA) sum_keys_AoS = 0;
    for (auto p = data; p != data + 2 * num_tuples; p += 2)
        sum_keys_AoS += *p;
    assert(sum_keys_SoA == sum_keys_AoS);
#endif
    assert(total_sum == compute_sum(data, 2 * num_tuples));
    run_benchmarks_AoS(file, max_num_threads, distribution_name, num_tuples, sorted, cracks, data);

    free(data);
    free(sorted);
    free(cracks);
}


/*-- CPUSET handler --------------------------------------------------------------------------------------------------*/
#ifdef LINUX
void cpuset_create(const unsigned num_threads)
{
    const unsigned threads_total = std::thread::hardware_concurrency();
    const bool is_numa_available = numa_available() != -1;
    std::ostringstream oss;
    std::cerr << "Create CPUSETs...\n";

    std::cerr << "` Assign cores 0-" << (num_threads - 1) << " and memory region 0 to cpuset \"hpc\".\n";
    oss << 0 << '-' << (num_threads - 1);
    cpuset_create_partition("hpc", oss.str().c_str(), "0");

    oss.str("");
    oss << num_threads << '-' << (threads_total - 1);
    std::string s_cores = oss.str();
    oss.str("");
    if (is_numa_available) {
        const auto max_node = numa_max_node();
        if (max_node)
            oss << 1 << '-' << max_node;
        else oss << 0;
    } else
        oss << 0;
    std::string s_memory = oss.str();
    std::cerr << "` Assign cores " << s_cores << " and memory regions " << s_memory << " to cpuset \"general\".\n";
    cpuset_create_partition("general", s_cores.c_str(), s_memory.c_str());

    const auto num_general = cpuset_move_all_tasks_to_partition("general");
    std::cerr << "` Moved " << num_general << " taks to \"general\".\n";

    const uint32_t pid = getpid();
    std::cerr << "` Move this task (" << pid << ") to \"hpc\" and set CPU affinity.\n";
    cpuset_move_task_to_partition("hpc", pid);
    std::vector<unsigned> cpuids;
    for (unsigned i = 0; i != num_threads; ++i)
        cpuids.push_back(i);
    set_cpu_affinity(cpuids);
    std::cerr << '\n';
}

void cpuset_remove()
{
    std::cerr << "Remove CPUSETs..." << std::flush;
    cpuset_remove_partition("general");
    cpuset_remove_partition("hpc");
    std::cerr << " DONE" << std::endl;
}

void my_signal_handler(int signum)
{
    switch (signum) {
        case 2: /* SIGINT */
        case 6: /* SIGABRT */
            cpuset_remove();
            signal(signum, SIG_DFL);
            raise(signum);

        default:;
    }
}
#endif

void usage(const char *prg)
{
    std::cout << "Usage:\n\n\t" << prg << " <SIZE> <DISTRIBUTION>\n\nwhere\n\tdistribution\tone of DENSE, SPARSE, GRID\n";
}


int main(int argc, char **argv)
{
    if (argc <= 1) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    if (strcmp(argv[1], "-h") == 0 or strcmp(argv[1], "--help") == 0) {
        usage(argv[0]);
        exit(EXIT_SUCCESS);
    }

    const unsigned num_threads = atoi(*++argv);
#ifdef LINUX
    cpuset_create(num_threads);
    struct sigaction sigint_handler;
    sigint_handler.sa_handler = my_signal_handler;
    sigaction(/* SIGINT  */ 2, &sigint_handler, nullptr);
    sigaction(/* SIGABRT */ 6, &sigint_handler, nullptr);
    atexit(cpuset_remove);
#endif

#ifndef NDEBUG
    /* Test rewiring. */
    {
        rewiring::memory_file<false> smallpages{"smallpages"};
        rewiring::memory_file<true> hugepages{"hugepages"};

        {
            auto vm_first = smallpages.map(2 * decltype(smallpages)::PAGE_SIZE, 0);
            auto vm_second = smallpages.map(2 * decltype(smallpages)::PAGE_SIZE, smallpages.size);

            int *p0 = (int*) vm_first.addr;
            int *p1 = (int*) vm_second.addr;

            *p0 = 0;
            *p1 = 1;

            swap(vm_first, 0, vm_second, 0);
            assert(*p0 == 1);
            assert(*p1 == 0);

            vm_first.unmap();
            vm_second.unmap();
        }
        {
            auto vm_first = hugepages.map(2 * decltype(hugepages)::PAGE_SIZE, 0);
            auto vm_second = hugepages.map(2 * decltype(hugepages)::PAGE_SIZE, hugepages.size);

            int *p0 = (int*) vm_first.addr;
            int *p1 = (int*) vm_second.addr;

            *p0 = 0;
            *p1 = 1;

            swap(vm_first, 0, vm_second, 0);
            assert(*p0 == 1);
            assert(*p1 == 0);

            vm_first.unmap();
            vm_second.unmap();
        }
    }
#endif

    /* Initialize PAPI. */
    {
        int ret = PAPI_library_init(PAPI_VER_CURRENT);
        if (ret != PAPI_VER_CURRENT && ret > 0) {
            std::cerr << "PAPI library version mismatch!" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (ret < 0) {
            std::cerr << "Initializing PAPI failed: " << PAPI_strerror(ret) << std::endl;
            exit(EXIT_FAILURE);
        }
        if (PAPI_is_initialized() == PAPI_NOT_INITED) {
            std::cerr << "Failed to initialize PAPI." << std::endl;
            exit(EXIT_FAILURE);
        }

        std::cout << "Initialized PAPI.  " << PAPI_num_counters() << " hardware counters available.\n";
    }

    /* Initialize PCM. */
    PCM::getInstance();

    /* Parse command line arguments. */
    std::vector<DISTRIBUTION> distributions;
    for (++argv; *argv; ++argv) {
        std::string s{*argv};
        auto it = STR_TO_DISTRIBUTION.find(s);
        if (it == STR_TO_DISTRIBUTION.end())
            std::cerr << "unknown distribution '" << s << "'\n";
        else
            distributions.push_back(it->second);
    }

    /* Create measurements file. */
#ifdef LINUX
    char hostname[256];
    gethostname(hostname, 256);
#endif
    std::ostringstream path_builder;
    path_builder
        << "partitioning"
#ifdef LINUX
        << '_' << hostname
#endif
        << ".csv";
    std::string path = path_builder.str();
    FILE *file = fopen(path.c_str(), "w");
    assert(file);
    print_title(file);

    for (DISTRIBUTION d : distributions) {
        run_benchmarks<uint32_t>(file, SIZE, num_threads, d);
        run_benchmarks<uint64_t>(file, SIZE, num_threads, d);
    }

    fclose(file);
    exit(EXIT_SUCCESS);
}
