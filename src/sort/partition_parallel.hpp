#pragma once

#include "sort/sort.hpp"
#include "util/monitor.hpp"
#include "util/rewiring.hpp"
#include "util/system.hpp"
#include "util/timer.hpp"
#include <cmath>
#include <thread>
#include <vector>


template<typename T>
struct partition_thread_pirk
{
    std::thread t;
    T *begin_left;
    T *end_left;
    T *begin_right;
    T *end_right;
    T *crack;
};

template<typename T>
void partition_MT_pirk(partition_thread_pirk<T> *t, const T pivot)
{
    partition_predicated_pirk p;
    t->crack = p.AoS(pivot, t->begin_left, t->end_left, t->begin_right, t->end_right);
    assert(t->begin_left <= t->crack);
    assert(t->crack <= t->end_right);
    assert(t->crack < t->end_left or t->begin_right <= t->crack);
}

template<typename T>
T * refined_partition_and_merge_AoS(const T pivot, T * const begin, T * const end, const std::size_t num_threads, uint64_t &tsc_merge)
{
    const std::size_t num_tuples = end - begin;
    const std::size_t slice_size = num_tuples / (2 * num_threads);

    partition_thread_pirk<T> *threads = new partition_thread_pirk<T>[num_threads];
    threads[0].begin_left = begin;
    threads[0].end_right = end;
    for (std::size_t i = 1; i != num_threads; ++i) {
        threads[i - 1].end_left    = threads[i].begin_left = begin + i * slice_size;
        threads[i - 1].begin_right = threads[i].end_right  = end   - i * slice_size;
        assert(begin <= threads[i].begin_left);
        assert(threads[i].end_right <= end);
        assert(threads[i].begin_left < threads[i].end_right);
        assert(threads[i - 1].begin_left  <= threads[i - 1].end_left);
        assert(threads[i - 1].end_left    <  threads[i - 1].begin_right);
        assert(threads[i - 1].begin_right <= threads[i - 1].end_right);
    }
    threads[num_threads - 1].end_left = threads[num_threads - 1].end_right;
    threads[num_threads - 1].begin_right = threads[num_threads - 1].begin_left;

    /* Run a thread per pair of slices. */
    auto f = partition_MT_pirk<T>;
    for (std::size_t i = 0; i != num_threads - 1; ++i)
        threads[i].t = std::thread(f, &threads[i], pivot);
    /* Run partitioning for the innermost slice. */
    threads[num_threads - 1].crack = partition_predicated_pirk{}.AoS(pivot, threads[num_threads - 1].begin_left,
                                                                            threads[num_threads - 1].end_right);
    /* Join the threads. */
    for (std::size_t i = 0; i != num_threads - 1; ++i)
        threads[i].t.join();

    const uint64_t tsc_merge_begin = rdtsc();

#ifndef NDEBUG
    for (std::size_t i = 0; i != num_threads - 1; ++i) {
        partition_thread_pirk<T> &t = threads[i];
        if (t.crack < t.begin_right) {
            assert(verify_partition(t.begin_left, t.crack, t.end_left));
            for (T *p = t.begin_right; p != t.end_right; ++p)
                assert(*p >= pivot);
        } else {
            assert(verify_partition(t.begin_right, t.crack, t.end_right));
            for (T *p = t.begin_left; p != t.end_left; ++p)
                assert(*p < pivot);
        }
    }
    assert(verify_partition(threads[num_threads - 1].begin_left,
                            threads[num_threads - 1].crack,
                            threads[num_threads - 1].end_right));
#endif

    /* Merge the slices. */
    struct step {
        T *crack;
        T *left;
        T *right;
    };
    step *steps = new step[num_threads];

    for (std::size_t i = 0, l = 0, r = num_threads; i != num_threads; ++i) {
        partition_thread_pirk<T> &t = threads[i];
        if (t.crack >= t.begin_right)
            steps[--r] = step{ t.crack, t.begin_right, t.end_right };
        else
            steps[l++] = step{ t.crack, t.begin_left, t.end_left };
    }

    std::size_t l = 0, r = num_threads - 1;
    T *runner_left  = steps[l].crack;
    T *runner_right = steps[r].crack;
    while (runner_left < runner_right) {
        if (runner_right == steps[r].left) {
            runner_right = steps[--r].crack;
        } else if (runner_left == steps[l].right) {
            runner_left = steps[++l].crack;
        } else {
            using std::swap;
            assert(*runner_left >= pivot);
            assert(runner_right[-1] < pivot);
            swap (*runner_left++, *--runner_right);
        }
    }
    assert(*runner_left >= pivot);
    assert(runner_right[-1] < pivot);

    T *crack = *runner_right >= pivot ? runner_right : runner_left;
    assert(begin <= crack);
    assert(crack <= end);
    assert(crack == begin or crack[-1] < pivot);
    assert(crack == end or *crack >= pivot);

    tsc_merge = rdtscp() - tsc_merge_begin;

    delete[] threads;
    delete[] steps;
    return crack;
}


template<bool UseHugepages, std::size_t CHUNK_SIZE, typename T>
struct partition_thread_rewired
{
    std::thread t;
    rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> buffer;
    T *begin;
    T *end;
    T *crack;
};

template<bool UseHugepages, std::size_t CHUNK_SIZE, typename T, bool SIMD = false>
void partition_MT_rewired(partition_thread_rewired<UseHugepages, CHUNK_SIZE, T> *t, const T pivot,
                          rewiring::memory_file<UseHugepages> *memfile,
                          rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> *data)
{
#ifndef NDEBUG
    const auto sum = compute_sum(t->begin, t->end);
#endif
    if (SIMD) {
        partition_rewired_simd<UseHugepages> p{*memfile};
        t->crack = p.AoS(pivot, t->begin, t->end, *data, t->buffer);
    } else {
        partition_rewired<UseHugepages, /*use non-temporal stores */ true> p{*memfile};
        t->crack = p.AoS(pivot, t->begin, t->end, *data, t->buffer);
    }
    assert(t->begin <= t->crack);
    assert(t->crack <= t->end);
    assert(t->crack == t->end or *t->crack >= pivot);
    assert(verify_partition(t->begin, t->crack, t->end));
    assert(sum == compute_sum(t->begin, t->end));
}

template<bool UseHugepages, std::size_t CHUNK_SIZE, typename T, bool SIMD = false>
T * rewired_partition_and_merge_AoS(const T pivot, T * const begin, T * const end,
                                    rewiring::memory_file<UseHugepages> &memfile,
                                    rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &data,
                                    const std::size_t num_threads,
                                    uint64_t &tsc_merge)
{
#ifndef NDEBUG
    const auto sum = compute_sum(begin, end);
#endif

    using thread_t = partition_thread_rewired<UseHugepages, CHUNK_SIZE, T>;
    const std::size_t num_tuples = end - begin;
    constexpr std::size_t ELEMENTS_PER_CHUNK = CHUNK_SIZE / sizeof(T);

    /* Round an index to the next chunk boundary. */
    auto round_to_chunksize = [&](std::size_t pos) {
        const std::size_t chunk = (pos + ELEMENTS_PER_CHUNK / 2) / ELEMENTS_PER_CHUNK;
        return chunk * ELEMENTS_PER_CHUNK;
    };

    /* Create the workload for the threads. */
    auto *threads = new thread_t[num_threads];
    threads[0].begin = begin;
    threads[0].buffer = memfile.template map<CHUNK_SIZE>(4 * CHUNK_SIZE, memfile.size);
    for (std::size_t i = 1; i != num_threads; ++i) {
        auto &t = threads[i];
        t.buffer = memfile.template map<CHUNK_SIZE>(4 * CHUNK_SIZE, memfile.size);
        threads[i - 1].end = t.begin = begin + round_to_chunksize(i * num_tuples / num_threads);
    }
    threads[num_threads - 1].end = end;

    /* Run a thread per slice. */
    auto f = partition_MT_rewired<UseHugepages, CHUNK_SIZE, T, SIMD>;
    for (std::size_t i = 0; i != num_threads; ++i)
        threads[i].t = std::thread(f, &threads[i], pivot, &memfile, &data);
    /* Join the threads. */
    for (std::size_t i = 0; i != num_threads; ++i)
        threads[i].t.join();
    /* Unmap the buffers. */
    for (std::size_t i = 0; i != num_threads; ++i)
        threads[i].buffer.unmap();

    assert(sum == compute_sum(begin, end));

    const uint64_t tsc_merge_begin = rdtsc();

    /* Classify chunks. */
    std::vector<std::size_t> chunks_less;
    std::vector<std::size_t> chunks_greater;
    std::vector<std::size_t> cracks;
    for (std::size_t tid = 0; tid != num_threads; ++tid) {
        auto &t = threads[tid];

        const std::size_t delta = t.crack - begin;
        if (delta & (ELEMENTS_PER_CHUNK - 1)) /* not a whole multiple */
            cracks.push_back(delta);

        for (std::size_t i = (t.begin - begin) / ELEMENTS_PER_CHUNK; i < (t.crack - begin) / ELEMENTS_PER_CHUNK; ++i)
            chunks_less.push_back(i);
        for (std::size_t i = (t.crack - begin + ELEMENTS_PER_CHUNK - 1) / ELEMENTS_PER_CHUNK; i < std::size_t(t.end - begin) / ELEMENTS_PER_CHUNK; ++i)
            chunks_greater.push_back(i);
    }
    assert(is_sorted(cracks.begin(), cracks.end()));

    std::size_t n_rewire = 0;
    while (n_rewire != chunks_less.size() and n_rewire != chunks_greater.size() and
            chunks_greater[n_rewire] < chunks_less[chunks_less.size() - n_rewire - 1])
        ++n_rewire;

    for (std::size_t i = 0; i != n_rewire; ++i) {
        const std::size_t the_chunk_less = chunks_less[chunks_less.size() - n_rewire + i];
        const std::size_t the_chunk_greater = chunks_greater[n_rewire - i - 1];
        using std::swap;
        swap(data, the_chunk_less, data, the_chunk_greater);

        /* check if we introduced a new crack */
        const T *p_left = begin + (the_chunk_greater + 1) * ELEMENTS_PER_CHUNK;
        const T *p_right = begin + the_chunk_less * ELEMENTS_PER_CHUNK;

        if (p_left != end and *p_left >= pivot) {
            const std::size_t crack = (the_chunk_greater + 1) * ELEMENTS_PER_CHUNK;
            cracks.insert(std::upper_bound(cracks.begin(), cracks.end(), crack), crack);
        }
        if (p_right != begin and p_right[-1] < pivot) {
            const std::size_t crack = the_chunk_less * ELEMENTS_PER_CHUNK;
            cracks.insert(std::upper_bound(cracks.begin(), cracks.end(), crack), crack);
        }
    }

    /* Rewire incorrectly placed chunks and swap incorrectly placed elements. */
    std::size_t l = 0, r = cracks.size() - 1;
    T *runner_left  = begin + cracks[l];
    T *runner_right = begin + cracks[r];
    while (runner_left < runner_right) {
        if (*runner_left < pivot) {
            runner_left = begin + cracks[++l];
        } else if (runner_right[-1] >= pivot) {
            runner_right = begin + cracks[--r];
        } else {
            using std::swap;
            assert(*runner_left >= pivot);
            assert(runner_right[-1] < pivot);
            swap(*runner_left++, *--runner_right);
        }
    }

    T *crack = *runner_right >= pivot ? runner_right : runner_left;
    assert(begin <= crack);
    assert(crack <= end);
    assert(crack == begin or crack[-1] < pivot);
    assert(crack == end or *crack >= pivot);
    assert(verify_partition(begin, crack, end));

    tsc_merge = rdtscp() - tsc_merge_begin;

    delete[] threads;
    return crack;
}


template<typename T>
struct partition_thread_hybrid
{
    std::thread t;
    T *begin;
    T *end;
    T *crack;
};

template<typename T>
void partition_MT_hybrid(partition_thread_hybrid<T> *t, const T pivot)
{
    partition_predicated_register p;
    t->crack = p.AoS(pivot, t->begin, t->end);
    assert(t->begin <= t->crack);
    assert(t->crack <= t->end);
    assert(t->crack == t->end or *t->crack >= pivot);
}

template<bool UseHugepages, std::size_t CHUNK_SIZE, typename T>
T * hybrid_partition_and_merge_AoS(const T pivot, T * const begin, T * const end,
                                   rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &data,
                                   const std::size_t num_threads,
                                   uint64_t &tsc_merge)
{
    using thread_t = partition_thread_hybrid<T>;
    const std::size_t num_tuples = end - begin;
    constexpr std::size_t ELEMENTS_PER_CHUNK = CHUNK_SIZE / sizeof(T);

    /* Round an index to the next chunk boundary. */
    auto round_to_chunksize = [&](std::size_t pos) {
        const std::size_t chunk = (pos + ELEMENTS_PER_CHUNK / 2) / ELEMENTS_PER_CHUNK;
        return chunk * ELEMENTS_PER_CHUNK;
    };

    /* Create the workload for the threads. */
    auto *threads = new thread_t[num_threads];
    threads[0].begin = begin;
    for (std::size_t i = 1; i != num_threads; ++i) {
        auto &t = threads[i];
        threads[i - 1].end = t.begin = begin + round_to_chunksize(i * num_tuples / num_threads);
    }
    threads[num_threads - 1].end = end;

    /* Run a thread per slice. */
    auto f = partition_MT_hybrid<T>;
    for (std::size_t i = 0; i != num_threads; ++i) {
        threads[i].t = std::thread(f, &threads[i], pivot);
    }
    /* Join the threads. */
    for (std::size_t i = 0; i != num_threads; ++i)
        threads[i].t.join();

    const uint64_t tsc_merge_begin = rdtsc();

    /* Classify chunks. */
    std::vector<std::size_t> chunks_less;
    std::vector<std::size_t> chunks_greater;
    std::vector<std::size_t> cracks;
    for (std::size_t tid = 0; tid != num_threads; ++tid) {
        auto &t = threads[tid];
        const std::size_t num_chunks = (t.end - t.begin) / ELEMENTS_PER_CHUNK;
        assert(num_chunks * ELEMENTS_PER_CHUNK == std::size_t(t.end - t.begin));

        const std::size_t delta = t.crack - begin;
        if (delta & (CHUNK_SIZE - 1)) /* not a whole multiple */
            cracks.push_back(delta);

        for (std::size_t i = (t.begin - begin) / ELEMENTS_PER_CHUNK; i != (t.crack - begin) / ELEMENTS_PER_CHUNK; ++i)
            chunks_less.push_back(i);
        for (std::size_t i = (t.crack - begin + ELEMENTS_PER_CHUNK - 1) / ELEMENTS_PER_CHUNK; i != std::size_t(t.end - begin) / ELEMENTS_PER_CHUNK; ++i)
            chunks_greater.push_back(i);
    }
    assert(is_sorted(cracks.begin(), cracks.end()));

    std::size_t n_rewire = 0;
    while (n_rewire != chunks_less.size() and n_rewire != chunks_greater.size() and
            chunks_greater[n_rewire] < chunks_less[chunks_less.size() - n_rewire - 1])
        ++n_rewire;

    for (std::size_t i = 0; i != n_rewire; ++i) {
        const std::size_t the_chunk_less = chunks_less[chunks_less.size() - n_rewire + i];
        const std::size_t the_chunk_greater = chunks_greater[n_rewire - i - 1];
        using std::swap;
        swap(data, the_chunk_less, data, the_chunk_greater);

        /* check if we introduced a new crack */
        const T *p_left = begin + (the_chunk_greater + 1) * ELEMENTS_PER_CHUNK;
        const T *p_right = begin + the_chunk_less * ELEMENTS_PER_CHUNK;

        if (p_left != end and *p_left >= pivot) {
            const std::size_t crack = (the_chunk_greater + 1) * ELEMENTS_PER_CHUNK;
            cracks.insert(std::upper_bound(cracks.begin(), cracks.end(), crack), crack);
        }
        if (p_right != begin and p_right[-1] < pivot) {
            const std::size_t crack = the_chunk_less * ELEMENTS_PER_CHUNK;
            cracks.insert(std::upper_bound(cracks.begin(), cracks.end(), crack), crack);
        }
    }

    /* Rewire incorrectly placed chunks and swap incorrectly placed elements. */
    std::size_t l = 0, r = cracks.size() - 1;
    T *runner_left  = begin + cracks[l];
    T *runner_right = begin + cracks[r];
    while (runner_left < runner_right) {
        if (*runner_left < pivot) {
            runner_left = begin + cracks[++l];
        } else if (runner_right[-1] >= pivot) {
            runner_right = begin + cracks[--r];
        } else {
            using std::swap;
            assert(*runner_left >= pivot);
            assert(runner_right[-1] < pivot);
            swap(*runner_left++, *--runner_right);
        }
    }

    T *crack = *runner_right >= pivot ? runner_right : runner_left;
    assert(begin <= crack);
    assert(crack <= end);
    assert(crack == begin or crack[-1] < pivot);
    assert(crack == end or *crack >= pivot);
    assert(verify_partition(begin, crack, end));

    tsc_merge = rdtscp() - tsc_merge_begin;

    delete[] threads;
    return crack;
}
