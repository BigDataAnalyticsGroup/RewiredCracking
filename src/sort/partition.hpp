#pragma once

#include "sort/permute.hpp"
#include "sort/tuple_type.hpp"
#include "util/assert.hpp"
#include "util/checks.hpp"
#include "util/minmax.hpp"
#include "util/rewiring.hpp"
#include "util/timer.hpp"
#include "util/vprint.hpp"
#include <algorithm>
#include <bitset>
#include <chrono>
#include <cstddef>
#include <cstring>
#include <emmintrin.h>
#include <immintrin.h>
#include <iostream>
#include <iterator>
#include <map>


struct partition_branching
{
    template<typename T>
    T * AoS(const T pivot, T *begin, T *end)
    {
        using std::swap;
        while (begin != end) {
            if (*begin < pivot) ++begin;
            else if (end[-1] >= pivot) --end;
            else {
                assert(*begin >= pivot);
                assert(end[-1] < pivot);
                swap(*begin, end[-1]);
                assert(*begin < pivot);
                assert(end[-1] >= pivot);
            }
        }
        assert(begin == end);
        return begin;
    }

    template<typename T>
    T * SoA(const T pivot, T * const begin, T * const end, std::size_t num_tuples)
    {
        using std::swap;

        T *keys_begin = begin;
        T *payloads_begin = begin + num_tuples;
        T *keys_end = end;
        T *payloads_end = end + num_tuples;

        while (keys_begin != keys_end) {
            if (*keys_begin < pivot) {
                ++keys_begin;
                ++payloads_begin;
            } else if (keys_end[-1] >= pivot) {
                --keys_end;
                --payloads_end;
            } else {
                assert(*keys_begin >= pivot);
                assert(keys_end[-1] < pivot);
                swap(*keys_begin, keys_end[-1]);
                swap(*payloads_begin, payloads_end[-1]);
                assert(*keys_begin < pivot);
                assert(keys_end[-1] >= pivot);
            }
        }
        assert(keys_begin == keys_end);
        return keys_begin;
    }
};

struct partition_predicated_naive
{
    template<typename T>
    T * AoS(const T pivot, T *begin, T *end)
    {
        while (begin < end) {
            const T left = *begin;
            const T right = end[-1];
            *begin = right;
            end[-1] = left;
            const ptrdiff_t adv_lo = right < pivot;
            const ptrdiff_t adv_hi = left >= pivot;
            begin += adv_lo;
            end   -= adv_hi;
        }
        return begin;
    }

    template<typename T>
    T * SoA(const T pivot, T * const begin, T * const end, std::size_t num_tuples)
    {
        T *keys_begin = begin;
        T *payloads_begin = begin + num_tuples;
        T *keys_end = end;
        T *payloads_end = end + num_tuples;

        while (keys_begin < keys_end) {
            const T key_left = *keys_begin;
            const T key_right = keys_end[-1];
            const T payload_left = *payloads_begin;
            const T payload_right = payloads_end[-1];

            *keys_begin = key_right;
            keys_end[-1] = key_left;
            *payloads_begin = payload_right;
            payloads_end[-1] = payload_left;

            const ptrdiff_t adv_lo = key_right < pivot;
            const ptrdiff_t adv_hi = key_left >= pivot;
            keys_begin += adv_lo;
            payloads_begin += adv_lo;
            keys_end -= adv_hi;
            payloads_end -= adv_hi;
        }
        return keys_begin;
    }
};

struct partition_predicated_pirk
{
    template<typename T>
    T * AoS(const T pivot, T *begin, T *end)
    {
        struct {
            unsigned which;
            T values[2];
        } localbuffer[2];

        localbuffer[0].which = 0;
        localbuffer[0].values[0] = *begin;
        localbuffer[0].values[1] = end[-1];

        localbuffer[1].which = 1;
        localbuffer[1].values[0] = *begin;
        localbuffer[1].values[1] = end[-1];

        unsigned i = 0;
        while (begin < end) {
            const T value = localbuffer[i].values[localbuffer[i].which];
            *begin = end[-1] = value;
            const unsigned advance_lower  = value <  pivot;
            const unsigned advance_higher = value >= pivot;
            begin += advance_lower;
            end   -= advance_higher;
            localbuffer[i].which = advance_higher;
            localbuffer[i].values[0] = *begin;
            localbuffer[i].values[1] = end[-1];
            i = 1 - i;
        }
        assert(begin == end);
        return begin;
    }

    template<typename T>
    T * SoA(const T pivot, T *begin, T *end, std::size_t num_tuples)
    {
        struct {
            struct {
                T key;
                T payload;
            } values[2];
            unsigned which;
        } localbuffer[2];

        localbuffer[0].values[0].key = *begin;
        localbuffer[0].values[0].payload = begin[num_tuples];
        localbuffer[0].values[1].key = end[-1];
        localbuffer[0].values[1].payload = end[num_tuples - 1];
        localbuffer[0].which = 0;

        localbuffer[1].values[0].key = *begin;
        localbuffer[1].values[0].payload = begin[num_tuples];
        localbuffer[1].values[1].key = end[-1];
        localbuffer[1].values[1].payload = end[num_tuples - 1];
        localbuffer[1].which = 1;

        unsigned i = 0;
        while (begin < end) {
            const T key = localbuffer[i].values[localbuffer[i].which].key;
            const T payload = localbuffer[i].values[localbuffer[i].which].payload;
            *begin = end[-1] = key;
            begin[num_tuples] = end[num_tuples - 1] = payload;
            const unsigned advance_lower  = key <  pivot;
            const unsigned advance_higher = key >= pivot;
            begin += advance_lower;
            end   -= advance_higher;
            localbuffer[i].values[0].key = *begin;
            localbuffer[i].values[0].payload = begin[num_tuples];
            localbuffer[i].values[1].key = end[-1];
            localbuffer[i].values[1].payload = end[num_tuples - 1];
            localbuffer[i].which = advance_higher;
            i = 1 - i;
        }
        assert(begin == end);
        return begin;
    }

    template<typename T>
    T * AoS(const T pivot, T *begin_left, T * end_left, T *begin_right, T *end_right)
    {
        assert(begin_left < end_left);
        assert(begin_right < end_right);
        assert(end_left <= begin_right);

        struct {
            unsigned which;
            T values[2];
        } localbuffer[2];

        localbuffer[0].which = 0;
        localbuffer[0].values[0] = *begin_left;
        localbuffer[0].values[1] = end_right[-1];

        localbuffer[1].which = 1;
        localbuffer[1].values[0] = *begin_left;
        localbuffer[1].values[1] = end_right[-1];

        T *p_left  = begin_left;
        T *p_right = end_right;
        unsigned i = 0;
        while (p_left < p_right) {
            const T value = localbuffer[i].values[localbuffer[i].which];
            *p_left = p_right[-1] = value;
            const unsigned advance_lower  = value <  pivot;
            const unsigned advance_higher = value >= pivot;
            p_left  += advance_lower;
            p_right -= advance_higher;
            if (p_left == end_left) p_left = begin_right;
            if (p_right == begin_right) p_right = end_left;
            localbuffer[i].which = advance_higher;
            localbuffer[i].values[0] = *p_left;
            localbuffer[i].values[1] = p_right[-1];
            i = 1 - i;
        }
        assert(p_left == p_right);
        assert(p_left != end_left);
        return p_left;
    }
};

struct partition_predicated_decoupled
{
    template<typename T>
    T * AoS(const T pivot, T *begin, T *end)
    {
        struct {
            ptrdiff_t which;
            T values[2];
        } buf[2] = {
            {0, {*begin, end[-1]}},
            {1, {*begin, end[-1]}},
        };
        ptrdiff_t i = 0;

        while (begin != end) {
            const T value = buf[i].values[buf[i].which];
            *begin = end[-1] = value;
            buf[i].values[1] = begin[1];
            buf[i].values[0] = end[-2];
            const ptrdiff_t advance_lower = value < pivot;
            begin += advance_lower;
            end   += ptrdiff_t(-1) + advance_lower;
            buf[i].which = advance_lower;
            i = ptrdiff_t(1) ^ i;
        }
        assert(begin == end);
        return begin;
    }

    template<typename T>
    T * SoA(const T pivot, T * const begin, T * const end, std::size_t num_tuples)
    {
        T *keys_begin = begin;
        T *payloads_begin = begin + num_tuples;
        T *keys_end = end;
        T *payloads_end = end + num_tuples;

        struct {
            ptrdiff_t which;
            T keys[2];
            T payloads[2];
        } buf[2] = {
            { 0, { *keys_begin, keys_end[-1] }, { *payloads_begin, payloads_end[-1] } },
            { 1, { *keys_begin, keys_end[-1] }, { *payloads_begin, payloads_end[-1] } },
        };
        ptrdiff_t i = 0;

        while (keys_begin != keys_end) {
            const T key = buf[i].keys[buf[i].which];
            const T payload = buf[i].payloads[buf[i].which];
            *keys_begin = keys_end[-1] = key;
            *payloads_begin = payloads_end[-1] = payload;
            buf[i].keys[1] = keys_begin[1];
            buf[i].keys[0] = keys_end[-2];
            buf[i].payloads[1] = payloads_begin[1];
            buf[i].payloads[0] = payloads_end[-2];
            const ptrdiff_t advance_lower = key < pivot;
            keys_begin += advance_lower;
            keys_end   += ptrdiff_t(-1) + advance_lower;
            payloads_begin += advance_lower;
            payloads_end   += ptrdiff_t(-1) + advance_lower;
            buf[i].which = advance_lower;
            i = ptrdiff_t(1) ^ i;
        }
        assert(keys_begin == keys_end);
        return keys_begin;
    }
};

struct partition_predicated_register
{
    template<typename T>
    T * AoS(const T pivot, T *begin, T *end)
    {
        /* Corner case handling for odd number of elements. */
        if ((end - begin) & 0x1) {
            if (end[-1] >= pivot)
                --end;
            else {
                using std::swap;
                swap(end[-1], *begin);
                ++begin;
            }
        }
        assert(not(end - begin & 0x1) && "not a multiple of 2");

        T first = *begin;
        T second = end[-1];

        while (begin < end) {
            {
                *begin = end[-1] = first;
                const T left = begin[1];
                const T right = end[-2];
                const ptrdiff_t advance_lower = first < pivot;
                begin += advance_lower;
                end   += ptrdiff_t(-1) + advance_lower;
                first = advance_lower ? left : right;
            }
            {
                *begin = end[-1] = second;
                const T left = begin[1];
                const T right = end[-2];
                const ptrdiff_t advance_lower = second < pivot;
                begin += advance_lower;
                end   += ptrdiff_t(-1) + advance_lower;
                second = advance_lower ? left : right;
            }
        }

        return begin;
    }

    template<typename T>
    T * SoA(const T pivot, T * const begin, T * const end, std::size_t num_tuples)
    {
        T *keys_begin = begin;
        T *payloads_begin = begin + num_tuples;
        T *keys_end = end;
        T *payloads_end = end + num_tuples;

        /* Corner case handling for odd number of elements. */
        if ((keys_end - keys_begin) & 0x1) {
            if (keys_end[-1] >= pivot) {
                --keys_end;
                --payloads_end;
            } else {
                using std::swap;
                swap(keys_end[-1], *keys_begin);
                swap(payloads_end[-1], *payloads_begin);
                ++keys_begin;
                ++payloads_begin;
            }
        }
        assert(not(keys_end - keys_begin & 0x1) && "not a multiple of 2");

        T key_first = *keys_begin;
        T payload_first = *payloads_begin;
        T key_second = keys_end[-1];
        T payload_second = payloads_end[-1];

        while (keys_begin < keys_end) {
            {
                *keys_begin = keys_end[-1] = key_first;
                *payloads_begin = payloads_end[-1] = payload_first;
                const T key_left = keys_begin[1];
                const T payload_left = payloads_begin[1];
                const T key_right = keys_end[-2];
                const T payload_right = payloads_end[-2];
                const ptrdiff_t advance_lower = key_first < pivot;
                keys_begin += advance_lower;
                keys_end   += ptrdiff_t(-1) + advance_lower;
                payloads_begin += advance_lower;
                payloads_end   += ptrdiff_t(-1) + advance_lower;
                key_first = advance_lower ? key_left : key_right;
                payload_first = advance_lower ? payload_left : payload_right;
            }
            {
                *keys_begin = keys_end[-1] = key_second;
                *payloads_begin = payloads_end[-1] = payload_second;
                const T key_left = keys_begin[1];
                const T payload_left = payloads_begin[1];
                const T key_right = keys_end[-2];
                const T payload_right = payloads_end[-2];
                const ptrdiff_t advance_lower = key_second < pivot;
                keys_begin += advance_lower;
                keys_end   += ptrdiff_t(-1) + advance_lower;
                payloads_begin += advance_lower;
                payloads_end   += ptrdiff_t(-1) + advance_lower;
                key_second = advance_lower ? key_left : key_right;
                payload_second = advance_lower ? payload_left : payload_right;
            }
        }

        return keys_begin;
    }
};

template<bool UseHugepages, bool UseNonTemporalStore = false>
struct partition_rewired
{
    partition_rewired(rewiring::memory_file<UseHugepages> &memfile) : memfile(memfile) { }

    template<typename T, std::size_t CHUNK_SIZE>
    T * AoS(const T pivot, T * const begin, T * const end,
            rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &data,
            rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &vm_buffer) const
    {
        assert(vm_buffer.size == 4 * CHUNK_SIZE && "buffer of wrong size");
        using namespace rewiring;
        constexpr std::size_t ELEMENTS_PER_CHUNK = CHUNK_SIZE / sizeof(T);

        T * const p_buffer = static_cast<T*>(vm_buffer.addr);
        T * const p_buffer_lo_begin = p_buffer + 0 * ELEMENTS_PER_CHUNK;
        T * const p_buffer_hi_end = p_buffer + 4 * ELEMENTS_PER_CHUNK;

        enum {
            BUFFER_LO_FIRST  = 0,
            BUFFER_LO_SECOND = 1,
            BUFFER_HI_SECOND = 2,
            BUFFER_HI_FIRST  = 3,
        };

        /* Define pointers. */
        T *p_buf_lo = p_buffer_lo_begin;
        T *p_buf_hi = p_buffer_hi_end;

        T * const p_data = static_cast<T*>(data.addr);

        const std::size_t chunk_lo_begin = (begin - p_data) / ELEMENTS_PER_CHUNK;
        const std::size_t chunk_hi_end   = (end - p_data - 1) / ELEMENTS_PER_CHUNK;
        assert(p_data + chunk_lo_begin * ELEMENTS_PER_CHUNK <= begin);
        assert(end <= p_data + (chunk_hi_end + 1) * ELEMENTS_PER_CHUNK);
        assert(p_data + chunk_hi_end * ELEMENTS_PER_CHUNK < end);
        std::size_t next_chunk_lo = chunk_lo_begin;
        std::size_t next_chunk_hi = chunk_hi_end;
        assert(p_data + next_chunk_lo * ELEMENTS_PER_CHUNK <= begin);
        assert(p_data + (next_chunk_lo + 1) * ELEMENTS_PER_CHUNK > begin);
        assert(p_data + next_chunk_hi * ELEMENTS_PER_CHUNK < end);
        assert(p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK >= end);

#ifndef NDEBUG
        const auto sum_before  = compute_sum(p_data + chunk_lo_begin * ELEMENTS_PER_CHUNK, begin);
        const auto sum_after   = compute_sum(end, p_data + (chunk_hi_end + 1) * ELEMENTS_PER_CHUNK);
        const auto sum_between = compute_sum(begin, end);
#endif

        T *p_src_lo = p_data + next_chunk_lo * ELEMENTS_PER_CHUNK;
        T *p_src_hi = p_data + next_chunk_hi * ELEMENTS_PER_CHUNK;

        auto store = [](void *dest, T &value) {
            if constexpr (UseNonTemporalStore) {
                /* non-temporal store */
                if constexpr (sizeof(T) == 8)
                    _mm_stream_si64((long long*) dest, *((long long*) &value));
                else
                    _mm_stream_si128((__m128i*) dest, *((__m128i*) &value));
            } else {
                *((T*) dest) = value;
            }
        };
        auto insert = [&](T &value) {
            if constexpr (UseNonTemporalStore) { /* use predicated execution, avoid speculation on mem ops */
                const ptrdiff_t to_left = value < pivot;
                T *dest = to_left ? p_buf_lo : p_buf_hi - 1;
                store(dest, value);
                p_buf_lo += to_left;
                p_buf_hi += ptrdiff_t(-1) + to_left;
            } else { /* use speculative execution */
                if (value < pivot)
                    *p_buf_lo++ = value;
                else
                    *--p_buf_hi = value;
            }
        };

        bool isFirst;
        bool fromLo;

        {
            /* memcpy portions outside [begin,end) to the buffers */
            T * const start_lo = p_data + next_chunk_lo * ELEMENTS_PER_CHUNK;
            T * const end_hi = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK;
            assert(begin - start_lo < long(ELEMENTS_PER_CHUNK));
            assert(end_hi - end < long(ELEMENTS_PER_CHUNK));

            if (start_lo != begin) {
                for (auto p = start_lo; p != begin; ++p) {
                    assert(*p < pivot);
                    store(p_buf_lo++, *p);
                }
#ifndef NDEBUG
                for (auto p = p_buffer_lo_begin; p != p_buf_lo; ++p)
                    assert(*p < pivot);
#endif
            }
            if (end_hi != end) {
                for (auto p = end_hi; p != end; --p) {
                    assert(p[-1] >= pivot);
                    store(--p_buf_hi, p[-1]);
                }
#ifndef NDEBUG
                for (auto p = p_buf_hi; p != p_buffer_hi_end; ++p)
                    assert(*p >= pivot);
#endif
            }

            /* partition the rest of the pages */
            if (next_chunk_lo == next_chunk_hi) {
                for (auto p = begin; p != end; ++p) {
                    insert(*p);
                }
                p_src_lo += ELEMENTS_PER_CHUNK;
            } else {
                if (start_lo != begin) {
                    for (auto p = begin; p != start_lo + ELEMENTS_PER_CHUNK; ++p) {
                        insert(*p);
                    }
                    p_src_lo += ELEMENTS_PER_CHUNK;
                }
                if (end_hi != end) {
                    for (auto p = end_hi - ELEMENTS_PER_CHUNK; p != end; ++p) {
                        insert(*p);
                    }
                    p_src_hi -= ELEMENTS_PER_CHUNK;
                }
            }

#ifndef NDEBUG
            for (auto p = p_buffer_lo_begin; p != p_buf_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_hi; p != p_buffer_hi_end; ++p)
                assert(*p >= pivot);
#endif

            /* If one page of a buffer has been filled, swap the virtual memory region of this buffer into the virtual
             * memory region of the input data. */
            const bool is_full_lo = p_buf_lo - p_buffer_lo_begin >= long(ELEMENTS_PER_CHUNK);
            const bool is_full_hi = p_buffer_hi_end - p_buf_hi >= long(ELEMENTS_PER_CHUNK);
            if (is_full_lo) {
                swap(vm_buffer, BUFFER_LO_FIRST, data, next_chunk_lo);
                swap(vm_buffer, BUFFER_LO_FIRST, vm_buffer, BUFFER_LO_SECOND);
                p_buf_lo -= ELEMENTS_PER_CHUNK;
                ++next_chunk_lo;
            }
            if (is_full_hi) {
                swap(vm_buffer, BUFFER_HI_FIRST, data, next_chunk_hi);
                swap(vm_buffer, BUFFER_HI_FIRST, vm_buffer, BUFFER_HI_SECOND);
                p_buf_hi += ELEMENTS_PER_CHUNK;
                --next_chunk_hi;
            }

#ifndef NDEBUG
            for (auto p = begin; p < p_data + next_chunk_lo * ELEMENTS_PER_CHUNK; ++p)
                assert(*p < pivot);
            for (auto p = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK; p < end; ++p)
                assert(*p >= pivot);
            for (auto p = p_buffer_lo_begin; p != p_buf_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_hi; p != p_buffer_hi_end; ++p)
                assert(*p >= pivot);
#endif

            if (is_full_lo and is_full_hi) {
                isFirst = true;
                fromLo = true;
            } else if (is_full_lo) {
                isFirst = false;
                fromLo = true;
            } else if (is_full_hi) {
                isFirst = false;
                fromLo = false;
            } else {
                isFirst = start_lo == begin and end_hi == end;
                fromLo = start_lo == begin;
            }
        }

        /* Repeatedly partiton a chunk at a time into the buffers.  Rewire buffers if they run full. */
        while (p_src_lo <= p_src_hi) {
            T *p_src = fromLo ? p_src_lo : p_src_hi;
            assert((p_src - p_data) % ELEMENTS_PER_CHUNK == 0);

            /* Partition one chunk from the input data into buffer_lo/buffer_hi. */
            for (const T *end = p_src + ELEMENTS_PER_CHUNK; p_src != end; ++p_src)
                insert(*p_src);

#ifndef NDEBUG
            for (auto p = p_buffer_lo_begin; p != p_buf_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_hi; p != p_buffer_hi_end; ++p)
                assert(*p >= pivot);
#endif

            /* Update pointers. */
            if (not fromLo)
                p_src -= 2lu * ELEMENTS_PER_CHUNK;
            (fromLo ? p_src_lo : p_src_hi) = p_src;

            /* If one page of a buffer has been filled, swap the virtual memory region of this buffer into the virtual
             * memory region of the input data. */
            const bool is_full_lo = p_buf_lo - p_buffer_lo_begin >= long(ELEMENTS_PER_CHUNK);
            const bool is_full_hi = p_buffer_hi_end - p_buf_hi >= long(ELEMENTS_PER_CHUNK);
            if (is_full_lo) {
                assert(p_src_lo >= p_data + (next_chunk_lo + 1) * ELEMENTS_PER_CHUNK);
                swap(vm_buffer, BUFFER_LO_FIRST, data, next_chunk_lo);
                swap(vm_buffer, BUFFER_LO_FIRST, vm_buffer, BUFFER_LO_SECOND);
                p_buf_lo -= ELEMENTS_PER_CHUNK;
                ++next_chunk_lo;
            }
            if (is_full_hi) {
                assert(p_src_hi + ELEMENTS_PER_CHUNK <= p_data + next_chunk_hi * ELEMENTS_PER_CHUNK);
                swap(vm_buffer, BUFFER_HI_FIRST, data, next_chunk_hi);
                swap(vm_buffer, BUFFER_HI_FIRST, vm_buffer, BUFFER_HI_SECOND);
                p_buf_hi += ELEMENTS_PER_CHUNK;
                --next_chunk_hi;
            }

#ifndef NDEBUG
            for (auto p = begin; p < p_data + next_chunk_lo * ELEMENTS_PER_CHUNK; ++p)
                assert(*p < pivot);
            for (auto p = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK; p < end; ++p)
                assert(*p >= pivot);
            for (auto p = p_buffer_lo_begin; p != p_buf_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_hi; p != p_buffer_hi_end; ++p)
                assert(*p >= pivot);
#endif

            if (isFirst and not is_full_lo and not is_full_hi) {
                isFirst = false;
                fromLo = not fromLo;
            } else if (is_full_lo and is_full_hi) {
                isFirst = true;
                fromLo = true;
            } else if (is_full_lo)
                fromLo = true;
            else if (is_full_hi)
                fromLo = false;
        }

#ifndef NDEBUG
        /* Verify partitioning so far was correct. */
        for (auto p = begin; p < p_data + next_chunk_lo * ELEMENTS_PER_CHUNK; ++p)
            assert(*p < pivot);
        for (auto p = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK; p < end; ++p)
            assert(*p >= pivot);

        /* Verify the remaining data inside the buffers was placed correctly. */
        for (auto p = p_buffer_lo_begin; p != p_buf_lo; ++p)
            assert(*p < pivot);
        for (auto p = p_buf_hi; p != p_buffer_hi_end; ++p)
            assert(*p >= pivot);
#endif

        const std::size_t num_elements_buffer_lo = p_buf_lo - p_buffer_lo_begin;
        const std::size_t num_elements_buffer_hi = p_buffer_hi_end - p_buf_hi;

        assert(num_elements_buffer_lo < ELEMENTS_PER_CHUNK);
        assert(num_elements_buffer_hi < ELEMENTS_PER_CHUNK);

        std::size_t next_chunk = fromLo ? next_chunk_lo : next_chunk_hi;
        T *partition = p_data + next_chunk * ELEMENTS_PER_CHUNK + num_elements_buffer_lo;
        assert(num_elements_buffer_lo != 0 or num_elements_buffer_hi == 0);
        if (num_elements_buffer_lo != 0 or num_elements_buffer_hi != 0) {
            assert(begin <= partition);
            assert(partition <= end);

            /* Merge buffer_hi into buffer_lo.  This will exactly fill the first page of buffer_lo. */
            assert(num_elements_buffer_lo + num_elements_buffer_hi == ELEMENTS_PER_CHUNK);
            std::memcpy(p_buf_lo, p_buf_hi, sizeof(T) * num_elements_buffer_hi);
            assert(verify_partition(p_buffer_lo_begin, p_buf_lo, p_buffer_lo_begin + ELEMENTS_PER_CHUNK));
            assert(*p_buf_lo >= pivot);
            assert(*(p_buf_lo - 1) < pivot);

            /* Swap the first page of buffer_lo into the virtual address space of the data. */
            assert(vm_buffer.addr == (void*) p_buffer_lo_begin);
            swap(vm_buffer, BUFFER_LO_FIRST, data, next_chunk);

            assert(*partition >= pivot);
            assert(*(partition - 1) < pivot);

#ifndef NDEBUG
            /* Verify partitioning so far was correct. */
            for (auto p = begin; p != partition; ++p)
                assert(*p < pivot);
            for (auto p = partition; p != end; ++p)
                assert(*p >= pivot);
#endif
        }

#ifndef NDEBUG
        assert(sum_before == compute_sum(p_data + chunk_lo_begin * ELEMENTS_PER_CHUNK, begin));
        assert(sum_after  == compute_sum(end, p_data + (chunk_hi_end + 1) * ELEMENTS_PER_CHUNK));
        assert(sum_between == compute_sum(begin, end));
#endif
        return partition;
    }

    template<typename T, std::size_t CHUNK_SIZE>
    T * SoA(const T pivot, T * const begin, T * const end,
            rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &data,
            rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &vm_buffer) const
    {
        assert(vm_buffer.size == 8 * CHUNK_SIZE, "buffer of wrong size");
        assert(data.size % CHUNK_SIZE == 0, "input data of wrong size");
        using namespace rewiring;
        constexpr std::size_t ELEMENTS_PER_CHUNK = CHUNK_SIZE / sizeof(T);
        const std::size_t num_tuples = data.size / 2 / sizeof(T);
        const std::size_t num_chunks_offset = data.size / CHUNK_SIZE / 2;

        T * const p_buffer = static_cast<T*>(vm_buffer.addr);
        T * const p_buffer_keys_lo_begin = p_buffer + 0 * ELEMENTS_PER_CHUNK;
        T * const p_buffer_keys_hi_end   = p_buffer + 4 * ELEMENTS_PER_CHUNK;
        T * const p_buffer_payloads_lo_begin = p_buffer + 4 * ELEMENTS_PER_CHUNK;
        T * const p_buffer_payloads_hi_end   = p_buffer + 8 * ELEMENTS_PER_CHUNK;

        enum {
            BUFFER_KEYS_LO_FIRST  = 0,
            BUFFER_KEYS_LO_SECOND = 1,
            BUFFER_KEYS_HI_SECOND = 2,
            BUFFER_KEYS_HI_FIRST  = 3,

            BUFFER_PAYLOADS_LO_FIRST  = 4,
            BUFFER_PAYLOADS_LO_SECOND = 5,
            BUFFER_PAYLOADS_HI_SECOND = 6,
            BUFFER_PAYLOADS_HI_FIRST  = 7,
        };

        /* Define pointers. */
        T *p_buf_keys_lo = p_buffer_keys_lo_begin;
        T *p_buf_keys_hi = p_buffer_keys_hi_end;
        T *p_buf_payloads_lo = p_buffer_payloads_lo_begin;
        T *p_buf_payloads_hi = p_buffer_payloads_hi_end;

        T * const p_data = static_cast<T*>(data.addr);
        T * const p_keys = p_data;
#ifndef NDEBUG
        T * const p_payloads = p_data + num_tuples;
#endif

#ifndef NDEBUG
        const auto sum_before  = compute_sum(p_data, begin);
        const auto sum_after   = compute_sum(end, p_data + num_tuples);
        const auto sum_between = compute_sum(begin, end);
#endif

        std::size_t next_chunk_lo = (begin - p_keys) / ELEMENTS_PER_CHUNK;
        std::size_t next_chunk_hi = (end - p_keys - 1) / ELEMENTS_PER_CHUNK;
        assert(p_keys + next_chunk_lo * ELEMENTS_PER_CHUNK <= begin);
        assert(p_keys + (next_chunk_lo + 1) * ELEMENTS_PER_CHUNK > begin);
        assert(p_keys + next_chunk_hi * ELEMENTS_PER_CHUNK < end);
        assert(p_keys + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK >= end);

        T *p_src_keys_lo = p_data + next_chunk_lo * ELEMENTS_PER_CHUNK;
        T *p_src_keys_hi = p_data + next_chunk_hi * ELEMENTS_PER_CHUNK;
        T *p_src_payloads_lo = p_src_keys_lo + num_tuples;
        T *p_src_payloads_hi = p_src_keys_hi + num_tuples;

#if 0
        auto insert = [&](T *p) {
            assert(begin <= p);
            assert(p < end);
            const T key = *p;
            const T payload = p[num_tuples];
            if (key < pivot) {
                *p_buf_keys_lo++ = key;
                *p_buf_payloads_lo++ = payload;
            } else {
                *--p_buf_keys_hi = key;
                *--p_buf_payloads_hi = payload;
            }
        };
#else
        auto insert = [&](T *p) {
            assert(begin <= p);
            assert(p < end);
            const T key = *p;
            const T payload = p[num_tuples];
            const ptrdiff_t to_lower = key < pivot;
            *p_buf_keys_lo = p_buf_keys_hi[-1] = key;
            *p_buf_payloads_lo = p_buf_payloads_hi[-1] = payload;
            p_buf_keys_lo += to_lower;
            p_buf_payloads_lo += to_lower;
            p_buf_keys_hi += ptrdiff_t(-1) + to_lower;
            p_buf_payloads_hi += ptrdiff_t(-1) + to_lower;
        };
#endif

        bool isFirst;
        bool fromLo;

        {
            /* memcpy portions outside [begin,end) to the buffers */
            T * const start_lo = p_keys + next_chunk_lo * ELEMENTS_PER_CHUNK;
            T * const end_hi = p_keys + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK;
            const auto num_lo = begin - start_lo;
            const auto num_hi = end_hi - end;
            assert(num_lo < long(ELEMENTS_PER_CHUNK));
            assert(num_hi < long(ELEMENTS_PER_CHUNK));

            if (num_lo) {
#ifndef NDEBUG
                for (auto p = start_lo; p != begin; ++p)
                    assert(*p < pivot);
#endif
                std::memcpy(p_buf_keys_lo,     start_lo,              num_lo * sizeof(T));
                std::memcpy(p_buf_payloads_lo, start_lo + num_tuples, num_lo * sizeof(T));
                p_buf_keys_lo += num_lo;
                p_buf_payloads_lo += num_lo;
#ifndef NDEBUG
                for (auto p = p_buffer_keys_lo_begin; p != p_buf_keys_lo; ++p)
                    assert(*p < pivot);
#endif
            }

            if (num_hi) {
#ifndef NDEBUG
                for (auto p = end; p != end_hi; ++p)
                    assert(*p >= pivot);
#endif
                std::memcpy(p_buf_keys_hi - num_hi,     end,              num_hi * sizeof(T));
                std::memcpy(p_buf_payloads_hi - num_hi, end + num_tuples, num_hi * sizeof(T));
                p_buf_keys_hi -= num_hi;
                p_buf_payloads_hi -= num_hi;
#ifndef NDEBUG
                for (auto p = p_buf_keys_hi; p != p_buffer_keys_hi_end; ++p)
                    assert(*p >= pivot);
#endif
            }

            /* partition the rest of the pages */
            if (next_chunk_lo == next_chunk_hi) {
                for (auto p = begin; p != end; ++p)
                    insert(p);
                p_src_keys_lo += ELEMENTS_PER_CHUNK;
                p_src_payloads_lo += ELEMENTS_PER_CHUNK;
            } else {
                if (num_lo) {
                    for (auto p = begin; p != start_lo + ELEMENTS_PER_CHUNK; ++p)
                        insert(p);
                    p_src_keys_lo += ELEMENTS_PER_CHUNK;
                    p_src_payloads_lo += ELEMENTS_PER_CHUNK;
                }
                if (num_hi) {
                    for (auto p = end_hi - ELEMENTS_PER_CHUNK; p != end; ++p)
                        insert(p);
                    p_src_keys_hi -= ELEMENTS_PER_CHUNK;
                    p_src_payloads_hi -= ELEMENTS_PER_CHUNK;
                }
            }

            for (auto p = p_buffer_keys_lo_begin; p != p_buf_keys_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_keys_hi; p != p_buffer_keys_hi_end; ++p)
                assert(*p >= pivot);

            /* If one page of a buffer has been filled, swap the virtual memory region of this buffer into the virtual
             * memory region of the input data. */
            const bool is_full_lo = p_buf_keys_lo - p_buffer_keys_lo_begin >= long(ELEMENTS_PER_CHUNK);
            const bool is_full_hi = p_buffer_keys_hi_end - p_buf_keys_hi >= long(ELEMENTS_PER_CHUNK);
            if (is_full_lo) {
                swap(vm_buffer, BUFFER_KEYS_LO_FIRST, data, next_chunk_lo);
                swap(vm_buffer, BUFFER_PAYLOADS_LO_FIRST, data, next_chunk_lo + num_chunks_offset);
                swap(vm_buffer, BUFFER_KEYS_LO_FIRST, vm_buffer, BUFFER_KEYS_LO_SECOND);
                swap(vm_buffer, BUFFER_PAYLOADS_LO_FIRST, vm_buffer, BUFFER_PAYLOADS_LO_SECOND);
                p_buf_keys_lo -= ELEMENTS_PER_CHUNK;
                p_buf_payloads_lo -= ELEMENTS_PER_CHUNK;
                ++next_chunk_lo;
            }
            if (is_full_hi) {
                swap(vm_buffer, BUFFER_KEYS_HI_FIRST, data, next_chunk_hi);
                swap(vm_buffer, BUFFER_PAYLOADS_HI_FIRST, data, next_chunk_hi + num_chunks_offset);
                swap(vm_buffer, BUFFER_KEYS_HI_FIRST, vm_buffer, BUFFER_KEYS_HI_SECOND);
                swap(vm_buffer, BUFFER_PAYLOADS_HI_FIRST, vm_buffer, BUFFER_PAYLOADS_HI_SECOND);
                p_buf_keys_hi += ELEMENTS_PER_CHUNK;
                p_buf_payloads_hi += ELEMENTS_PER_CHUNK;
                --next_chunk_hi;
            }

#ifndef NDEBUG
            for (auto p = p_data; p != p_data + next_chunk_lo * ELEMENTS_PER_CHUNK; ++p)
                assert(*p < pivot);
            for (auto p = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK; p != p_data + num_tuples; ++p)
                assert(*p >= pivot);
            for (auto p = p_buffer_keys_lo_begin; p != p_buf_keys_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_keys_hi; p != p_buffer_keys_hi_end; ++p)
                assert(*p >= pivot);
#endif

            if (is_full_lo and is_full_hi) {
                isFirst = true;
                fromLo = true;
            } else if (is_full_lo) {
                isFirst = false;
                fromLo = true;
            } else if (is_full_hi) {
                isFirst = false;
                fromLo = false;
            } else {
                isFirst = start_lo == begin and end_hi == end;
                fromLo = start_lo == begin;
            }
        }

        while (p_src_keys_lo <= p_src_keys_hi) {
            T *p_src = fromLo ? p_src_keys_lo : p_src_keys_hi;

            /* Partition one page from the input data into buffer_lo/buffer_hi. */
            for (const T *end = p_src + ELEMENTS_PER_CHUNK; p_src != end; ++p_src)
                insert(p_src);

#ifndef NDEBUG
            for (auto p = p_buffer_keys_lo_begin; p != p_buf_keys_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_keys_hi; p != p_buffer_keys_hi_end; ++p)
                assert(*p >= pivot);
#endif

            /* Update pointers. */
            if (not fromLo)
                p_src -= 2lu * ELEMENTS_PER_CHUNK;
            (fromLo ? p_src_keys_lo : p_src_keys_hi) = p_src;

            /* If one page of a buffer has been filled, swap the virtual memory region of this buffer into the virtual
             * memory region of the input data. */
            const bool is_full_lo = p_buf_keys_lo - p_buffer_keys_lo_begin >= long(ELEMENTS_PER_CHUNK);
            const bool is_full_hi = p_buffer_keys_hi_end - p_buf_keys_hi >= long(ELEMENTS_PER_CHUNK);
            if (is_full_lo) {
                swap(vm_buffer, BUFFER_KEYS_LO_FIRST, data, next_chunk_lo);
                swap(vm_buffer, BUFFER_KEYS_LO_FIRST, vm_buffer, BUFFER_KEYS_LO_SECOND);
                swap(vm_buffer, BUFFER_PAYLOADS_LO_FIRST, data, next_chunk_lo + num_chunks_offset);
                swap(vm_buffer, BUFFER_PAYLOADS_LO_FIRST, vm_buffer, BUFFER_PAYLOADS_LO_SECOND);
                p_buf_keys_lo -= ELEMENTS_PER_CHUNK;
                p_buf_payloads_lo -= ELEMENTS_PER_CHUNK;
                ++next_chunk_lo;
            }
            if (is_full_hi) {
                swap(vm_buffer, BUFFER_KEYS_HI_FIRST, data, next_chunk_hi);
                swap(vm_buffer, BUFFER_KEYS_HI_FIRST, vm_buffer, BUFFER_KEYS_HI_SECOND);
                swap(vm_buffer, BUFFER_PAYLOADS_HI_FIRST, data, next_chunk_hi + num_chunks_offset);
                swap(vm_buffer, BUFFER_PAYLOADS_HI_FIRST, vm_buffer, BUFFER_PAYLOADS_HI_SECOND);
                p_buf_keys_hi += ELEMENTS_PER_CHUNK;
                p_buf_payloads_hi += ELEMENTS_PER_CHUNK;
                --next_chunk_hi;
            }

            if (isFirst and not is_full_lo and not is_full_hi) {
                isFirst = false;
                fromLo = not fromLo;
            } else if (is_full_lo and is_full_hi) {
                isFirst = true;
                fromLo = true;
            } else if (is_full_lo)
                fromLo = true;
            else if (is_full_hi)
                fromLo = false;

#ifndef NDEBUG
            for (auto p = p_keys; p != p_data + next_chunk_lo * ELEMENTS_PER_CHUNK; ++p)
                assert(*p < pivot);
            for (auto p = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK; p != p_payloads; ++p)
                assert(*p >= pivot);
            for (auto p = p_buffer_keys_lo_begin; p != p_buf_keys_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_keys_hi; p != p_buffer_keys_hi_end; ++p)
                assert(*p >= pivot);
#endif
        }

#ifndef NDEBUG
        /* Verify partitioning so far was correct. */
        for (auto p = p_keys; p != p_data + next_chunk_lo * ELEMENTS_PER_CHUNK; ++p)
            assert(*p < pivot);
        for (auto p = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK; p != p_payloads; ++p)
            assert(*p >= pivot);

        /* Verify the remaining data inside the buffers was placed correctly. */
        for (auto p = p_buffer_keys_lo_begin; p != p_buf_keys_lo; ++p)
            assert(*p < pivot);
        for (auto p = p_buf_keys_hi; p != p_buffer_keys_hi_end; ++p)
            assert(*p >= pivot);
#endif

        const std::size_t num_elements_buffer_lo = p_buf_keys_lo - p_buffer_keys_lo_begin;
        const std::size_t num_elements_buffer_hi = p_buffer_keys_hi_end - p_buf_keys_hi;

        assert(num_elements_buffer_lo < ELEMENTS_PER_CHUNK);
        assert(num_elements_buffer_hi < ELEMENTS_PER_CHUNK);

        T *partition;

        /* If the unlikely situation occurs that both buffers were filled and swapped, we are done. */
        if (num_elements_buffer_lo == 0 and num_elements_buffer_hi == 0) {
            assert(next_chunk_lo > next_chunk_hi);
            partition = p_data + next_chunk_lo * ELEMENTS_PER_CHUNK;
            goto exit;
        }

        assert(next_chunk_lo >= next_chunk_hi);

        {
            /* Compute the position of the partition. */
            //std::size_t next_chunk = fromLo ? next_chunk_lo : next_chunk_hi;
            partition = p_data + next_chunk_lo * ELEMENTS_PER_CHUNK + num_elements_buffer_lo;

            /* Merge buffer_hi into buffer_lo.  This will exactly fill the first page of buffer_lo. */
            assert(num_elements_buffer_lo + num_elements_buffer_hi == ELEMENTS_PER_CHUNK);
            std::memcpy(p_buf_keys_lo, p_buf_keys_hi, sizeof(T) * num_elements_buffer_hi);
            std::memcpy(p_buf_payloads_lo, p_buf_payloads_hi, sizeof(T) * num_elements_buffer_hi);
            assert(verify_partition(p_buffer_keys_lo_begin, p_buf_keys_lo, p_buffer_keys_lo_begin + ELEMENTS_PER_CHUNK));
            assert(*p_buf_keys_lo >= pivot);
            assert(*(p_buf_keys_lo - 1) < pivot);

            /* Swap the first page of buffer_lo into the virtual address space of the data. */
            assert(vm_buffer.addr == (void*) p_buffer_keys_lo_begin);
            swap(vm_buffer, BUFFER_KEYS_LO_FIRST, data, next_chunk_lo);
            swap(vm_buffer, BUFFER_PAYLOADS_LO_FIRST, data, next_chunk_lo + num_chunks_offset);

            assert(partition == end or *partition >= pivot);
            assert(partition == begin or *(partition - 1) < pivot);
        }

exit:
#ifndef NDEBUG
        assert(sum_before  == compute_sum(p_data, begin));
        assert(sum_after   == compute_sum(end, p_data + num_tuples));
        assert(sum_between == compute_sum(begin, end));
#endif
        return partition;
    }

    private:
    rewiring::memory_file<UseHugepages> &memfile;
};

#if __AVX2__
template<bool UseHugepages>
struct partition_rewired_simd
{
    partition_rewired_simd(rewiring::memory_file<UseHugepages> &memfile) : memfile(memfile) { }

    template<typename T, std::size_t CHUNK_SIZE>
    T * AoS(const T pivot, T * const begin, T * const end,
            rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &data,
            rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &vm_buffer) const
    {
        assert(vm_buffer.size == 4 * CHUNK_SIZE && "buffer of wrong size");
        using namespace rewiring;
        constexpr std::size_t ELEMENTS_PER_CHUNK = CHUNK_SIZE / sizeof(T);
        const std::size_t num_tuples = data.size / sizeof(T);
#ifdef NDEBUG
        (void) num_tuples;
#endif

        T * const p_buffer = static_cast<T*>(vm_buffer.addr);
        T * const p_buffer_lo_begin = p_buffer + 0 * ELEMENTS_PER_CHUNK;
        T * const p_buffer_hi_end = p_buffer + 4 * ELEMENTS_PER_CHUNK;

        enum {
            BUFFER_LO_FIRST  = 0,
            BUFFER_LO_SECOND = 1,
            BUFFER_HI_SECOND = 2,
            BUFFER_HI_FIRST  = 3,
        };

        /* Define pointers. */
        T *p_buf_lo = p_buffer_lo_begin;
        T *p_buf_hi = p_buffer_hi_end;

        T * const p_data = static_cast<T*>(data.addr);

#ifndef NDEBUG
        const auto sum_before  = compute_sum(p_data, begin);
        const auto sum_after   = compute_sum(end, p_data + num_tuples);
        const auto sum_between = compute_sum(begin, end);
#endif

        std::size_t next_chunk_lo = (begin - p_data) / ELEMENTS_PER_CHUNK;
        std::size_t next_chunk_hi = (end - p_data - 1) / ELEMENTS_PER_CHUNK;
        assert(p_data + next_chunk_lo * ELEMENTS_PER_CHUNK <= begin);
        assert(p_data + (next_chunk_lo + 1) * ELEMENTS_PER_CHUNK > begin);
        assert(p_data + next_chunk_hi * ELEMENTS_PER_CHUNK < end);
        assert(p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK >= end);

        T *p_src_lo = p_data + next_chunk_lo * ELEMENTS_PER_CHUNK;
        T *p_src_hi = p_data + next_chunk_hi * ELEMENTS_PER_CHUNK;

        auto insert = [&](T value) {
            if (value < pivot)
                *p_buf_lo++ = value;
            else
                *--p_buf_hi = value;
        };

        bool isFirst;
        bool fromLo;

        {
            /* memcpy portions outside [begin,end) to the buffers */
            T * const start_lo = p_data + next_chunk_lo * ELEMENTS_PER_CHUNK;
            T * const end_hi = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK;
            assert(begin - start_lo < long(ELEMENTS_PER_CHUNK));
            assert(end_hi - end < long(ELEMENTS_PER_CHUNK));

            if (start_lo != begin) {
#ifndef NDEBUG
                for (auto p = start_lo; p != begin; ++p)
                    assert(*p < pivot);
#endif
                std::memcpy(p_buf_lo, start_lo, (begin - start_lo) * sizeof(T));
                p_buf_lo += (begin - start_lo);
#ifndef NDEBUG
                for (auto p = p_buffer_lo_begin; p != p_buf_lo; ++p)
                    assert(*p < pivot);
#endif
            }
            if (end_hi != end) {
#ifndef NDEBUG
                for (auto p = end; p != end_hi; ++p)
                    assert(*p >= pivot);
#endif
                std::memcpy(p_buf_hi - (end_hi - end), end, (end_hi - end) * sizeof(T));
                p_buf_hi -= (end_hi - end);
#ifndef NDEBUG
                for (auto p = p_buf_hi; p != p_buffer_hi_end; ++p)
                    assert(*p >= pivot);
#endif
            }

            /* partition the rest of the pages */
            if (next_chunk_lo == next_chunk_hi) {
                for (auto p = begin; p != end; ++p)
                    insert(*p);
                p_src_lo += ELEMENTS_PER_CHUNK;
            } else {
                if (start_lo != begin) {
                    for (auto p = begin; p != start_lo + ELEMENTS_PER_CHUNK; ++p)
                        insert(*p);
                    p_src_lo += ELEMENTS_PER_CHUNK;
                }
                if (end_hi != end) {
                    for (auto p = end_hi - ELEMENTS_PER_CHUNK; p != end; ++p)
                        insert(*p);
                    p_src_hi -= ELEMENTS_PER_CHUNK;
                }
            }

#ifndef NDEBUG
            for (auto p = p_buffer_lo_begin; p != p_buf_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_hi; p != p_buffer_hi_end; ++p)
                assert(*p >= pivot);
#endif

            /* If one page of a buffer has been filled, swap the virtual memory region of this buffer into the virtual
             * memory region of the input data. */
            const bool is_full_lo = p_buf_lo - p_buffer_lo_begin >= long(ELEMENTS_PER_CHUNK);
            const bool is_full_hi = p_buffer_hi_end - p_buf_hi >= long(ELEMENTS_PER_CHUNK);
            if (is_full_lo) {
                swap(vm_buffer, BUFFER_LO_FIRST, data, next_chunk_lo);
                swap(vm_buffer, BUFFER_LO_FIRST, vm_buffer, BUFFER_LO_SECOND);
                p_buf_lo -= ELEMENTS_PER_CHUNK;
                ++next_chunk_lo;
            }
            if (is_full_hi) {
                swap(vm_buffer, BUFFER_HI_FIRST, data, next_chunk_hi);
                swap(vm_buffer, BUFFER_HI_FIRST, vm_buffer, BUFFER_HI_SECOND);
                p_buf_hi += ELEMENTS_PER_CHUNK;
                --next_chunk_hi;
            }

#ifndef NDEBUG
            for (auto p = p_data; p != p_data + next_chunk_lo * ELEMENTS_PER_CHUNK; ++p)
                assert(*p < pivot);
            for (auto p = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK; p != p_data + num_tuples; ++p)
                assert(*p >= pivot);
            for (auto p = p_buffer_lo_begin; p != p_buf_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_hi; p != p_buffer_hi_end; ++p)
                assert(*p >= pivot);
#endif

            if (is_full_lo and is_full_hi) {
                isFirst = true;
                fromLo = true;
            } else if (is_full_lo) {
                isFirst = false;
                fromLo = true;
            } else if (is_full_hi) {
                isFirst = false;
                fromLo = false;
            } else {
                isFirst = start_lo == begin and end_hi == end;
                fromLo = start_lo == begin;
            }
        }

        /* SIMD constants */
        constexpr int LANES = 32 / sizeof(typename T::value_type);
        const __m256i vpivot = LANES == 8 ? _mm256_set1_epi32(pivot.key) : _mm256_set1_epi64x(pivot.key);

        /* Repeatedly partiton a chunk at a time into the buffers.  Rewire buffers if they run full. */
        while (p_src_lo <= p_src_hi) {
            T *p_src = fromLo ? p_src_lo : p_src_hi;
            assert((p_src - p_data) % ELEMENTS_PER_CHUNK == 0);

            /* Partition one chunk from the input data into buffer_lo/buffer_hi. */
            for (const T *end = p_src + ELEMENTS_PER_CHUNK; p_src != end; p_src += LANES) {
                const __m256i *p_vdata = (__m256i*) p_src;
                const __m256i vfirst  = _mm256_load_si256(p_vdata);
                const __m256i vsecond = _mm256_load_si256(p_vdata + 1);

                if constexpr (LANES == 8) {
                    const __m256i vkeys     = _mm256_shuffle_ps(vfirst, vsecond, _MM_SHUFFLE(2, 0, 2, 0));
                    const __m256i vpayloads = _mm256_shuffle_ps(vfirst, vsecond, _MM_SHUFFLE(3, 1, 3, 1));

                    const __m256i vcmpmask = _mm256_cmpgt_epi32(vpivot, vkeys); // 1 means less than pivot
                    const int mask = _mm256_movemask_ps(vcmpmask);
                    const int num_less = __builtin_popcount(mask);
                    assert(0 <= num_less);
                    assert(num_less <= LANES);

                    const uint64_t shuffle = *reinterpret_cast<const uint64_t*>(permutation8_packable_table + mask);
                    const __m256i vshuffle = _mm256_cvtepu8_epi32(_mm_cvtsi64_si128(shuffle));

                    const __m256i vpartitioned_keys = _mm256_permutevar8x32_epi32(vkeys, vshuffle);
                    const __m256i vpartitioned_payloads = _mm256_permutevar8x32_epi32(vpayloads, vshuffle);

                    const __m256i vunpacked_lo = _mm256_unpacklo_epi32(vpartitioned_keys, vpartitioned_payloads);
                    const __m256i vunpacked_hi = _mm256_unpackhi_epi32(vpartitioned_keys, vpartitioned_payloads);

#ifndef NDEBUG
                    T *p_lo = (T*) &vunpacked_lo;
                    T *p_hi = (T*) &vunpacked_hi;
                    switch (num_less) {
                        case 8: assert(p_hi[3] < pivot);
                        case 7: assert(p_hi[2] < pivot);
                        case 6: assert(p_hi[1] < pivot);
                        case 5: assert(p_hi[0] < pivot);
                        case 4: assert(p_lo[3] < pivot);
                        case 3: assert(p_lo[2] < pivot);
                        case 2: assert(p_lo[1] < pivot);
                        case 1: assert(p_lo[0] < pivot);
                    }
#endif

                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_lo),     vunpacked_lo);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_lo) + 1, vunpacked_hi);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_hi) - 1, vunpacked_hi);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_hi) - 2, vunpacked_lo);

                    p_buf_hi -= (LANES - num_less);
#ifndef NDEBUG
                    for (int i = 0; i != num_less; ++i)
                        assert(p_buf_lo[i] < pivot);
                    for (int i = 0; i != LANES - num_less; ++i)
                        assert(p_buf_hi[i] >= pivot);
#endif
                    p_buf_lo += num_less;
                }
                if constexpr (LANES == 4) {
                    const __m256i vkeys     = _mm256_shuffle_pd(vfirst, vsecond, _MM_SHUFFLE(0, 0, 0, 0));
                    const __m256i vpayloads = _mm256_shuffle_pd(vfirst, vsecond, _MM_SHUFFLE(1, 1, 1, 1));

                    const __m256i vcmpmask = _mm256_cmpgt_epi64(vpivot, vkeys); // 1 means less than pivot
                    const int mask = _mm256_movemask_pd(vcmpmask);
                    const int num_less = __builtin_popcount(mask);
                    assert(0 <= num_less);
                    assert(num_less <= LANES);

                    const uint64_t shuffle = *reinterpret_cast<const uint64_t*>(permutation4_packable_table + mask);
                    const __m256i vshuffle = _mm256_cvtepu8_epi32(_mm_cvtsi64_si128(shuffle));

                    const __m256i vpartitioned_keys = _mm256_permutevar8x32_epi32(vkeys, vshuffle);
                    const __m256i vpartitioned_payloads = _mm256_permutevar8x32_epi32(vpayloads, vshuffle);

                    const __m256i vunpacked_lo = _mm256_unpacklo_epi64(vpartitioned_keys, vpartitioned_payloads);
                    const __m256i vunpacked_hi = _mm256_unpackhi_epi64(vpartitioned_keys, vpartitioned_payloads);

#ifndef NDEBUG
                    T *p_lo = (T*) &vunpacked_lo;
                    T *p_hi = (T*) &vunpacked_hi;
                    switch (num_less) {
                        case 4: assert(p_hi[1] < pivot);
                        case 3: assert(p_hi[0] < pivot);
                        case 2: assert(p_lo[1] < pivot);
                        case 1: assert(p_lo[0] < pivot);
                    }
#endif

                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_lo),     vunpacked_lo);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_lo) + 1, vunpacked_hi);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_hi) - 1, vunpacked_hi);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_hi) - 2, vunpacked_lo);

                    p_buf_hi -= (LANES - num_less);
#ifndef NDEBUG
                    for (int i = 0; i != num_less; ++i)
                        assert(p_buf_lo[i] < pivot);
                    for (int i = 0; i != LANES - num_less; ++i)
                        assert(p_buf_hi[i] >= pivot);
#endif
                    p_buf_lo += num_less;
                }
            }

#ifndef NDEBUG
            for (auto p = p_buffer_lo_begin; p != p_buf_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_hi; p != p_buffer_hi_end; ++p)
                assert(*p >= pivot);
#endif

            /* Update pointers. */
            if (not fromLo)
                p_src -= 2lu * ELEMENTS_PER_CHUNK;
            (fromLo ? p_src_lo : p_src_hi) = p_src;

            /* If one page of a buffer has been filled, swap the virtual memory region of this buffer into the virtual
             * memory region of the input data. */
            const bool is_full_lo = p_buf_lo - p_buffer_lo_begin >= long(ELEMENTS_PER_CHUNK);
            const bool is_full_hi = p_buffer_hi_end - p_buf_hi >= long(ELEMENTS_PER_CHUNK);
            if (is_full_lo) {
                assert(p_src_lo >= p_data + (next_chunk_lo + 1) * ELEMENTS_PER_CHUNK);
                swap(vm_buffer, BUFFER_LO_FIRST, data, next_chunk_lo);
                swap(vm_buffer, BUFFER_LO_FIRST, vm_buffer, BUFFER_LO_SECOND);
                p_buf_lo -= ELEMENTS_PER_CHUNK;
                ++next_chunk_lo;
            }
            if (is_full_hi) {
                assert(p_src_hi + ELEMENTS_PER_CHUNK <= p_data + next_chunk_hi * ELEMENTS_PER_CHUNK);
                swap(vm_buffer, BUFFER_HI_FIRST, data, next_chunk_hi);
                swap(vm_buffer, BUFFER_HI_FIRST, vm_buffer, BUFFER_HI_SECOND);
                p_buf_hi += ELEMENTS_PER_CHUNK;
                --next_chunk_hi;
            }

#ifndef NDEBUG
            for (auto p = p_data; p != p_data + next_chunk_lo * ELEMENTS_PER_CHUNK; ++p)
                assert(*p < pivot);
            for (auto p = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK; p != p_data + num_tuples; ++p)
                assert(*p >= pivot);
            for (auto p = p_buffer_lo_begin; p != p_buf_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_hi; p != p_buffer_hi_end; ++p)
                assert(*p >= pivot);
#endif

            if (isFirst and not is_full_lo and not is_full_hi) {
                isFirst = false;
                fromLo = not fromLo;
            } else if (is_full_lo and is_full_hi) {
                isFirst = true;
                fromLo = true;
            } else if (is_full_lo)
                fromLo = true;
            else if (is_full_hi)
                fromLo = false;
        }

#ifndef NDEBUG
        /* Verify partitioning so far was correct. */
        for (auto p = p_data; p != p_data + next_chunk_lo * ELEMENTS_PER_CHUNK; ++p)
            assert(*p < pivot);
        for (auto p = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK; p != p_data + num_tuples; ++p)
            assert(*p >= pivot);

        /* Verify the remaining data inside the buffers was placed correctly. */
        for (auto p = p_buffer_lo_begin; p != p_buf_lo; ++p)
            assert(*p < pivot);
        for (auto p = p_buf_hi; p != p_buffer_hi_end; ++p)
            assert(*p >= pivot);
#endif

        const std::size_t num_elements_buffer_lo = p_buf_lo - p_buffer_lo_begin;
        const std::size_t num_elements_buffer_hi = p_buffer_hi_end - p_buf_hi;

        assert(num_elements_buffer_lo < ELEMENTS_PER_CHUNK);
        assert(num_elements_buffer_hi < ELEMENTS_PER_CHUNK);

        T *partition;

        /* If the unlikely situation occurs that both buffers were filled and swapped, we are done. */
        if (num_elements_buffer_lo == 0 and num_elements_buffer_hi == 0) {
            assert(next_chunk_lo > next_chunk_hi);
            partition = p_data + next_chunk_lo * ELEMENTS_PER_CHUNK;
            assert(begin <= partition);
            assert(partition <= end);
            goto exit;
        }

        assert(next_chunk_lo >= next_chunk_hi);
        {
            /* Compute the position of the partition. */
            std::size_t const next_chunk = fromLo ? next_chunk_lo : next_chunk_hi;
            partition = p_data + next_chunk * ELEMENTS_PER_CHUNK + num_elements_buffer_lo;
            assert(begin <= partition);
            assert(partition <= end);

            /* Merge buffer_hi into buffer_lo.  This will exactly fill the first page of buffer_lo. */
            assert(num_elements_buffer_lo + num_elements_buffer_hi == ELEMENTS_PER_CHUNK);
            std::memcpy(p_buf_lo, p_buf_hi, sizeof(T) * num_elements_buffer_hi);
            assert(verify_partition(p_buffer_lo_begin, p_buf_lo, p_buffer_lo_begin + ELEMENTS_PER_CHUNK));
            assert(*p_buf_lo >= pivot);
            assert(*(p_buf_lo - 1) < pivot);

            /* Swap the first page of buffer_lo into the virtual address space of the data. */
            assert(vm_buffer.addr == (void*) p_buffer_lo_begin);
            swap(vm_buffer, BUFFER_LO_FIRST, data, next_chunk);

            assert(*partition >= pivot);
            assert(*(partition - 1) < pivot);

#ifndef NDEBUG
            /* Verify partitioning so far was correct. */
            for (auto p = p_data; p != partition; ++p)
                assert(*p < pivot);
            for (auto p = partition; p != p_data + num_tuples; ++p)
                assert(*p >= pivot);
#endif
        }

exit:
#ifndef NDEBUG
        assert(sum_before == compute_sum(p_data, begin));
        assert(sum_after == compute_sum(end, p_data + num_tuples));
        assert(sum_between == compute_sum(begin, end));
#endif
        return partition;
    }

    template<typename T, std::size_t CHUNK_SIZE>
    T * SoA(const T pivot, T * const begin, T * const end,
            rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &data,
            rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &vm_buffer) const
    {
        assert(vm_buffer.size == 8 * CHUNK_SIZE, "buffer of wrong size");
        assert(data.size % CHUNK_SIZE == 0, "input data of wrong size");
        using namespace rewiring;
        constexpr std::size_t ELEMENTS_PER_CHUNK = CHUNK_SIZE / sizeof(T);
        const std::size_t num_tuples = data.size / 2 / sizeof(T);
        const std::size_t num_chunks_offset = data.size / CHUNK_SIZE / 2;

        T * const p_buffer = static_cast<T*>(vm_buffer.addr);
        T * const p_buffer_keys_lo_begin = p_buffer + 0 * ELEMENTS_PER_CHUNK;
        T * const p_buffer_keys_hi_end   = p_buffer + 4 * ELEMENTS_PER_CHUNK;
        T * const p_buffer_payloads_lo_begin = p_buffer + 4 * ELEMENTS_PER_CHUNK;
        T * const p_buffer_payloads_hi_end   = p_buffer + 8 * ELEMENTS_PER_CHUNK;

        enum {
            BUFFER_KEYS_LO_FIRST  = 0,
            BUFFER_KEYS_LO_SECOND = 1,
            BUFFER_KEYS_HI_SECOND = 2,
            BUFFER_KEYS_HI_FIRST  = 3,

            BUFFER_PAYLOADS_LO_FIRST  = 4,
            BUFFER_PAYLOADS_LO_SECOND = 5,
            BUFFER_PAYLOADS_HI_SECOND = 6,
            BUFFER_PAYLOADS_HI_FIRST  = 7,
        };

        /* Define pointers. */
        T *p_buf_keys_lo = p_buffer_keys_lo_begin;
        T *p_buf_keys_hi = p_buffer_keys_hi_end;
        T *p_buf_payloads_lo = p_buffer_payloads_lo_begin;
        T *p_buf_payloads_hi = p_buffer_payloads_hi_end;

        T * const p_data = static_cast<T*>(data.addr);
        T * const p_keys = p_data;
#ifndef NDEBUG
        T * const p_payloads = p_data + num_tuples;
#endif

#ifndef NDEBUG
        const auto sum_before  = compute_sum(p_data, begin);
        const auto sum_after   = compute_sum(end, p_data + num_tuples);
        const auto sum_between = compute_sum(begin, end);
#endif

        std::size_t next_chunk_lo = (begin - p_keys) / ELEMENTS_PER_CHUNK;
        std::size_t next_chunk_hi = (end - p_keys - 1) / ELEMENTS_PER_CHUNK;
        assert(p_keys + next_chunk_lo * ELEMENTS_PER_CHUNK <= begin);
        assert(p_keys + (next_chunk_lo + 1) * ELEMENTS_PER_CHUNK > begin);
        assert(p_keys + next_chunk_hi * ELEMENTS_PER_CHUNK < end);
        assert(p_keys + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK >= end);

        /* SIMD constants */
        constexpr int LANES = 32 / sizeof(T);
        const __m256i vpivot = LANES == 8 ? _mm256_set1_epi32(pivot) : _mm256_set1_epi64x(pivot);
        //const __m256i vreverse = LANES == 8 ? _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7) : _mm256_set_epi64x(0, 1, 2, 3);

        T *p_src_keys_lo = p_data + next_chunk_lo * ELEMENTS_PER_CHUNK;
        T *p_src_keys_hi = p_data + next_chunk_hi * ELEMENTS_PER_CHUNK;
        T *p_src_payloads_lo = p_src_keys_lo + num_tuples;
        T *p_src_payloads_hi = p_src_keys_hi + num_tuples;

        auto insert = [&](T *p) {
            const T key = *p;
            const T payload = p[num_tuples];
            if (key < pivot) {
                *p_buf_keys_lo++ = key;
                *p_buf_payloads_lo++ = payload;
            } else {
                *--p_buf_keys_hi = key;
                *--p_buf_payloads_hi = payload;
            }
        };

        bool isFirst;
        bool fromLo;

        {
            /* memcpy portions outside [begin,end) to the buffers */
            T * const start_lo = p_keys + next_chunk_lo * ELEMENTS_PER_CHUNK;
            T * const end_hi = p_keys + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK;
            const auto num_lo = begin - start_lo;
            const auto num_hi = end_hi - end;
            assert(num_lo < long(ELEMENTS_PER_CHUNK));
            assert(num_hi < long(ELEMENTS_PER_CHUNK));

            if (num_lo) {
#ifndef NDEBUG
                for (auto p = start_lo; p != begin; ++p)
                    assert(*p < pivot);
#endif
                std::memcpy(p_buf_keys_lo,     start_lo,              num_lo * sizeof(T));
                std::memcpy(p_buf_payloads_lo, start_lo + num_tuples, num_lo * sizeof(T));
                p_buf_keys_lo += num_lo;
                p_buf_payloads_lo += num_lo;
#ifndef NDEBUG
                for (auto p = p_buffer_keys_lo_begin; p != p_buf_keys_lo; ++p)
                    assert(*p < pivot);
#endif
            }

            if (num_hi) {
#ifndef NDEBUG
                for (auto p = end; p != end_hi; ++p)
                    assert(*p >= pivot);
#endif
                std::memcpy(p_buf_keys_hi - num_hi,     end,              num_hi * sizeof(T));
                std::memcpy(p_buf_payloads_hi - num_hi, end + num_tuples, num_hi * sizeof(T));
                p_buf_keys_hi -= num_hi;
                p_buf_payloads_hi -= num_hi;
#ifndef NDEBUG
                for (auto p = p_buf_keys_hi; p != p_buffer_keys_hi_end; ++p)
                    assert(*p >= pivot);
#endif
            }

            /* partition the rest of the pages */
            if (next_chunk_lo == next_chunk_hi) {
                for (auto p = begin; p != end; ++p)
                    insert(p);
                p_src_keys_lo += ELEMENTS_PER_CHUNK;
                p_src_payloads_lo += ELEMENTS_PER_CHUNK;
            } else {
                if (num_lo) {
                    for (auto p = begin; p != start_lo + ELEMENTS_PER_CHUNK; ++p)
                        insert(p);
                    p_src_keys_lo += ELEMENTS_PER_CHUNK;
                    p_src_payloads_lo += ELEMENTS_PER_CHUNK;
                }
                if (num_hi) {
                    for (auto p = end_hi - ELEMENTS_PER_CHUNK; p != end; ++p)
                        insert(p);
                    p_src_keys_hi -= ELEMENTS_PER_CHUNK;
                    p_src_payloads_hi -= ELEMENTS_PER_CHUNK;
                }
            }

            for (auto p = p_buffer_keys_lo_begin; p != p_buf_keys_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_keys_hi; p != p_buffer_keys_hi_end; ++p)
                assert(*p >= pivot);

            /* If one page of a buffer has been filled, swap the virtual memory region of this buffer into the virtual
             * memory region of the input data. */
            const bool is_full_lo = p_buf_keys_lo - p_buffer_keys_lo_begin >= long(ELEMENTS_PER_CHUNK);
            const bool is_full_hi = p_buffer_keys_hi_end - p_buf_keys_hi >= long(ELEMENTS_PER_CHUNK);
            if (is_full_lo) {
                swap(vm_buffer, BUFFER_KEYS_LO_FIRST, data, next_chunk_lo);
                swap(vm_buffer, BUFFER_PAYLOADS_LO_FIRST, data, next_chunk_lo + num_chunks_offset);
                swap(vm_buffer, BUFFER_KEYS_LO_FIRST, vm_buffer, BUFFER_KEYS_LO_SECOND);
                swap(vm_buffer, BUFFER_PAYLOADS_LO_FIRST, vm_buffer, BUFFER_PAYLOADS_LO_SECOND);
                p_buf_keys_lo -= ELEMENTS_PER_CHUNK;
                p_buf_payloads_lo -= ELEMENTS_PER_CHUNK;
                ++next_chunk_lo;
            }
            if (is_full_hi) {
                swap(vm_buffer, BUFFER_KEYS_HI_FIRST, data, next_chunk_hi);
                swap(vm_buffer, BUFFER_PAYLOADS_HI_FIRST, data, next_chunk_hi + num_chunks_offset);
                swap(vm_buffer, BUFFER_KEYS_HI_FIRST, vm_buffer, BUFFER_KEYS_HI_SECOND);
                swap(vm_buffer, BUFFER_PAYLOADS_HI_FIRST, vm_buffer, BUFFER_PAYLOADS_HI_SECOND);
                p_buf_keys_hi += ELEMENTS_PER_CHUNK;
                p_buf_payloads_hi += ELEMENTS_PER_CHUNK;
                --next_chunk_hi;
            }

#ifndef NDEBUG
            for (auto p = p_data; p != p_data + next_chunk_lo * ELEMENTS_PER_CHUNK; ++p)
                assert(*p < pivot);
            for (auto p = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK; p != p_data + num_tuples; ++p)
                assert(*p >= pivot);
            for (auto p = p_buffer_keys_lo_begin; p != p_buf_keys_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_keys_hi; p != p_buffer_keys_hi_end; ++p)
                assert(*p >= pivot);
#endif

            if (is_full_lo and is_full_hi) {
                isFirst = true;
                fromLo = true;
            } else if (is_full_lo) {
                isFirst = false;
                fromLo = true;
            } else if (is_full_hi) {
                isFirst = false;
                fromLo = false;
            } else {
                isFirst = start_lo == begin and end_hi == end;
                fromLo = start_lo == begin;
            }
        }

        while (p_src_keys_lo <= p_src_keys_hi) {
            T *p_src = fromLo ? p_src_keys_lo : p_src_keys_hi;

            /* Partition one page from the input data into buffer_lo/buffer_hi. */
            for (const T *end = p_src + ELEMENTS_PER_CHUNK; p_src != end; p_src += LANES) {
                const __m256i vkeys = _mm256_load_si256(reinterpret_cast<__m256i*>(p_src));
                const __m256i vpayloads = _mm256_load_si256(reinterpret_cast<__m256i*>(p_src + num_tuples));

                if constexpr (LANES == 8) {
                    const __m256i vcmpmask = _mm256_cmpgt_epi32(vpivot, vkeys);
                    const int mask = _mm256_movemask_ps(vcmpmask);
                    const int num_less = __builtin_popcount(mask);
                    assert(0 <= num_less);
                    assert(num_less <= LANES);

                    const uint64_t shuffle = *reinterpret_cast<const uint64_t*>(permutation8_table + mask);
                    const __m256i vshuffle = _mm256_cvtepu8_epi32(_mm_cvtsi64_si128(shuffle));

                    const __m256i vpartitioned_keys = _mm256_permutevar8x32_epi32(vkeys, vshuffle);
                    const __m256i vpartitioned_payloads = _mm256_permutevar8x32_epi32(vpayloads, vshuffle);

#ifndef NDEBUG
                    for (int i = 0; i != num_less; ++i)
                        assert(reinterpret_cast<const uint32_t*>(&vpartitioned_keys)[i] < pivot);
                    for (int i = num_less; i != LANES; ++i)
                        assert(reinterpret_cast<const uint32_t*>(&vpartitioned_keys)[i] >= pivot);
#endif

                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_keys_lo), vpartitioned_keys);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_payloads_lo), vpartitioned_payloads);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_keys_hi - LANES), vpartitioned_keys);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_payloads_hi - LANES), vpartitioned_payloads);

                    p_buf_keys_hi     -= (LANES - num_less);
                    p_buf_payloads_hi -= (LANES - num_less);
#ifndef NDEBUG
                    for (int i = 0; i != num_less; ++i)
                        assert(p_buf_keys_lo[i] < pivot);
                    for (int i = 0; i != LANES - num_less; ++i)
                        assert(p_buf_keys_hi[i] >= pivot);
#endif
                    p_buf_keys_lo     += num_less;
                    p_buf_payloads_lo += num_less;
                }
                if constexpr (LANES == 4) {
                    const __m256i vcmpmask = _mm256_cmpgt_epi64(vpivot, vkeys);
                    const int mask = _mm256_movemask_pd(vcmpmask);
                    const int num_less = __builtin_popcount(mask);
                    assert(0 <= num_less);
                    assert(num_less <= LANES);

                    const uint64_t shuffle = *reinterpret_cast<const uint64_t*>(permutation4_table + mask);
                    const __m256i vshuffle = _mm256_cvtepu8_epi32(_mm_cvtsi64_si128(shuffle));

                    const __m256i vpartitioned_keys = _mm256_permutevar8x32_epi32(vkeys, vshuffle);
                    const __m256i vpartitioned_payloads = _mm256_permutevar8x32_epi32(vpayloads, vshuffle);

#ifndef NDEBUG
                    for (int i = 0; i != num_less; ++i)
                        assert(reinterpret_cast<const uint64_t*>(&vpartitioned_keys)[i] < pivot);
                    for (int i = num_less; i != LANES; ++i)
                        assert(reinterpret_cast<const uint64_t*>(&vpartitioned_keys)[i] >= pivot);
#endif

                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_keys_lo), vpartitioned_keys);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_payloads_lo), vpartitioned_payloads);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_keys_hi - LANES), vpartitioned_keys);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(p_buf_payloads_hi - LANES), vpartitioned_payloads);

                    p_buf_keys_hi     -= (LANES - num_less);
                    p_buf_payloads_hi -= (LANES - num_less);
#ifndef NDEBUG
                    for (int i = 0; i != num_less; ++i)
                        assert(p_buf_keys_lo[i] < pivot);
                    for (int i = 0; i != LANES - num_less; ++i)
                        assert(p_buf_keys_hi[i] >= pivot);
#endif
                    p_buf_keys_lo     += num_less;
                    p_buf_payloads_lo += num_less;
                }
            }

#ifndef NDEBUG
            for (auto p = p_buffer_keys_lo_begin; p != p_buf_keys_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_keys_hi; p != p_buffer_keys_hi_end; ++p)
                assert(*p >= pivot);
#endif

            /* Update pointers. */
            if (not fromLo)
                p_src -= 2lu * ELEMENTS_PER_CHUNK;
            (fromLo ? p_src_keys_lo : p_src_keys_hi) = p_src;

            /* If one page of a buffer has been filled, swap the virtual memory region of this buffer into the virtual
             * memory region of the input data. */
            const bool is_full_lo = p_buf_keys_lo - p_buffer_keys_lo_begin >= long(ELEMENTS_PER_CHUNK);
            const bool is_full_hi = p_buffer_keys_hi_end - p_buf_keys_hi >= long(ELEMENTS_PER_CHUNK);
            if (is_full_lo) {
                swap(vm_buffer, BUFFER_KEYS_LO_FIRST, data, next_chunk_lo);
                swap(vm_buffer, BUFFER_KEYS_LO_FIRST, vm_buffer, BUFFER_KEYS_LO_SECOND);
                swap(vm_buffer, BUFFER_PAYLOADS_LO_FIRST, data, next_chunk_lo + num_chunks_offset);
                swap(vm_buffer, BUFFER_PAYLOADS_LO_FIRST, vm_buffer, BUFFER_PAYLOADS_LO_SECOND);
                p_buf_keys_lo -= ELEMENTS_PER_CHUNK;
                p_buf_payloads_lo -= ELEMENTS_PER_CHUNK;
                ++next_chunk_lo;
            }
            if (is_full_hi) {
                swap(vm_buffer, BUFFER_KEYS_HI_FIRST, data, next_chunk_hi);
                swap(vm_buffer, BUFFER_KEYS_HI_FIRST, vm_buffer, BUFFER_KEYS_HI_SECOND);
                swap(vm_buffer, BUFFER_PAYLOADS_HI_FIRST, data, next_chunk_hi + num_chunks_offset);
                swap(vm_buffer, BUFFER_PAYLOADS_HI_FIRST, vm_buffer, BUFFER_PAYLOADS_HI_SECOND);
                p_buf_keys_hi += ELEMENTS_PER_CHUNK;
                p_buf_payloads_hi += ELEMENTS_PER_CHUNK;
                --next_chunk_hi;
            }

            if (isFirst and not is_full_lo and not is_full_hi) {
                isFirst = false;
                fromLo = not fromLo;
            } else if (is_full_lo and is_full_hi) {
                isFirst = true;
                fromLo = true;
            } else if (is_full_lo)
                fromLo = true;
            else if (is_full_hi)
                fromLo = false;

#ifndef NDEBUG
            for (auto p = p_keys; p != p_data + next_chunk_lo * ELEMENTS_PER_CHUNK; ++p)
                assert(*p < pivot);
            for (auto p = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK; p != p_payloads; ++p)
                assert(*p >= pivot);
            for (auto p = p_buffer_keys_lo_begin; p != p_buf_keys_lo; ++p)
                assert(*p < pivot);
            for (auto p = p_buf_keys_hi; p != p_buffer_keys_hi_end; ++p)
                assert(*p >= pivot);
#endif
        }

#ifndef NDEBUG
        /* Verify partitioning so far was correct. */
        for (auto p = p_keys; p != p_data + next_chunk_lo * ELEMENTS_PER_CHUNK; ++p)
            assert(*p < pivot);
        for (auto p = p_data + (next_chunk_hi + 1) * ELEMENTS_PER_CHUNK; p != p_payloads; ++p)
            assert(*p >= pivot);

        /* Verify the remaining data inside the buffers was placed correctly. */
        for (auto p = p_buffer_keys_lo_begin; p != p_buf_keys_lo; ++p)
            assert(*p < pivot);
        for (auto p = p_buf_keys_hi; p != p_buffer_keys_hi_end; ++p)
            assert(*p >= pivot);
#endif

        const std::size_t num_elements_buffer_lo = p_buf_keys_lo - p_buffer_keys_lo_begin;
        const std::size_t num_elements_buffer_hi = p_buffer_keys_hi_end - p_buf_keys_hi;

        assert(num_elements_buffer_lo < ELEMENTS_PER_CHUNK);
        assert(num_elements_buffer_hi < ELEMENTS_PER_CHUNK);

        T *partition;

        /* If the unlikely situation occurs that both buffers were filled and swapped, we are done. */
        if (num_elements_buffer_lo == 0 and num_elements_buffer_hi == 0) {
            assert(next_chunk_lo > next_chunk_hi);
            partition = p_data + next_chunk_lo * ELEMENTS_PER_CHUNK;
            goto exit;
        }

        assert(next_chunk_lo >= next_chunk_hi);

        {
            /* Compute the position of the partition. */
            std::size_t next_chunk = fromLo ? next_chunk_lo : next_chunk_hi;
            partition = p_data + next_chunk * ELEMENTS_PER_CHUNK + num_elements_buffer_lo;

            /* Merge buffer_hi into buffer_lo.  This will exactly fill the first page of buffer_lo. */
            assert(num_elements_buffer_lo + num_elements_buffer_hi == ELEMENTS_PER_CHUNK);
            std::memcpy(p_buf_keys_lo, p_buf_keys_hi, sizeof(T) * num_elements_buffer_hi);
            std::memcpy(p_buf_payloads_lo, p_buf_payloads_hi, sizeof(T) * num_elements_buffer_hi);
            assert(verify_partition(p_buffer_keys_lo_begin, p_buf_keys_lo, p_buffer_keys_lo_begin + ELEMENTS_PER_CHUNK));
            assert(*p_buf_keys_lo >= pivot);
            assert(*(p_buf_keys_lo - 1) < pivot);

            /* Swap the first page of buffer_lo into the virtual address space of the data. */
            assert(vm_buffer.addr == (void*) p_buffer_keys_lo_begin);
            swap(vm_buffer, BUFFER_KEYS_LO_FIRST, data, next_chunk_lo);
            swap(vm_buffer, BUFFER_PAYLOADS_LO_FIRST, data, next_chunk_lo + num_chunks_offset);

            assert(*partition >= pivot);
            assert(*(partition - 1) < pivot);
        }

exit:
#ifndef NDEBUG
        assert(sum_before  == compute_sum(p_data, begin));
        assert(sum_after   == compute_sum(end, p_data + num_tuples));
        assert(sum_between == compute_sum(begin, end));
#endif
        return partition;
    }

    private:
    rewiring::memory_file<UseHugepages> &memfile;
};
#endif
