#pragma once

#include "sort/partition.hpp"
#include "sort/partition_parallel.hpp"
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


template<typename T>
using crackerindex_t = std::map<T, std::size_t>;


template<typename Partitioning, typename T>
void crack_AoS(crackerindex_t<T> &index, tuple_type<T> *data, const std::size_t num_tuples, const T pred,
               Partitioning &p, Timer &timer)
{
    std::size_t lower_bound = 0;
    std::size_t upper_bound = num_tuples;

#ifndef NDEBUG
    tuple_type<T> value_lower = 0;
    tuple_type<T> value_upper = T(-1);
#endif

    auto it = index.upper_bound(pred);
    if (it != index.end()) {
        assert(it->first > pred);
        upper_bound = it->second;
#ifndef NDEBUG
        value_upper = it->first;
        assert(upper_bound < num_tuples);
        for (auto p = data; p != data + it->second; ++p)
            assert(*p < it->first);
        for (auto p = data + it->second; p != data + num_tuples; ++p)
            assert(*p >= it->first);
#endif
    }
    if (it != index.begin()) {
        it = std::prev(it);
        assert(it->first <= pred);
        if (it->first == pred) return; // crack already exists
        lower_bound = it->second;
#ifndef NDEBUG
        value_lower = it->first;
        for (auto p = data; p != data + it->second; ++p)
            assert(*p < it->first);
        for (auto p = data + it->second; p != data + num_tuples; ++p)
            assert(*p >= it->first);
#endif
    }
    assert(lower_bound <= upper_bound);

    if (lower_bound + 1 < upper_bound) {
        timer.start();
        tuple_type<T> *crack = p.AoS(tuple_type<T>(pred), data + lower_bound, data + upper_bound);
        timer.stop();
        assert(data + lower_bound <= crack);
        assert(crack < data + upper_bound);
        assert(crack == data + upper_bound or crack->key >= pred);
        assert(crack == data + lower_bound or (crack - 1)->key < pred);
#ifndef NDEBUG
        for (auto p = data; p != data + lower_bound; ++p)
            assert(*p < value_lower);
        for (auto p = data + lower_bound; p != crack; ++p) {
            assert(*p >= value_lower);
            assert(*p < pred);
        }
        for (auto p = crack; p != data + upper_bound; ++p) {
            assert(*p >= pred);
            assert(*p < value_upper);
        }
        for (auto p = data + upper_bound; p != data + num_tuples; ++p)
            assert(*p >= value_upper);
#endif
        if (crack != data + lower_bound) // test if the predicate created a new crack
            index.emplace(pred, crack - data);
    }
}

template<typename Partitioning, typename T, bool UseHugepages, std::size_t CHUNK_SIZE>
void crack_AoS(crackerindex_t<T> &index, rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &data,
               const std::size_t num_tuples, const T pred, Partitioning &p,
               rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &vm_buffer, Timer &timer)
{
    tuple_type<T> *p_data = static_cast<tuple_type<T>*>(data.addr);
    std::size_t lower_bound = 0;
    std::size_t upper_bound = num_tuples;

#ifndef NDEBUG
    tuple_type<T> value_lower = 0;
    tuple_type<T> value_upper = T(-1);
#endif

    auto it = index.upper_bound(pred);
    if (it != index.end()) {
        assert(it->first > pred);
        upper_bound = it->second;
#ifndef NDEBUG
        value_upper = it->first;
        assert(upper_bound < num_tuples);
        for (auto p = p_data; p != p_data + it->second; ++p)
            assert(*p < it->first);
        for (auto p = p_data + it->second; p != p_data + num_tuples; ++p)
            assert(*p >= it->first);
#endif
    }
    if (it != index.begin()) {
        it = std::prev(it);
        assert(it->first <= pred);
        if (it->first == pred) return; // crack already exists
        lower_bound = it->second;
#ifndef NDEBUG
        value_lower = it->first;
        for (auto p = p_data; p != p_data + it->second; ++p)
            assert(*p < it->first);
        for (auto p = p_data + it->second; p != p_data + num_tuples; ++p)
            assert(*p >= it->first);
#endif
    }
    assert(lower_bound <= upper_bound);

    if (lower_bound + 1 < upper_bound) {
        timer.start();
        tuple_type<T> *crack = p.AoS(tuple_type<T>(pred), p_data + lower_bound, p_data + upper_bound, data, vm_buffer);
        timer.stop();
        assert(p_data + lower_bound <= crack);
        assert(crack < p_data + upper_bound);
        assert(crack == p_data + upper_bound or crack->key >= pred);
        assert(crack == p_data + lower_bound or (crack - 1)->key < pred);
#ifndef NDEBUG
        for (auto p = p_data; p != p_data + lower_bound; ++p)
            assert(*p < value_lower);
        for (auto p = p_data + lower_bound; p != crack; ++p) {
            assert(*p >= value_lower);
            assert(*p < pred);
        }
        for (auto p = crack; p != p_data + upper_bound; ++p) {
            assert(*p >= pred);
            assert(*p < value_upper);
        }
        for (auto p = p_data + upper_bound; p != p_data + num_tuples; ++p)
            assert(*p >= value_upper);
#endif
        if (crack != p_data + lower_bound) // test if the predicate created a new crack
            index.emplace(pred, crack - p_data);
    }
}

template<typename Partitioning, typename T, bool UseHugepages, std::size_t CHUNK_SIZE>
void mixed_crack_AoS(crackerindex_t<T> &index, rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &data,
                     const std::size_t num_tuples, const T pred, Partitioning &p,
                     rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &vm_buffer)
{
    tuple_type<T> *p_data = static_cast<tuple_type<T>*>(data.addr);
    std::size_t lower_bound = 0;
    std::size_t upper_bound = num_tuples;
    constexpr std::size_t MIN_DELTA_FOR_REWIRING = 56lu * 1024 * 1024 / sizeof(tuple_type<T>);

#ifndef NDEBUG
    tuple_type<T> value_lower = 0;
    tuple_type<T> value_upper = T(-1);
#endif

    auto it = index.upper_bound(pred);
    if (it != index.end()) {
        assert(it->first > pred);
        upper_bound = it->second;
#ifndef NDEBUG
        value_upper = it->first;
        assert(upper_bound < num_tuples);
        for (auto p = p_data; p != p_data + it->second; ++p)
            assert(*p < it->first);
        for (auto p = p_data + it->second; p != p_data + num_tuples; ++p)
            assert(*p >= it->first);
#endif
    }
    if (it != index.begin()) {
        it = std::prev(it);
        assert(it->first <= pred);
        if (it->first == pred) return; // crack already exists
        lower_bound = it->second;
#ifndef NDEBUG
        value_lower = it->first;
        for (auto p = p_data; p != p_data + it->second; ++p)
            assert(*p < it->first);
        for (auto p = p_data + it->second; p != p_data + num_tuples; ++p)
            assert(*p >= it->first);
#endif
    }
    assert(lower_bound <= upper_bound);

    if (lower_bound + 1 < upper_bound) {
        tuple_type<T> *crack;
        if (upper_bound - lower_bound > MIN_DELTA_FOR_REWIRING) {
            crack = p.AoS(tuple_type<T>(pred), p_data + lower_bound, p_data + upper_bound, data, vm_buffer);
        } else {
            partition_predicated_register p;
            crack = p.AoS(tuple_type<T>(pred), p_data + lower_bound, p_data + upper_bound);
        }
        assert(p_data + lower_bound <= crack);
        assert(crack < p_data + upper_bound);
        assert(crack == p_data + upper_bound or crack->key >= pred);
        assert(crack == p_data + lower_bound or (crack - 1)->key < pred);
#ifndef NDEBUG
        for (auto p = p_data; p != p_data + lower_bound; ++p)
            assert(*p < value_lower);
        for (auto p = p_data + lower_bound; p != crack; ++p) {
            assert(*p >= value_lower);
            assert(*p < pred);
        }
        for (auto p = crack; p != p_data + upper_bound; ++p) {
            assert(*p >= pred);
            assert(*p < value_upper);
        }
        for (auto p = p_data + upper_bound; p != p_data + num_tuples; ++p)
            assert(*p >= value_upper);
#endif
        if (crack != p_data + lower_bound and crack != p_data + upper_bound and
                (crack - 1)->key < pred and crack->key >= pred) // test if the predicate created a new crack
            index.emplace(pred, crack - p_data);
    }
}

template<typename Partitioning, typename T>
void crack_SoA(crackerindex_t<T> &index, T *data, const std::size_t num_tuples, const T pred,
               Partitioning &p, Timer &timer)
{
    std::size_t lower_bound = 0;
    std::size_t upper_bound = num_tuples;

#ifndef NDEBUG
    T value_lower = 0;
    T value_upper = T(-1);
#endif

    auto it = index.upper_bound(pred);
    if (it != index.end()) {
        assert(it->first > pred);
        upper_bound = it->second;
#ifndef NDEBUG
        value_upper = it->first;
        assert(upper_bound < num_tuples);
        for (auto p = data; p != data + it->second; ++p)
            assert(*p < it->first);
        for (auto p = data + it->second; p != data + num_tuples; ++p)
            assert(*p >= it->first);
#endif
    }
    if (it != index.begin()) {
        it = std::prev(it);
        assert(it->first <= pred);
        if (it->first == pred) return; // crack already exists
        lower_bound = it->second;
#ifndef NDEBUG
        value_lower = it->first;
        for (auto p = data; p != data + it->second; ++p)
            assert(*p < it->first);
        for (auto p = data + it->second; p != data + num_tuples; ++p)
            assert(*p >= it->first);
#endif
    }
    assert(lower_bound <= upper_bound);

    if (lower_bound + 1 < upper_bound) {
        timer.start();
        T *crack = p.SoA(pred, data + lower_bound, data + upper_bound, num_tuples);
        timer.stop();
        assert(data + lower_bound <= crack);
        assert(crack < data + upper_bound);
        assert(crack == data + upper_bound or *crack >= pred);
        assert(crack == data + lower_bound or crack[-1] < pred);
#ifndef NDEBUG
        for (auto p = data; p != data + lower_bound; ++p)
            assert(*p < value_lower);
        for (auto p = data + lower_bound; p != crack; ++p) {
            assert(*p >= value_lower);
            assert(*p < pred);
        }
        for (auto p = crack; p != data + upper_bound; ++p) {
            assert(*p >= pred);
            assert(*p < value_upper);
        }
        for (auto p = data + upper_bound; p != data + num_tuples; ++p)
            assert(*p >= value_upper);
#endif
        if (crack != data + lower_bound) // test if the predicate created a new crack
            index.emplace(pred, crack - data);
    }
}

template<typename Partitioning, typename T, bool UseHugepages, std::size_t CHUNK_SIZE>
void crack_SoA(crackerindex_t<T> &index, rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &data,
               const std::size_t num_tuples, const T pred, Partitioning &p,
               rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &vm_buffer, Timer &timer)
{
    T *p_data = static_cast<T*>(data.addr);
    std::size_t lower_bound = 0;
    std::size_t upper_bound = num_tuples;

#ifndef NDEBUG
    T value_lower = 0;
    T value_upper = T(-1);
#endif

    auto it = index.upper_bound(pred);
    if (it != index.end()) {
        assert(it->first > pred);
        upper_bound = it->second;
#ifndef NDEBUG
        value_upper = it->first;
        assert(upper_bound < num_tuples);
        for (auto p = p_data; p != p_data + it->second; ++p)
            assert(*p < it->first);
        for (auto p = p_data + it->second; p != p_data + num_tuples; ++p)
            assert(*p >= it->first);
#endif
    }
    if (it != index.begin()) {
        it = std::prev(it);
        assert(it->first <= pred);
        if (it->first == pred) return; // crack already exists
        lower_bound = it->second;
#ifndef NDEBUG
        value_lower = it->first;
        for (auto p = p_data; p != p_data + it->second; ++p)
            assert(*p < it->first);
        for (auto p = p_data + it->second; p != p_data + num_tuples; ++p)
            assert(*p >= it->first);
#endif
    }
    assert(lower_bound <= upper_bound);

    if (lower_bound + 1 < upper_bound) {
        T *crack;
        timer.start();
        crack = p.SoA(pred, p_data + lower_bound, p_data + upper_bound, data, vm_buffer);
        timer.stop();
        assert(p_data + lower_bound <= crack);
        assert(crack < p_data + upper_bound);
        assert(crack == p_data + upper_bound or *crack >= pred);
        assert(crack == p_data + lower_bound or crack[-1] < pred);
#ifndef NDEBUG
        for (auto p = p_data; p != p_data + lower_bound; ++p)
            assert(*p < value_lower);
        for (auto p = p_data + lower_bound; p != crack; ++p) {
            assert(*p >= value_lower);
            assert(*p < pred);
        }
        for (auto p = crack; p != p_data + upper_bound; ++p) {
            assert(*p >= pred);
            assert(*p < value_upper);
        }
        for (auto p = p_data + upper_bound; p != p_data + num_tuples; ++p)
            assert(*p >= value_upper);
#endif
        if (crack != p_data + lower_bound and crack != p_data + upper_bound and
                crack[-1] < pred and *crack >= pred) // test if the predicate created a new crack
            index.emplace(pred, crack - p_data);
    }
}

template<typename Partitioning, typename T, bool UseHugepages, std::size_t CHUNK_SIZE>
void mixed_crack_SoA(crackerindex_t<T> &index, rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &data,
                      const std::size_t num_tuples, const T pred, Partitioning &p,
                      rewiring::memory_mapping<UseHugepages, CHUNK_SIZE> &vm_buffer)
{
    T *p_data = static_cast<T*>(data.addr);
    std::size_t lower_bound = 0;
    std::size_t upper_bound = num_tuples;
    constexpr std::size_t MIN_DELTA_FOR_REWIRING = 56lu * 1024 * 1024 / sizeof(T);

#ifndef NDEBUG
    T value_lower = 0;
    T value_upper = T(-1);
#endif

    auto it = index.upper_bound(pred);
    if (it != index.end()) {
        assert(it->first > pred);
        upper_bound = it->second;
#ifndef NDEBUG
        value_upper = it->first;
        assert(upper_bound < num_tuples);
        for (auto p = p_data; p != p_data + it->second; ++p)
            assert(*p < it->first);
        for (auto p = p_data + it->second; p != p_data + num_tuples; ++p)
            assert(*p >= it->first);
#endif
    }
    if (it != index.begin()) {
        it = std::prev(it);
        assert(it->first <= pred);
        lower_bound = it->second;
#ifndef NDEBUG
        value_lower = it->first;
        for (auto p = p_data; p != p_data + it->second; ++p)
            assert(*p < it->first);
        for (auto p = p_data + it->second; p != p_data + num_tuples; ++p)
            assert(*p >= it->first);
#endif
    }
    assert(lower_bound <= upper_bound);

    if (lower_bound + 1 < upper_bound) {
        T *crack;
        if (upper_bound - lower_bound > MIN_DELTA_FOR_REWIRING) {
            crack = p.SoA(pred, p_data + lower_bound, p_data + upper_bound, data, vm_buffer);
        } else {
            partition_predicated_register p;
            crack = p.SoA(pred, p_data + lower_bound, p_data + upper_bound, num_tuples);
        }
        assert(p_data + lower_bound <= crack);
        assert(crack < p_data + upper_bound);
        assert(crack == p_data + upper_bound or *crack >= pred);
        assert(crack == p_data + lower_bound or crack[-1] < pred);
#ifndef NDEBUG
        for (auto p = p_data; p != p_data + lower_bound; ++p)
            assert(*p < value_lower);
        for (auto p = p_data + lower_bound; p != crack; ++p) {
            assert(*p >= value_lower);
            assert(*p < pred);
        }
        for (auto p = crack; p != p_data + upper_bound; ++p) {
            assert(*p >= pred);
            assert(*p < value_upper);
        }
        for (auto p = p_data + upper_bound; p != p_data + num_tuples; ++p)
            assert(*p >= value_upper);
#endif
        if (crack != p_data + lower_bound and crack != p_data + upper_bound and
                crack[-1] < pred and *crack >= pred) // test if the predicate created a new crack
            index.emplace(pred, crack - p_data);
    }
}
