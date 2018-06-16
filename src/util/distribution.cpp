#include "distribution.hpp"

#include "util/assert.hpp"
#include "util/int.hpp"
#include <algorithm>
#include <iostream>
#include <random>
#include <unordered_set>


/*-- uniform sparse distribution --{{{--------------------------------------------------------------------------------*/
template<typename T>
void create_uniform_sparse(T *data, std::size_t size, const T sentinel, const T max)
{
    static_assert(std::numeric_limits<T>::is_integer, "type T is not an integer type");
    const T min = std::numeric_limits<T>::min();

    std::cout << "Generating uniform sparse data in range [" << min << ", " << max << "]\n";

    T n = max - T(size);
    std::mt19937_64 g(42);
    std::uniform_int_distribution<T> d(min, n);
    std::unordered_set<T> used(2 * size);

    for (std::size_t pos = 0; pos != size; ++pos, ++n) {
        T v = d(g);
        if (v == sentinel || used.count(v))
            v = n;
        assert(min <= v && v <= max);
        data[pos] = v;
        used.insert(v);
    }
}
/*--}}}---------------------------------------------------------------------------------------------------------------*/


/*-- ascending distribution --{{{-------------------------------------------------------------------------------------*/
template<typename T>
void create_ascending(T *data, std::size_t size, const T sentinel)
{
    static_assert(std::numeric_limits<T>::is_integer, "type T is not an integer type");
    const T min = std::numeric_limits<T>::min();

    std::cout << "Generating uniform ascending data from " << min << " to " << (min + size) << "\n";

    T v = min;
    const auto end = data + size;
    while (data != end) {
        v += v == sentinel;
        *data++ = v++;
    }
}
/*--}}}---------------------------------------------------------------------------------------------------------------*/


/*-- uniform dense distribution --{{{--------------------------------------------------------------------------------*/
template<typename T>
void create_uniform_dense(T *data, std::size_t size, const T sentinel)
{
    static_assert(std::numeric_limits<T>::is_integer, "type T is not an integer type");
    using std::shuffle;
    const T min = std::numeric_limits<T>::min();

    std::cout << "Generating uniform dense data from " << min << " to " << (min + size) << "\n";

    T v = min;
    for (T *p = data, *end = data + size; p != end; ++p, ++v) {
        v += v == sentinel;
        *p = v;
    }

    std::mt19937_64 g(42);
    shuffle(data, data + size, g);
}
/*--}}}---------------------------------------------------------------------------------------------------------------*/


/*-- grid distribution --{{{------------------------------------------------------------------------------------------*/
template<typename T>
void create_grid(T *data, std::size_t size, uint8_t bits_per_byte)
{
    static_assert(std::numeric_limits<T>::is_integer, "type T is not an integer type");
    constexpr std::size_t NUM_BYTES = sizeof(T) / sizeof(uint8_t);
    const uint8_t MAX = uint8_t(-1) >> (8 - bits_per_byte);

    std::cout << "Generating grid data with " << bits_per_byte << " bits per byte\n";

    union counter_t {
        T v;
        uint8_t bytes[NUM_BYTES];
    };

    auto shuffle = [](counter_t &c) {
        std::size_t i = NUM_BYTES;
        uint8_t tmp = c.bytes[--i];
        while (i) {
            c.bytes[i] = c.bytes[i - 1];
            --i;
        }
        c.bytes[0] = tmp;
    };

    auto increment = [=](counter_t &c) {
        std::size_t pos = 0;
        for (;;) {
            if (c.bytes[pos] == MAX) {
                c.bytes[pos] = 1;
                ++pos;
            } else {
                ++c.bytes[pos];
                break;
            }
        }
    };

    counter_t counter;
    counter.v = 0;

    for(std::size_t i = 0; i != size; ++i) {
        if (i % NUM_BYTES == 0) increment(counter);
        *data++ = counter.v;
        shuffle(counter);
    }
}
/*--}}}---------------------------------------------------------------------------------------------------------------*/


/* 16 bit */
//template void create_uniform_sparse<uint16_t>(uint16_t*, std::size_t, const uint16_t, const uint16_t);
//template void create_ascending     <uint16_t>(uint16_t*, std::size_t, const uint16_t);
//template void create_uniform_dense <uint16_t>(uint16_t*, std::size_t, const uint16_t);
//template void create_grid          <uint16_t>(uint16_t*, std::size_t, uint8_t);

/* 32 bit */
template void create_uniform_sparse<uint32_t>(uint32_t*, std::size_t, const uint32_t, const uint32_t);
template void create_ascending     <uint32_t>(uint32_t*, std::size_t, const uint32_t);
template void create_uniform_dense <uint32_t>(uint32_t*, std::size_t, const uint32_t);
template void create_grid          <uint32_t>(uint32_t*, std::size_t, uint8_t);

/* 64 bit */
template void create_uniform_sparse<uint64_t>(uint64_t*, std::size_t, const uint64_t, const uint64_t);
template void create_ascending     <uint64_t>(uint64_t*, std::size_t, const uint64_t);
template void create_uniform_dense <uint64_t>(uint64_t*, std::size_t, const uint64_t);
template void create_grid          <uint64_t>(uint64_t*, std::size_t, uint8_t);

/* 128 bit */
//template void create_uniform_sparse<uint128_t>(uint128_t*, std::size_t, const uint128_t, const uint128_t);
//template void create_ascending     <uint128_t>(uint128_t*, std::size_t, const uint128_t);
//template void create_uniform_dense <uint128_t>(uint128_t*, std::size_t, const uint128_t);
//template void create_grid          <uint128_t>(uint128_t*, std::size_t, uint8_t);

/* 256 bit */
//template void create_uniform_sparse<uint256_t>(uint256_t*, std::size_t, const uint256_t, const uint256_t);
//template void create_ascending     <uint256_t>(uint256_t*, std::size_t, const uint256_t);
//template void create_uniform_dense <uint256_t>(uint256_t*, std::size_t, const uint256_t);
//template void create_grid          <uint256_t>(uint256_t*, std::size_t, uint8_t);
