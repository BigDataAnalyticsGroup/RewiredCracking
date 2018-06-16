#pragma once

#include "util/assert.hpp"
#include <cstdint>
#include <immintrin.h>


void * memcpy_fast(void * __restrict__ dest, const void * __restrict__ const src, std::size_t size)
{
    uint8_t * __restrict__ byte_dest = (uint8_t*) dest;
    uint8_t * __restrict__ byte_src  = (uint8_t*) src;

#ifndef NDEBUG
    constexpr std::size_t BULK_SIZE = 4 * 32; // 4x 256 bit
    assert((byte_dest - byte_src) % BULK_SIZE == 0, "not a multiple of the bulk size");
    assert(((uintptr_t) src) % BULK_SIZE == 0, "not alignned");
#endif

    __m256i *vdest = (__m256i*) byte_dest;
    __m256i *vsrc  = (__m256i*) byte_src;
    const std::size_t n = size / 32;

    for (const __m256i *end = vsrc + n; vsrc != end; vsrc += 4, vdest += 4) {
        /* Prefetch the next 128 bytes. */
        __builtin_prefetch(vsrc +  4, /* rw= */ 0, /* locality= */ 0);
        __builtin_prefetch(vsrc +  5, /* rw= */ 0, /* locality= */ 0);
        __builtin_prefetch(vsrc +  6, /* rw= */ 0, /* locality= */ 0);
        __builtin_prefetch(vsrc +  7, /* rw= */ 0, /* locality= */ 0);

        /* Load 128 bytes to registers (bypassing the caches). */
        __m256i r0 = _mm256_stream_load_si256(vsrc + 0);
        __m256i r1 = _mm256_stream_load_si256(vsrc + 1);
        __m256i r2 = _mm256_stream_load_si256(vsrc + 2);
        __m256i r3 = _mm256_stream_load_si256(vsrc + 3);

        /* Store 128 bytes through main memory (bypassing the caches). */
        _mm256_stream_si256(vdest + 0, r0);
        _mm256_stream_si256(vdest + 1, r1);
        _mm256_stream_si256(vdest + 2, r2);
        _mm256_stream_si256(vdest + 3, r3);
    }
    assert((uint8_t*) vsrc == ((uint8_t*) src) + size);

    return dest;
}
