#pragma once

#include "sort/permute.hpp"
#include "sort/tuple_type.hpp"
#include "util/assert.hpp"
#include "util/checks.hpp"
#include "util/minmax.hpp"
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


template<typename Partitioning, typename T>
void qsort(T *begin, T *end, Partitioning p)
{
    using std::swap;
    using std::min;
    using std::max;
    assert(begin < end);

    if (end - begin == 2) {
        if (begin[0] > begin[1])
            swap(begin[0], begin[1]);
        return;
    }

    /* Compute median of three. */
    T l = *begin;
    T r = end[-1];
    T m = *(begin + (end - begin) / 2);
    minmax(l, r);
    minmax(l, m);
    minmax(m, r);
    assert(l <= m);
    assert(m <= r);

    T *mid = p.AoS(m, begin, end); // mid points to an element in the hi part
    assert(*(mid - 1) <= m);
    assert(*mid >= m);
    assert(verify_partition(begin, mid, end));

    if (mid - begin >= 2) qsort<Partitioning>(begin, mid, p);
    if (end - mid   >= 2) qsort<Partitioning>(mid, end, p);
}
