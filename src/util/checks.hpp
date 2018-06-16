#pragma once

#include <array>


template<typename T>
void print(T *arr, std::size_t size)
{
    while (size--)
        std::cout << *arr++ << ", ";
    std::cout << std::endl;
}

template<typename T>
bool verify_partition(const T *begin, const T *mid, const T *end)
{
    T max_lo{std::numeric_limits<T>::lowest()};
    T min_hi{std::numeric_limits<T>::max()};
    assert(max_lo < min_hi);
    for (const T *p = begin; p != mid; ++p)
        max_lo = std::max(max_lo, *p);
    for (const T *p = mid; p != end; ++p)
        min_hi = std::min(min_hi, *p);
    return max_lo <= min_hi;
}

template<typename T>
bool is_sorted(const T *begin, const T *end)
{
    T last = *begin++;
    while (begin != end) {
        if (last > *begin) return false;
        last = *begin++;
    }
    return true;
}

template<typename T>
long long compute_sum(const T *begin, std::size_t n)
{
    long long sum(0);
    while (n--)
        sum += static_cast<long>(*begin++);
    return sum;
}

template<typename T>
long long compute_sum(const T *begin, const T *end) { return compute_sum(begin, end - begin); }

template<typename T>
T median_of_three(const T *begin, const T *end)
{
    const T *mid = begin + (end - begin) / 2;
    std::array<T, 3> samples = {{*begin, *mid, *(end - 1)}};
    std::sort(samples.begin(), samples.end());
    return samples[1];
}
