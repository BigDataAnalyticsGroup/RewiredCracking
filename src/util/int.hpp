#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <type_traits>


namespace
{
    /**
     * Returns -1 if first is smaller than second, +1 if first is bigger than second, and 0 if first and second are
     * equal.
     */
    template<typename T>
    int int_cmp(T first, T second)
    {
        static_assert(std::is_integral<T>::value, "not an integral type");
        return int(first > second) - int(first < second);
    }
}


template<unsigned N>
struct alignas(8) int_varying
{
    static constexpr unsigned LENGTH = N;

    uint64_t chunks[N];

    /* Constructor */
    int_varying() = default;

    int_varying(uint64_t n)
    {
        chunks[0] = n;
        for (unsigned i = 1; i != N; ++i)
            chunks[i] = 0;
    }

    int_varying(const int_varying &other) { std::copy(other.chunks, other.chunks + N, this->chunks); }

    int_varying & operator=(const int_varying &other) {
        std::copy(other.chunks, other.chunks + N, this->chunks);
        return *this;
    }

    /* Comparison */
    bool operator==(int_varying other) const {
        bool eq = true;
        for (unsigned i = 0; i != N; ++i)
            eq &= this->chunks[i] == other.chunks[i];
        return eq;
    }
    bool operator< (int_varying other) const { return cmp<true,  false>(other); }
    bool operator> (int_varying other) const { return cmp<false, false>(other); }
    bool operator<=(int_varying other) const { return cmp<true,  true >(other); }
    bool operator>=(int_varying other) const { return cmp<false, true >(other); }

    /* Bitwise */
#define BITWISE(OP) \
    int_varying res; \
    for (unsigned i = 0; i != N; ++i) \
        res.chunks[i] = this->chunks[i] OP other.chunks[i]; \
    return res

    int_varying operator&(int_varying other) const { BITWISE(&); }
    int_varying operator|(int_varying other) const { BITWISE(|); }
    int_varying operator^(int_varying other) const { BITWISE(^); }
    int_varying operator~() const {
        int_varying res;
        for (unsigned i = 0; i != N; ++i)
            res.chunks[i] = ~this->chunks[i];
        return res;
    }

    /* Arithmetic */
    int_varying & operator+=(uint64_t n) {
        bool carry;
        unsigned i = 0;
        do {
            carry = __builtin_uaddl_overflow(this->chunks[i], n, &this->chunks[i]);
        } while (carry && i != N);
        return *this;
    }
    int_varying & operator-=(uint64_t n) {
        bool carry;
        unsigned i = 0;
        do {
            carry = __builtin_usubl_overflow(this->chunks[i], n, &this->chunks[i]);
        } while (carry && i != N);
        return *this;
    }

    int_varying & operator+=(int_varying other) {
        bool carry = false;
        for (unsigned i = 0; i != N; ++i) {
            if (carry) {
               carry = __builtin_uaddl_overflow(this->chunks[i], 1lu, &this->chunks[i]);
               this->chunks[i] += other.chunks[i]; // can't overflow *again*
            } else {
                carry = __builtin_uaddl_overflow(this->chunks[i], other.chunks[i], &this->chunks[i]);
            }
        }
        return *this;
    }
    int_varying & operator-=(int_varying other) {
        bool carry = false;
        for (unsigned i = 0; i != N; ++i) {
            if (carry) {
               carry = __builtin_uaddl_overflow(this->chunks[i], -1lu, &this->chunks[i]);
               this->chunks[i] -= other.chunks[i]; // can't overflow *again*
            } else {
                carry = __builtin_usubl_overflow(this->chunks[i], other.chunks[i], &this->chunks[i]);
            }
        }
        return *this;
    }

    int_varying operator+(int_varying other) const {
        int_varying res{*this};
        return res += other;
    }
    int_varying operator-(int_varying other) const {
        int_varying res{*this};
        return res -= other;
    }

    int_varying & operator++() { return this->operator+=(1); }
    int_varying & operator--() { return this->operator-=(1); }
    int_varying operator++(int) {
        int_varying old{*this};
        this->operator++();
        return old;
    }
    int_varying operator--(int) {
        int_varying old{*this};
        this->operator--();
        return old;
    }

    private:
#if 0
    template<bool Less, bool Eq>
    bool cmp(int_varying other) const {
        for (unsigned i = N; i-- != 0;) {
            if (this->chunks[i] < other.chunks[i]) return Less;
            else if (this->chunks[i] > other.chunks[i]) return not Less;
        }
        return Eq;
    }
#else
    template<bool Less, bool Eq>
    bool cmp(int_varying other) const {
        int res = 0;
        for (unsigned i = N; i-- != 0;) {
            res <<= 1;
            const int c = int_cmp(this->chunks[i], other.chunks[i]); // this ? other
            res += c;
        }
        return ((res < 0) & Less) | ((res > 0) & not Less) | ((res == 0) & Eq);
    }
#endif
};


namespace std
{
    template<unsigned N>
    struct numeric_limits<int_varying<N>> : numeric_limits<uint64_t>
    {
        static int_varying<N> max() {
            int_varying<N> res;
            for (unsigned i = 0; i != N; ++i)
                res.chunks[i] = numeric_limits<uint64_t>::max();
            return res;
        }

        static int_varying<N> min() { return int_varying<N>{0}; }
    };

    template<unsigned N>
    struct is_integral<int_varying<N>> : is_integral<uint64_t> { };

    template<unsigned N>
    struct hash<int_varying<N>>
    {
        using result_type = uint64_t;
        result_type operator()(int_varying<N> v) const {
            result_type h = v.chunks[0];
            for (unsigned i = 1; i != N; ++i)
                h = h * 37 + v.chunks[i];
            return h;
        }
    };

    template<unsigned N>
    struct uniform_int_distribution<int_varying<N>> : uniform_int_distribution<uint64_t>
    {
        uniform_int_distribution(int_varying<N> a, int_varying<N> b)
            : uniform_int_distribution<uint64_t>(a.chunks[0], b.chunks[0])
        { }
    };
}


using uint128_t = int_varying<2>; // 16 B
using uint256_t = int_varying<4>; // 32 B
