#include "test/catch.hpp"
#include "util/bit_reversal.hpp"
#include <cstdint>


#define TEST(X, Y) { \
    auto t = reverse_bits(X); \
    REQUIRE(t == Y); \
    REQUIRE(reverse_bits(t) == X); \
}

TEST_CASE( "reverse_bits/32-bit" )
{
    using T = uint32_t;
    TEST(T(0), T(0));
    TEST(T(1), T(1) << 31);
    TEST(T(8), T(1) << 28);
    TEST(T(3), T(3) << 30);
}

TEST_CASE( "reverse_bits/64-bit" )
{
    using T = uint64_t;
    TEST(T(0), T(0));
    TEST(T(1), T(1) << 63);
    TEST(T(8), T(1) << 60);
    TEST(T(3), T(3) << 62);
}
