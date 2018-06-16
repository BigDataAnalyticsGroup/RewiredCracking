#include "test/catch.hpp"
#include "util/distribution.hpp"
#include <cstdint>


constexpr std::size_t SIZE = 1 << 16;

TEST_CASE( "distribution" )
{
    using T = uint64_t;
    T *data = reinterpret_cast<T*>(malloc(SIZE * sizeof(T)));

    SECTION( "sparse" ) {
        create_uniform_sparse(data, SIZE, T(-1));
    }
    SECTION( "dense" ) {
        create_uniform_dense(data, SIZE, T(-1));
    }
    SECTION( "grid" ) {
        create_grid(data, SIZE, /* bits_per_byte= */ 4);
    }

    std::sort(data, data + SIZE);
    for (auto p = data, end = data + SIZE - 1; p != end; ++p)
        REQUIRE(p[0] != p[1]);

    free(data);
}
