#include "hash/linear_hashing.hpp"
#include "test/catch.hpp"
#include "util/distribution.hpp"


/* General hash table parameters. */
#ifndef KEY_TYPE
#error "must define KEY_TYPE"
#endif

#ifndef VALUE_TYPE
#error "must define VALUE_TYPE"
#endif

#ifndef HASH
#error "must define HASH"
#endif

#ifndef ALLOCATOR
#error "must define ALLOCATOR"
#endif

#ifndef THE_SENTINEL
#error "must define THE_SENTINEL"
#endif

/* Linear hashing specific parameters. */
#ifndef THE_BUCKET_SIZE
#error "must define THE_BUCKET_SIZE"
#endif


#define STRINGIFY(X) #X
#define STR(X) STRINGIFY(X)

#define NAME "linear_hashing_" STR(KEY_TYPE) "_" STR(VALUE_TYPE) "_" STR(HASH) "_" STR(ALLOCATOR)
#define SIZE (std::size_t(1) <<10)


TEST_CASE( NAME )
{
    using HT = linear_hashing<KEY_TYPE, VALUE_TYPE, HASH, std::equal_to<KEY_TYPE>, ALLOCATOR, THE_BUCKET_SIZE, THE_SENTINEL>;
    HT ht{10};
    REQUIRE(ht.raw_capacity() == (1 << 10) * THE_BUCKET_SIZE);
}

TEST_CASE( NAME "/integration" )
{
    using HT = linear_hashing<KEY_TYPE, VALUE_TYPE, HASH, std::equal_to<KEY_TYPE>, ALLOCATOR, THE_BUCKET_SIZE, THE_SENTINEL>;
    KEY_TYPE *data = new KEY_TYPE[SIZE];

    {
        HT ht{8};
        create_uniform_sparse(data, SIZE, KEY_TYPE(0));
        for (std::size_t i = 0; i != SIZE; ++i)
            ht.put(data[i], i);

        for (std::size_t i = 0; i != SIZE; ++i) {
            VALUE_TYPE *p = ht.get(data[i]);
            REQUIRE(*p == i);
        }
    }

    {
        HT ht{8};
        create_uniform_dense(data, SIZE, KEY_TYPE(0));
        for (std::size_t i = 0; i != SIZE; ++i)
            ht.put(data[i], i);

        for (std::size_t i = 0; i != SIZE; ++i) {
            VALUE_TYPE *p = ht.get(data[i]);
            REQUIRE(*p == i);
        }
    }

    {
        HT ht{8};
        create_grid(data, SIZE, 4);
        for (std::size_t i = 0; i != SIZE; ++i)
            ht.put(data[i], i);

        for (std::size_t i = 0; i != SIZE; ++i) {
            VALUE_TYPE *p = ht.get(data[i]);
            REQUIRE(*p == i);
        }
    }

    delete[] data;
}


#undef NAME
#undef KEY_TYPE
#undef VALUE_TYPE
#undef HASH
#undef ALLOCATOR
#undef THE_SENTINEL
#undef THE_BUCKET_SIZE
