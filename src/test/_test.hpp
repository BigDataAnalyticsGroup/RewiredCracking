#ifndef TABLE
#error "#define TABLE to some hash table"
#endif

#ifndef HASH
#error "#define HASH to some hash function"
#endif


#include "test/catch.hpp"
#include <cstdint>


#define STRINGIFY(X) #X
#define STR(X) STRINGIFY(X)
#define PREFIX STR(TABLE) "/" STR(HASH)

using table_type = TABLE<uint64_t, uint64_t, HASH<uint64_t>>;


TEST_CASE( PREFIX )
{
    table_type table(-1, 32);
    REQUIRE(table.find(42) == table.end());
    REQUIRE(table.end() == table.end());
}

TEST_CASE( PREFIX "/insert" )
{
    table_type table(-1, 1024);
    for (uint64_t i = 0; i != 42; ++i)
        table.insert(i, i);
    for (uint64_t i = 0; i != 42; ++i) {
        auto it = table.find(i);
        REQUIRE(it != table.end());
        REQUIRE((*it).key   == i);
        REQUIRE((*it).value == i);
    }
    for (uint64_t i = 42; i != 1024; ++i)
        REQUIRE(table.find(i) == table.end());
}

TEST_CASE( PREFIX "/rehash" )
{
    table_type table(-1, 32);

    /* initial fill */
    for (uint64_t i = 0; i != 10; ++i)
        table.insert(i, i);
    for (uint64_t i = 0; i != 10; ++i) {
        auto it = table.find(i);
        REQUIRE(it != table.end());
        REQUIRE((*it).key == i);
    }

    /* fill to trigger rehash(es) */
    for (uint64_t i = 10; i != 1000; ++i)
        table.insert(i, i);

    /* verify state */
    REQUIRE(table.size() == 1000);
    REQUIRE(table.capacity() >= 1000);
    for (uint64_t i = 0; i != 1000; ++i) {
        auto it = table.find(i);
        REQUIRE(it != table.end());
        REQUIRE((*it).key == i);
    }
}


#undef TABLE
#undef HASH
#undef STRINGIFY
#undef STR
#undef PREFIX
