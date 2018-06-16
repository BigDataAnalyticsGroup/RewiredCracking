#pragma once

#include <algorithm>
#include <iostream>
#include <limits>


template<typename T>
struct tuple_type
{
    using value_type = T;

    T key;
    T payload;

    tuple_type() = default;
    tuple_type(T key, T payload) : key(key), payload(payload) { }
    tuple_type(T key) : key(key), payload(0) { }

    /* Comparison */
    bool operator==(tuple_type other) const { return this->key == other.key; }
    bool operator!=(tuple_type other) const { return this->key != other.key; }
    bool operator< (tuple_type other) const { return this->key <  other.key; }
    bool operator> (tuple_type other) const { return this->key >  other.key; }
    bool operator<=(tuple_type other) const { return this->key <= other.key; }
    bool operator>=(tuple_type other) const { return this->key >= other.key; }

    /* Bitwise */
    tuple_type operator& (tuple_type other) const { return { key_type(this->key & other.key), payload_type(this->payload & other.payload) }; }
    tuple_type operator| (tuple_type other) const { return { key_type(this->key | other.key), payload_type(this->payload | other.payload) }; }
    tuple_type operator^ (tuple_type other) const { return { key_type(this->key ^ other.key), payload_type(this->payload ^ other.payload) }; }
    tuple_type operator~ () const { return { key_type(~key), payload_type(~payload)}; }

    tuple_type & operator+=(tuple_type other) {
        this->key += other.key;
        this->payload += other.payload;
        return *this;
    }
    tuple_type & operator-=(tuple_type other) {
        this->key -= other.key;
        this->payload -= other.payload;
        return *this;
    }

    tuple_type operator+(tuple_type other) {
        tuple_type res{*this};
        return res.operator+=(other);
    }
    tuple_type operator-(tuple_type other) {
        tuple_type res{*this};
        return res.operator-=(other);
    }

    friend std::ostream & operator<<(std::ostream &out, tuple_type t) {
        return out << '(' << t.key << ", " << t.payload << ')';
    }

    explicit operator long() const { return long(key); }
};

static_assert(sizeof(tuple_type<uint16_t>) ==  4, "unexpected padding");
static_assert(sizeof(tuple_type<uint32_t>) ==  8, "unexpected padding");
static_assert(sizeof(tuple_type<uint64_t>) == 16, "unexpected padding");

namespace std
{
    template<typename T>
    struct numeric_limits<tuple_type<T>> : numeric_limits<T> { };
}

