#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <sys/mman.h>


struct default_allocator
{
    template<typename T>
    T * allocate(std::size_t n) const { return static_cast<T*>(malloc(sizeof(T) * n)); }

    template<typename T>
    void deallocate(T *p, std::size_t n) const { (void) n; free(p); }

    template<typename T>
    T * reallocate(T *p, std::size_t old_size, std::size_t new_size) const {
        (void) old_size;
        return static_cast<T*>(realloc(p, sizeof(T) * new_size));
    }
};

template<std::size_t Align>
struct aligned_allocator
{
    static_assert((Align & (Align - 1)) == 0, "not a power of 2");

    template<typename T>
    T * allocate(std::size_t n) const { return static_cast<T*>(aligned_alloc(Align, sizeof(T) * n)); }

    template<typename T>
    void deallocate(T *p, std::size_t n) const { (void) n; free(p); }

    template<typename T>
    T * reallocate(T *p, std::size_t old_size, std::size_t new_size) const {
        T *dst = allocate<T>(new_size);
        std::memcpy(dst, p, sizeof(T) * old_size);
        deallocate(p, old_size);
        return dst;
    }
};

struct mremap_allocator
{
    template<typename T>
    T * allocate(std::size_t n) const {
        return static_cast<T*>(mmap(nullptr, sizeof(T) * n, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
    }

    template<typename T>
    void deallocate(T *p, std::size_t n) const { munmap(p, sizeof(T) * n); }

    template<typename T>
    T * reallocate(T *p, std::size_t old_size, std::size_t new_size) const {
        return static_cast<T*>(mremap(p, sizeof(T) * old_size, sizeof(T) * new_size, MREMAP_MAYMOVE));
    }
};
