#pragma once

#ifndef NDEBUG // -- Debug ---------------------------------------------------------------------------------------------

#include <cstdio>
#include <cstdlib>

#define _ASSERT1(P) _ASSERT2(P, nullptr)
#define _ASSERT2(P, M) __assert(__FILE__, __LINE__, (P), #P, (M))
#define GET_MACRO(_1, _2, MACRO, ...) MACRO
#define ASSERT(...) GET_MACRO(__VA_ARGS__, _ASSERT2, _ASSERT1, 0)(__VA_ARGS__)

inline void __assert(const char *file, const unsigned line, int pred, const char *predstr, const char *comment)
{
    if (pred) return;
    fprintf(stderr, "%s:%u: ASSERTION: (%s) is false", file, line, predstr);
    if (comment) fprintf(stderr, ", \"%s\"", comment);
    fputc('\n', stderr);
    fflush(stderr);
    std::abort();
}

inline void __assert(const char *file, const unsigned line, void *pred, const char *predstr, const char *comment)
{
    if (pred) return;
    fprintf(stderr, "%s:%u: ASSERTION: %s is NULL", file, line, predstr);
    if (comment) fprintf(stderr, ", \"%s\"", comment);
    fputc('\n', stderr);
    fflush(stderr);
    std::abort();
}

#else // -- Release ----------------------------------------------------------------------------------------------------

#define ASSERT(...)

#endif
