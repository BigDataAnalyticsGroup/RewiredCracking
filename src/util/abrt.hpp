#pragma once

#include <cstdio>
#include <cstdlib>


#define abrt(MSG) __abort(__FILE__, __LINE__, (MSG))

namespace
{

[[noreturn]] void __abort(const char *filename, unsigned line, const char *msg)
{
    fflush(stdout);
    fprintf(stderr, "%s:%u: %s\n", filename, line, msg);
    fflush(stderr);
    abort();
}

}
