#include "util/system.hpp"


/*-- LINUX {{{--------------------------------------------------------------------------------------------------------*/
#ifdef LINUX

#include <sched.h>
#include <unistd.h>
#include <vector>

void set_cpu_affinity(unsigned cpuid)
{
    cpu_set_t cpus;
    CPU_ZERO(&cpus);
    CPU_SET(cpuid, &cpus);
    sched_setaffinity(getpid(), sizeof(cpus), &cpus);
}

void set_cpu_affinity(std::vector<unsigned> cpuids)
{
    cpu_set_t cpus;
    CPU_ZERO(&cpus);
    for (unsigned cpuid : cpuids)
        CPU_SET(cpuid, &cpus);
    sched_setaffinity(getpid(), sizeof(cpus), &cpus);
}

#endif
/*--}}}---------------------------------------------------------------------------------------------------------------*/
