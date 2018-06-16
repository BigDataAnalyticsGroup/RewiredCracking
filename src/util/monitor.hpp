#pragma once

#include <algorithm>
#include <chrono>
#include <ctime>
#include <iostream>
#include <papi.h>
#include <unordered_map>
#include <vector>
#include <x86intrin.h>


/** Return the value of the Time Stamp Counter (TSC) register. */
inline uint64_t rdtsc()
{
    return __rdtsc();
}

/** Return the value of the Time Stamp Counter (TSC) register after finishing all previous instructions. */
inline uint64_t rdtscp()
{
    uint32_t aux;
    return __rdtscp(&aux);
}

struct EventSet
{
    EventSet(std::initializer_list<int> papi_events)
        : events_(papi_events)
    {
        for (std::size_t i = 0, e = events_.size(); i != e; ++i)
            index_.emplace(events_[i], i);
    }

    const int * get() const { return &*events_.begin(); }
    std::size_t operator[](int papi_event) const { return index_.at(papi_event); }
    std::size_t size() const { return events_.size(); }

    private:
    std::vector<int> events_;
    std::unordered_map<int, std::size_t> index_;
};

/** Implements a simple interface to PAPI. */
struct Monitor
{
    Monitor(const EventSet &e)
        : e(e)
        , values(new long long[e.size()]())
    { }

    ~Monitor() {
        delete[] values;
    }

    void start() {
        asm volatile ("" : : : "memory"); // enforce instruction ordering
        int ret = PAPI_start_counters((int*) e.get(), e.size());
        if (ret != PAPI_OK) {
            std::cerr << "Error programming counters: " << PAPI_strerror(ret) << std::endl;
            std::terminate();
        }
        tsc_start = rdtsc();
    }

    void stop() {
        tsc_stop = rdtscp();
        PAPI_stop_counters(values, e.size());
    }

    long long operator[](int papi_event) const { return values[e[papi_event]]; }
    uint64_t getTSC() const { return tsc_stop - tsc_start; }

    const EventSet &e;
    long long *values;
    uint64_t tsc_start, tsc_stop;
};
