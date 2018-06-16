#pragma once

#include <chrono>


struct Timer
{
    using clk        = std::chrono::high_resolution_clock;
    using time_point = clk::time_point;
    using duration   = clk::duration;

    Timer() : dur_(0) { }

    time_point start() { return tp_ = clk::now(); }

    void stop() {
        const auto t = clk::now();
        dur_ += (t - tp_);
        tp_ = t;
    }

    duration reset() {
        const auto d = dur_;
        dur_ = duration(0);
        return d;
    }

    duration operator()() { return dur_; }

    private:
    time_point tp_;
    duration dur_;
};
