#pragma once


template<typename T>
void minmax(T &a, T &b)
{
    using std::min;
    using std::max;
    T _a = min(a, b);
    T _b = max(a, b);
    a = _a; b = _b;
}
