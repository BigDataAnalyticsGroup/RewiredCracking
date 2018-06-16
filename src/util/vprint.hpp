#pragma once

#include <immintrin.h>
#include <iomanip>
#include <iostream>


template<typename value_type, typename vector_type>
void vprint(std::ostream &out, vector_type vec)
{
    const int digits = (2.4f * sizeof(value_type) + .5f);
    unsigned num_values = sizeof(vector_type) / sizeof(value_type);
    const value_type *values = reinterpret_cast<value_type*>(&vec);

    while (num_values--) {
        out << std::setw(digits) << *values++;
        if (num_values)
            out << ' ';
    }
}
