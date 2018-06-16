#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>


enum DISTRIBUTION
{
    DIST_Uniform_Sparse,
    DIST_Ascending,
    DIST_Uniform_Dense,
    DIST_Grid,
};

static const std::unordered_map<std::string, DISTRIBUTION> STR_TO_DISTRIBUTION = {
    {"SPARSE",      DIST_Uniform_Sparse},
    {"ASCENDING",   DIST_Ascending},
    {"DENSE",       DIST_Uniform_Dense},
    {"GRID",        DIST_Grid},
};
static const std::unordered_map<DISTRIBUTION, std::string> DISTRIBUTION_TO_STR = {
    {DIST_Uniform_Sparse,  "uniform sparse"},
    {DIST_Ascending,       "ascending"},
    {DIST_Uniform_Dense,   "uniform dense"},
    {DIST_Grid,            "grid"},
};

template<typename T>
void create_uniform_sparse(T *data, std::size_t size, const T sentinel, const T max = std::numeric_limits<T>::max());

template<typename T>
void create_ascending(T *data, std::size_t size, const T sentinel);

template<typename T>
void create_uniform_dense(T *data, std::size_t size, const T sentinel);

template<typename T>
void create_grid(T *data, std::size_t size, uint8_t bits_per_byte);

template<typename T>
void create_distribution(DISTRIBUTION dist, T *data, std::size_t size,
                         const T sentinel,
                         uint8_t bits_per_byte,
                         const T max = std::numeric_limits<T>::max())
{
    switch (dist) {
        case DIST_Uniform_Sparse:
            create_uniform_sparse(data, size, sentinel, max);
            return;

        case DIST_Ascending:
            create_ascending(data, size, sentinel);
            return;

        case DIST_Uniform_Dense:
            create_uniform_dense(data, size, sentinel);
            return;

        case DIST_Grid:
            create_grid(data, size, bits_per_byte);
            return;
    }
}
