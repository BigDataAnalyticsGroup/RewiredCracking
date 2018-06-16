#pragma once

#include <cxxabi.h>
#include <typeinfo>


template<typename T>
const char * get_type_name()
{
    int status;
    const char *name = abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status);
    if (status != 0) name = typeid(T).name();
    return name;
}

