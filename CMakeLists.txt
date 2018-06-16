cmake_minimum_required(VERSION 3.0)

project("linear_rehash_table" C CXX)

# Threads
find_package(Threads REQUIRED)

include(ProcessorCount)
ProcessorCount(J)

# Processor Counter Monitor (PCM)
add_custom_target(
    pcm
    ALL
    COMMAND make -j${J} lib
    BYPRODUCTS "${CMAKE_SOURCE_DIR}/processorcountermonitor/libPCM.a"
    WORKING_DIRECTORY ../processorcountermonitor
    COMMENT "Building Processor Counter Monitor (PCM) library..."
    USES_TERMINAL
    )

# Build flags
SET(EXECUTABLE_OUTPUT_PATH      "${PROJECT_BINARY_DIR}/bin")
SET(LIBRARY_OUTPUT_PATH         "${PROJECT_BINARY_DIR}/lib")

# Compile flags
execute_process(COMMAND getconf LEVEL1_DCACHE_LINESIZE
    OUTPUT_VARIABLE LEVEL1_DCACHE_LINESIZE
    OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND getconf PAGESIZE
    OUTPUT_VARIABLE PAGESIZE
    OUTPUT_STRIP_TRAILING_WHITESPACE)
SET(CMAKE_COMPILE_FLAGS         "-W -Wall -pedantic -DLEVEL1_DCACHE_LINESIZE=${LEVEL1_DCACHE_LINESIZE} -DPAGESIZE=${PAGESIZE} -march=native -ggdb -Wno-nested-anon-types -Wno-gnu-statement-expression")
SET(CMAKE_C_FLAGS               "${CMAKE_C_FLAGS} ${CMAKE_COMPILE_FLAGS}")
SET(CMAKE_CXX_FLAGS             "-std=c++1z ${CMAKE_CXX_FLAGS} ${CMAKE_COMPILE_FLAGS}")
SET(CMAKE_CXX_FLAGS_DEBUG       "${CMAKE_CXX_FLAGS_DEBUG} -DVERBOSE")

# Directories
include_directories(src/ processorcountermonitor/)
add_subdirectory(src/)
