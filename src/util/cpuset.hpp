#pragma once

#include "util/system.hpp"
#include <cstdint>


#ifdef LINUX

/* Provide access to the CPUSET pseudo-filesystem to isolate CPUs for HPC benchmarks. */

/**
 * Creates a new CPUSET for the specified CPUs and memory banks with the given name.
 */
void cpuset_create_partition(const char *name, const char *CPUs, const char *mem_bank);

/**
 * Removes the specified CPUSET.  All tasks in that CPUSET are moved to the root CPUSET.
 */
void cpuset_remove_partition(const char *name);

/* Moves all tasks in the root CPUSET to the specified CPUSET.  Returns the number of tasks moved. */
std::size_t cpuset_move_all_tasks_to_partition(const char *name);

/**
 * Moves the task with the given PID from the root CPUSET to the specified CPUSET.
 */
void cpuset_move_task_to_partition(const char *name, uint32_t pid);

/**
 * Removes the task with the given PID from the CPUSET by moving it to the root CPUSET.
 */
void cpuset_move_task_to_root(uint32_t pid);

#endif
