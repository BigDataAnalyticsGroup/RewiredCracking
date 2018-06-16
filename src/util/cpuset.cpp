#include "util/cpuset.hpp"

#ifdef LINUX


#include "util/assert.hpp"
#include <cstdlib>
#include <cstring>
#include <err.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>


constexpr const char * CPUSET_DIR = "/dev/cpuset";

bool exists(const char *path)
{
    struct stat st;
    errno = 0;
    stat(path, &st);
    if (ENOENT == errno) return false;
    if (errno)
        err(EXIT_FAILURE, "error stat'ing '%s'", path);
    return true;
}

bool isdir(const char *path)
{
    struct stat st;
    errno = 0;
    stat(path, &st);
    if (ENOENT == errno) return false;
    if (errno)
        err(EXIT_FAILURE, "error stat'ing '%s'", path);
    return S_ISDIR(st.st_mode);
}

bool isfile(const char *path)
{
    struct stat st;
    errno = 0;
    stat(path, &st);
    if (ENOENT == errno) return false;
    if (errno) {
        warn("error stat'ing '%s'", path);
        std::terminate();
    }
    return S_ISREG(st.st_mode);
}

bool cpuset_exists() { return isdir(CPUSET_DIR); }

#define ASSERT_CPUSET_EXISTS \
    if (not cpuset_exists()) { \
        std::cerr << __FILE__ << ":" << __LINE__ << ": /dev/cpuset does not exist" << std::endl; \
        std::exit(EXIT_FAILURE); \
    }

void cpuset_create_partition(const char *name, const char *CPUs, const char *mem_banks)
{
    ASSERT_CPUSET_EXISTS;

    const std::string path = std::string{CPUSET_DIR} + '/' + name;
    const std::string p_cpus = path + "/cpuset.cpus";
    const std::string p_mems = path + "/cpuset.mems";

    if (isdir(path.c_str())) {
        std::cerr << path << " already exists" << std::endl;
        std::terminate();
    }

    if (mkdir(path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) {
        warn("error creating directory '%s'", path.c_str());
        std::terminate();
    }

    FILE *f_cpus = fopen(p_cpus.c_str(), "w");
    assert(f_cpus);
    fprintf(f_cpus, "%s", CPUs);
    fclose(f_cpus);

    FILE *f_mems = fopen(p_mems.c_str(), "w");
    assert(f_mems);
    fprintf(f_mems, "%s", mem_banks);
    fclose(f_mems);
}

void cpuset_remove_partition(const char *name)
{
    ASSERT_CPUSET_EXISTS;

    const std::string path = std::string{CPUSET_DIR} + '/' + name;
    const std::string p_tasks = path + "/tasks";
    const std::string p_root_tasks = std::string{CPUSET_DIR} + "/tasks";

    if (not isdir(path.c_str())) {
        std::cerr << path << " does not exist" << std::endl;
        std::terminate();
    }

    /* Migrate the tasks. */
    std::vector<uint32_t> PIDs;
    FILE *f_tasks = fopen(p_tasks.c_str(), "r");
    assert(f_tasks);
    while (not feof(f_tasks)) {
        uint32_t pid;
        fscanf(f_tasks, "%u", &pid);
        PIDs.push_back(pid);
    }
    fclose(f_tasks);

    FILE *f_root_tasks = fopen(p_root_tasks.c_str(), "a");
    for (uint32_t pid : PIDs) {
        assert(f_root_tasks);
        fprintf(f_root_tasks, "%u\n", pid);
        fflush(f_root_tasks);
    }
    fclose(f_root_tasks);

    rmdir(path.c_str());
}

std::size_t cpuset_move_all_tasks_to_partition(const char *name)
{
    ASSERT_CPUSET_EXISTS;

    const std::string p_root_tasks = std::string{CPUSET_DIR} + "/tasks";
    const std::string p_tasks = std::string{CPUSET_DIR} + '/' + name + "/tasks";

    std::vector<uint32_t> PIDs;
    FILE *f_root_tasks = fopen(p_root_tasks.c_str(), "r");
    assert(f_root_tasks);
    while (not feof(f_root_tasks)) {
        uint32_t pid;
        fscanf(f_root_tasks, "%u", &pid);
        PIDs.push_back(pid);
    }
    fclose(f_root_tasks);

    const std::size_t num_tasks_moved = PIDs.size();

    FILE *f_tasks = fopen(p_tasks.c_str(), "a");
    for (uint32_t pid : PIDs) {
        assert(f_tasks);
        fprintf(f_tasks, "%u\n", pid);
        fflush(f_tasks);
    }
    fclose(f_tasks);

    return num_tasks_moved;
}

void cpuset_move_task_to_partition(const char *name, uint32_t pid)
{
    ASSERT_CPUSET_EXISTS;

    const std::string p_tasks = std::string{CPUSET_DIR} + '/' + name + "/tasks";

    FILE *f_tasks = fopen(p_tasks.c_str(), "a");
    assert(f_tasks);
    fprintf(f_tasks, "%u\n", pid);
    fclose(f_tasks);
}

void cpuset_move_task_to_root(uint32_t pid)
{
    ASSERT_CPUSET_EXISTS;

    const std::string p_root_tasks = std::string{CPUSET_DIR} + "/tasks";

    FILE *f_root_tasks = fopen(p_root_tasks.c_str(), "a");
    assert(f_root_tasks);
    fprintf(f_root_tasks, "%u\n", pid);
    fclose(f_root_tasks);
}


#endif
