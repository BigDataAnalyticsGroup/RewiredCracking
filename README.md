# An Analysis and Evaluation of Database Cracking Kernels

## Instructions

### Download the Source Code

```plain
git clone https://github.com/BigDataAnalyticsGroup/RewiredCracking.git
cd RewiredCracking
git submodule update --init --recursive
```
or
```plain
wget https://github.com/BigDataAnalyticsGroup/RewiredCracking/archive/master.zip
unzip master.zip
rm master.zip
cd RewiredCracking
wget https://github.com/opcm/pcm/archive/master.zip
rmdir processorcountermonitor
unzip master.zip
mv pcm-master processorcountermonitor
```

### Build the Project

```plain
mkdir build_debug
cd build_debug
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release ..
make -j 42
```

### Mount the `hugetlbfs`

```plain
sudo ./mount_hugetlbfs.sh
```

### Mount the `cpuset` pseudo filesystem

The following commands require privileged access.

```plain
mkdir /dev/cpuset
mount -t cpuset cpuset /dev/cpuset
```

### Execute the Benchmarks

Use `build_debug/bin/partition -h` for help.


```plain
build_debug/bin/partition 4 DENSE
```

You find the results in a CSV file that is named after the name of the host, e.g.
`partitioning_MY_AWESOME_HOSTNAME.csv`.

### Generate the Plots

The plot script is written in R.  Make sure to have all required libraries installed.

```plain
R -q --no-save --file=partitioning.r --args partitioning_MY_AWESOME_HOSTNAME.csv
```

### Explore the Results

The plots are named `partitioning_MY_AWESOME_HOSTNAME_<NAME>.pdf`.
