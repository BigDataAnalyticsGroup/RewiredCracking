# An Analysis and Evaluation of Database Cracking Kernels

## Instructions

### Download the Source Code

```
git clone TODO
```
or
```
wget | unzip TODO
```

### Build the Project

```
mkdir build_debug
cd build_debug
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release ..
make -j 42
```

### Mount the `hugetlbfs`

```
sudo ./mount_hugetlbfs.sh
```

### Mount the `cpuset` pseudo filesystem

TODO

### Execute the Benchmarks

```
build_debug/bin/partition 100000 DENSE
```

You find the results in a CSV file that is named after the name of the host, e.g.
`partitioning_MY_AWESOME_HOSTNAME.csv`.

### Generate the Plots

The plot script is written in R.  Make sure to have all required libraries installed.

```bash
R -q --no-save --file=partitioning.r --args partitioning_MY_AWESOME_HOSTNAME.csv
```

### Explore the Results

The plots are named `partitioning_MY_AWESOME_HOSTNAME_<NAME>.pdf`.
