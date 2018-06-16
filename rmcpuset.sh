#!/bin/bash

for l in $(cat /dev/cpuset/general/tasks);
do
    echo "$l" >> /dev/cpuset/tasks
done

rmdir /dev/cpuset/general/
rmdir /dev/cpuset/hpc/
