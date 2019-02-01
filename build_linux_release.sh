#!/bin/sh

mkdir -p build/linux/release
cd build/linux/release
cmake -D CMAKE_BUILD_TYPE=Release ../../../
make