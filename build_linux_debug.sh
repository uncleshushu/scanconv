#!/bin/sh

mkdir -p build/linux/debug
cd build/linux/debug
cmake -D CMAKE_BUILD_TYPE=Debug ../../../
cmake --build .