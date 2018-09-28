#!/bin/bash

echo "Configuring and building DSO ..."

mkdir build
cd build
cmake ..
make -j4
