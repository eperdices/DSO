#!/bin/bash

echo "Ruilding ROS node"

cd ROS/DSO
mkdir build
cd build
cmake .. -DROS_BUILD_TYPE=Release
make -j4
