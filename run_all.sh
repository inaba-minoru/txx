#!/usr/bin/env bash

# If project not ready, generate cmake file.
if [[ ! -d build ]]; then
    mkdir -p build
    cd build
    cmake ..
    cd ..
fi

# Build project.
cd build
make -j
cd ..

# Run all testcases. 
# You can comment some lines to disable the run of specific examples.
mkdir -p output
# time bin/PA1 testcases/scene09_s.txt output/scene09_s.bmp


# time bin/PA1 mycase/case1.txt output/scene.bmp
time bin/PA1 mycase/case2.txt output/scene.bmp
# time bin/PA1 mycase/case3.txt output/scene.bmp

# time bin/PA1 testcases/scene.txt output/scene.bmp
# time bin/PA1 testcases/scene0X.txt output/scene0X.bmp
# time bin/PA1 testcases/scene00.txt output/scene00.bmp
# time bin/PA1 testcases/scene10_wineglass.txt output/scene10_wineglass.bmp
# time bin/PA1 testcases/scene06_bunny_1k.txt output/scene06.bmp
# time bin/PA1 testcases/scene01_basic.txt output/scene01.bmp
# time bin/PA1 testcases/scene02_cube.txt output/scene02.bmp
# time bin/PA1 testcases/scene03_sphere.txt output/scene03.bmp
# time bin/PA1 testcases/scene04_axes.txt output/scene04.bmp
# time bin/PA1 testcases/scene05_bunny_200.txt output/scene05.bmp
# time bin/PA1 testcases/scene07_shine.txt output/scene07.bmp
# time bin/PA1 testcases/scene08_dragon.txt output/scene08.bmp
