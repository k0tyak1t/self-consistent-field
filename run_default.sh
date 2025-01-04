#!/bin/sh
mkdir -p build geometry basis
rm -f build/scf.out
g++ src/*.cpp -Iinclude -llapack -lboost_math_c99 -L/usr/lib/x86_64-linux-gnu -o build/scf.out -O3
time ./build/scf.out geometry/h2o.xyz basis/6311.basis
exit 0
