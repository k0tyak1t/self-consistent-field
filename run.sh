#!/bin/sh
mkdir -p compiled geometry basis
rm -f compiled/scf.out
g++ src/*.cpp -Iinclude -llapack -lboost_math_c99 -L/usr/lib/x86_64-linux-gnu -o compiled/scf.out -Wall -pedantic -O3
time ./compiled/scf.out geometry/h2o.xyz basis/6311.basis
exit 0
