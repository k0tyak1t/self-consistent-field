#!/bin/sh

mkdir -p build geometry basis

mol=$1
basis=$2
diis=$3

compile_flags="-Iinclude -DDEFAULT_MAX_ITER=500 -llapack -llapacke -lboost_math_c99 -L/usr/lib/x86_64-linux-gnu -O3"
if [ "$diis" = "-DNDIIS" ]; then
  compile_flags="${compile_flags} ${diis}"
fi

rm -f build/scf.out

g++ src/*.cpp ${compile_flags} -o build/scf.out

time ./build/scf.out geometry/${mol}.xyz basis/${basis}.basis

exit 0
