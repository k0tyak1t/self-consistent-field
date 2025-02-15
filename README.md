# SCF

<!--toc:start-->
- [SCF](#scf)
  - [Dependencies](#dependencies)
  - [How to use](#how-to-use)
    - [Geometry example](#geometry-example)
    - [Basis example](#basis-example)
<!--toc:end-->

## Dependencies

1) lapack
2) boost math c99

## How to use

1. Save geometry file in .xyz format

### Geometry example

```
3
H2O (blah blah, comment string)
H    0.0   0.0  0.0
H    3.0   0.0  0.0
O    1.5   0.5  0.0
```

2. Save basis in GAMESS-US format

### Basis example

```
HYDROGEN
S   4
1         1.301000E+01           1.968500E-02
2         1.962000E+00           1.379770E-01
3         4.446000E-01           4.781480E-01
4         1.220000E-01           5.012400E-01
S   1
1         1.220000E-01           1.000000E+00
P   1
1         7.270000E-01           1.0000000

OXYGEN
S   9
1         1.172000E+04           7.100000E-04
2         1.759000E+03           5.470000E-03
3         4.008000E+02           2.783700E-02
4         1.137000E+02           1.048000E-01
5         3.703000E+01           2.830620E-01
6         1.327000E+01           4.487190E-01
7         5.025000E+00           2.709520E-01
8         1.013000E+00           1.545800E-02
9         3.023000E-01          -2.585000E-03
S   9
1         1.172000E+04          -1.600000E-04
2         1.759000E+03          -1.263000E-03
3         4.008000E+02          -6.267000E-03
4         1.137000E+02          -2.571600E-02
5         3.703000E+01          -7.092400E-02
6         1.327000E+01          -1.654110E-01
7         5.025000E+00          -1.169550E-01
8         1.013000E+00           5.573680E-01
9         3.023000E-01           5.727590E-01
S   1
1         3.023000E-01           1.000000E+00
P   4
1         1.770000E+01           4.301800E-02
2         3.854000E+00           2.289130E-01
3         1.046000E+00           5.087280E-01
4         2.753000E-01           4.605310E-01
P   1
1         2.753000E-01           1.000000E+00
D   1
1         1.185000E+00           1.0000000
```

3. Build project using Makefile

```
# use DIIS
make -j4

# without DIIS
make -j4 DEFS=NDIIS
```

4. How to run

```
./scf.x path/to/geometry.xyz path/to/basis.basis
```
