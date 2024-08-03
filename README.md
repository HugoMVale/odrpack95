# odrpack

[![CI](https://github.com/HugoMVale/odrpack95/actions/workflows/CI.yml/badge.svg)](https://github.com/HugoMVale/odrpack95/actions)
[![codecov](https://codecov.io/gh/HugoMVale/odrpack95/branch/main/graph/badge.svg?token=1XL5LQSO9P)](https://codecov.io/gh/HugoMVale/odrpack95)
[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)

## Description

`odrpack` is a package for weighted orthogonal distance regression (ODR), also known as [errors-in-variables regression](https://en.wikipedia.org/wiki/Errors-in-variables_models). 
It is designed primarily for instances when both the explanatory and response variables have significant errors. 
The package implements a highly efficient algorithm for minimizing the sum of the squares of the weighted orthogonal
distances between each data point and the curve described by the model equation, subject to parameter bounds. The nonlinear
model can be either explicit or implicit. Additionally, `odrpack` can be used to solve the ordinary least squares problem where all of
the errors are attributed to the observations of the dependent variable.

<p align="center">
  <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/8/81/Total_least_squares.svg/220px-Total_least_squares.svg.png" width="200" alt="Deming regression; special case of ODR.">
</p>

## History

The first version of the library, named ODRPACK, was originally released in 1989 [1], and the
last "official" update, named ODRPACK95, dates from 2004 [2]. The full revision history is
provided on page iv of the [User's Reference Guide for ODRPACK95](https://github.com/HugoMVale/odrpack95/blob/main/original/Doc/guide.pdf).

`odrpack` is a modernization of the ODRPACK95 code [3], intented to make the library easier to
use and maintain. The main changes include:

* [x] Conversion from fixed-form (`.f`) to free-form (`.f90`).
* [x] Conversion from upper case to lower case.
* [x] Modularization.
* [x] Removal of `DATA` statements, labeled do loops, and (most) `goto`s.
* [x] Addition of `intent(in/out)` to all procedures.
* [x] Addition of explicit interfaces to BLAS routines.
* [x] Implementation of a C API.
* [x] Automatic code documentation with FORD.

|    Version    | Year |   Standard   |
|:-------------:|:----:|:------------:|
|  odrpack      | 2024 | Fortran 2018 |
| ODRPACK95 1.0 | 2004 |  Fortran 95  |
|  ODRPACK 2.0  | 1992 |  FORTRAN 77  |
|  ODRPACK 1.0  | 1989 |  FORTRAN 77  |

**References**

[1] Paul T. Boggs, Janet R. Donaldson, Richaard h. Byrd, and Robert B. Schnabel. 1989. Algorithm 676: ODRPACK: software for weighted orthogonal distance regression. ACM Trans. Math. Softw. 15, 4 (Dec. 1989), 348–364. https://doi.org/10.1145/76909.76913

[2] Jason W. Zwolak, Paul T. Boggs, and Layne T. Watson. 2007. Algorithm 869: ODRPACK95: A weighted orthogonal distance regression code with bound constraints. ACM Trans. Math. Softw. 33, 4 (August 2007), 27–es. https://doi.org/10.1145/1268776.1268782

[3] Original source code from [Netlib](https://www.netlib.org/odrpack/).

## Build instructions

### Dependencies

`odrpack` depends on a small number of functions from BLAS and LINPACK.

* The build configuration files provided with the code (see further below) assume
OpenBLAS is locally installed. If another BLAS source is preferred, the configuration files
should be adjusted accordingly. Alternatively, the subset of required BLAS functions is
available in [src/blas.f_](/src/blas.f_); it suffices to remove the underscore and edit the
build configuration files.
* The subset of required LINPACK functions is available in [src/linpack.f](/src/linpack.f) and
is already included in the build configuration files. No action is required.

### With fpm

The easiest way to build/test the code and run the examples is by means of [`fpm`](https://fpm.fortran-lang.org/).

To build the library, do:

```sh
fpm build --profile release
```

To run the tests, do:

```sh
fpm test --profile release
```

To run the provided examples, do:

```sh
fpm run --example "example_name" --profile release
```

### With meson

First, setup the build:

```sh
meson setup _build
```

To build the libraries (static and dynamic), do:

```sh
meson compile -C _build
```

To run the tests, do:

```sh
meson test -C _build
```

## Licence

* The original ODERPACK95 code is [public domain](https://github.com/scipy/scipy/issues/7107#issuecomment-307378785).
* Any modications done in the course of this project are covered by BSD-3.
