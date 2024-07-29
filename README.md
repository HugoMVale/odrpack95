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
  <img src="https://en.wikipedia.org/wiki/Deming_regression#/media/File:Total_least_squares.svg" width="500" alt="Deming regression; special case of ODR.">
</p>

## Revision history

|    Version    | Year |   Standard   |
|:-------------:|:----:|:------------:|
|  odrpack      | 2024 | Fortran 2018 |
| ODRPACK95 1.0 | 2004 |  Fortran 95  |
|  ODRPACK 2.0  | 1992 |  FORTRAN 77  |
|  ODRPACK 1.0  | 1989 |  FORTRAN 77  |

The detailed revision history of ODRPACK and ODRPACK95 is provided on page iv of the 
[User's Reference Guide for ODRPACK95](https://github.com/HugoMVale/odrpack95/blob/main/original/Doc/guide.pdf).

### odrpack vs ODRPACK95

This project aims to modernize the ODRPACK95 code, namely:

* [x] Modify the tests so they can be automatically run in the CI.
* [x] Convert from fixed-form (`.f`) to free-form (`.f90`).
* [x] Convert from upper case to lower case.
* [x] Split the code in modules.
* [x] Add `intent(in/out)` to all procedures.
* [x] Remove labeled do loops, and (almost) all gotos.
* [x] Implement a C API.
* [x] Generate automatic code documentation with FORD.
* [ ] Implement python bindings to the C API.
* [ ] Replace LINPACK calls by LAPACK. Code profiling shows that LINPACK calls are not time
      consuming, so this change is motivated by code maintainability, not by performance.

## Build instructions

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

To run the provided example, do:

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

## References

* Original code from [Netlib](https://www.netlib.org/odrpack/).
* Paul T. Boggs, Janet R. Donaldson, Richaard h. Byrd, and Robert B. Schnabel. 1989. Algorithm 676: ODRPACK: software for weighted orthogonal distance regression. ACM Trans. Math. Softw. 15, 4 (Dec. 1989), 348–364. https://doi.org/10.1145/76909.76913
* Jason W. Zwolak, Paul T. Boggs, and Layne T. Watson. 2007. Algorithm 869: ODRPACK95: A weighted orthogonal distance regression code with bound constraints. ACM Trans. Math. Softw. 33, 4 (August 2007), 27–es. https://doi.org/10.1145/1268776.1268782
