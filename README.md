# odrpack95

[![CI](https://github.com/HugoMVale/odrpack95/actions/workflows/CI.yml/badge.svg)](https://github.com/HugoMVale/odrpack95/actions)
[![codecov](https://codecov.io/gh/HugoMVale/odrpack95/branch/main/graph/badge.svg?token=1XL5LQSO9P)](https://codecov.io/gh/HugoMVale/odrpack95)
[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)

## Description

ODRPACK95 is a portable collection of Fortran 95 subprograms for
fitting a model to data with bound constraints on the model
parameters.  It is designed primarily for instances when the
explanatory as well as the response variables have significant
errors, implementing a highly efficient algorithm for solving the
weighted orthogonal distance regression problem, i.e., for minimizing
the sum of the squares of the weighted orthogonal distances between
each data point and the curve described by the model equation.  It
can also be used to solve the ordinary least squares problem where
all of the errors are attributed to the observations of the dependent
variable.

This project aims to modernize the original code, namely:

* [x] Modify the tests so they can be automatically run in the CI.
* [ ] Fix some warnings, e.g. "Fortran runtime warning: An array temporary was created".
* [ ] Implement a C API.
* [ ] Implement python bindings to the C API.
* [ ] Conversion from fixed (.f) to free-form (.f90) source.

## Installation

The easiest way to build/test the code and run the examples is by means of [`fpm`](https://fpm.fortran-lang.org/). 

To build the library, do:
```
fpm build --profile release
```

To run the tests, do:
```
fpm test --profile release
```
To run the provided example, do:
```
fpm run --example "simple_example"
```
 
## Licence

* The original ODERPACK95 code is [public domain](https://github.com/scipy/scipy/issues/7107#issuecomment-307378785).
* Any modications done in the course of this project are also public domain.

## References

* Original code from [Netlib](https://www.netlib.org/odrpack/).
* Paul T. Boggs, Janet R. Donaldson, Richaard h. Byrd, and Robert B. Schnabel. 1989. Algorithm 676: ODRPACK: software for weighted orthogonal distance regression. ACM Trans. Math. Softw. 15, 4 (Dec. 1989), 348–364. https://doi.org/10.1145/76909.76913
* Jason W. Zwolak, Paul T. Boggs, and Layne T. Watson. 2007. Algorithm 869: ODRPACK95: A weighted orthogonal distance regression code with bound constraints. ACM Trans. Math. Softw. 33, 4 (August 2007), 27–es. https://doi.org/10.1145/1268776.1268782
