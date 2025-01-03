---
project: odrpack
license: mit
summary: odrpack is a package for weighted orthogonal distance regression (ODR), also known as errors-in-variables regression.
src_dir: src
         example
         c
exclude: src/linpack.f
output_dir: _site
page_dir: doc
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            iso_c_binding:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html#ISO_005fC_005fBINDING
source: true
proc_internals: true
graph: true
coloured_edges: true
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
project_github: https://github.com/HugoMVale/odrpack95
author: HugoMVale
github: https://github.com/HugoMVale/
email: 57530119+HugoMVale@users.noreply.github.com
dbg: true
predocmark: >
docmark_alt: #
predocmark_alt: <
md_extensions: markdown.extensions.toc
---

`odrpack` is a package for weighted orthogonal distance regression (ODR), also known as [errors-in-variables regression](https://en.wikipedia.org/wiki/Errors-in-variables_models). 
It is designed primarily for instances when both the explanatory and response variables have significant errors. 
The package implements a highly efficient algorithm for minimizing the sum of the squares of the weighted orthogonal
distances between each data point and the curve described by the model equation, subject to parameter bounds. The nonlinear
model can be either explicit or implicit. Additionally, `odrpack` can be used to solve the ordinary least squares problem where all of
the errors are attributed to the observations of the dependent variable.

<p align="center">
  <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/8/81/Total_least_squares.svg/220px-Total_least_squares.svg.png" width="200" alt="Deming regression; special case of ODR.">
</p>
