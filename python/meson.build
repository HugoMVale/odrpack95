project(
  'odrpack',
  ['fortran', 'c', 'cpp'],
  version : '0.1.0',
  meson_version: '>=1.5.0',
  default_options: [
    'buildtype=debugoptimized',
    'b_ndebug=if-release',
    'c_std=c17',
    'cpp_std=c++17'
  ],
)

# meson_path = meson.current_build_dir()
py = import('python').find_installation(pure: false)
py_dep = dependency('python3', required: true)
pybind11_dep = dependency('pybind11', required: true)
openblas_dep = dependency('openblas', required: true)

# odrpack fortran
odrpack_source = files(
  '../src/linpack.f',
  '../src/blas_interfaces.f90',
  '../src/odrpack_kinds.F90',
  '../src/odrpack_core.f90',
  '../src/odrpack.f90',
  '../src/odrpack_reports.f90',
  '../c/odrpack_capi.f90'
)

odrpack_lib = library(
  meson.project_name(),
  sources: odrpack_source,
  dependencies: [openblas_dep],
  install: true
)

# Extension fortran module
py.extension_module(
 '_odrpack',
'_odrpack.py.cpp',
 dependencies : [py_dep, pybind11_dep],
 link_with: odrpack_lib,
 include_directories: include_directories('../c/include/odrpack'),
 install : true,
)
