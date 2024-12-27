# meson.build

project(
  'odrpack',
  ['fortran', 'c'],
  version: '1.0.0',
  meson_version: '>=0.55',
  default_options: ['default_library=static', 'buildtype=debugoptimized', 'fortran_std=none'],
)

odrpack_lic = files(
  'LICENSE',
)

# Sources
srcs = files(
  'src/linpack.f',
  'src/blas_interfaces.f90',
  'src/odrpack_kinds.F90',
  'src/odrpack_core.f90',
  'src/odrpack.f90',
  'src/odrpack_reports.f90',
)

openblas = dependency('openblas', required: true)

odrpack_lib = library(
  meson.project_name(),
  sources: srcs,
  dependencies: [openblas],
  install: true
)

# C interface
subdir('c')

# Examples
subdir('example')

# Tests
subdir('test')
