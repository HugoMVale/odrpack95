# meson.build

project(
  'odrpack_c',
  'fortran', 'c',
  version: '1.0.0',
  meson_version: '>=0.55',
  default_options: ['default_library=static', 'buildtype=debugoptimized', 'fortran_std=none'],
)

odrpack_lic = files(
  '../LICENSE',
)

openblas = dependency('openblas', required: true)

# Sources
srcs = files(
  '../src/linpack.f',
  '../src/odrpack_kinds.F90',
  '../src/odrpack_core.f90',
  '../src/odrpack.f90',
  '../src/odrpack_reports.f90',
  'odrpack_capi.f90'
)

odrpack_c_lib = library(
  meson.project_name(),
  sources: srcs,
  dependencies: [openblas],
  install: false,
)

# Examples
subdir('example')

# # Tests
# subdir('test')
