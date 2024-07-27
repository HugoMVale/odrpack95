# meson.build

project(
  'odrpack',
  'fortran',
  version: '1.0.0',
  meson_version: '>=0.55',
  default_options: ['default_library=static', 'buildtype=debugoptimized', 'fortran_std=none'],
)

odrpack_lic = files(
  'LICENSE',
)

# Sources
srcs = files(
  'src/blas.f',
  'src/linpack.f',
  'src/odrpack_kinds.F90',
  'src/odrpack_core.f90',
  'src/odrpack.f90',
  'src/odrpack_reports.f90',
)

odrpack_lib = library(
  meson.project_name(),
  sources: srcs,
  install: false,
)

# Examples
subdir('example')

# Tests
subdir('test')
