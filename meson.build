# meson.build

project('odrpack', 'fortran',
  version: '1.0.0',
  meson_version: '>=0.55',
  default_options: [
    'default_library=both',
    'buildtype=debugoptimized',
    'fortran_std=none'
  ],
)

odrpack_lic = files(
  'LICENSE',
)

# Sources
srcs = files(
  'src/odrpack_kinds.F90',
  'src/odrpack.f90',
  'src/lpkbls.f'
)

odrpack_lib = library(
  meson.project_name(),
  sources: srcs,
  install: false,
)

odrpack_dep = declare_dependency(
  link_with: odrpack_lib
)

# minpack_inc = minpack_lib.private_dir_include()
# minpack_dep = declare_dependency(
#   link_with: minpack_lib,
#   include_directories: [minpack_inc, include_directories('include')],
#   variables: {'includedir': meson.current_source_dir() / 'include'},
# )

odrpack_header = files(
  'c/include/odrpack/odrpack.h',
)

# Examples
# subdir('example')

# Tests
subdir('test')