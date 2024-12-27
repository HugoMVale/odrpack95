# meson.build

odrpack_c_lib = library(
  meson.project_name() + '_c',
  sources: 'odrpack_capi.f90',
  link_with: odrpack_lib,
  install: true,
)

# Examples
subdir('example')

# # Tests
# subdir('test')
