odrpack_c_lib = library(
  meson.project_name() + '_c',
  sources: 'odrpack_capi.f90',
  link_with: odrpack_lib,
  install: true,
)

if get_option('build_examples')
  subdir('example')
endif

