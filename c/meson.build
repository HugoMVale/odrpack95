odrpack_c_lib = library(
  meson.project_name() + '_c',
  sources: 'odrpack_capi.f90',
  link_with: odrpack_lib,
  install: true,
)

incdir = include_directories('include')

odrpack_c_dep = declare_dependency(
  include_directories: [incdir],
  link_with: [odrpack_c_lib]
)

if get_option('build_examples')
  subdir('example')
endif

