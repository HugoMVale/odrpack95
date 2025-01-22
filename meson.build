project(
  'odrpack95',
  ['fortran', 'c'],
  version: '1.0.1',
  meson_version: '>=1.1',
  license: 'MIT',
  default_options: [
  'default_library=static',
  'buildtype=release',
  'fortran_std=none'
  ],
)

# fortran_compiler = meson.get_compiler('fortran')
c_compiler = meson.get_compiler('c')

# BLAS
if host_machine.system() == 'darwin'
  message('Using `Accelerate` framework as BLAS library.')
  add_project_link_arguments('-framework', 'Accelerate', language: ['c', 'fortran'])
  blas_dep = declare_dependency()
else
  blas_dep = dependency('blas', required: false)
  if not blas_dep.found()
    blas_dep = c_compiler.find_library('blas', required: false)
  endif
endif

if blas_dep.found()
  blas_source = []
else
  message('No specific BLAS library detected. Falling back to reference BLAS.')
  blas_source = join_paths(meson.current_source_dir(), 'extern/blas.f')
endif

# Sources
odrpack_source = files(
  'src/linpack.f',
  'src/blas_interfaces.f90',
  'src/odrpack.f90',
  'src/odrpack_kinds.F90',
  'src/odrpack_core.f90',
  'src/odrpack_reports.f90',
)

odrpack_lib = library(
  meson.project_name(),
  sources: [odrpack_source, blas_source],
  dependencies: [blas_dep],
  install: true
)

# C interface
subdir('c')

# Examples
if get_option('build_examples')
    subdir('example')
endif

# Tests
subdir('test')
