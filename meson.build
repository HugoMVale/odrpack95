project(
  'odrpack95',
  ['fortran', 'c'],
  version: '2.0.1',
  meson_version: '>=1.1',
  license: 'MIT',
  default_options: [
  'buildtype=release',
  'fortran_std=none'
  ],
)

cc = meson.get_compiler('c')
os = host_machine.system()

# Compiler flags
if cc.get_id() == 'gcc'
  add_project_arguments('-Wno-compare-reals', '-Wno-unused-dummy-argument', language: ['fortran'])
  if cc.version().version_compare('>=14') and os != 'windows'
    add_project_arguments('-ftrampoline-impl=heap', language: ['c', 'fortran'])
  endif
else
  warning('Compiler is @0@. If the build fails, try using the GCC compiler.'.format(cc.get_id()))
endif

# BLAS
if os == 'darwin'
  message('Using `Accelerate` framework as BLAS library.')
  add_project_link_arguments('-framework', 'Accelerate')
  blas_dep = declare_dependency()
else
  # static linking does not work in cibuildwheel
  blas_dep = dependency('openblas', required: false)
  if not blas_dep.found()
    blas_dep = cc.find_library('openblas', required: false)
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
  'src/odrpack_capi.f90'
)

incdir = include_directories('c/include')

if get_option('build_shared')
  odrpack_lib = shared_library(
    meson.project_name(),
    sources: [odrpack_source, blas_source],
    dependencies: [blas_dep],
    install: true
  )
else
  odrpack_lib = static_library(
    meson.project_name(),
    sources: [odrpack_source, blas_source],
    dependencies: [blas_dep]
  )
endif

odrpack_dep = declare_dependency(
  include_directories: [incdir],
  link_with: [odrpack_lib]
)

# Examples
if get_option('build_examples')
    subdir('example')
    subdir('c')
endif

# Tests
if get_option('build_tests')
  subdir('test')
endif
