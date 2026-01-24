project(
  'odrpack95',
  ['fortran', 'c'],
  version: '3.0.1',
  meson_version: '>=1.1',
  license: 'MIT',
  default_options: [
  'buildtype=release',
  'fortran_std=none'
  ],
)

cc = meson.get_compiler('c')
fc = meson.get_compiler('fortran')
os = host_machine.system()

# Compiler flags
if fc.get_id() == 'gcc'
  add_project_arguments('-Wno-compare-reals', '-Wno-unused-dummy-argument', language: ['fortran'])
else
  warning('Compiler is @0@. If the build fails, try using the GCC compiler.'.format(fc.get_id()))
endif

# BLAS
blas_option = get_option('blas_option')
blas_dep = disabler()
blas_source = files()

if blas_option == 'refblas'
  message('`Reference BLAS` selected.')
  blas_source = files('extern/blas.f')

elif blas_option == 'openblas'
  message('`OpenBLAS` selected.')
  blas_dep = dependency('openblas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('openblas', required: false)
  endif
  if not blas_dep.found()
    error('OpenBLAS library not found. Install OpenBLAS or select another BLAS option.')
  endif

elif blas_option == 'auto'
  message('`Auto` BLAS selection enabled.')
  
  if os == 'darwin'
    blas_dep = dependency('appleframeworks', modules: ['Accelerate'], required: false)
    if blas_dep.found()
      message('Using Apple Accelerate framework (BLAS/LAPACK).')
    endif
  endif
  
  if not blas_dep.found()
    foreach lib : ['blas', 'openblas']
      blas_dep = dependency(lib, required: false)
      if not blas_dep.found()
        blas_dep = fc.find_library(lib, required: false)
      endif
      if blas_dep.found()
        message('Detected BLAS library: @0@'.format(lib))
        break
      endif
    endforeach
  endif

  if not blas_dep.found()
    warning('No BLAS library detected. Falling back to reference BLAS.')
    blas_source = files('extern/blas.f')
  endif

else
   error('Invalid BLAS option: @0@'.format(blas_option))
endif

# Sources
odrpack_source = files(
  'src/linpack.f',
  'src/blas_interfaces.f90',
  'src/odrpack.f90',
  'src/odrpack_kinds.f90',
  'src/odrpack_core.f90',
  'src/odrpack_reports.f90',
  'src/odrpack_capi.f90'
)

all_sources = odrpack_source + blas_source


# Library
alls_deps = []
if blas_dep.found()
  alls_deps += blas_dep
endif

if get_option('build_shared')
  odrpack_lib = shared_library(
    meson.project_name(),
    sources: all_sources,
    dependencies: alls_deps,
    install: true
  )
else
  odrpack_lib = static_library(
    meson.project_name(),
    sources: all_sources,
    dependencies: alls_deps
  )
endif

incdir = include_directories('c/include')
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
