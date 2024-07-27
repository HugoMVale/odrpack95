examples = ['1', '2', '3', '4', '5']

foreach example : examples
    executable(
      'example@0@'.format(example),
      sources: files('example@0@.f90'.format(example.underscorify())),
      link_with: odrpack_lib,
    )
endforeach
