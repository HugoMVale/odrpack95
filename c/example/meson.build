examples = ['3', '5']

foreach example : examples
    executable(
      'example@0@'.format(example),
      sources: files('example@0@.c'.format(example.underscorify())),
      link_with: odrpack_c_lib,
    )
endforeach
