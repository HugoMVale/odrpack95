examples = ['3', '5']

foreach example : examples
    executable(
      'example@0@'.format(example),
      sources: files('example@0@.c'.format(example.underscorify())),
      link_with: odrpack_lib,
      link_args: host_machine.system() == 'windows' ? [] : ['-lm'],
      link_language: 'fortran',
    )
endforeach
