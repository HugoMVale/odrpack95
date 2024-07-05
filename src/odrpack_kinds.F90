module odrpack_kinds
!! Real kinds and common numeric constants.
    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none
    private

    public :: wp
    public :: negone, zero, one, two, three, ten
    
#ifdef REAL32
    integer, parameter :: wp = real32
#elif REAL64
    integer, parameter :: wp = real64
#else
    integer, parameter :: wp = real64
#endif

    real(wp), parameter :: negone = -1.0_wp
    real(wp), parameter :: zero = 0.0_wp
    real(wp), parameter :: one = 1.0_wp
    real(wp), parameter :: two = 2.0_wp
    real(wp), parameter :: three = 3.0_wp
    real(wp), parameter :: ten = 10.0_wp

end module odrpack_kinds