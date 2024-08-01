module odrpack_kinds
!! Real kinds and common numeric constants.
    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none
    private

    public :: sp, dp, wp
    public :: negone, zero, half, one, two, three, eight, ten, fiftn, hundred
    public :: pi
    
    integer, parameter :: sp = real32
    integer, parameter :: dp = real64

#ifdef REAL32
    integer, parameter :: wp = sp
#elif REAL64
    integer, parameter :: wp = dp
#else
    integer, parameter :: wp = dp
#endif

    real(wp), parameter :: negone = -1.0_wp
    real(wp), parameter :: zero = 0.0_wp
    real(wp), parameter :: half = 0.5_wp
    real(wp), parameter :: one = 1.0_wp
    real(wp), parameter :: two = 2.0_wp
    real(wp), parameter :: three = 3.0_wp
    real(wp), parameter :: eight = 8.0_wp
    real(wp), parameter :: ten = 10.0_wp
    real(wp), parameter :: fiftn = 15.0_wp
    real(wp), parameter :: hundred = 100.0_wp
    real(wp), parameter :: pi = 4*atan(one)

end module odrpack_kinds