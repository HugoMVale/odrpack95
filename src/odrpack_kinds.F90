module odrpack_kinds
!! Real kinds and common numeric constants.
    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none
    private

    public :: wp
    public :: ZERO, ONE, NEGONE
    
#ifdef REAL32
    integer, parameter :: wp = real32
#elif REAL64
    integer, parameter :: wp = real64
#else
    integer, parameter :: wp = real64
#endif
    
    real(wp), parameter :: ZERO = 0.0_wp
    real(wp), parameter :: ONE = 1.0_wp
    real(wp), parameter :: NEGONE = -1.0_wp
    
end module odrpack_kinds