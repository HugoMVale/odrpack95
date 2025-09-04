module odrpack_kinds
!! Real kinds and common numeric constants.

   use iso_fortran_env, only: real64
   implicit none
   private

   public :: dp
   public :: negone, zero, half, one, two, three, eight, ten, fiftn, hundred
   public :: pi

   integer, parameter :: dp = real64

   real(dp), parameter :: negone = -1.0_dp
   real(dp), parameter :: zero = 0.0_dp
   real(dp), parameter :: half = 0.5_dp
   real(dp), parameter :: one = 1.0_dp
   real(dp), parameter :: two = 2.0_dp
   real(dp), parameter :: three = 3.0_dp
   real(dp), parameter :: eight = 8.0_dp
   real(dp), parameter :: ten = 10.0_dp
   real(dp), parameter :: fiftn = 15.0_dp
   real(dp), parameter :: hundred = 100.0_dp
   real(dp), parameter :: pi = 4*atan(one)

end module odrpack_kinds
