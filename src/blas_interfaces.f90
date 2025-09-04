module blas_interfaces
  !! Double precision interfaces for the BLAS procedures used by odrpack.

   use odrpack_kinds, only: dp
   implicit none

   interface

      pure real(dp) function dasum(n, x, incx)
      !! Sum of magnitudes of vector components.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function dasum

      pure subroutine daxpy(n, a, x, incx, y, incy)
      !! Computation `Y = A*X + Y`.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: a
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine daxpy

      pure subroutine dcopy(n, x, incx, y, incy)
      !! Vector copy `Y = X`.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine dcopy

      pure real(dp) function ddot(n, x, incx, y, incy)
      !! Inner product of vectors.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(in) :: y(*)
         integer, intent(in) :: incy
      end function ddot

      pure real(dp) function dnrm2(n, x, incx)
      !! Euclidean length (L2 Norm) of vector.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function dnrm2

      pure subroutine drot(n, x, incx, y, incy, c, s)
      !! Apply Givens rotation.
         import :: dp
         integer, parameter :: wp = dp
         integer, intent(in) :: n
         real(dp), intent(inout) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
         real(dp), intent(in) :: c
         real(dp), intent(in) :: s
      end subroutine drot

      pure subroutine drotg(a, b, c, s)
      !! Construct plane Givens rotation.
         import :: dp
         real(dp), intent(in) :: a
         real(dp), intent(in) :: b
         real(dp), intent(out) :: c
         real(dp), intent(out) :: s
      end subroutine drotg

      pure subroutine dscal(n, a, x, incx)
      !! Vector scale `X = A*X`.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: a
         real(dp), intent(inout) :: x(*)
         integer, intent(in) :: incx
      end subroutine dscal

      pure subroutine dswap(n, x, incx, y, incy)
      !! Interchange vectors.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine dswap

      pure integer function idamax(n, x, incx)
      !! Find largest component of vector.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function idamax

   end interface

end module blas_interfaces
