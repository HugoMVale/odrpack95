module blas_interfaces
  !! Single and double precision interfaces for the BLAS procedures used by odrpack.
   use odrpack_kinds, only: sp, dp
   implicit none

   interface

      pure real(sp) function sasum(n, x, incx)
         import :: sp
         integer, intent(in) :: n
         real(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function sasum

      pure real(dp) function dasum(n, x, incx)
      !! Takes the sum of the absolute values.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function dasum

      pure subroutine saxpy(n, a, x, incx, y, incy)
         import :: sp
         integer, intent(in) :: n
         real(sp), intent(in) :: a
         real(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(sp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine saxpy

      pure subroutine daxpy(n, a, x, incx, y, incy)
      !! Constant times a vector plus a vector.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: a
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine daxpy

      pure subroutine scopy(n, x, incx, y, incy)
         import :: sp
         integer, intent(in) :: n
         real(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(sp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine scopy

      pure subroutine dcopy(n, x, incx, y, incy)
      !! Copies a vector, x, to a vector, y.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine dcopy

      pure real(sp) function sdot(n, x, incx, y, incy)
         import :: sp
         integer, intent(in) :: n
         real(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(sp), intent(in) :: y(*)
         integer, intent(in) :: incy
      end function sdot

      pure real(dp) function ddot(n, x, incx, y, incy)
      !! Forms the dot product of two vectors.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(in) :: y(*)
         integer, intent(in) :: incy
      end function ddot

      pure real(sp) function snrm2(n, x, incx)
         import :: sp
         integer, intent(in) :: n
         real(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function snrm2

      pure real(dp) function dnrm2(n, x, incx)
      !! Euclidean norm of a vector.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function dnrm2

      pure subroutine srot(n, x, incx, y, incy, c, s)
         import :: sp
         integer, intent(in) :: n
         real(sp), intent(inout) :: x(*)
         integer, intent(in) :: incx
         real(sp), intent(inout) :: y(*)
         integer, intent(in) :: incy
         real(sp), intent(in) :: c
         real(sp), intent(in) :: s
      end subroutine srot

      pure subroutine drot(n, x, incx, y, incy, c, s)
      !! Applies a plane rotation.
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

      pure subroutine srotg(a, b, c, s)
         import :: sp
         real(sp), intent(in) :: a
         real(sp), intent(in) :: b
         real(sp), intent(out) :: c
         real(sp), intent(out) :: s
      end subroutine srotg

      pure subroutine drotg(a, b, c, s)
         import :: dp
         real(dp), intent(in) :: a
         real(dp), intent(in) :: b
         real(dp), intent(out) :: c
         real(dp), intent(out) :: s
      end subroutine drotg

      pure subroutine sscal(n, a, x, incx)
         import :: sp
         integer, intent(in) :: n
         real(sp), intent(in) :: a
         real(sp), intent(inout) :: x(*)
         integer, intent(in) :: incx
      end subroutine sscal

      pure subroutine dscal(n, a, x, incx)
      !! Scales a vector by a constant.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: a
         real(dp), intent(inout) :: x(*)
         integer, intent(in) :: incx
      end subroutine dscal

      pure subroutine sswap(n, x, incx, y, incy)
         import :: sp
         integer, intent(in) :: n
         real(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(sp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine sswap

      pure subroutine dswap(n, x, incx, y, incy)
      !! Interchanges two vectors.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine dswap

      pure integer function isamax(n, x, incx)
         import :: sp
         integer, intent(in) :: n
         real(sp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function isamax

      pure integer function idamax(n, x, incx)
      !! Finds the index of the first element having maximum absolute value.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function idamax

   end interface

end module blas_interfaces
