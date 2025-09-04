module example2_model
!! Model for example2.

   use iso_c_binding, only: c_ptr
   use odrpack_kinds, only: dp, one, zero
   implicit none

contains

   pure subroutine fcn( &
      n, m, q, np, ldifx, beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop, data)
   !! User-supplied subroutine for evaluating the model.

      integer, intent(in) :: n, m, q, np, ldifx, ideval, ifixb(np), ifixx(ldifx, m)
      real(dp), intent(in) :: beta(np), xplusd(n, m)
      real(dp), intent(out) :: f(n, q), fjacb(n, np, q), fjacd(n, m, q)
      integer, intent(out) :: istop
      type(c_ptr), intent(in), value :: data

      ! Local variables
      integer :: i

      ! Check for unacceptable values for this problem
      if (beta(1) > zero) then
         istop = 1
         return
      else
         istop = 0
      end if

      ! Compute predicted values
      if (mod(ideval, 10) > 0) then
         do i = 1, ubound(f, 2)
            f(:, i) = beta(3)*(xplusd(:, 1) - beta(1))**2 + &
                      2*beta(4)*(xplusd(:, 1) - beta(1))*(xplusd(:, 2) - beta(2)) + &
                      beta(5)*(xplusd(:, 2) - beta(2))**2 - one
         end do
      end if

   end subroutine fcn

end module example2_model

program example2
   !! Implicit ODR job.

   use odrpack_kinds, only: dp
   use odrpack, only: odr, odrpack_model
   use example2_model, only: fcn
   implicit none

   ! Variable declarations
   type(odrpack_model) :: model
   integer :: i, iprint, j, job, lundata, lunrpt, m, n, np, q
   real(dp), allocatable :: beta(:), x(:, :), y(:, :)

   ! Set model procedure
   model%fcn => fcn

   ! Set up report files
   open (newunit=lunrpt, file='./example/report2.dat')

   ! Read problem dimensions
   open (newunit=lundata, file='./example/data2.dat')
   read (lundata, fmt=*) n, m, np, q

   ! Allocate arrays
   allocate (beta(np), x(n, m), y(n, q))

   ! Read problem data
   read (lundata, fmt=*) (beta(i), i=1, np)
   do i = 1, n
      read (lundata, fmt=*) (x(i, j), j=1, m)
   end do
   close (lundata)

   ! Specify task: Implicit orthogonal distance regression
   !       With forward finite difference derivatives
   !       Covariance matrix constructed with recomputed derivatives
   !       DELTA initialized to zero
   !       Not a restart
   job = 00001
   iprint = 2002

   ! Compute solution
   call odr(model, n, m, q, np, beta, y, x, &
            job=job, lunerr=lunrpt, lunrpt=lunrpt, iprint=iprint)

   close (lunrpt)

end program example2
