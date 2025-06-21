module example2_model
!! Model for example2.

   use odrpack_kinds, only: wp, one, zero
   implicit none

contains

   pure subroutine fcn( &
      n, m, q, np, beta, xplusd, ifixb, ifixx, ldifx, ideval, f, fjacb, fjacd, istop)
   !! User-supplied subroutine for evaluating the model.

      integer, intent(in) :: ideval, ldifx, m, n, np, q
      integer, intent(in) :: ifixb(np), ifixx(ldifx, m)
      real(kind=wp), intent(in) :: beta(np), xplusd(n, m)
      real(kind=wp), intent(out) :: f(n, q), fjacb(n, np, q), fjacd(n, m, q)
      integer, intent(out) :: istop

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
      if (mod(ideval, 10) >= 1) then
         do i = 1, q
            f(:, i) = beta(3)*(xplusd(:, 1) - beta(1))**2 + &
                      2*beta(4)*(xplusd(:, 1) - beta(1))*(xplusd(:, 2) - beta(2)) + &
                      beta(5)*(xplusd(:, 2) - beta(2))**2 - one
         end do
      end if

   end subroutine fcn

end module example2_model

program example2
   !! Implicit ODR job.

   use odrpack, only: odr
   use odrpack_kinds, only: wp
   use example2_model, only: fcn
   implicit none

   ! Variable declarations
   integer :: i, info, iprint, j, job, lundata, lunrpt, m, n, np, q
   real(kind=wp), allocatable :: beta(:), x(:, :), y(:, :)

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
   call odr(fcn=fcn, &
            n=n, m=m, q=q, np=np, &
            beta=beta, &
            y=y, x=x, &
            job=job, &
            lunerr=lunrpt, lunrpt=lunrpt, iprint=iprint, &
            info=info)

   close (lunrpt)

end program example2
