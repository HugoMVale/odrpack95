program example2
!! Implicit ODR job.
   use odrpack, only: odr
   use odrpack_kinds, only: wp
   implicit none

   ! Variable declarations
   integer :: i, info, j, job, lunerr, lunrpt, m, n, np, nq
   real(kind=wp), allocatable :: beta(:), x(:, :), y(:, :)
   external :: fcn

   ! Set up report files
   lunerr = 9
   lunrpt = 9
   open (unit=9, file='./example/report2.dat')

   ! Read problem dimensions
   open (unit=5, file='./example/data2.dat')
   read (5, fmt=*) n, m, np, nq

   ! Allocate arrays
   allocate (beta(np), x(n, m), y(n, nq))

   ! Read problem data
   read (5, fmt=*) (beta(i), i=1, np)
   do i = 1, n
      read (5, fmt=*) (x(i, j), j=1, m)
   end do
   close(5)

   ! Specify task: Implicit orthogonal distance regression
   !       With forward finite difference derivatives
   !       Covariance matrix constructed with recomputed derivatives
   !       DELTA initialized to zero
   !       Not a restart
   job = 00001

   ! Compute solution
   call odr(fcn=fcn, &
            n=n, m=m, np=np, nq=nq, &
            beta=beta, &
            y=y, x=x, &
            job=job, &
            lunerr=lunerr, lunrpt=lunrpt, &
            info=info)

end program example2

subroutine fcn(n, m, np, nq, ldn, ldm, ldnp, beta, xplusd, ifixb, ifixx, &
               ldifx, ideval, f, fjacb, fjacd, istop)

   use odrpack_kinds, only: wp, one, zero
   implicit none

   integer, intent(in) :: ideval, ldifx, ldm, ldn, ldnp, m, n, np, nq
   integer, intent(in) :: ifixb(np), ifixx(ldifx, m)
   real(kind=wp), intent(in) :: beta(np), xplusd(ldn, m)
   real(kind=wp), intent(out) :: f(ldn, nq), fjacb(ldn, ldnp, nq), fjacd(ldn, ldm, nq)
   integer, intent(out) :: istop

   ! Local variables
   integer :: i

   ! Check for unacceptable values for this problem
   if (beta(1) .gt. zero) then
      istop = 1
      return
   else
      istop = 0
   end if

   ! Compute predicted values
   if (mod(ideval, 10) .ge. 1) then
      do i = 1, nq
            f(:, i) = beta(3)*(xplusd(:, 1) - beta(1))**2 + &
                      2*beta(4)*(xplusd(:, 1) - beta(1))*(xplusd(:, 2) - beta(2)) + &
                      beta(5)*(xplusd(:, 2) - beta(2))**2 - one
      end do
   end if

end subroutine fcn
