module example3_model
!! Model for example3.

   use odrpack_kinds, only: wp, one, zero
   implicit none

contains

   pure subroutine fcn(beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop)
   !! User-supplied subroutine for evaluating the model.

      integer, intent(in) :: ideval, ifixb(:), ifixx(:, :)
      real(kind=wp), intent(in) :: beta(:), xplusd(:, :)
      real(kind=wp), intent(out) :: f(:, :), fjacb(:, :, :), fjacd(:, :, :)
      integer, intent(out) :: istop

      ! Local variables
      real(kind=wp) :: freq, omega, ctheta, stheta, theta, phi, r
      real(kind=wp), parameter :: pi = 4*atan(one)
      integer :: i

      ! Check for unacceptable values for this problem
      do i = 1, ubound(xplusd, 1)
         if (xplusd(i, 1) < zero) then
            istop = 1
            return
         end if
      end do
      istop = 0

      theta = pi*beta(4)*0.5E0_wp
      ctheta = cos(theta)
      stheta = sin(theta)

      ! Compute predicted values
      if (mod(ideval, 10) >= 1) then
         do i = 1, ubound(xplusd, 1)
            freq = xplusd(i, 1)
            omega = (2.0E0_wp*pi*freq*exp(-beta(3)))**beta(4)
            phi = atan2((omega*stheta), (1 + omega*ctheta))
            r = (beta(1) - beta(2))*sqrt((1 + omega*ctheta)**2 + (omega*stheta)**2)**(-beta(5))
            f(i, 1) = beta(2) + r*cos(beta(5)*phi)
            f(i, 2) = r*sin(beta(5)*phi)
         end do
      end if

   end subroutine fcn

end module example3_model

program example3
!! Explicit ODR job, with parameter bounds, weights, `delta` initialized by user, and central
!! difference derivatives.

   use odrpack_kinds, only: wp
   use odrpack, only: odr
   use example3_model, only: fcn
   implicit none

   ! Variable declarations
   integer :: i, iprint, j, job, lundata, lunrpt, m, n, np, q
   integer, allocatable :: ifixx(:, :)
   real(kind=wp), allocatable :: beta(:), x(:, :), y(:, :), wd(:, :, :), we(:, :, :), &
                                 delta(:, :)

   ! Set up report files
   open (newunit=lunrpt, file='./example/report3.dat')

   ! Read problem dimensions
   open (newunit=lundata, file='./example/data3.dat')
   read (lundata, fmt=*) n, m, np, q

   ! Allocate arrays
   allocate (beta(np), x(n, m), y(n, q), delta(n, m), we(n, q, q), wd(n, m, m), ifixx(n, m))

   ! Read problem data
   read (lundata, fmt=*) (beta(i), i=1, np)
   do i = 1, n
      read (lundata, fmt=*) (x(i, j), j=1, m), (y(i, j), j=1, q)
   end do
   close (lundata)

   ! Specify task as explicit orthogonal distance regression
   !       With central difference derivatives
   !       Covariance matrix constructed with recomputed derivatives
   !       `delta` initialized by user
   !       Not a restart
   ! And indicate long initial report
   !       No iteration reports
   !       Long final report
   job = 01010
   iprint = 2002

   ! Initialize `delta`, and specify first decade of frequencies as fixed
   do i = 1, n
      if (x(i, 1) < 100.0E0_wp) then
         delta(i, 1) = 0.0E0_wp
         ifixx(i, 1) = 0
      elseif (x(i, 1) <= 150.0E0_wp) then
         delta(i, 1) = 0.0E0_wp
         ifixx(i, 1) = 1
      elseif (x(i, 1) <= 1000.0E0_wp) then
         delta(i, 1) = 25.0E0_wp
         ifixx(i, 1) = 1
      elseif (x(i, 1) <= 10000.0E0_wp) then
         delta(i, 1) = 560.0E0_wp
         ifixx(i, 1) = 1
      elseif (x(i, 1) <= 100000.0E0_wp) then
         delta(i, 1) = 9500.0E0_wp
         ifixx(i, 1) = 1
      else
         delta(i, 1) = 144000.0E0_wp
         ifixx(i, 1) = 1
      end if
   end do

   ! Set weights
   do i = 1, n
      if (x(i, 1) == 100.0E0_wp .or. x(i, 1) == 150.0E0_wp) then
         we(i, 1, 1) = 0.0E0_wp
         we(i, 1, 2) = 0.0E0_wp
         we(i, 2, 1) = 0.0E0_wp
         we(i, 2, 2) = 0.0E0_wp
      else
         we(i, 1, 1) = 559.6E0_wp
         we(i, 1, 2) = -1634.0E0_wp
         we(i, 2, 1) = -1634.0E0_wp
         we(i, 2, 2) = 8397.0E0_wp
      end if
      wd(i, 1, 1) = (1.0E-4_wp)/(x(i, 1)**2)
   end do

   ! Compute solution
   call odr(fcn, n, m, q, np, beta, y, x, &
            delta=delta, we=we, wd=wd, ifixx=ifixx, &
            job=job, iprint=iprint, lunerr=lunrpt, lunrpt=lunrpt)

   close (lunrpt)

end program example3
