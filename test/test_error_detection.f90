module test_error_detection_m

   use odrpack_kinds, only: wp, zero, one
   implicit none

   real(wp), dimension(2) :: lower, upper

contains

   subroutine fcn( &
      n, m, np, nq, beta, xplusd, ifixb, ifixx, ldifx, ideval, f, fjacb, fjacd, istop)
   !! User-supplied subroutine for evaluating the model.

      integer, intent(in) :: ideval, ldifx, m, n, np, nq
      integer, intent(in) :: ifixb(np), ifixx(ldifx, m)
      real(kind=wp), intent(in) :: beta(np), xplusd(n, m)
      real(kind=wp), intent(out) :: f(n, nq), fjacb(n, np, nq), fjacd(n, m, nq)
      integer, intent(out) :: istop

      integer :: i

      ! Do something with FJACD, FJACB, IFIXB and IFIXX to avoid warnings that they
      ! are not being used.  This is simply not to worry users that the example
      ! program is failing.
      if (ifixb(1) > 0 .and. ifixx(1, 1) > 0 .and. fjacb(1, 1, 1) > zero &
         .and. fjacd(1, 1, 1) > zero) then
         ! Do nothing.
      end if

      if (any(lower > beta)) then
         write (0, *) "LOWER BOUNDS VIOLATED"
         do i = 1, np
            if (lower(i) > beta(i)) then
               write (0, *) "   IN THE ", i, " POSITION WITH ", beta(i), &
                  "<", lower(i)
            end if
         end do
      end if

      if (any(upper < beta)) then
         write (0, *) "UPPER BOUNDS VIOLATED"
         do i = 1, np
            if (upper(i) < beta(i)) then
               write (0, *) "   IN THE ", i, " POSITION WITH ", beta(i), &
                  ">", upper(i)
            end if
         end do
      end if

      istop = 0

      if (mod(ideval, 10) /= 0) then
         do i = 1, n
            f(i, 1) = beta(1)*exp(beta(2)*xplusd(i, 1))
         end do
      end if

      if (mod(ideval/10, 10) /= 0) then
         do i = 1, n
            fjacb(i, 1, 1) = exp(beta(2)*xplusd(i, 1))
            fjacb(i, 2, 1) = beta(1)*xplusd(i, 1)*exp(beta(2)*xplusd(i, 1))
         end do
      end if

      if (mod(ideval/100, 10) /= 0) then
         do i = 1, n
            fjacd(i, 1, 1) = beta(1)*beta(2)*exp(beta(2)*xplusd(i, 1))
         end do
      end if

   end subroutine fcn

end module test_error_detection_m

program test_error_detection
   !! Error detection tests for [[odrpack]].

   use odrpack_kinds, only: wp
   use odrpack, only: odr
   use test_error_detection_m, only: fcn, lower, upper
   implicit none

   integer :: n, m, nq, np, info, lunerr, lunrpt
   integer, allocatable :: iwork(:)
   real(wp), allocatable :: beta(:), y(:, :), x(:, :), delta(:, :), work(:)
   logical :: passed

   passed = .true.

   open (newunit=lunrpt, file="./test/test_error_detection_report.txt")
   open (newunit=lunerr, file="./test/test_error_detection_error.txt")

   ! Invalid dimensions
   n = 0
   m = 0
   nq = -1
   np = 0
   allocate (beta(np), y(n, nq), x(n, m))
   y = 0.0_wp
   x = 0.0_wp
   beta = 0.0_wp

   call odr(fcn, n, m, np, nq, beta, y, x, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 11111) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent dimensions of mandatory arrays X, Y, BETA
   n = 10
   m = 4
   nq = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np+1), y(n-1, nq+1), x(n+1, m-1))
   y = 0.0_wp
   x = 0.0_wp
   beta = 0.0_wp

   call odr(fcn, n, m, np, nq, beta, y, x, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 10**5) passed = .false.

   write (lunrpt, *) "INFO = ", info
   
   ! Inconsistent dimensions of optional WORK, IWORK
   n = 10
   m = 4
   nq = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, nq), x(n, m), iwork(42), work(69))
   y = 0.0_wp
   x = 0.0_wp
   beta = 0.0_wp

   call odr(fcn, n, m, np, nq, beta, y, x, iwork=iwork, work=work, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 10**5) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent dimensions of optional arrays DELTA, IFIXB, SCLB, STPB, IFIXX, SCLD, STPD
   n = 10
   m = 4
   nq = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, nq), x(n, m), delta(m, n))
   y = 0.0_wp
   x = 0.0_wp
   beta = 0.0_wp

   call odr(fcn, n, m, np, nq, beta, y, x, &
            delta=delta, ifixb=[1, 0], sclb=[1.0_wp, 1.0_wp], stpb=[1.0_wp, 1.0_wp], &
            ifixx=reshape([1, 0, 0, 1], [2, 2]), scld=reshape([1.0_wp, 1.0_wp], [2, 1]), &
            stpd=reshape([1.0_wp, 1.0_wp], [2, 1]), &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 10**5) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent dimensions of optional arrays WE, WD
   n = 10
   m = 4
   nq = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, nq), x(n, m))
   y = 0.0_wp
   x = 0.0_wp
   beta = 0.0_wp

   call odr(fcn, n, m, np, nq, beta, y, x, &
            wd=reshape([-1.0_wp], [1, 1, 1]), we=reshape([-1.0_wp], [1, 1, 1]), &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 10**5) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent LOWER and UPPER
   n = 10
   m = 4
   nq = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, nq), x(n, m))
   y = 0.0_wp
   x = 0.0_wp
   beta = 0.0_wp

   call odr(fcn, n, m, np, nq, beta, y, x, &
            lower=beta+1e-10_wp, upper=beta-1e-10_wp, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 91000) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! ERROR IN JOB SPECIFICATION WITH WORK AND IWORK
   n = 1
   m = 1
   nq = 1
   np = 1
   deallocate (beta, y, x)
   allocate (beta(np), y(n, nq), x(n, m))
   y = 0.0_wp
   x = 0.0_wp
   beta = 0.0_wp

   call odr(fcn, n, m, np, nq, beta, y, x, iprint=1, info=info, job=10000, &
            lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 70110) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! ERROR IN JOB SPECIFICATION WITH DELTA
   n = 1
   m = 1
   nq = 1
   np = 1
   deallocate (beta, y, x)
   allocate (beta(np), y(n, nq), x(n, m))
   y = 0.0_wp
   x = 0.0_wp
   beta = 0.0_wp

   call odr(fcn, n, m, np, nq, beta, y, x, iprint=1, info=info, job=1000, &
            lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 71000) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! BOUNDS TOO SMALL FOR DERIVATIVE CHECKER WHEN DERIVATIVES DON'T AGREE.
   n = 4
   m = 1
   nq = 1
   np = 2
   deallocate (beta, y, x)
   allocate (beta(np), y(n, nq), x(n, m))
   beta = [-200.0_wp, -5.0_wp]
   upper = [-200.0_wp, 0.0_wp]
   lower = [-200.000029802322_wp, -5.0_wp]
   y(:, 1) = [2.718281828459045_wp, 7.389056098930650_wp, &
              148.4131591025766_wp, 403.4287934927353_wp]
   x(:, 1) = [1.0_wp, 2.0_wp, 5.0_wp, 6.0_wp]

   call odr(fcn, n, m, np, nq, beta, y, x, iprint=1, info=info, job=0020, &
            lunrpt=lunrpt, lunerr=lunerr, lower=lower, upper=upper)

   if (info /= 1) passed = .false.

   write (lunrpt, *) "INFO = ", info

!  ERROR IN ARRAY ALLOCATION
!  The following code is intended to force memory allocation failure.  An
!  appropriate N for your machine must be chosen to ensure memory allocation
!  will fail within ODRPACK95.  A value of about 1/4 the total memory available
!  to a process should do the trick.  However, most modern operating systems and
!  Fortran compilers will not likely deny ODRPACK95 memory before they fail for
!  another reason.  Therefore, the memory allocation checks in ODRPACK95 are not
!  easy to provoke.  An operating system may return successfull memory
!  allocation but fail to guarantee the memory causing a segfault when some
!  memory locations are accessed.  A Fortran compiler or operating system may
!  allow limited sized stacks during subroutine invocation causing the ODRPACK95
!  call to fail before ODRPACK95 executes its first line.
!
!       N  = 032000000
!       M  = 1
!       NQ = 1
!       NP = 1
!       DEALLOCATE(BETA,Y,X)
!       ALLOCATE(BETA(NP),Y(N,NQ),X(N,M),STAT=STAT)
!       IF (STAT/=0) THEN
!       WRITE(0,*)
!       &       "SYSTEM ERROR: COULD NOT ALLOCATE MEMORY, TESTER ",
!       &       "FAILED TO RUN."
!       STOP
!       END IF
!       Y(:,:) = 0.0_wp
!       X(:,:) = 0.0_wp
!       BETA(:) = 0.0_wp
!
!       CALL ODR(FCN,N,M,NP,NQ,BETA,Y,X,IPRINT=1,INFO=INFO,
!       &   LUNRPT=LUN,LUNERR=LUN)
!
!       WRITE(LUN,*) "INFO = ", INFO
!
   close (lunrpt)
   close (lunerr)

   if (passed) then
      stop "Error detection tests passed. See report."
   else
      error stop "Error detection tests failed. See report."
   end if

end program test_error_detection
