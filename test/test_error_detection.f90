module test_error_detection_m

   use iso_c_binding, only: c_ptr
   use odrpack_kinds, only: dp, one, zero
   implicit none

   real(dp), allocatable :: lower(:), upper(:)

contains

   subroutine fcn( &
      n, m, q, np, ldifx, beta, xplusd, ifixb, ifixx, ideval, f, fjacb, fjacd, istop, data)
   !! User-supplied subroutine for evaluating the model.

      integer, intent(in) :: n, m, q, np, ldifx, ideval, ifixb(np), ifixx(ldifx, m)
      real(dp), intent(in) :: beta(np), xplusd(n, m)
      real(dp), intent(out) :: f(n, q), fjacb(n, np, q), fjacd(n, m, q)
      integer, intent(out) :: istop
      type(c_ptr), intent(in), value :: data

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
         do i = 1, size(beta)
            if (lower(i) > beta(i)) then
               write (0, *) "   IN THE ", i, " POSITION WITH ", beta(i), "<", lower(i)
            end if
         end do
      end if

      if (any(upper < beta)) then
         write (0, *) "UPPER BOUNDS VIOLATED"
         do i = 1, size(beta)
            if (upper(i) < beta(i)) then
               write (0, *) "   IN THE ", i, " POSITION WITH ", beta(i), ">", upper(i)
            end if
         end do
      end if

      istop = 0

      if (mod(ideval, 10) /= 0) then
         do i = 1, ubound(f, 1)
            f(i, 1) = beta(1)*exp(beta(2)*xplusd(i, 1))
         end do
      end if

      if (mod(ideval/10, 10) /= 0) then
         do i = 1, ubound(f, 1)
            fjacb(i, 1, 1) = exp(beta(2)*xplusd(i, 1))
            fjacb(i, 2, 1) = beta(1)*xplusd(i, 1)*exp(beta(2)*xplusd(i, 1))
         end do
      end if

      if (mod(ideval/100, 10) /= 0) then
         do i = 1, ubound(f, 1)
            fjacd(i, 1, 1) = beta(1)*beta(2)*exp(beta(2)*xplusd(i, 1))
         end do
      end if

   end subroutine fcn

end module test_error_detection_m

program test_error_detection
   !! Error detection tests for [[odrpack]].

   use odrpack_kinds, only: dp
   use odrpack, only: odr, odrpack_model
   use test_error_detection_m, only: fcn, lower, upper
   implicit none

   type(odrpack_model) :: model
   integer :: n, m, q, np, info, lunerr, lunrpt
   integer, allocatable :: iwork(:)
   real(dp), allocatable :: beta(:), y(:, :), x(:, :), delta(:, :), rwork(:), &
                            we(:, :, :), wd(:, :, :)
   logical :: passed

   passed = .true.

   open (newunit=lunrpt, file="./test/test_error_detection_report.txt")
   open (newunit=lunerr, file="./test/test_error_detection_error.txt")

   ! Invalid problem dimensions
   n = 0
   m = 0
   q = -1
   np = 0
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   ! Set model procedure
   model%fcn => fcn

   call odr(model, n, m, q, np, beta, y, x, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 11111) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent dimensions of mandatory arrays X, Y, BETA
   n = 10
   m = 4
   q = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np + 1), y(n - 1, q + 1), x(n + 1, m - 1))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 10**5) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent dimensions of optional RWORK, IWORK
   n = 10
   m = 4
   q = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m), iwork(42), rwork(69))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, iwork=iwork, rwork=rwork, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 10**5) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent dimensions of optional arrays DELTA, IFIXB, SCLB, STPB, IFIXX, SCLD, STPD
   n = 10
   m = 4
   q = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m), delta(m, n))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            delta=delta, ifixb=[1, 0], sclb=[1.0_dp, 1.0_dp], stpb=[1.0_dp, 1.0_dp], &
            ifixx=reshape([1, 0, 0, 1], [2, 2]), scld=reshape([1.0_dp, 1.0_dp], [2, 1]), &
            stpd=reshape([1.0_dp, 1.0_dp], [2, 1]), &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 10**5) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent dimensions of optional arrays WE, WD
   n = 10
   m = 4
   q = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            wd=reshape([-1.0_dp], [1, 1, 1]), we=reshape([-1.0_dp], [1, 1, 1]), &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 10**5) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent dimensions of optional arrays LOWER and UPPER
   n = 10
   m = 4
   q = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            lower=[1.0_dp, 1.0_dp], upper=[1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 10**5) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional arrays LOWER > UPPER
   n = 10
   m = 4
   q = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            lower=beta + 1e-10_dp, upper=beta - 1e-10_dp, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 91000) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional arrays BETA < LOWER
   n = 10
   m = 4
   q = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            lower=beta + 1e-10_dp, upper=beta + 1.0_dp, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 90100) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional arrays UPPER - LOWER ~ 0
   n = 10
   m = 4
   q = 2
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m), lower(np), upper(np))
   y = 0.0_dp
   x = 0.0_dp
   beta = 1.0_dp
   lower = beta - 1e-100_dp
   upper = beta + 1e-100_dp

   call odr(model, n, m, q, np, beta, y, x, &
            lower=lower, upper=upper, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 90010) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional arrays UPPER - LOWER ~ 0

   lower = beta - 1e-10_dp
   upper = beta + 1e-10_dp

   call odr(model, n, m, q, np, beta, y, x, &
            lower=lower, upper=upper, stpb=beta*1e-5_dp, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 90020) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array SCLB
   n = 10
   m = 1
   q = 1
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            sclb=[1.0_dp, -1.0_dp, 1.0_dp], &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 30100) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array SCLD(n,m)
   n = 3
   m = 1
   q = 1
   np = 2
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            scld=reshape([1.0_dp, -1.0_dp, 1.0_dp], [n, m]), &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 30200) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional arrays SCLD(1,m)
   n = 3
   m = 2
   q = 1
   np = 2
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            scld=reshape([1.0_dp, -1.0_dp], [1, m]), &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 30200) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array STPB
   n = 10
   m = 1
   q = 1
   np = 3
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            stpb=[1.0_dp, -1.0_dp, 1.0_dp], &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 31000) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array STPD(n,m)
   n = 3
   m = 1
   q = 1
   np = 2
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            stpd=reshape([1.0_dp, -1.0_dp, 1.0_dp], [n, m]), &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 32000) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array STPD(1,m)
   n = 3
   m = 2
   q = 1
   np = 2
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            stpd=reshape([1.0_dp, -1.0_dp], [1, m]), &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 32000) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array WE(1,1,q)
   n = 5
   m = 2
   q = 3
   np = 2
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m), we(1, 1, q))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp
   we = -1.0_dp
   we(1, 1, 1) = 1.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            we=we, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 30010) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array WE(1,q,q)
   n = 5
   m = 2
   q = 3
   np = 2
   deallocate (beta, y, x, we)
   allocate (beta(np), y(n, q), x(n, m), we(1, q, q))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp
   we = -1.0_dp
   we(1, 1, 1) = 1.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            we=we, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 30010) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array WE(n,1,q)
   n = 5
   m = 2
   q = 3
   np = 2
   deallocate (beta, y, x, we)
   allocate (beta(np), y(n, q), x(n, m), we(n, 1, q))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp
   we = -1.0_dp
   we(1, 1, 1) = 1.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            we=we, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 30010) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array WE(n,q,q)
   n = 5
   m = 2
   q = 3
   np = 2
   deallocate (beta, y, x, we)
   allocate (beta(np), y(n, q), x(n, m), we(n, q, q))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp
   we = -1.0_dp
   we(1, 1, 1) = 1.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            we=we, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 30010) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array WD(1,1,m)
   n = 5
   m = 2
   q = 3
   np = 2
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m), wd(1, 1, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp
   wd = -1.0_dp
   wd(1, 1, 1) = 1.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            wd=wd, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 30001) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array WD(1,m,m)
   n = 5
   m = 2
   q = 3
   np = 2
   deallocate (beta, y, x, wd)
   allocate (beta(np), y(n, q), x(n, m), wd(1, m, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp
   wd = -1.0_dp
   wd(1, 1, 1) = 1.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            wd=wd, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 30001) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array WD(n,1,m)
   n = 5
   m = 2
   q = 3
   np = 2
   deallocate (beta, y, x, wd)
   allocate (beta(np), y(n, q), x(n, m), wd(n, 1, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp
   wd = -1.0_dp
   wd(1, 1, 1) = 1.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            wd=wd, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 30001) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! Inconsistent values of optional array WD(n,m,m)
   n = 5
   m = 2
   q = 3
   np = 2
   deallocate (beta, y, x, wd)
   allocate (beta(np), y(n, q), x(n, m), wd(n, m, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp
   wd = -1.0_dp
   wd(1, 1, 1) = 1.0_dp

   call odr(model, n, m, q, np, beta, y, x, &
            wd=wd, &
            iprint=1, info=info, lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 30001) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! ERROR IN JOB SPECIFICATION WITH RWORK AND IWORK
   n = 1
   m = 1
   q = 1
   np = 1
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, iprint=1, info=info, job=10000, &
            lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 70110) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! ERROR IN JOB SPECIFICATION WITH DELTA
   n = 1
   m = 1
   q = 1
   np = 1
   deallocate (beta, y, x)
   allocate (beta(np), y(n, q), x(n, m))
   y = 0.0_dp
   x = 0.0_dp
   beta = 0.0_dp

   call odr(model, n, m, q, np, beta, y, x, iprint=1, info=info, job=1000, &
            lunrpt=lunrpt, lunerr=lunerr)

   if (info /= 71000) passed = .false.

   write (lunrpt, *) "INFO = ", info

   ! BOUNDS TOO SMALL FOR DERIVATIVE CHECKER WHEN DERIVATIVES DON'T AGREE.
   n = 4
   m = 1
   q = 1
   np = 2
   deallocate (beta, y, x, lower, upper)
   allocate (beta(np), y(n, q), x(n, m), lower(np), upper(np))
   beta = [-200.0_dp, -5.0_dp]
   upper = [-200.0_dp, 0.0_dp]
   lower = [-200.000029802322_dp, -5.0_dp]
   y(:, 1) = [2.718281828459045_dp, 7.389056098930650_dp, &
              148.4131591025766_dp, 403.4287934927353_dp]
   x(:, 1) = [1.0_dp, 2.0_dp, 5.0_dp, 6.0_dp]

   call odr(model, n, m, q, np, beta, y, x, iprint=1, info=info, job=0020, &
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
!       Q  = 1
!       NP = 1
!       DEALLOCATE(BETA,Y,X)
!       ALLOCATE(BETA(NP),Y(N,Q),X(N,M),STAT=STAT)
!       IF (STAT/=0) THEN
!       WRITE(0,*)
!       &       "SYSTEM ERROR: COULD NOT ALLOCATE MEMORY, TESTER ",
!       &       "FAILED TO RUN."
!       STOP
!       END IF
!       Y(:,:) = 0.0_dp
!       X(:,:) = 0.0_dp
!       BETA(:) = 0.0_dp
!
!       CALL ODR(FCN,N,M,Q,NP,BETA,Y,X,IPRINT=1,INFO=INFO,
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
