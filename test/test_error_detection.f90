program tester
!***BEGIN PROLOGUE  TESTER
!***REFER TO ODR
!***ROUTINES CALLED  ODR
!***DATE WRITTEN   20040322   (YYYYMMDD)
!***REVISION DATE  20040322   (YYYYMMDD)
!***PURPOSE  EXCERCISE ERROR REPORTING OF THE F90 VERSION OF ODRPACK95
!***END PROLOGUE  TESTER

!...USED MODULES
   use odrpack_kinds, only: wp
   use odrpack, only: odr
   implicit none

!...LOCAL SCALARS
   integer :: n, m, nq, np, info, lun
   logical :: passed
!  STAT

!...LOCAL ARRAYS
   real(wp) :: beta(:), y(:, :), x(:, :), upper(2), lower(2)

!...ALLOCATABLE ARRAYS
   allocatable :: beta, y, x

!...EXTERNAL SUBPROGRAMS
   external fcn

   common/bounds/upper, lower

!***FIRST EXECUTABLE STATEMENT  TESTER

   passed = .true.

   open (unit=8, file="./test/test_error_summary.txt")
   write (8, *) "NO SUMMARY AVAILABLE"
   close (8)

   lun = 9
   open (unit=lun, file="./test/test_error_report.txt")

!  ERROR IN PROBLEM SIZE

   n = 0
   m = 0
   nq = 0
   np = 0
   allocate (beta(np), y(n, nq), x(n, m))
   y(:, :) = 0.0_wp
   x(:, :) = 0.0_wp
   beta(:) = 0.0_wp

   call odr(fcn, n, m, np, nq, beta, y, x, iprint=1, info=info, &
            lunrpt=lun, lunerr=lun)

   if (info .ne. 11111) passed = .false.

   write (lun, *) "INFO = ", info

!  ERROR IN JOB SPECIFICATION WITH WORK AND IWORK

   n = 1
   m = 1
   nq = 1
   np = 1
   deallocate (beta, y, x)
   allocate (beta(np), y(n, nq), x(n, m))
   y(:, :) = 0.0_wp
   x(:, :) = 0.0_wp
   beta(:) = 0.0_wp

   call odr(fcn, n, m, np, nq, beta, y, x, iprint=1, info=info, job=10000, &
            lunrpt=lun, lunerr=lun)

   if (info .ne. 70110) passed = .false.

   write (lun, *) "INFO = ", info

!  ERROR IN JOB SPECIFICATION WITH DELTA

   n = 1
   m = 1
   nq = 1
   np = 1
   deallocate (beta, y, x)
   allocate (beta(np), y(n, nq), x(n, m))
   y(:, :) = 0.0_wp
   x(:, :) = 0.0_wp
   beta(:) = 0.0_wp

   call odr(fcn, n, m, np, nq, beta, y, x, iprint=1, info=info, job=1000, &
            lunrpt=lun, lunerr=lun)

   if (info .ne. 71000) passed = .false.

   write (lun, *) "INFO = ", info

!  BOUNDS TOO SMALL FOR DERIVATIVE CHECKER WHEN DERIVATIVES DON'T AGREE.

   n = 4
   m = 1
   nq = 1
   np = 2
   deallocate (beta, y, x)
   allocate (beta(np), y(n, nq), x(n, m))
   beta(:) = (/-200.0_wp, -5.0_wp/)
   upper(1:2) = (/-200.0_wp, 0.0_wp/)
   lower(1:2) = (/-200.000029802322_wp, -5.0_wp/)
   y(:, 1) = (/2.718281828459045_wp, 7.389056098930650_wp, &
               148.4131591025766_wp, 403.4287934927353_wp/)
   x(:, 1) = (/1.0_wp, 2.0_wp, 5.0_wp, 6.0_wp/)

   call odr(fcn, n, m, np, nq, beta, y, x, iprint=1, info=info, job=0020, &
            lunrpt=lun, lunerr=lun, lower=lower, upper=upper)

   if (info .ne. 1) passed = .false.

   write (lun, *) "INFO = ", info

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
!       IF (STAT.NE.0) THEN
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
   close (lun)

   if (passed) then
      stop "Error detection tests passed. See report."
   else
      error stop "Error detection tests failed. See report."
   end if

end program tester

subroutine fcn &
   (n, m, np, nq, &
    ldn, ldm, ldnp, &
    beta, xplusd, &
    ifixb, ifixx, ldifx, &
    ideval, f, fjacb, fjacd, &
    istop)
!***BEGIN PROLOGUE  FCN
!***REFER TO  ODR
!***ROUTINES CALLED  (NONE)
!***DATE WRITTEN   20040322 (YYYYMMDD)
!***REVISION DATE  20040322 (YYYYMMDD)
!***PURPOSE  DUMMY ROUTINE FOR ODRPACK95 ERROR EXERCISER
!***END PROLOGUE  FCN

!...USED MODULES
   use odrpack_kinds, only: wp, zero
   implicit none

!...SCALAR ARGUMENTS
   integer :: ideval, istop, ldifx, ldm, ldn, ldnp, m, n, np, nq

!...ARRAY ARGUMENTS
   real(wp) :: beta(np), f(ldn, nq), fjacb(ldn, ldnp, nq), fjacd(ldn, ldm, nq), &
                    xplusd(ldn, m)
   integer :: ifixb(np), ifixx(ldifx, m)

!...ARRAYS IN COMMON
   real(wp) :: lower(2), upper(2)

!...LOCAL SCALARS
   integer :: i

   common /bounds/ upper, lower

!***FIRST EXECUTABLE STATEMENT
!
!  Do something with FJACD, FJACB, IFIXB and IFIXX to avoid warnings that they
!  are not being used.  This is simply not to worry users that the example
!  program is failing.
   if (ifixb(1) .gt. 0 .and. ifixx(1, 1) .gt. 0 .and. fjacb(1, 1, 1) .gt. zero &
       .and. fjacd(1, 1, 1) .gt. zero) then
      ! Do nothing.
   end if

   if (any(lower(1:np) .gt. beta(1:np))) then
      write (0, *) "LOWER BOUNDS VIOLATED"
      do i = 1, np
         if (lower(i) .gt. beta(i)) then
            write (0, *) "   IN THE ", i, " POSITION WITH ", beta(i), &
               "<", lower(i)
         end if
      end do
   end if

   if (any(upper(1:np) .lt. beta(1:np))) then
      write (0, *) "UPPER BOUNDS VIOLATED"
      do i = 1, np
         if (upper(i) .lt. beta(i)) then
            write (0, *) "   IN THE ", i, " POSITION WITH ", beta(i), &
               ">", upper(i)
         end if
      end do
   end if

   istop = 0

   if (mod(ideval, 10) .ne. 0) then
      do i = 1, n
         f(i, 1) = beta(1)*exp(beta(2)*xplusd(i, 1))
      end do
   end if

   if (mod(ideval/10, 10) .ne. 0) then
      do i = 1, n
         fjacb(i, 1, 1) = exp(beta(2)*xplusd(i, 1))
         fjacb(i, 2, 1) = beta(1)*xplusd(i, 1)*exp(beta(2)*xplusd(i, 1))
      end do
   end if

   if (mod(ideval/100, 10) .ne. 0) then
      do i = 1, n
         fjacd(i, 1, 1) = beta(1)*beta(2)*exp(beta(2)*xplusd(i, 1))
      end do
   end if

end subroutine fcn
