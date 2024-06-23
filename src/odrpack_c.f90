module odrpack_c
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
   use odrpack, only: odr
   use odrpack_kinds, only: wp
   implicit none
   private

   public :: odr_short_c, open_file, close_file
   
   abstract interface
      subroutine fcn_tc(n, m, np, nq, ldn, ldm, ldnp, beta, xplusd, ifixb, ifixx, ldifx, &
                        ideval, f, fjacb, fjacd, istop) bind(C)
         !! User-supplied subroutine for evaluating the model.
         import :: c_int, c_double
         integer(c_int), intent(in) :: n
         integer(c_int), intent(in) :: m
         integer(c_int), intent(in) :: np
         integer(c_int), intent(in) :: nq
         integer(c_int), intent(in) :: ldn
         integer(c_int), intent(in) :: ldm
         integer(c_int), intent(in) :: ldnp
         real(c_double), intent(in) :: beta(np)
         real(c_double), intent(in) :: xplusd(ldn, m)
         integer(c_int), intent(in) :: ifixb(np)
         integer(c_int), intent(in) :: ifixx(ldifx, m)
         integer(c_int), intent(in) :: ldifx
         integer(c_int), intent(in) :: ideval
         real(c_double), intent(out) :: f(ldn, nq)
         real(c_double), intent(out) :: fjacb(ldn, ldnp, nq)
         real(c_double), intent(out) :: fjacd(ldn, ldm, nq)
         integer(c_int), intent(out) :: istop
      end subroutine
   end interface

contains

   subroutine odr_short_c(fcn, n, m, np, nq, beta, y, x, job, lower, upper) bind(C)
      !! "Short" wrapper for the ODR routine including only principal arguments.
      procedure(fcn_tc) :: fcn
         !! User-supplied subroutine for evaluating the model.
      integer(c_int), intent(in) :: n
         !! Number of observations.
      integer(c_int), intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer(c_int), intent(in) :: np
         !! Number of function parameters.
      integer(c_int), intent(in) :: nq
         !! Number of responses per observation.
      real(c_double), intent(inout) :: beta(np)
         !! Function parameters.
      real(c_double), intent(in) :: y(n, nq)
         !! Dependent variable. Unused when the model is implicit.
      real(c_double), intent(in) :: x(n, m)
         !! Explanatory variable.
      integer(c_int), intent(in), optional :: job
         !! Variable controlling problem initialization and computational method.
      real(c_double), intent(in), optional :: lower(np)
         !! Lower bound on `beta`.
      real(kind=wp), intent(in), optional :: upper(np)
         !! Upper bound on `beta`.

      call odr(fcn, n, m, np, nq, beta, y, x, job=job, lower=lower, upper=upper)

   end subroutine odr_short_c

   ! subroutine odr_c( &
   !    fcn, &
   !    n, m, np, nq, &
   !    beta, &
   !    y, x, &
   !    delta, &
   !    we, wd, &
   !    ifixb, ifixx, &
   !    job, ndigit, taufac, &
   !    sstol, partol, maxit, &
   !    iprint, lunerr, lunrpt, &
   !    stpb, stpd, &
   !    sclb, scld, &
   !    info, &
   !    lower, upper) bind(C)
   !    !! Driver routine for finding the weighted explicit or implicit orthogonal distance
   !    !! regression (ODR) or ordinary linear or nonlinear least squares (OLS) solution (long call
   !    !! statement).
   !    procedure(fcn_tc) :: fcn
   !       !! User-supplied subroutine for evaluating the model.
   !    integer(c_int), intent(in), value :: n
   !       !! Number of observations.
   !    integer(c_int), intent(in), value :: m
   !       !! Number of columns of data in the independent variable.
   !    integer(c_int), intent(in), value :: np
   !       !! Number of function parameters.
   !    integer(c_int), intent(in), value :: nq
   !       !! Number of responses per observation.
   !    real(c_double), intent(inout) :: beta(np)
   !       !! Function parameters.
   !    real(c_double), intent(in) :: y(n, nq)
   !       !! Dependent variable. Unused when the model is implicit.
   !    real(c_double), intent(in) :: x(n, m)
   !       !! Explanatory variable.
   !    real(c_double), intent(inout), optional :: delta(n, m)
   !       !! Initial error in the `x` data.
   !    real(c_double), intent(in), optional :: we(ldwe, ld2we, nq)
   !       !! `epsilon` weights.
   !    real(c_double), intent(in), optional :: wd(ldwd, ld2wd, m)
   !       !! `delta` weights.
   !    integer(c_int), intent(in), optional :: ifixb(np)
   !       !! Values designating whether the elements of `beta` are fixed at their input values or not.
   !    integer(c_int), intent(in), optional :: ifixx(ldifx, m)
   !    !! Values designating whether the elements of `x` are fixed at their input values or not.
   !    integer(c_int), intent(in), optional :: job
   !       !! Variable controlling problem initialization and computational method.
   !    integer(c_int), intent(in), optional :: ndigit
   !       !! Number of accurate digits in the function results, as supplied by the user.
   !    real(c_double), intent(in), optional :: taufac
   !       !! Factor used to compute the initial trust region diameter.
   !    real(c_double), intent(in), optional :: sstol
   !       !! Sum-of-squares convergence stopping tolerance.
   !    real(c_double), intent(in), optional :: partol
   !       !! Parameter convergence stopping tolerance.
   !    integer(c_int), intent(in), optional :: maxit
   !       !! Maximum number of iterations allowed.
   !    integer(c_int), intent(in), optional :: iprint
   !       !! Print control variable.
   !    integer(c_int), intent(in), optional :: lunerr
   !       !! Logical unit number for error messages.
   !    integer(c_int), intent(in), optional :: lunrpt
   !       !! Logical unit number for computation reports.
   !    real(c_double), intent(in), optional :: stpb(np)
   !       !! Relative step for computing finite difference derivatives with respect to `beta`.
   !    integer(c_int), intent(in) :: ldstpd
   !    real(c_double), intent(in), optional :: stpd(ldstpd, m)
   !       !! Relative step for computing finite difference derivatives with respect to `delta`.
   !    real(c_double), intent(in), optional :: sclb(np)
   !       !! Scaling values for `beta`.
   !    real(c_double), intent(in), optional :: scld(ldscld, m)
   !       !! Scaling values for `delta`.
   !    integer(c_int), intent(out), optional :: info
   !       !! Variable designating why the computations were stopped.
   !    real(c_double), intent(in), optional :: lower(np)
   !       !! Lower bound on `beta`.
   !    real(c_double), intent(in), optional :: upper(np)
   !       !! Upper bound on `beta`.

   !    call odr(fcn, n, m, np, nq, beta, y, x, &
   !             delta=delta, we=we, wd=wd, ifixb=ifixb, ifixx=ifixx, job=job, ndigit=ndigit, &
   !             taufac=taufac, sstol=sstol, partol=partol, maxit=maxit, iprint=iprint, &
   !             lunerr=lunerr, lunrpt=lunrpt, stpb=stpb, stpd=stpd, sclb=sclb, scld=scld, &
   !             info=info, lower=lower, upper=upper)

   ! end subroutine

   subroutine open_file(lun, fn, len, ierr) bind(C)
   !! Open a new file associated with a specified logical unit number.
      integer(c_int), intent(inout) :: lun
         !! Logical unit number. If `lun>0`, the user-supplied logical unit number is used.
         !! Otherwise, a new logical unit number is assigned.
      integer(c_int), intent(in) :: len
         !! Length of the string containing the file name, i.e. `strlen(fn)`.
      character(kind=c_char), intent(in) :: fn(len)
         !! String containing the file name.
      integer(c_int), intent(out) :: ierr
         !! Error code returned by `iostat`.

      character(len=len) :: fn_
      integer :: i

      do i = 1, len
         fn_(i:i) = fn(i)
      end do 

      if (lun > 0) then
         open (file=fn_, unit=lun, status='new', iostat=ierr)
      else
         open (file=fn_, newunit=lun, status='new', iostat=ierr)
      end if

   end subroutine

   subroutine close_file(lun) bind(C)
   !! Close a file associated with a specified logical unit number.
      integer(c_int), intent(in) :: lun
         !! Logical unit number.
      close (unit=lun)
   end subroutine

end module odrpack_c
