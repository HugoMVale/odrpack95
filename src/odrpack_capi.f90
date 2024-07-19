module odrpack_capi
   !! C bindings for 'odrpack'.

   use, intrinsic :: iso_c_binding, only: c_bool, c_char, c_double, c_f_pointer, c_int, c_ptr
   use odrpack_kinds, only: wp
   use odrpack, only: odr
   use odrpack_core, only: dwinf
   implicit none
   private

   public :: odr_basic_c
   public :: odr_short_c
   public :: odr_long_c
   public :: fcn_tc
   public :: dwinf_c
   public :: workidx_t
   public :: open_file
   public :: close_file

   abstract interface
      subroutine fcn_tc(n, m, np, nq, ldn, ldm, ldnp, beta, xplusd, ifixb, ifixx, ldifx, &
                        ideval, f, fjacb, fjacd, istop) bind(C)
      !! User-supplied subroutine for evaluating the model.
         import :: c_int, c_double
         integer(c_int), intent(in) :: n
            !! Number of observations.
         integer(c_int), intent(in) :: m
            !! Number of columns of data in the independent variable.
         integer(c_int), intent(in) :: np
            !! Number of function parameters.
         integer(c_int), intent(in) :: nq
            !! Number of responses per observation.
         integer(c_int), intent(in) :: ldn
            !! Leading dimension declarator equal or exceeding `n`.
         integer(c_int), intent(in) :: ldm
            !! Leading dimension declarator equal or exceeding `m`.
         integer(c_int), intent(in) :: ldnp
            !! Leading dimension declarator equal or exceeding `np`.
         real(c_double), intent(in) :: beta(np)
            !! Current values of parameters.
         real(c_double), intent(in) :: xplusd(ldn, m)
            !! Current value of explanatory variable, i.e., `x + delta`.
         integer(c_int), intent(in) :: ifixb(np)
            !! Indicators for "fixing" parameters (`beta`).
         integer(c_int), intent(in) :: ifixx(ldifx, m)
            !! Indicators for "fixing" explanatory variable (`x`).
         integer(c_int), intent(in) :: ldifx
            !! Leading dimension of array `ifixx`.
         integer(c_int), intent(in) :: ideval
            !! Indicator for selecting computation to be performed.
         real(c_double), intent(out) :: f(ldn, nq)
            !! Predicted function values.
         real(c_double), intent(out) :: fjacb(ldn, ldnp, nq)
            !! Jacobian with respect to `beta`.
         real(c_double), intent(out) :: fjacd(ldn, ldm, nq)
            !! Jacobian with respect to errors `delta`.
         integer(c_int), intent(out) :: istop
            !! Stopping condition.
      end subroutine
   end interface

   interface
      integer(c_int) function strlen(string) bind(C)
      !! Length of a C-string.
         import :: c_int, c_ptr
         type(c_ptr), intent(in), value :: string
      end function
   end interface

   type, bind(C) :: workidx_t
   !! 0-based indices of the variables stored in the real work array.
      integer(c_int) :: delta
      integer(c_int) :: eps
      integer(c_int) :: xplus
      integer(c_int) :: fn
      integer(c_int) :: sd
      integer(c_int) :: vcv
      integer(c_int) :: rvar
      integer(c_int) :: wss
      integer(c_int) :: wssde
      integer(c_int) :: wssep
      integer(c_int) :: rcond
      integer(c_int) :: eta
      integer(c_int) :: olmav
      integer(c_int) :: tau
      integer(c_int) :: alpha
      integer(c_int) :: actrs
      integer(c_int) :: pnorm
      integer(c_int) :: rnors
      integer(c_int) :: prers
      integer(c_int) :: partl
      integer(c_int) :: sstol
      integer(c_int) :: taufc
      integer(c_int) :: epsma
      integer(c_int) :: beta0
      integer(c_int) :: betac
      integer(c_int) :: betas
      integer(c_int) :: betan
      integer(c_int) :: s
      integer(c_int) :: ss
      integer(c_int) :: ssf
      integer(c_int) :: qraux
      integer(c_int) :: u
      integer(c_int) :: fs
      integer(c_int) :: fjacb
      integer(c_int) :: we1
      integer(c_int) :: diff
      integer(c_int) :: delts
      integer(c_int) :: deltn
      integer(c_int) :: t
      integer(c_int) :: tt
      integer(c_int) :: omega
      integer(c_int) :: fjacd
      integer(c_int) :: wrk1
      integer(c_int) :: wrk2
      integer(c_int) :: wrk3
      integer(c_int) :: wrk4
      integer(c_int) :: wrk5
      integer(c_int) :: wrk6
      integer(c_int) :: wrk7
      integer(c_int) :: lower
      integer(c_int) :: upper
      integer(c_int) :: lwkmn
   end type

contains

   subroutine odr_basic_c(fcn, n, m, np, nq, beta, y, x, lower, upper, job) bind(C)
   !! "Basic" wrapper for the `odr` routine including mandatory arguments and very few optional
   !! arguments.
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
      real(c_double), intent(in), optional :: lower(np)
         !! Lower bound on `beta`.
      real(wp), intent(in), optional :: upper(np)
         !! Upper bound on `beta`.
      integer(c_int), intent(in), optional :: job
         !! Variable controlling initialization and computational method.

      call odr(fcn, n, m, np, nq, beta, y, x, lower=lower, upper=upper, job=job)

   end subroutine odr_basic_c

   subroutine odr_short_c( &
      fcn, &
      n, m, np, nq, &
      beta, y, x, &
      we, ldwe, ld2we, &
      wd, ldwd, ld2wd, &
      lower, upper, &
      job, &
      iprint, lunerr, lunrpt, &
      info) bind(C)
   !! "Short" wrapper for the `odr` routine including mandatory arguments and most commonly used
   !! optional arguments. Similar to the short-call statement of the original ODRPACK `DODR`.
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
      real(c_double), intent(in) :: we(ldwe, ld2we, nq)
         !! `epsilon` weights.
      integer(c_int), intent(in) :: ldwe
         !! Leading dimension of array `we`, `ldwe ∈ {1, n}`.
      integer(c_int), intent(in) :: ld2we
         !! Second dimension of array `we`, `ld2we ∈ {1, n}`.
      real(c_double), intent(in) :: wd(ldwd, ld2wd, m)
         !! `delta` weights.
      integer(c_int), intent(in) :: ldwd
         !! Leading dimension of array `wd`, `ldwd ∈ {1, n}`.
      integer(c_int), intent(in) :: ld2wd
         !! Second dimension of array `wd`, `ld2wd ∈ {1, m}`.
      real(c_double), intent(in), optional :: lower(np)
         !! Lower bound on `beta`.
      real(wp), intent(in), optional :: upper(np)
         !! Upper bound on `beta`.
      integer(c_int), intent(in), optional :: job
         !! Variable controlling initialization and computational method.
      integer(c_int), intent(in), optional :: iprint
         !! Print control variable.
      integer(c_int), intent(in), optional :: lunerr
         !! Logical unit number for error messages.
      integer(c_int), intent(in), optional :: lunrpt
         !! Logical unit number for computation reports.
      integer(c_int), intent(out), optional :: info
         !! Logical unit number for computation reports.

      call odr(fcn, n, m, np, nq, beta, y, x, &
               we=we, wd=wd, &
               lower=lower, upper=upper, &
               job=job, &
               iprint=iprint, lunerr=lunerr, lunrpt=lunrpt, &
               info=info)

   end subroutine odr_short_c

   subroutine odr_long_c( &
      fcn, &
      n, m, np, nq, &
      beta, y, x, &
      we, ldwe, ld2we, &
      wd, ldwd, ld2wd, &
      ifixb, ifixx, ldifx, &
      stpb, stpd, ldstpd, &
      sclb, scld, ldscld, &
      lower, upper, &
      delta, &
      job, ndigit, taufac, &
      sstol, partol, maxit, &
      iprint, lunerr, lunrpt, &
      info) bind(C)
   !! "Long" wrapper for the `odr` routine including mandatory arguments and all optional
   !! arguments.
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
      real(c_double), intent(in) :: we(ldwe, ld2we, nq)
         !! `epsilon` weights.
      integer(c_int), intent(in) :: ldwe
         !! Leading dimension of array `we`, `ldwe ∈ {1, n}`.
      integer(c_int), intent(in) :: ld2we
         !! Second dimension of array `we`, `ld2we ∈ {1, n}`.
      real(c_double), intent(in) :: wd(ldwd, ld2wd, m)
         !! `delta` weights.
      integer(c_int), intent(in) :: ldwd
         !! Leading dimension of array `wd`, `ldwd ∈ {1, n}`.
      integer(c_int), intent(in) :: ld2wd
         !! Second dimension of array `wd`, `ld2wd ∈ {1, m}`.
      integer(c_int), intent(in) :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer(c_int), intent(in) :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      integer(c_int), intent(in) :: ldifx
         !! Leading dimension of array `ifixx`, `ldifx ∈ {1, n}`.
      real(c_double), intent(in) :: stpb(np)
         !! Relative step for computing finite difference derivatives with respect to `beta`.
      real(c_double), intent(in) :: stpd(ldstpd, m)
         !! Relative step for computing finite difference derivatives with respect to `delta`.
      integer(c_int), intent(in) :: ldstpd
         !! Leading dimension of array `stpd`, `ldstpd ∈ {1, n}`.
      real(c_double), intent(in) :: sclb(np)
         !! Scaling values for `beta`.
      real(c_double), intent(in) :: scld(ldscld, m)
         !! Scaling values for `delta`.
      integer(c_int), intent(in) :: ldscld
         !! Leading dimension of array `scld`, `ldscld ∈ {1, n}`.
      real(c_double), intent(in), optional :: lower(np)
         !! Lower bound on `beta`.
      real(wp), intent(in), optional :: upper(np)
         !! Upper bound on `beta`.
      real(c_double), intent(inout), optional :: delta(n, m)
         !! Initial error in the `x` data.
      integer(c_int), intent(in), optional :: job
         !! Variable controlling initialization and computational method.
      integer(c_int), intent(in), optional :: ndigit
         !! Number of accurate digits in the function results, as supplied by the user.
      real(c_double), intent(in), optional :: taufac
         !! Factor used to compute the initial trust region diameter.
      real(c_double), intent(in), optional :: sstol
         !! Sum-of-squares convergence stopping tolerance.
      real(c_double), intent(in), optional :: partol
         !! Parameter convergence stopping tolerance.
      integer(c_int), intent(in), optional :: maxit
         !! Maximum number of iterations allowed.
      integer(c_int), intent(in), optional :: iprint
         !! Print control variable.
      integer(c_int), intent(in), optional :: lunerr
         !! Logical unit number for error messages.
      integer(c_int), intent(in), optional :: lunrpt
         !! Logical unit number for computation reports.
      integer(c_int), intent(out), optional :: info
         !! Variable designating why the computations were stopped.

      call odr(fcn, n, m, np, nq, beta, y, x, &
               delta=delta, &
               we=we, wd=wd, &
               lower=lower, upper=upper, &
               ifixb=ifixb, &
               ifixx=ifixx, &
               job=job, ndigit=ndigit, taufac=taufac, &
               sstol=sstol, partol=partol, maxit=maxit, &
               iprint=iprint, lunerr=lunerr, lunrpt=lunrpt, &
               stpb=stpb, &
               stpd=stpd, &
               sclb=sclb, &
               scld=scld, &
               info=info)

   end subroutine odr_long_c

   subroutine open_file(lun, filename_cptr, ierr) bind(C)
   !! Open a new file associated with a specified logical unit number.
      integer(c_int), intent(inout) :: lun
         !! Logical unit number. If `lun > 0`, the user-supplied logical unit number is used.
         !! Otherwise, a new logical unit number is assigned.
      type(c_ptr), intent(in), value :: filename_cptr
         !! C-string containing the file name.
      integer(c_int), intent(out) :: ierr
         !! Error code returned by `iostat`.

      character(kind=c_char), pointer :: filename_fptr(:)
      character(len=:), allocatable :: filename
      integer :: length, i

      length = strlen(filename_cptr)
      allocate (filename_fptr(length))
      call c_f_pointer(cptr=filename_cptr, fptr=filename_fptr, shape=[length])

      allocate (character(len=length) :: filename)
      do i = 1, length
         filename(i:i) = filename_fptr(i)
      end do

      if (lun > 0) then
         open (file=filename, unit=lun, status='replace', iostat=ierr)
      else
         open (file=filename, newunit=lun, status='replace', iostat=ierr)
      end if

   end subroutine open_file

   subroutine close_file(lun, ierr) bind(C)
   !! Close a file associated with a specified logical unit number.
      integer(c_int), intent(in) :: lun
         !! Logical unit number.
      integer(c_int), intent(out) :: ierr
         !! Error code returned by `iostat`.

      close (unit=lun, iostat=ierr)

   end subroutine close_file

   subroutine dwinf_c(n, m, np, nq, ldwe, ld2we, isodr, workidx) bind(C)
   !! Get storage locations within real work space.
      integer(c_int), intent(in) :: n
         !! Number of observations.
      integer(c_int), intent(in) :: m
         !! Number of columns of data in the explanatory variable.
      integer(c_int), intent(in) :: np
         !! Number of function parameters.
      integer(c_int), intent(in) :: nq
         !! Number of responses per observation.
      integer(c_int), intent(in) :: ldwe
         !! Leading dimension of array `we`.
      integer(c_int), intent(in) :: ld2we
         !! Second dimension of array `we`.
      logical(c_bool), intent(in) :: isodr
         !! Variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      type(workidx_t), intent(out) :: workidx
         !! 0-based indexes of real work array.

      integer :: deltai, epsi, xplusi, fni, sdi, vcvi, rvari, wssi, wssdei, wssepi, &
                 rcondi, etai, olmavi, taui, alphai, actrsi, pnormi, rnorsi, prersi, partli, &
                 sstoli, taufci, epsmai, beta0i, betaci, betasi, betani, si, ssi, ssfi, &
                 qrauxi, ui, fsi, fjacbi, we1i, diffi, deltsi, deltni, ti, tti, omegai, &
                 fjacdi, wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, loweri, upperi, &
                 lwkmn

      call dwinf(n, m, np, nq, ldwe, ld2we, logical(isodr, kind=kind(.true.)), &
                 deltai, epsi, xplusi, fni, sdi, vcvi, &
                 rvari, wssi, wssdei, wssepi, rcondi, etai, &
                 olmavi, taui, alphai, actrsi, pnormi, rnorsi, prersi, &
                 partli, sstoli, taufci, epsmai, &
                 beta0i, betaci, betasi, betani, si, ssi, ssfi, qrauxi, ui, &
                 fsi, fjacbi, we1i, diffi, &
                 deltsi, deltni, ti, tti, omegai, fjacdi, &
                 wrk1i, wrk2i, wrk3i, wrk4i, wrk5i, wrk6i, wrk7i, &
                 loweri, upperi, &
                 lwkmn)

      workidx%eps = epsi - 1
      workidx%xplus = xplusi - 1
      workidx%fn = fni - 1
      workidx%sd = sdi - 1
      workidx%vcv = vcvi - 1
      workidx%rvar = rvari - 1
      workidx%wss = wssi - 1
      workidx%wssde = wssdei - 1
      workidx%wssep = wssepi - 1
      workidx%rcond = rcondi - 1
      workidx%eta = etai - 1
      workidx%olmav = olmavi - 1
      workidx%tau = taui - 1
      workidx%alpha = alphai - 1
      workidx%actrs = actrsi - 1
      workidx%pnorm = pnormi - 1
      workidx%rnors = rnorsi - 1
      workidx%prers = prersi - 1
      workidx%partl = partli - 1
      workidx%sstol = sstoli - 1
      workidx%taufc = taufci - 1
      workidx%epsma = epsmai - 1
      workidx%beta0 = beta0i - 1
      workidx%betac = betaci - 1
      workidx%betas = betasi - 1
      workidx%betan = betani - 1
      workidx%s = si - 1
      workidx%ss = ssi - 1
      workidx%ssf = ssfi - 1
      workidx%qraux = qrauxi - 1
      workidx%u = ui - 1
      workidx%fs = fsi - 1
      workidx%fjacb = fjacbi - 1
      workidx%we1 = we1i - 1
      workidx%diff = diffi - 1
      workidx%delts = deltsi - 1
      workidx%deltn = deltni - 1
      workidx%t = ti - 1
      workidx%tt = tti - 1
      workidx%omega = omegai - 1
      workidx%fjacd = fjacdi - 1
      workidx%wrk1 = wrk1i - 1
      workidx%wrk2 = wrk2i - 1
      workidx%wrk3 = wrk3i - 1
      workidx%wrk4 = wrk4i - 1
      workidx%wrk5 = wrk5i - 1
      workidx%wrk6 = wrk6i - 1
      workidx%wrk7 = wrk7i - 1
      workidx%lower = loweri - 1
      workidx%upper = upperi - 1
      workidx%lwkmn = lwkmn

   end subroutine dwinf_c

end module odrpack_capi
