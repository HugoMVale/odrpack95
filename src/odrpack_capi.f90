module odrpack_capi
   !! C-bindings for 'odrpack'.

   use, intrinsic :: iso_c_binding, only: c_bool, c_char, c_double, c_f_pointer, c_int, c_ptr
   use odrpack, only: odr, workspace_dimensions
   use odrpack_core, only: diwinf, dwinf
   implicit none

   abstract interface
      subroutine fcn_tc(n, m, np, nq, ldn, ldm, ldnp, beta, xplusd, ifixb, ifixx, ldifx, &
                        ideval, f, fjacb, fjacd, istop) bind(C)
      !! User-supplied subroutine for evaluating the model.
         import :: c_int, c_double
         implicit none
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
   end type workidx_t

   type, bind(C) :: iworkidx_t
   !! 0-based indices of the variables stored in the integer work array.
      integer(c_int) :: msgb
      integer(c_int) :: msgd
      integer(c_int) :: ifix2
      integer(c_int) :: istop
      integer(c_int) :: nnzw
      integer(c_int) :: npp
      integer(c_int) :: idf
      integer(c_int) :: job
      integer(c_int) :: iprin
      integer(c_int) :: luner
      integer(c_int) :: lunrp
      integer(c_int) :: nrow
      integer(c_int) :: ntol
      integer(c_int) :: neta
      integer(c_int) :: maxit
      integer(c_int) :: niter
      integer(c_int) :: nfev
      integer(c_int) :: njev
      integer(c_int) :: int2
      integer(c_int) :: irank
      integer(c_int) :: ldtt
      integer(c_int) :: bound
      integer(c_int) :: liwkmn
   end type iworkidx_t

contains

   subroutine odr_short_c( &
      fcn, &
      n, m, np, nq, &
      beta, y, x, &
      delta, &
      lower, upper, &
      job) bind(C)
   !! "Short-call" wrapper for the `odr` routine including very few optional arguments.
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
      real(c_double), intent(inout), optional :: delta(n, m)
         !! Error in the `x` data. Initial guess on input and estimated value on output.
      real(c_double), intent(in), optional :: lower(np)
         !! Lower bound on `beta`.
      real(c_double), intent(in), optional :: upper(np)
         !! Upper bound on `beta`.
      integer(c_int), intent(in), optional :: job
         !! Variable controlling initialization and computational method.

      call odr(fcn, n, m, np, nq, beta, y, x, &
               delta=delta, &
               lower=lower, upper=upper, &
               job=job)

   end subroutine odr_short_c

   subroutine odr_medium_c( &
      fcn, &
      n, m, np, nq, &
      ldwe, ld2we, &
      ldwd, ld2wd, &
      ldifx, &
      beta, y, x, &
      we, wd, &
      ifixb, ifixx, &
      delta, &
      lower, upper, &
      job, iprint, lunerr, lunrpt, info) bind(C)
   !! "Medium-call" wrapper for the `odr` routine including most commonly used optional arguments.
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
      integer(c_int), intent(in) :: ldwe
         !! Leading dimension of array `we`, `ldwe ∈ {1, n}`.
      integer(c_int), intent(in) :: ld2we
         !! Second dimension of array `we`, `ld2we ∈ {1, nq}`.
      integer(c_int), intent(in) :: ldwd
         !! Leading dimension of array `wd`, `ldwd ∈ {1, n}`.
      integer(c_int), intent(in) :: ld2wd
         !! Second dimension of array `wd`, `ld2wd ∈ {1, m}`.
      integer(c_int), intent(in) :: ldifx
         !! Leading dimension of array `ifixx`, `ldifx ∈ {1, n}`.
      real(c_double), intent(inout) :: beta(np)
         !! Function parameters.
      real(c_double), intent(in) :: y(n, nq)
         !! Dependent variable. Unused when the model is implicit.
      real(c_double), intent(in) :: x(n, m)
         !! Explanatory variable.
      real(c_double), intent(in), optional :: we(ldwe, ld2we, nq)
         !! `epsilon` weights.
      real(c_double), intent(in), optional :: wd(ldwd, ld2wd, m)
         !! `delta` weights.
      integer(c_int), intent(in), optional :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer(c_int), intent(in), optional :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      real(c_double), intent(inout), optional :: delta(n, m)
         !! Error in the `x` data. Initial guess on input and estimated value on output.
      real(c_double), intent(in), optional :: lower(np)
         !! Lower bound on `beta`.
      real(c_double), intent(in), optional :: upper(np)
         !! Upper bound on `beta`.
      integer(c_int), intent(in), optional :: job
         !! Variable controlling initialization and computational method.
      integer(c_int), intent(in), optional :: iprint
         !! Print control variable.
      integer(c_int), intent(in), optional :: lunerr
         !! Logical unit number for error messages. Available options are:
         !!   0 => no output.
         !!   6 => output to standard error.
         !!   other => output to logical unit number `lunerr`.
      integer(c_int), intent(in), optional :: lunrpt
         !! Logical unit number for computation reports. Available options are:
         !!   0 => no output.
         !!   6 => output to standard error.
         !!   other => output to logical unit number `lunrpt`.
      integer(c_int), intent(out), optional :: info
         !! Logical unit number for computation reports.

      call odr(fcn, n, m, np, nq, beta, y, x, &
               we=we, wd=wd, &
               ifixb=ifixb, ifixx=ifixx, &
               delta=delta, &
               lower=lower, upper=upper, &
               job=job, &
               iprint=iprint, lunerr=lunerr, lunrpt=lunrpt, &
               info=info)

   end subroutine odr_medium_c

   subroutine odr_long_c( &
      fcn, &
      n, m, np, nq, &
      ldwe, ld2we, &
      ldwd, ld2wd, &
      ldifx, &
      ldstpd, ldscld, &
      lwork, liwork, &
      beta, y, x, &
      we, wd, &
      ifixb, ifixx, &
      stpb, stpd, &
      sclb, scld, &
      delta, &
      lower, upper, &
      work, iwork, &
      job, ndigit, taufac, &
      sstol, partol, maxit, &
      iprint, lunerr, lunrpt, &
      info) bind(C)
   !! "Long-call" wrapper for the `odr` routine including all optional arguments.
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
      integer(c_int), intent(in) :: ldwe
         !! Leading dimension of array `we`, `ldwe ∈ {1, n}`.
      integer(c_int), intent(in) :: ld2we
         !! Second dimension of array `we`, `ld2we ∈ {1, nq}`.
      integer(c_int), intent(in) :: ldwd
         !! Leading dimension of array `wd`, `ldwd ∈ {1, n}`.
      integer(c_int), intent(in) :: ld2wd
         !! Second dimension of array `wd`, `ld2wd ∈ {1, m}`.
      integer(c_int), intent(in) :: ldifx
         !! Leading dimension of array `ifixx`, `ldifx ∈ {1, n}`.
      integer(c_int), intent(in) :: ldstpd
         !! Leading dimension of array `stpd`, `ldstpd ∈ {1, n}`.
      integer(c_int), intent(in) :: ldscld
         !! Leading dimension of array `scld`, `ldscld ∈ {1, n}`.
      integer(c_int), intent(in) :: lwork
         !! Length of array `work`.
      integer(c_int), intent(in) :: liwork
         !! Length of array `iwork`.
      real(c_double), intent(inout) :: beta(np)
         !! Function parameters.
      real(c_double), intent(in) :: y(n, nq)
         !! Dependent variable. Unused when the model is implicit.
      real(c_double), intent(in) :: x(n, m)
         !! Explanatory variable.
      real(c_double), intent(in), optional :: we(ldwe, ld2we, nq)
         !! `epsilon` weights.
      real(c_double), intent(in), optional :: wd(ldwd, ld2wd, m)
         !! `delta` weights.
      integer(c_int), intent(in), optional :: ifixb(np)
         !! Values designating whether the elements of `beta` are fixed at their input values or not.
      integer(c_int), intent(in), optional :: ifixx(ldifx, m)
         !! Values designating whether the elements of `x` are fixed at their input values or not.
      real(c_double), intent(in), optional :: stpb(np)
         !! Relative step for computing finite difference derivatives with respect to `beta`.
      real(c_double), intent(in), optional :: stpd(ldstpd, m)
         !! Relative step for computing finite difference derivatives with respect to `delta`.
      real(c_double), intent(in), optional :: sclb(np)
         !! Scaling values for `beta`.
      real(c_double), intent(in), optional :: scld(ldscld, m)
         !! Scaling values for `delta`.
      real(c_double), intent(inout), optional :: delta(n, m)
         !! Error in the `x` data. `Shape: (n, m)`. Initial guess on input and estimated value
         !! on output.
      real(c_double), intent(in), optional :: lower(np)
         !! Lower bound on `beta`.
      real(c_double), intent(in), optional :: upper(np)
         !! Upper bound on `beta`.
      real(c_double), intent(inout), optional :: work(lwork)
         !! Real work space.
      integer(c_int), intent(inout), optional :: iwork(liwork)
         !! Integer work space.
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
         !! Logical unit number for error messages. Available options are:
         !!   0 => no output.
         !!   6 => output to standard error.
         !!   other => output to logical unit number `lunerr`.
      integer(c_int), intent(in), optional :: lunrpt
         !! Logical unit number for computation reports. Available options are:
         !!   0 => no output.
         !!   6 => output to standard error.
         !!   other => output to logical unit number `lunrpt`.
      integer(c_int), intent(out), optional :: info
         !! Variable designating why the computations were stopped.

      call odr(fcn, n, m, np, nq, beta, y, x, &
               we=we, wd=wd, &
               ifixb=ifixb, ifixx=ifixx, &
               delta=delta, &
               lower=lower, upper=upper, &
               job=job, &
               ndigit=ndigit, taufac=taufac, sstol=sstol, partol=partol, maxit=maxit, &
               iprint=iprint, lunerr=lunerr, lunrpt=lunrpt, &
               stpb=stpb, stpd=stpd, &
               sclb=sclb, scld=scld, &
               info=info, &
               work=work, iwork=iwork)

   end subroutine odr_long_c

   subroutine open_file(filename_cptr, lun, ierr) bind(C)
   !! Open a new file associated with a specified logical unit number.
      type(c_ptr), intent(in), value :: filename_cptr
         !! C-string containing the file name.
      integer(c_int), intent(inout) :: lun
         !! Logical unit number. If `lun > 0`, the user-supplied logical unit number is used.
         !! Otherwise, a new logical unit number is assigned.
      integer(c_int), intent(out) :: ierr
         !! Error code returned by `iostat`.

      character(kind=c_char), pointer :: filename_fptr(:)
      character(len=:), allocatable :: filename
      character(len=256) :: errmsg
      integer :: length, i

      length = strlen(filename_cptr)
      allocate (filename_fptr(length))
      call c_f_pointer(cptr=filename_cptr, fptr=filename_fptr, shape=[length])

      allocate (character(len=length) :: filename)
      do i = 1, length
         filename(i:i) = filename_fptr(i)
      end do

      if (lun > 0) then
         open (file=filename, unit=lun, status='replace', iostat=ierr, iomsg=errmsg)
      else
         open (file=filename, newunit=lun, status='replace', iostat=ierr, iomsg=errmsg)
      end if

      if (ierr /= 0) then
         print *, "I/O error: ", trim(errmsg)
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

   pure subroutine diwinf_c(m, np, nq, iworkidx) bind(C)
   !! Get storage locations within integer work space.
      integer(c_int), intent(in) :: m
         !! The number of columns of data in the independent variable.
      integer(c_int), intent(in) :: np
         !! The number of function parameters.
      integer(c_int), intent(in) :: nq
         !! The number of responses per observation.
      type(iworkidx_t), intent(out) :: iworkidx
         !! 0-based indexes of integer work array.

      integer :: msgbi, msgdi, ifix2i, istopi, nnzwi, nppi, idfi, jobi, iprini, luneri, &
                 lunrpi, nrowi, ntoli, netai, maxiti, niteri, nfevi, njevi, int2i, iranki, &
                 ldtti, boundi, liwkmn

      call diwinf(m, np, nq, &
                  msgbi, msgdi, ifix2i, istopi, &
                  nnzwi, nppi, idfi, &
                  jobi, iprini, luneri, lunrpi, &
                  nrowi, ntoli, netai, &
                  maxiti, niteri, nfevi, njevi, int2i, iranki, ldtti, &
                  boundi, &
                  liwkmn)

      iworkidx%msgb = msgbi - 1
      iworkidx%msgd = msgdi - 1
      iworkidx%ifix2 = ifix2i - 1
      iworkidx%istop = istopi - 1
      iworkidx%nnzw = nnzwi - 1
      iworkidx%npp = nppi - 1
      iworkidx%idf = idfi - 1
      iworkidx%job = jobi - 1
      iworkidx%iprin = iprini - 1
      iworkidx%luner = luneri - 1
      iworkidx%lunrp = lunrpi - 1
      iworkidx%nrow = nrowi - 1
      iworkidx%ntol = ntoli - 1
      iworkidx%neta = netai - 1
      iworkidx%maxit = maxiti - 1
      iworkidx%niter = niteri - 1
      iworkidx%nfev = nfevi - 1
      iworkidx%njev = njevi - 1
      iworkidx%int2 = int2i - 1
      iworkidx%irank = iranki - 1
      iworkidx%ldtt = ldtti - 1
      iworkidx%bound = boundi - 1
      iworkidx%liwkmn = liwkmn

   end subroutine diwinf_c

   pure subroutine dwinf_c(n, m, np, nq, ldwe, ld2we, isodr, workidx) bind(C)
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

      workidx%delta = deltai - 1
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

   pure subroutine workspace_dimensions_c(n, m, np, nq, isodr, lwork, liwork) bind(C)
   !! Calculate the dimensions of the workspace arrays.
      integer(c_int), intent(in) :: n
         !! Number of observations.
      integer(c_int), intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer(c_int), intent(in) :: np
         !! Number of function parameters.
      integer(c_int), intent(in) :: nq
         !! Number of responses per observation.
      logical(c_bool), intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`)
         !! or by OLS (`isodr = .false.`).
      integer(c_int), intent(out) :: lwork
         !! Length of real `work` array.
      integer(c_int), intent(out) :: liwork
         !! Length of integer `iwork` array.

      call workspace_dimensions(n, m, np, nq, logical(isodr, kind=kind(.true.)), &
                                lwork, liwork)

   end subroutine workspace_dimensions_c

end module odrpack_capi
