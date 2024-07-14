module odrpack

   use odrpack_kinds, only: wp
   implicit none
   private

   public :: odr, tempret

   ! A temporary work array for holding return values before copying to a lower rank array.
   real(kind=wp), allocatable :: tempret(:, :)

   abstract interface
      subroutine fcn_t(n, m, np, nq, ldn, ldm, ldnp, beta, xplusd, ifixb, ifixx, ldifx, &
                       ideval, f, fjacb, fjacd, istop)
      !! User-supplied subroutine for evaluating the model.
         import :: wp
         integer, intent(in) :: n
            !! Number of observations.
         integer, intent(in) :: m
            !! Number of columns of data in the independent variable.
         integer, intent(in) :: np
            !! Number of function parameters.
         integer, intent(in) :: nq
            !! Number of responses per observation.
         integer, intent(in) :: ldn
            !! Leading dimension declarator equal or exceeding `n`.
         integer, intent(in) :: ldm
            !! Leading dimension declarator equal or exceeding `m`.
         integer, intent(in) :: ldnp
            !! Leading dimension declarator equal or exceeding `np`.
         real(kind=wp), intent(in) :: beta(np)
            !! Current values of parameters.
         real(kind=wp), intent(in) :: xplusd(ldn, m)
            !! Current value of explanatory variable, i.e., `x + delta`.
         integer, intent(in) :: ifixb(np)
            !! Indicators for "fixing" parameters (`beta`).
         integer, intent(in) :: ifixx(ldifx, m)
            !!  Indicators for "fixing" explanatory variable (`x`).	
         integer, intent(in) :: ldifx
            !! Leading dimension of array `ifixx`.
         integer, intent(in) :: ideval
            !! Indicator for selecting computation to be performed.
         real(kind=wp), intent(out) :: f(ldn, nq)
            !! Predicted function values.
         real(kind=wp), intent(out) :: fjacb(ldn, ldnp, nq)
            !! Jacobian with respect to `beta`.
         real(kind=wp), intent(out) :: fjacd(ldn, ldm, nq)
            !! Jacobian with respect to errors `delta`.
         integer, intent(out) :: istop
            !! Stopping condition, with meaning as follows. 0 means current `beta` and
            !! `x+delta` were acceptable and values were computed successfully. 1 means current
            !! `beta` and `x+delta` are not acceptable;  ODRPACK95 should select values closer
            !! to most recently used values if possible. -1 means current `beta` and `x+delta`
            !! are not acceptable; ODRPACK95 should stop.
      end subroutine fcn_t
   end interface

contains

   subroutine odr &
      (fcn, &
       n, m, np, nq, &
       beta, &
       y, x, &
       delta, &
       we, wd, &
       ifixb, ifixx, &
       job, ndigit, taufac, &
       sstol, partol, maxit, &
       iprint, lunerr, lunrpt, &
       stpb, stpd, &
       sclb, scld, &
       work, iwork, &
       info, &
       lower, upper)
   !! Driver routine for finding the weighted explicit or implicit orthogonal distance
   !! regression (ODR) or ordinary linear or nonlinear least squares (OLS) solution (long call
   !! statement).
   ! Date Written:   860529   (YYMMDD)
   ! Revision Date:  20040301 (YYYYMMDD)
   ! Category No.:  G2E,I1B1
   ! Keywords:  Orthogonal distance regression,
   !   Nonlinear least squares,
   !   Measurement error models,
   !   Errors in variables
   ! Authors:
   !   Boggs, Paul T.
   !     Applied and Computational Mathematics Division
   !     National Institute of Standards and Technology
   !     Gaithersburg, MD 20899
   !   Byrd, Richard H.
   !     Department of Computer Science
   !     University of Colorado, Boulder, CO 80309
   !   Rogers, Janet E.
   !     Applied and Computational Mathematics Division
   !     National Institute of Standards and Technology
   !     Boulder, CO 80303-3328
   !   Schnabel, Robert B.
   !     Department of Computer Science
   !     University of Colorado, Boulder, CO 80309
   !     and
   !     Applied and Computational Mathematics Division
   !     National Institute of Standards and Technology
   !     Boulder, CO 80303-3328
   ! Purpose:  REAL (KIND=wp) driver routine for finding the weighted explicit or implicit
   !   orthogonal distance regression (ODR) or ordinary linear or nonlinear least squares (OLS)
   !   solution (long call statement)
   ! Description:
   !   For details, see ODRPACK95 User's Reference Guide.
   ! References:
   !   Boggs, P. T., R. H. Byrd, J. R. Donaldson, and R. B. Schnabel (1989),
   !     "Algorithm 676 --- ODRPACK: Software for Weighted Orthogonal Distance Regression,"
   !     ACM Trans. Math. Software., 15(4):348-364.
   !       Boggs, P. T., R. H. Byrd, J. E. Rogers, and
   !   R. B. Schnabel (1992),
   !     "User's Reference Guide for ODRPACK Version 2.01,
   !     Software for Weighted Orthogonal Distance Regression,"
   !     National Institute of Standards and Technology
   !     Internal Report Number 92-4834.
   !  Boggs, P. T., R. H. Byrd, and R. B. Schnabel (1987),
   !    "A Stable and Efficient Algorithm for Nonlinear Orthogonal Distance Regression,"
   !    SIAM J. Sci. Stat. Comput., 8(6):1052-1078.

      use odrpack_kinds, only: negone, zero

      procedure(fcn_t) :: fcn
         !! User-supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! Number of observations.
      integer, intent(in) :: m
         !! Number of columns of data in the independent variable.
      integer, intent(in) :: np
         !! Number of function parameters.
      integer, intent(in) :: nq
         !! Number of responses per observation.
      real(kind=wp), intent(inout) :: beta(:)
         !! Function parameters. `Shape: (np)`.
      real(kind=wp), intent(in) :: y(:, :)
         !! Dependent variable. `Shape: (n, nq)`. Unused when the model is implicit.
      real(kind=wp), intent(in) :: x(:, :)
         !! Explanatory variable. `Shape: (n, m)`.
      real(kind=wp), intent(inout), optional :: delta(:, :)
         !! Initial error in the `x` data. `Shape: (n, m)`.
      real(kind=wp), intent(in), optional :: we(:, :, :)
         !! `epsilon` weights. `Shape: (1<=ldwe<=n, 1<=ld2we<=nq, nq)`. See p. 25.
      real(kind=wp), intent(in), optional :: wd(:, :, :)
         !! `delta` weights. `Shape: (1<=ldwd<=n, 1<=ld2wd<=m, m)`. See p. 26.
      integer, intent(in), optional :: ifixb(:)
         !! Values designating whether the elements of `beta` are fixed at their input values
         !! or not. `Shape: (np)`.
      integer, intent(in), optional :: ifixx(:, :)
         !! Values designating whether the elements of `x` are fixed at their input values
         !! or not. `Shape: (1<=ldifx<=n, m)`. See p. 27.
      integer, intent(in), optional :: job
         !! Variable controlling problem initialization and computational method.
      integer, intent(in), optional :: ndigit
         !! Number of accurate digits in the function results, as supplied by the user.
      real(kind=wp), intent(in), optional :: taufac
         !! Factor used to compute the initial trust region diameter.
      real(kind=wp), intent(in), optional :: sstol
         !! Sum-of-squares convergence stopping tolerance.
      real(kind=wp), intent(in), optional :: partol
      !! Parameter convergence stopping tolerance.
      integer, intent(in), optional :: maxit
         !! Maximum number of iterations allowed.
      integer, intent(in), optional :: iprint
         !! Print control variable.
      integer, intent(in), optional :: lunerr
         !! Logical unit number for error messages.
      integer, intent(in), optional :: lunrpt
         !! Logical unit number for computation reports.
      real(kind=wp), intent(in), optional :: stpb(:)
         !! Relative step for computing finite difference derivatives with respect to `beta`.
         !! `Shape: (np)`.
      real(kind=wp), intent(in), optional :: stpd(:, :)
         !! Relative step for computing finite difference derivatives with respect to `delta`.
         !! `Shape: (1<=ldstpd<=n, m)`. See p. 31.
      real(kind=wp), intent(in), optional :: sclb(:)
      !! Scaling values for `beta`. `Shape: (np)`.
      real(kind=wp), intent(in), optional :: scld(:, :)
         !! Scaling values for `delta`. `Shape: (1<=ldscld<=n, m)`. See p. 32.
      real(kind=wp), intent(inout), pointer, optional :: work(:)
         !! Real work space.
      integer, intent(out), pointer, optional :: iwork(:)
         !! Integer work space.
      integer, intent(out), optional :: info
         !! Variable designating why the computations were stopped.
      real(kind=wp), intent(in), optional :: lower(:)
         !! Lower bound on `beta`. `Shape: (np)`.
      real(kind=wp), intent(in), optional :: upper(:)
         !! Upper bound on `beta`. `Shape: (np)`.

      ! Local variables
      ! TODO: remove save? replace pointer with allocatable?
      integer :: ldwe, ld2we, ldwd, ld2wd, ldifx, ldscld, ldstpd, ljob, lndigit, lmaxit, &
                 liprint, llunerr, llunrpt, linfo, lenwork, leniwork, linfo1, linfo2, linfo3, &
                 linfo4, linfo5
      integer :: lifixb(np), lifixx(n, m)
      real(kind=wp) :: ltaufac, lsstol, lpartol
      real(kind=wp) :: llower(np), lwe(n, nq, nq), lwd(n, m, m), lstpb(np), lstpd(n, m), &
                       lsclb(np), lscld(n, m), lupper(np), wd1(1, 1, 1)

      real(kind=wp), pointer, save :: lwork(:)
      integer, pointer, save :: liwork(:)  
      logical :: head

      ! Set LINFO to zero indicating no errors have been found thus far
      linfo = 0
      linfo1 = 0
      linfo2 = 0
      linfo3 = 0
      linfo4 = 0
      linfo5 = 0

      ! Set all scalar variable defaults except JOB
      ldwe = 1
      ld2we = 1
      ldwd = 1
      ld2wd = 1
      ldifx = 1
      ldscld = 1
      ldstpd = 1
      liprint = -1
      llunerr = -1
      llunrpt = -1
      lmaxit = -1
      lndigit = -1
      lpartol = negone
      lsstol = negone
      ltaufac = negone
      head = .true.

      !  Check for the option arguments for printing (so error messages can be
      !  printed appropriately from here on out
      if (present(iprint)) then
         liprint = iprint
      end if

      if (present(lunrpt)) then
         llunrpt = lunrpt
      end if
      if (llunrpt .lt. 0) then
         llunrpt = 6
      end if

      if (present(lunerr)) then
         llunerr = lunerr
      end if
      if (llunerr .lt. 0) then
         llunerr = 6
      end if

      ! Ensure the problem size is valid
      if (n .le. 0) then
         linfo5 = 1
         linfo4 = 1
      end if

      if (m .le. 0) then
         linfo5 = 1
         linfo3 = 1
      end if

      if (np .le. 0) then
         linfo5 = 1
         linfo2 = 1
      end if

      if (nq .le. 0) then
         linfo5 = 1
         linfo1 = 1
      end if

      if (linfo5 .ne. 0) then
         linfo = 10000*linfo5 + 1000*linfo4 + 100*linfo3 + 10*linfo2 + linfo1
         if (llunerr .gt. 0 .and. liprint .ne. 0) then
            call dodphd(head, llunrpt)
            call dodpe1( &
               llunerr, linfo, linfo5, linfo4, linfo3, linfo2, linfo1, &
               n, m, nq, &
               ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
               lenwork, leniwork &
               )
         end if
         if (present(info)) then
            info = linfo
         end if
         return
      end if

      ! Define LJOB and check that necessary arguments are passed for JOB
      if (present(job)) then
         ljob = job
         if (mod(job, 10000)/1000 .ge. 1) then
            if (.not. present(delta)) then
               linfo5 = 7
               linfo4 = 1
            ! elseif (.not. allocated(delta)) then
            !    linfo5 = 7
            !    linfo4 = 1
            end if
         end if
         if (job .ge. 10000) then
            if (.not. present(iwork)) then
               linfo5 = 7
               linfo2 = 1
            elseif (.not. associated(iwork)) then
               linfo5 = 7
               linfo2 = 1
            end if
         end if
         if (job .ge. 10000) then
            if (.not. present(work)) then
               linfo5 = 7
               linfo3 = 1
            elseif (.not. associated(work)) then
               linfo5 = 7
               linfo3 = 1
            end if
         end if
      else
         ljob = -1
      end if

      if (linfo5 .ne. 0) then
         linfo = 10000*linfo5 + 1000*linfo4 + 100*linfo3 + 10*linfo2 + linfo1
         if (llunerr .gt. 0 .and. liprint .ne. 0) then
            call dodphd(head, llunrpt)
            call dodpe1( &
               llunerr, linfo, linfo5, linfo4, linfo3, linfo2, linfo1, &
               n, m, nq, &
               ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
               lenwork, leniwork &
               )
         end if
         if (present(info)) then
            info = linfo
         end if
         return
      end if

      ! Determine the size of WORK
      if (ljob .lt. 0 .or. mod(ljob, 10) .le. 1) then
         lenwork = 18 + 13*np + np**2 + m + m**2 + 4*n*nq + 6*n*m + 2*n*nq*np + &
                   2*n*nq*m + nq**2 + 5*nq + nq*(np + m) + n*nq*nq
      else
         lenwork = 18 + 13*np + np**2 + m + m**2 + 4*n*nq + 2*n*m + 2*n*nq*np + &
                   5*nq + nq*(np + m) + n*nq*nq
      end if

      ! Determine the size of IWORK
      leniwork = 20 + 2*np + nq*(np + m)

      ! Allocate the work arrays
      allocate (lwork(lenwork), tempret(max(n, np), max(nq, m)), stat=linfo3)
      allocate (liwork(leniwork), stat=linfo2)
      lwork = zero
      liwork = 0
      ! if (present(delta)) then
      !    if (.not. allocated(delta)) then
      !       allocate (delta(n, m), stat=linfo4)
      !    end if
      ! end if
      if (linfo4 .ne. 0 .or. linfo3 .ne. 0 .or. linfo2 .ne. 0) then
         linfo5 = 8
      end if

      if (linfo5 .ne. 0) then
         linfo = 10000*mod(linfo5, 10) + 1000*mod(linfo4, 10) + &
                 100*mod(linfo3, 10) + 10*mod(linfo2, 10) + mod(linfo1, 10)
         if (llunerr .gt. 0 .and. liprint .ne. 0) then
            call dodphd(head, llunrpt)
            call dodpe1( &
               llunerr, linfo, linfo5, linfo4, linfo3, linfo2, linfo1, &
               n, m, nq, &
               ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
               lenwork, leniwork)
         end if
         if (present(info)) then
            info = linfo
         end if
         return
      end if

      ! Set array variable defaults except IWORK
      lwork(1:n*m) = zero
      lifixb(1) = -1
      lifixx(1, 1) = -1
      llower(1:np) = -huge(zero)
      lsclb(1) = negone
      lscld(1, 1) = negone
      lstpb(1) = negone
      lstpd(1, 1) = negone
      lupper(1:np) = huge(zero)
      lwe(1, 1, 1) = negone
      lwd(1, 1, 1) = negone

      ! Check the size of required arguments and return errors if they are too small
      if (size(beta) .lt. np) then
         linfo1 = linfo1 + 1
      end if

      if (any(size(y) .lt. (/n, nq/))) then
         linfo1 = linfo1 + 2
      end if

      if (any(size(x) .lt. (/n, m/))) then
         linfo1 = linfo1 + 4
      end if

      ! Check the presence of optional arguments and copy their values internally or
      ! report errors as necessary
      if (present(ifixb)) then
         if (size(ifixb) .lt. np) then
            linfo1 = linfo1 + 64
         end if
         if (ifixb(1) .lt. 0) then
            lifixb(1) = ifixb(1)
         else
            lifixb(1:np) = ifixb(1:np)
         end if
      end if

      if (present(ifixx)) then
         ldifx = size(ifixx, 1)
         if (any(size(ifixx) .le. (/0, 0/))) then
            linfo1 = linfo1 + 128
         end if
         if (.not. (ifixx(1, 1) .lt. 0 .or. ldifx .eq. 1 .or. ldifx .ge. n) &
             .or. size(ifixx, 2) .lt. m) then
            linfo1 = linfo1 + 128
         end if
         if (ldifx .gt. n) then
            ldifx = n
         end if
         if (ifixx(1, 1) .lt. 0) then
            lifixx(1, 1) = ifixx(1, 1)
         else
            lifixx(1:ldifx, 1:m) = ifixx(1:ldifx, 1:m)
         end if
      end if

      if (present(iwork)) then
         if (associated(iwork)) then
            if (size(iwork) .lt. leniwork) then
               linfo1 = linfo1 + 8192
            end if
            !  This is a restart, copy IWORK.
            if (mod(ljob/10000, 10) .ge. 1) then
               liwork(1:leniwork) = iwork(1:leniwork)
            end if
         end if
      end if

      if (present(maxit)) then
         lmaxit = maxit
      end if

      if (present(ndigit)) then
         lndigit = ndigit
      end if

      if (present(partol)) then
         lpartol = partol
      end if

      if (present(sclb)) then
         if (size(sclb) .lt. np) then
            linfo1 = linfo1 + 1024
         end if
         if (sclb(1) .le. zero) then
            lsclb(1) = sclb(1)
         else
            lsclb(1:np) = sclb(1:np)
         end if
      end if

      if (present(scld)) then
         ldscld = size(scld, 1)
         if (any(size(scld) .le. (/0, 0/))) then
            linfo1 = linfo1 + 2048
         end if
         if (.not. (scld(1, 1) .le. zero .or. ldscld .eq. 1 .or. ldscld .ge. n) &
             .or. size(scld, 2) .lt. m) then
            linfo1 = linfo1 + 2048
         end if
         if (ldscld .gt. n) then
            ldscld = n
         end if
         if (scld(1, 1) .le. zero) then
            lscld(1, 1) = scld(1, 1)
         else
            lscld(1:ldscld, 1:m) = scld(1:ldscld, 1:m)
         end if
      end if

      if (present(sstol)) then
         lsstol = sstol
      end if

      if (present(stpb)) then
         if (size(stpb) .lt. np) then
            linfo1 = linfo1 + 256
         end if
         if (stpb(1) .le. zero) then
            lstpb(1) = stpb(1)
         else
            lstpb(1:np) = stpb(1:np)
         end if
      end if

      if (present(stpd)) then
         ldstpd = size(stpd, 1)
         if (any(size(stpd) .le. (/0, 0/))) then
            linfo1 = linfo1 + 512
         end if
         if (.not. (stpd(1, 1) .le. zero .or. ldstpd .eq. 1 .or. ldstpd .ge. n) &
             .or. size(stpd, 2) .lt. m) then
            linfo1 = linfo1 + 512
         end if
         if (ldstpd .gt. n) then
            ldstpd = n
         end if
         if (stpd(1, 1) .le. zero) then
            lstpd(1, 1) = stpd(1, 1)
         else
            lstpd(1:ldstpd, 1:m) = stpd(1:ldstpd, 1:m)
         end if
      end if

      if (present(taufac)) then
         ltaufac = taufac
      end if

      if (present(we)) then
         ldwe = size(we, 1)
         ld2we = size(we, 2)
         if (any(size(we) .le. (/0, 0, 0/))) then
            linfo1 = linfo1 + 16
         end if
         if (.not. (we(1, 1, 1) .lt. zero .or. &
                    ((ldwe .eq. 1 .or. ldwe .ge. n) .and. &
                     (ld2we .eq. 1 .or. ld2we .ge. nq))) .or. size(we, 3) .lt. nq) then
            linfo1 = linfo1 + 16
         end if
         if (ldwe .gt. n) then
            ldwe = n
         end if
         if (ld2we .gt. nq) then
            ld2we = nq
         end if
         if (we(1, 1, 1) .lt. zero) then
            lwe(1, 1, 1) = we(1, 1, 1)
         else
            lwe(1:ldwe, 1:ld2we, 1:nq) = we(1:ldwe, 1:ld2we, 1:nq)
         end if
      end if

      if (present(wd)) then
         ldwd = size(wd, 1)
         ld2wd = size(wd, 2)
         if (any(size(wd) .le. (/0, 0, 0/))) then
            linfo1 = linfo1 + 32
         end if
         if (.not. (wd(1, 1, 1) .lt. zero .or. &
                    ((ldwd .eq. 1 .or. ldwd .ge. n) .and. &
                     (ld2wd .eq. 1 .or. ld2wd .ge. m))) .or. size(wd, 3) .lt. m) then
            linfo1 = linfo1 + 32
         end if
         if (ldwd .gt. n) then
            ldwd = n
         end if
         if (ld2wd .gt. m) then
            ld2wd = m
         end if
         if (wd(1, 1, 1) .le. 0.0_wp) then
            lwd(1, 1, 1) = wd(1, 1, 1)
         else
            lwd(1:ldwd, 1:ld2wd, 1:m) = wd(1:ldwd, 1:ld2wd, 1:m)
         end if
      end if

      if (present(work)) then
         if (associated(work)) then
            if (size(work) .lt. lenwork) then
               linfo1 = linfo1 + 4096
            end if
            !  Deltas are in WORK, copy them.
            if (mod(ljob/1000, 10) .ge. 1 .and. .not. present(delta)) then
               lwork(1:n*m) = work(1:n*m)
            end if
            !  This is a restart, copy WORK.
            if (mod(ljob/10000, 10) .ge. 1) then
               lwork(1:lenwork) = work(1:lenwork)
            end if
         end if
      end if

      if (present(delta)) then
         !if (allocated(delta)) then
            if (any(shape(delta) .lt. (/n, m/))) then
               linfo1 = linfo1 + 8
            end if
            lwork(1:n*m) = reshape(delta(1:n, 1:m), (/n*m/))
         !end if
      end if

      if (present(lower)) then
         if (size(lower) .lt. np) then
            linfo1 = linfo1 + 32768
         end if
         llower(1:np) = lower(1:np)
      end if

      if (present(upper)) then
         if (size(upper) .lt. np) then
            linfo1 = linfo1 + 16384
         end if
         lupper(1:np) = upper(1:np)
      end if

      ! Report an error if any of the array sizes didn't match.
      if (linfo1 .ne. 0) then
         linfo = 100000 + linfo1
         linfo1 = 0
         if (llunerr .gt. 0 .and. liprint .ne. 0) then
            call dodphd(head, llunrpt)
            call dodpe1( &
               llunerr, linfo, linfo5, linfo4, linfo3, linfo2, linfo1, &
               n, m, nq, &
               ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
               lenwork, leniwork &
               )
         end if
         if (present(info)) then
            info = linfo
         end if
         return
      end if

      if (lwd(1, 1, 1) .ne. 0) then
         call dodcnt &
            (fcn, &
             n, m, np, nq, &
             beta(1:np), &
             y(1:n, 1:nq), n, x(1:n, 1:m), n, &
             lwe(1:ldwe, 1:ld2we, 1:nq), ldwe, ld2we, &
             lwd(1:ldwd, 1:ld2wd, 1:m), ldwd, ld2wd, &
             lifixb, lifixx(1:ldifx, 1:m), ldifx, &
             ljob, lndigit, ltaufac, &
             lsstol, lpartol, lmaxit, &
             liprint, llunerr, llunrpt, &
             lstpb, lstpd(1:ldstpd, 1:m), ldstpd, &
             lsclb, lscld(1:ldscld, 1:m), ldscld, &
             lwork, lenwork, liwork, leniwork, &
             linfo, &
             llower, lupper)
      else
         wd1(1, 1, 1) = negone
         call dodcnt &
            (fcn, &
             n, m, np, nq, &
             beta(1:np), &
             y(1:n, 1:nq), n, x(1:n, 1:m), n, &
             lwe(1:ldwe, 1:ld2we, 1:nq), ldwe, ld2we, &
             wd1, 1, 1, &
             lifixb, lifixx(1:ldifx, 1:m), ldifx, &
             ljob, lndigit, ltaufac, &
             lsstol, lpartol, lmaxit, &
             liprint, llunerr, llunrpt, &
             lstpb, lstpd(1:ldstpd, 1:m), ldstpd, &
             lsclb, lscld(1:ldscld, 1:m), ldscld, &
             lwork, lenwork, liwork, leniwork, &
             linfo, &
             llower, lupper)
      end if

      if (present(delta)) then
         !if (allocated(delta)) then
            delta(1:n, 1:m) = reshape(lwork(1:n*m), (/n, m/))
         !end if
      end if

      if (present(info)) then
         info = linfo
      end if

      if (present(iwork)) then
         if (.not. associated(iwork)) then
            iwork => liwork
         else
            iwork(1:leniwork) = liwork(1:leniwork)
            deallocate (liwork)
         end if
      else
         deallocate (liwork)
      end if

      if (present(work)) then
         if (.not. associated(work)) then
            work => lwork
         else
            work(1:lenwork) = lwork(1:lenwork)
            deallocate (lwork)
         end if
      else
         deallocate (lwork)
      end if

      deallocate (tempret)

   end subroutine odr

   subroutine dodcnt &
      (fcn, n, m, np, nq, beta, y, ldy, x, ldx, &
       we, ldwe, ld2we, wd, ldwd, ld2wd, ifixb, ifixx, ldifx, &
       job, ndigit, taufac, sstol, partol, maxit, iprint, lunerr, lunrpt, &
       stpb, stpd, ldstpd, sclb, scld, ldscld, &
       work, lwork, iwork, liwork, &
       info, &
       lower, upper)
   !! Driver routine for finding the weighted explicit or implicit orthogonal distance
   !! regression (ODR) or ordinary linear or nonlinear least squares (OLS) solution.
   ! Routines Called  DODDRV
   ! Date Written   860529   (YYMMDD)
   ! Revision Date  920304   (YYMMDD)
   
      use odrpack_kinds, only: zero, one, three
   
      procedure(fcn_t) :: fcn
         !! User-supplied subroutine for evaluating the model.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the independent variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(kind=wp), intent(inout) :: beta(np)
         !! The function parameters.
      real(kind=wp), intent(in) :: y(ldy, nq)
         !! The dependent variable. Unused when the model is implicit.
      integer, intent(in) :: ldy
         !! The leading dimension of array `y`.
      real(kind=wp), intent(in) :: x(ldx, m)
         !! The independent variable.
      integer, intent(in) :: ldx
         !! The leading dimension of array `x`.
      real(kind=wp), intent(in) :: we(ldwe, ld2we, nq)
         !! The `epsilon` weights.
      integer, intent(in) :: ldwe
         !! The leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! The second dimension of array `we`.
      real(kind=wp), intent(in) :: wd(ldwd, ld2wd, m)
         !! The `delta` weights.
      integer, intent(in) :: ldwd
         !! The leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! The second dimension of array `wd`.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input
         !! values or not.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values
         !! or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      integer, intent(in) :: job
         !! The variable controlling problem initialization and computational method.
      integer, intent(in) :: ndigit
         !! The number of accurate digits in the function results, as supplied by the user.
      real(kind=wp), intent(in) :: taufac
         !! The factor used to compute the initial trust region diameter.
      real(kind=wp), intent(in) :: sstol
         !! The sum-of-squares convergence stopping tolerance.
      real(kind=wp), intent(in) :: partol
         !! The user-supplied parameter convergence stopping tolerance.
      integer, intent(in) :: maxit
         !! The maximum number of iterations allowed.
      integer, intent(in) :: iprint
         !! The print control variables.
      integer, intent(in) :: lunerr
         !! The logical unit number used for error messages.
      integer, intent(in) :: lunrpt
         !! The logical unit number used for computation reports.
      real(kind=wp), intent(in) :: stpb(np)
         !! The relative step for computing finite difference derivatives with respect to `beta`.
      real(kind=wp), intent(in) :: stpd(ldstpd, m)
         !! The relative step for computing finite difference derivatives with respect to `delta`.
      integer, intent(in) :: ldstpd
         !! The leading dimension of array `stpd`.
      real(kind=wp), intent(in) :: sclb(np)
         !! The scaling values for `beta`.
      real(kind=wp), intent(in) :: scld(ldscld, m)
         !! The scaling value for `delta`.
      integer, intent(in) :: ldscld
         !! The leading dimension of array `scld`.
      real(kind=wp), intent(inout) :: work(lwork)
         !! The real work space.
      integer, intent(in) :: lwork
         !! The length of vector `work`.
      integer, intent(inout) :: iwork(liwork)
         !! The integer work space.
      integer, intent(in) :: liwork
         !! The length of vector `iwork`.
      integer, intent(out) :: info
         !! The variable designating why the computations were stopped.
      real(kind=wp), intent(in) :: lower(np)
         !! The lower bound on `beta`.
      real(kind=wp), intent(in) :: upper(np)
         !! The upper bound on `beta`.
   
      ! Local scalars
      real(kind=wp), parameter :: pcheck = 1.0E3_wp, pstart = 1.0E1_wp, pfac = 1.0E1_wp
      real(kind=wp) :: cnvtol, tstimp
      integer :: iprnti, ipr1, ipr2, ipr2f, ipr3, jobi, job1, job2, job3, job4, job5, &
                 maxiti, maxit1
      logical :: done, fstitr, head, implct, prtpen
   
      ! Local arrays
      real(kind=wp) :: pnlty(1, 1, 1)
   
      ! External subroutines
      external :: doddrv
   
      ! Variable Definitions (alphabetically)
      !  BETA:    The function parameters.
      !  CNVTOL:  The convergence tolerance for implicit models.
      !  DONE:    The variable designating whether the inplicit solution has been found (DONE=TRUE)
      !           or not (DONE=FALSE).
      !  FCN:     The user-supplied subroutine for evaluating the model.
      !  FSTITR:  The variable designating whether this is the first iteration (FSTITR=TRUE)
      !           or not (FSTITR=FALSE).
      !  HEAD:    The variable designating whether the heading is to be printed (HEAD=TRUE)
      !           or not (HEAD=FALSE).
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input
      !           values or not.
      !  IFIXX:   The values designating whether the elements of X are fixed at their input
      !           values or not.
      !  IMPLCT:  The variable designating whether the solution is by implicit ODR (IMPLCT=TRUE)
      !           or explicit ODR (IMPLCT=FALSE).
      !  INFO:    The variable designating why the computations were stopped.
      !  IPRINT:  The print control variables.
      !  IPRNTI:  The print control variables.
      !  IPR1:    The 1st digit of the print control variable.
      !  IPR2:    The 2nd digit of the print control variable.
      !  IPR3:    The 3rd digit of the print control variable.
      !  IPR4:    The 4th digit of the print control variable.
      !  IWORK:   The integer work space.
      !  JOB:     The variable controling problem initialization and computational method.
      !  JOBI:    The variable controling problem initialization and computational method.
      !  JOB1:    The 1st digit of the variable JOB.
      !  JOB2:    The 2nd digit of the variable JOB.
      !  JOB3:    The 3rd digit of the variable JOB.
      !  JOB4:    The 4th digit of the variable JOB.
      !  JOB5:    The 5th digit of the variable JOB.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LDSCLD:  The leading dimension of array SCLD.
      !  LDSTPD:  The leading dimension of array STPD.
      !  LDWD:    The leading dimension of array WD.
      !  LDWE:    The leading dimension of array WE.
      !  LDX:     The leading dimension of array X.
      !  LDY:     The leading dimension of array Y.
      !  LD2WD:   The second dimension of array WD.
      !  LD2WE:   The second dimension of array WE.
      !  LIWORK:  The length of vector IWORK.
      !  LOWER:   The lower bound for BETA.
      !  LUNERR:  The logical unit number used for error messages.
      !  LUNRPT:  The logical unit number used for computation reports.
      !  LWORK:   The length of vector work.
      !  M:       The number of columns of data in the independent variable.
      !  MAXIT:   The maximum number of iterations allowed.
      !  MAXITI:  For implicit models, the number of iterations allowed for the current penalty
      !           parameter value.
      !  MAXIT1:  For implicit models, the number of iterations allowed for the next penalty
      !           parameter value.
      !  N:       The number of observations.
      !  NDIGIT:  The number of accurate digits in the function results, as supplied by the user.
      !  NP:      The number of function parameters.
      !  NQ:      The number of responses per observation.
      !  PARTOL:  The user supplied parameter convergence stopping tolerance.
      !  PCHECK:  The value designating the minimum penalty parameter allowed before the implicit
      !           problem can be considered solved.
      !  PFAC:    The factor for increasing the penalty parameter.
      !  PNLTY:   The penalty parameter for an implicit model.
      !  PRTPEN:  The value designating whether the penalty parameter is to be printed in the
      !           iteration report (PRTPEN=TRUE) or not (PRTPEN=FALSE).
      !  PSTART:  The factor for increasing the penalty parameter.
      !  SCLB:    The scaling values for BETA.
      !  SCLD:    The scaling values for DELTA.
      !  STPB:    The relative step for computing finite difference derivatives with respect to BETA.
      !  STPD:    The relative step for computing finite difference derivatives with respect to DELTA.
      !  SSTOL:   The sum-of-squares convergence stopping tolerance.
      !  TAUFAC:  The factor used to compute the initial trust region diameter.
      !  TSTIMP:  The relative change in the parameters between the initial values and the solution.
      !  UPPER:   The upper bound for BETA.
      !  WD:      The DELTA weights.
      !  WE:      The EPSILON weights.
      !  WORK:    The real work space.
      !  X:       The independent variable.
      !  Y:       The dependent variable. Unused when the model is implicit.
   
      implct = mod(job, 10) .eq. 1
      fstitr = .true.
      head = .true.
      prtpen = .false.
   
      if (implct) then
         !  Set up for implicit problem
         if (iprint .ge. 0) then
            ipr1 = mod(iprint, 10000)/1000
            ipr2 = mod(iprint, 1000)/100
            ipr2f = mod(iprint, 100)/10
            ipr3 = mod(iprint, 10)
         else
            ipr1 = 2
            ipr2 = 0
            ipr2f = 0
            ipr3 = 1
         end if
         iprnti = ipr1*1000 + ipr2*100 + ipr2f*10
   
         job5 = mod(job, 100000)/10000
         job4 = mod(job, 10000)/1000
         job3 = mod(job, 1000)/100
         job2 = mod(job, 100)/10
         job1 = mod(job, 10)
         jobi = job5*10000 + job4*1000 + job3*100 + job2*10 + job1
   
         if (we(1, 1, 1) .le. zero) then
            pnlty(1, 1, 1) = -pstart
         else
            pnlty(1, 1, 1) = -we(1, 1, 1)
         end if
   
         if (partol .lt. zero) then
            cnvtol = epsilon(zero)**(one/three)
         else
            cnvtol = min(partol, one)
         end if
   
         if (maxit .ge. 1) then
            maxiti = maxit
         else
            maxiti = 100
         end if
   
         done = maxiti .eq. 0
         prtpen = .true.
   
         do while (.true.)
            call doddrv &
               (head, fstitr, prtpen, &
                fcn, n, m, np, nq, beta, y, ldy, x, ldx, &
                pnlty, 1, 1, wd, ldwd, ld2wd, ifixb, ifixx, ldifx, &
                jobi, ndigit, taufac, sstol, cnvtol, maxiti, &
                iprnti, lunerr, lunrpt, &
                stpb, stpd, ldstpd, sclb, scld, ldscld, &
                work, lwork, iwork, liwork, &
                maxit1, tstimp, info, lower, upper)
   
            if (done) then
               return
            else
               done = maxit1 .le. 0 .or. (abs(pnlty(1, 1, 1)) .ge. pcheck .and. tstimp .le. cnvtol)
            end if
   
            if (done) then
               if (tstimp .le. cnvtol) then
                  info = (info/10)*10 + 2
               else
                  info = (info/10)*10 + 4
               end if
               jobi = 10000 + 1000 + job3*100 + job2*10 + job1
               maxiti = 0
               iprnti = ipr3
            else
               prtpen = .true.
               pnlty(1, 1, 1) = pfac*pnlty(1, 1, 1)
               jobi = 10000 + 1000 + 000 + job2*10 + job1
               maxiti = maxit1
               iprnti = 0000 + ipr2*100 + ipr2f*10
            end if
         end do
      else
         ! Explicit problem
         call doddrv &
            (head, fstitr, prtpen, &
             fcn, n, m, np, nq, beta, y, ldy, x, ldx, &
             we, ldwe, ld2we, wd, ldwd, ld2wd, ifixb, ifixx, ldifx, &
             job, ndigit, taufac, sstol, partol, maxit, &
             iprint, lunerr, lunrpt, &
             stpb, stpd, ldstpd, sclb, scld, ldscld, &
             work, lwork, iwork, liwork, &
             maxit1, tstimp, info, lower, upper)
      end if
   
   end subroutine dodcnt

end module odrpack
