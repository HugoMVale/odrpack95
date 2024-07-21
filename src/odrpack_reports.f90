module odrpack_reports
   !! Reporting routines.

   use odrpack_kinds, only: wp
   implicit none

contains

   impure subroutine dodpc1 &
      (ipr, lunrpt, &
      anajac, cdjac, chkjac, initd, restrt, isodr, implct, dovcv, redoj, &
      msgb1, msgb, msgd1, msgd, &
      n, m, np, nq, npp, nnzw, &
      x, ldx, ifixx, ldifx, delta, wd, ldwd, ld2wd, tt, ldtt, stpd, ldstpd &
      , y, ldy, we, ldwe, ld2we, pnlty, beta, ifixb, ssf, stpb, lower, &
      upper, job, neta, taufac, sstol, partol, maxit, wss, wssdel, wsseps)
   !! Generate initial summary report.
   ! Routines Called  DHSTEP
   ! Date Written   860529   (YYMMDD)
   ! Revision Date  920619   (YYMMDD)

      use odrpack_kinds, only: zero
      use odrpack_core, only: dhstep

      integer, intent(in) :: ipr
         !! The value indicating the report to be printed.
      integer, intent(in) :: lunrpt
         !! The logical unit number for the computation reports.
      logical, intent(in) :: anajac
         !! The variable designating whether the Jacobians are computed by finite differences
         !! (`anajac = .false.`) or not (`anajac = .true.`).
      logical, intent(in) :: cdjac
         !! The variable designating whether the Jacobians are computed by central differences
         !! (`cdjac = .true.`) or forward differences (`cdjac = .false.`).
      logical, intent(in) :: chkjac
         !! The variable designating whether the user-supplied Jacobians are to be checked
         !! (`chkjac = .true.`) or not (`chkjac = .false.`).
      logical, intent(in) :: initd
         !! The variable designating whether `delta` is initialized to zero (`initd = .true.`) 
         !! or to the values in the first `n` by `m` elements of array `work` (`initd = .false.`).
      logical, intent(in) :: restrt
         !! The variable designating whether the call is a restart (`restrt = .true.`) or
         !! not (`restrt = .false.`).
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      logical, intent(in) :: implct
         !! The variable designating whether the solution is by implicit ODR (`implct = .true.`)
         !! or explicit ODR (`implct = .false.`).
      logical, intent(in) :: dovcv
         !! The variable designating whether the covariance matrix is to be computed
         !! (`dovcv = .true.`) or not (`dovcv = .false.`).
      logical, intent(in) :: redoj
         !! The variable designating whether the Jacobian matrix is to be recomputed for the
         !! computation of the covariance matrix (`redoj = .true.`) or not (`redoj = .false.`).
      integer, intent(in) :: msgb1
         !! The error checking results for the Jacobian with respect to `beta`.
      integer, intent(in) :: msgb(nq, np)
         !! The error checking results for the Jacobian with respect to `beta`.
      integer, intent(in) :: msgd1
         !! The error checking results for the Jacobian with respect to `delta`.
      integer, intent(in) :: msgd(nq, m)
         !! The error checking results for the Jacobian with respect to `delta`.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(in) :: npp
         !! The number of function parameters being estimated.
      integer, intent(in) :: nnzw
         !! The number of nonzero observational error weights.
      real(wp), intent(in) :: x(ldx, m)
         !! The explanatory variable.
      integer, intent(in) :: ldx
         !! The leading dimension of array `x`.
      integer, intent(in) :: ifixx(ldifx, m)
         !! The values designating whether the elements of `x` are fixed at their input values or not.
      integer, intent(in) :: ldifx
         !! The leading dimension of array `ifixx`.
      real(wp), intent(in) :: delta(n, m)
         !! The estimated errors in the explanatory variables.
      real(wp), intent(in) :: wd(ldwd, ld2wd, m)
         !! The `delta` weights.
      integer, intent(in) :: ldwd
         !! The leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! The second dimension of array `wd`.
      real(wp), intent(in) :: tt(ldtt, m)
         !! The scaling values for `delta`.
      integer, intent(in) :: ldtt
         !! The leading dimension of array `tt`.
      real(wp), intent(in) :: stpd(ldstpd, m)
         !! The relative step used for computing finite difference derivatives with respect to `delta`.
      integer, intent(in) :: ldstpd
         !! The leading dimension of array `stpd`.
      real(wp), intent(in) :: y(ldy, nq)
         !! The response variable. Unused when the model is implicit.
      integer, intent(in) :: ldy
         !! The leading dimension of array `y`.
      real(wp), intent(in) :: we(ldwe, ld2we, nq)
         !! The `epsilon` weights.
      integer, intent(in) :: ldwe
         !! The leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! The second dimension of array `we`.
      real(wp), intent(in) :: pnlty
         !! The penalty parameter for an implicit model.
      real(wp), intent(in) :: beta(np)
         !! The function parameters.
      integer, intent(in) :: ifixb(np)
         !! The values designating whether the elements of `beta` are fixed at their input values
         !! or not.
      real(wp), intent(in) :: ssf(np)
         !! The scaling values for `beta`.
      real(wp), intent(in) :: stpb(np)
         !! The relative step used for computing finite difference derivatives with respect to `beta`.
      real(wp), intent(in) :: lower(np)
         !! The lower bounds for `beta`.
      real(wp), intent(in) :: upper(np)
         !! The upper bounds for `beta`.
      integer, intent(in) :: job
         !! The variable controlling problem initialization and computational method.
      integer, intent(in) :: neta
         !! The number of accurate digits in the function results. A negative value indicates
         !! that `neta` was estimated by 'odrpack'. A positive value indicates the value was
         !! supplied by the user.
      real(wp), intent(in) :: taufac
         !! The factor used to compute the initial trust region diameter.
      real(wp), intent(in) :: sstol
         !! The sum-of-squares convergence stopping tolerance.
      real(wp), intent(in) :: partol
         !! The parameter convergence stopping tolerance.
      integer, intent(in) :: maxit
         !! The maximum number of iterations allowed.
      real(wp), intent(in) :: wss
         !! The sum-of-squares of the weighted `epsilon`s and `delta`s.
      real(wp), intent(in) :: wssdel
         !! The sum-of-squares of the weighted `delta`s.
      real(wp), intent(in) :: wsseps
         !! The sum-of-squares of the weighted `epsilon`s.

      ! Local scalars
      real(wp) :: temp1, temp2, temp3
      integer :: i, itemp, j, job1, job2, job3, job4, job5, l

      ! Local arrays
      character(len=2) :: tempc0
      character(len=5) :: tempc1
      character(len=13) :: tempc2

      ! Variable Definitions (alphabetically)
      !  ANAJAC:  The variable designating whether the Jacobians are computed by finite differences
      !           (ANAJAC=FALSE) or not (ANAJAC=TRUE).
      !  BETA:    The function parameters.
      !  CDJAC:   The variable designating whether the Jacobians are computed by central differences
      !           (CDJAC=TRUE) or forward differences (CDJAC=FALSE).
      !  CHKJAC:  The variable designating whether the user supplied Jacobians are to be checked
      !           (CHKJAC=TRUE) or not (CHKJAC=FALSE).
      !  DELTA:   The estimated errors in the explanatory variables.
      !  DOVCV:   The variable designating whether the covariance matrix is to be computed
      !           (DOVCV=TRUE) or not (DOVCV=FALSE).
      !  I:       An indexing variable.
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input
      !           values or not.
      !  IFIXX:   The values designating whether the elements of X are fixed at their input
      !           values or not.
      !  IMPLCT:  The variable designating whether the solution is by implicit ODR (IMPLCT=TRUE)
      !           or explicit ODR (IMPLCT=FALSE).
      !  INITD:   The variable designating whether DELTA is initialized to zero (INITD=TRUE) or
      !           to the values in the first N by M elements of array WORK (INITD=FALSE).
      !  IPR:     The value indicating the report to be printed.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or by
      !           OLS (ISODR=FALSE).
      !  ITEMP:   A temporary integer value.
      !  J:       An indexing variable.
      !  JOB:     The variable controling problem initialization and computational method.
      !  JOB1:    The 1st digit (from the left) of variable JOB.
      !  JOB2:    The 2nd digit (from the left) of variable JOB.
      !  JOB3:    The 3rd digit (from the left) of variable JOB.
      !  JOB4:    The 4th digit (from the left) of variable JOB.
      !  JOB5:    The 5th digit (from the left) of variable JOB.
      !  L:       An indexing variable.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LDTT:    The leading dimension of array TT.
      !  LDWD:    The leading dimension of array WD.
      !  LDWE:    The leading dimension of array WE.
      !  LDX:     The leading dimension of array X.
      !  LDY:     The leading dimension of array Y.
      !  LD2WD:   The second dimension of array WD.
      !  LD2WE:   The second dimension of array WE.
      !  LUNRPT:  The logical unit number for the computation reports.
      !  M:       The number of columns of data in the explanatory variable.
      !  MAXIT:   The maximum number of iterations allowed.
      !  MSGB:    The error checking results for the Jacobian wrt beta.
      !  MSGB1:   The error checking results for the Jacobian wrt BETA.
      !  MSGD:    The error checking results for the Jacobian wrt DELTA.
      !  MSGD1:   The error checking results for the Jacobian wrt DELTA.
      !  N:       The number of observations.
      !  NETA:    The number of accurate digits in the function results. A negative value
      !           indicates that NETA was estimated by ODRPACK. A positive value indictes the
      !           value was supplied by the user.
      !  NNZW:    The number of nonzero observational error weights.
      !  NP:      The number of function parameters.
      !  NPP:     The number of function parameters being estimated.
      !  NQ:      The number of responses per observation.
      !  PARTOL:  The parameter convergence stopping tolerance.
      !  PNLTY:   The penalty parameter for an implicit model.
      !  REDOJ:   The variable designating whether the Jacobian matrix is to be recomputed for 
      !           the computation of the covariance matrix (REDOJ=TRUE) or not (REDOJ=FALSE).
      !  RESTRT:  The variable designating whether the call is a restart (RESTRT=TRUE) or
      !           not (RESTRT=FALSE).
      !  SSF:     The scaling values for BETA.
      !  SSTOL:   The sum-of-squares convergence stopping tolerance.
      !  STPB:    The relative step used for computing finite difference derivatives with respect
      !           to BETA.
      !  STPD:    The relative step used for computing finite difference derivatives with respect
      !           to DELTA.
      !  TAUFAC:  The factor used to compute the initial trust region diameter.
      !  TEMPC0:  A temporary CHARACTER*2 value.
      !  TEMPC1:  A temporary CHARACTER*5 value.
      !  TEMPC2:  A temporary CHARACTER*13 value.
      !  TEMP1:   A temporary REAL (wp) value.
      !  TEMP2:   A temporary REAL (wp) value.
      !  TEMP3:   A temporary REAL (wp) value.
      !  TT:      The scaling values for DELTA.
      !  WD:      The DELTA weights.
      !  WE:      The EPSILON weights.
      !  WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS.
      !  WSSDEL:  The sum-of-squares of the weighted DELTAS.
      !  WSSEPS:  The sum-of-squares of the weighted EPSILONS.
      !  X:       The explanatory variable.
      !  Y:       The response variable.  Unused when the model is implicit.

      ! Print problem size specification
      write (lunrpt, 1000) n, nnzw, nq, m, np, npp

      ! Print control values
      job1 = job/10000
      job2 = mod(job, 10000)/1000
      job3 = mod(job, 1000)/100
      job4 = mod(job, 100)/10
      job5 = mod(job, 10)
      write (lunrpt, 1100) job
      if (restrt) then
         write (lunrpt, 1110) job1
      else
         write (lunrpt, 1111) job1
      end if
      if (isodr) then
         if (initd) then
            write (lunrpt, 1120) job2
         else
            write (lunrpt, 1121) job2
         end if
      else
         write (lunrpt, 1122) job2, job5
      end if
      if (dovcv) then
         write (lunrpt, 1130) job3
         if (redoj) then
            write (lunrpt, 1131)
         else
            write (lunrpt, 1132)
         end if
      else
         write (lunrpt, 1133) job3
      end if
      if (anajac) then
         write (lunrpt, 1140) job4
         if (chkjac) then
            if (msgb1 >= 1 .or. msgd1 >= 1) then
               write (lunrpt, 1141)
            else
               write (lunrpt, 1142)
            end if
         else
            write (lunrpt, 1143)
         end if
      elseif (cdjac) then
         write (lunrpt, 1144) job4
      else
         write (lunrpt, 1145) job4
      end if
      if (isodr) then
         if (implct) then
            write (lunrpt, 1150) job5
         else
            write (lunrpt, 1151) job5
         end if
      else
         write (lunrpt, 1152) job5
      end if
      if (neta < 0) then
         write (lunrpt, 1200) - neta
      else
         write (lunrpt, 1210) neta
      end if
      write (lunrpt, 1300) taufac

      ! Print stopping criteria
      write (lunrpt, 1400) sstol, partol, maxit

      ! Print initial sum of squares
      if (implct) then
         write (lunrpt, 1500) wssdel
         if (isodr) then
            write (lunrpt, 1510) wss, wsseps, pnlty
         end if
      else
         write (lunrpt, 1600) wss
         if (isodr) then
            write (lunrpt, 1610) wssdel, wsseps
         end if
      end if

      if (ipr >= 2) then

         ! Print function parameter data
         write (lunrpt, 4000)
         if (chkjac .and. &
            ((msgb1 >= 1) .or. &
            (msgd1 >= 1))) then
            write (lunrpt, 4110)
         elseif (anajac) then
            write (lunrpt, 4120)
         else
            write (lunrpt, 4200)
         end if
         do j = 1, np
            if (ifixb(1) < 0) then
               tempc1 = '   NO'
            else
               if (ifixb(j) /= 0) then
                  tempc1 = '   NO'
               else
                  tempc1 = '  YES'
               end if
            end if
            if (anajac) then
               if (chkjac .and. &
                  ((msgb1 >= 1) .or. &
                  (msgd1 >= 1))) then
                  itemp = -1
                  do l = 1, nq
                     itemp = max(itemp, msgb(l, j))
                  end do
                  if (itemp <= -1) then
                     tempc2 = '    UNCHECKED'
                  elseif (itemp == 0) then
                     tempc2 = '     VERIFIED'
                  elseif (itemp >= 1) then
                     tempc2 = ' QUESTIONABLE'
                  end if
               else
                  tempc2 = '             '
               end if
            else
               tempc2 = '             '
            end if
            if (ssf(1) < zero) then
               temp1 = abs(ssf(1))
            else
               temp1 = ssf(j)
            end if
            if (anajac) then
               write (lunrpt, 4310) j, beta(j), tempc1, temp1, lower(j), &
                  upper(j), tempc2
            else
               if (cdjac) then
                  temp2 = dhstep(1, neta, 1, j, stpb, 1)
               else
                  temp2 = dhstep(0, neta, 1, j, stpb, 1)
               end if
               write (lunrpt, 4320) j, beta(j), tempc1, temp1, &
                  lower(j), upper(j), temp2
            end if
         end do

         ! Print explanatory variable data
         if (isodr) then
            write (lunrpt, 2010)
            if (chkjac .and. &
               ((msgb1 >= 1) .or. &
               (msgd1 >= 1))) then
               write (lunrpt, 2110)
            elseif (anajac) then
               write (lunrpt, 2120)
            else
               write (lunrpt, 2130)
            end if
         else
            write (lunrpt, 2020)
            write (lunrpt, 2140)
         end if
         if (isodr) then
            do j = 1, m
               tempc0 = '1,'
               do i = 1, n, n - 1

                  if (ifixx(1, 1) < 0) then
                     tempc1 = '   NO'
                  else
                     if (ldifx == 1) then
                        if (ifixx(1, j) == 0) then
                           tempc1 = '  YES'
                        else
                           tempc1 = '   NO'
                        end if
                     else
                        if (ifixx(i, j) == 0) then
                           tempc1 = '  YES'
                        else
                           tempc1 = '   NO'
                        end if
                     end if
                  end if

                  if (tt(1, 1) < zero) then
                     temp1 = abs(tt(1, 1))
                  else
                     if (ldtt == 1) then
                        temp1 = tt(1, j)
                     else
                        temp1 = tt(i, j)
                     end if
                  end if

                  if (wd(1, 1, 1) < zero) then
                     temp2 = abs(wd(1, 1, 1))
                  else
                     if (ldwd == 1) then
                        if (ld2wd == 1) then
                           temp2 = wd(1, 1, j)
                        else
                           temp2 = wd(1, j, j)
                        end if
                     else
                        if (ld2wd == 1) then
                           temp2 = wd(i, 1, j)
                        else
                           temp2 = wd(i, j, j)
                        end if
                     end if
                  end if

                  if (anajac) then
                     if (chkjac .and. &
                        (((msgb1 >= 1) .or. &
                           (msgd1 >= 1)) .and. &
                        (i == 1))) then
                        itemp = -1
                        do l = 1, nq
                           itemp = max(itemp, msgd(l, j))
                        end do
                        if (itemp <= -1) then
                           tempc2 = '    UNCHECKED'
                        elseif (itemp == 0) then
                           tempc2 = '     VERIFIED'
                        elseif (itemp >= 1) then
                           tempc2 = ' QUESTIONABLE'
                        end if
                     else
                        tempc2 = '             '
                     end if
                     if (m <= 9) then
                        write (lunrpt, 5110) &
                           tempc0, j, x(i, j), &
                           delta(i, j), tempc1, temp1, temp2, tempc2
                     else
                        write (lunrpt, 5120) &
                           tempc0, j, x(i, j), &
                           delta(i, j), tempc1, temp1, temp2, tempc2
                     end if
                  else
                     tempc2 = '             '
                     if (cdjac) then
                        temp3 = dhstep(1, neta, i, j, stpd, ldstpd)
                     else
                        temp3 = dhstep(0, neta, i, j, stpd, ldstpd)
                     end if
                     if (m <= 9) then
                        write (lunrpt, 5210) &
                           tempc0, j, x(i, j), &
                           delta(i, j), tempc1, temp1, temp2, temp3
                     else
                        write (lunrpt, 5220) &
                           tempc0, j, x(i, j), &
                           delta(i, j), tempc1, temp1, temp2, temp3
                     end if
                  end if

                  tempc0 = 'N,'

               end do
               if (j < m) write (lunrpt, 6000)
            end do
         else

            do j = 1, m
               tempc0 = '1,'
               do i = 1, n, n - 1
                  if (m <= 9) then
                     write (lunrpt, 5110) &
                        tempc0, j, x(i, j)
                  else
                     write (lunrpt, 5120) &
                        tempc0, j, x(i, j)
                  end if
                  tempc0 = 'N,'
               end do
               if (j < m) write (lunrpt, 6000)
            end do
         end if

         ! Print response variable data and observation error weights
         if (.not. implct) then
            write (lunrpt, 3000)
            write (lunrpt, 3100)
            do l = 1, nq
               tempc0 = '1,'
               do i = 1, n, n - 1
                  if (we(1, 1, 1) < zero) then
                     temp1 = abs(we(1, 1, 1))
                  elseif (ldwe == 1) then
                     if (ld2we == 1) then
                        temp1 = we(1, 1, l)
                     else
                        temp1 = we(1, l, l)
                     end if
                  else
                     if (ld2we == 1) then
                        temp1 = we(i, 1, l)
                     else
                        temp1 = we(i, l, l)
                     end if
                  end if
                  if (nq <= 9) then
                     write (lunrpt, 5110) &
                        tempc0, l, y(i, l), temp1
                  else
                     write (lunrpt, 5120) &
                        tempc0, l, y(i, l), temp1
                  end if
                  tempc0 = 'N,'
               end do
               if (l < nq) write (lunrpt, 6000)
            end do
         end if
      end if

      ! Format statements

   1000 format &
         (/' --- Problem Size:'/ &
         '            N = ', I5, &
         '          (number with nonzero weight = ', I5, ')'/ &
         '           NQ = ', I5/ &
         '            M = ', I5/ &
         '           NP = ', I5, &
         '          (number unfixed = ', I5, ')')
   1100 format &
         (/' --- Control Values:'/ &
         '          JOB = ', I5.5/ &
         '              = ABCDE, where')
   1110 format &
         ('                       A=', I1, ' ==> fit is a restart.')
   1111 format &
         ('                       A=', I1, ' ==> fit is not a restart.')
   1120 format &
         ('                       B=', I1, ' ==> deltas are initialized', &
         ' to zero.')
   1121 format &
         ('                       B=', I1, ' ==> deltas are initialized', &
         ' by user.')
   1122 format &
         ('                       B=', I1, ' ==> deltas are fixed at', &
         ' zero since E=', I1, '.')
   1130 format &
         ('                       C=', I1, ' ==> covariance matrix will', &
         ' be computed using')
   1131 format &
         ('                               derivatives re-', &
         'evaluated at the solution.')
   1132 format &
         ('                               derivatives from the', &
         ' last iteration.')
   1133 format &
         ('                       C=', I1, ' ==> covariance matrix will', &
         ' not be computed.')
   1140 format &
         ('                       D=', I1, ' ==> derivatives are', &
         ' supplied by user.')
   1141 format &
         ('                               derivatives were checked.'/ &
         '                               results appear questionable.')
   1142 format &
         ('                               derivatives were checked.'/ &
         '                               results appear correct.')
   1143 format &
         ('                               derivatives were not', &
         ' checked.')
   1144 format &
         ('                       D=', I1, ' ==> derivatives are', &
         ' estimated by central', &
         ' differences.')
   1145 format &
         ('                       D=', I1, ' ==> derivatives are', &
         ' estimated by forward', &
         ' differences.')
   1150 format &
         ('                       E=', I1, ' ==> method is implicit ODR.')
   1151 format &
         ('                       E=', I1, ' ==> method is explicit ODR.')
   1152 format &
         ('                       E=', I1, ' ==> method is explicit OLS.')
   1200 format &
         ('       NDIGIT = ', I5, '          (estimated by ODRPACK)')
   1210 format &
         ('       NDIGIT = ', I5, '          (supplied by user)')
   1300 format &
         ('       TAUFAC = ', 1P, E12.2)
   1400 format &
         (/' --- Stopping Criteria:'/ &
         '        SSTOL = ', 1P, E12.2, &
         '   (sum of squares stopping tolerance)'/ &
         '       PARTOL = ', 1P, E12.2, &
         '   (parameter stopping tolerance)'/ &
         '        MAXIT = ', I5, &
         '          (maximum number of iterations)')
   1500 format &
         (/' --- Initial Sum of Squared Weighted Deltas =', &
         17X, 1P, E17.8)
   1510 format &
         ('         Initial Penalty Function Value     =', 1P, E17.8/ &
         '                 Penalty Term               =', 1P, E17.8/ &
         '                 Penalty Parameter          =', 1P, E10.1)
   1600 format &
         (/' --- Initial Weighted Sum of Squares        =', &
         17X, 1P, E17.8)
   1610 format &
         ('         Sum of Squared Weighted Deltas     =', 1P, E17.8/ &
         '         Sum of Squared Weighted Epsilons   =', 1P, E17.8)
   2010 format &
         (/' --- Explanatory Variable and Delta Weight Summary:')
   2020 format &
         (/' --- Explanatory Variable Summary:')
   2110 format &
         (/'       Index      X(I,J)  DELTA(I,J)    Fixed', &
         '     Scale    Weight    Derivative'/ &
         '                                             ', &
         '                        Assessment'/, &
         '       (I,J)                          (IFIXX)', &
         '    (SCLD)      (WD)              '/)
   2120 format &
         (/'       Index      X(I,J)  DELTA(I,J)    Fixed', &
         '     Scale    Weight              '/ &
         '                                             ', &
         '                                  '/, &
         '       (I,J)                          (IFIXX)', &
         '    (SCLD)      (WD)              '/)
   2130 format &
         (/'       Index      X(I,J)  DELTA(I,J)    Fixed', &
         '     Scale    Weight    Derivative'/ &
         '                                             ', &
         '                         Step Size'/, &
         '       (I,J)                          (IFIXX)', &
         '    (SCLD)      (WD)        (STPD)'/)
   2140 format &
         (/'       Index      X(I,J)'/ &
         '       (I,J)            '/)
   3000 format &
         (/' --- Response Variable and Epsilon Error Weight', &
         ' Summary:')
   3100 format &
         (/'       Index      Y(I,L)      Weight'/ &
         '       (I,L)                    (WE)'/)
   4000 format &
         (/' --- Function Parameter Summary:')
   4110 format &
         (/'       Index   BETA(K)    Fixed     Scale   LOWER(K)', &
         '   UPPER(K)    Derivative'/ &
         '                                                    ', &
         '               Assessment'/, &
         '         (K)            (IFIXB)    (SCLB)           ', &
         '                         '/)
   4120 format &
         (/'       Index   BETA(K)    Fixed     Scale   LOWER(K)', &
         '   UPPER(K)              '/ &
         '                                                    ', &
         '                         '/, &
         '         (K)            (IFIXB)    (SCLB)           ', &
         '                         '/)
   4200 format &
         (/'       Index   BETA(K)    Fixed     Scale   LOWER(K)', &
         '   UPPER(K)    Derivative'/ &
         '                                                    ', &
         '                Step Size'/, &
         '         (K)            (IFIXB)    (SCLB)           ', &
         '                   (STPB)'/)
   4310 format &
         (7X, I5, 1P, E10.2, 4X, A5, E10.2, E11.2E3, E11.2E3, 1X, A13)
   4320 format &
         (7X, I5, 1P, E10.2, 4X, A5, E10.2, E11.2E3, E11.2E3, 1X, E13.5)
   5110 format &
         (9X, A2, I1, 1P, 2E12.3, 4X, A5, 2E10.2, 1X, A13)
   5120 format &
         (8X, A2, I2, 1P, 2E12.3, 4X, A5, 2E10.2, 1X, A13)
   5210 format &
         (9X, A2, I1, 1P, 2E12.3, 4X, A5, 2E10.2, 1X, E13.5)
   5220 format &
         (8X, A2, I2, 1P, 2E12.3, 4X, A5, 2E10.2, 1X, E13.5)
   6000 format &
         (' ')
   end subroutine dodpc1

   impure subroutine dodpc2 &
      (ipr, lunrpt, fstitr, implct, prtpen, &
      pnlty, &
      niter, nfev, wss, actred, prered, alpha, tau, pnorm, np, beta)
   !! Generate iteration reports.
   ! Routines Called  (NONE)
   ! Date Written   860529   (YYMMDD)
   ! Revision Date  920304   (YYMMDD)

      use odrpack_kinds, only: zero

      integer, intent(in) :: ipr
         !! The value indicating the report to be printed.
      integer, intent(in) :: lunrpt
         !! The logical unit number used for computation reports.
      logical, intent(in) :: fstitr
         !! The variable designating whether this is the first iteration (`fstitr = .true.`) or
         !! not (`fstitr = .false.`).
      logical, intent(in) :: implct
         !! The variable designating whether the solution is by implicit ODR (`implct = .true.`)
         !! or explicit ODR (`implct = .false.`).
      logical, intent(in) :: prtpen
         !! The variable designating whether the penalty parameter is to be printed in the
         !! iteration report (`prtpen = .true.`) or not (`prtpen = .false.`).
      real(wp), intent(in) :: pnlty
         !! The penalty parameter for an implicit model.
      integer, intent(in) :: niter
         !! The number of iterations.
      integer, intent(in) :: nfev
         !! The number of function evaluations.
      real(wp), intent(in) :: wss
         !! The sum-of-squares of the weighted `epsilon`s and `deltas`.
      real(wp), intent(in) :: actred
         !! The actual relative reduction in the sum-of-squares.
      real(wp), intent(in) :: prered
         !! The predicted relative reduction in the sum-of-squares.
      real(wp), intent(in) :: alpha
         !! The Levenberg-Marquardt parameter.
      real(wp), intent(in) :: tau
         !! The trust region diameter.
      real(wp), intent(in) :: pnorm
         !! The norm of the scaled estimated parameters.
      integer, intent(in) :: np
         !! The number of function parameters.
      real(wp), intent(in) :: beta(np)
         !! The function parameters.

      ! Local scalars
      real(wp) :: ratio
      integer :: j, k, l
      character(len=3) :: gn

      ! Variable Definitions (alphabetically)
      !  ACTRED:  The actual relative reduction in the sum-of-squares.
      !  ALPHA:   The Levenberg-Marquardt parameter.
      !  BETA:    The function parameters.
      !  FSTITR:  The variable designating whether this is the first iteration (FSTITR=.TRUE.)
      !           or not (FSTITR=.FALSE.).
      !  GN:      The CHARACTER*3 variable indicating whether a Gauss-Newton step was taken.
      !  IMPLCT:  The variable designating whether the solution is by implicit ODR (IMPLCT=TRUE)
      !           or explicit ODR (IMPLCT=FALSE).
      !  IPR:     The value indicating the report to be printed.
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  L:       An indexing variable.
      !  LUNRPT:  The logical unit number used for computation reports.
      !  NFEV:    The number of function evaluations.
      !  NITER:   The number of iterations.
      !  NP:      The number of function parameters.
      !  PNLTY:   The penalty parameter for an implicit model.
      !  PNORM:   The norm of the scaled estimated parameters.
      !  PRERED:  The predicted relative reduction in the sum-of-squares.
      !  PRTPEN:  The variable designating whether the penalty parameter is to be printed in
      !           the iteration report (PRTPEN=TRUE) or not (PRTPEN=FALSE).
      !  RATIO:   The ratio of TAU to PNORM.
      !  TAU:     The trust region diameter.
      !  WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS.

      if (fstitr) then
         if (ipr == 1) then
            if (implct) then
               write (lunrpt, 1121)
            else
               write (lunrpt, 1122)
            end if
         else
            if (implct) then
               write (lunrpt, 1131)
            else
               write (lunrpt, 1132)
            end if
         end if
      end if
      if (prtpen) then
         write (lunrpt, 1133) pnlty
      end if

      if (alpha == zero) then
         gn = 'YES'
      else
         gn = ' NO'
      end if
      if (pnorm /= zero) then
         ratio = tau/pnorm
      else
         ratio = zero
      end if
      if (ipr == 1) then
         write (lunrpt, 1141) niter, nfev, wss, actred, prered, &
            ratio, gn
      else
         j = 1
         k = min(3, np)
         if (j == k) then
            write (lunrpt, 1141) niter, nfev, wss, actred, prered, &
               ratio, gn, j, beta(j)
         else
            write (lunrpt, 1142) niter, nfev, wss, actred, prered, &
               ratio, gn, j, k, (beta(l), l=j, k)
         end if
         if (np > 3) then
            do j = 4, np, 3
               k = min(j + 2, np)
               if (j == k) then
                  write (lunrpt, 1151) j, beta(j)
               else
                  write (lunrpt, 1152) j, k, (beta(l), l=j, k)
               end if
            end do
         end if
      end if

      ! Format statements

   1121 format &
         (// &
         '         Cum.      Penalty    Act. Rel.   Pred. Rel.'/ &
         '  It.  No. FN     Function   Sum-of-Sqs   Sum-of-Sqs', &
         '              G-N'/ &
         ' Num.   Evals        Value    Reduction    Reduction', &
         '  TAU/PNORM  Step'/ &
         ' ----  ------  -----------  -----------  -----------', &
         '  ---------  ----')
   1122 format &
         (// &
         '         Cum.                 Act. Rel.   Pred. Rel.'/ &
         '  It.  No. FN     Weighted   Sum-of-Sqs   Sum-of-Sqs', &
         '              G-N'/ &
         ' Num.   Evals   Sum-of-Sqs    Reduction    Reduction', &
         '  TAU/PNORM  Step'/ &
         ' ----  ------  -----------  -----------  -----------', &
         '  ---------  ----'/)
   1131 format &
         (// &
         '         Cum.      Penalty    Act. Rel.   Pred. Rel.'/ &
         '  It.  No. FN     Function   Sum-of-Sqs   Sum-of-Sqs', &
         '              G-N      BETA -------------->'/ &
         ' Num.   Evals        Value    Reduction    Reduction', &
         '  TAU/PNORM  Step     Index           Value'/ &
         ' ----  ------  -----------  -----------  -----------', &
         '  ---------  ----     -----           -----')
   1132 format &
         (// &
         '         Cum.                 Act. Rel.   Pred. Rel.'/ &
         '  It.  No. FN     Weighted   Sum-of-Sqs   Sum-of-Sqs', &
         '              G-N      BETA -------------->'/ &
         ' Num.   Evals   Sum-of-Sqs    Reduction    Reduction', &
         '  TAU/PNORM  Step     Index           Value'/ &
         ' ----  ------  -----------  -----------  -----------', &
         '  ---------  ----     -----           -----'/)
   1133 format &
         (/' Penalty Parameter Value = ', 1P, E10.1)
   1141 format &
         (1X, I4, I8, 1X, 1P, E12.5, 2E13.4, E11.3, 3X, A3, 7X, I3, 3E16.8)
   1142 format &
         (1X, I4, I8, 1X, 1P, E12.5, 2E13.4, E11.3, 3X, A3, 1X, I3, ' To', I3, 3E16.8)
   1151 format &
         (76X, I3, 1P, E16.8)
   1152 format &
         (70X, I3, ' To', I3, 1P, 3E16.8)

   end subroutine dodpc2

   impure subroutine dodpc3 &
      (ipr, lunrpt, &
      isodr, implct, didvcv, dovcv, redoj, anajac, &
      n, m, np, nq, npp, &
      info, niter, nfev, njev, irank, rcond, istop, &
      wss, wssdel, wsseps, pnlty, rvar, idf, &
      beta, sdbeta, ifixb2, f, delta, &
      lower, upper)
   !! Generate final summary report.
   ! Routines Called  DPPT
   ! Date Written   860529   (YYMMDD)
   ! REvision Date  920619   (YYMMDD)

      use odrpack_core, only: dppt

      integer, intent(in) :: ipr
         !! The variable indicating what is to be printed.
      integer, intent(in) :: lunrpt
         !! The logical unit number used for computation reports.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      logical, intent(in) :: implct
         !! The variable designating whether the solution is by implicit ODR (`implct = .true.`)
         !! or explicit ODR (`implct = .false.`).
      logical, intent(in) :: didvcv
         !! The variable designating whether the covariance matrix was computed (`didvcv = .true.`)
         !! or not (`didvcv = .false.`).
      logical, intent(in) :: dovcv
         !! The variable designating whether the covariance matrix was to be computed
         !! (`dovcv = .true.`) or not (`dovcv = .false.`).
      logical, intent(in) :: redoj
         !! The variable designating whether the Jacobian matrix is to be recomputed for the
         !! computation of the covariance matrix (`redoj = .true.`) or not (`redoj = .false.`).
      logical, intent(in) :: anajac
         !! The variable designating whether the Jacobians are computed by finite differences
         !! (`anajac = .false.`) or not (`anajac = .true.`).
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(in) :: npp
         !! The number of function parameters being estimated.
      integer, intent(in) :: info
         !! The variable designating why the computations were stopped.
      integer, intent(in) :: niter
         !! The number of iterations.
      integer, intent(in) :: nfev
         !! The number of function evaluations.
      integer, intent(in) :: njev
         !! The number of Jacobian evaluations.
      integer, intent(in) :: irank
         !! The rank deficiency of the Jacobian with respect to `beta`.
      real(wp), intent(in) :: rcond
         !! The approximate reciprocal condition of `tfjacb`.
      integer, intent(in) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.
      real(wp), intent(in) :: wss
         !! The sum-of-squares of the weighted `epsilon`s and `delta`s.
      real(wp), intent(in) :: wssdel
         !! The sum-of-squares of the weighted `delta`s.
      real(wp), intent(in) :: wsseps
         !! The sum-of-squares of the weighted `epsilon`s.
      real(wp), intent(in) :: pnlty
         !! The penalty parameter for an implicit model.
      real(wp), intent(in) :: rvar
         !! The residual variance.
      integer, intent(in) :: idf
         !! The degrees of freedom of the fit, equal to the number of observations with nonzero
         !! weighted derivatives minus the number of parameters being estimated.
      real(wp), intent(in) :: beta(np)
         !! The function parameters.
      real(wp), intent(in) :: sdbeta(np)
         !! The standard errors of the estimated parameters.
      integer, intent(in) :: ifixb2(np)
         !! The values designating whether the elements of `beta` were estimated, fixed, or
         !! dropped because they caused rank deficiency.
      real(wp), intent(in) :: f(n, nq)
         !! The estimated values of `epsilon`.
      real(wp), intent(in) :: delta(n, m)
         !! The estimated errors in the explanatory variables.
      real(wp), intent(in) :: lower(np)
         !! Lower bound on `beta`.
      real(wp), intent(in) :: upper(np)
         !! Upper bound on `beta`.

      ! Local scalars
      real(wp) :: tval
      integer :: d1, d2, d3, d4, d5, i, j, k, l, nplm1
      character(len=90) :: fmt1

      ! Variable Definitions (alphabetically)
      !  ANAJAC:  The variable designating whether the JACOBIANS are computed by finite
      !           differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
      !  BETA:    The function parameters.
      !  D1:      The first digit of INFO.
      !  D2:      The second digit of INFO.
      !  D3:      The third digit of INFO.
      !  D4:      The fourth digit of INFO.
      !  D5:      The fifth digit of INFO.
      !  DELTA:   The estimated errors in the explanatory variables.
      !  DIDVCV:  The variable designating whether the covariance matrix was computed
      !           (DIDVCV=TRUE) or not (DIDVCV=FALSE).
      !  DOVCV:   The variable designating whether the covariance matrix was to be computed
      !           (DOVCV=TRUE) or not (DOVCV=FALSE).
      !  F:       The estimated values of EPSILON.
      !  FMT1:    A CHARACTER*90 variable used for formats.
      !  I:       An indexing variable.
      !  IDF:     The degrees of freedom of the fit, equal to the number of observations with
      !           nonzero weighted derivatives minus the number of parameters being estimated.
      !  IFIXB2:  The values designating whether the elements of BETA were estimated, fixed,
      !           or dropped because they caused rank deficiency, corresponding to values of
      !           IFIXB2 equaling 1, 0, and -1, respectively.  If IFIXB2 is -2, then no attempt
      !           was made to estimate the parameters because MAXIT = 0.
      !  IMPLCT:  The variable designating whether the solution is by implicit ODR (IMPLCT=TRUE)
      !           or explicit ODR (IMPLCT=FALSE).
      !  INFO:    The variable designating why the computations were stopped.
      !  IPR:     The variable indicating what is to be printed.
      !  IRANK:   The rank deficiency of the Jacobian wrt BETA.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or by
      !           OLS (ISODR=FALSE).
      !  ISTOP:   The variable designating whether there are problems computing the function at
      !           the current BETA and DELTA.
      !  J:       An indexing variable.
      !  K:       An indexing variable.
      !  L:       An indexing variable.
      !  LOWER:   Lower bound on BETA.
      !  LUNRPT:  The logical unit number used for computation reports.
      !  M:       The number of columns of data in the explanatory variable.
      !  N:       The number of observations.
      !  NFEV:    The number of function evaluations.
      !  NITER:   The number of iterations.
      !  NJEV:    The number of Jacobian evaluations.
      !  NP:      The number of function parameters.
      !  NPLM1:   The number of items to be printed per line, minus one.
      !  NPP:     The number of function parameters being estimated.
      !  NQ:      The number of responses per observation.
      !  PNLTY:   The penalty parameter for an implicit model.
      !  RCOND:   The approximate reciprocal condition of TFJACB.
      !  REDOJ:   The variable designating whether the Jacobian matrix is to be recomputed for
      !           the computation of the covariance matrix (REDOJ=TRUE) or not (REDOJ=FALSE).
      !  RVAR:    The residual variance.
      !  SDBETA:  The standard errors of the estimated parameters.
      !  TVAL:    The value of the 97.5 percent point function for the T distribution.
      !  UPPER:   Upper bound on BETA.
      !  WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS.
      !  WSSDEL:  The sum-of-squares of the weighted DELTAS.
      !  WSSEPS:  The sum-of-squares of the weighted EPSILONS.

      d1 = info/10000
      d2 = mod(info, 10000)/1000
      d3 = mod(info, 1000)/100
      d4 = mod(info, 100)/10
      d5 = mod(info, 10)

      ! Print stopping conditions
      write (lunrpt, 1000)
      if (info <= 9) then
         if (info == 1) then
            write (lunrpt, 1011) info
         elseif (info == 2) then
            write (lunrpt, 1012) info
         elseif (info == 3) then
            write (lunrpt, 1013) info
         elseif (info == 4) then
            write (lunrpt, 1014) info
         elseif (info <= 9) then
            write (lunrpt, 1015) info
         end if
      elseif (info <= 9999) then

         ! Print warning diagnostics
         write (lunrpt, 1020) info
         if (d2 == 1) write (lunrpt, 1021)
         if (d3 == 1) write (lunrpt, 1022)
         if (d4 == 1) write (lunrpt, 1023)
         if (d4 == 2) write (lunrpt, 1024)
         if (d5 == 1) then
            write (lunrpt, 1031)
         elseif (d5 == 2) then
            write (lunrpt, 1032)
         elseif (d5 == 3) then
            write (lunrpt, 1033)
         elseif (d5 == 4) then
            write (lunrpt, 1034)
         elseif (d5 <= 9) then
            write (lunrpt, 1035) d5
         end if
      else

         ! Print error messages
         write (lunrpt, 1040) info
         if (d1 == 5) then
            write (lunrpt, 1042)
            if (d2 /= 0) write (lunrpt, 1043) d2
            if (d3 == 3) then
               write (lunrpt, 1044) d3
            elseif (d3 /= 0) then
               write (lunrpt, 1045) d3
            end if
         elseif (d1 == 6) then
            write (lunrpt, 1050)
         else
            write (lunrpt, 1060) d1
         end if
      end if

      ! Print misc. stopping info
      write (lunrpt, 1300) niter
      write (lunrpt, 1310) nfev
      if (anajac) write (lunrpt, 1320) njev
      write (lunrpt, 1330) irank
      write (lunrpt, 1340) rcond
      write (lunrpt, 1350) istop

      ! Print final sum of squares
      if (implct) then
         write (lunrpt, 2000) wssdel
         if (isodr) then
            write (lunrpt, 2010) wss, wsseps, pnlty
         end if
      else
         write (lunrpt, 2100) wss
         if (isodr) then
            write (lunrpt, 2110) wssdel, wsseps
         end if
      end if
      if (didvcv) then
         write (lunrpt, 2200) sqrt(rvar), idf
      end if

      nplm1 = 3

      ! Print estimated BETA's, and, if, full rank, their standard errors
      write (lunrpt, 3000)
      if (didvcv) then
         write (lunrpt, 7300)
         tval = dppt(0.975E0_wp, idf)
         do j = 1, np
            if (ifixb2(j) >= 1) then
               write (lunrpt, 8400) j, beta(j), &
                  lower(j), upper(j), &
                  sdbeta(j), &
                  beta(j) - tval*sdbeta(j), &
                  beta(j) + tval*sdbeta(j)
            elseif (ifixb2(j) == 0) then
               write (lunrpt, 8600) j, beta(j), lower(j), upper(j)
            else
               write (lunrpt, 8700) j, beta(j), lower(j), upper(j)
            end if
         end do
         if (.not. redoj) write (lunrpt, 7310)
      else
         if (dovcv) then
            if (d1 <= 5) then
               write (lunrpt, 7410)
            else
               write (lunrpt, 7420)
            end if
         end if

         if ((irank == 0 .and. npp == np) .or. niter == 0) then
            if (np == 1) then
               write (lunrpt, 7100)
            else
               write (lunrpt, 7200)
            end if
            do j = 1, np, nplm1 + 1
               k = min(j + nplm1, np)
               if (k == j) then
                  write (lunrpt, 8100) j, beta(j)
               else
                  write (lunrpt, 8200) j, k, (beta(l), l=j, k)
               end if
            end do
            if (niter >= 1) then
               write (lunrpt, 8800)
            else
               write (lunrpt, 8900)
            end if
         else
            write (lunrpt, 7500)
            do j = 1, np
               if (ifixb2(j) >= 1) then
                  write (lunrpt, 8500) j, beta(j), lower(j), upper(j)
               elseif (ifixb2(j) == 0) then
                  write (lunrpt, 8600) j, beta(j), lower(j), upper(j)
               else
                  write (lunrpt, 8700) j, beta(j), lower(j), upper(j)
               end if
            end do
         end if
      end if

      if (ipr == 1) return

      ! Print EPSILON's and DELTA's together in a column if the number of
      ! columns of data in EPSILON and DELTA is less than or equal to three.
      if (implct .and. &
         (m <= 4)) then
         write (lunrpt, 4100)
         write (fmt1, 9110) m
         write (lunrpt, fmt1) (j, j=1, m)
         do i = 1, n
            write (lunrpt, 4130) i, (delta(i, j), j=1, m)
         end do

      elseif (isodr .and. &
            (nq + m <= 4)) then
         write (lunrpt, 4110)
         write (fmt1, 9120) nq, m
         write (lunrpt, fmt1) (l, l=1, nq), (j, j=1, m)
         do i = 1, n
            write (lunrpt, 4130) i, (f(i, l), l=1, nq), (delta(i, j), j=1, m)
         end do

      elseif (.not. isodr .and. &
            ((nq >= 2) .and. &
               (nq <= 4))) then
         write (lunrpt, 4120)
         write (fmt1, 9130) nq
         write (lunrpt, fmt1) (l, l=1, nq)
         do i = 1, n
            write (lunrpt, 4130) i, (f(i, l), l=1, nq)
         end do
      else

         ! Print EPSILON's and DELTA's separately
         if (.not. implct) then

            ! Print EPSILON'S
            do j = 1, nq
               write (lunrpt, 4200) j
               if (n == 1) then
                  write (lunrpt, 7100)
               else
                  write (lunrpt, 7200)
               end if
               do i = 1, n, nplm1 + 1
                  k = min(i + nplm1, n)
                  if (i == k) then
                     write (lunrpt, 8100) i, f(i, j)
                  else
                     write (lunrpt, 8200) i, k, (f(l, j), l=i, k)
                  end if
               end do
            end do
         end if

         ! Print DELTA'S
         if (isodr) then
            do j = 1, m
               write (lunrpt, 4300) j
               if (n == 1) then
                  write (lunrpt, 7100)
               else
                  write (lunrpt, 7200)
               end if
               do i = 1, n, nplm1 + 1
                  k = min(i + nplm1, n)
                  if (i == k) then
                     write (lunrpt, 8100) i, delta(i, j)
                  else
                     write (lunrpt, 8200) i, k, (delta(l, j), l=i, k)
                  end if
               end do
            end do
         end if
      end if

      ! Format statements

   1000 format &
         (/' --- Stopping Conditions:')
   1011 format &
         ('         INFO = ', I5, ' ==> sum of squares convergence.')
   1012 format &
         ('         INFO = ', I5, ' ==> parameter convergence.')
   1013 format &
         ('         INFO = ', I5, ' ==> sum of squares convergence and', &
         ' parameter convergence.')
   1014 format &
         ('         INFO = ', I5, ' ==> iteration limit reached.')
   1015 format &
         ('         INFO = ', I5, ' ==> unexpected value,', &
         ' probably indicating'/ &
         '                           incorrectly specified', &
         ' user input.')
   1020 format &
         ('         INFO = ', I5.4/ &
         '              =  ABCD, where a nonzero value for digit A,', &
         ' B, or C indicates why'/ &
         '                       the results might be questionable,', &
         ' and digit D indicates'/ &
         '                       the actual stopping condition.')
   1021 format &
         ('                       A=1 ==> derivatives are', &
         ' questionable.')
   1022 format &
         ('                       B=1 ==> user set ISTOP to', &
         ' nonzero value during last'/ &
         '                               call to subroutine FCN.')
   1023 format &
         ('                       C=1 ==> derivatives are not', &
         ' full rank at the solution.')
   1024 format &
         ('                       C=2 ==> derivatives are zero', &
         ' rank at the solution.')
   1031 format &
         ('                       D=1 ==> sum of squares convergence.')
   1032 format &
         ('                       D=2 ==> parameter convergence.')
   1033 format &
         ('                       D=3 ==> sum of squares convergence', &
         ' and parameter convergence.')
   1034 format &
         ('                       D=4 ==> iteration limit reached.')
   1035 format &
         ('                       D=', I1, ' ==> unexpected value,', &
         ' probably indicating'/ &
         '                               incorrectly specified', &
         ' user input.')
   1040 format &
         ('         INFO = ', I5.5/ &
         '              = ABCDE, where a nonzero value for a given', &
         ' digit indicates an'/ &
         '                       abnormal stopping condition.')
   1042 format &
         ('                       A=5 ==> user stopped computations', &
         ' in subroutine FCN.')
   1043 format &
         ('                       B=', I1, ' ==> computations were', &
         ' stopped during the'/ &
         '                                    function evaluation.')
   1044 format &
         ('                       C=', I1, ' ==> computations were', &
         ' stopped because'/ &
         '                                    derivatives with', &
         ' respect to delta were'/ &
         '                                    computed by', &
         ' subroutine FCN when'/ &
         '                                    fit is OLS.')
   1045 format &
         ('                       C=', I1, ' ==> computations were', &
         ' stopped during the'/ &
         '                                    jacobian evaluation.')
   1050 format &
         ('                       A=6 ==> numerical instabilities', &
         ' have been detected,'/ &
         '                               possibly indicating', &
         ' a discontinuity in the'/ &
         '                               derivatives or a poor', &
         ' poor choice of problem'/ &
         '                               scale or weights.')
   1060 format &
         ('                       A=', I1, ' ==> unexpected value,', &
         ' probably indicating'/ &
         '                               incorrectly specified', &
         ' user input.')
   1300 format &
         ('        NITER = ', I5, &
         '          (number of iterations)')
   1310 format &
         ('         NFEV = ', I5, &
         '          (number of function evaluations)')
   1320 format &
         ('         NJEV = ', I5, &
         '          (number of jacobian evaluations)')
   1330 format &
         ('        IRANK = ', I5, &
         '          (rank deficiency)')
   1340 format &
         ('        RCOND = ', 1P, E12.2, &
         '   (inverse condition number)')
   !1341 FORMAT
   !       +  ('                      ==> POSSIBLY FEWER THAN 2 SIGNIFICANT',
   !       +                        ' DIGITS IN RESULTS;'/
   !       +   '                          SEE ODRPACK95 REFERENCE',
   !       +                        ' GUIDE, SECTION 4.C.')
   1350 format &
         ('        ISTOP = ', I5, &
         '          (returned by user from', &
         ' subroutine FCN)')
   2000 format &
         (/' --- Final Sum of Squared Weighted Deltas = ', &
         17X, 1P, E17.8)
   2010 format &
         ('         Final Penalty Function Value     = ', 1P, E17.8/ &
         '               Penalty Term               = ', 1P, E17.8/ &
         '               Penalty Parameter          = ', 1P, E10.1)
   2100 format &
         (/' --- Final Weighted Sums of Squares       = ', 17X, 1P, E17.8)
   2110 format &
         ('         Sum of Squared Weighted Deltas   = ', 1P, E17.8/ &
         '         Sum of Squared Weighted Epsilons = ', 1P, E17.8)
   2200 format &
         (/' --- Residual Standard Deviation          = ', &
         17X, 1P, E17.8/ &
         '         Degrees of Freedom               =', I5)
   3000 format &
         (/' --- Estimated BETA(J), J = 1, ..., NP:')
   4100 format &
         (/' --- Estimated DELTA(I,*), I = 1, ..., N:')
   4110 format &
         (/' --- Estimated EPSILON(I) and DELTA(I,*), I = 1, ..., N:')
   4120 format &
         (/' --- Estimated EPSILON(I), I = 1, ..., N:')
   4130 format(5X, I5, 1P, 5E16.8)
   4200 format &
         (/' --- Estimated EPSILON(I,', I3, '), I = 1, ..., N:')
   4300 format &
         (/' --- Estimated DELTA(I,', I3, '), I = 1, ..., N:')
   7100 format &
         (/'           Index           Value'/)
   7200 format &
         (/'           Index           Value -------------->'/)
   7300 format &
         (/'                     BETA      LOWER     UPPER      S.D. ', &
         ' ___ 95% Confidence ___'/ &
         '                                                    BETA ', &
         '        Interval'/)
   7310 format &
         (/'     N.B. standard errors and confidence intervals are', &
         ' computed using'/ &
         '          derivatives calculated at the beginning', &
         ' of the last iteration,'/ &
         '          and not using derivatives re-evaluated at the', &
         ' final solution.')
   7410 format &
         (/'     N.B. the standard errors of the estimated betas were', &
         ' not computed because'/ &
         '          the derivatives were not available.  Either MAXIT', &
         ' is 0 and the third'/ &
         '          digit of JOB is greater than 1, or the most', &
         ' recently tried values of'/ &
         '          BETA and/or X+DELTA were identified as', &
         ' unacceptable by user supplied'/ &
         '          subroutine FCN.')
   7420 format &
         (/'     N.B. the standard errors of the estimated betas were', &
         ' not computed.'/ &
         '          (see info above.)')
   7500 format &
         (/'                     BETA         Status')
   8100 format &
         (11X, I5, 1P, E16.8)
   8200 format &
         (3X, I5, ' to', I5, 1P, 7E16.8)
   8400 format &
         (3X, I5, 1X, 1P, E16.8, 1X, E10.2, E10.2, E10.2, 1X, E10.2, 1X, 'to', E10.2)
   8500 format &
         (3X, I5, 1X, 1P, E16.8, 1X, E10.2, E10.2, 4X, 'Estimated')
   8600 format &
         (3X, I5, 1X, 1P, E16.8, 1X, E10.2, E10.2, 4X, '    Fixed')
   8700 format &
         (3X, I5, 1X, 1P, E16.8, 1X, E10.2, E10.2, 4X, '  Dropped')
   8800 format &
         (/'     N.B. no parameters were fixed by the user or', &
         ' dropped at the last'/ &
         '          iteration because they caused the model to be', &
         ' rank deficient.')
   8900 format &
         (/'     N.B. no change was made to the user supplied parameter', &
         ' values because'/ &
         '          MAXIT=0.')
   9110 format &
         ('(/''         I'',', &
         I2, '(''      DELTA(I,'',I1,'')'')/)')
   9120 format &
         ('(/''         I'',', &
         I2, '(''    EPSILON(I,'',I1,'')''),', &
         I2, '(''      DELTA(I,'',I1,'')'')/)')
   9130 format &
         ('(/''         I'',', &
         I2, '(''    EPSILON(I,'',I1,'')'')/)')

   end subroutine dodpc3

   impure subroutine dodpcr &
      (ipr, lunrpt, &
      head, prtpen, fstitr, didvcv, iflag, &
      n, m, np, nq, npp, nnzw, &
      msgb, msgd, beta, y, ldy, x, ldx, delta, &
      we, ldwe, ld2we, wd, ldwd, ld2wd, &
      ifixb, ifixx, ldifx, &
      lower, upper, &
      ssf, tt, ldtt, stpb, stpd, ldstpd, &
      job, neta, taufac, sstol, partol, maxit, &
      wss, rvar, idf, sdbeta, &
      niter, nfev, njev, actred, prered, &
      tau, pnorm, alpha, f, rcond, irank, info, istop)
   !! Generate computation reports.
   ! Routines Called  DFLAGS, DODPC1, DODPC2, DODPC3, DODPHD
   ! Date Written   860529   (YYMMDD)
   ! Revision Date  920619   (YYMMDD)

      use odrpack_core, only: dflags

      integer, intent(in) :: ipr
         !! The variable indicating what is to be printed.
      integer, intent(in) :: lunrpt
         !! The logical unit number used for computation reports.
      logical, intent(inout) :: head
         !! The variable designating whether the heading is to be printed (`head = .true.`)
         !! or not (`head = .false.`).
      logical, intent(in) :: prtpen
         !! The variable designating whether the penalty parameter is to be printed in the
         !! iteration report (`prtpen = .true.`) or not (`prtpen = .false.`).
      logical, intent(in) :: fstitr
         !! The variable designating whether this is the first iteration (`fstitr = .true.`)
         !! or not (`fstitr = .false.`).
      logical, intent(in) :: didvcv
         !! The variable designating whether the covariance matrix was computed
         !! (`didvcv = .true.`) or not (`didvcv = .false.`).
      integer, intent(in) :: iflag
         !! The variable designating what is to be printed.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(in) :: npp
         !! The number of function parameters being estimated.
      integer, intent(in) :: nnzw
         !! The number of nonzero weighted observations.
      integer, intent(in) :: msgb(nq*np + 1)
         !! The error checking results for the Jacobian with respect to `beta`.
      integer, intent(in) :: msgd(nq*m + 1)
         !! The error checking results for the Jacobian with respect to `delta`.
      real(wp), intent(in) :: beta(np)
         !! The function parameters.
      real(wp), intent(in) :: y(ldy, nq)
         !! The response variable. Unused when the model is implicit.
      integer, intent(in) :: ldy
         !! The leading dimension of array `y`.
      real(wp), intent(in) :: x(ldx, m)
         !! The explanatory variable.
      integer, intent(in) :: ldx
         !! The leading dimension of array `x`.
      real(wp), intent(in) :: delta(n, m)
         !! The estimated errors in the explanatory variables.
      real(wp), intent(in) :: we(ldwe, ld2we, nq)
         !! The `epsilon` weights.
      integer, intent(in) :: ldwe
         !! The leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! The second dimension of array `we`.
      real(wp), intent(in) :: wd(ldwd, ld2wd, m)
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
      real(wp), intent(in) :: lower(np)
         !! The lower bounds for `beta`.
      real(wp), intent(in) :: upper(np)
         !! The upper bounds for `beta`.
      real(wp), intent(in) :: ssf(np)
         !! The scaling values for `beta`.
      real(wp), intent(in) :: tt(ldtt, m)
         !! The scaling values for `delta`.
      integer, intent(in) :: ldtt
         !! The leading dimension of array `tt`.
      real(wp), intent(in) :: stpb(np)
         !! The relative step used for computing finite difference derivatives with respect
         !! to `beta`.
      real(wp), intent(in) :: stpd(ldstpd, m)
         !! The relative step used for computing finite difference derivatives with respect
         !! to `delta`.
      integer, intent(in) :: ldstpd
         !! The leading dimension of array `stpd`.
      integer, intent(in) :: job
         !! The variable controlling problem initialization and computational method.
      integer, intent(in) :: neta
         !! The number of accurate digits in the function results.
      real(wp), intent(in) :: taufac
         !! The factor used to compute the initial trust region diameter.
      real(wp), intent(in) :: sstol
         !! The sum-of-squares convergence stopping tolerance.
      real(wp), intent(in) :: partol
         !! The parameter convergence stopping tolerance.
      integer, intent(in) :: maxit
         !! The maximum number of iterations allowed.
      real(wp), intent(in) :: wss(3)
         !! The sum-of-squares of the weighted `epsilon`s and `delta`s, the sum-of-squares of 
         !! the weighted `delta`s, and the sum-of-squares of the weighted `epsilon`s.
      real(wp), intent(in) :: rvar
         !! The residual variance.
      integer, intent(in) :: idf
         !! The degrees of freedom of the fit, equal to the number of observations with nonzero
         !! weighted derivatives minus the number of parameters being estimated.
      real(wp), intent(in) :: sdbeta(np)
         !! The standard errors of the estimated parameters.
      integer, intent(in) :: niter
         !! The number of iterations.
      integer, intent(in) :: nfev
         !! The number of function evaluations.
      integer, intent(in) :: njev
         !! The number of Jacobian evaluations.
      real(wp), intent(in) :: actred
         !! The actual relative reduction in the sum-of-squares.
      real(wp), intent(in) :: prered
         !! The predicted relative reduction in the sum-of-squares.
      real(wp), intent(in) :: tau
         !! The trust region diameter.
      real(wp), intent(in) :: pnorm
         !! The norm of the scaled estimated parameters.
      real(wp), intent(in) :: alpha
         !! The Levenberg-Marquardt parameter.
      real(wp), intent(in) :: f(n, nq)
         !! The estimated values of `epsilon`.
      real(wp), intent(in) :: rcond
         !! The approximate reciprocal condition of `tfjacb`.
      integer, intent(in) :: irank
         !! The rank deficiency of the Jacobian with respect to `beta`.
      integer, intent(in) :: info
         !! The variable designating why the computations were stopped.
      integer, intent(in) :: istop
         !! The variable designating whether there are problems computing the function at the
         !! current `beta` and `delta`.

      ! Local scalars
      real(wp) :: pnlty
      logical :: anajac, cdjac, chkjac, dovcv, implct, initd, isodr, redoj, restrt
      character(len=3) :: typ

      ! Variable Definitions (alphabetically)
      !  ACTRED:  The actual relative reduction in the sum-of-squares.
      !  ALPHA:   The Levenberg-Marquardt parameter.
      !  ANAJAC:  The variable designating whether the Jacobians are computed by finite
      !           differences (ANAJAC=FALSE) or not (ANAJAC=TRUE).
      !  BETA:    The function parameters.
      !  CDJAC:   The variable designating whether the jacobians are computed by central
      !           differences (CDJAC=TRUE) or by forward differences (CDJAC=FALSE).
      !  CHKJAC:  The variable designating whether the user supplied Jacobians are to be checked
      !           (CHKJAC=TRUE) or not (CHKJAC=FALSE).
      !  DELTA:   The estimated errors in the explanatory variables.
      !  DIDVCV:  The variable designating whether the covariance matrix was computed
      !           (DIDVCV=TRUE) or not (DIDVCV=FALSE).
      !  DOVCV:   The variable designating whether the covariance matrix was to be computed
      !           (DOVCV=TRUE) or not (DOVCV=FALSE).
      !  F:       The (weighted) estimated values of EPSILON.
      !  FSTITR:  The variable designating whether this is the first iteration (FSTITR=TRUE) or
      !           not (FSTITR=FALSE).
      !  HEAD:    The variable designating whether the heading is to be printed (HEAD=TRUE) or
      !           not (HEAD=FALSE).
      !  IDF:     The degrees of freedom of the fit, equal to the number of observations with
      !           nonzero weighted derivatives minus the number of parameters being estimated.
      !  IFIXB:   The values designating whether the elements of BETA are fixed at their input
      !           values or not.
      !  IFIXX:   The values designating whether the elements of X are fixed at their input
      !           values or not.
      !  IFLAG:   The variable designating what is to be printed.
      !  IMPLCT:  The variable designating whether the solution is by implicit ODR (IMPLCT=TRUE)
      !           or explicit ODR (IMPLCT=FALSE).
      !  INFO:    The variable designating why the computations were stopped.
      !  INITD:   The variable designating whether DELTA is initialized to zero (INITD=TRUE) or
      !           to the values in the first N by M elements of array WORK (INITD=FALSE).
      !  IPR:     The value indicating the report to be printed.
      !  IRANK:   The rank deficiency of the Jacobian wrt BETA.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or by
      !           OLS (ISODR=FALSE).
      !  ISTOP:   The variable designating whether there are problems computing the function at
      !           the current BETA and DELTA.
      !  JOB:     The variable controling problem initialization and computational method.
      !  LDIFX:   The leading dimension of array IFIXX.
      !  LDSTPD:  The leading dimension of array STPD.
      !  LDTT:    The leading dimension of array TT.
      !  LDWD:    The leading dimension of array WD.
      !  LDWE:    The leading dimension of array WE.
      !  LDX:     The leading dimension of array X.
      !  LDY:     The leading dimension of array Y.
      !  LD2WD:   The second dimension of array WD.
      !  LD2WE:   The second dimension of array WE.
      !  LOWER:   Lower bound on BETA.
      !  LUNRPT:  The logical unit number for computation reports.
      !  M:       The number of columns of data in the explanatory variable.
      !  MAXIT:   The maximum number of iterations allowed.
      !  MSGB:    The error checking results for the Jacobian wrt BETA.
      !  MSGD:    The error checking results for the Jacobian wrt DELTA.
      !  N:       The number of observations.
      !  NETA:    The number of accurate digits in the function results.
      !  NFEV:    The number of function evaluations.
      !  NITER:   The number of iterations.
      !  NJEV:    The number of Jacobian evaluations.
      !  NNZW:    The number of nonzero weighted observations.
      !  NP:      The number of function parameters.
      !  NQ:      The number of responses per observation.
      !  NPP:     The number of function parameters being estimated.
      !  PARTOL:  The parameter convergence stopping tolerance.
      !  PNLTY:   The penalty parameter for an implicit model.
      !  PNORM:   The norm of the scaled estimated parameters.
      !  PRERED:  The predicted relative reduction in the sum-of-squares.
      !  PRTPEN:  The variable designating whether the penalty parameter is  to be printed in 
      !           the iteration report (PRTPEN=TRUE) or not (PRTPEN=FALSE).
      !  RCOND:   The approximate reciprocal condition number of TFJACB.
      !  REDOJ:   The variable designating whether the Jacobian matrix is to be recomputed for 
      !           the computation of the covariance matrix (REDOJ=TRUE) or not (REDOJ=FALSE).
      !  RESTRT:  The variable designating whether the call is a restart (RESTRT=TRUE) or not
      !           (RESTRT=FALSE).
      !  RVAR:    The residual variance.
      !  SDBETA:  The standard deviations of the estimated BETA'S.
      !  SSF:     The scaling values for BETA.
      !  SSTOL:   The sum-of-squares convergence stopping tolerance.
      !  STPB:    The relative step for computing finite difference derivatives with respect
      !           to BETA.
      !  STPD:    The relative step for computing finite difference derivatives with respect
      !           to DELTA.
      !  TAU:     The trust region diameter.
      !  TAUFAC:  The factor used to compute the initial trust region diameter.
      !  TT:      The scaling values for DELTA.
      !  TYP:     The CHARACTER*3 string "ODR" or "OLS".
      !  UPPER:   Upper bound on BETA.
      !  WE:      The EPSILON weights.
      !  WD:      The DELTA weights.
      !  WSS:     The sum-of-squares of the weighted EPSILONS and DELTAS, the sum-of-squares of
      !           the weighted DELTAS, and the sum-of-squares of the weighted EPSILONS.
      !  X:       The explanatory variable.
      !  Y:       The dependent variable.  Unused when the model is implicit.

      call dflags(job, restrt, initd, dovcv, redoj, anajac, cdjac, chkjac, isodr, implct)
      pnlty = abs(we(1, 1, 1))

      if (head) then
         call dodphd(head, lunrpt)
      end if
      if (isodr) then
         typ = 'ODR'
      else
         typ = 'OLS'
      end if

      ! Print initial summary
      if (iflag == 1) then
         write (lunrpt, 1200) typ
         call dodpc1 &
            (ipr, lunrpt, &
            anajac, cdjac, chkjac, initd, restrt, isodr, implct, dovcv, redoj, &
            msgb(1), msgb(2), msgd(1), msgd(2), n, m, np, nq, npp, nnzw, x, &
            ldx, ifixx, ldifx, delta, wd, ldwd, ld2wd, tt, ldtt, stpd, ldstpd, &
            y, ldy, we, ldwe, ld2we, pnlty, beta, ifixb, ssf, stpb, lower, &
            upper, job, neta, taufac, sstol, partol, maxit, wss(1), wss(2), &
            wss(3))

         ! Print iteration reports
      elseif (iflag == 2) then

         if (fstitr) then
            write (lunrpt, 1300) typ
         end if
         call dodpc2 &
            (ipr, lunrpt, fstitr, implct, prtpen, &
            pnlty, &
            niter, nfev, wss(1), actred, prered, alpha, tau, pnorm, np, beta)

         ! Print final summary
      elseif (iflag == 3) then

         write (lunrpt, 1400) typ
         call dodpc3 &
            (ipr, lunrpt, &
            isodr, implct, didvcv, dovcv, redoj, anajac, &
            n, m, np, nq, npp, &
            info, niter, nfev, njev, irank, rcond, istop, &
            wss(1), wss(2), wss(3), pnlty, rvar, idf, &
            beta, sdbeta, ifixb, f, delta, lower, upper)
      end if

      ! Format statements

   1200 format &
         (/' *** Initial summary for fit by method of ', A3, ' ***')
   1300 format &
         (/' *** Iteration reports for fit by method of ', A3, ' ***')
   1400 format &
         (/' *** Final summary for fit by method of ', A3, ' ***')

   end subroutine dodpcr

   impure subroutine dodpe1 &
      (lunerr, info, d1, d2, d3, d4, d5, &
      n, m, nq, &
      ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
      lwkmn, liwkmn)
   !! Print error reports.
   ! Routines Called  (NONE)
   ! Date Written   860529   (YYMMDD)
   ! Revision Date  920619   (YYMMDD)

      integer, intent(in) :: lunerr
         !! The logical unit number used for error messages.
      integer, intent(inout) :: info
         !! The variable designating why the computations were stopped.
      integer, intent(in) :: d1
         !! The 1st digit (from the left) of `info`.
      integer, intent(in) :: d2
         !! The 2nd digit (from the left) of `info`.
      integer, intent(in) :: d3
         !! The 3rd digit (from the left) of `info`.
      integer, intent(in) :: d4
         !! The 4th digit (from the left) of `info`.
      integer, intent(in) :: d5
         !! The 5th digit (from the left) of `info`.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(in) :: ldscld
         !! The leading dimension of array `scld`.
      integer, intent(in) :: ldstpd
         !! The leading dimension of array `stpd`.
      integer, intent(in) :: ldwe
         !! The leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! The second dimension of array `we`.
      integer, intent(in) :: ldwd
         !! The leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! The second dimension of array `wd`.
      integer, intent(in) :: lwkmn
         !! The minimum acceptable length of array `work`.
      integer, intent(in) :: liwkmn
         !! The minimum acceptable length of array `iwork`.

      ! Variable Definitions (alphabetically)
      !  D1:      The 1st digit (from the left) of INFO.
      !  D2:      The 2nd digit (from the left) of INFO.
      !  D3:      The 3rd digit (from the left) of INFO.
      !  D4:      The 4th digit (from the left) of INFO.
      !  D5:      The 5th digit (from the left) of INFO.
      !  INFO:    The variable designating why the computations were stopped.
      !  LDSCLD:  The leading dimension of array SCLD.
      !  LDSTPD:  The leading dimension of array STPD.
      !  LDWD:    The leading dimension of array WD.
      !  LDWE:    The leading dimension of array WE.
      !  LIWKMN:  The minimum acceptable length of array IWORK.
      !  LWKMN:   The minimum acceptable length of array WORK.
      !  LD2WD:   The second dimension of array WD.
      !  LD2WE:   The second dimension of array WE.
      !  LUNERR:  The logical unit number used for error messages.
      !  M:       The number of columns of data in the explanatory variable.
      !  N:       The number of observations.
      !  NQ:      The number of responses per observation.

      ! Print appropriate messages for errors in problem specification parameters

      if (d1 == 1) then
         if (d2 /= 0) then
            write (lunerr, 1100)
         end if
         if (d3 /= 0) then
            write (lunerr, 1200)
         end if
         if (d4 /= 0) then
            write (lunerr, 1300)
         end if
         if (d5 /= 0) then
            write (lunerr, 1400)
         end if

         ! Print appropriate messages for errors in dimension specification parameters

      elseif (d1 == 2) then

         if (d2 /= 0) then
            if (d2 == 1 .or. d2 == 3) then
               write (lunerr, 2110)
            end if
            if (d2 == 2 .or. d2 == 3) then
               write (lunerr, 2120)
            end if
         end if

         if (d3 /= 0) then
            if (d3 == 1 .or. d3 == 3 .or. d3 == 5 .or. d3 == 7) &
               then
               write (lunerr, 2210)
            end if
            if (d3 == 2 .or. d3 == 3 .or. d3 == 6 .or. d3 == 7) &
               then
               write (lunerr, 2220)
            end if
            if (d3 == 4 .or. d3 == 5 .or. d3 == 6 .or. d3 == 7) &
               then
               write (lunerr, 2230)
            end if
         end if

         if (d4 /= 0) then
            if (d4 == 1 .or. d4 == 3) then
               write (lunerr, 2310)
            end if
            if (d4 == 2 .or. d4 == 3) then
               write (lunerr, 2320)
            end if
         end if

         if (d5 /= 0) then
            if (d5 == 1 .or. d5 == 3) then
               write (lunerr, 2410) lwkmn
            end if
            if (d5 == 2 .or. d5 == 3) then
               write (lunerr, 2420) liwkmn
            end if
         end if

      elseif (d1 == 3) then

         ! Print appropriate messages for errors in scale values

         if (d3 /= 0) then
            if (d3 == 2 .or. d3 == 3) then
               if (ldscld >= n) then
                  write (lunerr, 3110)
               else
                  write (lunerr, 3120)
               end if
            end if
            if (d3 == 1 .or. d3 == 3) then
               write (lunerr, 3130)
            end if
         end if

         ! Print appropriate messages for errors in derivative step values

         if (d2 /= 0) then
            if (d2 == 2 .or. d2 == 3) then
               if (ldstpd >= n) then
                  write (lunerr, 3210)
               else
                  write (lunerr, 3220)
               end if
            end if
            if (d2 == 1 .or. d2 == 3) then
               write (lunerr, 3230)
            end if
         end if

         ! Print appropriate messages for errors in observational error weights

         if (d4 /= 0) then
            if (d4 == 1) then
               if (ldwe >= n) then
                  if (ld2we >= nq) then
                     write (lunerr, 3310)
                  else
                     write (lunerr, 3320)
                  end if
               else
                  if (ld2we >= nq) then
                     write (lunerr, 3410)
                  else
                     write (lunerr, 3420)
                  end if
               end if
            end if
            if (d4 == 2) then
               write (lunerr, 3500)
            end if
         end if

         ! Print appropriate messages for errors in DELTA weights

         if (d5 /= 0) then
            if (ldwd >= n) then
               if (ld2wd >= m) then
                  write (lunerr, 4310)
               else
                  write (lunerr, 4320)
               end if
            else
               if (ld2wd >= m) then
                  write (lunerr, 4410)
               else
                  write (lunerr, 4420)
               end if
            end if
         end if

      elseif (d1 == 7) then

         ! Print the appropriate messages for errors in JOB

         if (d2 /= 0) then
            write (lunerr, 5000)
         end if

         if (d3 /= 0) then
            write (lunerr, 5100)
         end if

         if (d4 /= 0) then
            write (lunerr, 5200)
         end if

      elseif (d1 == 8) then

         ! Print the appropriate messages for errors in array allocation

         if (d2 /= 0) then
            write (lunerr, 7200)
         end if

         if (d3 /= 0) then
            write (lunerr, 7300)
         end if

         if (d4 /= 0) then
            write (lunerr, 7400)
         end if

      elseif (d1 == 9) then

         ! Print the appropriate messages for errors in bounds

         if (d2 /= 0) then
            write (lunerr, 6000)
         end if

         if (d3 /= 0) then
            write (lunerr, 6100)
         end if

         if (d4 == 1) then
            write (lunerr, 6210)
         end if

         if (d4 == 2) then
            write (lunerr, 6220)
         end if

      end if

      ! Print error messages for array sizes incorrect

      if (info/100000 == 1) then
         info = info - 100000
         if (info >= 32768) then
            info = info - 32768
            write (lunerr, 8015)
         end if
         if (info >= 16384) then
            info = info - 16384
            write (lunerr, 8014)
         end if
         if (info >= 8192) then
            info = info - 8192
            write (lunerr, 8013)
         end if
         if (info >= 4096) then
            info = info - 4096
            write (lunerr, 8012)
         end if
         if (info >= 2048) then
            info = info - 2048
            write (lunerr, 8011)
         end if
         if (info >= 1024) then
            info = info - 1024
            write (lunerr, 8010)
         end if
         if (info >= 512) then
            info = info - 512
            write (lunerr, 8009)
         end if
         if (info >= 256) then
            info = info - 256
            write (lunerr, 8008)
         end if
         if (info >= 128) then
            info = info - 128
            write (lunerr, 8007)
         end if
         if (info >= 64) then
            info = info - 64
            write (lunerr, 8006)
         end if
         if (info >= 32) then
            info = info - 32
            write (lunerr, 8005)
         end if
         if (info >= 16) then
            info = info - 16
            write (lunerr, 8004)
         end if
         if (info >= 8) then
            info = info - 8
            write (lunerr, 8003)
         end if
         if (info >= 4) then
            info = info - 4
            write (lunerr, 8002)
         end if
         if (info >= 2) then
            info = info - 2
            write (lunerr, 8001)
         end if
         if (info >= 1) then
            info = info - 1
            write (lunerr, 8000)
         end if
      end if

      ! Format statements

   1100 format &
         (/' ERROR :  N is less than one.')
   1200 format &
         (/' ERROR :  M is less than one.')
   1300 format &
         (/' ERROR :  NP is less than one'/ &
         '          or NP is greater than N.')
   1400 format &
         (/' ERROR :  NQ is less than one.')
   2110 format &
         (/' ERROR :  LDX is less than N.')
   2120 format &
         (/' ERROR :  LDY is less than N.')
   2210 format &
         (/' ERROR :  LDIFX is less than N'/ &
         '          and LDIFX is not equal to one.')
   2220 format &
         (/' ERROR :  LDSCLD is less than N'/ &
         '          and LDSCLD is not equal to one.')
   2230 format &
         (/' ERROR :  LDSTPD is less than N'/ &
         '          and LDSTPD is not equal to one.')
   2310 format &
         (/' ERROR :  LDWE is less than N'/ &
         '          and LDWE is not equal to one or'/ &
         '          or'/ &
         '          LD2WE is less than NQ'/ &
         '          and LD2WE is not equal to one.')
   2320 format &
         (/' ERROR :  LDWD is less than N'/ &
         '          and LDWD is not equal to one.')
   2410 format &
         (/' ERROR :  LWORK is less than ', I7, ','/ &
         '          the smallest acceptable dimension of array WORK.')
   2420 format &
         (/' ERROR :  LIWORK is less than ', I7, ','/ &
         '          the smallest acceptable dimension of array', &
         ' IWORK.')
   3110 format &
         (/' ERROR :  SCLD(I,J) is less than or equal to zero'/ &
         '          for some I = 1, ..., N and J = 1, ..., M.'// &
         '          when SCLD(1,1) is greater than zero'/ &
         '          and LDSCLD is greater than or equal to N then'/ &
         '          each of the N by M elements of'/ &
         '          SCLD must be greater than zero.')
   3120 format &
         (/' ERROR :  SCLD(1,J) is less than or equal to zero'/ &
         '          for some J = 1, ..., M.'// &
         '          when SCLD(1,1) is greater than zero'/ &
         '          and LDSCLD is equal to one then'/ &
         '          each of the 1 by M elements of'/ &
         '          SCLD must be greater than zero.')
   3130 format &
         (/' ERROR :  SCLB(K) is less than or equal to zero'/ &
         '          for some K = 1, ..., NP.'// &
         '          all NP elements of', &
         '          SCLB must be greater than zero.')
   3210 format &
         (/' ERROR :  STPD(I,J) is less than or equal to zero'/ &
         '          for some I = 1, ..., N and J = 1, ..., M.'// &
         '          when STPD(1,1) is greater than zero'/ &
         '          and LDSTPD is greater than or equal to N then'/ &
         '          each of the N by M elements of'/ &
         '          STPD must be greater than zero.')
   3220 format &
         (/' ERROR :  STPD(1,J) is less than or equal to zero'/ &
         '          for some J = 1, ..., M.'// &
         '          when STPD(1,1) is greater than zero'/ &
         '          and LDSTPD is equal to one then'/ &
         '          each of the 1 by M elements of'/ &
         '          STPD must be greater than zero.')
   3230 format &
         (/' ERROR :  STPB(K) is less than or equal to zero'/ &
         '          for some K = 1, ..., NP.'// &
         '          all NP elements of', &
         ' STPB must be greater than zero.')
   3310 format &
         (/' ERROR :  At least one of the (NQ by NQ) arrays starting'/ &
         '          in WE(I,1,1), I = 1, ..., N, is not positive'/ &
         '          semidefinite.  When WE(1,1,1) is greater than'/ &
         '          or equal to zero, and LDWE is greater than or'/ &
         '          equal to N, and LD2WE is greater than or equal'/ &
         '          to NQ, then each of the (NQ by NQ) arrays in WE'/ &
         '          must be positive semidefinite.')
   3320 format &
         (/' ERROR :  At least one of the (1 by NQ) arrays starting'/ &
         '          in WE(I,1,1), I = 1, ..., N, has a negative'/ &
         '          element.  When WE(1,1,1) is greater than or'/ &
         '          equal to zero, and LDWE is greater than or equal'/ &
         '          to N, and LD2WE is equal to 1, then each of the'/ &
         '          (1 by NQ) arrays in WE must have only non-'/ &
         '          negative elements.')
   3410 format &
         (/' ERROR :  The (NQ by NQ) array starting in WE(1,1,1) is'/ &
         '          not positive semidefinite.  When WE(1,1,1) is'/ &
         '          greater than or equal to zero, and LDWE is equal'/ &
         '          to 1, and LD2WE is greater than or equal to NQ,'/ &
         '          then the (NQ by NQ) array in WE must be positive'/ &
         '          semidefinite.')
   3420 format &
         (/' ERROR :  The (1 by NQ) array starting in WE(1,1,1) has'/ &
         '          a negative element.  When WE(1,1,1) is greater'/ &
         '          than or equal to zero, and LDWE is equal to 1,'/ &
         '          and LD2WE is equal to 1, then the (1 by NQ)'/ &
         '          array in WE must have only nonnegative elements.')
   3500 format &
         (/' ERROR :  The number of nonzero arrays in array WE is'/ &
         '          less than NP.')
   4310 format &
         (/' ERROR :  At least one of the (M by M) arrays starting'/ &
         '          in WD(I,1,1), I = 1, ..., N, is not positive'/ &
         '          definite.  When WD(1,1,1) is greater than zero,'/ &
         '          and LDWD is greater than or equal to N, and'/ &
         '          LD2WD is greater than or equal to M, then each'/ &
         '          of the (M by M) arrays in WD must be positive'/ &
         '          definite.')
   4320 format &
         (/' ERROR :  At least one of the (1 by M) arrays starting'/ &
         '          in WD(I,1,1), I = 1, ..., N, has a nonpositive'/ &
         '          element.  When WD(1,1,1) is greater than zero,'/ &
         '          and LDWD is greater than or equal to N, and'/ &
         '          LD2WD is equal to 1, then each of the (1 by M)'/ &
         '          arrays in WD must have only positive elements.')
   4410 format &
         (/' ERROR :  The (M by M) array starting in WD(1,1,1) is'/ &
         '          not positive definite.  When WD(1,1,1) is'/ &
         '          greater than zero, and LDWD is equal to 1, and'/ &
         '          LD2WD is greater than or equal to M, then the'/ &
         '          (M by M) array in WD must be positive definite.')
   4420 format &
         (/' ERROR :  The (1 by M) array starting in WD(1,1,1) has a'/ &
         '          nonpositive element.  When WD(1,1,1) is greater'/ &
         '          than zero, and LDWD is equal to 1, and LD2WD is'/ &
         '          equal to 1, then the (1 by M) array in WD must'/ &
         '          have only positive elements.')
   5000 format &
         (/' ERROR :  JOB requires the optional argument DELTA and'/ &
         '          DELTA is not present or not associated.')
   5100 format &
         (/' ERROR :  JOB requires the optional argument WORK and'/ &
         '          WORK is not present or not associated.')
   5200 format &
         (/' ERROR :  JOB requires the optional argument IWORK and'/ &
         '          IWORK is not present or not associated.')
   6000 format &
         (/' ERROR :  LOWER(K)>UPPER(K) for some K.  Adjust the'/ &
         '          the bounds so that LOWER(K)<=UPPER(K) holds'/ &
         '          for all K.')
   6100 format &
         (/' ERROR :  BETA(K)>UPPER(K) or BETA(K)<LOWER(K) '/ &
         '          for some K.  Adjust the bounds or BETA so '/ &
         '          that LOWER(K)<=BETA(K)<=UPPER(K) holds'/ &
         '          for all K.')
   6210 format &
         (/' ERROR :  UPPER(K)-LOWER(K) < 400*BETA(K)*EPSMAC  '/ &
         '          for some K and EPSMAC having the largest '/ &
         '          value such that 1+EPSMAC/=1.  This '/ &
         '          constraint on UPPER and LOWER is necessary'/ &
         '          for the calculation of NDIGIT.  Increase the'/ &
         '          range of the bounds or specify NDIGIT '/ &
         '          explicitly.')
   6220 format &
         (/' ERROR :  UPPER(K)-LOWER(K) < ABS(STEP) for some'/ &
         '          K where step is the step size for numeric'/ &
         '          derivatives.  Increase the bounds or supply'/ &
         '          an analytic jacobian.')
   7200 format &
         (/' ERROR :  DELTA could not be allocated. ')
   7300 format &
         (/' ERROR :  WORK could not be allocated. ')
   7400 format &
         (/' ERROR :  IWORK could not be allocated. ')
   8000 format &
         (/' ERROR :  BETA has incorrect size. ')
   8001 format &
         (/' ERROR :  Y has incorrect size. ')
   8002 format &
         (/' ERROR :  X has incorrect size. ')
   8003 format &
         (/' ERROR :  DELTA has incorrect size. ')
   8004 format &
         (/' ERROR :  WE has incorrect size. ')
   8005 format &
         (/' ERROR :  WD has incorrect size. ')
   8006 format &
         (/' ERROR :  IFIXB has incorrect size. ')
   8007 format &
         (/' ERROR :  IFIXX has incorrect size. ')
   8008 format &
         (/' ERROR :  STPB has incorrect size. ')
   8009 format &
         (/' ERROR :  STPD has incorrect size. ')
   8010 format &
         (/' ERROR :  SCLB has incorrect size. ')
   8011 format &
         (/' ERROR :  SCLD has incorrect size. ')
   8012 format &
         (/' ERROR :  WORK has incorrect size. ')
   8013 format &
         (/' ERROR :  IWORK has incorrect size. ')
   8014 format &
         (/' ERROR :  UPPER has incorrect size. ')
   8015 format &
         (/' ERROR :  LOWER has incorrect size. ')

   end subroutine dodpe1

   impure subroutine dodpe2 &
      (lunerr, &
      n, m, np, nq, &
      fjacb, fjacd, &
      diff, msgb1, msgb, isodr, msgd1, msgd, &
      xplusd, nrow, neta, ntol)
   !! Generate the derivative checking report.
   ! Routines Called  (NONE)
   ! Date Written   860529   (YYMMDD)
   ! Revision Date  920619   (YYMMDD)

      integer, intent(in) :: lunerr
         !! The logical unit number used for error messages.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      real(wp), intent(in) :: fjacb(n, np, nq)
         !! The Jacobian with respect to `beta`.
      real(wp), intent(in) :: fjacd(n, m, nq)
         !! The Jacobian with respect to `delta`.
      real(wp), intent(in) :: diff(nq, np + m)
         !! The relative differences between the user-supplied and finite difference derivatives
         !! for each derivative checked.
      integer, intent(in) :: msgb1
         !! The error checking results for the Jacobian with respect to `beta`.
      integer, intent(in) :: msgb(nq, np)
         !! The error checking results for the Jacobian with respect to `beta`.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`) or by
         !! OLS (`isodr = .false.`).
      integer, intent(in) :: msgd1
         !! The error checking results for the Jacobian with respect to `delta`.
      integer, intent(in) :: msgd(nq, m)
         !! The error checking results for the Jacobian with respect to `delta`.
      real(wp), intent(in) :: xplusd(n, m)
         !! The values of `x` + `delta`.
      integer, intent(in) :: nrow
         !! The row number of the explanatory variable array at which the derivative is to be
         !! checked.
      integer, intent(in) :: neta
         !! The number of reliable digits in the model.
      integer, intent(in) :: ntol
         !! The number of digits of agreement required between the finite difference and the
         !! user-supplied derivatives.

      ! Local scalars
      integer :: i, j, k, l
      character(len=1) :: flag
      character(len=3) :: typ

      ! Local arrays
      logical :: ftnote(0:9)

      ! Variable Definitions (alphabetically)
      !  DIFF:    The relative differences between the user-supplied and finite difference
      !           derivatives for each derivative checked.
      !  FJACB:   The Jacobian with respect to BETA.
      !  FJACD:   The Jacobian with respect to DELTA.
      !  FLAG:    The character string indicating highly questionable results.
      !  FTNOTE:  The array controling footnotes.
      !  I:       An index variable.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=.TRUE.) or
      !           by OLS (ISODR=.FALSE.).
      !  J:       An index variable.
      !  K:       An index variable.
      !  L:       An index variable.
      !  LUNERR:  The logical unit number used for error messages.
      !  M:       The number of columns of data in the explanatory variable.
      !  MSGB:    The error checking results for the Jacobian wrt BETA.
      !  MSGB1:   The error checking results for the Jacobian wrt BETA.
      !  MSGD:    The error checking results for the Jacobian wrt DELTA.
      !  MSGD1:   The error checking results for the Jacobian wrt DELTA.
      !  N:       The number of observations.
      !  NETA:    The number of reliable digits in the model.
      !  NP:      The number of function parameters.
      !  NQ:      The number of responses per observation.
      !  NROW:    The row number of the explanatory variable array at which the derivative is to
      !           be checked.
      !  NTOL:    The number of digits of agreement required between the finite difference and
      !           the user supplied derivatives.
      !  TYP:     The character string indicating solution type, ODR or OLS.
      !  XPLUSD:  The values of X + DELTA.

      ! Set up for footnotes

      ftnote = .false.

      do l = 1, nq
         if (msgb1 >= 1) then
            do i = 1, np
               if (msgb(l, i) >= 1) then
                  ftnote(0) = .true.
                  ftnote(msgb(l, i)) = .true.
               end if
            end do
         end if

         if (msgd1 >= 1) then
            do i = 1, m
               if (msgd(l, i) >= 1) then
                  ftnote(0) = .true.
                  ftnote(msgd(l, i)) = .true.
               end if
            end do
         end if
      end do

      ! Print report

      if (isodr) then
         typ = 'ODR'
      else
         typ = 'OLS'
      end if
      write (lunerr, 1000) typ

      do l = 1, nq

         write (lunerr, 2100) l, nrow
         write (lunerr, 2200)

         do i = 1, np
            k = msgb(l, i)
            if (k == 7) then
               flag = '*'
            else
               flag = ' '
            end if
            if (k <= -1) then
               write (lunerr, 3100) i
            elseif (k == 0) then
               write (lunerr, 3200) i, fjacb(nrow, i, l), diff(l, i), flag
            elseif (k == 8) then
               write (lunerr, 3400) i, fjacb(nrow, i, l), flag, k
            elseif (k == 9) then
               write (lunerr, 3500) i, flag, k
            elseif (k >= 1) then
               write (lunerr, 3300) i, fjacb(nrow, i, l), diff(l, i), flag, &
                  k
            end if
         end do
         if (isodr) then
            do i = 1, m
               k = msgd(l, i)
               if (k == 7) then
                  flag = '*'
               else
                  flag = ' '
               end if
               if (k <= -1) then
                  write (lunerr, 4100) nrow, i
               elseif (k == 0) then
                  write (lunerr, 4200) nrow, i, fjacd(nrow, i, l), diff(l, np + i), flag
               elseif (k >= 1) then
                  write (lunerr, 4300) nrow, i, fjacd(nrow, i, l), diff(l, np + i), flag, k
               end if
            end do
         end if
      end do

      ! Print footnotes

      if (ftnote(0)) then

         write (lunerr, 5000)
         if (ftnote(1)) write (lunerr, 5100)
         if (ftnote(2)) write (lunerr, 5200)
         if (ftnote(3)) write (lunerr, 5300)
         if (ftnote(4)) write (lunerr, 5400)
         if (ftnote(5)) write (lunerr, 5500)
         if (ftnote(6)) write (lunerr, 5600)
         if (ftnote(7)) write (lunerr, 5700)
         if (ftnote(8)) write (lunerr, 5800)
         if (ftnote(9)) write (lunerr, 5900)
      end if

      if (neta < 0) then
         write (lunerr, 6000) - neta
      else
         write (lunerr, 6100) neta
      end if
      write (lunerr, 7000) ntol

      ! Print out row of explanatory variable which was checked.

      write (lunerr, 8100) nrow

      do j = 1, m
         write (lunerr, 8110) nrow, j, xplusd(nrow, j)
      end do

      ! Format statements

   1000 format &
         (//' *** Derivative checking report for fit by method of ', A3, &
         ' ***'/)
   2100 format(/'     For response ', I2, ' of observation ', I5/)
   2200 format('                      ', '         User', &
            '               ', '                '/ &
            '                      ', '     Supplied', &
            '     Relative', '    Derivative '/ &
            '        Derivative WRT', '        Value', &
            '   Difference', '    Assessment '/)
   3100 format('             BETA(', I3, ')', '       ---   ', &
            '       ---   ', '    Unchecked')
   3200 format('             BETA(', I3, ')', 1P, 2E13.2, 3X, A1, &
            'Verified')
   3300 format('             BETA(', I3, ')', 1P, 2E13.2, 3X, A1, &
            'Questionable (see note ', I1, ')')
   3400 format('             BETA(', I3, ')', 1P, 1E13.2, 13X, 3X, A1, &
            'Questionable (see note ', I1, ')')
   3500 format('             BETA(', I3, ')', 1P, 13X, 13X, 3X, A1, &
            'Small bounds (see note ', I1, ')')
   4100 format('          DELTA(', I2, ',', I2, ')', '       ---   ', &
            '       ---   ', '    Unchecked')
   4200 format('          DELTA(', I2, ',', I2, ')', 1P, 2E13.2, 3X, A1, &
            'Verified')
   4300 format('          DELTA(', I2, ',', I2, ')', 1P, 2E13.2, 3X, A1, &
            'Questionable (see note ', I1, ')')
   5000 format &
         (/'     NOTES:')
   5100 format &
         (/'      (1) User supplied and finite difference derivatives', &
         ' agree, but'/ &
         '          results are questionable because both are zero.')
   5200 format &
         (/'      (2) User supplied and finite difference derivatives', &
         ' agree, but'/ &
         '          results are questionable because one is', &
         ' identically zero'/ &
         '          and the other is only approximately zero.')
   5300 format &
         (/'      (3) User supplied and finite difference derivatives', &
         ' disagree, but'/ &
         '          results are questionable because one is', &
         ' identically zero'/ &
         '          and the other is not.')
   5400 format &
         (/'      (4) User supplied and finite difference derivatives', &
         ' disagree, but'/ &
         '          finite difference derivative is questionable', &
         ' because either'/ &
         '          the ratio of relative curvature to relative', &
         ' slope is too high'/ &
         '          or the scale is wrong.')
   5500 format &
         (/'      (5) User supplied and finite difference derivatives', &
         ' disagree, but'/ &
         '          finite difference derivative is questionable', &
         ' because the'/ &
         '          ratio of relative curvature to relative slope is', &
         ' too high.')
   5600 format &
         (/'      (6) User supplied and finite difference derivatives', &
         ' disagree, but'/ &
         '          have at least 2 digits in common.')
   5700 format &
         (/'      (7) User supplied and finite difference derivatives', &
         ' disagree, and'/ &
         '          have fewer than 2 digits in common.  derivative', &
         ' checking must'/ &
         '          be turned off in order to proceed.')
   5800 format &
         (/'      (8) User supplied and finite difference derivatives', &
         ' disagree, and'/ &
         '          bound constraints are too small to calculate', &
         ' further'/ &
         '          information.')
   5900 format &
         (/'      (9) Bound constraints too small to check derivative.')
   6000 format &
         (/'     Number of reliable digits in function results       ', &
         I5/ &
         '        (estimated by ODRPACK)')
   6100 format &
         (/'     Number of reliable digits in function results       ', &
         I5/ &
         '        (supplied by user)')
   7000 format &
         (/'     Number of digits of agreement required between      '/ &
         '     user supplied and finite difference derivative for  '/ &
         '     user supplied derivative to be considered verified  ', &
         I5)
   8100 format &
         (/'     Row number at which derivatives were checked        ', &
         I5// &
         '       -values of the explanatory variables at this row'/)
   8110 format &
         (10X, 'X(', I2, ',', I2, ')', 1X, 1P, 3E16.8)

   end subroutine dodpe2

   impure subroutine dodpe3(lunerr, d2, d3)
   !! Print error reports indicating that computations were stopped in user-supplied 
   !! subroutine `fcn`.
   ! Routines Called  (NONE)
   ! Date Written   860529   (YYMMDD)
   ! Revision Date  920619   (YYMMDD)

      integer :: lunerr
         !! The logical unit number used for error messages.
      integer :: d2
         !! The 2nd digit (from the left) of `info`.
      integer :: d3
         !! The 3rd digit (from the left) of `info`.

      ! Variable Definitions (alphabetically)
      !  D2:      The 2nd digit (from the left) of INFO.
      !  D3:      The 3rd digit (from the left) of INFO.
      !  LUNERR:  The logical unit number used for error messages.

      ! Print appropriate messages to indicate where computations were stopped

      if (d2 == 2) then
         write (lunerr, 1100)
      elseif (d2 == 3) then
         write (lunerr, 1200)
      elseif (d2 == 4) then
         write (lunerr, 1300)
      end if
      if (d3 == 2) then
         write (lunerr, 1400)
      end if

      ! Format statements

   1100 format &
         (//' Variable ISTOP has been returned with a nonzero value  '/ &
         ' from user supplied subroutine FCN when invoked using the'/ &
         ' initial estimates of BETA and DELTA supplied by the     '/ &
         ' user.  The initial estimates must be adjusted to allow  '/ &
         ' proper evaluation of subroutine FCN before the          '/ &
         ' regression procedure can continue.')
   1200 format &
         (//' Variable ISTOP has been returned with a nonzero value  '/ &
         ' from user supplied subroutine FCN.  This occurred during'/ &
         ' the computation of the number of reliable digits in the '/ &
         ' predicted values (F) returned from subroutine FCN, indi-'/ &
         ' cating that changes in the initial estimates of BETA(K),'/ &
         ' K=1,NP, as small as 2*BETA(K)*SQRT(MACHINE PRECISION),  '/ &
         ' where MACHINE PRECISION is defined as the smallest value'/ &
         ' E such that 1+E>1 on the computer being used, prevent   '/ &
         ' subroutine FCN from being properly evaluated.  The      '/ &
         ' initial estimates must be adjusted to allow proper      '/ &
         ' evaluation of subroutine FCN during these computations  '/ &
         ' before the regression procedure can continue.')
   1300 format &
         (//' Variable ISTOP has been returned with a nonzero value  '/ &
         ' from user supplied subroutine FCN.  This occurred during'/ &
         ' the derivative checking procedure, indicating that      '/ &
         ' changes in the initial estimates of BETA(K), K=1,NP, as '/ &
         ' small as MAX[BETA(K),1/SCLB(K)]*10**(-NETA/2), and/or   '/ &
         ' of DELTA(I,J), I=1,N and J=1,M, as small as             '/ &
         ' MAX[DELTA(I,J),1/SCLD(I,J)]*10**(-NETA/2), where NETA   '/ &
         ' is defined to be the number of reliable digits in       '/ &
         ' predicted values (F) returned from subroutine FCN,      '/ &
         ' prevent subroutine FCN from being properly evaluated.   '/ &
         ' the initial estimates must be adjusted to allow proper  '/ &
         ' evaluation of subroutine FCN during these computations  '/ &
         ' before the regression procedure can continue.')
   1400 format &
         (//' Variable ISTOP has been returned with a nonzero value  '/ &
         ' from user supplied subroutine FCN when invoked for '/ &
         ' derivative evaluations using the initial estimates of '/ &
         ' BETA and DELTA supplied by the user.  The initial '/ &
         ' estimates must be adjusted to allow proper evaluation '/ &
         ' of subroutine FCN before the regression procedure can '/ &
         ' continue.')

   end subroutine dodpe3

   impure subroutine dodper &
      (info, lunerr, &
      n, m, np, nq, &
      ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
      lwkmn, liwkmn, &
      fjacb, fjacd, &
      diff, msgb, isodr, msgd, &
      xplusd, nrow, neta, ntol)
   !! Controlling routine for printing error reports.
   ! Routines Called  DODPE1, DODPE2, DODPE3, DODPHD
   ! Date Written   860529   (YYMMDD)
   ! Revision Date  920619   (YYMMDD)

      integer, intent(inout) :: info
         !! The variable designating why the computations were stopped.
      integer, intent(in) :: lunerr
         !! The logical unit number used for error messages.
      integer, intent(in) :: n
         !! The number of observations.
      integer, intent(in) :: m
         !! The number of columns of data in the explanatory variable.
      integer, intent(in) :: np
         !! The number of function parameters.
      integer, intent(in) :: nq
         !! The number of responses per observation.
      integer, intent(in) :: ldscld
         !! The leading dimension of array `scld`.
      integer, intent(in) :: ldstpd
         !! The leading dimension of array `stpd`.
      integer, intent(in) :: ldwe
         !! The leading dimension of array `we`.
      integer, intent(in) :: ld2we
         !! The second dimension of array `we`.
      integer, intent(in) :: ldwd
         !! The leading dimension of array `wd`.
      integer, intent(in) :: ld2wd
         !! The second dimension of array `wd`.
      integer, intent(in) :: lwkmn
         !! The minimum acceptable length of array `work`.
      integer, intent(in) :: liwkmn
         !! The minimum acceptable length of array `iwork`.
      real(wp), intent(in) :: fjacb(n, np, nq)
         !! The Jacobian with respect to `beta`.
      real(wp), intent(in) :: fjacd(n, m, nq)
         !! The Jacobian with respect to `delta`.
      real(wp), intent(in) :: diff(nq, np + m)
         !! The relative differences between the user-supplied and finite difference 
         !! derivatives for each derivative checked.
      integer, intent(in) :: msgb(nq*np + 1)
         !! The error checking results for the Jacobian with respect to `beta`.
      logical, intent(in) :: isodr
         !! The variable designating whether the solution is by ODR (`isodr = .true.`) or
         !! by OLS (`isodr = .false.`).
      integer, intent(in) :: msgd(nq*m + 1)
         !! The error checking results for the Jacobian with respect to `delta`.
      real(wp), intent(in) :: xplusd(n, m)
         !! The values `x + delta`.
      integer, intent(in) :: nrow
         !! The row number of the explanatory variable array at which the derivative is to be
         !! checked.
      integer, intent(in) :: neta
         !! The number of reliable digits in the model.
      integer, intent(in) :: ntol
         !! The number of digits of agreement required between the finite difference and the
         !! user-supplied derivatives.

      ! Local scalars
      integer :: d1, d2, d3, d4, d5
      logical :: head

      ! Variable Definitions (alphabetically)
      !  D1:      The 1st digit (from the left) of INFO.
      !  D2:      The 2nd digit (from the left) of INFO.
      !  D3:      The 3rd digit (from the left) of INFO.
      !  D4:      The 4th digit (from the left) of INFO.
      !  D5:      The 5th digit (from the left) of INFO.
      !  DIFF:    The relative differences between the user supplied and finite difference
      !           derivatives for each derivative checked.
      !  FJACB:   The Jacobian with respect to BETA.
      !  FJACD:   The Jacobian with respect to DELTA.
      !  HEAD:    The variable designating whether the heading is to be printed (HEAD=.TRUE.)
      !           or not (HEAD=.FALSE.).
      !  INFO:    The variable designating why the computations were stopped.
      !  ISODR:   The variable designating whether the solution is by ODR (ISODR=.TRUE.) or
      !           by OLS (ISODR=.FALSE.).
      !  LDSCLD:  The leading dimension of array SCLD.
      !  LDSTPD:  The leading dimension of array STPD.
      !  LDWD:    The leading dimension of array WD.
      !  LDWE:    The leading dimension of array WE.
      !  LD2WD:   The second dimension of array WD.
      !  LD2WE:   The second dimension of array WE.
      !  LIWKMN:  The minimum acceptable length of array IWORK.
      !  LUNERR:  The logical unit number used for error messages.
      !  LWKMN:   The minimum acceptable length of array WORK.
      !  M:       The number of columns of data in the explanatory variable.
      !  MSGB:    The error checking results for the Jacobian wrt BETA.
      !  MSGD:    The error checking results for the Jacobian wrt DELTA.
      !  N:       The number of observations.
      !  NETA:    The number of reliable digits in the model.
      !  NP:      The number of function parameters.
      !  NQ:      The number of responses per observation.
      !  NROW:    The row number of the explanatory variable array at which the derivative is
      !           to be checked.
      !  NTOL:    The number of digits of agreement required between the finite difference and
      !           the user supplied derivatives.
      !  XPLUSD:  The values X + DELTA.

      ! Print heading
      head = .true.
      call dodphd(head, lunerr)

      ! Extract individual digits from variable INFO
      d1 = mod(info, 100000)/10000
      d2 = mod(info, 10000)/1000
      d3 = mod(info, 1000)/100
      d4 = mod(info, 100)/10
      d5 = mod(info, 10)

      ! Print appropriate error messages for ODRPACK invoked stop
      if ((d1 >= 1 .and. d1 <= 3) .or. (d1 == 7 .or. d1 == 9)) then

         ! Print appropriate messages for errors in
         !         problem specification parameters
         !         dimension specification parameters
         !         number of good digits in X
         !         weights
         call dodpe1(lunerr, info, d1, d2, d3, d4, d5, &
                     n, m, nq, &
                     ldscld, ldstpd, ldwe, ld2we, ldwd, ld2wd, &
                     lwkmn, liwkmn)

      elseif ((d1 == 4) .or. (msgb(1) >= 0)) then

         ! Print appropriate messages for derivative checking
         call dodpe2(lunerr, &
                     n, m, np, nq, &
                     fjacb, fjacd, &
                     diff, msgb(1), msgb(2), isodr, msgd(1), msgd(2), &
                     xplusd, nrow, neta, ntol)

      elseif (d1 == 5) then

         ! Print appropriate error message for user invoked stop from FCN
         call dodpe3(lunerr, d2, d3)

      end if

      ! Print correct form of call statement
      if ((d1 >= 1 .and. d1 <= 3) .or. &
         (d1 == 4 .and. (d2 == 2 .or. d3 == 2)) .or. &
         (d1 == 5)) then
         write (lunerr, 1100)
      end if

      ! Format statements

   1100 format &
         (//' The correct form of the call statement is '// &
         '       call odr'/ &
         '      +     (fcn,'/ &
         '      +     n, m, np, nq,'/ &
         '      +     beta,'/ &
         '      +     y, x,'/ &
         '      +     delta*,'/ &
         '      +     we*, wd*,'/ &
         '      +     ifixb*, ifixx*,'/ &
         '      +     job*, ndigit*, taufac*,'/ &
         '      +     sstol*, partol*, maxit*,'/ &
         '      +     iprint*, lunerr*, lunrpt*,'/ &
         '      +     stpb*, stpd*,'/ &
         '      +     sclb*, scld*,'/ &
         '      +     work*, iwork*,'/ &
         '      +     info*,'/ &
         '      +     lower*, upper*)'/ &
         ' * optional argument')
   
   end subroutine dodper

   impure subroutine dodphd(head, lunit)
   !! Print report heading.
   ! Routines Called  (NONE)
   ! Date Written   860529   (YYMMDD)
   ! Revision Date  920619   (YYMMDD)

      logical, intent(inout) :: head
         !! The variable designating whether the heading is to be printed (`head = .true.`) 
         !! or not (`head = .false.`).
      integer, intent(in) :: lunit
         !! The logical unit number to which the heading is written.

      ! Variable Definitions (alphabetically)
      !  HEAD:    The variable designating whether the heading is to be printed (HEAD=.TRUE.) 
      !           or not (HEAD=.FALSE.).
      !  LUNIT:    The logical unit number to which the heading is written.

      if (head) then
         write (lunit, 1000)
         head = .false.
      end if

      ! Format statements

   1000 format( &
         ' ********************************************************* '/ &
         ' * ODRPACK95 version 1.00 of 12-27-2005 (REAL (wp)) * '/ &
         ' ********************************************************* '/ &
         )

   end subroutine dodphd

end module odrpack_reports