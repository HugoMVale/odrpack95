
pure subroutine dwinf &
   (n, m, np, nq, ldwe, ld2we, isodr, &
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
!! Set storage locations within REAL (KIND=wp) work space
   ! Routines Called  (NONE)
   ! Date Written   860529   (YYMMDD)
   ! Revision Date  920619   (YYMMDD)
   integer, intent(in) :: n
      !! The number of observations.
   integer, intent(in) :: m
      !! The number of columns of data in the explanatory variable.
   integer, intent(in) :: np
      !! The number of function parameters.
   integer, intent(in) :: nq
      !! The number of responses per observation.
   integer, intent(in) :: ldwe
      !! The leading dimension of array `we`.
   integer, intent(in) :: ld2we
      !! The second dimension of array `we`.
   logical, intent(in) :: isodr
      !! The variable designating whether the solution is by ODR (`isodr`=.true.) or by OLS (`isodr`=.false.).
   integer, intent(out) :: deltai
      !! The starting location in array `work` of array `delta`.
   integer, intent(out) :: epsi
      !! The starting location in array `work` of array `eps`.
   integer, intent(out) :: xplusi
      !! The starting location in array `work` of array `xplusd`.
   integer, intent(out) :: fni
      !! The starting location in array `work` of array `fn`.
   integer, intent(out) :: sdi
      !! The starting location in array `work` of array `sd`.
   integer, intent(out) :: vcvi
      !! The starting location in array `work` of array `vcv`.
   integer, intent(out) :: rvari
      !! The location in array `work` of variable `rvar`.
   integer, intent(out) :: wssi
      !! The location in array `work` of variable `wss`.
   integer, intent(out) :: wssdei
      !! The location in array `work` of variable `wssdel`.
   integer, intent(out) :: wssepi
      !! The location in array `work` of variable `wsseps`.
   integer, intent(out) :: rcondi
      !! The location in array `work` of variable `rcondi`.
   integer, intent(out) :: etai
      !! The location in array `work` of variable `eta`.
   integer, intent(out) :: olmavi
      !! The location in array `work` of variable `olmavg`.
   integer, intent(out) :: taui
      !! The location in array `work` of variable `tau`.
   integer, intent(out) :: alphai
      !! The location in array `work` of variable `alpha`.
   integer, intent(out) :: actrsi
      !! The location in array `work` of variable `actrs`.
   integer, intent(out) :: pnormi
      !! The location in array `work` of variable `pnorm`.
   integer, intent(out) :: rnorsi
      !! The location in array `work` of variable `rnorms`.
   integer, intent(out) :: prersi
      !! The location in array `work` of variable `prers`.
   integer, intent(out) :: partli
      !! The location in array `work` of variable `partol`.
   integer, intent(out) :: sstoli
      !! The location in array `work` of variable `sstol`.
   integer, intent(out) :: taufci
      !! The location in array `work` of variable `taufac`.
   integer, intent(out) :: epsmai
      !! The location in array `work` of variable `epsmac`.
   integer, intent(out) :: beta0i
      !! The starting location in array `work` of array `beta0`.
   integer, intent(out) :: betaci
      !! The starting location in array `work` of array `betac`.
   integer, intent(out) :: betasi
      !! The starting location in array `work` of array `betas`.
   integer, intent(out) :: betani
      !! The starting location in array `work` of array `betan`.
   integer, intent(out) :: si
      !! The starting location in array `work` of array `s`.
   integer, intent(out) :: ssi
      !! The starting location in array `work` of array `ss`.
   integer, intent(out) :: ssfi
      !! The starting location in array `work` of array `ssf`.
   integer, intent(out) :: qrauxi
      !! The starting location in array `work` of array `qraux`.
   integer, intent(out) :: ui
      !! The starting location in array `work` of array `u`.
   integer, intent(out) :: fsi
      !! The starting location in array `work` of array `fs`.
   integer, intent(out) :: fjacbi
      !! The starting location in array `work` of array `fjacb`.
   integer, intent(out) :: we1i
      !! The starting location in array `work` of array `we1`.
   integer, intent(out) :: diffi
      !! The starting location in array `work` of array `diff`.
   integer, intent(out) :: deltsi
      !! The starting location in array `work` of array `deltas`.
   integer, intent(out) :: deltni
      !! The starting location in array `work` of array `deltan`.
   integer, intent(out) :: ti
      !! The starting location in array `work` of array `t`.
   integer, intent(out) :: tti
      !! The starting location in array `work` of array `tt`.
   integer, intent(out) :: omegai
      !! The starting location in array `work` of array `omega`.
   integer, intent(out) :: fjacdi
      !! The starting location in array `work` of array `fjacd`.
   integer, intent(out) :: wrk1i
      !! The starting location in array `work` of array `wrk1`.
   integer, intent(out) :: wrk2i
      !! The starting location in array `work` of array `wrk2`.
   integer, intent(out) :: wrk3i
      !! The starting location in array `work` of array `wrk3`.
   integer, intent(out) :: wrk4i
      !! The starting location in array `work` of array `wrk4`.
   integer, intent(out) :: wrk5i
      !! The starting location in array `work` of array `wrk5`.
   integer, intent(out) :: wrk6i
      !! The starting location in array `work` of array `wrk6`.
   integer, intent(out) :: wrk7i
      !! The starting location in array `work` of array `wrk7`.
   integer, intent(out) :: loweri
      !! The starting location in array `work` of array `lower`.
   integer, intent(out) :: upperi
      !! The starting location in array `work` of array `upper`.
   integer, intent(out) :: lwkmn
      !! The minimum acceptable length of vector `work`.

   ! Local scalars
   integer :: next

   ! Variable Definitions (alphabetically)
   !  ACTRSI:  The location in array WORK of variable ACTRS.
   !  ALPHAI:  The location in array WORK of variable ALPHA.
   !  BETACI:  The starting location in array WORK of array BETAC.
   !  BETANI:  The starting location in array WORK of array BETAN.
   !  BETASI:  The starting location in array WORK of array BETAS.
   !  BETA0I:  The starting location in array WORK of array BETA0.
   !  DELTAI:  The starting location in array WORK of array DELTA.
   !  DELTNI:  The starting location in array WORK of array DELTAN.
   !  DELTSI:  The starting location in array WORK of array DELTAS.
   !  DIFFI:   The starting location in array WORK of array DIFF.
   !  EPSI:    The starting location in array WORK of array EPS.
   !  EPSMAI:  The location in array WORK of variable EPSMAC.
   !  ETAI:    The location in array WORK of variable ETA.
   !  FJACBI:  The starting location in array WORK of array FJACB.
   !  FJACDI:  The starting location in array WORK of array FJACD.
   !  FNI:     The starting location in array WORK of array FN.
   !  FSI:     The starting location in array WORK of array FS.
   !  ISODR:   The variable designating whether the solution is by ODR (ISODR=TRUE) or by OLS (ISODR=FALSE).
   !  LDWE:    The leading dimension of array WE.
   !  LD2WE:   The second dimension of array WE.
   !  LWKMN:   The minimum acceptable length of vector work.
   !  M:       The number of columns of data in the explanatory variable.
   !  N:       The number of observations.
   !  NEXT:    The next available location with WORK.
   !  NP:      The number of function parameters.
   !  NQ:      The number of responses per observation.
   !  OLMAVI:  The location in array WORK of variable OLMAVG.
   !  OMEGAI:  The starting location in array WORK of array OMEGA.
   !  PARTLI:  The location in array WORK of variable PARTOL.
   !  PNORMI:  The location in array WORK of variable PNORM.
   !  PRERSI:  The location in array WORK of variable PRERS.
   !  QRAUXI:  The starting location in array WORK of array QRAUX.
   !  RCONDI:  The location in array WORK of variable RCONDI.
   !  RNORSI:  The location in array WORK of variable RNORMS.
   !  RVARI:   The location in array WORK of variable RVAR.
   !  SDI:     The starting location in array WORK of array SD.
   !  SI:      The starting location in array WORK of array S.
   !  SSFI:    The starting location in array WORK of array SSF.
   !  SSI:     The starting location in array WORK of array SS.
   !  SSTOLI:  The location in array WORK of variable SSTOL.
   !  TAUFCI:  The location in array WORK of variable TAUFAC.
   !  TAUI:    The location in array WORK of variable TAU.
   !  TI:      The starting location in array WORK of array T.
   !  TTI:     The starting location in array WORK of array TT.
   !  UI:      The starting location in array WORK of array U.
   !  VCVI:    The starting location in array WORK of array VCV.
   !  WE1I:    The starting location in array WORK of array WE1.
   !  WRK1I:   The starting location in array WORK of array WRK1.
   !  WRK2I:   The starting location in array WORK of array WRK2.
   !  WRK3I:   The starting location in array WORK of array WRK3.
   !  WRK4I:   The starting location in array WORK of array WRK4.
   !  WRK5I:   The starting location in array WORK of array WRK5.
   !  WRK6I:   The starting location in array WORK of array WRK6.
   !  WRK7I:   The starting location in array WORK of array WRK7.
   !  WSSI:    The location in array WORK of variable WSS.
   !  WSSDEI:  The location in array WORK of variable WSSDEL.
   !  WSSEPI:  The location in array work of variable WSSEPS.
   !  XPLUSI:  The starting location in array WORK of array XPLUSD.

   if (n .ge. 1 .and. m .ge. 1 .and. np .ge. 1 .and. nq .ge. 1 .and. &
       ldwe .ge. 1 .and. ld2we .ge. 1) then

      deltai = 1
      epsi = deltai + n*m
      xplusi = epsi + n*nq
      fni = xplusi + n*m
      sdi = fni + n*nq
      vcvi = sdi + np
      rvari = vcvi + np*np

      wssi = rvari + 1
      wssdei = wssi + 1
      wssepi = wssdei + 1
      rcondi = wssepi + 1
      etai = rcondi + 1
      olmavi = etai + 1

      taui = olmavi + 1
      alphai = taui + 1
      actrsi = alphai + 1
      pnormi = actrsi + 1
      rnorsi = pnormi + 1
      prersi = rnorsi + 1
      partli = prersi + 1
      sstoli = partli + 1
      taufci = sstoli + 1
      epsmai = taufci + 1
      beta0i = epsmai + 1

      betaci = beta0i + np
      betasi = betaci + np
      betani = betasi + np
      si = betani + np
      ssi = si + np
      ssfi = ssi + np
      qrauxi = ssfi + np
      ui = qrauxi + np
      fsi = ui + np

      fjacbi = fsi + n*nq

      we1i = fjacbi + n*np*nq

      diffi = we1i + ldwe*ld2we*nq

      next = diffi + nq*(np + m)

      if (isodr) then
         deltsi = next
         deltni = deltsi + n*m
         ti = deltni + n*m
         tti = ti + n*m
         omegai = tti + n*m
         fjacdi = omegai + nq*nq
         wrk1i = fjacdi + n*m*nq
         next = wrk1i + n*m*nq
      else
         deltsi = deltai
         deltni = deltai
         ti = deltai
         tti = deltai
         omegai = deltai
         fjacdi = deltai
         wrk1i = deltai
      end if

      wrk2i = next
      wrk3i = wrk2i + n*nq
      wrk4i = wrk3i + np
      wrk5i = wrk4i + m*m
      wrk6i = wrk5i + m
      wrk7i = wrk6i + n*nq*np
      loweri = wrk7i + 5*nq
      upperi = loweri + np
      next = upperi + np

      lwkmn = next
   else
      deltai = 1
      epsi = 1
      xplusi = 1
      fni = 1
      sdi = 1
      vcvi = 1
      rvari = 1
      wssi = 1
      wssdei = 1
      wssepi = 1
      rcondi = 1
      etai = 1
      olmavi = 1
      taui = 1
      alphai = 1
      actrsi = 1
      pnormi = 1
      rnorsi = 1
      prersi = 1
      partli = 1
      sstoli = 1
      taufci = 1
      epsmai = 1
      beta0i = 1
      betaci = 1
      betasi = 1
      betani = 1
      si = 1
      ssi = 1
      ssfi = 1
      qrauxi = 1
      fsi = 1
      ui = 1
      fjacbi = 1
      we1i = 1
      diffi = 1
      deltsi = 1
      deltni = 1
      ti = 1
      tti = 1
      fjacdi = 1
      omegai = 1
      wrk1i = 1
      wrk2i = 1
      wrk3i = 1
      wrk4i = 1
      wrk5i = 1
      wrk6i = 1
      wrk7i = 1
      loweri = 1
      upperi = 1
      lwkmn = 1
   end if

end subroutine dwinf

pure subroutine dwght(n, m, wt, ldwt, ld2wt, t, wtt)
   !! Scale matrix T using WT, i.e., compute WTT = WT*T
   ! Routines Called  (NONE)
   ! Date Written   860529   (YYMMDD)
   ! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp, zero

   integer, intent(in) :: n
      !! The number of rows of data in `t`.
   integer, intent(in) :: m
      !! The number of columns of data in `t`.
   integer, intent(in) :: ldwt
      !! The leading dimension of array `wt`.
   integer, intent(in) :: ld2wt
      !! The second dimension of array `wt`.
   real(kind=wp), intent(in) :: wt(:, :, :)
      !! The weights.
   real(kind=wp), intent(in) :: t(:, :)
      !! The array being scaled by `wt`.
   real(kind=wp), intent(out) :: wtt(:, :)
      !! The results of weighting array `t` by `wt`. Array `wtt` can be the same as `t` only if
      !! the arrays in `wt` are upper triangular with zeros below the diagonal.

   ! Local scalars
   real(kind=wp) :: temp
   integer :: i, j, k

   ! Variable Definitions (alphabetically)
   !  I:       An indexing variable.
   !  J:       An indexing variable.
   !  K:       An indexing variable.
   !  LDWT:    The leading dimension of array WT.
   !  LD2WT:   The second dimension of array WT.
   !  M:       The number of columns of data in T.
   !  N:       The number of rows of data in T.
   !  T:       The array being scaled by WT.
   !  TEMP:    A temporary scalar.
   !  WT:      The weights.
   !  WTT:     The results of weighting array T by WT. Array WTT can be the same as T only if
   !           the arrays in WT are upper triangular with zeros below the diagonal.

   if (n .eq. 0 .or. m .eq. 0) return

   if (wt(1, 1, 1) .ge. zero) then
      if (ldwt .ge. n) then
         if (ld2wt .ge. m) then
            ! WT is an N-array of M by M matrices
            do i = 1, n
               do j = 1, m
                  temp = zero
                  do k = 1, m
                     temp = temp + wt(i, j, k)*t(i, k)
                  end do
                  wtt(i, j) = temp
               end do
            end do
         else
            ! WT is an N-array of diagonal matrices
            do i = 1, n
               do j = 1, m
                  wtt(i, j) = wt(i, 1, j)*t(i, j)
               end do
            end do
         end if
      else
         if (ld2wt .ge. m) then
            ! WT is an M by M matrix
            do i = 1, n
               do j = 1, m
                  temp = zero
                  do k = 1, m
                     temp = temp + wt(1, j, k)*t(i, k)
                  end do
                  wtt(i, j) = temp
               end do
            end do
         else
            ! WT is a diagonal matrice
            do i = 1, n
               do j = 1, m
                  wtt(i, j) = wt(1, 1, j)*t(i, j)
               end do
            end do
         end if
      end if
   else
      ! WT is a scalar
      do j = 1, m
         do i = 1, n
            wtt(i, j) = abs(wt(1, 1, 1))*t(i, j)
         end do
      end do
   end if

end subroutine dwght

subroutine dvevtr &
   (m, nq, indx, &
    v, ldv, ld2v, e, lde, ve, ldve, ld2ve, vev, ldvev, &
    wrk5)
!! Compute  V*E*trans(V) for the (INDX)TH M by NQ array in V
! Routines Called  DSOLVE
! Date Written   910613   (YYMMDD)
! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp, zero

   integer, intent(in) :: m
      !! The number of columns of data in the independent variable.
   integer, intent(in) :: nq
      !! The number of responses per observation.
   integer, intent(in) :: indx
      !! The row in `v` in which the `m` by `nq` array is stored.
   integer, intent(in) :: ldv
      !! The leading dimension of array `v`.
   integer, intent(in) :: ld2v
      !! The second dimension of array `v`.
   integer, intent(in) :: lde
      !! The leading dimension of array `e`.
   integer, intent(in) :: ldve
      !! The leading dimension of array `ve`.
   integer, intent(in) :: ldvev
      !! The leading dimension of array `vev`.
   integer, intent(in) :: ld2ve
      !! The second dimension of array `ve`.
   real(kind=wp), intent(in) :: v(ldv, ld2v, nq)
      !! An array of `nq` by `m` matrices.
   real(kind=wp), intent(in) :: e(lde, m)
      !! The `m` by `m` matrix of the factors so `ete = (d**2 + alpha*t**2)`.
   real(kind=wp), intent(out) :: ve(ldve, ld2ve, m)
      !! The `nq` by `m` array `ve = v * inv(e)`.
   real(kind=wp), intent(out) :: vev(ldvev, nq)
      !! The `nq` by `nq` array `vev = v * inv(ete) * trans(v)`.
   real(kind=wp), intent(out) :: wrk5(m)
      !! An `m` work vector.

   ! Local scalars
   integer :: j, l1, l2

   ! External subroutines
   external :: dsolve

   ! Variable Definitions (alphabetically)
   !  INDX:    The row in V in which the M by NQ array is stored.
   !  J:       An indexing variable.
   !  LDE:     The leading dimension of array E.
   !  LDV:     The leading dimension of array V.
   !  LDVE:    The leading dimension of array VE.
   !  LDVEV:   The leading dimension of array VEV.
   !  LD2V:    The second dimension of array V.
   !  L1:      An indexing variable.
   !  L2:      An indexing variable.
   !  M:       The number of columns of data in the independent variable.
   !  NQ:      The number of responses per observation.
   !  E:       The M by M matrix of the factors so ETE = (D**2 + ALPHA*T**2).
   !  V:       An array of NQ by M matrices.
   !  VE:      The NQ by M array VE = V * inv(E)
   !  VEV:     The NQ by NQ array VEV = V * inv(ETE) * trans(V).
   !  WRK5:    An M work vector.

   if (nq .eq. 0 .or. m .eq. 0) return

   do l1 = 1, nq
      do j = 1, m
         wrk5(j) = v(indx, j, l1)
      end do
      call dsolve(m, e, lde, wrk5, 4)
      do j = 1, m
         ve(indx, l1, j) = wrk5(j)
      end do
   end do

   do l1 = 1, nq
      do l2 = 1, l1
         vev(l1, l2) = zero
         do j = 1, m
            vev(l1, l2) = vev(l1, l2) + ve(indx, l1, j)*ve(indx, l2, j)
         end do
         vev(l2, l1) = vev(l1, l2)
      end do
   end do

end subroutine dvevtr

subroutine dunpac(n2, v1, v2, ifix)
!! Copy the elements of `v1` into the locations of `v2` which are unfixed
! Routines Called  DCOPY
! Date Written   860529   (YYMMDD)
! Revision Date  920304   (YYMMDD)

   use odrpack_kinds, only: wp

   integer, intent(in) :: n2
      !! The number of items in `v2`.
   real(kind=wp), intent(in) :: v1(n2)
      !! The vector of the unfixed items.
   real(kind=wp), intent(out) :: v2(n2)
      !! The vector of the fixed and unfixed items into which the elements of `v1` are to be inserted.
   integer, intent(in) :: ifix(n2)
      !! The values designating whether the elements of `v2` are fixed at their input values or not.

   ! Local scalars
   integer :: i, n1

   ! External subroutines
   external :: dcopy

   ! Variable Definitions (alphabetically)
   !  I:       An indexing variable.
   !  IFIX:    The values designating whether the elements of V2 are fixed at their input values or not.
   !  N1:      The number of items in V1.
   !  N2:      The number of items in V2.
   !  V1:      The vector of the unfixed items.
   !  V2:      The vector of the fixed and unfixed items into which the elements of V1 are to be inserted.

   n1 = 0
   if (ifix(1) .ge. 0) then
      do i = 1, n2
         if (ifix(i) .ne. 0) then
            n1 = n1 + 1
            v2(i) = v1(n1)
         end if
      end do
   else
      n1 = n2
      call dcopy(n2, v1, 1, v2, 1)
   end if

end subroutine dunpac

subroutine dsolve(n, t, ldt, b, job)
!! Solve systems of the form:
!!  T * X = B  or  trans(T) * X = B
!! where T is an upper or lower triangular matrix of order N, and the solution X overwrites
!! the RHS B. (adapted from LINPACK subroutine DTRSL)
!! References:
!!  Dongarra J.J., Bunch J.R., Moler C.B., Stewart G.W., *LINPACK Users Guide*, SIAM, 1979.
! Routines Called  DAXPY,DDOT
! Date Written   920220   (YYMMDD)
! Revision Date  920619   (YYMMDD)

   use odrpack_kinds, only: wp, zero

   integer, intent(in) :: n
      !! The number of rows and columns of data in array `t`.
   real(kind=wp), intent(in) :: t(ldt, n)
      !! The upper or lower tridiagonal system.
   integer, intent(in) :: ldt
      !! The leading dimension of array `t`.
   real(kind=wp), intent(inout) :: b(n)
      !! On input: the right hand side; On exit: the solution.
   integer, intent(in) :: job
      !! What kind of system is to be solved:
      !!   1   Solve T*X=B, T lower triangular,
      !!   2   Solve T*X=B, T upper triangular,
      !!   3   Solve trans(T)*X=B, T lower triangular,
      !!   4   Solve trans(T)*X=B, T upper triangular.

   ! Local scalars
   real(kind=wp) :: temp
   integer :: j1, j, jn

   ! External functions
   real(kind=wp), external :: ddot

   ! External subroutines
   external :: daxpy

   ! Variable Definitions (alphabetically)
   !  B:       On input:  the right hand side;  On exit:  the solution
   !  J1:      The first nonzero entry in T.
   !  J:       An indexing variable.
   !  JN:      The last nonzero entry in T.
   !  JOB:     What kind of system is to be solved, where if JOB is
   !           1   Solve T*X=B, T lower triangular,
   !           2   Solve T*X=B, T upper triangular,
   !           3   Solve trans(T)*X=B, T lower triangular,
   !           4   Solve trans(T)*X=B, T upper triangular.
   !  LDT:     The leading dimension of array T.
   !  N:       The number of rows and columns of data in array T.
   !  T:       The upper or lower tridiagonal system.

   !  Find first nonzero diagonal entry in T
   j1 = 0
   do j = 1, n
      if (j1 .eq. 0 .and. t(j, j) .ne. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
         j1 = j
      elseif (t(j, j) .eq. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
         b(j) = zero
      end if
   end do
   if (j1 .eq. 0) return
   !  Find last nonzero diagonal entry in T
   jn = 0
   do j = n, j1, -1
      if (jn .eq. 0 .and. t(j, j) .ne. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
         jn = j
      elseif (t(j, j) .eq. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
         b(j) = zero
      end if
   end do

   if (job .eq. 1) then
      !  Solve T*X=B for T lower triangular
      b(j1) = b(j1)/t(j1, j1)
      do j = j1 + 1, jn
         temp = -b(j - 1)
         call daxpy(jn - j + 1, temp, t(j, j - 1), 1, b(j), 1)
         if (t(j, j) .ne. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
            b(j) = b(j)/t(j, j)
         else
            b(j) = zero
         end if
      end do

   elseif (job .eq. 2) then
      !  Solve T*X=B for T upper triangular.
      b(jn) = b(jn)/t(jn, jn)
      do j = jn - 1, j1, -1
         temp = -b(j + 1)
         call daxpy(j, temp, t(1, j + 1), 1, b(1), 1)
         if (t(j, j) .ne. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
            b(j) = b(j)/t(j, j)
         else
            b(j) = zero
         end if
      end do

   elseif (job .eq. 3) then
      !  Solve trans(T)*X=B for T lower triangular.
      b(jn) = b(jn)/t(jn, jn)
      do j = jn - 1, j1, -1
         b(j) = b(j) - ddot(jn - j + 1, t(j + 1, j), 1, b(j + 1), 1)
         if (t(j, j) .ne. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
            b(j) = b(j)/t(j, j)
         else
            b(j) = zero
         end if
      end do

   elseif (job .eq. 4) then
      !  Solve trans(T)*X=B for T upper triangular.
      b(j1) = b(j1)/t(j1, j1)
      do j = j1 + 1, jn
         b(j) = b(j) - ddot(j - 1, t(1, j), 1, b(1), 1)
         if (t(j, j) .ne. zero) then
   !!! FPT - 3087 REAL or COMPLEX quantity tested for exact equality/inequality
            b(j) = b(j)/t(j, j)
         else
            b(j) = zero
         end if
      end do
   end if

end subroutine dsolve
