! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Alperen Yavuz
! chain_general.inc -- chi(Disc A) mod q for lanes whose Euclid chain is not
! generic (some degree drops by more than one).  Included twice by
! pdfilter.f90, for the small primes (16-bit, R = 2^16) and the large primes
! (32-bit, R = 2^32), with GNAME, GK (integer kind), GQ, GQINV, GONE, GNEGR,
! GEXP, GVMAX and the arithmetic macros
!     GMUL(a, b)        Montgomery product a*b*R^-1 mod q
!     GDIFF(a, b, c, d) (a*b - c*d)*R^-1 mod q
!     GZERO(x)          x is 0 mod q
! set by pdfilter.f90 (see there for the ranges of the values).
!
! Polynomials are stored left-aligned (index 0 = leading coefficient, zero
! padded), so one elementary pseudo-division step
!     U' = lc(V)*U - lc(U)*x^(deg U - deg V)*V
! is the same formula  U'(i) = lc(V)*U(i+1) - lc(U)*V(i+1)  in every lane,
! whatever the degrees.  Leading zeros of U' are shifted out lane by lane.
!
! Bookkeeping, with Res(V,P) = lc(V)^deg(P) * prod_{V(b)=0} P(b):
!     Res(V,U') = lc(V)^(deg U' - deg U + deg V) * Res(V,U),
! and swapping U and V only changes the sign (chi(-1) = 1 for q = 1 mod 4).
! The odd powers of lc(V) are multiplied into mm; at the end (V constant c)
! Res = c^deg(U), so chi(Disc A) = chi(mm * c^(deg U mod 2)).  A zero U'
! means Res = 0 mod q: the lane is kept.  The factors R^-1 of the Montgomery
! products are squares and do not change chi.
!
! U, V and the new remainder W live in three slots of p.  While all lanes move
! in lockstep (same degrees), "U = W" and "swap U, V" just exchange slot
! indices; lanes that diverge are handled with masked merges.

  subroutine GNAME(a, j, keep)
    integer(GK), intent(in) :: a(BLK, 0:10)
    integer, intent(in) :: j
    logical, intent(out) :: keep(BLK)
    integer(GK) :: p(BLK, 0:11, 0:2)
    integer(GK) :: du(BLK), dv(BLK), dw(BLK), mm(BLK), lu(BLK), lv(BLK), ev(BLK)
    integer(GK) :: base(BLK), acc(BLK), tmp(BLK)
    ! Lane flags as 0/1 integers: gfortran does not vectorise logical arrays.
    integer(GK) :: act(BLK), done(BLK), resz(BLK), zf(BLK), sw(BLK)
    integer(GK) :: q, qinv, negr
    integer(i8) :: e
    integer :: i, it, imax, iu, iv, iw, t

    q = GQ(j); qinv = GQINV(j); negr = GNEGR(j)
    iu = 0; iv = 1; iw = 2
    ! U = A (degree 10, monic), V = A' (degree 9, leading coefficient 10).
    do i = 0, 10
      p(:, i, iu) = a(:, 10-i)
    end do
    p(:, 11, iu) = 0
    do i = 0, 9
      p(:, i, iv) = GMUL(int(10-i, GK), a(:, 10-i))
    end do
    p(:, 10:11, iv) = 0
    du = 10; dv = 9
    mm = 1
    done = 0; resz = 0

    do it = 1, 40
      act = int(1, GK) - done
      if (GVMAX(act) == 0) exit
      imax = int(GVMAX(du*act))
      lu = p(:, 0, iu); lv = p(:, 0, iv)
      do i = 0, imax - 1
        p(:, i, iw) = GDIFF(lv, p(:, i+1, iu), lu, p(:, i+1, iv))
      end do
      p(:, imax:11, iw) = 0
      dw = du - int(1, GK)
      ! Shift out leading zeros (a zero polynomial ends at degree -1).
      do
        zf = merge(act, int(0, GK), GZERO(p(:, 0, iw)) .and. dw >= 0)
        if (GVMAX(zf) == 0) exit
        do i = 0, imax - 1
          p(:, i, iw) = merge(p(:, i+1, iw), p(:, i, iw), zf > 0)
        end do
        p(:, imax, iw) = merge(int(0, GK), p(:, imax, iw), zf > 0)
        dw = dw - zf
      end do
      ! U' = 0: resultant vanishes mod q.
      tmp = merge(act, int(0, GK), dw < 0)
      resz = max(resz, tmp)
      done = max(done, tmp)
      act = act - tmp
      ! chi bookkeeping for lc(V)^(deg U' - deg U + deg V).
      ev = abs(dw - du + dv)
      tmp = GMUL(mm, lv)
      mm = merge(tmp, mm, act > 0 .and. iand(ev, int(1, GK)) /= 0)
      ! U = U'; swap U and V when deg U < deg V.
      sw = merge(act, int(0, GK), dw < dv)
      ! Lockstep: every lane active with the same deg U' and deg V.
      tmp = (int(1, GK) - act) + abs(dw - dw(1)) + abs(dv - dv(1))
      if (GVMAX(tmp) == 0) then
        t = iu; iu = iw; iw = t
        du = dw
        if (sw(1) > 0) then
          t = iu; iu = iv; iv = t
          tmp = du; du = dv; dv = tmp
        end if
      else
        do i = 0, imax
          p(:, i, iu) = merge(p(:, i, iw), p(:, i, iu), act > 0)
        end do
        du = merge(dw, du, act > 0)
        do i = 0, imax
          tmp = p(:, i, iu)
          p(:, i, iu) = merge(p(:, i, iv), p(:, i, iu), sw > 0)
          p(:, i, iv) = merge(tmp, p(:, i, iv), sw > 0)
        end do
        tmp = du
        du = merge(dv, du, sw > 0)
        dv = merge(tmp, dv, sw > 0)
      end if
      ! V constant c: Res = c^deg U.
      tmp = GMUL(mm, p(:, 0, iv))
      mm = merge(tmp, mm, act > 0 .and. dv == 0 .and. iand(du, int(1, GK)) /= 0)
      done = max(done, merge(act, int(0, GK), dv == 0))
    end do

    ! Euler criterion on mm in the Montgomery domain: acc ends at chi*R mod q,
    ! and chi = -1 is the residue -R (negr or negr - q).
    base = mm
    acc = GONE(j)
    e = GEXP(j)
    do while (e > 0)
      if (iand(e, 1_i8) == 1) acc = GMUL(acc, base)
      base = GMUL(base, base)
      e = shiftr(e, 1)
    end do
    keep = .not. (done > 0 .and. resz == 0 .and. (acc == negr .or. acc == negr - q))
  end subroutine GNAME
