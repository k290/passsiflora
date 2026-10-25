! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Alperen Yavuz
! pdfilter.f90 -- exact modular pre-filter for the NumberFields polDisc test.
!
! The upstream test asks: after removing every factor of the primes in S from
! |polDisc|, is what remains a perfect square?  Up to rational square factors
! (B's content c gives c^10, the last sub-resultant step gives (B0/h)^(2k)),
! polDisc equals Disc(A), so the question is the same for Disc(A).
!
! For a prime q with q == 1 (mod 4) and every p in S a quadratic residue mod q,
! a passing polynomial must have Legendre(Disc(A) mod q) in {0, 1}.  If the
! Legendre symbol is -1 the polynomial can never pass, so it is rejected here.
! Everything else is kept and handed to the original GMP test unchanged,
! which keeps the final flags bit-identical to the stock application.
!
! Disc(A) mod q is evaluated with a generic Euclid chain over F_q using
! pseudo-remainders (two rounds per step, degree drops by exactly one).  In
! that generic case every leading-coefficient power that appears has an even
! exponent, so chi(Disc(A)) = chi(final constant) and no inverses are needed.
!
! All arithmetic is integer arithmetic with signed Montgomery multiplication:
! no floating point is used, so the results do not depend on the rounding
! mode, on x87 or SSE arithmetic, or on flags such as -ffast-math.
!  * Small primes 2^8 < q < 2^10: 16-bit lanes (sixteen per AVX2 register),
!    R = 2^16 (mont).  Values are kept in (-2q, 2q), a difference of two
!    Montgomery products; a value is 0 mod q iff it is 0 or +-q.
!  * Large primes 2^24 < q < 2^25: 32-bit lanes, R = 2^32 (sredc), products
!    formed as 64-bit widening products.  A difference a*b - c*d is reduced
!    once; values are in (-q, q), where 0 is the only multiple of q.
!  * Every Montgomery product carries a factor R^-1, and R = 2^16 = (2^8)^2
!    and R = 2^32 = (2^16)^2 are squares, so the factors never change a
!    Legendre symbol.  Each polynomial of a chain is its exact counterpart
!    times a power of R^-1 (both terms of every difference carry the same
!    power), so the zero tests are exact.  Euler's criterion runs in the
!    Montgomery domain: acc starts at R mod q and ends at chi*R mod q, so -1
!    is the residue -R, i.e. negr or negr - q with negr = -R mod q in [0, q).
!  * Loading: every coefficient c is split once per buffer as
!    c = p0 + p1*2^21 + p2*2^42 (p0, p1 in [0, 2^21), p2 = floor(c/2^42) in
!    [-2^21, 2^21)); per prime, p0*d0 + p1*d1 + p2*d2 with d_i = 2^(21i)*2^32
!    mod q is reduced once by sredc, which gives c mod q itself, in (-q, q).
! The first stage uses the small primes; each prime rejects about half of the
! candidates, and the few survivors continue with the large primes.  All inner
! loops vectorise across polynomials (BLK lanes).

! Built twice (see the build scripts, -cpp): once per instruction-set level,
! with a distinct module name and C symbols; the C++ side picks one at run time.
#ifndef PDF_MOD
#  define PDF_MOD pdfilter
#  define PDF_INIT_NAME "pdf_init"
#  define PDF_FILTER_NAME "pdf_filter"
#endif

module PDF_MOD
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: pdf_init, pdf_filter

  integer, parameter :: i2 = c_int16_t, i4 = c_int32_t, i8 = c_int64_t
  integer, parameter :: MAXQ = 24, MAXS = 6, BLK = 64

  ! Large primes (32-bit stage, R = 2^32): qinv = q^-1 mod 2^32 as a signed
  ! 32-bit value, one = R mod q, negr = -R mod q, ld(i) = 2^(21i)*R mod q.
  integer :: nq = 0
  integer(i4) :: qv(MAXQ), qinvv(MAXQ), onev(MAXQ), negrv(MAXQ), ldv(0:2, MAXQ)
  integer(i8) :: qexp(MAXQ)
  ! Small primes (16-bit stage, R = 2^16): the same constants for R = 2^16;
  ! qsl, qinvsl, lds are the 32-bit constants for loading (sredc, R = 2^32).
  integer :: ns = 0
  integer(i2) :: qs(MAXS), qinvs(MAXS), ones(MAXS), negrs(MAXS)
  integer(i4) :: qsl(MAXS), qinvsl(MAXS), lds(0:2, MAXS)
  integer(i8) :: qsexp(MAXS)

  integer, allocatable, save :: idx(:), nidx(:), uidx(:)
  integer(i4), allocatable, save :: p0A(:, :), p1A(:, :), p2A(:, :)
  logical, save :: debug = .false., allgen = .false.
  ! Per prime stage: go straight to chain_general when most lanes of the last
  ! probe left the generic chain (structurally degenerate families).
  logical, save :: prefgen(MAXS + MAXQ) = .false.
  integer, save :: ncall = 0

contains

  pure integer(i8) function powmod_i(b, e, m) result(r)
    integer(i8), intent(in) :: b, e, m
    integer(i8) :: x, ee
    r = 1
    x = modulo(b, m)
    ee = e
    do while (ee > 0)
      if (iand(ee, 1_i8) == 1) r = modulo(r*x, m)
      x = modulo(x*x, m)
      ee = shiftr(ee, 1)
    end do
  end function powmod_i

  pure logical function is_prime(n)
    integer(i8), intent(in) :: n
    integer(i8) :: d
    is_prime = .false.
    if (n < 2) return
    if (mod(n, 2_i8) == 0) then
      is_prime = (n == 2)
      return
    end if
    d = 3
    do while (d*d <= n)
      if (mod(n, d) == 0) return
      d = d + 2
    end do
    is_prime = .true.
  end function is_prime

  ! q is usable when q == 1 (mod 4) is prime and every p in S is a nonzero
  ! quadratic residue mod q.
  logical function usable(q, numP, pSet)
    integer(i8), intent(in) :: q
    integer(c_int), intent(in) :: numP
    integer(c_int), intent(in) :: pSet(*)
    integer(i8) :: p
    integer :: k
    usable = .false.
    if (mod(q, 4_i8) /= 1 .or. .not. is_prime(q)) return
    do k = 1, numP
      p = int(pSet(k), i8)
      if (mod(p, q) == 0 .or. powmod_i(p, (q-1)/2, q) /= 1) return
    end do
    usable = .true.
  end function usable

  ! Inverse of an odd q modulo 2^32 (extended Euclid; every intermediate stays
  ! below 2^33).
  pure integer(i8) function inv32(q)
    integer(i8), intent(in) :: q
    integer(i8) :: r0, r1, s0, s1, t, d
    r0 = 2_i8**32; r1 = q; s0 = 0; s1 = 1
    do while (r1 /= 0)
      d = r0 / r1
      t = r0 - d*r1; r0 = r1; r1 = t
      t = s0 - d*s1; s0 = s1; s1 = t
    end do
    inv32 = modulo(s0, 2_i8**32)
  end function inv32

  ! x mod 2^b as a signed b-bit value.
  pure integer(i8) function signedmod(x, b)
    integer(i8), intent(in) :: x
    integer, intent(in) :: b
    signedmod = modulo(x, 2_i8**b)
    if (signedmod >= 2_i8**(b-1)) signedmod = signedmod - 2_i8**b
  end function signedmod

  ! Choose the filter primes for the set S of primes divided out by the test.
  subroutine pdf_init(numP, pSet) bind(C, name=PDF_INIT_NAME)
    integer(c_int), value :: numP
    integer(c_int), intent(in) :: pSet(*)
    integer(i8) :: q, qi, r
    character(len=8) :: ev
    integer :: evlen, i

    call get_environment_variable('NF_PDFDEBUG', ev, evlen)
    debug = (evlen > 0)
    ! NF_GENERAL_ALL=1 sends every lane through chain_general (testing aid).
    call get_environment_variable('NF_GENERAL_ALL', ev, evlen)
    allgen = (evlen > 0)

    ns = 0
    q = 2_i8**10 - 1
    do while (ns < MAXS .and. q > 256)
      if (usable(q, numP, pSet)) then
        ns = ns + 1
        qi = inv32(q)
        qs(ns)     = int(q, i2)
        qinvs(ns)  = int(signedmod(qi, 16), i2)
        ones(ns)   = int(modulo(2_i8**16, q), i2)
        negrs(ns)  = int(modulo(-2_i8**16, q), i2)
        qsl(ns)    = int(q, i4)
        qinvsl(ns) = int(signedmod(qi, 32), i4)
        r = modulo(2_i8**32, q)
        do i = 0, 2
          lds(i, ns) = int(r, i4)
          r = modulo(r*2_i8**21, q)
        end do
        qsexp(ns) = (q - 1) / 2
      end if
      q = q - 2
    end do

    nq = 0
    q = 2_i8**25 - 1
    do while (nq < MAXQ .and. q > 2_i8**24)
      if (usable(q, numP, pSet)) then
        nq = nq + 1
        qv(nq)    = int(q, i4)
        qinvv(nq) = int(signedmod(inv32(q), 32), i4)
        onev(nq)  = int(modulo(2_i8**32, q), i4)
        negrv(nq) = int(modulo(-2_i8**32, q), i4)
        r = modulo(2_i8**32, q)
        do i = 0, 2
          ldv(i, nq) = int(r, i4)
          r = modulo(r*2_i8**21, q)
        end do
        qexp(nq) = (q - 1) / 2
      end if
      q = q - 2
    end do
  end subroutine pdf_init

  ! Signed Montgomery multiplication, R = 2^16, for odd 2^8 < q < 2^10 and
  ! qinv = q^-1 mod 2^16 as a signed 16-bit value.  For |a|, |b| < 2q: with
  ! x = a*b and m = (x*qinv) mod 2^16 in [-2^15, 2^15), x - m*q is divisible by
  ! 2^16, so (x - m*q)/2^16 = hi(x) - hi(m*q) exactly (hi = 16-bit multiply-
  ! high, the low halves of x and m*q are equal).  It is congruent to a*b*2^-16
  ! and lies in (-q, q), as |x| + 2^15*q < 4q^2 + 2^15*q < 2^16*q for q < 2^13.
  ! The truncations to 16 bits are conversions (modulo 2^16), not overflows.
  elemental integer(i2) function mont(a, b, q, qinv)
    integer(i2), intent(in) :: a, b, q, qinv
    integer(i2) :: lo, m, h1, h2
    lo = int(int(a, i4)*int(b, i4), i2)
    m  = int(int(lo, i4)*int(qinv, i4), i2)
    h1 = int(shifta(int(a, i4)*int(b, i4), 16), i2)
    h2 = int(shifta(int(m, i4)*int(q, i4), 16), i2)
    mont = h1 - h2
  end function mont

  ! Signed Montgomery reduction, R = 2^32, for odd q < 2^31 and qinv = q^-1 mod
  ! 2^32 as a signed 32-bit value.  For |x| < 2^31*q: m = (x*qinv) mod 2^32 as a
  ! signed 32-bit value (the 64-bit product of the signed low half of x and
  ! qinv is exact, below 2^62, and is truncated), x - m*q is divisible by 2^32,
  ! and u = (x - m*q)/2^32 lies in (-q, q), since |x| + 2^31*q < 2^32*q; u is
  ! congruent to x*2^-32 mod q.  Callers form x from 64-bit products of 32-bit
  ! values, which gfortran vectorises as widening multiplies.
  elemental integer(i4) function sredc(x, q, qinv)
    integer(i8), intent(in) :: x
    integer(i4), intent(in) :: q, qinv
    integer(i4) :: m
    m = int(int(int(x, i4), i8)*int(qinv, i8), i4)
    sredc = int(shifta(x - int(m, i8)*int(q, i8), 32), i4)
  end function sredc

  ! c mod q in (-q, q) for c = p0 + p1*2^21 + p2*2^42, d_i = 2^(21i)*2^32 mod q:
  ! the sum is congruent to c*2^32 and below 3*2^21*q < 2^31*q in absolute value.
  elemental integer(i4) function ldres(p0, p1, p2, d0, d1, d2, q, qinv)
    integer(i4), intent(in) :: p0, p1, p2, d0, d1, d2, q, qinv
    ldres = sredc(int(p0, i8)*int(d0, i8) + int(p1, i8)*int(d1, i8) + int(p2, i8)*int(d2, i8), q, qinv)
  end function ldres

  ! Small-prime versions of load_block / chain_block (same algorithm).
  subroutine load_block_s(ids, c, j, a)
    integer, intent(in) :: c, j
    integer, intent(in) :: ids(c)
    integer(i2), intent(out) :: a(BLK, 0:10)
    integer(i4) :: q, qinv, d0, d1, d2
    integer :: l, k

    q = qsl(j); qinv = qinvsl(j); d0 = lds(0, j); d1 = lds(1, j); d2 = lds(2, j)
    do k = 0, 9
      do l = 1, c
        a(l, k) = int(ldres(p0A(ids(l), k), p1A(ids(l), k), p2A(ids(l), k), d0, d1, d2, q, qinv), i2)
      end do
      a(c+1:BLK, k) = 0
    end do
    a(:, 10) = 1
  end subroutine load_block_s

  subroutine chain_block_s(a, j, keep, ok)
    integer(i2), intent(in) :: a(BLK, 0:10)
    integer, intent(in) :: j
    logical, intent(out) :: keep(BLK), ok(BLK)
    integer(i2) :: p(BLK, 0:10, 0:2)
    integer(i2) :: lc(BLK), top(BLK), t(BLK), base(BLK), acc(BLK)
    integer(i2) :: q, qinv, negr
    integer(i8) :: e
    integer :: m, k, l, ia, ib, ir, tmp

    q = qs(j); qinv = qinvs(j); negr = negrs(j)
    ia = 0; ib = 1; ir = 2
    p(:, :, ia) = a
    do k = 1, 10
      p(:, k-1, ib) = mont(a(:, k), int(k, i2), q, qinv)
    end do
    ok = .true.

    do m = 9, 1, -1
      lc  = p(:, m, ib)
      top = p(:, m+1, ia)
      p(:, 0, ir) = mont(lc, p(:, 0, ia), q, qinv)
      do k = 1, m
        do l = 1, BLK
          p(l, k, ir) = mont(lc(l), p(l, k, ia), q, qinv) - mont(top(l), p(l, k-1, ib), q, qinv)
        end do
      end do
      t = p(:, m, ir)
      do k = 0, m-1
        do l = 1, BLK
          p(l, k, ir) = mont(lc(l), p(l, k, ir), q, qinv) - mont(t(l), p(l, k, ib), q, qinv)
        end do
      end do
      if (m >= 2) ok = ok .and. (p(:, m-1, ir) /= 0 .and. abs(p(:, m-1, ir)) /= q)
      tmp = ia; ia = ib; ib = ir; ir = tmp
    end do

    base = p(:, 0, ib)
    acc = ones(j)
    e = qsexp(j)
    do while (e > 0)
      if (iand(e, 1_i8) == 1) acc = mont(acc, base, q, qinv)
      base = mont(base, base, q, qinv)
      e = shiftr(e, 1)
    end do

    keep = .not. (ok .and. (acc == negr .or. acc == negr - q))
  end subroutine chain_block_s

  ! Load and reduce a block of polynomials into F_q.
  ! a(l,k) = coefficient of x^k of polynomial idx(l), k = 0..10, in (-q, q).
  ! p0A/p1A/p2A hold every coefficient split into 21-bit pieces, computed once
  ! per buffer so each prime only gathers and reduces.
  subroutine load_block(ids, c, j, a)
    integer, intent(in) :: c, j
    integer, intent(in) :: ids(c)
    integer(i4), intent(out) :: a(BLK, 0:10)
    integer(i4) :: q, qinv, d0, d1, d2
    integer :: l, k

    q = qv(j); qinv = qinvv(j); d0 = ldv(0, j); d1 = ldv(1, j); d2 = ldv(2, j)
    do k = 0, 9
      do l = 1, c
        a(l, k) = ldres(p0A(ids(l), k), p1A(ids(l), k), p2A(ids(l), k), d0, d1, d2, q, qinv)
      end do
      a(c+1:BLK, k) = 0
    end do
    a(:, 10) = 1
  end subroutine load_block

  ! keep(l) = .false. only when Disc(A_l) is proven to be a non-residue mod q.
  ! ok(l) = .false. when the chain left the generic path (some leading
  ! coefficient vanished mod q); pdf_filter then reruns that lane through
  ! chain_general, which follows any degree sequence.
  subroutine chain_block(a, j, keep, ok)
    integer(i4), intent(in) :: a(BLK, 0:10)
    integer, intent(in) :: j
    logical, intent(out) :: keep(BLK), ok(BLK)
    ! Three polynomial slots rotate roles (A, B, R) instead of being copied.
    integer(i4) :: p(BLK, 0:10, 0:2)
    integer(i4) :: lc(BLK), top(BLK), t(BLK), base(BLK), acc(BLK)
    integer(i4) :: q, qinv, negr
    integer(i8) :: e
    integer :: m, k, l, ia, ib, ir, tmp

    q = qv(j); qinv = qinvv(j); negr = negrv(j)
    ia = 0; ib = 1; ir = 2
    p(:, :, ia) = a
    do k = 1, 10
      p(:, k-1, ib) = sredc(int(k, i8)*int(a(:, k), i8), q, qinv)
    end do
    ok = .true.

    do m = 9, 1, -1
      ! A = slot ia (degree m+1), B = slot ib (degree m).
      lc  = p(:, m, ib)
      top = p(:, m+1, ia)
      ! Round 1: R = lc*A - top*x*B, degree <= m.
      p(:, 0, ir) = sredc(int(lc, i8)*int(p(:, 0, ia), i8), q, qinv)
      do k = 1, m
        do l = 1, BLK
          p(l, k, ir) = sredc(int(lc(l), i8)*int(p(l, k, ia), i8) - int(top(l), i8)*int(p(l, k-1, ib), i8), q, qinv)
        end do
      end do
      ! Round 2: R2 = lc*R - R_m*B, degree <= m-1 (in place).
      t = p(:, m, ir)
      do k = 0, m-1
        do l = 1, BLK
          p(l, k, ir) = sredc(int(lc(l), i8)*int(p(l, k, ir), i8) - int(t(l), i8)*int(p(l, k, ib), i8), q, qinv)
        end do
      end do
      if (m >= 2) ok = ok .and. (p(:, m-1, ir) /= 0)
      ! Next step: A = B, B = R2; the old A slot becomes scratch.
      tmp = ia; ia = ib; ib = ir; ir = tmp
    end do

    ! Euler criterion on the final constant, in the Montgomery domain.
    base = p(:, 0, ib)
    acc = onev(j)
    e = qexp(j)
    do while (e > 0)
      if (iand(e, 1_i8) == 1) acc = sredc(int(acc, i8)*int(base, i8), q, qinv)
      base = sredc(int(base, i8)*int(base, i8), q, qinv)
      e = shiftr(e, 1)
    end do

    keep = .not. (ok .and. (acc == negr .or. acc == negr - q))
  end subroutine chain_block

  ! max over the BLK lanes, 8 lanes at a time.
  pure integer(i2) function vmax_s(x)
    integer(i2), intent(in) :: x(BLK)
    integer(i2) :: m(8)
    integer :: i
    m = x(1:8)
    do i = 9, BLK, 8
      m = max(m, x(i:i+7))
    end do
    vmax_s = max(max(max(m(1), m(2)), max(m(3), m(4))), max(max(m(5), m(6)), max(m(7), m(8))))
  end function vmax_s

  pure integer(i4) function vmax_l(x)
    integer(i4), intent(in) :: x(BLK)
    integer(i4) :: m(8)
    integer :: i
    m = x(1:8)
    do i = 9, BLK, 8
      m = max(m, x(i:i+7))
    end do
    vmax_l = max(max(max(m(1), m(2)), max(m(3), m(4))), max(max(m(5), m(6)), max(m(7), m(8))))
  end function vmax_l

  ! General (non-generic degree sequence) chain, once per stage.
#define GNAME chain_general_s
#define GK i2
#define GQ qs
#define GQINV qinvs
#define GONE ones
#define GNEGR negrs
#define GEXP qsexp
#define GVMAX vmax_s
#define GMUL(a, b) mont(a, b, q, qinv)
#define GDIFF(a, b, c, d) mont(a, b, q, qinv) - mont(c, d, q, qinv)
#define GZERO(x) (x == 0 .or. abs(x) == q)
#include "chain_general.inc"
#undef GNAME
#undef GK
#undef GQ
#undef GQINV
#undef GONE
#undef GNEGR
#undef GEXP
#undef GVMAX
#undef GMUL
#undef GDIFF
#undef GZERO
#define GNAME chain_general
#define GK i4
#define GQ qv
#define GQINV qinvv
#define GONE onev
#define GNEGR negrv
#define GEXP qexp
#define GVMAX vmax_l
#define GMUL(a, b) sredc(int(a, i8)*int(b, i8), q, qinv)
#define GDIFF(a, b, c, d) sredc(int(a, i8)*int(b, i8) - int(c, i8)*int(d, i8), q, qinv)
#define GZERO(x) (x == 0)
#include "chain_general.inc"

  ! One block through prime stage j (small primes first, then large ones),
  ! with the fast generic chain or, for lanes that left it, the general one.
  subroutine run_block(ids, c, j, general, keep, ok)
    integer, intent(in) :: c, j
    integer, intent(in) :: ids(c)
    logical, intent(in) :: general
    logical, intent(out) :: keep(BLK), ok(BLK)
    integer(i4) :: a(BLK, 0:10)
    integer(i2) :: as(BLK, 0:10)
    if (j <= ns) then
      call load_block_s(ids, c, j, as)
      if (general) then
        call chain_general_s(as, j, keep)
      else
        call chain_block_s(as, j, keep, ok)
      end if
    else
      call load_block(ids, c, j - ns, a)
      if (general) then
        call chain_general(a, j - ns, keep)
      else
        call chain_block(a, j - ns, keep, ok)
      end if
    end if
  end subroutine run_block

  ! flag(i) = 0 for proven failures, 1 for everything that still needs the
  ! exact test.  pol is the row-major C buffer polBuf[i*11 + col].
  subroutine pdf_filter(pol, n, flag) bind(C, name=PDF_FILTER_NAME)
    integer(c_int), value :: n
    integer(i8), intent(in) :: pol(11, n)
    integer(c_signed_char), intent(out) :: flag(n)
    logical :: keep(BLK), ok(BLK), usegen, probe
    integer :: na, nn, nu, j, s, c, l, k
    integer(i8) :: x

    if (n <= 0) return
    if (allocated(idx)) then
      if (size(idx) < n) deallocate(idx, nidx, uidx, p0A, p1A, p2A)
    end if
    if (.not. allocated(idx)) allocate(idx(n), nidx(n), uidx(n), p0A(n, 0:9), p1A(n, 0:9), p2A(n, 0:9))

    do k = 0, 9
      do l = 1, n
        x = pol(11-k, l)
        p0A(l, k) = int(iand(x, 2_i8**21 - 1), i4)
        p1A(l, k) = int(iand(shifta(x, 21), 2_i8**21 - 1), i4)
        p2A(l, k) = int(shifta(x, 42), i4)
      end do
    end do
    do l = 1, n
      idx(l) = l
    end do
    na = n
    ! Stages 1..ns use the small primes, ns+1..ns+nq the large ones.
    ncall = ncall + 1
    probe = (mod(ncall, 64) == 1)
    do j = 1, ns + nq
      nn = 0
      nu = 0
      usegen = allgen .or. (prefgen(j) .and. .not. probe)
      ! Pass 1: every candidate; lanes that left the generic chain go to uidx.
      do s = 1, na, BLK
        c = min(BLK, na - s + 1)
        call run_block(idx(s:s+c-1), c, j, usegen, keep, ok)
        if (usegen) ok = .true.
        do l = 1, c
          if (.not. ok(l)) then
            nu = nu + 1
            uidx(nu) = idx(s + l - 1)
          else if (keep(l)) then
            nn = nn + 1
            nidx(nn) = idx(s + l - 1)
          end if
        end do
      end do
      if (.not. usegen) prefgen(j) = (2*nu > na)
      ! Pass 2: those lanes through the general chain.
      do s = 1, nu, BLK
        c = min(BLK, nu - s + 1)
        call run_block(uidx(s:s+c-1), c, j, .true., keep, ok)
        do l = 1, c
          if (keep(l)) then
            nn = nn + 1
            nidx(nn) = uidx(s + l - 1)
          end if
        end do
      end do
      idx(1:nn) = nidx(1:nn)
      if (debug) write(0, '(a,i3,a,i9,a,i9,a,i9)') 'pdf prime', j, ': in', na, ' retried', nu, ' kept', nn
      na = nn
      if (na == 0) exit
    end do
    debug = .false.

    flag = 0_c_signed_char
    do l = 1, na
      flag(idx(l)) = 1_c_signed_char
    end do
  end subroutine pdf_filter

end module PDF_MOD
