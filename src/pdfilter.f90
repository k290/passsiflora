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
!  * Small primes 2^8 < q < 2^13: 16-bit lanes (sixteen per AVX2 register),
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
!  * Loading: per prime and block, the coefficients of the block's polynomials
!    are gathered from pol (gather), and every coefficient c is split as
!    c = p0 + p1*2^21 + p2*2^42 (p0, p1 in [0, 2^21), p2 = floor(c/2^42) in
!    [-2^21, 2^21)); p0*d0 + p1*d1 + p2*d2 with d_i = 2^(21i)*2^32 mod q is
!    reduced once by sredc, which gives c mod q itself, in (-q, q) (ldcoef).
! The first 24 stages use the small primes; each prime rejects about half of
! the candidates, so almost only polynomials whose discriminant is a square
! (mod those primes) continue with the six large primes.  All inner loops
! vectorise across polynomials (BLK lanes).
!
! Contract (pdfilter_chains.tex, Section "The interface contract"):
!  * pdf_init(numP, pSet): pSet(1:numP) is the set S of primes that the exact
!    test divides out (upstream's testPolys allows at most 3).  It chooses the
!    filter primes.  Before the first call no prime is chosen: every flag is 1.
!  * pdf_filter(pol, n, flag): pol is the row-major C buffer polBuf[i*11 + col]
!    of n polynomials with int64 coefficients, col 0 the leading coefficient.
!    n <= 0: nothing is done; otherwise flag(1:n) is written.  flag = 0 only
!    for a monic polynomial whose discriminant is a non-square modulo a filter
!    prime; every other polynomial, among them every one with a leading
!    coefficient other than 1, gets flag = 1.
!  * Not thread-safe: the module keeps state between calls (the work arrays,
!    prefgen, ncall, the constants of pdf_init).  The application calls it from
!    one thread.
! The branch name barrett-reduction is historic: the arithmetic is Montgomery
! reduction.

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
#ifdef PDF_STAGES_NAME
  public :: pdf_stages
#endif

  integer, parameter :: i2 = c_int16_t, i4 = c_int32_t, i8 = c_int64_t
  ! At most MAXS small and MAXQ large primes, 30 stages in all.  A 16-bit stage has twice the lanes per register of a
  ! 32-bit one and cheaper reductions, so the small primes do the work and the large ones see few polynomials.
  integer, parameter :: MAXQ = 6, MAXS = 24, BLK = 64
  ! The limits the proofs need (pdfilter_chains.tex, Proposition "Every value is in range" and Theorem "The limits of the
  ! whole filter"): QMIN <= q <= QS_MAX for a small prime, QMIN <= q <= QL_MAX for a large one.  pdf_init skips any other
  ! prime (function chosen).  pdf_init searches the usable primes below QS_START down to QS_STOP (exclusive), then below
  ! QL_START down to QL_STOP; these ranges lie inside the limits.
  integer(i8), parameter :: QMIN = 11, QS_MAX = 8385, QL_MAX = 1097010265
  integer(i8), parameter :: QS_START = 2_i8**13 - 1, QS_STOP = 256
  integer(i8), parameter :: QL_START = 2_i8**25 - 1, QL_STOP = 2_i8**24

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
  logical, save :: debug = .false., allgen = .false.
  ! Per prime stage: go straight to chain_general when most lanes of the last
  ! probe left the generic chain (structurally degenerate families).
  logical, save :: prefgen(MAXS + MAXQ) = .false.
  integer, save :: ncall = 0
#ifdef PDF_STAGES_NAME
  ! Test aid, compiled only when the build defines PDF_STAGES_NAME (the test build scripts do, the
  ! application's do not): the stage at which each polynomial of the last call was rejected (0: kept),
  ! and per stage, summed over calls: lanes in, through chain_block in pass 1, through chain_general
  ! in pass 1, retried in pass 2, kept, and iterations of chain_general's lockstep and masked branches.
  integer(c_int), allocatable, save :: rstage(:)
  integer(c_long_long), save :: cnt(7, MAXS + MAXQ) = 0
#endif

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

  ! The guard of the proofs: a usable prime is chosen only if QMIN <= q <= qmax; otherwise it is skipped, with one line
  ! on stderr.  Test builds may define PDF_NO_LIMITS to switch the guard off (never the application's build).
  logical function chosen(q, qmax, numP, pSet)
    integer(i8), intent(in) :: q, qmax
    integer(c_int), intent(in) :: numP
    integer(c_int), intent(in) :: pSet(*)
    chosen = usable(q, numP, pSet)
#ifndef PDF_NO_LIMITS
    if (chosen .and. (q < QMIN .or. q > qmax)) then
      write(0, '(a,i0,a)') 'pdf_init: usable prime ', q, ' skipped: outside the proved limits of the integer arithmetic'
      chosen = .false.
    end if
#endif
  end function chosen

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
    ! A new initialisation starts a new routing history: the calls are counted from here.
    prefgen = .false.
    ncall = 0

    ns = 0
    q = QS_START
    do while (ns < MAXS .and. q > QS_STOP)
      if (chosen(q, QS_MAX, numP, pSet)) then
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
    q = QL_START
    do while (nq < MAXQ .and. q > QL_STOP)
      if (chosen(q, QL_MAX, numP, pSet)) then
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

    ! The start-up check: if the compiled arithmetic disagrees with the exact values, use no prime (every flag is 1).
    if (.not. arith_ok()) then
      write(0, '(a)') 'pdf_init: the integer arithmetic failed its start-up check; the pre-filter is off (every flag 1)'
      ns = 0
      nq = 0
    end if
  end subroutine pdf_init

  ! Signed Montgomery multiplication, R = 2^16, for the small primes and
  ! qinv = q^-1 mod 2^16 as a signed 16-bit value.  With x = a*b and
  ! m = (x*qinv) mod 2^16 in [-2^15, 2^15), x - m*q is divisible by 2^16, so
  ! (x - m*q)/2^16 = hi(x) - hi(m*q) exactly (hi = 16-bit multiply-high, the
  ! low halves of x and m*q are equal).  It is congruent to a*b*2^-16 and lies
  ! in (-q, q) when |x| < 2^15*q: for all |a|, |b| < 2q if q < 2^13, and for the
  ! values of the chains up to q = 8385 (pdfilter_chains.tex, Theorem "The
  ! limits of the whole filter").  The truncations to 16 bits are conversions
  ! (gfortran reduces them modulo 2^16; checked at start-up by arith_ok), not
  ! overflows.
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

  ! c mod q in (-q, q) for an int64 coefficient c: its pieces p0, p1, p2 (see the header) through ldres.
  elemental integer(i4) function ldcoef(c, d0, d1, d2, q, qinv)
    integer(i8), intent(in) :: c
    integer(i4), intent(in) :: d0, d1, d2, q, qinv
    ldcoef = ldres(int(iand(c, 2_i8**21 - 1), i4), int(iand(shifta(c, 21), 2_i8**21 - 1), i4), int(shifta(c, 42), i4), &
                   d0, d1, d2, q, qinv)
  end function ldcoef

  ! Test value k of a known-answer table with |value| <= bound: the extremes first, then pseudo-random values.
  integer(i8) function kat_val(k, bound, rng)
    integer, intent(in) :: k
    integer(i8), intent(in) :: bound
    integer(i8), intent(inout) :: rng
    select case (k)
    case (1); kat_val = bound
    case (2); kat_val = -bound
    case (3); kat_val = bound - 1
    case (4); kat_val = 1 - bound
    case (5); kat_val = 0
    case (6); kat_val = 1
    case (7); kat_val = -1
    case (8); kat_val = bound / 2
    case default
      rng = ieor(rng, shiftl(rng, 13)); rng = ieor(rng, shiftr(rng, 7)); rng = ieor(rng, shiftl(rng, 17))
      if (bound == huge(0_i8)) then
        kat_val = rng
      else
        kat_val = modulo(rng, 2*bound + 1) - bound
      end if
    end select
  end function kat_val

  ! Start-up check of the compiled arithmetic (pdfilter_chains.tex, Proposition "The start-up check").  For every chosen
  ! prime: mont (small primes), sredc and ldres (both, with the stored constants) and the 21-bit split, called
  ! elementally on arrays as in the chains, and the loads themselves (load_block_s or load_block: gather, split, ldres) on a
  ! block of the same coefficients, on extreme and pseudo-random inputs, against values computed with MODULO and exact
  ! division only (no conversion out of range, no bit operation on a negative value).  .false. on any mismatch.
  logical function arith_ok()
    integer(i2) :: a2(BLK), b2(BLK), u2(BLK), as(BLK, 0:10)
    integer(i4) :: u4(BLK), p0(BLK), p1(BLK), p2(BLK), q4, qi4, d(0:2), al(BLK, 0:10)
    integer(i8) :: x(BLK), c(BLK), q, lo, m, t, rng, pk(11, BLK), v
    integer :: j, l, k, stage, ids(BLK)
    arith_ok = .true.
    rng = 88172645463325252_i8
    do stage = 1, ns + nq
      if (stage <= ns) then
        j = stage
        q = qs(j)
        do l = 1, BLK
          a2(l) = int(kat_val(l, 2*q - 1, rng), i2)
          b2(l) = int(kat_val(BLK + 1 - l, 2*q - 1, rng), i2)
        end do
        u2 = mont(a2, b2, qs(j), qinvs(j))
        do l = 1, BLK
          x(l) = int(a2(l), i8)*int(b2(l), i8)
          lo = modulo(x(l) + 32768, 65536_i8) - 32768
          m = modulo(lo*int(qinvs(j), i8) + 32768, 65536_i8) - 32768
          if (int(u2(l), i8)*65536 /= x(l) - m*q) arith_ok = .false.
        end do
        if (modulo(q*int(qinvs(j), i8) - 1, 65536_i8) /= 0) arith_ok = .false.
        q4 = qsl(j); qi4 = qinvsl(j); d = lds(:, j)
      else
        j = stage - ns
        q = qv(j)
        q4 = qv(j); qi4 = qinvv(j); d = ldv(:, j)
      end if
      if (modulo(q*int(qi4, i8) - 1, 2_i8**32) /= 0) arith_ok = .false.
      do l = 1, BLK                                    ! sredc on |x| < 2^31 q
        x(l) = kat_val(l, 2_i8**31*q - 1, rng)
      end do
      u4 = sredc(x, q4, qi4)
      do l = 1, BLK
        lo = modulo(x(l) + 2_i8**31, 2_i8**32) - 2_i8**31
        m = modulo(lo*int(qi4, i8) + 2_i8**31, 2_i8**32) - 2_i8**31
        if (int(u4(l), i8)*2_i8**32 /= x(l) - m*q) arith_ok = .false.
      end do
      do l = 1, BLK                                    ! the split and ldres on every int64 coefficient class
        c(l) = kat_val(l, huge(0_i8), rng)
      end do
      c(1:4) = [-huge(0_i8) - 1, 2_i8**42, -2_i8**42, 2_i8**21 - 1]
      p0 = int(iand(c, 2_i8**21 - 1), i4)
      p1 = int(iand(shifta(c, 21), 2_i8**21 - 1), i4)
      p2 = int(shifta(c, 42), i4)
      u4 = ldres(p0, p1, p2, d(0), d(1), d(2), q4, qi4)
      do l = 1, BLK
        t = (c(l) - modulo(c(l), 2_i8**21)) / 2_i8**21
        if (p0(l) /= modulo(c(l), 2_i8**21) .or. p1(l) /= modulo(t, 2_i8**21) .or. &
            p2(l) /= (t - modulo(t, 2_i8**21)) / 2_i8**21) arith_ok = .false.
        if (abs(int(u4(l), i8)) >= q .or. modulo(int(u4(l), i8), q) /= modulo(c(l), q)) arith_ok = .false.
      end do
      ! The loads: polynomial l has the coefficients c rotated by l - 1 + 7k, in pol's layout; the block lists them in
      ! reverse order.
      do l = 1, BLK
        pk(1, l) = 1
        do k = 0, 9
          pk(11-k, l) = c(1 + mod(l - 1 + 7*k, BLK))
        end do
        ids(l) = BLK + 1 - l
      end do
      if (stage <= ns) then
        call load_block_s(pk, ids, BLK, j, as)
        al = int(as, i4)
      else
        call load_block(pk, ids, BLK, j, al)
      end if
      do l = 1, BLK
        do k = 0, 9
          v = int(al(l, k), i8)
          if (abs(v) >= q .or. modulo(v, q) /= modulo(pk(11-k, ids(l)), q)) arith_ok = .false.
        end do
        if (al(l, 10) /= 1) arith_ok = .false.
      end do
    end do
  end function arith_ok

  ! The coefficients of a block: x(l, k) = coefficient of x^k of polynomial ids(l), l = 1..c, k = 0..9 (pol(11-k, ids(l)):
  ! each polynomial is 11 consecutive values of pol); the lanes c+1..BLK are 0.  Gathering whole polynomials first keeps the
  ! arithmetic of the loads on contiguous lanes, where gfortran vectorises it.
  subroutine gather(pol, ids, c, x)
    integer(i8), intent(in) :: pol(11, *)
    integer, intent(in) :: c
    integer, intent(in) :: ids(c)
    integer(i8), intent(out) :: x(BLK, 0:9)
    integer :: l, k
    do l = 1, c
      do k = 0, 9
        x(l, k) = pol(11-k, ids(l))
      end do
    end do
    x(c+1:BLK, :) = 0
  end subroutine gather

  ! Small-prime versions of load_block / chain_block (same algorithm).
  subroutine load_block_s(pol, ids, c, j, a)
    integer(i8), intent(in) :: pol(11, *)
    integer, intent(in) :: c, j
    integer, intent(in) :: ids(c)
    integer(i2), intent(out) :: a(BLK, 0:10)
    integer(i4) :: q, qinv, d0, d1, d2
    integer(i8) :: x(BLK, 0:9)
    integer :: k

    q = qsl(j); qinv = qinvsl(j); d0 = lds(0, j); d1 = lds(1, j); d2 = lds(2, j)
    call gather(pol, ids, c, x)
    do k = 0, 9
      a(:, k) = int(ldcoef(x(:, k), d0, d1, d2, q, qinv), i2)
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
  ! a(l,k) = coefficient of x^k of polynomial ids(l), k = 0..10, in (-q, q);
  ! the lanes c+1..BLK are 0 (and 1 at k = 10).
  subroutine load_block(pol, ids, c, j, a)
    integer(i8), intent(in) :: pol(11, *)
    integer, intent(in) :: c, j
    integer, intent(in) :: ids(c)
    integer(i4), intent(out) :: a(BLK, 0:10)
    integer(i4) :: q, qinv, d0, d1, d2
    integer(i8) :: x(BLK, 0:9)
    integer :: k

    q = qv(j); qinv = qinvv(j); d0 = ldv(0, j); d1 = ldv(1, j); d2 = ldv(2, j)
    call gather(pol, ids, c, x)
    do k = 0, 9
      a(:, k) = ldcoef(x(:, k), d0, d1, d2, q, qinv)
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
#define GSTAGE j
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
#undef GSTAGE
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
#define GSTAGE (ns + j)
#define GMUL(a, b) sredc(int(a, i8)*int(b, i8), q, qinv)
#define GDIFF(a, b, c, d) sredc(int(a, i8)*int(b, i8) - int(c, i8)*int(d, i8), q, qinv)
#define GZERO(x) (x == 0)
#include "chain_general.inc"

  ! One block through prime stage j (small primes first, then large ones),
  ! with the fast generic chain or, for lanes that left it, the general one.
  subroutine run_block(pol, ids, c, j, general, keep, ok)
    integer(i8), intent(in) :: pol(11, *)
    integer, intent(in) :: c, j
    integer, intent(in) :: ids(c)
    logical, intent(in) :: general
    logical, intent(out) :: keep(BLK), ok(BLK)
    integer(i4) :: a(BLK, 0:10)
    integer(i2) :: as(BLK, 0:10)
    if (j <= ns) then
      call load_block_s(pol, ids, c, j, as)
      if (general) then
        call chain_general_s(as, j, keep)
      else
        call chain_block_s(as, j, keep, ok)
      end if
    else
      call load_block(pol, ids, c, j - ns, a)
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
    integer :: na, nn, nu, j, s, c, l

    if (n <= 0) return
    if (allocated(idx)) then
      if (size(idx) < n) deallocate(idx, nidx, uidx)
    end if
    if (.not. allocated(idx)) allocate(idx(n), nidx(n), uidx(n))

    do l = 1, n
      idx(l) = l
    end do
#ifdef PDF_STAGES_NAME
    if (allocated(rstage)) then
      if (size(rstage) < n) deallocate(rstage)
    end if
    if (.not. allocated(rstage)) allocate(rstage(n))
    rstage(1:n) = 0
#endif
    na = n
    ! Stages 1..ns use the small primes, ns+1..ns+nq the large ones.  Calls are counted modulo 64 (1, ..., 64, 1, ...),
    ! so the counter cannot overflow; the probe is calls 1, 65, 129, ...
    ncall = mod(ncall, 64) + 1
    probe = (ncall == 1)
    do j = 1, ns + nq
      nn = 0
      nu = 0
      usegen = allgen .or. (prefgen(j) .and. .not. probe)
      ! Pass 1: every candidate; lanes that left the generic chain go to uidx.
      do s = 1, na, BLK
        c = min(BLK, na - s + 1)
        call run_block(pol, idx(s:s+c-1), c, j, usegen, keep, ok)
        if (usegen) ok = .true.
        do l = 1, c
          if (.not. ok(l)) then
            nu = nu + 1
            uidx(nu) = idx(s + l - 1)
          else if (keep(l)) then
            nn = nn + 1
            nidx(nn) = idx(s + l - 1)
#ifdef PDF_STAGES_NAME
          else
            rstage(idx(s + l - 1)) = j
#endif
          end if
        end do
      end do
      if (.not. usegen) prefgen(j) = (2*nu > na)
      ! Pass 2: those lanes through the general chain.
      do s = 1, nu, BLK
        c = min(BLK, nu - s + 1)
        call run_block(pol, uidx(s:s+c-1), c, j, .true., keep, ok)
        do l = 1, c
          if (keep(l)) then
            nn = nn + 1
            nidx(nn) = uidx(s + l - 1)
#ifdef PDF_STAGES_NAME
          else
            rstage(uidx(s + l - 1)) = j
#endif
          end if
        end do
      end do
#ifdef PDF_STAGES_NAME
      cnt(1, j) = cnt(1, j) + na
      if (usegen) then
        cnt(3, j) = cnt(3, j) + na
      else
        cnt(2, j) = cnt(2, j) + na
      end if
      cnt(4, j) = cnt(4, j) + nu
      cnt(5, j) = cnt(5, j) + nn
#endif
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
    ! The chains treat every polynomial as monic (a(:, 10) = 1; pol(1, l) is never read): any other one is kept.
    do l = 1, n
      if (pol(1, l) /= 1) then
        flag(l) = 1_c_signed_char
#ifdef PDF_STAGES_NAME
        rstage(l) = 0
#endif
      end if
    end do
  end subroutine pdf_filter

#ifdef PDF_STAGES_NAME
  ! Test aid: st(i) = stage at which polynomial i of the last pdf_filter call was rejected (0: kept;
  ! stages 1..ns are the small primes, ns+1..ns+nq the large ones, in the order pdf_init chose them);
  ! counts(:, j) = the counters of stage j (see cnt above), zeroed afterwards when reset /= 0;
  ! primes(j) = the prime of stage j; nstages = ns + nq.  n = 0 reads only the counters and primes.
  subroutine pdf_stages(n, st, counts, primes, nstages, reset) bind(C, name=PDF_STAGES_NAME)
    integer(c_int), value :: n, reset
    integer(c_int), intent(out) :: st(n), nstages
    integer(c_long_long), intent(out) :: counts(7, MAXS + MAXQ), primes(MAXS + MAXQ)
    if (n > 0) st = rstage(1:n)
    counts = cnt
    primes = 0
    primes(1:ns) = qs(1:ns)
    primes(ns+1:ns+nq) = qv(1:nq)
    nstages = ns + nq
    if (reset /= 0) cnt = 0
  end subroutine pdf_stages
#endif

end module PDF_MOD
