! arith_checks.f90 -- checks the compiled integer arithmetic of the pre-filter against exact models.
! Built and run by arith_checks.sh, which extracts the text of inv32, signedmod, mont, sredc and ldres (and the line that
! defines the kinds) from ../src/pdfilter.f90 into $OUT/arith_kinds.inc and $OUT/arith_funcs.inc; they are included below,
! in this module, so that gfortran inlines and vectorises them as it does in pdfilter.f90.  The calls are elemental calls
! on arrays, as in the filter.
!
! Usage: arith_checks PART [arguments]                                    (each part prints one result line per check)
!   mont-pairs Q0 Q1 SEL   every pair |a|,|b| <= 2q-1 for every q in [Q0,Q1] (SEL = odd: all odd q; prime: primes only):
!                          the result equals the definition (x - m q)/2^16, m = s16(x q'), and has the residue a b 2^-16;
!                          counts the pairs with |u| >= q (Lemma mont: none for q < 2^13)
!   mont-scan Q0 Q1        for every odd q: the smallest-x witness of |u| >= q from the characterisation, each side,
!                          checked with the compiled mont; q without a witness are checked by brute force
!   mont-extremes Q0 Q1    for every odd q: the largest and smallest result over |a|,|b| <= 2q-1, and the first q at which
!                          their difference leaves 16 bits
!   sredc                 windows around x = +-2^31 q, the counterexamples, the extreme classes, random x, differences
!                          a b - c d with a,b,c,d at and near +-(q-1); model in 128-bit integers
!   ldres                  the loading on extreme and random coefficients for 2,585 moduli (with 200 random odd q < 2^31)
!   inv32                  q q^-1 = 1 mod 2^32 and the signed 16/32-bit forms, for many q
!   bounds-small Q0 Q1     the bound propagation of the small-prime stage (Section 8.8) for every odd q in [Q0,Q1]
!   bounds-large Q0 Q1     the same for the large-prime stage
! The exact models use integer(16) (128 bits) where a value can exceed 64 bits, and integer(8) where it provably cannot.

module arith_checks_mod
  use, intrinsic :: iso_c_binding
  implicit none
#include "arith_kinds.inc"
  integer, parameter :: i16 = 16
contains
#include "arith_funcs.inc"

  pure integer(i8) function s16v(x)          ! x mod 2^16 as a signed 16-bit value
    integer(i8), intent(in) :: x
    s16v = modulo(x, 65536_i8)
    if (s16v >= 32768) s16v = s16v - 65536
  end function s16v

  pure integer(i16) function s32v(x)          ! x mod 2^32 as a signed 32-bit value
    integer(i16), intent(in) :: x
    s32v = modulo(x, 2_i16**32)
    if (s32v >= 2_i16**31) s32v = s32v - 2_i16**32
  end function s32v

  pure logical function isprime(n)
    integer(i8), intent(in) :: n
    integer(i8) :: d
    isprime = .false.
    if (n < 2) return
    if (mod(n, 2_i8) == 0) then; isprime = (n == 2); return; end if
    d = 3
    do while (d*d <= n)
      if (mod(n, d) == 0) return
      d = d + 2
    end do
    isprime = .true.
  end function isprime

  pure integer(i8) function powm(b, e, m)
    integer(i8), intent(in) :: b, e, m
    integer(i16) :: r, x
    integer(i8) :: ee
    r = 1; x = modulo(int(b, i16), int(m, i16)); ee = e
    do while (ee > 0)
      if (iand(ee, 1_i8) == 1) r = modulo(r*x, int(m, i16))
      x = modulo(x*x, int(m, i16)); ee = shiftr(ee, 1)
    end do
    powm = int(r, i8)
  end function powm

  pure logical function usable25(q)            ! pdf_init's rule for S = {2, 5}
    integer(i8), intent(in) :: q
    usable25 = mod(q, 4_i8) == 1 .and. isprime(q)
    if (usable25) usable25 = powm(2_i8, (q-1)/2, q) == 1 .and. powm(5_i8, (q-1)/2, q) == 1
  end function usable25

  ! the 16-bit constants of a small prime, as pdf_init computes them
  subroutine small_consts(q, q2, qinv2)
    integer(i8), intent(in) :: q
    integer(i2), intent(out) :: q2, qinv2
    q2 = int(q, i2)
    qinv2 = int(signedmod(inv32(q), 16), i2)
  end subroutine small_consts

  ! every pair |a|,|b| <= 2q-1: returns the number of model mismatches, residue errors and pairs with |u| >= q
  subroutine mont_pairs(q, nmis, nres, nout, umax, umin)
    integer(i8), intent(in) :: q
    integer(i8), intent(out) :: nmis, nres, nout, umax, umin
    integer(i2), allocatable :: bv(:), uv(:)
    integer(i2) :: q2, qinv2, a2
    integer(i8) :: a, n, j, x, mm, um
    call small_consts(q, q2, qinv2)
    n = 4*q - 1
    allocate(bv(n), uv(n))
    do j = 1, n
      bv(j) = int(j - 2*q, i2)
    end do
    nmis = 0; nres = 0; nout = 0; umax = -huge(0_i8); umin = huge(0_i8)
    do a = -(2*q - 1), 2*q - 1
      a2 = int(a, i2)
      uv = mont(a2, bv, q2, qinv2)
      nout = nout + count(abs(int(uv, i8)) >= q)
      umax = max(umax, int(maxval(uv), i8)); umin = min(umin, int(minval(uv), i8))
      ! the model, exact in 64 bits (|x| <= 2^30, |x q'| <= 2^45): m = s16(x q'), and u 2^16 = x - m q must hold
      ! exactly; then u = a b 2^-16 mod q follows.  nres counts x - m q not divisible by 2^16 (impossible by Lemma 8.1).
      do j = 1, n
        x = a*int(bv(j), i8)
        mm = iand(x*int(qinv2, i8), 65535_i8)
        mm = mm - 2*iand(mm, 32768_i8)
        um = x - mm*q
        if (iand(um, 65535_i8) /= 0) nres = nres + 1
        if (um /= int(uv(j), i8)*65536) nmis = nmis + 1
      end do
    end do
  end subroutine mont_pairs

  ! characterisation search (Lemma 8.2 of the text): smallest x >= (R/2) q with x = (R+m) q + R t (positive side), or
  ! |x| = (R-m) q + R t (negative side), that is a product of two numbers <= 2q-1.  Returns 0 if there is none.
  subroutine mont_witness(q, neg, a, b)
    integer(i8), intent(in) :: q
    logical, intent(in) :: neg
    integer(i8), intent(out) :: a, b
    integer(i8), parameter :: R = 65536
    integer(i8) :: m, x, xx, lo, bnd, m0, m1, ms
    bnd = 2*q - 1
    a = 0; b = 0
    if (.not. neg) then; m0 = -R/2; m1 = R/2 - 1; ms = 1
    else; m0 = R/2 - 1; m1 = -R/2; ms = -1; end if
    do m = m0, m1, ms
      if (.not. neg) then; x = (R + m)*q; else; x = (R - m)*q; end if
      if (x > bnd*bnd) exit
      do xx = x, bnd*bnd, R
        lo = (xx + bnd - 1) / bnd
        do a = lo, bnd
          if (a*a > xx) exit
          if (mod(xx, a) == 0 .and. xx / a <= bnd) then
            b = xx / a
            if (neg) b = -b
            return
          end if
        end do
      end do
    end do
    a = 0; b = 0
  end subroutine mont_witness

  ! largest and smallest result of the compiled mont over |a|,|b| <= 2q-1.  mont depends on a b only, so a > 0 suffices;
  ! b > 0 gives the maximum and b < 0 the minimum (a result of the other sign has |u| <= q/2, Lemma 8.1).  For each a the
  ! b with |a b| too small to beat the current extremes (|u| <= |a b|/2^16 + q/2) are skipped; the rest is one array call.
  subroutine mont_extremes(q, umax, umin)
    integer(i8), intent(in) :: q
    integer(i8), intent(out) :: umax, umin
    integer(i2), allocatable :: bv(:), uv(:)
    integer(i2) :: q2, qinv2
    integer(i8) :: a, j, bnd, need, nb
    call small_consts(q, q2, qinv2)
    bnd = 2*q - 1
    allocate(bv(bnd), uv(bnd))
    do j = 1, bnd
      bv(j) = int(bnd + 1 - j, i2)                    ! descending
    end do
    umax = 0; umin = 0
    do a = bnd, 1, -1
      need = 65536*min(umax, -umin) - q*32768          ! |a b| must reach this to give a new extreme
      if (a*bnd < need) exit
      nb = bnd - max(1_i8, (need + a - 1)/a) + 1       ! the b >= need/a
      uv(1:nb) = mont(int(a, i2), bv(1:nb), q2, qinv2)
      umax = max(umax, int(maxval(uv(1:nb)), i8))
      uv(1:nb) = mont(int(a, i2), -bv(1:nb), q2, qinv2)
      umin = min(umin, int(minval(uv(1:nb)), i8))
    end do
  end subroutine mont_extremes

  ! sredc on a vector, compared with the model (128 bits): m = s32(x q'), u = (x - m q)/2^32.  Counts: model mismatches
  ! (where the model's x - m q and u fit 64 and 32 bits), cases where they do not fit (the compiled code then wraps),
  ! results outside (-q, q) with |x| < 2^31 q (Lemma 8.3 says none), and results outside (-q, q) in all.
  subroutine sredc_vec(xv, q, nmis, novf, nout_in, nout_all)
    integer(i8), intent(in) :: xv(:), q
    integer(i8), intent(inout) :: nmis, novf, nout_in, nout_all
    integer(i4) :: uv(size(xv)), q4, qinv4
    integer(i16) :: m, um, y
    integer :: j
    q4 = int(q, i4); qinv4 = int(signedmod(inv32(q), 32), i4)
    uv = sredc(xv, q4, qinv4)
    do j = 1, size(xv)
      m = s32v(int(xv(j), i16)*int(qinv4, i16))
      y = int(xv(j), i16) - m*int(q, i16)
      um = y / 2_i16**32
      if (um*2_i16**32 /= y) nmis = nmis + 1                      ! cannot happen (Lemma 8.1)
      if (abs(y) >= 2_i16**63 .or. abs(um) >= 2_i16**31) then
        novf = novf + 1
      else if (um /= int(uv(j), i16)) then
        nmis = nmis + 1
      end if
      if (abs(int(uv(j), i8)) >= q) then
        nout_all = nout_all + 1
        if (abs(int(xv(j), i16)) < 2_i16**31*int(q, i16)) nout_in = nout_in + 1
      end if
    end do
  end subroutine sredc_vec

  ! bound propagation (Section 8.8): returns 0 if every reduction of the four chains satisfies its lemma for prime size q
  ! and every 16-bit value fits; otherwise a code naming the first call site that may fail
  integer function bounds_small(q) result(site)
    integer(i8), intent(in) :: q
    integer(i16), parameter :: R = 65536
    integer(i16) :: A(0:10), B(0:10), Rr(0:10), P(0:10), L, lc, top, t, bu, bv, bw, mm, qq
    integer :: k, m, it
    qq = q; site = 0
    L = (2*(2_i16**21 - 1)*(qq - 1) + 2_i16**21*(qq - 1) + qq*2_i16**31) / 2_i16**32
    A(0:9) = L; A(10) = 1
    B = 0
    do k = 1, 10
      B(k-1) = mo(A(k), int(k, i16), 1); if (site /= 0) return
    end do
    do m = 9, 1, -1
      lc = B(m); top = A(m+1)
      Rr = 0; P = 0
      Rr(0) = mo(lc, A(0), 10 + m); if (site /= 0) return
      do k = 1, m
        Rr(k) = df(mo(lc, A(k), 10 + m), mo(top, B(k-1), 10 + m), 10 + m); if (site /= 0) return
      end do
      t = Rr(m)
      do k = 0, m - 1
        P(k) = df(mo(lc, Rr(k), 20 + m), mo(t, B(k), 20 + m), 20 + m); if (site /= 0) return
      end do
      A = B; B = P
    end do
    t = mo(B(0), B(0), 30); t = mo(qq - 1, B(0), 30); t = mo(qq - 1, qq - 1, 30); if (site /= 0) return
    bu = max(L, 1_i16); bv = mo(10_i16, L, 40); mm = 1
    do it = 0, 17
      bw = df(mo(bv, bu, 41 + it), mo(bu, bv, 41 + it), 41 + it); if (site /= 0) return
      mm = mo(max(mm, qq - 1), max(bu, bv, bw), 41 + it); if (site /= 0) return
      bu = max(bu, bv, bw); bv = bu
    end do
    t = mo(qq - 1, bu, 60); t = mo(qq - 1, qq - 1, 60)
  contains
    integer(i16) function mo(x, y, where)
      integer(i16), intent(in) :: x, y
      integer, intent(in) :: where
      if (x*y >= qq*(R/2) .and. site == 0) site = where
      mo = (x*y + qq*(R/2)) / R
    end function mo
    integer(i16) function df(x, y, where)
      integer(i16), intent(in) :: x, y
      integer, intent(in) :: where
      if (x + y > 32767 .and. site == 0) site = 100 + where
      df = x + y
    end function df
  end function bounds_small

  integer function bounds_large(q) result(site)
    integer(i8), intent(in) :: q
    integer(i16), parameter :: R = 2_i16**32
    integer(i16) :: A(0:10), B(0:10), Rr(0:10), P(0:10), L, lc, top, t, bu, bv, bw, qq
    integer :: k, m, it
    qq = q; site = 0
    L = (2*(2_i16**21 - 1)*(qq - 1) + 2_i16**21*(qq - 1) + qq*2_i16**31) / 2_i16**32
    A(0:9) = L; A(10) = 1
    B = 0
    do k = 1, 10
      B(k-1) = sr(int(k, i16)*A(k), 1)
    end do
    do m = 9, 1, -1
      lc = B(m); top = A(m+1)
      Rr = 0; P = 0
      Rr(0) = sr(lc*A(0), 10 + m)
      do k = 1, m
        Rr(k) = sr(lc*A(k) + top*B(k-1), 10 + m)
      end do
      t = Rr(m)
      do k = 0, m - 1
        P(k) = sr(lc*Rr(k) + t*B(k), 20 + m)
      end do
      A = B; B = P
    end do
    t = sr(B(0)*B(0), 30); t = sr((qq - 1)*(qq - 1), 30)
    bu = max(L, 1_i16); bv = sr(10*L, 40)
    do it = 0, 17
      bw = sr(2*bu*bv, 41 + it)
      t = sr((qq - 1)*max(bu, bv, bw), 41 + it)
      bu = max(bu, bv, bw); bv = bu
    end do
  contains
    integer(i16) function sr(x, where)
      integer(i16), intent(in) :: x
      integer, intent(in) :: where
      if (x >= qq*(R/2) .and. site == 0) site = where
      sr = (x + qq*(R/2)) / R
    end function sr
  end function bounds_large
end module arith_checks_mod

program arith_checks
  use arith_checks_mod
  implicit none
  character(len=32) :: part, s1, s2, s3
  integer(i8) :: q0, q1, q, nq, nmis, nres, nout, umax, umin, tmis, tres, tout, a, b, a2, b2, nfailq, nlemma
  integer(i8) :: firstodd, firstprime, first14, firstuse, nclean, nposw, nnegw, npairs
  integer(i2) :: q2, qinv2, u2
  integer :: site, nargs
  logical :: primes_only

  nargs = command_argument_count()
  call get_command_argument(1, part)
  s1 = ''; s2 = ''; s3 = ''
  if (nargs >= 2) call get_command_argument(2, s1)
  if (nargs >= 3) call get_command_argument(3, s2)
  if (nargs >= 4) call get_command_argument(4, s3)

  select case (trim(part))
  case ('mont-pairs')
    read(s1, *) q0; read(s2, *) q1; primes_only = trim(s3) == 'prime'
    nq = 0; tmis = 0; tres = 0; tout = 0; nfailq = 0; nlemma = 0; npairs = 0
    do q = ior(q0, 1_i8), q1, 2
      if (primes_only .and. .not. isprime(q)) cycle
      call mont_pairs(q, nmis, nres, nout, umax, umin)
      nq = nq + 1; tmis = tmis + nmis; tres = tres + nres; tout = tout + nout; npairs = npairs + (4*q - 1)**2
      if (nout > 0) nfailq = nfailq + 1
      if (q < 8192 .and. nout > 0) nlemma = nlemma + 1
      if (nout > 0 .or. q > 8192) write(*, '(a,i0,a,i0,a,i0,a,i0)') 'mont q=', q, ': pairs outside (-q,q) ', nout, &
        ', max ', umax, ', min ', umin
    end do
    write(*, '(a,i0,a,i0,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') 'mont-pairs [', q0, ',', q1, '] ', trim(s3), ': q tested ', nq, &
      ', pairs ', npairs, &
      ', model mismatches ', tmis, ', residue errors ', tres, ', q with |u| >= q ', nfailq, ', of them below 2^13 ', nlemma
    if (tmis + tres + nlemma > 0) error stop 1

  case ('mont-count')
    read(s1, *) q
    call mont_pairs(q, nmis, nres, nout, umax, umin)
    write(*, '(a,i0,a,i0,a,i0,a,i0,a,f6.3,a,i0,a,i0)') 'mont-count q=', q, ': pairs |a|,|b| < 2q: ', (4*q-1)**2, &
      ', outside (-q,q): ', nout, ', max |u| = ', max(umax, -umin), ' (', real(max(umax, -umin), 8)/q, ' q), model mismatches ', &
      nmis, ', residue errors ', nres
    if (nmis + nres > 0) error stop 1

  case ('mont-scan')
    read(s1, *) q0; read(s2, *) q1
    firstodd = 0; firstprime = 0; first14 = 0; firstuse = 0; nclean = 0; nposw = 0; nnegw = 0; tmis = 0
    do q = ior(q0, 1_i8), q1, 2
      call small_consts(q, q2, qinv2)
      call mont_witness(q, .false., a, b)
      call mont_witness(q, .true., a2, b2)
      if (a /= 0) then
        u2 = mont(int(a, i2), int(b, i2), q2, qinv2)
        if (u2 < q) tmis = tmis + 1            ! the witness must fail in the compiled code
        nposw = nposw + 1
      end if
      if (a2 /= 0) then
        u2 = mont(int(a2, i2), int(b2, i2), q2, qinv2)
        if (u2 > -q) tmis = tmis + 1
        nnegw = nnegw + 1
      end if
      if (a == 0 .and. a2 == 0) then          ! no witness: check every pair with the compiled mont
        call mont_pairs(q, nmis, nres, nout, umax, umin)
        tmis = tmis + nmis + nres + nout
        nclean = nclean + 1
        write(*, '(a,i0,a,l1)') 'mont-scan: no counterexample for q = ', q, ' (all pairs checked), prime: ', isprime(q)
      else
        if (firstodd == 0) then
          firstodd = q
          write(*, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') 'mont-scan: first failing odd q = ', q, ': a = ', a, ', b = ', b, &
            ', mont = ', mont(int(a, i2), int(b, i2), q2, qinv2), '; negative side a = ', a2, ', b = ', b2
        end if
        if (firstprime == 0 .and. isprime(q)) firstprime = q
        if (first14 == 0 .and. isprime(q) .and. mod(q, 4_i8) == 1) first14 = q
        if (firstuse == 0 .and. usable25(q)) firstuse = q
      end if
    end do
    write(*, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') 'mont-scan [', q0, ',', q1, ']: witnesses positive ', nposw, &
      ', negative ', nnegw, ', no counterexample ', nclean, '; first failing odd q ', firstodd, ', prime ', firstprime, &
      ', prime = 1 mod 4 ', first14, ', usable for S={2,5} ', firstuse, '; disagreements ', tmis
    if (tmis > 0) error stop 1

  case ('mont-extremes')
    ! a difference of two results of mont leaves 16 bits when max - min >= 2^15.  At the two q next to the first such q
    ! the pruned search is compared with mont_pairs (every pair).
    read(s1, *) q0; read(s2, *) q1
    firstodd = 0; firstprime = 0; nfailq = 0; tmis = 0
    do q = ior(q0, 1_i8), q1, 2
      call mont_extremes(q, umax, umin)
      if (umax - umin >= 32768) then
        if (firstodd == 0) firstodd = q
        if (firstprime == 0 .and. isprime(q)) firstprime = q
      else if (firstodd /= 0) then
        nfailq = nfailq + 1                    ! below 2^15 again after the first q above
      end if
      if (q == 12799 .or. q == 12801) then
        call mont_pairs(q, nmis, nres, nout, a, b)
        if (a /= umax .or. b /= umin) tmis = tmis + 1
        write(*, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') 'mont-extremes q = ', q, ': max ', umax, ', min ', umin, &
          ', max - min ', umax - umin, '; every pair: max ', a, ', min ', b
      end if
    end do
    write(*, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') 'mont-extremes [', q0, ',', q1, ']: first odd q with max - min >= 2^15 ', &
      firstodd, ', first prime ', firstprime, ', odd q above it with max - min < 2^15 ', nfailq, '; disagreements ', tmis
    if (tmis > 0) error stop 1

  case ('sredc')
    call sredc_checks()

  case ('ldres')
    call ldres_checks()

  case ('inv32')
    call inv32_checks()

  case ('bounds-small', 'bounds-large')
    read(s1, *) q0; read(s2, *) q1
    firstodd = 0; nq = 0; nfailq = 0
    do q = ior(q0, 1_i8), q1, 2
      if (trim(part) == 'bounds-small') then; site = bounds_small(q); else; site = bounds_large(q); end if
      if (site == 0) then
        nq = nq + 1
        if (firstodd /= 0) nfailq = nfailq + 1              ! passes again above a failure
      else if (firstodd == 0) then
        firstodd = q
        write(*, '(a,a,i0,a,i0)') trim(part), ': first odd q that fails the propagation: ', q, ', call site code ', site
      end if
    end do
    write(*, '(a,a,i0,a,i0,a,i0,a,i0,a,i0)') trim(part), ' [', q0, ',', q1, ']: odd q passing ', nq, ', first failing ', firstodd, &
      ', passing again above it ', nfailq

  case default
    write(*, '(a)') 'unknown part'
    error stop 2
  end select

contains

  subroutine sredc_checks()
    integer(i8), allocatable :: xv(:)
    integer(i8) :: qs(9), x0, wnd, j, k, mis, ovf, oin, oall, nwin, rng, qr, lim, nqd
    integer(i8) :: ca, cb, cc, cd, dd(5)
    integer(i4) :: q4, qinv4, u
    integer(i16) :: xw
    integer :: iq, ia, ib, ic, id
    ! moduli: below 2^30, 2^30 + 1 (the end of the sufficient bound), the smallest failing odd q, prime and usable prime,
    ! the family x = 2^31 q, the end of the filter-level bound, the first q where an i8 overflow is possible, 2^31 - 1
    qs = [1073741789_i8, 1073741825_i8, 1073765103_i8, 1073765117_i8, 1073765201_i8, 1073774593_i8, 1097010265_i8, &
          1676704567_i8, 2147483647_i8]
    ! (a) windows of 2*2^16 consecutive x around +2^31 q and -2^31 q
    wnd = 65536
    allocate(xv(2*wnd))
    mis = 0; ovf = 0; oin = 0; oall = 0; nwin = 0
    do iq = 1, size(qs)
      q = qs(iq)
      do k = -1, 1, 2
        x0 = k*(2_i8**31)*q
        do j = 1, 2*wnd
          xv(j) = x0 - wnd + j - 1
        end do
        call sredc_vec(xv, q, mis, ovf, oin, oall)
        nwin = nwin + 2*wnd
      end do
    end do
    write(*, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'sredc windows around +-2^31 q (9 moduli): ', nwin, ' x, model mismatches ', mis, &
      ', not representable ', ovf, ', |u| >= q with |x| < 2^31 q ', oin, ', |u| >= q in all ', oall
    if (mis + oin > 0) error stop 1
    ! (b) the counterexamples: x = ab - cd formed as in chain_block and GDIFF
    call replay(1073765103_i8, 1073765102_i8, 1073765101_i8, -1073764896_i8, 1073764888_i8)    ! smallest odd q
    call replay(1073765117_i8, 1073765111_i8, 1073765111_i8, -1073764914_i8, 1073764866_i8)    ! smallest prime
    call replay(1073765201_i8, 1073765200_i8, 1073765198_i8, -1073765108_i8, 1073765063_i8)    ! smallest usable prime
    call replay(1073774593_i8, 1073741825_i8, 1073741824_i8, -1073774592_i8, 1073774592_i8)    ! x = 2^31 q exactly
    ! (c) the extreme classes of Lemma 8.2: x = (2^31 + mu) q + 2^32 t gives u = q + t, x = -(2^31 + 1 + mu) q - 2^32 t
    !     gives u = -(q + t); just below them (x = 2^31 q - 1 and -(2^31 + 1) q + 1) |u| < q.  Where q + t does not fit 32 bits
    !     (q = 2^31 - 1) the conversion to i4 wraps: the result is q + t - 2^32 (counted as "wrapped")
    mis = 0; ovf = 0
    do iq = 1, size(qs)
      q = qs(iq)
      q4 = int(q, i4); qinv4 = int(signedmod(inv32(q), 32), i4)
      do j = 0, 3
        do k = 0, 3
          xw = (2_i16**31 + j)*int(q, i16) + 2_i16**32*k
          u = sredc(int(xw, i8), q4, qinv4)
          if (q + k <= 2_i8**31 - 1) then
            if (int(u, i8) /= q + k) mis = mis + 1
          else
            ovf = ovf + 1
            if (int(u, i8) /= q + k - 2_i8**32) mis = mis + 1
          end if
          u = sredc(-int(xw, i8) - q, q4, qinv4)
          if (q + k <= 2_i8**31) then
            if (int(u, i8) /= -(q + k)) mis = mis + 1
          else
            ovf = ovf + 1
            if (int(u, i8) /= -(q + k) + 2_i8**32) mis = mis + 1
          end if
        end do
      end do
      u = sredc((2_i8**31)*q - 1, q4, qinv4)
      if (abs(int(u, i8)) >= q) mis = mis + 1
      u = sredc(-(2_i8**31 + 1)*q + 1, q4, qinv4)
      if (abs(int(u, i8)) >= q) mis = mis + 1
    end do
    write(*, '(a,i0,a,i0)') 'sredc extreme classes x = +-((2^31 + mu) q + 2^32 t), mu,t = 0..3, and the x just below, ' // &
      '9 moduli: wrong ', mis, &
      ', results that wrap to 32 bits as predicted ', ovf
    if (mis > 0) error stop 1
    ! (d) random x with |x| < 2^31 q (Lemma 8.3: |u| < q), and with 2^31 q <= |x| <= 2 (q-1)^2, for 2000 random odd q < 2^31
    mis = 0; ovf = 0; oin = 0; oall = 0; rng = 88172645463325252_i8; nwin = 0
    do iq = 1, 2000
      qr = ior(modulo(xs(rng), 2_i8**31 - 3) + 3, 1_i8)
      lim = (2_i8**31)*qr
      do j = 1, 2*wnd
        if (j <= wnd .or. 2*(qr - 1)**2 <= lim) then
          xv(j) = modulo(xs(rng), 2*lim - 1) - (lim - 1)
        else
          xv(j) = modulo(xs(rng), 2*(qr - 1)**2 - lim + 1) + lim
          if (iand(j, 1_i8) == 0) xv(j) = -xv(j)
        end if
      end do
      call sredc_vec(xv, qr, mis, ovf, oin, oall)
      nwin = nwin + 2*wnd
    end do
    write(*, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'sredc random x for 2000 random odd q < 2^31: ', nwin, ' x, model mismatches ', mis, &
      ', not representable ', ovf, ', |u| >= q with |x| < 2^31 q ', oin, ', |u| >= q in all ', oall
    if (mis + oin > 0) error stop 1
    ! (e) differences a b - c d with a, b, c, d in {+-(q-1), ..., +-(q-5)}, all 10^4 choices, formed as in the code
    mis = 0; ovf = 0; oin = 0; oall = 0; nwin = 0
    deallocate(xv); allocate(xv(10**4))
    do iq = 1, size(qs)
      q = qs(iq)
      dd = [q - 1, q - 2, q - 3, q - 4, q - 5]
      j = 0
      do ia = -5, 5
        do ib = -5, 5
          do ic = -5, 5
            do id = -5, 5
              if (ia*ib*ic*id == 0) cycle
              ca = sign(dd(abs(ia)), int(ia, i8)); cb = sign(dd(abs(ib)), int(ib, i8))
              cc = sign(dd(abs(ic)), int(ic, i8)); cd = sign(dd(abs(id)), int(id, i8))
              j = j + 1
              xv(j) = int(int(ca, i4), i8)*int(int(cb, i4), i8) - int(int(cc, i4), i8)*int(int(cd, i4), i8)
            end do
          end do
        end do
      end do
      nqd = oall
      call sredc_vec(xv(1:j), q, mis, ovf, oin, oall)
      nwin = nwin + j
      write(*, '(a,i0,a,i0)') 'sredc differences at +-(q-1)..+-(q-5), q = ', q, ': |u| >= q in ', oall - nqd
    end do
    write(*, '(a,i0,a,i0,a,i0,a,i0)') 'sredc differences of products at +-(q-1)..+-(q-5), 9 moduli: ', nwin, &
      ' x, model mismatches ', mis, &
      ', not representable ', ovf, ', |u| >= q with |x| < 2^31 q ', oin
    if (mis + oin > 0) error stop 1
  end subroutine sredc_checks

  subroutine replay(qq, ca, cb, cc, cd)
    integer(i8), intent(in) :: qq, ca, cb, cc, cd
    integer(i8) :: xv(1), m0, v0, i0, o0
    integer(i4) :: u
    xv(1) = int(int(ca, i4), i8)*int(int(cb, i4), i8) - int(int(cc, i4), i8)*int(int(cd, i4), i8)
    m0 = 0; v0 = 0; i0 = 0; o0 = 0
    call sredc_vec(xv, qq, m0, v0, i0, o0)
    u = sredc(xv(1), int(qq, i4), int(signedmod(inv32(qq), 32), i4))
    write(*, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'sredc counterexample q = ', qq, ': a = ', ca, ', b = ', cb, &
      ', c = ', cc, &
      ', d = ', cd, ', x = a b - c d = ', xv(1), ', sredc = q + ', int(u, i8) - qq, '  (model mismatches ', m0, ')'
    if (o0 /= 1 .or. m0 /= 0) error stop 1
  end subroutine replay

  integer(i8) function xs(r)                   ! xorshift
    integer(i8), intent(inout) :: r
    r = ieor(r, shiftl(r, 13)); r = ieor(r, shiftr(r, 7)); r = ieor(r, shiftl(r, 17))
    xs = r
  end function xs

  subroutine ldres_checks()
    integer(i8) :: c(16384), rng, nmod, n2, bad2, r, qq, iq
    integer(i4) :: p0(16384), p1(16384), p2(16384), d(0:2), q4, qinv4, res(16384)
    integer :: i, k, bb
    rng = 2463534242_i8; nmod = 0; n2 = 0; bad2 = 0
    qq = 257
    do while (qq < 2_i8**31)
      nmod = nmod + 1
      q4 = int(qq, i4); qinv4 = int(signedmod(inv32(qq), 32), i4)
      r = modulo(2_i8**32, qq)
      do k = 0, 2
        d(k) = int(r, i4); r = modulo(r*2_i8**21, qq)
      end do
      do i = 1, 16384
        rng = ieor(rng, shiftl(rng, 13)); rng = ieor(rng, shiftr(rng, 7)); rng = ieor(rng, shiftl(rng, 17))
        if (mod(i, 2) == 0) then
          c(i) = rng
        else
          bb = int(modulo(shiftr(rng, 3), 63_i8)); c(i) = shifta(rng, 63 - bb)
        end if
      end do
      c(1:14) = [-huge(0_i8) - 1, huge(0_i8), 0_i8, 1_i8, -1_i8, qq, -qq, 2_i8**42, -2_i8**42, 2_i8**21, -2_i8**21, &
                 2_i8**42 - 1, -2_i8**42 - 1, qq*qq]
      if (qq < 2_i8**30) c(14) = 3*qq*qq
      p0 = int(iand(c, 2_i8**21 - 1), i4)
      p1 = int(iand(shifta(c, 21), 2_i8**21 - 1), i4)
      p2 = int(shifta(c, 42), i4)
      res = ldres(p0, p1, p2, d(0), d(1), d(2), q4, qinv4)
      n2 = n2 + 16384
      do i = 1, 16384
        if (abs(int(res(i), i8)) >= qq .or. modulo(int(res(i), i16) - int(c(i), i16), int(qq, i16)) /= 0) bad2 = bad2 + 1
      end do
      ! moduli: every odd q in (2^8, 2^10); 2^24 + 1 and every 8th odd number in the top 32000 below 2^25;
      ! then 200 random odd q < 2^31
      if (qq < 1023) then
        qq = qq + 2
      else if (qq == 1023) then
        qq = 2_i8**24 + 1
      else if (qq == 2_i8**24 + 1) then
        qq = 2_i8**25 - 1 - 2*1999*8
      else if (qq < 2_i8**25) then
        qq = qq + 16
        if (qq >= 2_i8**25) qq = 2_i8**25 + 1
      else
        iq = nmod - 2385
        if (iq >= 200) exit
        qq = ior(modulo(xs(rng), 2_i8**31 - 2_i8**25) + 2_i8**25, 1_i8)
      end if
    end do
    write(*, '(a,i0,a,i0,a,i0)') 'ldres: ', nmod, ' moduli, ', n2, ' coefficients, wrong (outside (-q,q) or wrong residue) ', bad2
    if (bad2 > 0) error stop 1
  end subroutine ldres_checks

  subroutine inv32_checks()
    integer(i8) :: qq, inv, nq, bad, rng, j, ns, nl
    nq = 0; bad = 0; rng = 7_i8
    ! the filter primes of pdf_init's rule for S = {2,5}, all odd q < 2^20, and 10^6 random odd q < 2^31
    ns = 0; qq = 1023
    do while (ns < 6 .and. qq > 256)
      if (usable25(qq)) then; ns = ns + 1; call inv_one(qq, nq, bad); end if
      qq = qq - 2
    end do
    nl = 0; qq = 2_i8**25 - 1
    do while (nl < 24 .and. qq > 2_i8**24)
      if (usable25(qq)) then; nl = nl + 1; call inv_one(qq, nq, bad); end if
      qq = qq - 2
    end do
    do qq = 3, 2_i8**20, 2
      call inv_one(qq, nq, bad)
    end do
    do j = 1, 1000000
      call inv_one(ior(modulo(xs(rng), 2_i8**31 - 3) + 3, 1_i8), nq, bad)
    end do
    inv = 0
    write(*, '(a,i0,a,i0)') 'inv32 / signedmod: 30 filter primes, every odd q < 2^20 and 10^6 random odd q < 2^31: ', nq, &
      ' q, wrong ', bad
    if (bad > 0) error stop 1
  end subroutine inv32_checks

  subroutine inv_one(q, nq, bad)
    integer(i8), intent(in) :: q
    integer(i8), intent(inout) :: nq, bad
    integer(i8) :: inv, s16_, s32_
    nq = nq + 1
    inv = inv32(q)
    if (inv < 0 .or. inv >= 2_i8**32 .or. modulo(int(inv, i16)*q, 2_i16**32) /= 1) bad = bad + 1
    s32_ = signedmod(inv, 32); s16_ = signedmod(inv, 16)
    if (s32_ < -2_i8**31 .or. s32_ >= 2_i8**31 .or. modulo(s32_*q, 2_i8**32) /= 1) bad = bad + 1
    if (s16_ < -2_i8**15 .or. s16_ >= 2_i8**15 .or. modulo(s16_*q, 2_i8**16) /= 1) bad = bad + 1
  end subroutine inv_one
end program arith_checks
