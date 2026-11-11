// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 Alperen Yavuz
// nf_kernel.h -- portable source of the GPU kernel (OpenCL C and CUDA C).  GPU version of the innermost GetDecics enumeration (m2 / m1
// loops of Mart52Engine_Tgt) fused with the small-prime Legendre pre-filter of
// pdfilter.f90.
//
// One work item = one (node, m2) pair.  A node is one fixed (a1..a5, m-range)
// state of the CPU loops, uploaded as 18 doubles + 35 integers.  The work item
// runs the complete m1 loop of the original code, replays its bb4..bb9 tests
// in the original order, counts every polynomial that passes them (the
// application reports that number) and writes the few that survive the
// filter to the output buffer.  Everything that can influence a branch is
// computed exactly as on the CPU:
//   * doubles: same operation order, no FMA contraction, IEEE sqrt / divide; a device without double precision does the
//     same operations in software (NF_SOFT_DOUBLE, below), with the same results bit for bit;
//   * integers: unsigned wrap-around arithmetic, like the CPU's long long;
//   * double -> int64 uses the x86 rule (out of range or NaN gives INT64_MIN).
// The pre-filter itself is exact integer arithmetic (redc below), with no floating-point operation.

// Backend specifics (address spaces, types, atomics, ...) come from nf_port_opencl.h or
// nf_port_cuda.h, which must be included first.

#define ND 18
#define NI 35
#define REC 12        // output record: node, m2, m1, bb1..bb9

// ---- doubles --------------------------------------------------------------
// The enumeration's bounds are IEEE binary64 arithmetic, exactly as on the CPU.  The kernel writes every double operation as
// one of the macros below.  Without NF_SOFT_DOUBLE they are the hardware's operators.  With NF_SOFT_DOUBLE (the OpenCL host
// defines it for a device without cl_khr_fp64: Intel Gen11 / Xe-LP / Arc A-series, Apple) nfdbl holds the bit pattern, the
// node doubles are read as such, and each operation runs in software on 64-bit integers, correctly rounded to nearest even
// like the hardware: the same results bit for bit, including signed zeros, infinities, NaN (any NaN: only comparisons and
// d2l see them) and subnormals.  Structure after Berkeley SoftFloat 3: significands with the leading bit at bit 62 and 10
// rounding bits; exponents biased, minus one where noted.  The port header supplies NF_CLZ64 and NF_MULHI64.
// pdfilter_Test/test_soft_double.cpp compares every operation with the CPU's hardware doubles.
#ifndef NF_SOFT_DOUBLE
typedef double nfdbl;
NF_INLINE i64 d2l(double x)
{
  if (!(x > -9223372036854775808.0 && x < 9223372036854775808.0)) return NF_I64_MIN;
  return (i64)x;
}
#define DADD(a, b) ((a) + (b))
#define DSUB(a, b) ((a) - (b))
#define DMUL(a, b) ((a) * (b))
#define DDIV(a, b) ((a) / (b))
#define DNEG(a) (-(a))
#define DSQRT(a) sqrt(a)
#define DCEIL(a) ceil(a)
#define DFLOOR(a) floor(a)
#define DFROMI(x) ((double)(i64)(x))
#define DGE(a, b) ((a) >= (b))
#define DLE(a, b) ((a) <= (b))
#define D2L(a) d2l(a)
#define DK(bits, val) (val)
#else
typedef u64 nfdbl;
#define SD_SIGN(a) ((a) >> 63)
#define SD_EXP(a) ((int)((a) >> 52) & 0x7FF)
#define SD_FRAC(a) ((a) & (u64)0x000FFFFFFFFFFFFFUL)
#define SD_PACK(s, e, sig) (((u64)(s) << 63) + ((u64)(e) << 52) + (sig))
#define SD_NAN ((u64)0x7FF8000000000000UL)
#define SD_HID ((u64)0x0010000000000000UL)

NF_INLINE u64 sd_srj(u64 a, int d)                // a >> d, with every bit shifted out ORed into bit 0
{
  if (d <= 0) return a;
  return d < 63 ? (a >> d) | (u64)((a << ((64 - d) & 63)) != 0) : (u64)(a != 0);
}

NF_INLINE u64 sd_roundpack(u64 s, int e, u64 sig) // e: biased exponent - 1; sig: leading bit at 62
{
  u64 rb = sig & 0x3FF;
  if ((unsigned)e >= 0x7FDu) {
    if (e < 0) { sig = sd_srj(sig, -e); e = 0; rb = sig & 0x3FF; }
    else if (e > 0x7FD || sig + 0x200 >= (u64)0x8000000000000000UL) return SD_PACK(s, 0x7FF, 0);
  }
  sig = (sig + 0x200) >> 10;
  if (rb == 0x200) sig &= ~(u64)1;
  if (!sig) e = 0;
  return SD_PACK(s, e, sig);
}

NF_INLINE u64 sd_normroundpack(u64 s, int e, u64 sig)
{
  int sh = NF_CLZ64(sig) - 1;
  e -= sh;
  if (sh >= 10 && (unsigned)e < 0x7FDu) return SD_PACK(s, sig ? e : 0, sig << (sh - 10));
  return sd_roundpack(s, e, sig << sh);
}

NF_INLINE u64 sd_addmags(u64 a, u64 b, u64 s)
{
  int ea = SD_EXP(a), eb = SD_EXP(b), d = ea - eb, ez;
  u64 fa = SD_FRAC(a), fb = SD_FRAC(b), z;
  if (d == 0) {
    if (ea == 0) return a + fb;
    if (ea == 0x7FF) return (fa | fb) ? SD_NAN : a;
    ez = ea;
    z = (2 * SD_HID + fa + fb) << 9;
  } else {
    fa <<= 9; fb <<= 9;
    if (d < 0) {
      if (eb == 0x7FF) return fb ? SD_NAN : SD_PACK(s, 0x7FF, 0);
      ez = eb;
      fa = ea ? fa + (u64)0x2000000000000000UL : fa << 1;
      fa = sd_srj(fa, -d);
    } else {
      if (ea == 0x7FF) return fa ? SD_NAN : a;
      ez = ea;
      fb = eb ? fb + (u64)0x2000000000000000UL : fb << 1;
      fb = sd_srj(fb, d);
    }
    z = (u64)0x2000000000000000UL + fa + fb;
    if (z < (u64)0x4000000000000000UL) { ez--; z <<= 1; }
  }
  return sd_roundpack(s, ez, z);
}

NF_INLINE u64 sd_submags(u64 a, u64 b, u64 s)
{
  int ea = SD_EXP(a), eb = SD_EXP(b), d = ea - eb, ez;
  u64 fa = SD_FRAC(a), fb = SD_FRAC(b), z;
  if (d == 0) {
    if (ea == 0x7FF) return SD_NAN;
    i64 df = (i64)fa - (i64)fb;
    if (df == 0) return 0;
    if (ea) ea--;
    if (df < 0) { s ^= 1; df = -df; }
    int sh = NF_CLZ64((u64)df) - 11;
    ez = ea - sh;
    if (ez < 0) { sh = ea; ez = 0; }
    return SD_PACK(s, ez, (u64)df << sh);
  }
  fa <<= 10; fb <<= 10;
  if (d < 0) {
    s ^= 1;
    if (eb == 0x7FF) return fb ? SD_NAN : SD_PACK(s, 0x7FF, 0);
    fa += ea ? (u64)0x4000000000000000UL : fa;
    fa = sd_srj(fa, -d);
    fb |= (u64)0x4000000000000000UL;
    ez = eb;
    z = fb - fa;
  } else {
    if (ea == 0x7FF) return fa ? SD_NAN : a;
    fb += eb ? (u64)0x4000000000000000UL : fb;
    fb = sd_srj(fb, d);
    fa |= (u64)0x4000000000000000UL;
    ez = ea;
    z = fa - fb;
  }
  return sd_normroundpack(s, ez - 1, z);
}

NF_INLINE u64 sd_add(u64 a, u64 b) { return SD_SIGN(a) == SD_SIGN(b) ? sd_addmags(a, b, SD_SIGN(a)) : sd_submags(a, b, SD_SIGN(a)); }
NF_INLINE u64 sd_sub(u64 a, u64 b) { return SD_SIGN(a) == SD_SIGN(b) ? sd_submags(a, b, SD_SIGN(a)) : sd_addmags(a, b, SD_SIGN(a)); }

NF_INLINE u64 sd_mul(u64 a, u64 b)
{
  u64 s = SD_SIGN(a) ^ SD_SIGN(b);
  int ea = SD_EXP(a), eb = SD_EXP(b);
  u64 fa = SD_FRAC(a), fb = SD_FRAC(b);
  if (ea == 0x7FF) { if (fa || (eb == 0x7FF && fb)) return SD_NAN; return (eb || fb) ? SD_PACK(s, 0x7FF, 0) : SD_NAN; }
  if (eb == 0x7FF) { if (fb) return SD_NAN; return (ea || fa) ? SD_PACK(s, 0x7FF, 0) : SD_NAN; }
  if (!ea) { if (!fa) return SD_PACK(s, 0, 0); int sh = NF_CLZ64(fa) - 11; ea = 1 - sh; fa <<= sh; }
  if (!eb) { if (!fb) return SD_PACK(s, 0, 0); int sh = NF_CLZ64(fb) - 11; eb = 1 - sh; fb <<= sh; }
  int ez = ea + eb - 0x3FF;
  fa = (fa | SD_HID) << 10;
  fb = (fb | SD_HID) << 11;
  u64 z = NF_MULHI64(fa, fb) | (u64)((fa * fb) != 0);
  if (z < (u64)0x4000000000000000UL) { ez--; z <<= 1; }
  return sd_roundpack(s, ez, z);
}

NF_INLINE u64 sd_div(u64 a, u64 b)
{
  u64 s = SD_SIGN(a) ^ SD_SIGN(b);
  int ea = SD_EXP(a), eb = SD_EXP(b);
  u64 fa = SD_FRAC(a), fb = SD_FRAC(b);
  if (ea == 0x7FF) { if (fa || eb == 0x7FF) return SD_NAN; return SD_PACK(s, 0x7FF, 0); }
  if (eb == 0x7FF) { if (fb) return SD_NAN; return SD_PACK(s, 0, 0); }
  if (!eb) {
    if (!fb) return (ea || fa) ? SD_PACK(s, 0x7FF, 0) : SD_NAN;
    int sh = NF_CLZ64(fb) - 11; eb = 1 - sh; fb <<= sh;
  }
  if (!ea) { if (!fa) return SD_PACK(s, 0, 0); int sh = NF_CLZ64(fa) - 11; ea = 1 - sh; fa <<= sh; }
  int ez = ea - eb + 0x3FE;
  u64 n = fa | SD_HID, dd = fb | SD_HID, q = 0;
  if (n < dd) { n <<= 1; ez--; }
  for (int i = 0; i < 55; i++) {                  // n in [dd, 2 dd): 55 quotient bits, the first one 1
    q <<= 1;
    if (n >= dd) { n -= dd; q |= 1; }
    n <<= 1;
  }
  return sd_roundpack(s, ez, (q << 8) | (u64)(n != 0));
}

NF_INLINE u64 sd_sqrt(u64 a)
{
  int ea = SD_EXP(a);
  u64 fa = SD_FRAC(a);
  if (ea == 0x7FF) return (fa || SD_SIGN(a)) ? SD_NAN : a;
  if (SD_SIGN(a)) return (ea || fa) ? SD_NAN : a;   // -0 stays -0
  if (!ea) { if (!fa) return a; int sh = NF_CLZ64(fa) - 11; ea = 1 - sh; fa <<= sh; }
  u64 m = fa | SD_HID;                            // a = m * 2^e
  int e = ea - 1075;
  if (e & 1) { m <<= 1; e--; }
  u64 hi = m >> 8, lo = m << 56, r = 0, rem = 0;  // the root of m * 2^56 (110 bits) to 55 bits, two bits at a time
  for (int i = 0; i < 55; i++) {
    int pos = 108 - 2 * i;
    u64 two = pos >= 64 ? (hi >> (pos - 64)) & 3 : (lo >> pos) & 3;
    rem = (rem << 2) | two;
    u64 t = (r << 2) | 1;
    if (rem >= t) { rem -= t; r = (r << 1) | 1; } else r <<= 1;
  }
  return sd_roundpack(0, 54 + (e - 56) / 2 + 0x3FE, (r << 8) | (u64)(rem != 0));
}

NF_INLINE u64 sd_from_i64(i64 x)
{
  u64 u = (u64)x;
  if (!(u & (u64)0x7FFFFFFFFFFFFFFFUL)) return x < 0 ? SD_PACK(1, 0x43E, 0) : 0;
  u64 s = x < 0;
  return sd_normroundpack(s, 0x43C, s ? (u64)0 - u : u);
}

NF_INLINE u64 sd_round_int(u64 a, int up)         // up = 1: ceil, 0: floor
{
  int e = SD_EXP(a);
  u64 s = SD_SIGN(a);
  if (e >= 0x433) return a;                       // |a| >= 2^52, inf, NaN
  if (e < 0x3FF) {                                // |a| < 1
    if (!(a << 1)) return a;
    if (up) return s ? SD_PACK(1, 0, 0) : (u64)0x3FF0000000000000UL;
    return s ? (u64)0xBFF0000000000000UL : 0;
  }
  u64 mask = ((u64)1 << (0x433 - e)) - 1;
  if (!(a & mask)) return a;
  if (s != (u64)up) a += mask;                    // away from zero: floor of a negative, ceil of a positive
  return a & ~mask;
}

NF_INLINE i64 sd_d2l(u64 a)                       // the x86 rule: out of range or NaN gives INT64_MIN
{
  int e = SD_EXP(a);
  if (e >= 0x43E) return NF_I64_MIN;
  if (e < 0x3FF) return 0;
  u64 m = SD_FRAC(a) | SD_HID;
  int sh = e - 0x433;
  u64 v = sh >= 0 ? m << sh : m >> (-sh);
  return SD_SIGN(a) ? -(i64)v : (i64)v;
}

NF_INLINE int sd_isnan(u64 a) { return (a << 1) > ((u64)0x7FF0000000000000UL << 1); }
NF_INLINE int sd_le(u64 a, u64 b)                 // a <= b; false when either is NaN
{
  if (sd_isnan(a) || sd_isnan(b)) return 0;
  u64 sa = SD_SIGN(a), sb = SD_SIGN(b);
  if (sa != sb) return sa || !((a | b) << 1);
  return a == b || (int)(sa ^ (u64)(a < b));
}

#define DADD(a, b) sd_add(a, b)
#define DSUB(a, b) sd_sub(a, b)
#define DMUL(a, b) sd_mul(a, b)
#define DDIV(a, b) sd_div(a, b)
#define DNEG(a) ((a) ^ (u64)0x8000000000000000UL)
#define DSQRT(a) sd_sqrt(a)
#define DCEIL(a) sd_round_int(a, 1)
#define DFLOOR(a) sd_round_int(a, 0)
#define DFROMI(x) sd_from_i64((i64)(x))
#define DGE(a, b) sd_le(b, a)
#define DLE(a, b) sd_le(a, b)
#define D2L(a) sd_d2l(a)
#define DK(bits, val) ((u64)bits##UL)
#endif

// ---- helpers -------------------------------------------------------------

// The filter is integer-only (docs/whitepaper.tex, Theorem E): unsigned 32-bit arithmetic, which wraps modulo 2^32 by the definition
// of OpenCL C and CUDA, for an odd prime q < 8192 (the table's primes are below 8192).  Residues are kept in [0, 2q): a value is 0 mod q
// exactly when it is 0 or q.  Montgomery reduction with R = 2^16: for 0 <= t < 8q^2, m = t*qneg mod 2^16 (qneg = -1/q mod 2^32, whose
// low 16 bits are -1/q mod 2^16) makes t + m*q (< 2^30, no wrap) divisible by 2^16, and (t + m*q) / 2^16 < t/2^16 + q < 2q is
// congruent to t/R mod q.
NF_INLINE nfuint redc(nfuint t, nfuint q, nfuint qneg)
{
  nfuint m = (t * qneg) & 0xFFFFu;
  return (t + m * q) >> 16;
}

// x mod q in [0, q).  Limbs of 16 bits, weights c1 = 2^16, c2 = 2^32, c3 = 2^48 mod q (all < q): the sum is at most
// (2^16 - 1)(1 + 3(q - 1)) < 2^31 for q < 8192 (it holds up to q = 10923).
NF_INLINE nfuint coef_mod(i64 x, int q, int c1, int c2, int c3, int w64)
{
  u64 u = (u64)x;
  int r = (int)(u & 0xFFFF) + (int)((u >> 16) & 0xFFFF) * c1
        + (int)((u >> 32) & 0xFFFF) * c2 + (int)(u >> 48) * c3;
  if (x < 0) r -= w64;          // x = u - 2^64
  r %= q;
  if (r < 0) r += q;
  return (nfuint)r;
}

// chi(Disc A) = -1 mod q proven?  Same algorithm as chain_general of pdfilter.f90 for one lane.  Every product of two residues and every
// lv*U(i) + (2q - lu)*V(i) is below 8q^2 and is reduced by redc, so each polynomial is its exact counterpart times a power of
// R^-1, a square; Euler's criterion runs in the Montgomery domain from one = R mod q, and -1 is negr = -R mod q or negr + q.
NF_INLINE int chi_reject(const nfuint *a, nfuint q, nfuint qneg, nfuint one, nfuint negr, int e)
{
  nfuint u[12], v[12], w[12];
  #pragma unroll
  for (int i = 0; i <= 10; i++) u[i] = a[10 - i];
  u[11] = 0u;
  #pragma unroll
  for (int i = 0; i <= 9; i++) v[i] = redc((nfuint)(10 - i) * a[10 - i], q, qneg);
  v[10] = 0u; v[11] = 0u;
  int du = 10, dv = 9, done = 0;
  nfuint mm = 1u;

  for (int it = 0; it < 40; it++) {
    nfuint lu = u[0], lv = v[0];
    #pragma unroll
    for (int i = 0; i < 11; i++)
      w[i] = (i < du) ? redc(lv * u[i + 1] + (2u * q - lu) * v[i + 1], q, qneg) : 0u;
    w[11] = 0u;
    int dw = du - 1;
    while (dw >= 0 && (w[0] == 0u || w[0] == q)) {
      #pragma unroll
      for (int i = 0; i < 11; i++) w[i] = w[i + 1];
      w[11] = 0u;
      dw--;
    }
    if (dw < 0) return 0;                       // resultant vanishes mod q: keep
    int ev = dw - du + dv;  if (ev < 0) ev = -ev;
    if (ev & 1) mm = redc(mm * lv, q, qneg);
    #pragma unroll
    for (int i = 0; i < 12; i++) u[i] = w[i];
    du = dw;
    if (du < dv) {
      #pragma unroll
      for (int i = 0; i < 12; i++) { nfuint t = u[i]; u[i] = v[i]; v[i] = t; }
      int t = du; du = dv; dv = t;
    }
    if (dv == 0) {
      if (du & 1) mm = redc(mm * v[0], q, qneg);
      done = 1;
      break;
    }
  }
  if (!done) return 0;
  nfuint acc = one, base = mm;
  for (int ee = e; ee > 0; ee >>= 1) {
    if (ee & 1) acc = redc(acc * base, q, qneg);
    base = redc(base * base, q, qneg);
  }
  return acc == negr || acc == negr + q;
}

// bb[1..10] as c[1..10] (leading coefficient 1 implied).  1 = survives all primes.  Row of qt: see primeRow in nf_gpu_common.cpp.
NF_INLINE int survives(const i64 *c, NF_C int *qt, nfuint np)
{
  for (nfuint p = 0; p < np; p++) {
    int q = qt[p * 8 + 0], c1 = qt[p * 8 + 1], c2 = qt[p * 8 + 2], c3 = qt[p * 8 + 3];
    int w64 = qt[p * 8 + 4], e = qt[p * 8 + 5];
    nfuint qneg = (nfuint)qt[p * 8 + 6], negr = (nfuint)qt[p * 8 + 7];
    nfuint a[11];
    #pragma unroll
    for (int k = 0; k <= 9; k++) a[k] = coef_mod(c[10 - k], q, c1, c2, c3, w64);
    a[10] = 1u;
    if (chi_reject(a, (nfuint)q, qneg, (nfuint)c1, negr, e)) return 0;
  }
  return 1;
}

// ---- start-up self-test --------------------------------------------------

// The filter's functions on inputs the host chose, with the device's compiler and arithmetic; the host compares the output with
// values it computes independently (nfSelfTestVerify, nf_gpu_common.cpp).  Same functions as nf_enum_t2, only the surrounding
// loop differs.  Work items, in this order:
//   np * nt       redc(ts[p * nt + i]) for prime p (0 <= t < 8q^2)
//   np * nx       coef_mod(xs[i]) for prime p
//   np * npoly    chi_reject of polynomial j at prime p (1 = reject), loaded as in survives
//   npoly         survives() of polynomial j over the whole table
// polys holds c1..c10 of each polynomial (leading coefficient 1 implied).
NF_KERNEL void nf_selftest_t2(NF_G const nfuint *ts, nfuint nt, NF_G const i64 *xs, nfuint nx, NF_G const i64 *polys, nfuint npoly,
                              NF_C int *qt, nfuint np, NF_G nfuint *out, nfuint total)
{
  nfuint g = (nfuint)NF_GID();
  if (g >= total) return;
  nfuint n1 = np * nt, n2 = np * nx, n3 = np * npoly;
  if (g < n1) {
    nfuint p = g / nt;
    out[g] = redc(ts[g], (nfuint)qt[p * 8 + 0], (nfuint)qt[p * 8 + 6]);
    return;
  }
  if (g < n1 + n2) {
    nfuint h = g - n1, p = h / nx, i = h % nx;
    out[g] = coef_mod(xs[i], qt[p * 8 + 0], qt[p * 8 + 1], qt[p * 8 + 2], qt[p * 8 + 3], qt[p * 8 + 4]);
    return;
  }
  if (g < n1 + n2 + n3) {
    nfuint h = g - n1 - n2, p = h / npoly, j = h % npoly;
    int q = qt[p * 8 + 0], c1 = qt[p * 8 + 1], c2 = qt[p * 8 + 2], c3 = qt[p * 8 + 3], w64 = qt[p * 8 + 4], e = qt[p * 8 + 5];
    nfuint a[11];
    for (int k = 0; k <= 9; k++) a[k] = coef_mod(polys[j * 10 + 9 - k], q, c1, c2, c3, w64);
    a[10] = 1u;
    out[g] = (nfuint)chi_reject(a, (nfuint)q, (nfuint)qt[p * 8 + 6], (nfuint)c1, (nfuint)qt[p * 8 + 7], e);
    return;
  }
  nfuint j = g - n1 - n2 - n3;
  i64 c[11];
  c[0] = 1;
  for (int k = 1; k <= 10; k++) c[k] = polys[j * 10 + k - 1];
  out[g] = (nfuint)survives(c, qt, np);
}

// ---- kernel ---------------------------------------------------------------

// The layout of the prime table qt (primeRow in nf_gpu_common.cpp).  Format 2: the integer filter above; format 1 was the
// floating-point filter (row[6] the bits of 1.0f/q).  The entry point carries the format in its name, so that a program never
// runs device code built for another layout (a CUDA fat binary or MUSA library left over from an older build): the host looks
// up NF_KERNEL_NAME (nf_gpu_backend.h), and the MUSA library reports NF_TABLE_FORMAT.  Change all of them together.
#define NF_TABLE_FORMAT 2

NF_KERNEL void nf_enum_t2(NF_G const nfdbl *nd, NF_G const i64 *ni,
                      NF_G const u64 *cum, nfuint nnodes, nfuint nodebase, u64 total,
                      NF_C int *qt, nfuint nprime,
                      NF_G i64 *out, NF_G nfuint *outcnt, nfuint outcap,
                      NF_G nfuint *cnt)
{
  u64 gid = NF_GID();
  if (gid >= total) return;

  // Node = last k with cum[k] <= gid.
  nfuint lo = 0, hi = nnodes;
  while (hi - lo > 1) {
    nfuint mid = (lo + hi) >> 1;
    if (cum[mid] <= gid) lo = mid; else hi = mid;
  }
  nfuint node = lo;
  NF_G const nfdbl *D = nd + (size_t)(nodebase + node) * ND;
  NF_G const i64 *I = ni + (size_t)(nodebase + node) * NI;

  i64 m2 = I[0] + (i64)(gid - cum[node]);

  // m1 range exactly as in the CPU code.
  nfdbl m2p = DSUB(D[1], DFROMI(m2));
  nfdbl B1 = DSQRT(DSUB(D[2], DMUL(DMUL(D[3], m2p), m2p)));
  i64 m1_L = D2L(DCEIL(DADD(D[0], DDIV(DADD(DNEG(B1), DMUL(D[4], m2p)), D[5]))));
  i64 m1_U = D2L(DFLOOR(DADD(D[0], DDIV(DADD(B1, DMUL(D[4], m2p)), D[5]))));

  const nfdbl lam1 = D[6], lam2 = D[7], lam1c = D[8], lam2c = D[9];
  const nfdbl T1 = D[10], T2 = D[11];
  const nfdbl t5 = D[12], t6 = D[13], t7 = D[14], t8 = D[15], t9 = D[16], bb8_Upre = D[17];

  const u64 M11 = (u64)I[2], M12 = (u64)I[3], M21 = (u64)I[4], M22 = (u64)I[5];
  const u64 cv41 = (u64)I[6], cv42 = (u64)I[7];
  const u64 a11 = (u64)I[8], a12 = (u64)I[9], a21 = (u64)I[10], a22 = (u64)I[11];
  const u64 a31 = (u64)I[12], a32 = (u64)I[13], a51 = (u64)I[14], a52 = (u64)I[15];
  const u64 Trw = (u64)I[16], Nmw = (u64)I[17];
  const i64 bb9L0 = I[18], bb9U0 = I[19];
  const u64 bb4pre2 = (u64)I[20];
  const i64 bb4L = I[21], bb4U = I[22];
  const u64 bb5pre2 = (u64)I[23], bb6pre2 = (u64)I[24], bb7pre = (u64)I[25], bb8pre = (u64)I[26];
  const u64 bb1 = (u64)I[27], bb2 = (u64)I[28], bb3 = (u64)I[29];
  const u64 S1 = (u64)I[30], S2 = (u64)I[31], S3 = (u64)I[32], Sum4 = (u64)I[33];
  const i64 bb10 = I[34];

  const u64 base41 = M12 * (u64)m2 + cv41;
  const u64 base42 = M22 * (u64)m2 + cv42;

  nfuint count = 0;

  for (i64 m1 = m1_L; m1 <= m1_U; ++m1) {
    u64 a41 = M11 * (u64)m1 + base41;
    u64 a42 = M21 * (u64)m1 + base42;
    u64 a41sq = a41 * a41, a42sq = a42 * a42, a41a42 = a41 * a42;
    nfdbl x1 = DFROMI(a41sq), x2 = DFROMI(a41a42), x3 = DFROMI(a42sq);
    nfdbl sd1 = DADD(DADD(x1, DMUL(lam1, x2)), DMUL(lam2, x3));
    nfdbl sd2 = DADD(DADD(x1, DMUL(lam1c, x2)), DMUL(lam2c, x3));
    if (!(DGE(sd1, DK(0, 0.0)) && DLE(sd1, T1) && DGE(sd2, DK(0, 0.0)) && DLE(sd2, T2))) continue;

    u64 bb9 = 2 * a41 * a51 + (a41 * a52 + a42 * a51) * Trw + 2 * a42 * a52 * Nmw;
    if (!((i64)bb9 >= bb9L0 && (i64)bb9 <= bb9U0)) continue;

    u64 bb4 = bb4pre2 + 2 * a41 + a42 * Trw;
    if (!((i64)bb4 >= bb4L && (i64)bb4 <= bb4U)) continue;

    u64 bb5 = bb5pre2 + 2 * a11 * a41 + (a11 * a42 + a12 * a41) * Trw + 2 * a12 * a42 * Nmw;
    u64 s4 = (u64)(-4L) * bb4 - Sum4;
    u64 Sum = bb4 * S1 + bb3 * S2 + bb2 * S3 + bb1 * s4;
    if (!((i64)bb5 >= D2L(DCEIL(DDIV(DSUB(DNEG(t5), DFROMI(Sum)), DK(0x4014000000000000, 5.0)))) &&
          (i64)bb5 <= D2L(DFLOOR(DDIV(DSUB(t5, DFROMI(Sum)), DK(0x4014000000000000, 5.0)))))) continue;

    u64 bb6 = bb6pre2 + 2 * a21 * a41 + (a21 * a42 + a22 * a41) * Trw + 2 * a22 * a42 * Nmw;
    u64 s5 = (u64)(-5L) * bb5 - Sum;
    Sum = bb5 * S1 + bb4 * S2 + bb3 * S3 + bb2 * s4 + bb1 * s5;
    if (!((i64)bb6 >= D2L(DCEIL(DDIV(DSUB(DNEG(t6), DFROMI(Sum)), DK(0x4018000000000000, 6.0)))) &&
          (i64)bb6 <= D2L(DFLOOR(DDIV(DSUB(t6, DFROMI(Sum)), DK(0x4018000000000000, 6.0)))))) continue;

    u64 bb7 = bb7pre + 2 * a31 * a41 + (a31 * a42 + a32 * a41) * Trw + 2 * a32 * a42 * Nmw;
    u64 s6 = (u64)(-6L) * bb6 - Sum;
    Sum = bb6 * S1 + bb5 * S2 + bb4 * S3 + bb3 * s4 + bb2 * s5 + bb1 * s6;
    if (!((i64)bb7 >= D2L(DCEIL(DDIV(DSUB(DNEG(t7), DFROMI(Sum)), DK(0x401C000000000000, 7.0)))) &&
          (i64)bb7 <= D2L(DFLOOR(DDIV(DSUB(t7, DFROMI(Sum)), DK(0x401C000000000000, 7.0)))))) continue;

    u64 bb8 = bb8pre + a41sq + a41a42 * Trw + a42sq * Nmw;
    u64 s7 = (u64)(-7L) * bb7 - Sum;
    Sum = bb7 * S1 + bb6 * S2 + bb5 * S3 + bb4 * s4 + bb3 * s5 + bb2 * s6 + bb1 * s7;
    if (!((i64)bb8 >= D2L(DCEIL(DDIV(DSUB(DNEG(t8), DFROMI(Sum)), DK(0x4020000000000000, 8.0)))) &&
          (i64)bb8 <= D2L(DFLOOR(DDIV(DSUB(t8, DFROMI(Sum)), DK(0x4020000000000000, 8.0)))))) continue;

    u64 s8 = (u64)(-8L) * bb8 - Sum;
    Sum = bb8 * S1 + bb7 * S2 + bb6 * S3 + bb5 * s4 + bb4 * s5 + bb3 * s6 + bb2 * s7 + bb1 * s8;
    if (!((i64)bb9 >= D2L(DCEIL(DDIV(DSUB(DNEG(t9), DFROMI(Sum)), DK(0x4022000000000000, 9.0)))) &&
          (i64)bb9 <= D2L(DFLOOR(DDIV(DSUB(t9, DFROMI(Sum)), DK(0x4022000000000000, 9.0)))))) continue;

    nfdbl b9 = DFROMI(bb9);
    nfdbl tmpdbl = DDIV(DMUL(DMUL(DK(0x3FE0000000000000, 0.5), b9), b9), DFROMI(bb10));
    if (!((i64)bb8 >= D2L(DCEIL(DSUB(tmpdbl, bb8_Upre))) &&
          (i64)bb8 <= D2L(DFLOOR(DADD(tmpdbl, bb8_Upre))))) continue;

    // A polynomial of the original enumeration.
    count++;

    i64 c[11];
    c[1] = (i64)bb1; c[2] = (i64)bb2; c[3] = (i64)bb3; c[4] = (i64)bb4; c[5] = (i64)bb5;
    c[6] = (i64)bb6; c[7] = (i64)bb7; c[8] = (i64)bb8; c[9] = (i64)bb9; c[10] = bb10;
    if (!survives(c, qt, nprime)) continue;

    nfuint idx = NF_ATOMIC_INC(outcnt);
    if (idx < outcap) {
      NF_G i64 *o = out + (size_t)idx * REC;
      o[0] = (i64)(nodebase + node);  o[1] = m2;  o[2] = m1;
      o[3] = c[1]; o[4] = c[2]; o[5] = c[3]; o[6] = c[4]; o[7] = c[5];
      o[8] = c[6]; o[9] = c[7]; o[10] = c[8]; o[11] = c[9];
    }
  }
  cnt[gid] = count;
}
