// test_pdfilter_oracle.cpp -- tests the pre-filter (pdfilter.f90) against an exact oracle.
//
// WHAT IS CHECKED.  For every polynomial f and every prime q of the filter (30, below 2^13, in the order of the filter's
// stages) an independent oracle computes Disc(f) mod q -- the classical Euclidean algorithm over F_q with inverses, which
// shares nothing with the filter's pseudo-remainder chains, its Montgomery arithmetic or its coefficient loading -- and the
// Legendre symbol.  The filter must reject f exactly when some symbol is -1, and (with the test hook, below) at the first
// such stage.  So both kinds of error are seen: a valid polynomial rejected (all that the "must pass" tests could see) and
// a polynomial that should be rejected kept (rounding down, -Ofast, a filter switched off).
// The oracle is checked too: against PARI's poldisc on a sample, against the square class that the constructed families
// are known to have, and for invariance under the transformations used below.  A disagreement there is a harness bug.
//
// WHY STRATA.  A random polynomial is rejected at the first prime with probability 1/2, so a fault at stage j stays masked
// unless the input reaches stage j; and random polynomials almost never meet the degenerate Euclid chains that the filter's
// general path (chain_general) exists for.  Each generator is aimed at a part of the filter that random input misses:
//   rand     dense coefficients of random bit length 0..63 and the int64 extremes (the coefficient loading)
//   small    small, sparse coefficients, like the application's candidates
//   even     g(x^2) with g random (Eric's second test without the restriction to valid polynomials)
//   valid    polynomials kept at every stage (factor products, symmetric constructions, depressed h(x^2), Tschirnhaus
//            transforms, then x -> c x + t, x -> -x): Disc = +-(S-smooth) * square
//   target   rejected at a chosen stage J = 1..30, or kept (J = 0): Disc lies in the square class of an integer D whose
//            Legendre symbols are +1 at the stages before J and -1 at J (D = D1*D2, matched on its vector of symbols);
//            built as g(x^2) with g(0) = D, a depressed h(x^2) with constant D, or x^10 + b x^5 + c with b^2 - 4c = D,
//            then transformed (Tschirnhaus, x -> c x + t, x -> -x)
//   crt      coefficients built by the Chinese remainder theorem over the first small primes (as many as keep the modulus
//            below 2^62: four of the 13-bit ones) from residues that pass there (coefficients up to 63 bits)
//   crtdeg   the same with a degenerate (non-generic) Euclid chain planted at one of these primes
//   crtdegL  a degenerate chain planted at one of the first large primes (CRT over that prime and as many small ones as
//            fit, at most 3, then sampled until the polynomial reaches that stage; only the large stages that sampling
//            reaches in about 2^12 tries are planted): for filters with large primes (-oldrule); the filter described
//            here has none, so the stratum is off
//   rep      a repeated factor over Z (Disc = 0: kept everywhere)
//   real     the application's passing polynomials (-real FILE) and transformed copies
//   cand     the application's candidates as the filter received them (-cand FILE; capshim.c captures them)
//
// OTHER CHECKS.  Every block goes through both builds (AVX2/FMA and baseline) and then once more transformed (per
// polynomial f(x+t), f(-x) or c^10 f(x/c), which leave every Legendre symbol and every degree sequence unchanged),
// permuted and cut into other block sizes: verdicts and stages must not change (lane, block and neighbour independence,
// a check that needs no oracle).  With the test hook of pdfilter.f90 (built when PDF_STAGES_NAME is defined, as the test
// build scripts do) the stage of every rejection and the filter's counters per stage are read: each stage must receive
// exactly the lanes the oracle says reach it, chain_block must hand exactly the lanes with a non-generic chain to
// chain_general, and the counters give a coverage report per stage.  The hook also reports the filter's primes, which
// the oracle then uses (they are checked to be usable and compared with pdf_init's rule).  With REF=<directory> in
// build-linux.sh a second filter is linked and compared with the first and with the oracle (differential test; it
// must have the first build's primes, which the floating-point filter of 327b04f no longer has: its are those of -oldrule).  -rm up|down|zero runs every filter call under another
// rounding mode.
//
// Output: a coverage table per stratum (from the oracle), a table per stage (from the hook), the failures per build, the
// oracle checks, and as the last lines "Total polynomials tested", "False rejections", "Missed rejections" and "Total
// failures" for the first build (AVX2/FMA if the CPU has it, the baseline build with NF_FORCE_BASE=1; fault_checks.sh reads these),
// then the failures of all builds.  Exit status 0: no failure in any build; 1: failures; 2: a harness bug (oracle or generator).
//
// Build and run (this directory):  TEST=test_pdfilter_oracle.cpp bash build-linux.sh && build/test_pdfilter_oracle 100 [seed]
//   100 = 100 blocks' worth of 10,000 polynomials (block sizes vary); bash fault_checks.sh 100 test_pdfilter_oracle.cpp
//   gives the fault table.  build/test_pdfilter_oracle -h lists the options.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <ctime>
#include <cmath>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <fenv.h>
#if defined(__SSE__)
#include <xmmintrin.h>
#endif

typedef long long i64;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef unsigned short u16;
#ifdef __SIZEOF_INT128__
typedef __int128 i128;
#else
#error "this test needs __int128 (gcc or clang on a 64-bit target)"
#endif

// The filter, as in TgtMartinet.h: two builds on x86-64 (AVX2+FMA and baseline), one elsewhere.
#if defined(__x86_64__) || defined(_M_X64)
  #define NF_PDF_DISPATCH 1
  extern "C" void pdf_init_v3(int, const int*);
  extern "C" void pdf_filter_v3(const long long*, int, signed char*);
  extern "C" void pdf_init_base(int, const int*);
  extern "C" void pdf_filter_base(const long long*, int, signed char*);
#else
  extern "C" void pdf_init(int, const int*);
  extern "C" void pdf_filter(const long long*, int, signed char*);
#endif

// Optional entry points.  These weak definitions stand in when the real ones are not linked:
//   pdf_stages_*   the filter's test hook; the stand-in reports 0 stages and the test then checks the verdicts only
//   pdf_*_ref      a second filter for the differential check (build-linux.sh REF=<directory>)
static int refAbsent = 0;
extern "C" {
__attribute__((weak)) void pdf_stages_v3(int, int *, long long *, long long *, int *nst, int) { *nst = 0; }
__attribute__((weak)) void pdf_stages_base(int, int *, long long *, long long *, int *nst, int) { *nst = 0; }
__attribute__((weak)) void pdf_stages(int, int *, long long *, long long *, int *nst, int) { *nst = 0; }
__attribute__((weak)) void pdf_init_ref(int, const int *) { refAbsent = 1; }
__attribute__((weak)) void pdf_filter_ref(const long long *, int, signed char *) {}
}

#include "pari.h"      // NB: on mingw64 this turns `long` into `long long`; below, use i64 / u64 and never the word `long`

namespace pdt {

struct Poly { i64 c[11]; };          // c[0] = 1, c[k] = coefficient of x^(10-k): the layout of the filter's buffer
const int MAXST = 30;                // MAXS + MAXQ of pdfilter.f90: the size of the hook's arrays

// ---------------------------------------------------------------- random numbers (one splitmix64 stream, seeded from the seed)
static u64 rs = 1;
static u64 rnd64() { rs += 0x9E3779B97F4A7C15ull; u64 z = rs; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }
static int rnd_int(int n) { return (int)(rnd64() % (u64)n); }                           // 0 .. n-1
static i64 rnd_range(i64 a, i64 b) { return a + (i64)(rnd64() % (u64)(b - a + 1)); }
static double rnd_unif() { return (rnd64() >> 11) * (1.0 / 9007199254740992.0); }
static i64 rnd_bits(int bits)        // uniform magnitude below 2^bits (bits = 0..63), random sign
{
  if (bits <= 0) return 0;
  i64 m = (i64)(rnd64() >> (64 - bits));
  return (rnd64() & 1) ? -m : m;
}

// ---------------------------------------------------------------- another rounding mode, switched on only around the filter calls
static int hostile = 0;              // 0 nearest, 1 up, 2 down, 3 toward zero
static void roundOn()
{
  if (!hostile) return;
  fesetround(hostile == 1 ? FE_UPWARD : hostile == 2 ? FE_DOWNWARD : FE_TOWARDZERO);
#if defined(__SSE__)
  _mm_setcsr((_mm_getcsr() & ~0x6000u) | (hostile == 1 ? 0x4000u : hostile == 2 ? 0x2000u : 0x6000u));
#endif
}
static void roundOff()
{
  if (!hostile) return;
  fesetround(FE_TONEAREST);
#if defined(__SSE__)
  _mm_setcsr(_mm_getcsr() & ~0x6000u);
#endif
}

// ---------------------------------------------------------------- arithmetic mod q (q < 2^31, so products fit in 64 bits)
static u64 powm(u64 b, u64 e, u64 q) { u64 r = 1; b %= q; while (e) { if (e & 1) r = r * b % q; b = b * b % q; e >>= 1; } return r; }
static u64 invm(u64 a, u64 q)        // a^-1 mod q (extended Euclid), a != 0 mod q, q prime
{
  i64 r0 = (i64)q, r1 = (i64)(a % q), s0 = 0, s1 = 1;
  while (r1) { i64 d = r0 / r1, t = r0 - d * r1; r0 = r1; r1 = t; t = s0 - d * s1; s0 = s1; s1 = t; }
  return (u64)(s0 < 0 ? s0 + (i64)q : s0);
}
static u64 modq(i64 c, u64 q) { i64 r = c % (i64)q; return (u64)(r < 0 ? r + (i64)q : r); }
// Barrett reduction for the oracle's inner loops: x < 2^64, m = floor((2^64 - 1)/q) (= floor(2^64/q), q odd), so
// floor(x m / 2^64) is floor(x/q) or one less, and one correction gives x mod q.
struct Fq { u64 q, m; explicit Fq(u64 qq) : q(qq), m(~0ull / qq) {} };
static inline u64 red(u64 x, const Fq &F) { u64 t = (u64)(((unsigned __int128)x * F.m) >> 64), r = x - t * F.q; return r >= F.q ? r - F.q : r; }
static u64 powf(u64 b, u64 e, const Fq &F) { u64 r = 1; b = red(b, F); while (e) { if (e & 1) r = red(r * b, F); b = red(b * b, F); e >>= 1; } return r; }

// ---------------------------------------------------------------- the filter's primes
static std::vector<int> Sset;
static std::vector<u64> FPr;         // the prime of every stage, in stage order
static int NST = 0, NS = 0;          // number of stages, number of small-prime stages
static bool primeU(u64 n) { if (n < 2) return false; for (u64 d = 2; d * d <= n; d++) if (n % d == 0) return false; return true; }
static bool usablePrime(u64 q)       // q = 1 mod 4 prime and every p in S a nonzero square mod q (pdfilter.f90's usable)
{
  if (q % 4 != 1 || !primeU(q)) return false;
  for (size_t k = 0; k < Sset.size(); k++) { if (q % (u64)Sset[k] == 0 || powm((u64)Sset[k], (q - 1) / 2, q) != 1) return false; }
  return true;
}
// pdf_init's rule: 30 below 2^13 (and above 256) from the top.  oldRule: the rule of the filters before 2026-10-29, among
// them the floating-point filter of 327b04f: 6 below 2^10, then 24 "large" ones below 2^25 (above 2^24).  (In between, until
// 2026-11-06, the integer filter took 24 below 2^13 and 6 below 2^25.)
static bool oldRule = false;
static std::vector<u64> specPrimes()
{
  std::vector<u64> pr;
  int ns = 0, nq = 0, maxs = oldRule ? 6 : 30, maxq = oldRule ? 24 : 0;
  for (u64 q = oldRule ? 1023 : 8191; ns < maxs && q > 256; q -= 2) if (usablePrime(q)) { pr.push_back(q); ns++; }
  for (u64 q = (1ull << 25) - 1; nq < maxq && q > (1ull << 24); q -= 2) if (usablePrime(q)) { pr.push_back(q); nq++; }
  return pr;
}

// ---------------------------------------------------------------- the oracle
// Res(A, B) over F_q by the Euclidean algorithm with remainders:  Res(A, B) = (-1)^(deg A deg B) Res(B, A) and
// Res(B, A) = lc(B)^(deg A - deg r) Res(B, r) for r = A mod B;  Res(A, c) = c^deg A for a constant c.  For monic f of
// degree 10, Disc(f) = (-1)^45 Res(f, f').  Also the degrees of the remainders: the chain is "generic" when they are
// 8, 7, ..., 1 (a zero remainder in place of the constant is allowed: Res = 0), which is exactly when chain_block keeps
// the lane; pat has bit d set for every remainder of degree d (d <= 8) and bit 9 when a zero remainder ends the chain.
struct ModInfo { u64 disc; int chi; bool generic; u16 pat; };
static ModInfo oracleRes(u64 *a, u64 *b, const Fq &F)   // a: monic degree 10 (a[i] = coefficient of x^i, in [0, q)), b = a'
{
  const u64 q = F.q;
  int n = 10, m = 9;
  u64 res = 1;
  bool neg = false, generic = true;
  u16 pat = 0;
  for (;;) {
    if (m == 0) { res = red(res * powf(b[0], (u64)n, F), F); break; }
    u64 il = invm(b[m], q);
    for (int i = n; i >= m; i--) {
      u64 t = red(a[i] * il, F);
      if (t) for (int k = 0; k <= m; k++) { u64 s = red(t * b[k], F), &x = a[i - m + k]; x = x >= s ? x - s : x + q - s; }
    }
    int k = m - 1;
    while (k >= 0 && a[k] == 0) k--;
    if (k != m - 1 && m - 1 >= 1) generic = false;
    if (k < 0) { pat |= 1u << 9; res = 0; break; }
    pat |= (u16)(1u << k);
    if ((n & 1) && (m & 1)) neg = !neg;
    res = red(res * powf(b[m], (u64)(n - k), F), F);
    u64 *t = a; a = b; b = t; n = m; m = k;
  }
  if (neg) res = (q - res) % q;
  ModInfo r;
  r.disc = (q - res) % q;
  u64 e = r.disc ? powf(r.disc, (q - 1) / 2, F) : 0;
  r.chi = r.disc == 0 ? 0 : (e == 1 ? 1 : -1);
  r.generic = generic;
  r.pat = pat;
  return r;
}
static ModInfo oracleRes(const u64 *g, const Fq &F)     // g: monic degree 10 residues
{
  const u64 q = F.q;
  u64 a[11], b[11];
  for (int i = 0; i <= 10; i++) a[i] = g[i];
  for (int i = 0; i <= 9; i++) b[i] = (u64)(i + 1) * a[i + 1] % q;
  b[10] = 0;
  return oracleRes(a, b, F);
}
static ModInfo oracleRes(const u64 *g, u64 q) { return oracleRes(g, Fq(q)); }
static std::vector<Fq> FF;          // the stages' primes, for the Barrett reduction
static ModInfo oracleMod(const Poly &p, const Fq &F)
{
  u64 g[11];
  for (int i = 0; i <= 10; i++) g[i] = modq(p.c[10 - i], F.q);
  return oracleRes(g, F);
}
// The verdict: the first stage whose symbol is -1 (0: kept).  The stages are evaluated up to that one.
struct OracleOut {
  int stage, reached;                // reached = number of stages evaluated
  u32 nongen, zero, neg;             // bit j-1: at stage j the chain is non-generic / Disc = 0 mod q / symbol -1
  u16 pat[MAXST];
};
static OracleOut oracle(const Poly &p, int upto = 0)   // upto > 0: evaluate that many stages even past the first -1
{
  OracleOut o;
  o.stage = 0; o.reached = 0; o.nongen = o.zero = o.neg = 0;
  for (int j = 1; j <= NST; j++) {
    ModInfo m = oracleMod(p, FF[j - 1]);
    o.reached = j;
    o.pat[j - 1] = m.pat;
    if (!m.generic) o.nongen |= 1u << (j - 1);
    if (m.chi == 0) o.zero |= 1u << (j - 1);
    if (m.chi < 0) { o.neg |= 1u << (j - 1); if (!o.stage) o.stage = j; }
    if (o.stage && j >= upto) break;
  }
  return o;
}

// ---------------------------------------------------------------- PARI helpers
static GEN LIM, NLIM;                // +-2^63
static bool fitsGEN(GEN a) { return cmpii(a, LIM) < 0 && cmpii(a, NLIM) >= 0; }
static GEN polyToGEN(const Poly &p)
{
  GEN f = cgetg(13, t_POL);
  f[1] = evalsigne(1) | evalvarn(0);
  for (int i = 0; i <= 10; i++) gel(f, 2 + i) = stoi((i64)p.c[10 - i]);
  return f;
}
static bool genToPoly(GEN f, Poly &p)
{
  if (f == NULL || typ(f) != t_POL || degpol(f) != 10) return false;
  for (int i = 0; i <= 10; i++) {
    GEN a = gel(f, 2 + i);
    if (typ(a) != t_INT || !fitsGEN(a)) return false;       // coefficients must fit the int64 the filter takes
    p.c[10 - i] = (i64)itos(a);
  }
  return p.c[0] == 1;
}
// Disc = +-(S-smooth) * square, as upstream's test decides
static bool passesUpstream(GEN D)
{
  if (gequal0(D)) return false;
  pari_sp av = avma;
  GEN m = absi(D), r;
  for (size_t k = 0; k < Sset.size(); k++) { (void)Z_pvalrem(m, utoipos((u64)Sset[k]), &r); m = r; }
  bool ok = Z_issquare(m) != 0;
  avma = av;
  return ok;
}

// ---------------------------------------------------------------- transformations that keep every Legendre symbol (from test_pdfilter_full.cpp)
// Tschirnhaus: the characteristic polynomial of y(alpha), y of random degree 1..dmax with coefficients in [-B, B];
// Disc(new) = det(M)^2 Disc(f) for an integer matrix M.
static GEN tsOwn(GEN f, int dmax, int B)
{
  for (int attempt = 0; attempt < 100; attempt++) {
    pari_sp av = avma;
    int d = 1 + rnd_int(dmax);
    GEN g = cgetg(d + 3, t_POL);
    g[1] = evalsigne(1) | evalvarn(0);
    for (int i = 0; i <= d; i++) {
      i64 v = rnd_range(-B, B);
      while (i == d && v == 0) v = rnd_range(-B, B);
      gel(g, 2 + i) = stoi(v);
    }
    GEN u = RgXQ_charpoly(g, f, 0);
    if (!gequal0(poldisc0(u, -1))) return u;
    avma = av;
  }
  return NULL;
}
// exact affine maps on the coefficients (128 bit): x -> x + t, x -> c x (monic form c^10 f(x/c)), x -> -x.
// Disc(f(x+t)) = Disc(f(-x)) = Disc(f), Disc(c^10 f(x/c)) = c^90 Disc(f); the Euclid chains mod q (q not dividing c)
// have the same degrees.
static bool affine(Poly &p, i64 t, int c, bool flip)
{
  i128 e[11];
  for (int i = 0; i <= 10; i++) e[i] = p.c[10 - i];
  if (flip) for (int i = 1; i <= 10; i += 2) e[i] = -e[i];
  if (c != 1) { i128 pw = 1; for (int i = 10; i >= 0; i--) { e[i] *= pw; pw *= c; } }
  if (t != 0) {
    const i128 big = (i128)1 << 120;                       // |t| < 2^7, so |t * e| < 2^127: the guard keeps every step inside 128 bits
    for (int i = 0; i < 10; i++) for (int j = 9; j >= i; j--) {
      e[j] += (i128)t * e[j + 1];
      if (e[j] > big || e[j] < -big) return false;
    }
  }
  const i128 i64max = (i128)9223372036854775807ll, i64min = -i64max - 1;
  for (int i = 0; i <= 10; i++) if (e[i] > i64max || e[i] < i64min) return false;
  for (int i = 0; i <= 10; i++) p.c[10 - i] = (i64)e[i];
  return p.c[0] == 1;
}
static bool scramble(Poly &p)
{
  bool flip = rnd_int(2) == 0;
  int c = (rnd_int(3) == 0) ? 2 + rnd_int(2) : 1;
  i64 t = 0;
  if (rnd_unif() < 0.8) { double mx = std::exp(rnd_unif() * std::log(90.0)); i64 m = (i64)std::ceil(mx); t = rnd_range(-m, m); }
  return affine(p, t, c, flip);
}
static void scrambleOrKeep(Poly &p, double prob) { Poly q = p; if (rnd_unif() < prob && scramble(q)) p = q; }

// ---------------------------------------------------------------- polynomials over Z (vectors, a[i] = coefficient of x^i)
static bool polyMul(const std::vector<i64> &a, const std::vector<i64> &b, std::vector<i64> &r)
{
  std::vector<i128> h(a.size() + b.size() - 1, 0);
  for (size_t i = 0; i < a.size(); i++) for (size_t j = 0; j < b.size(); j++) h[i + j] += (i128)a[i] * b[j];
  r.resize(h.size());
  const i128 lim = (i128)4000000000000000000ll;
  for (size_t i = 0; i < h.size(); i++) { if (h[i] > lim || h[i] < -lim) return false; r[i] = (i64)h[i]; }
  return true;
}
static void fromVec(const std::vector<i64> &f, Poly &p) { for (int i = 0; i <= 10; i++) p.c[10 - i] = f[i]; }
static i64 smoothVal(int maxExp)     // a random product of the primes of S
{
  i64 v = 1;
  for (size_t k = 0; k < Sset.size(); k++) { int e = rnd_int(maxExp + 1); for (int i = 0; i < e; i++) v *= Sset[k]; }
  return v;
}

// ---------------------------------------------------------------- valid families (from test_pdfilter_full.cpp)
static GEN seed0;                    // x^10 - 12 x^5 - 4, the seed of Eric's first test
static bool seed0ok = false;
struct Factor { std::vector<i64> a; int deg; };
static std::vector<Factor> pool;
static void addFactor(const std::vector<i64> &a)          // keep it if its own discriminant has the right class
{
  Factor F; F.a = a; F.deg = (int)a.size() - 1;
  pari_sp av = avma;
  GEN f = cgetg(F.deg + 3, t_POL); f[1] = evalsigne(1) | evalvarn(0);
  for (int i = 0; i <= F.deg; i++) gel(f, 2 + i) = stoi(a[i]);
  bool ok = (F.deg == 1) || passesUpstream(poldisc0(f, -1));
  avma = av;
  if (ok) pool.push_back(F);
}
static void buildPool()
{
  pool.clear();
  pari_sp av0 = avma;
  for (int r = -9; r <= 9; r++) addFactor({(i64)-r, 1});
  for (int b = -6; b <= 6; b++) for (int k = 1; k <= 8; k++) for (int sg = -1; sg <= 1; sg += 2) for (int rep = 0; rep < 6; rep++) {
    i64 d = (i64)sg * smoothVal(2) * k * k;
    if (((i64)b * b - d) % 4 != 0) continue;
    addFactor({((i64)b * b - d) / 4, (i64)b, 1});
  }
  for (int k = 2; k <= 10; k++) for (int rep = 0; rep < 25; rep++) {
    std::vector<i64> a(k + 1, 0); a[k] = 1; a[0] = -smoothVal(3) * (rnd_int(2) ? 1 : -1);
    addFactor(a);
  }
  for (int n = 3; n <= 30; n++) {
    GEN c = polcyclo(n, 0);
    if (degpol(c) > 10) continue;
    std::vector<i64> a(degpol(c) + 1);
    for (int i = 0; i <= degpol(c); i++) a[i] = (i64)itos(gel(c, 2 + i));
    addFactor(a);
  }
  avma = av0;
}
static bool mixedProduct(Poly &p)
{
  if (pool.empty()) return false;
  for (int tries = 0; tries < 200; tries++) {
    std::vector<i64> f(1, 1);
    int deg = 0, guard = 0;
    while (deg < 10 && guard++ < 50) {
      const Factor &F = pool[rnd_int((int)pool.size())];
      if (deg + F.deg > 10) continue;
      std::vector<i64> r;
      if (!polyMul(f, F.a, r)) { deg = -1; break; }
      f = r; deg += F.deg;
    }
    if (deg != 10) continue;
    fromVec(f, p);
    return p.c[0] == 1;
  }
  return false;
}
static bool symmetric(Poly &p)       // -g(x)g(-x);  h(x^5);  g(x^2) with g(0) = +-(S-smooth) k^2
{
  std::vector<i64> f(11, 0);
  i64 sgn = rnd_int(2) ? 1 : -1;
  int which = rnd_int(3);
  if (which == 0) {
    std::vector<i64> g(6), h(6), r;
    g[5] = 1;
    int B = 1 + rnd_int(9);
    for (int i = 0; i < 5; i++) g[i] = rnd_range(-B, B);
    for (int i = 0; i <= 5; i++) h[i] = (i % 2) ? -g[i] : g[i];
    if (!polyMul(g, h, r)) return false;
    for (int i = 0; i <= 10; i++) f[i] = -r[i];
  } else if (which == 1) {
    int b = rnd_int(41) - 20, k = 1 + rnd_int(12);
    i64 d = sgn * smoothVal(2) * k * k;
    if (((i64)b * b - d) % 4 != 0) return false;
    f[10] = 1; f[5] = b; f[0] = ((i64)b * b - d) / 4;
  } else {
    int k = 1 + rnd_int(12), B = 1 + rnd_int(30);
    f[10] = 1; f[0] = sgn * smoothVal(2) * k * k;
    for (int i = 1; i <= 4; i++) f[2 * i] = rnd_range(-B, B);
  }
  fromVec(f, p);
  return p.c[0] == 1;
}
static bool depressedEven(Poly &p)   // h(x^2), h a depressed quintic, constant +-(S-smooth) k^2
{
  int k = 1 + rnd_int(12), B = 1 + rnd_int(60);
  i64 e = (rnd_int(2) ? 1 : -1) * smoothVal(3) * k * k;
  std::vector<i64> f(11, 0);
  f[10] = 1; f[6] = rnd_range(-B, B); f[4] = rnd_range(-B, B); f[2] = rnd_range(-B, B); f[0] = e;
  fromVec(f, p);
  return true;
}
static bool tsBase(const Poly &q, Poly &p, int tsPercent)  // with probability tsPercent a Tschirnhaus transform of q
{
  pari_sp av = avma;
  GEN f = polyToGEN(q);
  bool ok = true;
  if (gequal0(poldisc0(f, -1))) ok = false;
  else if (rnd_int(100) < tsPercent) { GEN u = tsOwn(f, 9, 1 + rnd_int(3)); ok = u && genToPoly(u, p); }
  else p = q;
  avma = av;
  return ok;
}
static bool makeValid(Poly &p)
{
  Poly q;
  int w = rnd_int(seed0ok ? 5 : 4);
  switch (w) {
    case 0: return mixedProduct(q) && tsBase(q, p, 67);
    case 1: return symmetric(q) && tsBase(q, p, 67);
    case 2: case 3: return depressedEven(q) && tsBase(q, p, 30);
    default: { pari_sp av = avma; GEN u = tsOwn(seed0, 9, 1 + rnd_int(3)); bool ok = u && genToPoly(u, p); avma = av; return ok; }
  }
}

// ---------------------------------------------------------------- the target stratum: D with prescribed Legendre symbols
// chi_j(D) for the stages j = 1..NST as a bit vector, stage 1 in the most significant place (bit NST-j set: symbol -1).
static int jacobi(u64 a, u64 n)      // Jacobi symbol (a/n), n odd
{
  a %= n;
  int t = 1;
  while (a) {
    while (!(a & 1)) { a >>= 1; u64 r = n & 7; if (r == 3 || r == 5) t = -t; }
    std::swap(a, n);
    if ((a & 3) == 3 && (n & 3) == 3) t = -t;
    a %= n;
  }
  return n == 1 ? t : 0;
}
static bool chiVec(u64 D, int upto, u32 &v)   // the first `upto` stages; false if some prime divides D
{
  v = 0;
  for (int j = 1; j <= upto; j++) {
    int s = jacobi(D % FPr[j - 1], FPr[j - 1]);
    if (s == 0) return false;
    if (s < 0) v |= 1u << (NST - j);
  }
  return true;
}
// D = D1 * D2 with D1, D2 < 2^24: the symbols multiply, so the vectors add (xor).  A table of 2^17 vectors of D1, sorted,
// and random D2: the D1 that complete the wanted prefix of J stages are a range of the table.
struct DEntry { i64 D; u32 v; };
static std::vector<std::vector<DEntry> > dpool;           // dpool[J]: D rejected first at stage J (J = 0: never)
static void buildDPool()
{
  const int K = 1 << 17;
  std::vector<std::pair<u32, u32> > tab;
  tab.reserve(K);
  while ((int)tab.size() < K) { u32 d = 3 + (u32)(rnd64() % ((1u << 24) - 3)); u32 v; if (chiVec(d, NST, v)) tab.push_back(std::make_pair(v, d)); }
  std::sort(tab.begin(), tab.end());
  dpool.assign(NST + 1, std::vector<DEntry>());
  for (int J = 0; J <= NST; J++) {
    int L = J ? J : NST, want = J >= 24 ? 8 : 24;
    u32 tp = J ? 1u : 0u;
    for (i64 tries = 0; (int)dpool[J].size() < want && tries < 3000000; tries++) {
      u32 d2 = 3 + (u32)(rnd64() % ((1u << 24) - 3)), v2;
      if (!chiVec(d2, NST, v2)) continue;
      u32 need = tp ^ (v2 >> (NST - L));
      u32 lo = need << (NST - L), hi = (need + 1) << (NST - L);
      std::vector<std::pair<u32, u32> >::iterator a = std::lower_bound(tab.begin(), tab.end(), std::make_pair(lo, 0u));
      std::vector<std::pair<u32, u32> >::iterator b = std::lower_bound(a, tab.end(), std::make_pair(hi, 0u));
      if (a == b) continue;
      std::vector<std::pair<u32, u32> >::iterator e = a + (i64)(rnd64() % (u64)(b - a));
      DEntry de; de.D = (i64)e->second * (i64)d2; de.v = e->first ^ v2;
      dpool[J].push_back(de);
    }
  }
}
// Disc(g(x^2)) = (-4)^5 g(0) Disc(g)^2;  Disc(x^10 + b x^5 + c) = +-5^10 c^4 (b^2 - 4c)^5: the class of g(0), of b^2 - 4c
static bool makeTarget(Poly &p, int &J, u32 &v)
{
  J = rnd_int(10) == 0 ? 0 : 1 + rnd_int(NST);
  if (dpool[J].empty()) return false;
  const DEntry &de = dpool[J][rnd_int((int)dpool[J].size())];
  i64 D = (rnd64() & 1) ? de.D : -de.D;
  v = de.v;
  std::vector<i64> f(11, 0);
  f[10] = 1;
  i64 B = (i64)1 << (1 + rnd_int(12));
  int fam = rnd_int(10) < 6 ? 0 : rnd_int(2) + 1;   // mostly g(x^2), whose chains are mostly generic
  if (fam < 2) {                     // g(x^2), fam 1 depressed (no x^8): constant D
    for (int i = 1; i <= 4; i++) f[2 * i] = rnd_range(-B, B);
    if (fam == 1) f[8] = 0;
    f[0] = D;
  } else {                           // x^10 + b x^5 + c, b^2 - 4c = D or 4D
    i64 r = ((D % 4) + 4) % 4, Dm = (r == 2 || r == 3) ? 4 * D : D;
    i64 b = rnd_range(-B, B);
    if ((b & 1) != (Dm & 1)) b += 1;
    f[5] = b; f[0] = (b * b - Dm) / 4;
  }
  fromVec(f, p);
  if (rnd_unif() < 0.25) {
    pari_sp av = avma;
    GEN u = tsOwn(polyToGEN(p), 9, 1 + rnd_int(3));
    Poly q;
    if (u && genToPoly(u, q)) p = q;
    avma = av;
  }
  scrambleOrKeep(p, 0.7);
  return true;
}

// ---------------------------------------------------------------- residue polynomials mod q and the CRT strata
typedef std::array<u64, 11> Res;     // monic degree 10, r[i] = coefficient of x^i
static void resRandom(u64 q, Res &g) { for (int i = 0; i < 10; i++) g[i] = rnd64() % q; g[10] = 1; }
static void resMul(const std::vector<u64> &a, const std::vector<u64> &b, std::vector<u64> &r, u64 q)
{
  r.assign(a.size() + b.size() - 1, 0);
  for (size_t i = 0; i < a.size(); i++) for (size_t j = 0; j < b.size(); j++) r[i + j] = (r[i + j] + a[i] * b[j]) % q;
}
static void resShift(Res &g, u64 q, u64 t)                // g(x + t)
{
  for (int i = 0; i < 10; i++) for (int j = 9; j >= i; j--) g[j] = (g[j] + t * g[j + 1]) % q;
}
static void resWithRoots(u64 q, Res &g, int mult)          // (x - r)^mult * (x - s)^mult' * h: a repeated root mod q
{
  std::vector<u64> f(1, 1), r;
  int deg = 0, nsq = mult == 2 && rnd_int(2) ? 2 : 1;
  for (int k = 0; k < nsq; k++) {
    u64 root = rnd64() % q;
    for (int m = 0; m < mult; m++) { resMul(f, std::vector<u64>{(q - root) % q, 1}, r, q); f = r; deg++; }
  }
  std::vector<u64> h(11 - deg);
  for (int i = 0; i < 10 - deg; i++) h[i] = rnd64() % q;
  h[10 - deg] = 1;
  resMul(f, h, r, q);
  for (int i = 0; i <= 10; i++) g[i] = r[i];
}
// a degenerate (non-generic) chain mod q: trinomials x^10 + a x^k + b (first remainder of degree k), depressed
// h(x^2), two double roots or a triple root (gcd(f, f') of degree 2), then x -> x + t; or, for small q, a random
// polynomial that happens to be degenerate
static bool resDegenerate(u64 q, Res &g)
{
  for (int attempt = 0; attempt < 1000; attempt++) {
    int kind = rnd_int(q < 16384 ? 5 : 4);
    g.fill(0); g[10] = 1;
    if (kind == 0) { int k = 1 + rnd_int(7); g[k] = 1 + rnd64() % (q - 1); g[0] = rnd64() % q; }
    else if (kind == 1) { for (int i = 0; i <= 6; i += 2) g[i] = rnd64() % q; }
    else if (kind == 2) resWithRoots(q, g, 2);
    else if (kind == 3) resWithRoots(q, g, 3);
    else { for (int t = 0; t < 4000; t++) { resRandom(q, g); if (!oracleRes(g.data(), q).generic) return true; } continue; }
    resShift(g, q, rnd64() % q);
    if (!oracleRes(g.data(), q).generic) return true;
  }
  return false;
}
static bool resPass(u64 q, Res &g, double pZero)          // symbol 0 or +1 mod q
{
  if (rnd_unif() < pZero) { resWithRoots(q, g, 2); return true; }
  for (int t = 0; t < 1000; t++) { resRandom(q, g); if (oracleRes(g.data(), q).chi >= 0) return true; }
  return false;
}
// the coefficients from their residues (Garner): x[k] in [0, M); then crtPick adds to each a random multiple of M that
// keeps it in int64 (wide), or 0 or -M
struct Crt { i128 M; i128 x[10]; };
static void crtCombine(const std::vector<u64> &mods, const std::vector<Res> &res, Crt &c)
{
  c.M = 1;
  for (size_t i = 0; i < mods.size(); i++) c.M *= mods[i];
  for (int k = 0; k <= 9; k++) {
    i128 x = 0, m = 1;
    for (size_t i = 0; i < mods.size(); i++) {
      u64 q = mods[i], xm = (u64)(x % (i128)q), d = (res[i][k] + q - xm) % q;
      d = d * invm((u64)(m % (i128)q), q) % q;
      x += m * (i128)d; m *= q;
    }
    c.x[k] = x;
  }
}
static void crtPick(const Crt &c, Poly &p, bool wide)
{
  const i128 lo64 = -((i128)1 << 63), hi64 = ((i128)1 << 63) - 1;
  for (int k = 0; k <= 9; k++) {
    i128 x = c.x[k], t = 0;
    i128 tlo = -((x - lo64) / c.M), thi = (hi64 - x) / c.M;  // 0 <= x < M: the smallest and largest t with x + t*M in int64
    if (wide) t = tlo + (i128)(rnd64() % (u64)(thi - tlo + 1));
    else t = rnd_int(2) ? 0 : -1;
    p.c[10 - k] = (i64)(x + t * c.M);
  }
  p.c[0] = 1;
}

// ---------------------------------------------------------------- strata
enum St { S_RAND, S_SMALL, S_EVEN, S_VALID, S_TARGET, S_CRT, S_CRTDEG, S_CRTDEGL, S_REP, S_REAL, S_CAND, S_N };
static const char *stName[S_N] = { "rand", "small", "even", "valid", "target", "crt", "crtdeg", "crtdegL", "rep", "real", "cand" };
static std::vector<Poly> realSeeds, cands;
static std::vector<int> candFlag;    // the verdict recorded with a candidate (-1: none)
static size_t candPos = 0;
// What a construction guarantees, for the consistency check of the oracle:
//   1  kept at every stage (Disc = +-(S-smooth) * square, or 0)
//   2  Disc in the square class of D: wherever Disc != 0 mod q, the symbol is that of D (vector v)
struct Item { Poly p; int st, kind, J, cand; u32 v; };

static bool makeRand(Poly &p)
{
  static const i64 ext[] = { (i64)(-9223372036854775807ll - 1), 9223372036854775807ll, 0, 1, -1, 1ll << 21, -(1ll << 21), (1ll << 21) - 1,
                             1ll << 42, -(1ll << 42), (1ll << 42) - 1, -(1ll << 42) - 1, 1ll << 62, -(1ll << 62) };
  p.c[0] = 1;
  for (int k = 1; k <= 10; k++) p.c[k] = rnd_int(32) == 0 ? ext[rnd_int(14)] : rnd_bits(rnd_int(64));
  return true;
}
static bool makeSmall(Poly &p)
{
  int B = 1 + rnd_int(16);
  p.c[0] = 1;
  for (int k = 1; k <= 10; k++) p.c[k] = rnd_unif() < 0.3 ? 0 : rnd_bits(B);
  return true;
}
static bool makeEven(Poly &p)
{
  int B = 1 + rnd_int(20);
  p.c[0] = 1;
  for (int k = 1; k <= 10; k++) p.c[k] = (k % 2) ? 0 : rnd_bits(B);
  return true;
}
static bool makeRep(Poly &p)         // (x - a)^2 h or (x^2 + u x + w)^2 h over Z
{
  std::vector<i64> s, h, f;
  if (rnd_int(2)) { i64 a = rnd_range(-20, 20); s = {a * a, -2 * a, 1}; }
  else { i64 u = rnd_range(-9, 9), w = rnd_range(-30, 30); s = {w * w, 2 * u * w, u * u + 2 * w, 2 * u, 1}; }
  int dh = 10 - (int)(s.size() - 1), B = 1 + rnd_int(8);
  h.assign(dh + 1, 0);
  for (int i = 0; i < dh; i++) h[i] = rnd_bits(B);
  h[dh] = 1;
  if (!polyMul(s, h, f)) return false;
  fromVec(f, p);
  scrambleOrKeep(p, 0.7);
  return true;
}
// The CRT strata use as many primes as keep the modulus below 2^62 (four of the 30 filter primes, each below 2^13),
// so that every residue class has representatives in int64.
static const i128 CRTMAX = (i128)1 << 62;
static bool makeCrt(Poly &p, bool degen)
{
  std::vector<u64> mods; std::vector<Res> res;
  int nuse = 0;
  for (i128 M = 1; nuse < NS && M * FPr[nuse] < CRTMAX; nuse++) M *= FPr[nuse];
  if (nuse < 1) return false;
  int i0 = degen ? rnd_int(nuse) : -1;
  for (int i = 0; i < nuse; i++) {
    Res g;
    bool ok = (i == i0) ? resDegenerate(FPr[i], g) : resPass(FPr[i], g, 0.1);
    if (!ok) return false;
    mods.push_back(FPr[i]); res.push_back(g);
  }
  Crt c;
  crtCombine(mods, res, c);
  crtPick(c, p, rnd_unif() < 0.7);
  return true;
}
// The large stages that crtdegL can plant at: nsm small primes fit into the modulus next to a large one; the other NS - nsm
// small stages and the large stages before the planted one are passed by sampling, each with probability about 1/2, and
// the planted stage L is kept to at most 2^12 expected tries.  0: none (the stratum is then off).
static int crtdegLMax()
{
  int NL = NST - NS;
  if (NL < 1 || NS < 1) return 0;
  int nsm = 0;
  for (i128 M = FPr[NS]; nsm < std::min(3, NS) && M * FPr[nsm] < CRTMAX; nsm++) M *= FPr[nsm];
  return std::max(0, std::min(std::min(6, NL), 13 - (NS - nsm)));
}
static bool makeCrtDegL(Poly &p)
{
  int Lmax = crtdegLMax();
  if (Lmax < 1) return false;
  int L = 1 + rnd_int(Lmax), target = NS + L;
  std::vector<int> small(NS);
  std::iota(small.begin(), small.end(), 0);
  for (int i = NS - 1; i > 0; i--) std::swap(small[i], small[rnd_int(i + 1)]);
  std::vector<u64> mods; std::vector<Res> res;
  i128 M = FPr[target - 1];
  for (int i = 0; i < std::min(3, NS) && M * FPr[small[i]] < CRTMAX; i++) {
    Res g; if (!resPass(FPr[small[i]], g, 0.0)) return false;
    mods.push_back(FPr[small[i]]); res.push_back(g); M *= FPr[small[i]];
  }
  { Res g; if (!resDegenerate(FPr[target - 1], g)) return false; mods.push_back(FPr[target - 1]); res.push_back(g); }
  Crt c;
  crtCombine(mods, res, c);
  const int maxTries = 16 << (NS - (int)mods.size() + L);
  for (int tries = 0; tries < maxTries; tries++) {
    crtPick(c, p, true);
    OracleOut o = oracle(p);
    if (o.stage == 0 || o.stage >= target) return true;
  }
  return false;
}

// ---------------------------------------------------------------- coverage, from the oracle
struct Cov {
  i64 polys = 0, rejS = 0, rejL = 0, kept = 0, reachL = 0;
  i64 visS = 0, visL = 0, ngS = 0, ngL = 0, rejNgS = 0, rejNgL = 0, zS = 0, zL = 0;
  std::set<u16> patS, patL;
  i64 hist[6] = {0, 0, 0, 0, 0, 0};  // rejected at stages 1-6, 7-12, 13-18, 19-24, 25-30, kept
  i64 bits[8] = {0, 0, 0, 0, 0, 0, 0, 0};
};
static int bitLen(const Poly &p)
{
  int m = 0;
  for (int i = 1; i <= 10; i++) { u64 v = p.c[i] < 0 ? (u64)(-(p.c[i] + 1)) + 1 : (u64)p.c[i]; int b = 0; while (v) { b++; v >>= 1; } if (b > m) m = b; }
  return m;
}
static void cover(Cov &c, const Poly &p, const OracleOut &o)
{
  c.polys++;
  int last = o.stage ? o.stage : NST;
  for (int j = 1; j <= last; j++) {
    bool small = j <= NS, ng = (o.nongen >> (j - 1)) & 1, z = (o.zero >> (j - 1)) & 1;
    (small ? c.visS : c.visL)++;
    if (ng) { (small ? c.ngS : c.ngL)++; (small ? c.patS : c.patL).insert(o.pat[j - 1]); if (j == o.stage) (small ? c.rejNgS : c.rejNgL)++; }
    if (z) (small ? c.zS : c.zL)++;
  }
  if (o.stage == 0) c.kept++; else if (o.stage <= NS) c.rejS++; else c.rejL++;
  if (o.stage == 0 || o.stage > NS) c.reachL++;
  c.hist[o.stage ? std::min(4, (o.stage - 1) / 6) : 5]++;
  int b = bitLen(p);
  c.bits[b <= 8 ? 0 : std::min(7, (b - 1) / 8)]++;
}

// ---------------------------------------------------------------- the filter builds under test
struct Build {
  const char *name;
  void (*init)(int, const int *);
  void (*filter)(const i64 *, int, signed char *);
  void (*stages)(int, int *, long long *, long long *, int *, int);
  bool hook;
  i64 falseRej = 0, missed = 0, wrongStage = 0, route = 0, meta = 0, stageFalse = 0;
  i64 failSt[S_N] = {0};             // false, missed, wrong stage and transformed-copy failures by stratum
};
static i64 shownFail = 0;
static void showPoly(const char *what, const char *build, const Item &it, int bs, int got, int want)
{
  if (shownFail >= 12) return;
  shownFail++;
  printf("  %s by %s (%s, block of %d, filter stage %d, oracle stage %d): x^10", what, build, stName[it.st], bs, got, want);
  for (int k = 1; k <= 10; k++) printf(" %+lld", (i64)it.p.c[k]);
  printf("\n");
}
// the counters of one filter call against the oracle: stage j receives the lanes that reach it, chain_block (when used)
// retries exactly the non-generic ones, and keeps the lanes that pass
static long long hookTot[7 * MAXST];
static i64 routeCheck(const std::vector<OracleOut> &orc, const int *ids, int n, const long long *cnt)
{
  i64 bad = 0;
  for (int j = 1; j <= NST; j++) {
    i64 in = 0, ng = 0, out = 0;
    for (int i = 0; i < n; i++) {
      const OracleOut &o = orc[ids[i]];
      if (o.stage != 0 && o.stage < j) continue;
      in++;
      if ((o.nongen >> (j - 1)) & 1) ng++;
      if (o.stage != j) out++;
    }
    const long long *c = cnt + 7 * (j - 1);
    bool ok = c[0] == in && c[1] + c[2] == c[0] && (c[1] == 0 || c[2] == 0) && c[3] == (c[1] > 0 ? ng : 0) && c[4] == out;
    if (!ok) bad++;
  }
  return bad;
}

static void usage(const char *prog)
{
  printf("Usage: %s [N] [seed] [options]\n", prog);
  printf("  N              blocks' worth of 10,000 polynomials (default 100); -n TOTAL gives the number of polynomials\n");
  printf("  seed           random seed (default: the time; it is printed); also -seed X\n");
  printf("  -g LIST        strata with weights, e.g. target=3,crt,rand (default: all; weight 1 if omitted)\n");
  printf("                 rand small even valid target crt crtdeg crtdegL rep real cand\n");
  printf("  -S 2,5         the primes divided out by the upstream test (default 2,5; at most 3)\n");
  printf("  -block B       always call the filter with B polynomials (default: sizes vary)\n");
  printf("  -real FILE     the application's passing polynomials (bb1..bb10 or 11 integers per line)\n");
  printf("  -cand FILE     the application's candidates (11 integers per line, optionally the recorded verdict 0/1 as a 12th)\n");
  printf("  -xcheck F      fraction of the polynomials whose discriminant is also computed by PARI (default 0.01; the first 100 of each stratum always)\n");
  printf("  -nometa        no transformed second run\n");
  printf("  -rm MODE       up | down | zero: run every filter call under another rounding mode\n");
  printf("  -stages        print the table per stage (the counters of the test hook)\n");
  printf("  -oldrule       without a test hook, take the primes of the rule before 2026-10-29 (6 below 2^10, 24 below 2^25), as\n"
         "                 the floating-point filter of 327b04f has; a filter linked with REF must have the primes of the first build\n");
}

static int run(int argc, char **argv)
{
  i64 N = 100 * 10000; int fixedBlock = 0; u64 seed = (u64)time(NULL);
  double xfrac = 0.01; bool doMeta = true, showStages = false;
  std::string realFile, candFile, Sarg = "2,5", gArg = "all";
  int npos = 0;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    auto nxt = [&]() -> std::string { return i + 1 < argc ? std::string(argv[++i]) : std::string(); };
    if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
    else if (a == "-g") gArg = nxt();
    else if (a == "-n") { N = atoll(nxt().c_str()); if (npos == 0) npos = 1; }
    else if (a == "-S") Sarg = nxt();
    else if (a == "-seed") seed = (u64)atoll(nxt().c_str());
    else if (a == "-block") fixedBlock = atoi(nxt().c_str());
    else if (a == "-real") realFile = nxt();
    else if (a == "-cand") candFile = nxt();
    else if (a == "-xcheck") xfrac = atof(nxt().c_str());
    else if (a == "-nometa") doMeta = false;
    else if (a == "-stages") showStages = true;
    else if (a == "-oldrule") oldRule = true;
    else if (a == "-rm") { std::string v = nxt(); hostile = v == "up" ? 1 : v == "down" ? 2 : v == "zero" ? 3 : 0; if (!hostile) { usage(argv[0]); return 2; } }
    else if (a[0] != '-' && npos == 0 && atoll(a.c_str()) > 0) { N = atoll(a.c_str()) * 10000; npos++; }
    else if (a[0] != '-' && npos == 1) { seed = (u64)atoll(a.c_str()); npos++; }
    else { usage(argv[0]); return 2; }
  }
  rs = seed;

  pari_init(200000000, 10000);
  setrand(utoi((u64)seed));
  LIM = shifti(gen_1, 63); NLIM = negi(LIM);
  { size_t pos = 0; while (pos < Sarg.size()) { size_t e = Sarg.find(',', pos); if (e == std::string::npos) e = Sarg.size(); Sset.push_back(atoi(Sarg.substr(pos, e - pos).c_str())); pos = e + 1; } }
  if (Sset.empty() || Sset.size() > 3) { printf("-S: one to three primes\n"); return 2; }
  seed0 = gclone(mkpoln(11, gen_1, gen_0, gen_0, gen_0, gen_0, stoi(-12), gen_0, gen_0, gen_0, gen_0, stoi(-4)));
  seed0ok = passesUpstream(poldisc0(seed0, -1));
  buildPool();
  // polynomials from files: 10 integers (bb1..bb10, leading 1 implied) or 11 (leading 1 first), and for -cand an optional verdict
  auto readPolys = [](const std::string &file, std::vector<Poly> &out, std::vector<int> *flags) -> bool {
    FILE *f = fopen(file.c_str(), "r"); char line[4096];
    if (!f) { printf("cannot read %s\n", file.c_str()); return false; }
    while (fgets(line, sizeof line, f)) {
      if (line[0] == '#') continue;
      i64 v[13]; int n = 0; char *s = line;
      while (n < 13) { char *e; i64 x = strtoll(s, &e, 10); if (e == s) break; v[n++] = x; s = e; }
      Poly q; q.c[0] = 1; int fl = -1;
      if (n == 10) for (int k = 0; k < 10; k++) q.c[1 + k] = v[k];
      else if (n >= 11 && v[0] == 1) { for (int k = 1; k <= 10; k++) q.c[k] = v[k]; if (n >= 12) fl = (int)v[11]; }
      else continue;
      out.push_back(q);
      if (flags) flags->push_back(fl);
    }
    fclose(f);
    return true;
  };
  if (!realFile.empty() && !readPolys(realFile, realSeeds, NULL)) return 2;
  if (!candFile.empty() && !readPolys(candFile, cands, &candFlag)) return 2;

  // builds
  std::vector<Build> builds;
  int pS[3] = { 0, 0, 0 }; for (size_t k = 0; k < Sset.size(); k++) pS[k] = Sset[k];
#ifdef NF_PDF_DISPATCH
  __builtin_cpu_init();
  const char *fb = getenv("NF_FORCE_BASE");
  if (!(fb && *fb == '1') && __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")) { Build b; b.name = "AVX2/FMA"; b.init = pdf_init_v3; b.filter = pdf_filter_v3; b.stages = pdf_stages_v3; builds.push_back(b); }
  { Build b; b.name = "baseline"; b.init = pdf_init_base; b.filter = pdf_filter_base; b.stages = pdf_stages_base; builds.push_back(b); }
#else
  { Build b; b.name = "single build"; b.init = pdf_init; b.filter = pdf_filter; b.stages = pdf_stages; builds.push_back(b); }
#endif
  std::vector<int> stbuf(10240 + 64);
  long long cntbuf[7 * MAXST], prbuf[MAXST];
  std::vector<u64> specPr = specPrimes();
  std::vector<u64> hookPr;
  bool primesDiffer = false, primesUnusable = false;
  for (size_t b = 0; b < builds.size(); b++) {
    builds[b].init((int)Sset.size(), pS);
    int nst = 0;
    builds[b].stages(0, stbuf.data(), cntbuf, prbuf, &nst, 1);
    builds[b].hook = nst > 0;
    if (builds[b].hook) {
      std::vector<u64> pr(prbuf, prbuf + nst);
      if (hookPr.empty()) hookPr = pr; else if (pr != hookPr) primesDiffer = true;
    }
  }
  pdf_init_ref((int)Sset.size(), pS);
  bool haveRef = !refAbsent;
  FPr = hookPr.empty() ? specPr : hookPr;
  NST = (int)FPr.size();
  NS = 0;
  for (int j = 0; j < NST; j++) if (FPr[j] < 65536) NS++;
  for (int j = 0; j < NST; j++) FF.push_back(Fq(FPr[j]));
  for (int j = 0; j < NST; j++) if (!usablePrime(FPr[j])) primesUnusable = true;
  bool asSpec = FPr == specPr;
  buildDPool();

  // strata and weights
  std::vector<double> wt(S_N, 0.0);
  if (gArg == "all") {
    double d[S_N] = { 1, 1, 1, 2, 3, 2, 2, 2, 1, 1, 2 };
    for (int s = 0; s < S_N; s++) wt[s] = d[s];
    if (realSeeds.empty()) wt[S_REAL] = 0;
    if (cands.empty()) wt[S_CAND] = 0;
  } else {
    size_t pos = 0;
    while (pos < gArg.size()) {
      size_t e = gArg.find(',', pos); if (e == std::string::npos) e = gArg.size();
      std::string tok = gArg.substr(pos, e - pos), nm = tok; double w = 1;
      size_t eq = tok.find('='); if (eq != std::string::npos) { nm = tok.substr(0, eq); w = atof(tok.substr(eq + 1).c_str()); }
      int found = -1; for (int s = 0; s < S_N; s++) if (nm == stName[s]) found = s;
      if (found < 0) { printf("-g: unknown stratum %s\n", nm.c_str()); return 2; }
      wt[found] = w; pos = e + 1;
    }
    if (wt[S_REAL] > 0 && realSeeds.empty()) { printf("stratum real needs -real FILE\n"); return 2; }
    if (wt[S_CAND] > 0 && cands.empty()) { printf("stratum cand needs -cand FILE\n"); return 2; }
  }
  if (wt[S_CRTDEGL] > 0 && crtdegLMax() < 1) {
    if (gArg != "all") { printf("stratum crtdegL: with %d small stages no sampled polynomial reaches a large one\n", NS); return 2; }
    wt[S_CRTDEGL] = 0;
    printf("stratum crtdegL off: with %d small stages no sampled polynomial reaches a large one\n", NS);
  }
  double wsum = 0; for (int s = 0; s < S_N; s++) wsum += wt[s];
  if (wsum <= 0) { printf("no stratum\n"); return 2; }

  printf("pdfilter oracle test: S=%s, %lld polynomials, seed=%llu, builds:", Sarg.c_str(), (i64)N, (u64)seed);
  for (size_t b = 0; b < builds.size(); b++) printf(" %s%s", builds[b].name, builds[b].hook ? "" : " (no test hook)");
  if (haveRef) printf(" + reference filter");
  printf("%s\n", hostile ? "  [filters run under a non-default rounding mode]" : "");
  printf("strata:"); for (int s = 0; s < S_N; s++) if (wt[s] > 0) printf(" %s=%g", stName[s], wt[s]); printf("\n");
  printf("filter primes (%s): %d small %llu..%llu", hookPr.empty() ? "pdf_init's rule; no test hook" : "reported by the filter",
         NS, NS ? FPr[0] : 0ull, NS ? FPr[NS - 1] : 0ull);
  if (NST > NS) printf(", %d large %llu..%llu", NST - NS, FPr[NS], FPr[NST - 1]);
  printf("%s%s\n", asSpec ? ", as pdf_init's rule gives" : ", NOT the primes of pdf_init's rule", primesUnusable ? ", SOME NOT USABLE" : "");
  if (primesDiffer) printf("THE BUILDS REPORT DIFFERENT PRIMES\n");
  if (!realSeeds.empty()) printf("real polynomials: %zu\n", realSeeds.size());
  if (!cands.empty()) printf("candidates: %zu\n", cands.size());

  // ---- the main loop
  const int blockSizes[] = { 1, 7, 63, 64, 65, 127, 255, 1000, 4096, 10000, 10240 };
  std::vector<Cov> cov(S_N);
  std::vector<i64> xchecked(S_N, 0);
  i64 tested = 0, xcheckPolys = 0, oracleBug = 0, genBug = 0, classChecks = 0, invChecks = 0, genFail = 0;
  i64 refFalse = 0, refMissed = 0, refDiff = 0, candChecked = 0, candDiff = 0;
  std::vector<Item> queue, blk;
  std::vector<OracleOut> orc;
  std::vector<i64> buf, buf2;
  std::vector<signed char> flag, flag2;
  std::vector<std::vector<signed char> > flagB(builds.size());
  std::vector<std::vector<int> > stageB(builds.size());
  memset(hookTot, 0, sizeof hookTot);
  auto T0 = std::chrono::steady_clock::now();
  pari_sp ltop = avma;
  i64 blockNo = 0;

  // Half of the blocks draw from all strata, half from a single one (as the application's blocks come from one work unit):
  // the filter's routing (prefgen, decided on every 64th call) then sees both mostly generic and mostly degenerate blocks.
  int blockStratum = -1;
  auto drawStratum = [&]() -> int { double r = rnd_unif() * wsum; for (int s = 0; s < S_N; s++) { if (r < wt[s]) return s; r -= wt[s]; } return S_N - 1; };
  auto pickStratum = [&]() -> int { return blockStratum >= 0 ? blockStratum : drawStratum(); };
  auto nextItem = [&](Item &it) -> bool {
    for (;;) {
      avma = ltop;
      if (!queue.empty()) { it = queue.back(); queue.pop_back(); return true; }
      int s = pickStratum();
      Item b; b.st = s; b.kind = 0; b.J = -1; b.v = 0; b.cand = -1;
      bool ok = false; int reuse = 1;
      switch (s) {
        case S_RAND: ok = makeRand(b.p); break;
        case S_SMALL: ok = makeSmall(b.p); break;
        case S_EVEN: ok = makeEven(b.p); break;
        case S_VALID: ok = makeValid(b.p); b.kind = 1; reuse = 4; break;
        case S_TARGET: ok = makeTarget(b.p, b.J, b.v); b.kind = 2; reuse = 2; break;
        case S_CRT: ok = makeCrt(b.p, false); break;
        case S_CRTDEG: ok = makeCrt(b.p, true); break;
        case S_CRTDEGL: ok = makeCrtDegL(b.p); break;
        case S_REP: ok = makeRep(b.p); b.kind = 1; break;
        case S_REAL: b.p = realSeeds[rnd_int((int)realSeeds.size())]; b.kind = 1; ok = true; reuse = 2; break;
        case S_CAND: b.cand = (int)(candPos % cands.size()); b.p = cands[b.cand]; candPos++; ok = true; break;
      }
      if (!ok) { if (++genFail > 1000 + N) { printf("a generator keeps failing (%s)\n", stName[s]); exit(3); } continue; }
      for (int v = 0; v < reuse; v++) {
        Item c = b;
        if (v > 0 || s == S_VALID || s == S_REAL) scrambleOrKeep(c.p, v > 0 ? 0.95 : 0.8);
        if (v > 0) c.cand = -1;
        queue.push_back(c);
      }
    }
  };

  while (tested < N) {
    int bs = fixedBlock ? fixedBlock : blockSizes[rnd_int(11)];
    if ((i64)bs > N - tested) bs = (int)(N - tested);
    blockStratum = rnd_int(2) ? drawStratum() : -1;
    queue.clear();
    blk.resize(bs); orc.resize(bs);
    if ((int)stbuf.size() < bs) stbuf.resize(bs);
    buf.assign((size_t)bs * 11, 0);
    for (int i = 0; i < bs; i++) {
      nextItem(blk[i]);
      const Item &it = blk[i];
      bool xc = xchecked[it.st] < 100 || rnd_unif() < xfrac;
      orc[i] = oracle(it.p, it.kind || xc ? NST : 0);
      const OracleOut &o = orc[i];
      // the oracle against what the construction guarantees
      if (it.kind) {
        classChecks++;
        bool bad = false;
        if (it.kind == 1) bad = o.stage != 0;
        else for (int j = 1; j <= o.reached; j++) {
          int want = (it.v >> (NST - j)) & 1;
          bool z = (o.zero >> (j - 1)) & 1, ng = (o.neg >> (j - 1)) & 1;
          if (!z && (int)ng != want) bad = true;
        }
        if (bad) {
          pari_sp av = avma;
          bool genIsWrong = it.kind == 1 && !gequal0(poldisc0(polyToGEN(it.p), -1)) && !passesUpstream(poldisc0(polyToGEN(it.p), -1));
          avma = av;
          (genIsWrong ? genBug : oracleBug)++;
          if (genBug + oracleBug <= 5) { printf("%s (%s): the oracle gives stage %d:", genIsWrong ? "GENERATOR BUG" : "ORACLE INCONSISTENT WITH THE CONSTRUCTION", stName[it.st], o.stage); for (int k = 1; k <= 10; k++) printf(" %lld", (i64)it.p.c[k]); printf("\n"); }
        }
      }
      // the oracle against PARI: Disc mod q at every stage
      if (xc) {
        xchecked[it.st]++; xcheckPolys++;
        pari_sp av = avma;
        GEN D = poldisc0(polyToGEN(it.p), -1);
        for (int j = 1; j <= NST; j++) {
          u64 dq = umodiu(D, (u64)FPr[j - 1]);
          if (dq != oracleMod(it.p, FF[j - 1]).disc) { oracleBug++; if (oracleBug <= 5) { printf("ORACLE DIFFERS FROM PARI at q=%llu:", FPr[j - 1]); for (int k = 1; k <= 10; k++) printf(" %lld", (i64)it.p.c[k]); printf("\n"); } break; }
        }
        avma = av;
      }
      cover(cov[it.st], it.p, o);
      if (it.cand >= 0 && candFlag[it.cand] >= 0) { candChecked++; if (candFlag[it.cand] != (o.stage == 0 ? 1 : 0)) candDiff++; }
      for (int k = 0; k <= 10; k++) buf[(size_t)i * 11 + k] = it.p.c[k];
    }
    // every build, the whole block in one call
    std::vector<int> ids(bs);
    std::iota(ids.begin(), ids.end(), 0);
    for (size_t b = 0; b < builds.size(); b++) {
      Build &B = builds[b];
      flag.assign(bs, 99);
      roundOn();
      B.filter(buf.data(), bs, flag.data());
      roundOff();
      flagB[b] = flag;
      stageB[b].assign(bs, -1);
      if (B.hook) {
        int nst = 0;
        B.stages(bs, stbuf.data(), cntbuf, prbuf, &nst, 1);
        for (int i = 0; i < bs; i++) stageB[b][i] = stbuf[i];
        B.route += routeCheck(orc, ids.data(), bs, cntbuf);
        if (b == 0) for (int k = 0; k < 7 * MAXST; k++) hookTot[k] += cntbuf[k];
      }
      for (int i = 0; i < bs; i++) {
        int want = orc[i].stage, got = B.hook ? stageB[b][i] : (flag[i] ? 0 : -1);
        bool kept = flag[i] != 0;
        if (B.hook && (stageB[b][i] == 0) != kept) B.stageFalse++;           // the hook contradicts the verdict
        if (want == 0 && !kept) { B.falseRej++; B.failSt[blk[i].st]++; showPoly("FALSE REJECTION", B.name, blk[i], bs, got, want); }
        else if (want != 0 && kept) { B.missed++; B.failSt[blk[i].st]++; showPoly("MISSED REJECTION", B.name, blk[i], bs, got, want); }
        else if (want != 0 && B.hook && got != want) { B.wrongStage++; B.failSt[blk[i].st]++; showPoly("WRONG STAGE", B.name, blk[i], bs, got, want); }
      }
    }
    // the reference filter (differential)
    if (haveRef) {
      flag2.assign(bs, 99);
      roundOn();
      pdf_filter_ref(buf.data(), bs, flag2.data());
      roundOff();
      for (int i = 0; i < bs; i++) {
        bool kept = flag2[i] != 0, want = orc[i].stage == 0;
        if (want && !kept) refFalse++;
        if (!want && kept) refMissed++;
        if (kept != (flagB[0][i] != 0)) refDiff++;
      }
    }
    // transformed, permuted and re-cut: the same verdicts and stages (one build per block, alternating)
    if (doMeta) {
      size_t b = (size_t)(blockNo % (i64)builds.size());
      Build &B = builds[b];
      std::vector<int> perm(bs);
      std::iota(perm.begin(), perm.end(), 0);
      for (int i = bs - 1; i > 0; i--) std::swap(perm[i], perm[rnd_int(i + 1)]);
      std::vector<Poly> tp(bs);
      for (int i = 0; i < bs; i++) {
        tp[i] = blk[i].p;
        Poly q = tp[i];
        int kind = rnd_int(4);
        bool ok = kind == 0 ? affine(q, rnd_range(1, 8) * (rnd_int(2) ? 1 : -1), 1, false)
                : kind == 1 ? affine(q, 0, 1, true)
                : kind == 2 ? affine(q, 0, 2 + rnd_int(2), false)
                : affine(q, rnd_range(1, 8), 1, true);
        if (!ok) { q = tp[i]; ok = affine(q, rnd_int(2) ? 1 : -1, 1, false); }
        if (ok) tp[i] = q;
        if (rnd_int(20) == 0) {          // the oracle is invariant too (verdict, degree sequences, zeros)
          invChecks++;
          OracleOut o2 = oracle(tp[i], orc[i].reached);
          const OracleOut &o = orc[i];
          bool same = o2.stage == o.stage && o2.nongen == o.nongen && o2.zero == o.zero;
          for (int j = 0; same && j < o.reached; j++) same = o2.pat[j] == o.pat[j];
          if (!same) { oracleBug++; if (oracleBug <= 5) printf("ORACLE NOT INVARIANT under a transformation (%s)\n", stName[blk[i].st]); }
        }
      }
      int pos = 0;
      while (pos < bs) {
        int cs = blockSizes[rnd_int(11)];
        if (cs > bs - pos) cs = bs - pos;
        buf2.assign((size_t)cs * 11, 0);
        for (int i = 0; i < cs; i++) for (int k = 0; k <= 10; k++) buf2[(size_t)i * 11 + k] = tp[perm[pos + i]].c[k];
        flag2.assign(cs, 99);
        roundOn();
        B.filter(buf2.data(), cs, flag2.data());
        roundOff();
        if (B.hook) {
          int nst = 0;
          B.stages(cs, stbuf.data(), cntbuf, prbuf, &nst, 1);
          B.route += routeCheck(orc, perm.data() + pos, cs, cntbuf);
        }
        for (int i = 0; i < cs; i++) {
          int o = perm[pos + i];
          bool differ = (flag2[i] != 0) != (flagB[b][o] != 0) || (B.hook && stbuf[i] != stageB[b][o]);
          if (differ) { B.meta++; B.failSt[blk[o].st]++; if (shownFail < 12) { shownFail++; printf("  TRANSFORMED COPY DIFFERS in %s (%s): verdict %d -> %d\n", B.name, stName[blk[o].st], flagB[b][o] != 0, flag2[i] != 0); } }
        }
        pos += cs;
      }
    }
    tested += bs;
    blockNo++;
  }

  // ---- report
  double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
  printf("---- coverage by stratum, from the oracle (%.1f s) ----\n", wall);
  printf("%-8s %9s %8s %8s %8s %8s | %9s %9s | %8s %8s | %7s %7s | %8s %8s | %5s %5s\n", "stratum", "polys", "rej S", "rej L", "kept", "reach L",
         "visits S", "visits L", "nongen S", "nongen L", "rejNG S", "rejNG L", "Disc0 S", "Disc0 L", "pat S", "pat L");
  Cov tot;
  for (int s = 0; s < S_N; s++) {
    const Cov &c = cov[s];
    if (!c.polys) continue;
    printf("%-8s %9lld %8lld %8lld %8lld %8lld | %9lld %9lld | %8lld %8lld | %7lld %7lld | %8lld %8lld | %5zu %5zu\n", stName[s], c.polys, c.rejS, c.rejL, c.kept, c.reachL,
           c.visS, c.visL, c.ngS, c.ngL, c.rejNgS, c.rejNgL, c.zS, c.zL, c.patS.size(), c.patL.size());
    tot.polys += c.polys; tot.rejS += c.rejS; tot.rejL += c.rejL; tot.kept += c.kept; tot.reachL += c.reachL; tot.visS += c.visS; tot.visL += c.visL;
    tot.ngS += c.ngS; tot.ngL += c.ngL; tot.rejNgS += c.rejNgS; tot.rejNgL += c.rejNgL; tot.zS += c.zS; tot.zL += c.zL;
    tot.patS.insert(c.patS.begin(), c.patS.end()); tot.patL.insert(c.patL.begin(), c.patL.end());
    for (int k = 0; k < 6; k++) tot.hist[k] += c.hist[k];
    for (int k = 0; k < 8; k++) tot.bits[k] += c.bits[k];
  }
  printf("%-8s %9lld %8lld %8lld %8lld %8lld | %9lld %9lld | %8lld %8lld | %7lld %7lld | %8lld %8lld | %5zu %5zu\n", "all", tot.polys, tot.rejS, tot.rejL, tot.kept, tot.reachL,
         tot.visS, tot.visL, tot.ngS, tot.ngL, tot.rejNgS, tot.rejNgL, tot.zS, tot.zL, tot.patS.size(), tot.patL.size());
  printf("(rej S / L: rejected at a small (1-%d) / large stage; reach L: kept or rejected at a large stage; visits: (polynomial, stage) pairs\n"
         " evaluated at small (S) / large (L) stages; nongen: non-generic Euclid chains among them; rejNG: rejections with a\n"
         " non-generic chain (by chain_general); Disc0: Disc = 0 mod q; pat: distinct degree sequences of the non-generic chains)\n", NS);
  printf("rejection stage, all strata: 1-6 %lld, 7-12 %lld, 13-18 %lld, 19-24 %lld, 25-30 %lld, kept %lld\n", tot.hist[0], tot.hist[1], tot.hist[2], tot.hist[3], tot.hist[4], tot.hist[5]);
  {
    static const char *nm[] = { "<=8", "9-16", "17-24", "25-32", "33-40", "41-48", "49-56", "57-63" };
    printf("largest coefficient, bits:"); for (int k = 0; k < 8; k++) printf(" %s:%lld", nm[k], tot.bits[k]); printf("\n");
  }
  if (cov[S_TARGET].polys) {
    const Cov &c = cov[S_TARGET];
    printf("target stratum, rejection stage: 1-6 %lld, 7-12 %lld, 13-18 %lld, 19-24 %lld, 25-30 %lld, kept %lld\n", c.hist[0], c.hist[1], c.hist[2], c.hist[3], c.hist[4], c.hist[5]);
  }
  if (builds[0].hook) {
    printf("---- per stage, from the test hook of the %s build (first runs) ----\n", builds[0].name);
    printf("%5s %9s %10s %10s %11s %13s %9s %10s %11s %10s\n", "stage", "prime", "lanes in", "rejected", "chain_block", "chain_general", "retried", "kept", "lockstep it", "masked it");
    for (int j = 1; j <= NST; j++) {
      const long long *c = hookTot + 7 * (j - 1);
      if (!showStages && j > NS + 2 && j < NST) continue;
      printf("%5d %9llu %10lld %10lld %11lld %13lld %9lld %10lld %11lld %10lld\n", j, FPr[j - 1], c[0], c[0] - c[4], c[1], c[2], c[3], c[4], c[5], c[6]);
    }
    if (!showStages) printf("(stages %d..%d left out; -stages prints them)\n", NS + 3, NST - 1);
  }
  printf("---- results ----\n");
  // The totals that fault_checks.sh reads are those of the first build (AVX2/FMA when the CPU has it, as the application would pick;
  // the baseline build with NF_FORCE_BASE=1), so that its two columns are one build each; every build counts for the exit status.
  i64 fails = 0, first = 0;
  for (size_t b = 0; b < builds.size(); b++) {
    const Build &B = builds[b];
    printf("[%s] false rejections %lld, missed rejections %lld", B.name, B.falseRej, B.missed);
    if (B.hook) printf(", wrong stage %lld, routing mismatches %lld (block-stages), hook contradicts verdict %lld", B.wrongStage, B.route, B.stageFalse);
    if (doMeta) printf(", transformed copies differing %lld", B.meta);
    printf("\n");
    i64 fs = 0; for (int s = 0; s < S_N; s++) fs += B.failSt[s];
    if (fs) { printf("    by stratum:"); for (int s = 0; s < S_N; s++) if (B.failSt[s]) printf(" %s=%lld", stName[s], B.failSt[s]); printf("\n"); }
    i64 fb = B.falseRej + B.missed + B.wrongStage + B.route + B.meta + B.stageFalse;
    fails += fb;
    if (b == 0) first = fb;
  }
  if (haveRef) {
    printf("[reference] false rejections %lld, missed rejections %lld, verdicts differing from %s %lld\n", refFalse, refMissed, builds[0].name, refDiff);
    fails += refDiff; first += refDiff;
  }
  if (candChecked) printf("candidates with a recorded verdict: %lld, verdicts differing from the oracle: %lld\n", candChecked, candDiff);
  fails += candDiff; first += candDiff;
  printf("oracle checked: against PARI's poldisc on %lld polynomials x %d primes; against the construction on %lld; invariance on %lld: %lld disagreements\n",
         xcheckPolys, NST, classChecks, invChecks, oracleBug);
  printf("Total polynomials tested = %lld\n", (i64)tested);
  printf("False rejections = %lld   (%s build)\n", builds[0].falseRej, builds[0].name);
  printf("Missed rejections = %lld   (%s build)\n", builds[0].missed, builds[0].name);
  printf("Total failures = %lld   (%s build)\n", first, builds[0].name);
  if (builds.size() > 1) printf("All builds: failures %lld\n", fails);
  if (!asSpec || primesUnusable || primesDiffer) printf("The filter's primes are %s\n", primesUnusable ? "NOT ALL USABLE" : primesDiffer ? "DIFFERENT IN THE BUILDS" : "not those of pdf_init's rule");
  if (oracleBug || genBug) { printf("HARNESS BUG: %lld oracle disagreements, %lld generator errors\n", oracleBug, genBug); return 2; }
  return (fails || primesUnusable || primesDiffer) ? 1 : 0;
}

}  // namespace pdt

int main(int argc, char **argv) { return pdt::run(argc, argv); }
