// test_pdfilter.cpp -- stress test of the pre-filter (pdfilter.f90) with polynomials that MUST be kept.
//
// The idea is Eric's.  In the real search so few polynomials pass the discriminant test that a filter which
// wrongly rejected some of them would go unnoticed.  So we manufacture polynomials that are known to pass --
// their discriminant is +-(S-smooth) * (a perfect square) -- and check that the filter keeps every one of them.
// A transformation that changes the polynomial but keeps its discriminant in the same square class does this:
//   * Tschirnhaus: f(x) -> the characteristic polynomial of g(alpha), alpha a root of f, g in Z[x]:
//     Disc(new) = det(M)^2 * Disc(f) with an integer matrix M, so Disc(new)/Disc(f) is an integer square;
//   * x -> x + t (Disc unchanged), x -> -x (unchanged), x -> c*x in monic form c^10 f(x/c) (Disc * c^90).
// Only polynomials with such a discriminant are used here (no polynomial that should be rejected).
//
// ---------------------------------------------------------------------------------------------------------
// HOW THIS VERSION IMPROVES ON ERIC'S (his original is in git history; `-g orig` still produces its polynomials)
//
//  1. Diversity.  PARI's tschirnhaus() draws only three small numbers (g = a*x^2 + b*x + c, a in {1,2,3}, b and c in
//     {-4..3}): at most 192 different transformations.  So the 10,000*N polynomials of the original are at most 192
//     DISTINCT polynomials, repeated; all have coefficients below 2^28, all lie in one number field, and modulo the 30
//     filter primes they show 3 of the 22 possible factorisation patterns (a square discriminant allows the 22
//     partitions of 10 into an even number of parts).  A failure that hits one polynomial in 10,000 is invisible.
//     Here g has random degree 1..9 and coefficients up to +-3, followed by exact shift/scale/flip: millions of distinct
//     polynomials, coefficients from a few bits up to 63 bits.
//  2. Seeds.  Not one seed polynomial but families (all with a discriminant of the right square class, constructed
//     and then verified): products of small factors, "symmetric" constructions (g(x)g(-x), h(x^5), g(x^2)) and
//     h(x^2) with h a depressed quintic (-g dep).  The last one has a structurally degenerate Euclid chain (the
//     power sums p1=p2=p3=0), which is what the filter's general path and its zero-pivot test exist for, and which
//     polynomials from the real search have often.  Random dense polynomials almost never reach that path.
//     Optionally the application's own passing polynomials can be used as seeds (-real FILE).
//  3. Block sizes.  The original always calls the filter with 10,000 polynomials, so the lanes that do not fill a SIMD
//     block (the "tails") have the same size in every call.  Here the sizes vary (1, 7, 63, 64, 65, 127, 255, 1000,
//     4096, 10000).
//  4. Both builds, AVX2/FMA and baseline, are tested in the same run (the original tests whichever one the CPU picks).
//  5. The generators are checked, not trusted: the exact discriminant of a sample of the polynomials is computed
//     (the first 500 of every generator, then a fraction) and must be S-smooth times a square, and on a sample of the
//     transformation steps Disc(new)/Disc(old) must be an integer square.  A violation is reported as a generator bug.
//  6. The run can be repeated (the seed is printed, -seed), the failing polynomials are printed with the name of the
//     generator that made them, the exit status is non-zero on any failure (the original always returned 0), and
//     the numbers of the original are still printed ("Total polynomials tested", "Total failures").
//
// WHAT THAT BUYS.  Measured on Windows x86-64 with PARI 2.8, the unmodified pdfilter.f90 of this directory, and
// deliberately broken variants of it (one run each, seed 5; numbers are polynomials the filter wrongly REJECTED,
// AVX2 build / baseline build; NOT run on Linux):
//
//    fault put into the filter                         Eric's test as it was      this version, -g all, 1,000,000
//                                                      (his polys, block 10,000)  polynomials
//    gfortran -mfpmath=387 (x87, as a 32-bit build)    0 / 0 of 400,000           198 / 196   (same 192 polys, varied blocks: 82 / 70)
//    rounding mode "up" around the filter (-rm up)     0 of 400,000               37,746 (3.8 %)   [ts 271, mixed 1,565, sym 3,557, dep 113,743]
//    large primes up to 2^29 (silent: ~1 in 500,000)   0 / 0 of 400,000           1 / 2
//    small primes up to 2^14 (loud)                    32,777 / 3,073 (18 distinct polynomials)    94,877 / 8,540 (78,945 distinct)
//    rounding down / toward zero, -Ofast               0 (filter rejects nothing)  0   <- see the limit below
//
// Why: the 192 polynomials never reach the filter's degenerate (general) path, and their tails always have one size;
// the rounding-up fault only shows where a residue is zero mod q without being literally 0 (a cancellation
// lc*a - top*b = k*q with k != 0, which rounding up turns into -q, so the zero-pivot test misses it): in the dep
// family the Euclid chain is degenerate at all 30 primes and, after a shift, those zeros are such cancellations;
// and a rare fault needs a million different polynomials.  Among the 123,274 distinct polynomials that the
// application finds passing the discriminant test in upstream's 31 test work units, every one has a zero coefficient,
// 55 % are even polynomials and 25 % have bb1 = bb2 = bb3 = 0 (as dep does): dense random polynomials are not what
// the filter sees in practice.
//
// Reproduce on Linux (this directory):   bash build-linux.sh && ./test_pdfilter 100                    -> 0 failures
//                                        FFLAGS_EXTRA=-mfpmath=387 bash build-linux.sh && ./test_pdfilter 100   -> x87 failures
//                                        bash build-linux.sh && ./test_pdfilter 100 -rm up                     -> rounding-up failures
//                                        ./test_pdfilter 100 -g orig -block 10000   is Eric's configuration (-g orig = his generator)
//
// WHAT IT CANNOT SEE.  Because it only uses polynomials that must be kept, it only detects a filter that rejects
// TOO MUCH.  A filter that rejects too little (rounding down or toward zero, -Ofast, a self-test that switched the
// filter off) keeps everything and passes.  Detecting that needs polynomials that must be rejected; they were left out.
//
// Build: bash build-linux.sh.  The PARI calls used were checked against the declarations of PARI 2.17.3; the program was
// compiled and run against PARI 2.8 only.
// ---------------------------------------------------------------------------------------------------------

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <ctime>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <chrono>
#include <unordered_set>
#include <fenv.h>
#if defined(__SSE__)
#include <xmmintrin.h>
#endif

typedef long long i64;
typedef unsigned long long u64;
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

#include "pari.h"      // NB: on mingw64 this turns `long` into `long long`; below, use i64 / u64 and never the word `long`

namespace pdt {

struct Poly { i64 c[11]; };          // c[0] = 1, c[k] = coefficient of x^(10-k)

// ---------------------------------------------------------------- random numbers (one splitmix64 stream, seeded from -seed)
static u64 rs = 1;
static u64 rnd64() { rs += 0x9E3779B97F4A7C15ull; u64 z = rs; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }
static int rnd_int(int n) { return (int)(rnd64() % (u64)n); }                           // 0 .. n-1
static i64 rnd_range(i64 a, i64 b) { return a + (i64)(rnd64() % (u64)(b - a + 1)); }
static double rnd_unif() { return (rnd64() >> 11) * (1.0 / 9007199254740992.0); }

// ---------------------------------------------------------------- another rounding mode, switched on only around the filter call
// (PARI's own arithmetic, used here to generate and to check, is not exact under a non-default rounding mode.)
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

// ---------------------------------------------------------------- the square class: Disc = +-(S-smooth) * square, as upstream's test decides
static std::vector<int> Sset;
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
// Eric's invariant: a transformation keeps the discriminant in its square class, Disc(new)/Disc(old) = integer square.
static double ratioFrac = 0.02;
static i64 ratioChecked = 0, ratioBad = 0;
static void ratioGEN(GEN fold, GEN fnew)
{
  if (rnd_unif() >= ratioFrac) return;
  pari_sp av = avma;
  GEN Do = poldisc0(fold, -1), Dn = poldisc0(fnew, -1), r;
  ratioChecked++;
  if (gequal0(Do) || gequal0(Dn)) ratioBad++;
  else { GEN q = dvmdii(Dn, Do, &r); if (signe(r) != 0 || Z_issquare(q) == 0) ratioBad++; }
  avma = av;
}
static void ratioPoly(const Poly &a, const Poly &b) { pari_sp av = avma; ratioGEN(polyToGEN(a), polyToGEN(b)); avma = av; }

// ---------------------------------------------------------------- transformations that keep the square class
// (1) wider Tschirnhaus: y of random degree 1..dmax, coefficients in [-B,B]; result = characteristic polynomial of y(alpha)
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
    if (!gequal0(poldisc0(u, -1))) return u;                 // separable: a polynomial with a nonzero discriminant
    avma = av;
  }
  return NULL;
}
// (2) exact affine maps on the coefficient array (128 bit): x -> x + t, x -> c x (monic form c^10 f(x/c)), x -> -x
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

// ---------------------------------------------------------------- seed polynomials (all of them have the right square class, and are verified)
static GEN seed0;                    // x^10 - 12 x^5 - 4, the seed of the original test
static bool polyMul(const std::vector<i64> &a, const std::vector<i64> &b, std::vector<i64> &r)
{
  std::vector<i128> h(a.size() + b.size() - 1, 0);
  for (size_t i = 0; i < a.size(); i++) for (size_t j = 0; j < b.size(); j++) h[i + j] += (i128)a[i] * b[j];
  r.resize(h.size());
  const i128 lim = (i128)4000000000000000000ll;
  for (size_t i = 0; i < h.size(); i++) { if (h[i] > lim || h[i] < -lim) return false; r[i] = (i64)h[i]; }
  return true;
}
static i64 smoothVal(int maxExp)     // a random product of the primes of S
{
  i64 v = 1;
  for (size_t k = 0; k < Sset.size(); k++) { int e = rnd_int(maxExp + 1); for (int i = 0; i < e; i++) v *= Sset[k]; }
  return v;
}

// (a) products of small factors whose own discriminant has the right class (Disc of a product = prod Disc * prod Res^2)
struct Factor { std::vector<i64> a; int deg; };           // a[i] = coefficient of x^i
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
  for (int r = -9; r <= 9; r++) addFactor({(i64)-r, 1});                                      // linear
  for (int b = -6; b <= 6; b++) for (int k = 1; k <= 8; k++) for (int sg = -1; sg <= 1; sg += 2) for (int rep = 0; rep < 6; rep++) {
    i64 d = (i64)sg * smoothVal(2) * k * k;                                                   // x^2 + b x + c with b^2 - 4c = d
    if (((i64)b * b - d) % 4 != 0) continue;
    addFactor({((i64)b * b - d) / 4, (i64)b, 1});
  }
  for (int k = 2; k <= 10; k++) for (int rep = 0; rep < 25; rep++) {                          // x^k - a
    std::vector<i64> a(k + 1, 0); a[k] = 1; a[0] = -smoothVal(3) * (rnd_int(2) ? 1 : -1);
    addFactor(a);
  }
  for (int n = 3; n <= 30; n++) {                                                             // cyclotomic polynomials of degree <= 10
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
    for (int i = 0; i <= 10; i++) p.c[10 - i] = f[i];
    return p.c[0] == 1;
  }
  return false;
}

// (b) symmetric constructions
//   -g(x) g(-x) for a random monic quintic g: even, and Disc = Disc(g)^2 * Res(g(x), g(-x))^2 is a perfect square;
//   h(x^5), h quadratic with b^2 - 4c = +-(S-smooth) k^2:  Disc = +-5^10 c^4 Disc(h)^5, class of Disc(h);
//   g(x^2), g quintic with g(0) = +-(S-smooth) k^2:        Disc = +-2^10 g(0) Disc(g)^2, class of g(0).
static bool symmetric(Poly &p)
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
    for (int i = 0; i <= 10; i++) f[i] = -r[i];                // monic: the leading coefficient of g(x)g(-x) is -1
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
  for (int i = 0; i <= 10; i++) p.c[10 - i] = f[i];
  return p.c[0] == 1;
}

// (c) h(x^2), h a DEPRESSED quintic:  f = x^10 + b x^6 + c x^4 + d x^2 + e,  e = +-(S-smooth) k^2.
//     Disc(f) = +-2^10 e Disc(h)^2, class of e.  The power sums p1 = p2 = p3 = 0: the Euclid chain of (f, f') is degenerate
//     at the very first step, and stays so under x -> c x + t.
static bool depressedEven(Poly &p)
{
  int k = 1 + rnd_int(12), B = 1 + rnd_int(60);
  i64 e = (rnd_int(2) ? 1 : -1) * smoothVal(3) * k * k;
  std::vector<i64> f(11, 0);
  f[10] = 1; f[6] = rnd_range(-B, B); f[4] = rnd_range(-B, B); f[2] = rnd_range(-B, B); f[0] = e;
  for (int i = 0; i <= 10; i++) p.c[10 - i] = f[i];
  return true;
}

// (d) the application's own passing polynomials (-real FILE: bb1 .. bb10 per line, the leading 1 implied)
static std::vector<Poly> realSeeds;

// ---------------------------------------------------------------- generators
enum Gen { G_ORIG, G_TS, G_MIXED, G_SYM, G_DEP, G_REAL, G_ALL };
static const char *genName[] = { "orig", "ts", "mixed", "sym", "dep", "real", "all" };
static int tsB = 3, tsD = 9;

// a base polynomial: a seed, followed (with probability tsPercent) by a wider Tschirnhaus transformation
static bool finishBase(const Poly &q, Poly &p, int tsPercent)
{
  GEN f = polyToGEN(q);
  if (gequal0(poldisc0(f, -1))) return false;                // repeated factor: not separable, discriminant 0
  if (rnd_int(100) < tsPercent) {
    GEN u = tsOwn(f, tsD, 1 + rnd_int(tsB));
    if (u == NULL || !genToPoly(u, p)) return false;
    ratioGEN(f, u);
    return true;
  }
  p = q;
  return true;
}
static bool makeBase(Gen g, Poly &p)
{
  pari_sp av = avma;
  bool ok = false;
  Poly q;
  switch (g) {
    case G_ORIG: { GEN u = tschirnhaus(seed0); ok = genToPoly(u, p); if (ok) ratioGEN(seed0, u); break; }   // Eric's
    case G_TS: {
      GEN u = tsOwn(seed0, tsD, 1 + rnd_int(tsB));
      ok = u && genToPoly(u, p);
      if (ok) ratioGEN(seed0, u);
      break; }
    case G_MIXED: ok = mixedProduct(q) && finishBase(q, p, 67); break;
    case G_SYM:   ok = symmetric(q) && finishBase(q, p, 67); break;
    case G_DEP:   ok = depressedEven(q) && finishBase(q, p, 30); break;
    case G_REAL:  if (!realSeeds.empty()) ok = finishBase(realSeeds[rnd_int((int)realSeeds.size())], p, 67); break;
    default: break;
  }
  avma = av;
  return ok;
}

// ---------------------------------------------------------------- statistics
static u64 hashPoly(const Poly &p) { u64 h = 1469598103934665603ull; for (int i = 1; i <= 10; i++) { h ^= (u64)p.c[i]; h *= 1099511628211ull; h ^= h >> 29; } return h; }
static int bitLen(const Poly &p)
{
  int m = 0;
  for (int i = 1; i <= 10; i++) { u64 v = p.c[i] < 0 ? (u64)(-(p.c[i] + 1)) + 1 : (u64)p.c[i]; int b = 0; while (v) { b++; v >>= 1; } if (b > m) m = b; }
  return m;
}

// ---------------------------------------------------------------- the filter's primes (only for -pat): same selection as pdf_init
static bool primeU(u64 n) { if (n < 2) return false; for (u64 d = 2; d * d <= n; d++) if (n % d == 0) return false; return true; }
static u64 powmodU(u64 b, u64 e, u64 m) { u64 r = 1; b %= m; while (e) { if (e & 1) r = (u64)((i128)r * b % m); b = (u64)((i128)b * b % m); e >>= 1; } return r; }
static bool usablePrime(u64 q)
{
  if (q % 4 != 1 || !primeU(q)) return false;
  for (size_t k = 0; k < Sset.size(); k++) { if (q % (u64)Sset[k] == 0 || powmodU((u64)Sset[k], (q - 1) / 2, q) != 1) return false; }
  return true;
}
static std::vector<u64> filterPrimes()
{
  std::vector<u64> pr;
  int ns = 0, nq = 0;
  for (u64 q = 1023; ns < 6 && q > 256; q -= 2) if (usablePrime(q)) { pr.push_back(q); ns++; }
  for (u64 q = (1ull << 25) - 1; nq < 24 && q > (1ull << 24); q -= 2) if (usablePrime(q)) { pr.push_back(q); nq++; }
  return pr;
}

// ---------------------------------------------------------------- the filter builds under test
struct Build { const char *name; void (*init)(int, const int *); void (*filter)(const i64 *, int, signed char *); };

static void usage(const char *prog)
{
  printf("Usage: %s [N] [options]\n", prog);
  printf("  N             number of blocks of 10,000 polynomials (default 100); -n TOTAL gives the number of polynomials\n");
  printf("  -g GEN        orig | ts | mixed | sym | dep | real | all   (default all; real needs -real)\n");
  printf("  -S 2,5        the primes divided out by the upstream test (default 2,5; at most 3)\n");
  printf("  -seed X       random seed (default: the time; it is printed)\n");
  printf("  -block B      always call the filter with B polynomials (default: sizes vary)\n");
  printf("  -validate F   fraction of polynomials whose discriminant is checked exactly (default 0.05; the first 500 of each generator always)\n");
  printf("  -reuse K      exact affine variants made from every expensive base polynomial (default 8)\n");
  printf("  -pscr P       probability that a variant is scrambled by an affine map (default 0.85)\n");
  printf("  -ratio F      fraction of transformation steps whose Disc(new)/Disc(old) is checked to be a square (default 0.02)\n");
  printf("  -tsB B -tsD D coefficient bound (default 3) and maximal degree (default 9) of the Tschirnhaus polynomial\n");
  printf("  -real FILE    extra seeds, 10 integers per line (bb1..bb10, leading 1 implied)\n");
  printf("  -rm MODE      up | down | zero: run the filter under another rounding mode (test aid)\n");
  printf("  -pat K        print the factorisation patterns of K polynomials modulo the filter primes\n");
}

static int run(int argc, char **argv)
{
  Gen gen = G_ALL; i64 N = 100 * 10000; int fixedBlock = 0, patN = 0; u64 seed = (u64)time(NULL);
  double validateFrac = 0.05, pscr = 0.85; int reuse = 8;
  std::string realFile, Sarg = "2,5";
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    auto nxt = [&]() -> std::string { return i + 1 < argc ? std::string(argv[++i]) : std::string(); };
    if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
    else if (a == "-g") { std::string v = nxt(); bool found = false; for (int k = 0; k <= G_ALL; k++) if (v == genName[k]) { gen = (Gen)k; found = true; } if (!found) { usage(argv[0]); return 2; } }
    else if (a == "-n") N = atoll(nxt().c_str());
    else if (a == "-S") Sarg = nxt();
    else if (a == "-seed") seed = (u64)atoll(nxt().c_str());
    else if (a == "-block") fixedBlock = atoi(nxt().c_str());
    else if (a == "-validate") validateFrac = atof(nxt().c_str());
    else if (a == "-reuse") reuse = atoi(nxt().c_str());
    else if (a == "-pscr") pscr = atof(nxt().c_str());
    else if (a == "-ratio") ratioFrac = atof(nxt().c_str());
    else if (a == "-tsB") tsB = atoi(nxt().c_str());
    else if (a == "-tsD") tsD = atoi(nxt().c_str());
    else if (a == "-real") realFile = nxt();
    else if (a == "-pat") patN = atoi(nxt().c_str());
    else if (a == "-rm") { std::string v = nxt(); hostile = v == "up" ? 1 : v == "down" ? 2 : v == "zero" ? 3 : 0; if (!hostile) { usage(argv[0]); return 2; } }
    else if (a[0] != '-' && atoll(a.c_str()) > 0) N = atoll(a.c_str()) * 10000;          // Eric's usage: ./test_pdfilter N
    else { usage(argv[0]); return 2; }
  }
  if (reuse < 1) reuse = 1;
  rs = seed;

  pari_init(100000000, 10000);                                 // 100 MB, as in the original
  setrand(utoi((u64)seed));
  LIM = shifti(gen_1, 63); NLIM = negi(LIM);
  { size_t pos = 0; while (pos < Sarg.size()) { size_t e = Sarg.find(',', pos); if (e == std::string::npos) e = Sarg.size(); Sset.push_back(atoi(Sarg.substr(pos, e - pos).c_str())); pos = e + 1; } }
  if (Sset.empty() || Sset.size() > 3) { printf("-S: one to three primes\n"); return 2; }
  seed0 = gclone(mkpoln(11, gen_1, gen_0, gen_0, gen_0, gen_0, stoi(-12), gen_0, gen_0, gen_0, gen_0, stoi(-4)));
  bool seed0ok = passesUpstream(poldisc0(seed0, -1));
  buildPool();
  if (!realFile.empty()) {
    FILE *f = fopen(realFile.c_str(), "r"); char line[2048]; std::set<u64> seen;
    if (!f) { printf("cannot read %s\n", realFile.c_str()); return 2; }
    while (fgets(line, sizeof line, f)) {
      if (line[0] == '#') continue;
      Poly q; q.c[0] = 1; int n = 0; char *s = line;
      for (; n < 10; n++) { char *e; i64 v = strtoll(s, &e, 10); if (e == s) break; q.c[1 + n] = v; s = e; }
      if (n == 10 && seen.insert(hashPoly(q)).second) realSeeds.push_back(q);
    }
    fclose(f);
  }

  // which generators, with what weight
  std::vector<Gen> mix;
  if (gen == G_ALL) {
    if (seed0ok) { mix.push_back(G_ORIG); mix.push_back(G_TS); mix.push_back(G_TS); }
    mix.push_back(G_MIXED); mix.push_back(G_MIXED); mix.push_back(G_SYM); mix.push_back(G_DEP); mix.push_back(G_DEP);
    if (!realSeeds.empty()) { mix.push_back(G_REAL); mix.push_back(G_REAL); }
  } else {
    if ((gen == G_ORIG || gen == G_TS) && !seed0ok) { printf("the seed x^10-12x^5-4 needs 2 and 5 in S; use -g mixed, sym or dep\n"); return 2; }
    if (gen == G_REAL && realSeeds.empty()) { printf("-g real needs -real FILE\n"); return 2; }
    mix.push_back(gen);
  }

  // builds
  std::vector<Build> builds;
  int pS[3] = { 0, 0, 0 }; for (size_t k = 0; k < Sset.size(); k++) pS[k] = Sset[k];
#ifdef NF_PDF_DISPATCH
  __builtin_cpu_init();
  const char *fb = getenv("NF_FORCE_BASE");
  if (!(fb && *fb == '1') && __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")) { Build b = { "AVX2/FMA", pdf_init_v3, pdf_filter_v3 }; builds.push_back(b); }
  { Build b = { "baseline", pdf_init_base, pdf_filter_base }; builds.push_back(b); }
#else
  { Build b = { "single build", pdf_init, pdf_filter }; builds.push_back(b); }
#endif
  for (size_t b = 0; b < builds.size(); b++) builds[b].init((int)Sset.size(), pS);

  printf("pdfilter test: generator=%s, S=%s, %lld polynomials, seed=%llu, builds:", genName[gen], Sarg.c_str(), (i64)N, (u64)seed);
  for (size_t b = 0; b < builds.size(); b++) printf(" %s", builds[b].name);
  printf("%s\n", hostile ? "  [filter runs under a non-default rounding mode]" : "");
  if (!realSeeds.empty()) printf("real seeds: %zu\n", realSeeds.size());
  if (!seed0ok) printf("note: the seed x^10-12x^5-4 does not pass for this S; generators orig/ts are off\n");

  const int blockSizes[] = { 1, 7, 63, 64, 65, 127, 255, 1000, 4096, 10000 };
  std::unordered_set<u64> distinct; const size_t distinctCap = 5000000; bool capped = false;
  std::map<int, i64> bitBucket;
  i64 tested = 0, invalid = 0, genFail = 0;
  std::vector<i64> genTested(G_ALL, 0), genChecked(G_ALL, 0);
  std::vector<std::vector<i64> > genFail_b(G_ALL, std::vector<i64>(builds.size(), 0));
  std::vector<i64> failTotal(builds.size(), 0);
  std::set<u64> failDistinct; int shown = 0;
  std::vector<u64> fprimes; std::map<std::string, i64> patGlobal; std::set<std::pair<int, std::string> > patPairs; i64 patSampled = 0;
  if (patN > 0) fprimes = filterPrimes();
  std::vector<Poly> queue; std::vector<Gen> queueGen;
  std::vector<i64> buf; std::vector<signed char> flag; std::vector<Gen> blockGen; std::vector<Poly> blockPoly;
  auto T0 = std::chrono::steady_clock::now();
  pari_sp ltop = avma;

  while (tested < N) {
    int bs = fixedBlock ? fixedBlock : blockSizes[rnd_int(10)];
    if ((i64)bs > N - tested) bs = (int)(N - tested);
    buf.assign((size_t)bs * 11, 0); blockGen.assign(bs, G_ORIG); blockPoly.assign(bs, Poly());
    int got = 0;
    while (got < bs) {
      avma = ltop;
      if (queue.empty()) {
        Gen gg = mix[rnd_int((int)mix.size())];
        Poly base;
        if (makeBase(gg, base)) {
          int K = (gg == G_ORIG) ? 1 : reuse;
          for (int v = 0; v < K; v++) {
            Poly q = base;
            if (gg != G_ORIG && (v > 0 || rnd_unif() < pscr) && !scramble(q)) continue;
            if (gg != G_ORIG && rnd_unif() < ratioFrac) ratioPoly(base, q);
            queue.push_back(q); queueGen.push_back(gg);
          }
        } else if (++genFail > 100 * (N + 1000)) { printf("a generator keeps failing\n"); return 3; }
        continue;
      }
      Poly p = queue.back(); Gen gg = queueGen.back(); queue.pop_back(); queueGen.pop_back();
      bool check = genChecked[gg] < 500 || validateFrac >= 1.0 || rnd_unif() < validateFrac || patSampled < patN;
      GEN f = polyToGEN(p);
      if (check) {
        genChecked[gg]++;
        if (!passesUpstream(poldisc0(f, -1))) {          // the polynomial left the square class: a bug in a generator
          invalid++;
          if (invalid <= 5) { printf("GENERATOR BUG (%s): discriminant is not S-smooth * square:", genName[gg]); for (int k = 1; k <= 10; k++) printf(" %lld", (i64)p.c[k]); printf("\n"); }
          continue;
        }
      }
      for (int k = 0; k <= 10; k++) buf[(size_t)got * 11 + k] = p.c[k];
      blockGen[got] = gg; blockPoly[got] = p;
      if (distinct.size() < distinctCap) distinct.insert(hashPoly(p)); else capped = true;
      bitBucket[(bitLen(p) + 7) / 8]++;
      if (patSampled < patN) {                            // factorisation pattern modulo every filter prime
        patSampled++;
        for (size_t k = 0; k < fprimes.size(); k++) {
          pari_sp av = avma;
          GEN fa = factormod0(f, utoipos(fprimes[k]), 1);
          std::vector<int> degs;
          for (int r = 1; r < lg(gel(fa, 1)); r++) { int d = (int)itos(gcoeff(fa, r, 1)), e = (int)itos(gcoeff(fa, r, 2)); for (int z = 0; z < e; z++) degs.push_back(d); }
          std::sort(degs.begin(), degs.end());
          std::string s; for (size_t z = 0; z < degs.size(); z++) { s += std::to_string(degs[z]); s += ' '; }
          patPairs.insert(std::make_pair((int)k, s)); patGlobal[s]++;
          avma = av;
        }
      }
      genTested[gg]++;
      got++;
    }
    for (size_t b = 0; b < builds.size(); b++) {
      flag.assign(bs, 99);
      roundOn();
      builds[b].filter(buf.data(), bs, flag.data());
      roundOff();
      for (int i = 0; i < bs; i++) {
        if (flag[i] != 0) continue;                        // kept: as it must be
        failTotal[b]++; genFail_b[blockGen[i]][b]++;
        failDistinct.insert(hashPoly(blockPoly[i]));
        if (shown < 10) { shown++; printf("  REJECTED by %s (%s, block of %d): x^10", builds[b].name, genName[blockGen[i]], bs); for (int k = 1; k <= 10; k++) printf(" %+lld", (i64)blockPoly[i].c[k]); printf("\n"); }
      }
    }
    tested += bs;
  }

  double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
  i64 fails = 0; for (size_t b = 0; b < builds.size(); b++) fails += failTotal[b];
  printf("---- results (%.1f s) ----\n", wall);
  printf("distinct polynomials      : %s%zu of %lld tested\n", capped ? ">= " : "", distinct.size(), (i64)tested);
  printf("coefficient size (bits of the largest coefficient), count per range:");
  { static const char *nm[] = { "<=8", "9-16", "17-24", "25-32", "33-40", "41-48", "49-56", "57-63" };
    for (std::map<int, i64>::iterator it = bitBucket.begin(); it != bitBucket.end(); ++it) printf(" %s:%lld", nm[std::min(7, std::max(0, it->first - (it->first > 0 ? 1 : 0)))], (i64)it->second);
    printf("\n"); }
  printf("polynomials per generator :"); for (int g = 0; g < G_ALL; g++) if (genTested[g]) printf(" %s=%lld", genName[g], (i64)genTested[g]); printf("\n");
  { i64 chk = 0; for (int g = 0; g < G_ALL; g++) chk += genChecked[g];
    printf("discriminants verified exactly: %lld polynomials; transformation steps checked (Disc ratio is a square): %lld\n", (i64)chk, (i64)ratioChecked); }
  if (patN > 0) {
    printf("factorisation patterns mod the %zu filter primes, %lld polynomials: %zu distinct (prime,pattern) pairs, %zu distinct patterns (of 22 possible)\n",
           fprimes.size(), (i64)patSampled, patPairs.size(), patGlobal.size());
  }
  for (size_t b = 0; b < builds.size(); b++) {
    printf("[%s] valid polynomials REJECTED by the filter: %lld", builds[b].name, (i64)failTotal[b]);
    if (failTotal[b]) { printf("   (by generator:"); for (int g = 0; g < G_ALL; g++) if (genFail_b[g][b]) printf(" %s=%lld", genName[g], (i64)genFail_b[g][b]); printf(")"); }
    printf("\n");
  }
  if (!failDistinct.empty()) printf("distinct failing polynomials: %zu\n", failDistinct.size());
  printf("Total polynomials tested = %lld\n", (i64)tested);
  printf("Total failures = %lld\n", (i64)fails);
  if (invalid || ratioBad) { printf("GENERATOR BUG: %lld polynomials left the square class of the discriminant, %lld Disc-ratio violations\n", (i64)invalid, (i64)ratioBad); return 2; }
  return fails ? 1 : 0;
}

}  // namespace pdt

int main(int argc, char **argv) { return pdt::run(argc, argv); }
