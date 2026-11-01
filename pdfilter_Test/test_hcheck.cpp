// test_hcheck.cpp -- could hypothesis (H) be checked at run time instead of assumed?  (H_CHECK.md has the answer; this program
// produces its numbers.)  (H): no int64 expression of upstream's first two sub-resultant iterations overflows (pdfilter_chains.tex,
// Proposition "What the exact test computes").  Built by hcheck_checks.sh with upstream's polDiscTest_cpu.cpp twice (as shipped,
// and with -fwrapv, linked as polDiscTest_wrapv) and the two builds of the filter; no PARI.
//   build/hcheck/test_hcheck [COUNT] [SEED] [-hits FILE]
// 1. hcheck.h's nf_h_ok and nf_h_ok_fast against the 128-bit replay of test_poldisc.cpp (upstream_overflows), COUNT polynomials per
//    family (default 400000): the same verdict on every polynomial?
// 2. how far the triangle-inequality bound of the int64 values is above the values (could a cheap screen replace the exact check?);
// 3. polynomials outside (H) that upstream passes and the filter rejects: f = g + 2^62 (constant term), g a small polynomial that
//    passes, with B9 = 10/cont(f') even;
// 4. the final flags of the stock call (exact test alone), of today's application (filter, then the exact test) and of the
//    application with a per-polynomial bypass (filter, then flag 1 again where (H) fails, then the exact test);
// 5. -fwrapv: does compiling upstream's test with defined signed overflow change any flag?
// 6. the cost of the check next to the filter, per polynomial.
// -hits FILE adds the results (one polynomial per line, leading coefficient first, 10 or 11 numbers) with x -> -x and shifts |t| <= 40.
// Exit status 1 if a check disagrees with the replay, the bypass gives a flag the stock call does not, or -fwrapv changes a flag.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include <algorithm>
#include "hcheck.h"

extern "C" int polDiscTest_cpu(long long *polBuf, int numPolys, char *polGoodFlag, char *polyMask, int numP, int *pSet);
extern "C" int polDiscTest_wrapv(long long *polBuf, int numPolys, char *polGoodFlag, char *polyMask, int numP, int *pSet);
extern "C" void pdf_init_v3(int, const int*);
extern "C" void pdf_filter_v3(const long long*, int, signed char*);
extern "C" void pdf_init_base(int, const int*);
extern "C" void pdf_filter_base(const long long*, int, signed char*);

typedef long long i64;
typedef __int128 L;
static uint64_t rs;
static uint64_t rnd() { rs += 0x9E3779B97F4A7C15ull; uint64_t z = rs; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }
static i64 urand(i64 b) { return (i64)(rnd() % (uint64_t)(2 * b + 1)) - b; }
static double now() { timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
static int pSet[3] = {2, 5, 0};
static const int nS = 2;

// test_poldisc.cpp's replay of upstream's int64 steps in 128 bits (verbatim), extended to return the largest |value| (vmax) and the
// same expressions with every term taken in absolute value (tri: the triangle-inequality bound, in long double).
static bool upstream_overflows(const i64 *c, double *vmax = 0, double *tri = 0)
{
  bool ov = false;
  L mx = 0; long double mt = 0;
  auto chk = [&](L v) { if (v > (L)INT64_MAX || v < (L)INT64_MIN) ov = true; L x = v < 0 ? -v : v; if (x > mx) mx = x; return v; };
  auto tb = [&](long double w) { if (w > mt) mt = w; return w; };
  L A[11], B[10], R[10]; long double a[11], b[10], r[10];
  for (int col = 0; col < 11; col++) { A[10 - col] = c[col]; a[10 - col] = fabsl((long double)c[col]); }
  for (int k = 1; k <= 10; k++) { B[k - 1] = chk((L)k * A[k]); b[k - 1] = tb(k * a[k]); }
  int cont2 = 2; for (int k = 0; k < 10; k += 2) if ((B[k] & 1) != 0) { cont2 = 1; break; }
  int cont5 = 5; for (int k = 0; k < 9; k++) if (B[k] % 5 != 0) { cont5 = 1; break; }
  for (int k = 0; k < 10; k++) { B[k] = B[k] / (cont2 * cont5); b[k] /= (cont2 * cont5); }
  R[0] = chk(B[9] * A[0]); r[0] = tb(b[9] * a[0]);
  for (int k = 1; k <= 9; k++) { R[k] = chk(chk(B[9] * A[k]) - B[k - 1]); r[k] = tb(b[9] * a[k] + b[k - 1]); }
  for (int k = 0; k <= 8; k++) { R[k] = chk(chk(B[9] * R[k]) - chk(R[9] * B[k])); r[k] = tb(b[9] * r[k] + r[9] * b[k]); }
  int degR = 0; for (int k = 8; k > 0; k--) if (R[k] != 0) { degR = k; break; }
  if (!ov && degR == 8) {
    for (int k = 0; k < 10; k++) { A[k] = B[k]; a[k] = b[k]; }
    for (int k = 0; k <= degR; k++) { B[k] = R[k]; b[k] = r[k]; }
    R[0] = chk(B[8] * A[0]); r[0] = tb(b[8] * a[0]);
    for (int k = 1; k <= 8; k++) { R[k] = chk(chk(B[8] * A[k]) - chk(A[9] * B[k - 1])); r[k] = tb(b[8] * a[k] + a[9] * b[k - 1]); }
    for (int k = 0; k <= 7; k++) { R[k] = chk(chk(B[8] * R[k]) - chk(R[8] * B[k])); r[k] = tb(b[8] * r[k] + r[8] * b[k]); }
    chk(B[9] * B[9]);
  } else if (!ov && degR == 7) {
    for (int k = 0; k < 10; k++) { A[k] = B[k]; a[k] = b[k]; }
    for (int k = 0; k <= degR; k++) { B[k] = R[k]; b[k] = r[k]; }
    R[0] = chk(B[7] * A[0]); R[1] = chk(B[7] * A[1]); r[0] = tb(b[7] * a[0]); r[1] = tb(b[7] * a[1]);
    for (int k = 2; k <= 8; k++) { R[k] = chk(chk(B[7] * A[k]) - chk(A[9] * B[k - 2])); r[k] = tb(b[7] * a[k] + a[9] * b[k - 2]); }
    R[0] = chk(B[7] * R[0]); r[0] = tb(b[7] * r[0]);
    for (int k = 1; k <= 7; k++) { R[k] = chk(chk(B[7] * R[k]) - chk(R[8] * B[k - 1])); r[k] = tb(b[7] * r[k] + r[8] * b[k - 1]); }
    for (int k = 0; k <= 6; k++) { R[k] = chk(chk(B[7] * R[k]) - chk(R[7] * B[k])); r[k] = tb(b[7] * r[k] + r[7] * b[k]); }
    chk(chk(B[9] * B[9]) * B[9]);
    chk(B[7] * B[7]); tb(b[7] * b[7]);
  }
  if (vmax) *vmax = mx > 0 ? log2((double)mx) : 0;
  if (tri) *tri = mt > 0 ? (double)log2l(mt) : 0;
  return ov;
}

static bool shift(const i64 *c, i64 t, i64 *out)          // f(x + t), leading first; false if a coefficient leaves int64
{
  L g[11];
  for (int k = 0; k <= 10; k++) g[k] = c[10 - k];
  for (int i = 0; i < 10; i++) for (int k = 9; k >= i; k--) g[k] = g[k] + (L)t * g[k + 1];
  for (int k = 0; k <= 10; k++) { if (g[k] > (L)INT64_MAX || g[k] < (L)INT64_MIN) return false; out[10 - k] = (i64)g[k]; }
  return true;
}

static std::vector<std::vector<i64>> hits;
static void readPolys(const char *file)
{
  FILE *f = fopen(file, "r"); char line[4096];
  if (!f) { fprintf(stderr, "cannot read %s\n", file); exit(2); }
  while (fgets(line, sizeof line, f)) {
    if (line[0] == '#') continue;
    i64 v[13]; int n = 0; char *s = line;
    while (n < 13) { char *e; i64 x = strtoll(s, &e, 10); if (e == s) break; v[n++] = x; s = e; }
    std::vector<i64> p(11);
    if (n == 10) { p[0] = 1; for (int k = 0; k < 10; k++) p[1 + k] = v[k]; }
    else if (n >= 11 && v[0] == 1) { for (int k = 0; k <= 10; k++) p[k] = v[k]; }
    else continue;
    hits.push_back(p);
  }
  fclose(f);
}

// Families.  randB: coefficients uniform in [-B, B].  passbig: test_poldisc.cpp's polynomials built to pass with large coefficients.
// wR: the coefficient of x^(10-j) uniform in +-C(10,j) R^j / 4, the size of a polynomial whose roots have modulus about R (the real
// candidates have coefficients up to 2.1e9, between w8 and w9).  hits: the results of -hits, x or -x, shifted by |t| <= 40.
static void make(const std::string &fam, i64 *c)
{
  c[0] = 1;
  if (fam.compare(0, 4, "rand") == 0) {
    i64 b = atoll(fam.c_str() + 4);
    for (int k = 1; k <= 10; k++) c[k] = urand(b);
  } else if (fam[0] == 'w') {
    static const double binom[11] = {1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1};
    double r = atof(fam.c_str() + 1), p = 1;
    for (int j = 1; j <= 10; j++) { p *= r; double b = binom[j] * p * 0.25; c[j] = urand((i64)(b < 1 ? 1 : b)); }
  } else if (fam == "passbig") {
    for (;;) {
      i64 G[6]; G[5] = 1;
      for (int k = 1; k <= 4; k++) G[k] = urand(1000000);
      i64 m = 1 + (i64)(rnd() % 100000);
      static const i64 u[8] = {1, -1, 2, -2, 5, -5, 10, -10};
      G[0] = u[rnd() % 8] * m * m;
      i64 f[11]; for (int k = 0; k <= 10; k++) f[k] = 0;
      for (int k = 0; k <= 5; k++) f[10 - 2 * k] = G[k];
      if (shift(f, urand(2), c)) return;
    }
  } else if (fam == "hits") {
    for (;;) {
      const std::vector<i64> &p = hits[rnd() % hits.size()];
      i64 g[11]; int s = rnd() & 1;
      for (int k = 0; k <= 10; k++) g[k] = (s && (k % 2)) ? -p[k] : p[k];
      if (shift(g, urand(40), c)) return;
    }
  } else { fprintf(stderr, "family %s?\n", fam.c_str()); exit(2); }
}

// upstream's B9 = 10/cont(f') is even exactly when an even-index coefficient (k+1)*a_(k+1) of f' is odd (cont2 = 1).  Then upstream's
// int64 steps see a0 only through B9^2*a0 with 4 | B9^2, so a0 + 2^62 gives, after the wrap-around, the same numbers as a0.
static bool b9_even(const i64 *g) { for (int k = 0; k < 10; k += 2) if ((((i64)(k + 1) * g[10 - (k + 1)]) & 1) == 1) return true; return false; }
static void plus262(const i64 *g, i64 *f) { for (int k = 0; k <= 10; k++) f[k] = g[k]; f[10] = g[10] + (g[10] < 0 ? 1 : -1) * ((i64)1 << 62); }

int main(int argc, char **argv)
{
  i64 count = argc > 1 && argv[1][0] != '-' ? atoll(argv[1]) : 400000;
  rs = (argc > 2 && argv[2][0] != '-' ? strtoull(argv[2], 0, 10) : 1) * 0x9E3779B97F4A7C15ull + 7;
  for (int i = 1; i + 1 < argc; i++) if (!strcmp(argv[i], "-hits")) readPolys(argv[i + 1]);
  pdf_init_v3(nS, pSet); pdf_init_base(nS, pSet);
  int bad = 0;

  // ---- 1 ----
  printf("1. the check against the 128-bit replay, %lld polynomials per family\n", count);
  printf("   %-9s %10s %9s   largest |int64 value| where (H) holds, log2: median / max\n", "family", "(H) fails", "disagree");
  std::vector<std::string> fams = {"rand3", "rand100", "rand3000", "rand1000000", "passbig", "w4", "w6", "w8", "w9", "w10", "w11", "w12"};
  if (!hits.empty()) fams.push_back("hits");
  for (const std::string &fam : fams) {
    i64 fails = 0, dis = 0; i64 c[11]; std::vector<double> mb;
    for (i64 i = 0; i < count; i++) {
      make(fam, c);
      double vm; bool ov = upstream_overflows(c, &vm);
      int ok = nf_h_ok(c);
      fails += ov; dis += (ok == (int)ov) + (nf_h_ok_fast(c) != ok);
      if (!ov) mb.push_back(vm);
    }
    std::sort(mb.begin(), mb.end());
    printf("   %-9s %10lld %9lld", fam.c_str(), fails, dis);
    if (!mb.empty()) printf("   %.2f / %.2f", mb[mb.size() / 2], mb.back());
    printf("\n");
    bad |= dis != 0;
  }

  // ---- 2 ----
  printf("\n2. the triangle-inequality bound against the true largest value, polynomials where (H) holds\n");
  printf("   %-6s %24s %26s\n", "family", "bound - true (bits): median / 99% / max", "share with bound >= 2^62");
  for (const char *fam : {"w6", "w8", "w9", "w10"}) {
    std::vector<double> gap; i64 n = 0, over = 0; i64 c[11];
    for (i64 i = 0; i < count; i++) {
      make(fam, c);
      double vm, tr; if (upstream_overflows(c, &vm, &tr)) continue;
      n++; gap.push_back(tr - vm); over += tr >= 62;
    }
    std::sort(gap.begin(), gap.end());
    printf("   %-6s %14.2f / %5.2f / %5.2f %22.2f%%\n", fam, gap[n / 2], gap[(size_t)(n * 0.99)], gap.back(), 100.0 * over / n);
  }

  // ---- 3 ----
  printf("\n3. f = g +- 2^62 (constant term), g monic with |a_i| <= 3, B9 even, passed by upstream\n");
  std::vector<i64> gs, fs;                         // the pairs, for part 4
  {
    const int n = 4096; std::vector<i64> buf(n * 11); std::vector<char> up(n); std::vector<signed char> fl(n);
    int found = 0, upf = 0, upw = 0, rej = 0, hok = 0; i64 tried = 0;
    while (found < 100 && tried < 100000000) {
      for (int i = 0; i < n; i++) { i64 *c = &buf[i * 11]; c[0] = 1; for (int k = 1; k <= 10; k++) c[k] = urand(3); }
      pdf_filter_v3(buf.data(), n, fl.data());
      for (int i = 0; i < n; i++) up[i] = fl[i];
      polDiscTest_cpu(buf.data(), n, up.data(), up.data(), nS, pSet);
      tried += n;
      for (int i = 0; i < n && found < 100; i++) {
        const i64 *g = &buf[i * 11];
        if (!up[i] || !b9_even(g)) continue;
        i64 f[11]; plus262(g, f);
        char uf, uw; signed char ff;
        polDiscTest_cpu(f, 1, &uf, NULL, nS, pSet);
        polDiscTest_wrapv(f, 1, &uw, NULL, nS, pSet);
        pdf_filter_v3(f, 1, &ff);
        if (found < 2) {
          printf("   g ="); for (int k = 0; k <= 10; k++) printf(" %lld", g[k]);
          printf("\n   f ="); for (int k = 0; k <= 10; k++) printf(" %lld", f[k]);
          printf("\n      upstream(f) = %d (-fwrapv: %d), filter(f) = %d, (H) holds for f: %d\n", uf, uw, ff, nf_h_ok(f));
        }
        found++; upf += uf; upw += uw; rej += ff == 0; hok += nf_h_ok(f);
        gs.insert(gs.end(), g, g + 11); fs.insert(fs.end(), f, f + 11);
      }
    }
    printf("   %d pairs (%lld polynomials searched): upstream passes %d of the f (-fwrapv: %d), the filter rejects %d, (H) holds for %d\n",
           found, tried, upf, upw, rej, hok);
  }

  // ---- 4 ----
  printf("\n4. final flags: stock call, today's application, application with the per-polynomial bypass\n");
  printf("   %-12s %9s %10s %11s %14s %15s\n", "family", "polys", "(H) fails", "stock pass", "today != stock", "bypass != stock");
  std::vector<std::string> fams4 = {"rand1000000", "passbig", "w12", "constructed"};
  if (!hits.empty()) fams4.insert(fams4.begin() + 3, "hits");
  for (const std::string &fam : fams4) {
    const int n = 4096; std::vector<i64> buf(n * 11); std::vector<char> st(n), a(n), b(n); std::vector<signed char> fl(n);
    i64 tot = 0, fails = 0, sp = 0, d1 = 0, d2 = 0;
    for (i64 done = 0; done < (fam == "constructed" ? 1 : count / 4); ) {
      int m = n;
      if (fam == "constructed") {                   // the f of part 3 and their g
        m = 0;
        for (size_t i = 0; i < gs.size() / 11 && m + 2 <= n; i++) {
          memcpy(&buf[m++ * 11], &gs[i * 11], 11 * sizeof(i64));
          memcpy(&buf[m++ * 11], &fs[i * 11], 11 * sizeof(i64));
        }
      } else for (int i = 0; i < n; i++) make(fam, &buf[i * 11]);
      done += m;
      polDiscTest_cpu(buf.data(), m, st.data(), NULL, nS, pSet);
      pdf_filter_v3(buf.data(), m, fl.data());
      for (int i = 0; i < m; i++) a[i] = fl[i];
      polDiscTest_cpu(buf.data(), m, a.data(), a.data(), nS, pSet);
      for (int i = 0; i < m; i++) b[i] = (fl[i] || !nf_h_ok(&buf[i * 11])) ? 1 : 0;
      polDiscTest_cpu(buf.data(), m, b.data(), b.data(), nS, pSet);
      for (int i = 0; i < m; i++) { tot++; fails += !nf_h_ok(&buf[i * 11]); sp += st[i]; d1 += a[i] != st[i]; d2 += b[i] != st[i]; }
    }
    printf("   %-12s %9lld %10lld %11lld %14lld %15lld\n", fam.c_str(), tot, fails, sp, d1, d2);
    bad |= d2 != 0;
  }

  // ---- 5 ----
  printf("\n5. upstream's test compiled with -fwrapv against the test as shipped\n");
  for (const char *fam : {"rand3", "rand100", "rand3000", "rand1000000", "passbig", "w12"}) {
    const int n = 4096; std::vector<i64> buf(n * 11); std::vector<char> x(n), y(n);
    i64 tot = 0, fails = 0, d = 0;
    for (i64 done = 0; done < count / 8; done += n) {
      for (int i = 0; i < n; i++) make(fam, &buf[i * 11]);
      polDiscTest_cpu(buf.data(), n, x.data(), NULL, nS, pSet);
      polDiscTest_wrapv(buf.data(), n, y.data(), NULL, nS, pSet);
      for (int i = 0; i < n; i++) { tot++; fails += !nf_h_ok(&buf[i * 11]); d += x[i] != y[i]; }
    }
    printf("   %-12s %9lld polynomials, (H) fails on %7lld, flags differ on %lld\n", fam, tot, fails, d);
    bad |= d != 0;
  }

  // ---- 6 ----
  printf("\n6. cost per polynomial: blocks of 10240 (the application's buffer), family w8 where (H) holds, best of 5\n");
  {
    const int n = 10240, nb = 200;
    std::vector<i64> buf((size_t)n * 11 * nb); std::vector<signed char> fl(n);
    for (size_t i = 0; i < (size_t)n * nb; i++) do make("w8", &buf[i * 11]); while (!nf_h_ok(&buf[i * 11]));
    double best[5] = {1e9, 1e9, 1e9, 1e9, 1e9}; volatile i64 sink = 0;
    for (int rep = 0; rep < 5; rep++) {
      double t[6]; i64 s = 0;
      t[0] = now();
      for (int b = 0; b < nb; b++) pdf_filter_v3(&buf[(size_t)b * n * 11], n, fl.data());
      t[1] = now();
      for (int b = 0; b < nb; b++) pdf_filter_base(&buf[(size_t)b * n * 11], n, fl.data());
      t[2] = now();
      for (size_t i = 0; i < (size_t)n * nb; i++) s += nf_h_ok(&buf[i * 11]);
      t[3] = now();
      for (size_t i = 0; i < (size_t)n * nb; i++) s += nf_h_ok_fast(&buf[i * 11]);
      t[4] = now();
      for (int b = 0; b < nb; b++) {               // the bypass as the glue would do it
        i64 *p = &buf[(size_t)b * n * 11];
        pdf_filter_v3(p, n, fl.data());
        for (int i = 0; i < n; i++) if (!fl[i] && !nf_h_ok_fast(&p[i * 11])) fl[i] = 1;
        s += fl[0];
      }
      t[5] = now();
      sink = sink + s;
      for (int k = 0; k < 5; k++) best[k] = std::min(best[k], t[k + 1] - t[k]);
    }
    double N = (double)n * nb;
    printf("   filter, AVX2 build            %6.1f ns\n", 1e9 * best[0] / N);
    printf("   filter, baseline build        %6.1f ns\n", 1e9 * best[1] / N);
    printf("   check nf_h_ok                 %6.1f ns\n", 1e9 * best[2] / N);
    printf("   check nf_h_ok_fast            %6.1f ns   (%.0f%% of the AVX2 filter, %.0f%% of the baseline filter)\n",
           1e9 * best[3] / N, 100 * best[3] / best[0], 100 * best[3] / best[1]);
    printf("   AVX2 filter + bypass          %6.1f ns\n", 1e9 * best[4] / N);
    std::vector<char> up(2000);
    double t0 = now();
    polDiscTest_cpu(buf.data(), 2000, up.data(), NULL, nS, pSet);
    printf("   exact test (GMP)              %6.2f us per polynomial it runs on\n", 1e6 * (now() - t0) / 2000);
  }
  printf("\n%s\n", bad ? "FAILED" : "ok");
  return bad;
}
