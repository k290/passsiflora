// test_poldisc.cpp -- the link between the filter and the exact test (pdfilter_chains.tex, Proposition "What the exact test
// computes"): upstream's compiled polDiscTest_cpu against Definition "Passing" computed with PARI (poldisc, the primes of S
// removed, issquare), on families of polynomials; and the filter against both.
//   EXTRA_CXX=<GetDecicsSrc>/polDiscTest_cpu.cpp LIBS_EXTRA=-lgmp TEST=test_poldisc.cpp bash build-linux.sh
//   build/test_poldisc FAMILY COUNT SEED [-cand FILE] [-hits FILE] [-part i/N] [-S 2,5]
// FAMILY: rand3, rand100, rand1e6 (random monic, coefficients uniform in [-B, B]); pass (G(x^2) with G(0) = +-(1, 2, 5, 10)*m^2, and
// shifts x+t, |t| <= 3: polynomials built to pass); passbig (the same with |G_k| <= 10^6, m <= 10^5, |t| <= 2, so that upstream's int64 steps
// overflow on polynomials that pass); wrap (g + 2^62 at the constant term, g with |a_i| <= 3 passing and B9 = 10/cont(f') even: upstream's
// int64 steps see a0 only through B9^2*a0 with 4 | B9^2, so the wrapped numbers are those of g and upstream passes a polynomial that does
// not pass); zero (disc = 0: a square factor); cand (the captured candidates of the
// application, -cand FILE, with the shifts x+t, |t| <= 2); hits (the 83 real results, -hits FILE, with x -> -x and the shifts
// |t| <= COUNT).  -S gives the set S (one to three primes, default 2,5; for pass and passbig G(0) is then +-(a product of primes of S)
// times m^2).  -part i/N runs every N-th polynomial of cand or hits, starting with the i-th.
// For every polynomial: upstream's flag, the definition, whether one of upstream's int64 steps (the first two sub-resultant
// iterations, done in int64 before GMP takes over) would overflow (replayed in 128 bits), and the filter's flags (AVX2 and
// baseline builds, which must agree).
// Counted: disagreements of upstream and the definition (split by overflow), and filter flags 0 where upstream or the
// definition says 1 (where upstream says 1, split by overflow: outside the int64 range upstream can pass a polynomial that does
// not pass, as in wrap).  Exit status 1 if upstream disagrees with the definition without an overflow, or the filter rejects a
// polynomial that the definition passes, or one that upstream passes without an overflow.
#include <pari/pari.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <string>

extern "C" int polDiscTest_cpu(long long *polBuf, int numPolys, char *polGoodFlag, char *polyMask, int numP, int *pSet);
extern "C" void pdf_init_v3(int, const int*);
extern "C" void pdf_filter_v3(const long long*, int, signed char*);
extern "C" void pdf_init_base(int, const int*);
extern "C" void pdf_filter_base(const long long*, int, signed char*);

typedef long long i64;
typedef __int128 L;
static uint64_t rs;
static uint64_t rnd() { rs += 0x9E3779B97F4A7C15ull; uint64_t z = rs; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }
static i64 urand(i64 b) { return (i64)(rnd() % (uint64_t)(2 * b + 1)) - b; }

// The int64 part of upstream's polDiscTest_cpu (polBuf row c[0..10], leading first), replayed in 128 bits: true if any product,
// difference or quotient that upstream stores in an int64 leaves the int64 range.
static bool upstream_overflows(const i64 *c)
{
  bool ov = false;
  auto chk = [&](L v) { if (v > (L)INT64_MAX || v < (L)INT64_MIN) ov = true; return v; };
  L A[11], B[10], R[10];
  for (int col = 0; col < 11; col++) A[10 - col] = c[col];
  for (int k = 1; k <= 10; k++) B[k - 1] = chk((L)k * A[k]);
  int cont2 = 2; for (int k = 0; k < 10; k += 2) if ((B[k] & 1) != 0) { cont2 = 1; break; }
  int cont5 = 5; for (int k = 0; k < 9; k++) if (B[k] % 5 != 0) { cont5 = 1; break; }
  for (int k = 0; k < 10; k++) B[k] = B[k] / (cont2 * cont5);
  R[0] = chk(B[9] * A[0]);
  for (int k = 1; k <= 9; k++) R[k] = chk(chk(B[9] * A[k]) - B[k - 1]);
  for (int k = 0; k <= 8; k++) R[k] = chk(chk(B[9] * R[k]) - chk(R[9] * B[k]));
  int degR = 0; for (int k = 8; k > 0; k--) if (R[k] != 0) { degR = k; break; }
  if (ov) return true;
  if (degR == 8) {
    for (int k = 0; k < 10; k++) A[k] = B[k];
    for (int k = 0; k <= degR; k++) B[k] = R[k];
    R[0] = chk(B[8] * A[0]);
    for (int k = 1; k <= 8; k++) R[k] = chk(chk(B[8] * A[k]) - chk(A[9] * B[k - 1]));
    for (int k = 0; k <= 7; k++) R[k] = chk(chk(B[8] * R[k]) - chk(R[8] * B[k]));
    chk(B[9] * B[9]);
  } else if (degR == 7) {
    for (int k = 0; k < 10; k++) A[k] = B[k];
    for (int k = 0; k <= degR; k++) B[k] = R[k];
    R[0] = chk(B[7] * A[0]); R[1] = chk(B[7] * A[1]);
    for (int k = 2; k <= 8; k++) R[k] = chk(chk(B[7] * A[k]) - chk(A[9] * B[k - 2]));
    R[0] = chk(B[7] * R[0]);
    for (int k = 1; k <= 7; k++) R[k] = chk(chk(B[7] * R[k]) - chk(R[8] * B[k - 1]));
    for (int k = 0; k <= 6; k++) R[k] = chk(chk(B[7] * R[k]) - chk(R[7] * B[k]));
    chk(chk(B[9] * B[9]) * B[9]);
    chk(B[7] * B[7]);
  }
  return ov;
}

static std::vector<int> Sset = {2, 5};

// Definition "Passing": disc f != 0 and |disc f| with the factors of the primes of S removed a square
static int passes(const i64 *c)
{
  pari_sp av = avma;
  GEN f = cgetg(13, t_POL);
  f[1] = evalsigne(1) | evalvarn(0);
  for (int k = 0; k <= 10; k++) gel(f, 2 + k) = stoi(c[10 - k]);
  f = normalizepol(f);
  GEN D = ZX_disc(f);
  int r = 0;
  if (signe(D) != 0) {
    D = absi(D);
    for (int p : Sset) (void)Z_lvalrem(D, p, &D);
    r = Z_issquare(D) ? 1 : 0;
  }
  set_avma(av);
  return r;
}

static void shift(const i64 *c, i64 t, i64 *out)          // f(x + t), c leading first, exact in 128 bits
{
  L g[11];
  for (int k = 0; k <= 10; k++) g[k] = c[10 - k];
  for (int i = 0; i < 10; i++) for (int k = 9; k >= i; k--) g[k] = g[k] + (L)t * g[k + 1];
  for (int k = 0; k <= 10; k++) {
    if (g[k] > (L)INT64_MAX || g[k] < (L)INT64_MIN) { fprintf(stderr, "shift by %lld leaves the int64 range\n", t); exit(3); }
    out[10 - k] = (i64)g[k];
  }
}

static std::vector<std::vector<i64>> readPolys(const char *file)
{
  std::vector<std::vector<i64>> out;
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
    out.push_back(p);
  }
  fclose(f);
  return out;
}

int main(int argc, char **argv)
{
  if (argc < 4) { fprintf(stderr, "usage: test_poldisc FAMILY COUNT SEED [-cand FILE] [-hits FILE] [-part i/N] [-S 2,5]\n"); return 2; }
  std::string fam = argv[1]; i64 count = atoll(argv[2]); rs = strtoull(argv[3], 0, 10) * 0x9E3779B97F4A7C15ull + 7;
  const char *candf = 0, *hitsf = 0;
  i64 part_i = 1, part_n = 1;           // -part i/N: of the cand and hits polynomials only those with index = i-1 mod N
  for (int i = 4; i + 1 < argc; i++) {
    if (!strcmp(argv[i], "-cand")) candf = argv[i + 1];
    if (!strcmp(argv[i], "-hits")) hitsf = argv[i + 1];
    if (!strcmp(argv[i], "-part") && sscanf(argv[i + 1], "%lld/%lld", &part_i, &part_n) != 2) { fprintf(stderr, "-part i/N\n"); return 2; }
    if (!strcmp(argv[i], "-S")) {
      Sset.clear();
      for (const char *s = argv[i + 1]; *s; ) { Sset.push_back(atoi(s)); while (*s && *s != ',') s++; if (*s) s++; }
      if (Sset.empty() || Sset.size() > 3) { fprintf(stderr, "-S: one to three primes\n"); return 2; }
    }
  }
  pari_init(64000000, 0);
  int nS = (int)Sset.size(), pSet[3] = {0, 0, 0};
  for (int k = 0; k < nS; k++) pSet[k] = Sset[k];
  pdf_init_v3(nS, pSet); pdf_init_base(nS, pSet);
  std::vector<i64> units;                  // +- the products of subsets of S: the S-units used by pass and passbig
  for (int mask = 0; mask < (1 << nS); mask++) {
    i64 u = 1;
    for (int k = 0; k < nS; k++) if (mask >> k & 1) u *= Sset[k];
    units.push_back(u); units.push_back(-u);
  }
  if (nS != 2 || Sset[0] != 2 || Sset[1] != 5) { printf("S"); for (int p : Sset) printf(" %d", p); printf("\n"); }
  std::vector<std::vector<i64>> base;
  if (fam == "cand") { if (!candf) { fprintf(stderr, "cand needs -cand FILE\n"); return 2; } base = readPolys(candf); }
  if (fam == "hits") { if (!hitsf) { fprintf(stderr, "hits needs -hits FILE\n"); return 2; } base = readPolys(hitsf); }
  // the polynomials of the family, generated in blocks of up to 4096
  i64 made = 0, idx = part_i - 1, tshift = 0, sign = 1;
  i64 n_poly = 0, n_ov = 0, n_up = 0, n_def = 0, dis_noov = 0, dis_ov = 0, filt_bad_def = 0, filt_bad_up = 0, filt_bad_up_ov = 0, filt_rej = 0, base_disagree = 0;
  std::vector<i64> buf; std::vector<char> up; std::vector<signed char> fl, fb;
  i64 total = (fam == "cand") ? (i64)base.size() * 5 : (fam == "hits") ? (i64)base.size() * 2 * (2 * count + 1) : count;
  if (fam == "cand" || fam == "hits") total = total > part_i - 1 ? (total - part_i + part_n) / part_n : 0;
  while (made < total) {
    int n = (int)std::min<i64>(4096, total - made);
    buf.assign((size_t)n * 11, 0);
    for (int i = 0; i < n; i++) {
      i64 *c = &buf[(size_t)i * 11];
      c[0] = 1;
      if (fam == "rand3" || fam == "rand100" || fam == "rand1e6") {
        i64 b = fam == "rand3" ? 3 : fam == "rand100" ? 100 : 1000000;
        for (int k = 1; k <= 10; k++) c[k] = urand(b);
      } else if (fam == "pass" || fam == "passbig") {          // G(x^2) shifted
        bool big = fam == "passbig";
        i64 G[6]; G[5] = 1;
        for (int k = 1; k <= 4; k++) G[k] = urand(big ? 1000000 : 20);
        i64 m = 1 + (i64)(rnd() % (big ? 100000 : 30));
        G[0] = units[rnd() % units.size()] * m * m;
        i64 f[11]; for (int k = 0; k <= 10; k++) f[k] = 0;
        for (int k = 0; k <= 5; k++) f[10 - 2 * k] = G[k];             // f[col] is the coefficient of x^(10-col), G[k] of x^k
        shift(f, urand(big ? 2 : 3), c);               // |coefficients| <= 10^11 * 3^10 < 2^63 for passbig
      } else if (fam == "wrap") {
        for (;;) {
          for (int k = 1; k <= 10; k++) c[k] = urand(3);
          bool b9even = false;
          for (int k = 0; k < 10; k += 2) if ((((i64)(k + 1) * c[10 - (k + 1)]) & 1) == 1) b9even = true;
          signed char fl;
          if (!b9even) continue;
          pdf_filter_v3(c, 1, &fl);                       // the filter first: it rejects almost every one that does not pass
          if (!fl || !passes(c)) continue;
          c[10] += (c[10] < 0 ? 1 : -1) * ((i64)1 << 62);
          break;
        }
      } else if (fam == "zero") {          // (x^2 + a x + b)^2 * (x^6 + ...)
        i64 a = urand(20), b = urand(20), h[7]; h[0] = 1; for (int k = 1; k <= 6; k++) h[k] = urand(20);
        i64 q2[5] = {1, 2 * a, a * a + 2 * b, 2 * a * b, b * b};    // (x^2 + a x + b)^2, leading first
        for (int k = 0; k <= 10; k++) c[k] = 0;
        for (int u = 0; u <= 4; u++) for (int v = 0; v <= 6; v++) c[u + v] += q2[u] * h[v];
      } else if (fam == "cand") {
        const std::vector<i64> &p = base[(size_t)(idx % base.size())];
        i64 t = (idx / (i64)base.size()) - 2;      // t = -2 .. 2
        shift(p.data(), t, c); idx += part_n;
      } else if (fam == "hits") {
        // each hit, for x and -x, with the shifts t = -count .. count
        i64 per = 2 * (2 * count + 1);
        const std::vector<i64> &p = base[(size_t)(idx / per)];
        i64 r = idx % per; sign = r < (2 * count + 1) ? 1 : -1; tshift = (r % (2 * count + 1)) - count;
        i64 g[11]; for (int k = 0; k <= 10; k++) g[k] = (sign < 0 && (k % 2)) ? -p[k] : p[k];   // f(-x) for degree 10
        shift(g, tshift, c); idx += part_n;
      } else { fprintf(stderr, "unknown family %s\n", fam.c_str()); return 2; }
    }
    up.assign(n, 0); fl.assign(n, 0); fb.assign(n, 0);
    polDiscTest_cpu(buf.data(), n, up.data(), NULL, nS, pSet);
    pdf_filter_v3(buf.data(), n, fl.data());
    pdf_filter_base(buf.data(), n, fb.data());
    for (int i = 0; i < n; i++) {
      const i64 *c = &buf[(size_t)i * 11];
      int def = passes(c);
      bool ov = upstream_overflows(c);
      n_poly++; n_ov += ov; n_up += up[i] != 0; n_def += def;
      if ((up[i] != 0) != (def != 0)) {
        if (ov) dis_ov++;
        else { dis_noov++; if (dis_noov <= 5) { fprintf(stderr, "DISAGREE upstream %d definition %d:", up[i], def); for (int k = 0; k <= 10; k++) fprintf(stderr, " %lld", c[k]); fprintf(stderr, "\n"); } }
      }
      if (fl[i] != fb[i]) base_disagree++;
      if (fl[i] == 0) { filt_rej++; if (def) filt_bad_def++; if (up[i]) (ov ? filt_bad_up_ov : filt_bad_up)++; }
    }
    made += n;
  }
  // machine line: family polynomials upstream-passes definition-passes overflows disagree-without-overflow disagree-with-overflow
  // filter-rejects rejects-passing-by-definition rejects-passing-by-upstream (without / with overflow) AVX2-baseline-differences
  printf("POLDISC %s %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld\n", fam.c_str(), (long long)n_poly, (long long)n_up, (long long)n_def,
         (long long)n_ov, (long long)dis_noov, (long long)dis_ov, (long long)filt_rej, (long long)filt_bad_def, (long long)filt_bad_up,
         (long long)base_disagree, (long long)filt_bad_up_ov);
  return (dis_noov || filt_bad_def || filt_bad_up || base_disagree) ? 1 : 0;
}
