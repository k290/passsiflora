/* limit_probe.c -- what happens in the chains of the pre-filter when a prime is larger than the proved limits
   (docs/whitepaper.tex, Section 4.4 "The size limits, and the guard").  A one-lane model of load_block + chain_block +
   chain_general (16-bit values, mont with R = 2^16; the loading's sredc with R = 2^32), written from src/pdfilter.f90 with
   the same integer semantics (conversions to a narrower kind and the 16-bit subtraction wrap modulo 2^k, as gfortran does).
   It is a model, not the filter: the filter-level tests are fault_checks.sh rows SF..S15; this program measures how close
   the values come to the limits and searches for inputs that change a verdict.
     gcc -O2 -o limit_probe limit_probe.c
     ./limit_probe Q SEED N           N random polynomials (residues in (-q/2, q/2)) at the prime Q
     ./limit_probe Q SEED N climb ITERS [x|u]   N random starts, each improved by ITERS one-coefficient changes that push
                                      up the largest |x| / (q R/2) (x, the default) or the largest |u| / q and
                                      |difference| / 2^15 (u): a search for directed inputs, not a sample
   Prints, for the N polynomials: those with a reduction result outside (-q, q), 16-bit or 64-bit overflows, zeros that a
   zero test misses, verdicts that differ from the exact ones (exact: Euclid with inverses mod Q, Q must be prime), the
   largest |x| / (q R/2) met by a reduction (the lemmas need < 1), its histogram, and the first polynomials whose verdict
   differs (11 coefficients, leading one first: directed inputs for the filter). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef int16_t i2; typedef int32_t i4; typedef int64_t i8; typedef __int128 i128;

static i8 Q;
static i2 QI2, ONE2, NEGR2; static i4 QL, QIL, LD[3];
static long ev_range, ev_ovf, ev_zero;
static double polymax, polyu;     /* largest |x| / (q R/2); largest |u| / q and |difference| / 2^15 (an event iff >= 1) */
static uint64_t rs;
static uint64_t rnd(void) { rs += 0x9E3779B97F4A7C15ull; uint64_t z = rs; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }
static void ratio(i128 x, double half) { double r = (double)(x < 0 ? -x : x) / half; if (r > polymax) polymax = r; }
static void uratio(i128 v, double lim) { double r = (double)(v < 0 ? -v : v) / lim; if (r > polyu) polyu = r; }

/* ---- mont exactly as the Fortran (lo, m, h1, h2), 16-bit subtraction */
static i2 mont(i2 a, i2 b) {
  i4 x = (i4)a * b; ratio(x, (double)Q * 32768.0);
  i2 lo = (i2)x, m = (i2)((i4)lo * QI2), h1 = (i2)(x >> 16), h2 = (i2)(((i4)m * (i4)Q) >> 16);
  i2 u = (i2)(h1 - h2); uratio(u, (double)Q); if (u >= Q || u <= -Q) ev_range++; return u;
}
static i2 sub2(i2 a, i2 b) { i4 d = (i4)a - b; uratio(d, 32768.0); if (d > 32767 || d < -32768) ev_ovf++; return (i2)d; }
/* loading: 32-bit sredc of p0 d0 + p1 d1 + p2 d2 */
static i4 sredc32(i8 x) { i4 m = (i4)((i8)(i4)x * QIL); return (i4)((x - (i8)m * QL) >> 32); }
static i4 load(i8 c) { i4 p0 = (i4)(c & ((1 << 21) - 1)), p1 = (i4)((c >> 21) & ((1 << 21) - 1)), p2 = (i4)(c >> 42);
  return sredc32((i8)p0 * LD[0] + (i8)p1 * LD[1] + (i8)p2 * LD[2]); }

/* ---- exact answer: chi(disc) and normality, by Euclid with inverses mod Q.  Q < 2^31 (it fits i4), so every product of two
   residues is below 2^62 and the arithmetic is done in i8. */
static i8 pw(i8 b, i8 e, i8 m) { i8 r = 1, x = ((b % m) + m) % m; while (e) { if (e & 1) r = r * x % m; x = x * x % m; e >>= 1; } return r; }
static void exact(const i8 *c, int *chi, int *normal) {
  i8 a[11], b[11]; int n = 10, m = 9; i8 q = Q;
  for (int i = 0; i <= 10; i++) { a[i] = c[10 - i] % q; if (a[i] < 0) a[i] += q; }
  for (int i = 0; i <= 9; i++) b[i] = (i + 1) * a[i + 1] % q;
  i8 res = 1; int neg = 0; *normal = 1; i8 *A = a, *B = b;
  for (;;) {
    if (m == 0) { res = res * pw(B[0], n, q) % q; break; }
    i8 il = pw(B[m], q - 2, q);
    for (int i = n; i >= m; i--) { i8 t = A[i] * il % q; if (t) for (int k = 0; k <= m; k++) A[i - m + k] = ((A[i - m + k] - t * B[k]) % q + q) % q; }
    int k = m - 1; while (k >= 0 && A[k] == 0) k--;
    if (k != m - 1 && m - 1 >= 1) *normal = 0;
    if (k < 0) { res = 0; break; }
    if ((n & 1) && (m & 1)) neg = !neg;
    res = res * pw(B[m], n - k, q) % q;
    i8 *t = A; A = B; B = t; n = m; m = k;
  }
  if (neg) res = (q - res) % q;
  i8 d = (q - res) % q;
  *chi = d == 0 ? 0 : (pw(d, (q - 1) / 2, q) == 1 ? 1 : -1);
}

/* ---- the chains, one lane; values held in 64-bit variables but converted exactly as the code's kinds */
static i8 MUL(i8 a, i8 b) { return (i8)mont((i2)a, (i2)b); }
static i8 DIFF(i8 a, i8 b, i8 c, i8 d) { return (i8)sub2(mont((i2)a, (i2)b), mont((i2)c, (i2)d)); }
static int iszero(i8 v) { return v == 0 || v == Q || v == -Q; }
static int euler(i8 base) {
  i8 acc = ONE2, negr = NEGR2, e = (Q - 1) / 2;
  while (e > 0) { if (e & 1) acc = MUL(acc, base); base = MUL(base, base); e >>= 1; }
  return acc == negr || acc == negr - Q;
}
static int chain_block(const i8 *a, int *ok) {
  i8 p[3][11]; memset(p, 0, sizeof p); int ia = 0, ib = 1, ir = 2;
  for (int k = 0; k <= 10; k++) p[ia][k] = a[k];
  for (int k = 1; k <= 10; k++) p[ib][k - 1] = MUL(a[k], k);
  *ok = 1;
  for (int m = 9; m >= 1; m--) {
    i8 lc = p[ib][m], top = p[ia][m + 1];
    p[ir][0] = MUL(lc, p[ia][0]);
    for (int k = 1; k <= m; k++) p[ir][k] = DIFF(lc, p[ia][k], top, p[ib][k - 1]);
    i8 t = p[ir][m];
    for (int k = 0; k <= m - 1; k++) p[ir][k] = DIFF(lc, p[ir][k], t, p[ib][k]);
    if (m >= 2) { i8 v = p[ir][m - 1]; if (iszero(v)) *ok = 0; else if (v % Q == 0) ev_zero++; }
    int tmp = ia; ia = ib; ib = ir; ir = tmp;
  }
  return !(*ok && euler(p[ib][0]));
}
static int chain_general(const i8 *a) {
  i8 p[3][12]; memset(p, 0, sizeof p); int iu = 0, iv = 1, iw = 2;
  for (int i = 0; i <= 10; i++) p[iu][i] = a[10 - i];
  for (int i = 0; i <= 9; i++) p[iv][i] = MUL(10 - i, a[10 - i]);
  int du = 10, dv = 9, done = 0, resz = 0; i8 mm = 1;
  for (int it = 1; it <= 40 && !done; it++) {
    int imax = du; i8 lu = p[iu][0], lv = p[iv][0];
    for (int i = 0; i <= imax - 1; i++) p[iw][i] = DIFF(lv, p[iu][i + 1], lu, p[iv][i + 1]);
    for (int i = imax; i <= 11; i++) p[iw][i] = 0;
    int dw = du - 1;
    while (dw >= 0 && iszero(p[iw][0])) { for (int i = 0; i <= imax - 1; i++) p[iw][i] = p[iw][i + 1]; p[iw][imax] = 0; dw--; }
    if (dw >= 0 && p[iw][0] % Q == 0) ev_zero++;
    if (dw < 0) { resz = 1; done = 1; break; }
    int ev = abs(dw - du + dv); i8 tmp = MUL(mm, lv); if (ev & 1) mm = tmp;
    int sw = dw < dv;
    { int t = iu; iu = iw; iw = t; du = dw; if (sw) { t = iu; iu = iv; iv = t; t = du; du = dv; dv = t; } }
    tmp = MUL(mm, p[iv][0]); if (dv == 0 && (du & 1)) mm = tmp;
    if (dv == 0) done = 1;
  }
  return !(done && !resz && euler(mm));
}
static i4 qinv32(i8 q) { uint64_t x = 1; for (int i = 0; i < 6; i++) x *= 2 - (uint64_t)q * x; return (i4)(uint32_t)x; }

/* one polynomial (c[0] = 1, c[k] the coefficient of x^(10-k)): loading and both chains; returns the objective of climb
   (OBJU: the largest |u| / q and |difference| / 2^15; else the largest |x| / (q R/2)), and whether a
   reduction left (-q, q) or a subtraction overflowed */
static double allmax, allu;         /* the largest |x| / (q R/2) and |u| / q over every polynomial evaluated by climb */
static int OBJU;
static double run_poly(const i8 *c, int *event, int *ok, int *kb, int *kg) {
  i8 a[11]; for (int k = 0; k <= 9; k++) a[k] = (i8)(i2)load(c[10 - k]); a[10] = 1;
  long r0 = ev_range, o0 = ev_ovf; polymax = 0; polyu = 0;
  *kb = chain_block(a, ok); *kg = chain_general(a);
  *event = ev_range > r0 || ev_ovf > o0;
  if (polymax > allmax) allmax = polymax;
  if (polyu > allu) allu = polyu;
  return OBJU ? polyu : polymax;
}

/* climb: from each of N random polynomials, ITERS steps that change one coefficient (to a random residue, or by a small
   amount) and keep the change when the objective (see run_poly) does not decrease.  Every
   polynomial met with an event is compared with the exact answer; those whose verdict differs are printed (directed
   inputs). */
static int climb(long N, long iters) {
  double best = 0; long nev = 0, ndiff = 0, nevals = 0, shown = 0; i8 bestc[11];
  for (long it = 0; it < N; it++) {
    i8 c[11], d[11]; c[0] = 1;
    for (int k = 1; k <= 10; k++) c[k] = (i8)(rnd() % (uint64_t)Q) - Q / 2;
    int ev, ok, kb, kg; double cur = run_poly(c, &ev, &ok, &kb, &kg); nevals++;
    for (long s = 0; s <= iters; s++) {
      if (ev) {
        nev++;
        int chi, normal; exact(c, &chi, &normal); int want = chi != -1;
        if ((ok && kb != want) || kg != want) {
          ndiff++;
          if (shown < 20) {
            shown++; printf("directed input (exact chi %d; chain_block ok %d keep %d; chain_general keep %d; |u|/q %.4f):", chi, ok, kb, kg, cur);
            for (int k = 0; k <= 10; k++) printf(" %lld", (long long)c[k]);
            printf("\n");
          }
        }
        break;                                 /* an event: this start is done */
      }
      if (s == iters) break;
      memcpy(d, c, sizeof d);
      int k = 1 + (int)(rnd() % 10);
      if (rnd() & 1) d[k] = (i8)(rnd() % (uint64_t)Q) - Q / 2;
      else { d[k] += (i8)(rnd() % 201) - 100; if (d[k] > Q / 2) d[k] -= Q; if (d[k] < -Q / 2) d[k] += Q; }
      int ev2, ok2, kb2, kg2; double v = run_poly(d, &ev2, &ok2, &kb2, &kg2); nevals++;
      if (v >= cur || ev2) { memcpy(c, d, sizeof c); cur = v; ev = ev2; ok = ok2; kb = kb2; kg = kg2; }
    }
    if (cur > best) { best = cur; memcpy(bestc, c, sizeof bestc); }
  }
  printf("q=%lld climb (objective %s): %ld starts, %ld steps each, %ld polynomials evaluated; largest |u| / q %.4f, "
         "largest |x| / (q R/2) %.4f; starts that reached an event %ld, of them with a different verdict %ld\n",
         (long long)Q, OBJU ? "u" : "x", N, iters, nevals, allu, allmax, nev, ndiff);
  printf("best polynomial:"); for (int k = 0; k <= 10; k++) printf(" %lld", (long long)bestc[k]); printf("\n");
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: limit_probe Q SEED N [climb ITERS [x|u]]\n"); return 2; }
  Q = atoll(argv[1]); rs = strtoull(argv[2], 0, 10); long N = atol(argv[3]);
  if (Q <= 5 || Q >= 32768 || !(Q & 1)) { fprintf(stderr, "Q must be odd, > 5 and below 2^15\n"); return 2; }
  i8 R16 = 1 << 16, R32 = (i8)1 << 32;
  QI2 = (i2)(qinv32(Q) & 0xFFFF); ONE2 = (i2)(R16 % Q); NEGR2 = (i2)((Q - R16 % Q) % Q);
  QL = (i4)Q; QIL = qinv32(Q);
  i8 r = R32 % Q; for (int i = 0; i < 3; i++) { LD[i] = (i4)r; r = (i8)((i128)r * ((i8)1 << 21) % Q); }
  if (argc >= 6 && !strcmp(argv[4], "climb")) { OBJU = argc >= 7 && argv[6][0] == 'u'; return climb(N, atol(argv[5])); }
  long nrange = 0, nerr_b = 0, nerr_g = 0, nok = 0, shown = 0, hist[60] = {0}; double maxr = 0;
  for (long it = 0; it < N; it++) {
    i8 c[11]; c[0] = 1;
    for (int k = 1; k <= 10; k++) c[k] = (i8)(rnd() % (uint64_t)Q) - Q / 2;
    i8 a[11]; for (int k = 0; k <= 9; k++) a[k] = (i8)(i2)load(c[10 - k]); a[10] = 1;
    int chi, normal; exact(c, &chi, &normal);
    long r0 = ev_range; polymax = 0;
    int ok, kb = chain_block(a, &ok), kg = chain_general(a);
    if (polymax > maxr) maxr = polymax;
    { int bin = (int)(polymax * 20); if (bin > 59) bin = 59; hist[bin]++; }
    if (ev_range > r0) nrange++;
    int want = chi != -1;
    if (ok != normal) nok++;
    if (ok && kb != want) nerr_b++;
    if (kg != want) nerr_g++;
    if (((ok && kb != want) || kg != want) && shown < 3) {
      shown++; printf("verdict differs (exact chi %d, normal %d; chain_block ok %d keep %d; chain_general keep %d):", chi, normal, ok, kb, kg);
      for (int k = 0; k <= 10; k++) printf(" %lld", (long long)c[k]);
      printf("\n");
    }
  }
  printf("q=%lld: %ld polynomials; with a result outside (-q,q): %ld (%ld results); overflows %ld; zeros missed %ld; "
         "ok differs from normality %ld; verdict differs: chain_block %ld, chain_general %ld\n",
         (long long)Q, N, nrange, ev_range, ev_ovf, ev_zero, nok, nerr_b, nerr_g);
  printf("largest |x| / (q R/2): %.3f; per-polynomial maximum, bins of 0.05:", maxr);
  for (int i = 0; i < 60; i++) if (hist[i]) printf(" %.2f:%ld", i * 0.05, hist[i]);
  printf("\n");
  return 0;
}
