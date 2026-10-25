/* signtest: the round-up lemma at all 30 filter primes.
   Under round-up, red(k*q) = x - q*((x*qi + MAGIC) - MAGIC) (pdfilter.f90 lines 162-172) returns -q instead of 0 exactly when
   (k > 0 and qi > 1/q) or (k < 0 and qi < 1/q).  For every filter prime (chosen as pdf_init does) and both reciprocals the
   filter can hold - qi rounded to nearest (pdf_init in round-to-nearest) and qi rounded up (pdf_init after a library has set
   round-up) - this program decides EXACTLY whether qi is above 1/q (integer arithmetic on the mantissa), then reduces exact
   multiples k*q under round-up, with separate operations and with a fused multiply-add for x*qi + MAGIC, and counts the
   results that differ from the lemma.  Stages 1-6 use single precision (q < 2^10, |x| <= 2^19), stages 7-30 double precision
   (q < 2^25, |x| < 2^51), as the filter does.
   Then the rest of Lemma "Directed rounding in red" of pdfilter_chains.tex for random |x| <= 2^(p-3): (a) red(x) = r - q under round-up
   when q does not divide x; (b) red(x) in [0, q], congruent to x, and = r for non-multiples, under round-down and toward zero.  Finally
   Lemma "red in 80-bit registers": red evaluated in long double (the x87 format on x86-64 Linux), and the two-step load
   red(red(hi)*s + lo) of load_block, against L_q (the symmetric residue of 2048*x modulo q, divided by 2048).  Exit status 1 on a difference.
     gcc -O0 -o signtest signtest.c -lm && ./signtest                                                                       */
#include <stdio.h>
#include <fenv.h>
#include <math.h>

#define MAGIC  6755399441055744.0       /* 1.5*2^52 */
#define MAGICS 12582912.0f              /* 1.5*2^23 */

static int is_prime(long n) { if (n < 2) return 0; if (n % 2 == 0) return n == 2; for (long d = 3; d * d <= n; d += 2) if (n % d == 0) return 0; return 1; }
static long powmod(long b, long e, long m) { long r = 1; b %= m; while (e > 0) { if (e & 1) r = (__int128)r * b % m; b = (__int128)b * b % m; e >>= 1; } return r; }
static int usable(long q) { return q % 4 == 1 && is_prime(q) && powmod(2, (q - 1) / 2, q) == 1 && powmod(5, (q - 1) / 2, q) == 1; }

/* exact test: is the binary number v (> 0) above 1/q ?  v = m * 2^e with an integer mantissa m: compare q*m with 2^-e */
static int above(double v, long q)
{
  int e; double fr = frexp(v, &e);                  /* v = fr * 2^e, 0.5 <= fr < 1 */
  unsigned long long m = (unsigned long long)ldexp(fr, 53); e -= 53;   /* v = m * 2^e, exact */
  __int128 lhs = (__int128)q * m, rhs = (__int128)1 << (-e);
  return lhs > rhs;                                  /* never equal: 1/q is not a binary fraction (q odd > 1) */
}

static volatile double vd;                            /* volatile: one rounding per operation, no contraction */
static volatile float vf;
static double red_up_d(double x, double q, double qi, int fused)
{
  fesetround(FE_UPWARD);
  double t;
  if (fused) t = fma(x, qi, MAGIC); else { vd = x * qi; vd = vd + MAGIC; t = vd; }
  vd = t - MAGIC; t = vd;
  double r = x - q * t;
  fesetround(FE_TONEAREST);
  return r;
}
static float red_up_f(float x, float q, float qi, int fused)
{
  fesetround(FE_UPWARD);
  float t;
  if (fused) t = fmaf(x, qi, MAGICS); else { vf = x * qi; vf = vf + MAGICS; t = vf; }
  vf = t - MAGICS; t = vf;
  float r = x - q * t;
  fesetround(FE_TONEAREST);
  return r;
}

int main(void)
{
  long primes[30]; int np = 0;
  for (long q = 1023; np < 6 && q > 256; q -= 2) if (usable(q)) primes[np++] = q;
  for (long q = (1L << 25) - 1; np < 30 && q > (1L << 24); q -= 2) if (usable(q)) primes[np++] = q;
  unsigned long lcg = 12345;
  long tested = 0, wrong = 0; int above_near = 0;
  char row[30][96];                                    /* printed in two columns: stages 1-15 | 16-30 */
  for (int s = 0; s < 30; s++) {
    long q = primes[s]; int single = s < 6;
    long kmax = single ? (1L << 19) / q : ((1L << 51) - 1) / q;
    double qn, qu; int an, au;
    if (single) { volatile float a = (float)q; fesetround(FE_TONEAREST); vf = 1.0f / a; qn = vf; fesetround(FE_UPWARD); vf = 1.0f / a; qu = vf; fesetround(FE_TONEAREST); }
    else        { volatile double a = (double)q; fesetround(FE_TONEAREST); vd = 1.0 / a; qn = vd; fesetround(FE_UPWARD); vd = 1.0 / a; qu = vd; fesetround(FE_TONEAREST); }
    an = above(qn, q); au = above(qu, q); above_near += an;
    long n = 0, bad = 0;
    for (int which = 0; which < 2; which++) {
      double qi = which ? qu : qn; int ab = which ? au : an;
      for (long i = 0; i < 4000; i++) {
        long k;
        if (i < 1000 || kmax <= 4000) k = (i % 1000) + 1;
        else { lcg = lcg * 6364136223846793005UL + 1442695040888963407UL; k = 1 + (long)((lcg >> 11) % (unsigned long)kmax); }
        if (k > kmax) continue;
        for (int sg = -1; sg <= 1; sg += 2) {
          long kk = sg * k;
          double expect = ((kk > 0 && ab) || (kk < 0 && !ab)) ? -(double)q : 0.0;
          for (int fused = 0; fused < 2; fused++) {
            double r = single ? (double)red_up_f((float)(kk * q), (float)q, (float)qi, fused) : red_up_d((double)(kk * q), (double)q, qi, fused);
            n++; if (r != expect) bad++;
          }
        }
      }
    }
    tested += n; wrong += bad;
    snprintf(row[s], sizeof row[s], "%5d  %-8ld  %-6s  %-7s  %-5s  %6ld  %ld", s + 1, q, single ? "single" : "double",
             an ? "above" : "below", au ? "above" : "below", n, bad);
  }
  char head[96];
  snprintf(head, sizeof head, "%5s  %-8s  %-6s  %-7s  %-5s  %6s  %s", "stage", "q", "prec", "qi.near", "qi.up", "mult.", "differ");
  printf("%-50s  |  %s\n", head, head);
  for (int s = 0; s < 15; s++) printf("%-50s  |  %s\n", row[s], row[s + 15]);
  printf("30 primes, %ld reductions of exact multiples under round-up, %ld differ from the lemma\n", tested, wrong);
  printf("qi rounded to nearest is above 1/q at %d of 30 stages, qi rounded up at 30\n", above_near);

  /* the rest of the lemma: (a) for integers that are not multiples, (b) rounding down and toward zero, random |x| <= 2^(p-3) */
  long na = 0, nb = 0, bada = 0, badb = 0;
  for (int s = 0; s < 30; s++) {
    long q = primes[s]; int single = s < 6;
    long xmax = single ? (1L << 21) : (1L << 50);
    double qi;
    if (single) { volatile float a = (float)q; vf = 1.0f / a; qi = vf; } else { volatile double a = (double)q; vd = 1.0 / a; qi = vd; }
    for (int i = 0; i < 20000; i++) {
      lcg = lcg * 6364136223846793005UL + 1442695040888963407UL;
      long x = (long)((lcg >> 11) % (unsigned long)(2 * xmax + 1)) - xmax;
      if (i % 4 == 0) x = (x / q) * q;                                      /* a quarter exact multiples */
      long r = ((x % q) + q) % q;                                           /* 0 <= r < q */
      int fused = i & 1;
      double up = single ? (double)red_up_f((float)x, (float)q, (float)qi, fused) : red_up_d((double)x, (double)q, qi, fused);
      if (r != 0) { na++; if (up != (double)(r - q)) bada++; }
      for (int mode = 0; mode < 2; mode++) {                                /* 0: down, 1: toward zero */
        int rm = mode ? FE_TOWARDZERO : FE_DOWNWARD;
        double v;
        fesetround(rm);
        if (single) { float t; if (fused) t = fmaf((float)x, (float)qi, MAGICS); else { vf = (float)x * (float)qi; vf = vf + MAGICS; t = vf; }
                      vf = t - MAGICS; t = vf; v = (float)x - (float)q * t; }
        else        { double t; if (fused) t = fma((double)x, qi, MAGIC); else { vd = (double)x * qi; vd = vd + MAGIC; t = vd; }
                      vd = t - MAGIC; t = vd; v = (double)x - (double)q * t; }
        fesetround(FE_TONEAREST);
        nb++;
        if (v < 0 || v > q || v != (double)(long)v || ((long)v - x) % q != 0 || (r != 0 && v != (double)r)) badb++;
      }
    }
  }
  printf("part (a), not multiples: %ld reductions under round-up, %ld differ\n", na, bada);
  printf("part (b): %ld reductions under round-down and toward zero, %ld differ\n", nb, badb);

  /* Lemma (red in 80-bit registers): every operation in long double (x87, 64-bit significand, round to nearest), result stored as a double,
     for x = X/2048 with |x| < 2^39, and the two-step load red(red(hi)*s + lo) of a coefficient c, |c| < 2^51.  Expected: L_q(c), which is the
     symmetric residue of 2048*c (resp. of X) modulo q, divided by 2048. */
  long nx = 0, badx = 0;
  for (int s = 0; s < 30; s++) {
    long q = primes[s];
    volatile double qa = (double)q; volatile double qi = 1.0 / qa;          /* the double reciprocal of pdf_init (load_block uses doubles) */
    long double M = 6755399441055744.0L, Q = (long double)q, QI = (long double)qi;
    long s25 = (1L << 25) % q;
    for (int i = 0; i < 20000; i++) {
      lcg = lcg * 6364136223846793005UL + 1442695040888963407UL;
      long c = (long)((lcg >> 11) % (1UL << 52)) - (1L << 51) + 1;         /* |c| < 2^51 */
      long X = c >> 1;                                                      /* |X| < 2^50, x = X/2048 */
      long double x = (long double)X / 2048.0L;
      long double t = x * QI, u = t + M, w = u - M, red1 = x - Q * w;
      long qr = (long)(((__int128)X % q + q) % q); if (qr > q / 2) qr -= q;
      nx++; if ((double)red1 != (double)qr / 2048.0) badx++;
      long hi = c >> 25, lo = c & ((1L << 25) - 1);
      long double h = (long double)hi;
      t = h * QI; u = t + M; w = u - M; long double rh = h - Q * w;
      long double x2 = rh * (long double)s25 + (long double)lo;
      t = x2 * QI; u = t + M; w = u - M; long double load = x2 - Q * w;
      long cr = (long)((((__int128)c * 2048) % q + q) % q); if (cr > q / 2) cr -= q;
      nx++; if ((double)load != (double)cr / 2048.0) badx++;
    }
  }
  printf("x87 lemma: %ld reductions and loads in 80-bit long double, %ld differ from L_q\n", nx, badx);
  return wrong != 0 || bada != 0 || badb != 0 || badx != 0;
}
