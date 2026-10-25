/* What does red(x,q,qi) = x - q*((x*qi + MAGIC) - MAGIC) return when the arithmetic is not IEEE double with round-to-nearest?
   (red is the integer reduction of src/pdfilter.f90 lines 162-165.)   gcc -O0 -o red_demo red_demo.c -lm && ./red_demo        */
#include <stdio.h>
#include <fenv.h>
#include <math.h>

#define MAGIC 6755399441055744.0                      /* 1.5 * 2^52 */

static double red_double(double x, double q, double qi)           /* SSE2 double, round to nearest: the intended behaviour */
{
  volatile double t = x * qi;  t = t + MAGIC;  t = t - MAGIC;
  return x - q * t;
}
static long double red_x87(double x, double q, double qi)         /* the same expression kept in 80-bit x87 registers */
{
  volatile long double t = (long double)x * qi;  t = t + (long double)MAGIC;  t = t - (long double)MAGIC;
  return (long double)x - (long double)q * t;
}
static double red_mode(double x, double q, double qi, int mode)   /* SSE2 double with another rounding mode */
{
  int old = fegetround();
  fesetround(mode);
  volatile double t = x * qi;  t = t + MAGIC;  t = t - MAGIC;
  double r = x - q * t;
  fesetround(old);
  return r;
}

int main(void)
{
  double q = 1009.0, qi = 1.0 / q;                 /* a small filter prime */
  double xs[] = { 123456.0, 5.0 * 1009.0, 7.0 * 1009.0, 200.0, 1008.0, 1009.0 * 100.0 + 504.0 };
  printf("q = %.0f, qi = fl(1/q) = %.20g  (qi %s 1/q)\n", q, qi, (long double)qi > 1.0L / 1009.0L ? ">" : "<");
  printf("%12s %14s %26s %12s %12s %12s   x mod q\n", "x", "round-nearest", "x87 (80-bit registers)", "round-up", "round-down", "round-zero");
  for (unsigned i = 0; i < sizeof xs / sizeof *xs; i++) {
    double x = xs[i];
    printf("%12.0f %14.0f %26.17Lg %12.0f %12.0f %12.0f   %ld\n", x, red_double(x, q, qi), red_x87(x, q, qi),
           red_mode(x, q, qi, FE_UPWARD), red_mode(x, q, qi, FE_DOWNWARD), red_mode(x, q, qi, FE_TOWARDZERO), (long)fmod(x, q));
  }
  /* how many of the 30 filter primes have fl(1/q) > 1/q?  (these turn an exact multiple of q into -q under round-up) */
  return 0;
}
