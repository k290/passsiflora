/* Used by fault_checks.sh: calls the real filter (its entry points renamed pdf_filter_*_in) with the floating point
   rounding mode set to "up", only for the duration of the call, so that PARI's own arithmetic in the test program
   is not disturbed.  This imitates a machine whose rounding mode is wrong while the filter runs. */
#include <fenv.h>
#include <xmmintrin.h>

typedef long long i64;
typedef void (*filterfn)(const i64*, int, signed char*);
extern void pdf_filter_v3_in(const i64*, int, signed char*);
extern void pdf_filter_base_in(const i64*, int, signed char*);

static void run_up(filterfn f, const i64 *p, int n, signed char *o)
{
  unsigned csr = _mm_getcsr();
  int old = fegetround();
  fesetround(FE_UPWARD);
  _mm_setcsr((_mm_getcsr() & ~0x6000u) | 0x4000u);
  f(p, n, o);
  _mm_setcsr(csr);
  fesetround(old);
}

void pdf_filter_v3(const i64 *p, int n, signed char *o)   { run_up(pdf_filter_v3_in, p, n, o); }
void pdf_filter_base(const i64 *p, int n, signed char *o) { run_up(pdf_filter_base_in, p, n, o); }
