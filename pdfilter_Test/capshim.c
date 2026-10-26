/* Captures the polynomials that the application hands to the pre-filter, with the filter's verdict, for
   test_pdfilter_oracle -cand FILE.  Like rmshim.c it is linked in place of the filter's entry points (the filter compiled with
   them renamed pdf_filter_*_in, -DPDF_FILTER_NAME='pdf_filter_v3_in' and '..._base_in'); every call is passed on unchanged.
   With NF_CAPTURE=<file>, one polynomial in NF_CAPTURE_EVERY (default 1000, chosen by a hash of its coefficients, so the
   same polynomials are chosen in every run) is appended to <file>: its 11 coefficients, leading one first, and the verdict
   (1 kept, 0 rejected).  Without NF_CAPTURE the shim does nothing.  Link it with the application's objects as the fault guide
   links app_feed and cpu_good, e.g. (in its work directory W, filter objects compiled into cap/):
     g++ -O3 -o cpu_cap obj/main.o $APP cap/pdfilter_v3.o cap/pdfilter_base.o cap/capshim.o $LIBS
     NF_CAPTURE=$W/cand.txt NF_CAPTURE_EVERY=500 ./cpu_cap        (in a directory holding the work unit as "in") */
#include <stdio.h>
#include <stdlib.h>

typedef long long i64;
typedef void (*filterfn)(const i64*, int, signed char*);
extern void pdf_filter_v3_in(const i64*, int, signed char*);
extern void pdf_filter_base_in(const i64*, int, signed char*);

static FILE *cap;
static unsigned long long every = 0;      /* 0: not yet initialised */

static void capture(filterfn f, const i64 *p, int n, signed char *o)
{
  f(p, n, o);
  if (!every) {
    const char *file = getenv("NF_CAPTURE"), *e = getenv("NF_CAPTURE_EVERY");
    every = e ? strtoull(e, NULL, 10) : 1000;
    if (!every) every = 1;
    cap = file ? fopen(file, "a") : NULL;
  }
  if (!cap) return;
  for (int i = 0; i < n; i++) {
    unsigned long long h = 1469598103934665603ull;
    for (int k = 0; k < 11; k++) { h ^= (unsigned long long)p[i * 11 + k]; h *= 1099511628211ull; h ^= h >> 29; }
    if (h % every) continue;
    for (int k = 0; k < 11; k++) fprintf(cap, "%lld ", p[i * 11 + k]);
    fprintf(cap, "%d\n", o[i] ? 1 : 0);
  }
  fflush(cap);
}

void pdf_filter_v3(const i64 *p, int n, signed char *o)   { capture(pdf_filter_v3_in, p, n, o); }
void pdf_filter_base(const i64 *p, int n, signed char *o) { capture(pdf_filter_base_in, p, n, o); }
