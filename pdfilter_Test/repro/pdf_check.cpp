// pdf_check: runs the pre-filter (pdfilter.f90) ALONE on polynomials read from standard input and prints its verdicts.
// Used for the comparisons that the application cannot make: the order in which the rounding mode is changed.
//   input : one monic degree-10 polynomial per line, 11 integers, leading coefficient first
//   usage : pdf_check [alone|block] [up|down|zero] [ocl|ocl-late]
//             alone     one filter call per polynomial (a block of 1)          block  all polynomials in ONE call
//             up|down|zero  set the rounding mode around every filter call; the filter was initialised in round-to-nearest
//                       (the order of rmshim.c in fault_checks.sh)
//             ocl       call clGetPlatformIDs of $OCL_LIB (default ./libOpenCL.so.1) BEFORE pdf_init, as the application does
//             ocl-late  the same call AFTER pdf_init
//   output: "<index> <flag> <polynomial>" (flag 1 kept, 0 REJECTED); stderr: "<n> polynomials, <r> rejected (...)"
//   env   : NF_FORCE_BASE=1 (baseline build), NF_GENERAL_ALL=1, NF_PDFDEBUG=1 as for the application
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <fenv.h>
#include <xmmintrin.h>
#include <dlfcn.h>

extern "C" void pdf_init_v3(int, const int *);
extern "C" void pdf_filter_v3(const long long *, int, signed char *);
extern "C" void pdf_init_base(int, const int *);
extern "C" void pdf_filter_base(const long long *, int, signed char *);

static void (*filter)(const long long *, int, signed char *);

static void call(const long long *p, int n, signed char *flag, const char *mode)
{
  if (!mode) { filter(p, n, flag); return; }
  unsigned csr = _mm_getcsr();
  int old = fegetround();
  bool up = !strcmp(mode, "up"), down = !strcmp(mode, "down");
  fesetround(up ? FE_UPWARD : down ? FE_DOWNWARD : FE_TOWARDZERO);
  _mm_setcsr((_mm_getcsr() & ~0x6000u) | (up ? 0x4000u : down ? 0x2000u : 0x6000u));
  filter(p, n, flag);
  _mm_setcsr(csr);
  fesetround(old);
}

static void call_opencl()
{
  const char *lib = getenv("OCL_LIB");
  void *h = dlopen(lib ? lib : "./libOpenCL.so.1", RTLD_NOW);
  typedef int (*Fn)(unsigned, void *, unsigned *);
  Fn f = h ? (Fn)dlsym(h, "clGetPlatformIDs") : 0;
  if (!f) { fprintf(stderr, "cannot load clGetPlatformIDs from the OpenCL library\n"); exit(2); }
  unsigned np = 0;
  f(0, 0, &np);
}

static std::string show(const long long *c)
{
  std::string s = "x^10";
  char buf[64];
  for (int k = 1; k <= 10; k++) {
    if (c[k] == 0) continue;
    snprintf(buf, sizeof buf, " %c %lld", c[k] < 0 ? '-' : '+', c[k] < 0 ? -c[k] : c[k]);
    s += buf;
    if (k < 10) { snprintf(buf, sizeof buf, "*x^%d", 10 - k); s += buf; }
  }
  return s;
}

int main(int argc, char **argv)
{
  const char *ctx = "alone", *mode = 0;
  bool ocl = false, ocl_late = false;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "alone") || !strcmp(argv[i], "block")) ctx = argv[i];
    else if (!strcmp(argv[i], "up") || !strcmp(argv[i], "down") || !strcmp(argv[i], "zero")) mode = argv[i];
    else if (!strcmp(argv[i], "ocl")) ocl = true;
    else if (!strcmp(argv[i], "ocl-late")) ocl_late = true;
    else { fprintf(stderr, "usage: %s [alone|block] [up|down|zero] [ocl|ocl-late] < polynomials\n", argv[0]); return 2; }
  }
  void (*init)(int, const int *);
  const char *fb = getenv("NF_FORCE_BASE");
  __builtin_cpu_init();
  if (!(fb && *fb == '1') && __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")) {
    init = pdf_init_v3; filter = pdf_filter_v3;
  } else {
    init = pdf_init_base; filter = pdf_filter_base;
  }
  if (ocl) call_opencl();
  int pSet[2] = {2, 5};
  init(2, pSet);
  if (ocl_late) call_opencl();

  std::vector<long long> polys;
  long long v;
  while (scanf("%lld", &v) == 1) polys.push_back(v);
  if (polys.size() % 11) { fprintf(stderr, "input must be a multiple of 11 integers\n"); return 2; }
  int n = (int)(polys.size() / 11);
  for (int i = 0; i < n; i++) if (polys[i * 11] != 1) { fprintf(stderr, "polynomial %d is not monic\n", i + 1); return 2; }

  std::vector<signed char> flag(n, 1);
  if (!strcmp(ctx, "alone")) { for (int i = 0; i < n; i++) call(&polys[i * 11], 1, &flag[i], mode); }
  else call(polys.data(), n, flag.data(), mode);

  int rejected = 0;
  for (int i = 0; i < n; i++) {
    printf("%d %d %s\n", i + 1, (int)flag[i], show(&polys[i * 11]).c_str());
    if (!flag[i]) rejected++;
  }
  fprintf(stderr, "%d polynomials, %d rejected\n", n, rejected);
  return 0;
}
