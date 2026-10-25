/* A stand-in for an OpenCL library whose FIRST call changes the floating-point rounding mode of the calling thread, as the
   proxy ("patched") OpenCL.dll of the 2026-10-20 experiments did.  PROXY_MODE=up|down|zero chooses the mode; any other value
   (or none) leaves it alone.  The mode is set at the first API call, as in the proxy: a value set in a library constructor
   (DllMain on Windows) does not survive, because the loader restores the thread's floating-point state afterwards.
   Only clGetPlatformIDs exists, and it reports no platform, which is all the experiments need.
     gcc -shared -fPIC -O2 -o libOpenCL.so.1 stub_opencl.c                                                                 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fenv.h>
#include <xmmintrin.h>

int clGetPlatformIDs(unsigned num_entries, void *platforms, unsigned *num_platforms)
{
  static int done = 0;
  if (!done) {
    done = 1;
    const char *m = getenv("PROXY_MODE");
    if (m && (!strcmp(m, "up") || !strcmp(m, "down") || !strcmp(m, "zero"))) {
      int up = !strcmp(m, "up"), down = !strcmp(m, "down");
      fesetround(up ? FE_UPWARD : down ? FE_DOWNWARD : FE_TOWARDZERO);                       /* x87 control word and MXCSR */
      _mm_setcsr((_mm_getcsr() & ~0x6000u) | (up ? 0x4000u : down ? 0x2000u : 0x6000u));   /* MXCSR, explicitly */
      /* only the rounding-control field (0 nearest, 1 down, 2 up, 3 toward zero) is printed: the rest of MXCSR holds sticky
         exception flags that depend on what the process computed before, so the whole register is not a stable line */
      fprintf(stderr, "[stub OpenCL: rounding mode now %s, MXCSR rounding control = %u]\n", m, (_mm_getcsr() >> 13) & 3u);
    }
  }
  (void)num_entries; (void)platforms;
  if (num_platforms) *num_platforms = 0;
  return -1;                                   /* CL_PLATFORM_NOT_FOUND_KHR */
}
