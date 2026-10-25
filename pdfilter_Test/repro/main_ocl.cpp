// main_ocl: the BOINC-free main of the CPU application (stubs/main.cpp: the same call into TgtMartinet(), the same "StatVec"
// line) that first calls clGetPlatformIDs of the library named by $OCL_LIB, as nfGpuInit does in the GPU application.  A
// library that changes the rounding mode at that call then leaves the thread that runs the CPU pre-filter in that mode.
#include <dlfcn.h>
#include <cstdio>
#include <cstdlib>
#include <xmmintrin.h>
#include "GetDecics.h"
#include "TgtMartinet.h"

int numBlocks = 10240;
int threadsPerBlock = 1;
int polyBufferSize = 10240;
extern int nfArgc;
extern char **nfArgv;

int do_checkpoint(pari_long*, pari_long*) { return 0; }

static void call_opencl(const char *lib)
{
  void *h = dlopen(lib, RTLD_NOW);
  typedef int (*Fn)(unsigned, void *, unsigned *);
  Fn f = h ? (Fn)dlsym(h, "clGetPlatformIDs") : 0;
  if (!f) { fprintf(stderr, "cannot load clGetPlatformIDs from %s\n", lib); exit(2); }
  unsigned np = 0;
  f(0, 0, &np);
  fprintf(stderr, "[clGetPlatformIDs of %s called; MXCSR rounding bits now %u]\n", lib, (_mm_getcsr() >> 13) & 3);
}

int main(int argc, char** argv)
{
  nfArgc = argc;  nfArgv = argv;
  const char* in  = argc > 1 ? argv[1] : "in";
  const char* out = argc > 2 ? argv[2] : "out";
  if (getenv("OCL_LIB")) call_opencl(getenv("OCL_LIB"));
  pari_long SV[5] = {0, 0, 0, 0, 0};
  TgtMartinet((char*)in, (char*)out, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, SV);
  printf("StatVec: %lld %lld %lld %lld\n", (long long)SV[0], (long long)SV[1], (long long)SV[2], (long long)SV[3]);
  return 0;
}
