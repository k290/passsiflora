// app_feed: feeds a file of polynomials to the REAL polynomial pipeline of GetDecics with the Passiflora patch (CPU version).
// The pipeline is testPolys() of the patched TgtMartinet.cpp: pre-filter (pdfilter.f90) -> exact polynomial-discriminant test
// (polDiscTest_cpu) -> irreducibility test -> field-discriminant test -> result file.  Only the enumeration of the search space
// is replaced by a file reader.  The order of TgtMartinet() is kept: pari_init, read S, then (as nfGpuInit does in the
// application) the optional OpenCL-library call, then testPolys(); its first call runs pdf_init.
//   usage : app_feed polys.txt [result file]       polys.txt: one polynomial per line, 11 integers, leading coefficient first
//   output: one line "StatVec: a b c d" on stdout (polynomials fed, passing the polynomial-discriminant test, irreducible,
//           passing the field-discriminant test); the result lines go to the result file (default "out")
//   env   : OCL_LIB=<OpenCL library>  call its clGetPlatformIDs before the first testPolys, as nfGpuInit does in the application
//           FEED_ALONE=1               one testPolys call (a batch of 1) per polynomial instead of one call for the whole file
//           NF_VERIFY=1, NF_PDFDEBUG=1, NF_FORCE_BASE=1, NF_GENERAL_ALL=1: read by the application and the filter themselves
// Linked with the application's own TgtMartinet.o, polDiscTest_cpu.o and nf_gpu_stub.o and the filter objects (see repro.sh).
#include <dlfcn.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <xmmintrin.h>
#include "GetDecics.h"
#include "TgtMartinet.h"

using std::ofstream;

extern GEN Sg;
extern char *polGoodFlag;
extern int nfArgc;
extern char **nfArgv;
int numBlocks = 10240;
int threadsPerBlock = 1;
int polyBufferSize = 10240;                       // the CPU batch size of the application
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

int main(int argc, char **argv)
{
  const char *file = argc > 1 ? argv[1] : "polys.txt";
  const char *outname = argc > 2 ? argv[2] : "out";
  nfArgc = argc;  nfArgv = argv;
  std::ifstream in(file);
  if (!in) { fprintf(stderr, "cannot read %s\n", file); return 2; }
  std::vector<pari_long> buf;
  long long v;
  while (in >> v) buf.push_back((pari_long)v);
  if (buf.size() % 11) { fprintf(stderr, "%s: the number of integers must be a multiple of 11\n", file); return 2; }
  int n = (int)(buf.size() / 11);
  for (int i = 0; i < n; i++) if (buf[i * 11] != 1) { fprintf(stderr, "polynomial %d is not monic\n", i + 1); return 2; }

  pari_init(300000000, 10000);                    // as TgtMartinet()
  Sg = gp_read_str((char *)"[2, 5]");              // the set S of the work units
  polGoodFlag = (char *)malloc(polyBufferSize > n ? polyBufferSize : n);
  if (getenv("OCL_LIB")) call_opencl(getenv("OCL_LIB"));   // nfGpuInit runs before the first pdf_init in the application

  ofstream out(outname);
  pari_long SV[5] = {0, 0, 0, 0, 0};
  pari_long loopIdx[4] = {0, 0, 0, 0};
  if (getenv("FEED_ALONE")) {
    for (int i = 0; i < n; i++) testPolys((long long *)&buf[(size_t)i * 11], 1, SV, loopIdx, out);
  } else {
    testPolys((long long *)buf.data(), n, SV, loopIdx, out);
  }
  printf("StatVec: %lld %lld %lld %lld\n", (long long)SV[0], (long long)SV[1], (long long)SV[2], (long long)SV[3]);
  return 0;
}
