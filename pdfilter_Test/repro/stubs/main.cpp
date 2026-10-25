// Replaces GetDecics.cpp (the BOINC wrapper) for builds without the BOINC libraries: the same call into TgtMartinet(),
// which reads the work unit "in" and writes the result file "out", and the same final "StatVec: a b c d" line.
#include <cstdio>
#include "GetDecics.h"
#include "TgtMartinet.h"

int numBlocks = 10240;
int threadsPerBlock = 1;
int polyBufferSize = 10240;
extern int nfArgc;
extern char **nfArgv;

int do_checkpoint(pari_long*, pari_long*) { return 0; }

int main(int argc, char** argv)
{
  nfArgc = argc;  nfArgv = argv;
  const char* in  = argc > 1 ? argv[1] : "in";
  const char* out = argc > 2 ? argv[2] : "out";
  pari_long SV[5] = {0, 0, 0, 0, 0};
  TgtMartinet((char*)in, (char*)out, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, SV);
  printf("StatVec: %lld %lld %lld %lld\n", (long long)SV[0], (long long)SV[1], (long long)SV[2], (long long)SV[3]);
  return 0;
}
