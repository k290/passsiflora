// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 Alperen Yavuz
// nf_gpu_musa.cpp -- back end for GPUs that BOINC does not know about (Moore Threads MUSA).
//
// The kernel and the runtime calls live in nf_musa_dev.cu, built with the vendor's compiler into the
// shared library libnf_musa.so (which links the vendor's runtime).  This program loads that library at
// run time from its own directory.
// BOINC cannot schedule such a GPU: a project server ignores coprocessor types it does not know (an app
// version declaring one counts as a CPU version, and requests for its work are skipped), so a <coproc> in
// cc_config.xml does not work with an unmodified project.  The program is therefore a CPU app version, as
// BOINC sees it, that also uses the card: without the library, the runtime, a device or a passed
// self-test it runs on the CPU, which is exactly what it declares to BOINC.  One task at a time uses each
// card (a lock file next to the program); the other tasks run on the CPU.  --device N picks the card
// (default 0); NF_GPU=0 switches the card off.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "boinc_api.h"
#include "app_ipc.h"
#include "nf_gpu_backend.h"

const bool nfBackendCpuApp = true;            // a CPU app version that may use a card

#ifdef __linux__
#include <dlfcn.h>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

static int (*p_count)();
static int (*p_format)();
static int (*p_init)(int, const int *, int, char *, int, int *, int *, unsigned long long *, int *);
static int (*p_set_nodes)(const double *, const long long *, size_t);
static int (*p_run)(unsigned, unsigned, const unsigned long long *, unsigned long long, long long *, unsigned *,
                    unsigned long long *);
static int (*p_selftest)(const unsigned *, unsigned, const long long *, unsigned, const long long *, unsigned, unsigned *,
                         unsigned);

// The directory of the program file itself (under BOINC the project directory, not the task's slot directory).
static std::string exeDir()
{
  char exe[4096];
  ssize_t n = readlink("/proc/self/exe", exe, sizeof(exe) - 1);
  if (n <= 0) return "";
  exe[n] = 0;
  std::string dir = exe;
  size_t p = dir.rfind('/');
  return p == std::string::npos ? std::string(".") : dir.substr(0, p);
}

static bool loadLibrary()
{
  std::string dir = exeDir();
  void *h = 0;
  if (!dir.empty()) h = dlopen((dir + "/libnf_musa.so").c_str(), RTLD_NOW);
  if (!h) h = dlopen("libnf_musa.so", RTLD_NOW);
  if (!h) { fprintf(stderr, "GPU mode: libnf_musa.so or the MUSA runtime could not be loaded; running on the CPU.\n"); return false; }
  p_count = (int (*)())dlsym(h, "nfm_device_count");
  // A library built for another prime-table format (an older libnf_musa.so, which has no nfm_table_format) must not run.
  p_format = (int (*)())dlsym(h, "nfm_table_format");
  if (!p_format || p_format() != NF_TABLE_FORMAT) {
    fprintf(stderr, "GPU mode: libnf_musa.so was built for another prime-table format (%d, this program %d); rebuild it\n"
                    "with build/build-musa.sh.  Running on the CPU.\n", p_format ? p_format() : 1, NF_TABLE_FORMAT);
    return false;
  }
  p_init = (int (*)(int, const int *, int, char *, int, int *, int *, unsigned long long *, int *))dlsym(h, "nfm_init");
  p_set_nodes = (int (*)(const double *, const long long *, size_t))dlsym(h, "nfm_set_nodes");
  p_run = (int (*)(unsigned, unsigned, const unsigned long long *, unsigned long long, long long *, unsigned *,
                   unsigned long long *))dlsym(h, "nfm_run");
  p_selftest = (int (*)(const unsigned *, unsigned, const long long *, unsigned, const long long *, unsigned, unsigned *,
                        unsigned))dlsym(h, "nfm_selftest");
  return p_count && p_init && p_set_nodes && p_run;
}

// One task per card: an exclusive lock on nf_musa_gpu<N>.lock next to the program, held until the process ends (the
// system releases it then, also when BOINC kills the task).  false: another task has the card.
static bool lockDevice(int dev)
{
  std::string dir = exeDir();
  std::string path = (dir.empty() ? std::string(".") : dir) + "/nf_musa_gpu" + std::to_string(dev) + ".lock";
  int fd = open(path.c_str(), O_RDWR | O_CREAT, 0666);
  if (fd < 0) {
    fprintf(stderr, "GPU mode: cannot open %s; using the card without a lock.\n", path.c_str());
    return true;
  }
  if (flock(fd, LOCK_EX | LOCK_NB) != 0) { close(fd); return false; }
  return true;                                 // fd stays open: the lock lasts as long as the process
}

// The start-up self-test of the filter on the device (nfm_selftest, with the table of nfm_init): false = the library has none,
// or the device does not reproduce the exact values, and the task runs on the CPU (it is a CPU task to BOINC).
static bool selfTest(const std::vector<int> &primes)
{
  if (!p_selftest) {
    fprintf(stderr, "GPU mode: libnf_musa.so has no self-test; rebuild it with build/build-musa.sh.  Running on the CPU.\n");
    return false;
  }
  NfSelfTest t;
  nfSelfTestBuild(primes, t);
  std::vector<unsigned> out(t.expect.size(), 0u);
  boinc_begin_critical_section();
  int rc = p_selftest(t.ts.data(), t.nt, t.xs.data(), t.nx, t.polys.data(), t.npoly, out.data(), (unsigned)out.size());
  boinc_end_critical_section();
  if (rc != 0) nfGpuFail(false, 1, "the MUSA self-test could not run (%d).", rc);
  if (!nfSelfTestVerify(t, out)) {
    fprintf(stderr, "GPU mode: the card is not used; running on the CPU.\n");
    return false;
  }
  return true;
}

bool nfBackendInit(int argc, char **argv, const std::vector<int> &primes)
{
  APP_INIT_DATA aid;
  boinc_get_init_data(aid);
  std::string gt = aid.gpu_type;
  if (!gt.empty() && gt != "CPU") {            // listed in app_info.xml as a GPU app version of another kind
    fprintf(stderr, "GPU mode: BOINC assigned a GPU of type \"%s\" to this CPU program; running on the CPU.\n", gt.c_str());
    return false;
  }
  if (!loadLibrary()) return false;

  int n = p_count();
  if (n <= 0) { fprintf(stderr, "GPU mode: no MUSA device found; running on the CPU.\n"); return false; }
  int dev = nfGpuDeviceArg(argc, argv);
  if (dev < 0) dev = 0;
  if (dev >= n) { fprintf(stderr, "GPU mode: no MUSA device %d (found %d); running on the CPU.\n", dev, n); return false; }
  if (!lockDevice(dev)) {
    fprintf(stderr, "GPU mode: MUSA device %d is in use by another task; running on the CPU.\n", dev);
    return false;
  }

  char name[256] = "?";
  int major = 0, minor = 0, mps = 0;
  unsigned long long mem = 0;
  if (p_init(dev, primes.data(), (int)primes.size(), name, sizeof(name), &major, &minor, &mem, &mps) != 0)
    nfGpuFail(false, 1, "the MUSA device %d could not be initialised.", dev);
  fprintf(stderr, "GPU mode: MUSA device %d \"%s\", %d multiprocessors, compute capability %d.%d, %llu MB, %u filter primes.\n",
          dev, name, mps, major, minor, (unsigned long long)(mem >> 20), (unsigned)(primes.size() / 8));
  return selfTest(primes);
}

void nfBackendShutdown() {}

void nfBackendSetNodes(const double *nd, const long long *ni, size_t n)
{
  if (p_set_nodes(nd, ni, n) != 0) nfGpuFail(false, 1, "MUSA: the nodes could not be copied to the device.");
}

long long nfBackendRun(size_t a, size_t b, const unsigned long long *cum, std::vector<NfSurv> &out, bool *overflow)
{
  unsigned long long total = cum[b - a];
  *overflow = false;
  std::vector<long long> raw((size_t)NF_OUT_CAP * 12);
  unsigned outcnt = 0;
  unsigned long long inspected = 0;
  boinc_begin_critical_section();
  int rc = p_run((unsigned)(b - a), (unsigned)a, cum, total, raw.data(), &outcnt, &inspected);
  boinc_end_critical_section();
  if (rc != 0) nfGpuFail(false, 1, "MUSA: the kernel failed (%d).", rc);
  if (outcnt > NF_OUT_CAP) { *overflow = true; return 0; }
  size_t base = out.size();
  out.resize(base + outcnt);
  for (unsigned k = 0; k < outcnt; k++) {
    NfSurv &s = out[base + k];
    s.node = raw[k * 12 + 0];  s.m2 = raw[k * 12 + 1];  s.m1 = raw[k * 12 + 2];
    for (int j = 0; j < 9; j++) s.b[j] = raw[k * 12 + 3 + j];
  }
  return (long long)inspected;
}

#else   // not Linux: no MUSA back end

bool nfBackendInit(int, char **, const std::vector<int> &) { return false; }
void nfBackendShutdown() {}
void nfBackendSetNodes(const double *, const long long *, size_t) {}
long long nfBackendRun(size_t, size_t, const unsigned long long *, std::vector<NfSurv> &, bool *) { return 0; }

#endif
