// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 Alperen Yavuz
// nf_gpu_opencl.cpp -- OpenCL back end (NVIDIA, AMD, Intel, ...).

#define CL_TARGET_OPENCL_VERSION 120
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <CL/cl.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#ifdef _WIN32
#  include <windows.h>
#  ifndef CREATE_WAITABLE_TIMER_HIGH_RESOLUTION
#    define CREATE_WAITABLE_TIMER_HIGH_RESOLUTION 0x00000002
#  endif
#else
#  include <time.h>
#endif

#include "boinc_api.h"
#include "boinc_opencl.h"
#include "app_ipc.h"
#include "nf_gpu_backend.h"

static const char *kernelSource =
#include "nf_enum_cl.h"
;

static cl_context g_ctx;
static cl_command_queue g_queue;
static cl_program g_prog;
static cl_kernel g_kern;
static cl_mem g_bNodeD, g_bNodeI, g_bCum, g_bPrimes, g_bOut, g_bOutCnt, g_bCnt;
static size_t g_cntCap = 0;
static unsigned g_np = 0;

const bool nfBackendCpuApp = false;           // a GPU app version: it uses the GPU BOINC assigned, or stops

// OpenCL errors stop the task.  CLCHK: the device not available counts as temporary (in use by another program); CLALLOC
// (memory allocations): running out of device or host memory also counts as temporary (BOINC: retry in a minute).
#define CLCHK(x) do { cl_int e_ = (x); if (e_ != CL_SUCCESS) \
    nfGpuFail(e_ == CL_DEVICE_NOT_AVAILABLE, 1, "OpenCL error %d at %s:%d (%s)", (int)e_, __FILE__, __LINE__, #x); } while (0)
#define CLALLOC(x) do { cl_int e_ = (x); if (e_ != CL_SUCCESS) \
    nfGpuFail(e_ == CL_DEVICE_NOT_AVAILABLE || e_ == CL_MEM_OBJECT_ALLOCATION_FAILURE || e_ == CL_OUT_OF_RESOURCES || \
              e_ == CL_OUT_OF_HOST_MEMORY, 1, "OpenCL error %d at %s:%d (%s): could not allocate memory", (int)e_, __FILE__, \
              __LINE__, #x); } while (0)

// One millisecond of sleep.  Windows' Sleep(1) can last a whole 15.6 ms timer tick; a high-resolution waitable timer
// (Windows 10 1803 and newer) sleeps about 1 ms.
static void sleep1ms()
{
#ifdef _WIN32
  static HANDLE timer = CreateWaitableTimerExW(NULL, NULL, CREATE_WAITABLE_TIMER_HIGH_RESOLUTION, TIMER_ALL_ACCESS);
  LARGE_INTEGER due;
  due.QuadPart = -10000;                     // 1 ms, relative, in 100 ns units
  if (timer && SetWaitableTimer(timer, &due, 0, NULL, NULL, FALSE)) WaitForSingleObject(timer, INFINITE);
  else Sleep(1);
#else
  struct timespec ts = { 0, 1000000 };
  nanosleep(&ts, NULL);
#endif
}

// Waits for a kernel by polling its event every millisecond.  NVIDIA's OpenCL spins a CPU core in its blocking calls
// (clFinish, a blocking read) for as long as the kernel runs; a GPU task then took a whole core.  Launches last about 100 ms.
static void waitKernel(cl_event ev)
{
  CLCHK(clFlush(g_queue));
  for (;;) {
    cl_int s = CL_QUEUED;
    CLCHK(clGetEventInfo(ev, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(s), &s, NULL));
    if (s == CL_COMPLETE) break;
    if (s < 0) nfGpuFail(false, 1, "OpenCL error %d: the kernel failed", (int)s);
    sleep1ms();
  }
  clReleaseEvent(ev);
}

// The start-up self-test of the filter on this device (nf_selftest_t2 with the table of g_bPrimes).  A device that does not
// reproduce the exact values stops the task (exit status 6).
static void selfTest(const std::vector<int> &primes)
{
  cl_int st;
  cl_kernel kern = clCreateKernel(g_prog, NF_SELFTEST_NAME, &st);
  if (st != CL_SUCCESS) nfGpuFail(false, 1, "the kernel has no self-test (%d).", (int)st);
  NfSelfTest t;
  nfSelfTestBuild(primes, t);
  cl_uint nt = t.nt, nx = t.nx, npoly = t.npoly, np = t.np, total = (cl_uint)t.expect.size();
  cl_mem bTs = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, t.ts.size() * sizeof(cl_uint), (void *)t.ts.data(), &st);  CLALLOC(st);
  cl_mem bXs = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, t.xs.size() * sizeof(cl_long), (void *)t.xs.data(), &st);  CLALLOC(st);
  cl_mem bPol = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, t.polys.size() * sizeof(cl_long), (void *)t.polys.data(), &st);  CLALLOC(st);
  cl_mem bOut = clCreateBuffer(g_ctx, CL_MEM_WRITE_ONLY, (size_t)total * sizeof(cl_uint), NULL, &st);  CLALLOC(st);
  int ai = 0;
  CLCHK(clSetKernelArg(kern, ai++, sizeof(cl_mem), &bTs));
  CLCHK(clSetKernelArg(kern, ai++, sizeof(cl_uint), &nt));
  CLCHK(clSetKernelArg(kern, ai++, sizeof(cl_mem), &bXs));
  CLCHK(clSetKernelArg(kern, ai++, sizeof(cl_uint), &nx));
  CLCHK(clSetKernelArg(kern, ai++, sizeof(cl_mem), &bPol));
  CLCHK(clSetKernelArg(kern, ai++, sizeof(cl_uint), &npoly));
  CLCHK(clSetKernelArg(kern, ai++, sizeof(cl_mem), &g_bPrimes));
  CLCHK(clSetKernelArg(kern, ai++, sizeof(cl_uint), &np));
  CLCHK(clSetKernelArg(kern, ai++, sizeof(cl_mem), &bOut));
  CLCHK(clSetKernelArg(kern, ai++, sizeof(cl_uint), &total));
  size_t local = 64, global = ((size_t)total + local - 1) / local * local;
  std::vector<unsigned> out(total, 0u);
  boinc_begin_critical_section();
  CLCHK(clEnqueueNDRangeKernel(g_queue, kern, 1, NULL, &global, &local, 0, NULL, NULL));
  CLCHK(clEnqueueReadBuffer(g_queue, bOut, CL_TRUE, 0, (size_t)total * sizeof(cl_uint), out.data(), 0, NULL, NULL));
  boinc_end_critical_section();
  clReleaseMemObject(bTs);  clReleaseMemObject(bXs);  clReleaseMemObject(bPol);  clReleaseMemObject(bOut);
  // kern is kept, like g_kern (build-windows-gpu.sh imports no clReleaseKernel)
  if (!nfSelfTestVerify(t, out)) nfGpuFail(false, 6, "this GPU or driver does not reproduce the filter's exact values.");
}

// Run by hand without an init_data.xml: the OpenCL GPUs of all platforms, numbered in order; --device N picks one (default 0).
static void pickByHand(int idx, cl_platform_id *plat, cl_device_id *dev)
{
  cl_platform_id plats[16];  cl_uint np = 0;
  if (clGetPlatformIDs(16, plats, &np) != CL_SUCCESS || np == 0)
    nfGpuFail(false, 1, "no OpenCL platform found (is a graphics driver with OpenCL installed?).");
  std::vector<cl_platform_id> gp;
  std::vector<cl_device_id> gd;
  for (cl_uint p = 0; p < np && p < 16; p++) {
    cl_device_id devs[16];  cl_uint nd = 0;
    if (clGetDeviceIDs(plats[p], CL_DEVICE_TYPE_GPU, 16, devs, &nd) != CL_SUCCESS) continue;
    char pname[256] = "?";
    clGetPlatformInfo(plats[p], CL_PLATFORM_NAME, sizeof(pname), pname, NULL);
    for (cl_uint d = 0; d < nd && d < 16; d++) {
      char name[256] = "?";
      clGetDeviceInfo(devs[d], CL_DEVICE_NAME, sizeof(name), name, NULL);
      fprintf(stderr, "OpenCL GPU %u: \"%s\" (%s)\n", (unsigned)gd.size(), name, pname);
      gp.push_back(plats[p]);  gd.push_back(devs[d]);
    }
  }
  if (gd.empty()) nfGpuFail(false, 1, "no OpenCL GPU found.");
  if (idx < 0 || idx >= (int)gd.size()) nfGpuFail(false, 1, "no OpenCL GPU %d (found %u).", idx, (unsigned)gd.size());
  *plat = gp[idx];  *dev = gd[idx];
}

bool nfBackendInit(int argc, char **argv, const std::vector<int> &primes)
{
  NfGpuChoice c = nfGpuChoose(argc, argv);
  g_np = (unsigned)(primes.size() / 8);

  cl_platform_id plat = NULL;
  cl_device_id dev = NULL;
  const char *how;
  if (!c.type.empty()) {                       // the GPU BOINC assigned (or init_data.xml names): BOINC finds it
    int ptype = c.type == "NVIDIA" ? PROC_TYPE_NVIDIA_GPU : c.type == "ATI" ? PROC_TYPE_AMD_GPU :
                c.type == "intel_gpu" ? PROC_TYPE_INTEL_GPU : -1;
    if (ptype < 0) nfGpuFail(false, 1, "BOINC assigned a GPU of type \"%s\", which this OpenCL program does not use.", c.type.c_str());
    int rv = boinc_get_opencl_ids(argc, argv, ptype, &dev, &plat);
    if (rv) nfGpuFail(false, 1, "BOINC could not find the OpenCL device it assigned (%s device %d, error %d).", c.type.c_str(),
                      c.device, rv);
    how = c.byHand ? "named by init_data.xml" : "assigned by BOINC";
  } else {
    pickByHand(c.device, &plat, &dev);
    how = "run by hand; --device N picks another";
  }

  char name[256] = "?";
  clGetDeviceInfo(dev, CL_DEVICE_NAME, sizeof(name), name, NULL);
  cl_uint cus = 0;  clGetDeviceInfo(dev, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cus), &cus, NULL);
  cl_ulong gmem = 0;  clGetDeviceInfo(dev, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(gmem), &gmem, NULL);
  char ext[8192] = "";  clGetDeviceInfo(dev, CL_DEVICE_EXTENSIONS, sizeof(ext) - 1, ext, NULL);
  fprintf(stderr, "GPU mode: OpenCL device %d \"%s\" (%s), %u compute units, %llu MB, %u filter primes.\n",
          c.device, name, how, cus, (unsigned long long)(gmem >> 20), g_np);
  if (!strstr(ext, "cl_khr_fp64")) nfGpuFail(false, 1, "the device has no double precision (cl_khr_fp64).");

  cl_int st;
  cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)plat, 0 };
  g_ctx = clCreateContext(cps, 1, &dev, NULL, NULL, &st);  CLALLOC(st);
  g_queue = clCreateCommandQueue(g_ctx, dev, 0, &st);  CLALLOC(st);

  const char *src = kernelSource;
  size_t len = strlen(src);
  g_prog = clCreateProgramWithSource(g_ctx, 1, &src, &len, &st);  CLCHK(st);
  st = clBuildProgram(g_prog, 1, &dev, "-cl-std=CL1.2", NULL, NULL);
  if (st != CL_SUCCESS) {
    size_t n = 0;
    clGetProgramBuildInfo(g_prog, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &n);
    std::vector<char> log(n + 1, 0);
    clGetProgramBuildInfo(g_prog, dev, CL_PROGRAM_BUILD_LOG, n, log.data(), NULL);
    fprintf(stderr, "OpenCL build log:\n%s\n", log.data());
    nfGpuFail(false, 1, "the driver could not compile the kernel (OpenCL error %d).", (int)st);
  }
  g_kern = clCreateKernel(g_prog, NF_KERNEL_NAME, &st);  CLCHK(st);

  g_bNodeD = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY, NF_MAX_NODES * NF_ND * sizeof(double), NULL, &st);  CLALLOC(st);
  g_bNodeI = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY, NF_MAX_NODES * NF_NI * sizeof(long long), NULL, &st);  CLALLOC(st);
  g_bCum = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY, (NF_MAX_NODES + 1) * sizeof(unsigned long long), NULL, &st);  CLALLOC(st);
  g_bPrimes = clCreateBuffer(g_ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, primes.size() * sizeof(int),
                             (void *)primes.data(), &st);  CLALLOC(st);
  g_bOut = clCreateBuffer(g_ctx, CL_MEM_WRITE_ONLY, (size_t)NF_OUT_CAP * 12 * sizeof(long long), NULL, &st);  CLALLOC(st);
  g_bOutCnt = clCreateBuffer(g_ctx, CL_MEM_READ_WRITE, sizeof(unsigned), NULL, &st);  CLALLOC(st);
  g_bCnt = NULL;  g_cntCap = 0;
  selfTest(primes);
  return true;
}

void nfBackendShutdown()
{
  clFinish(g_queue);
}

void nfBackendSetNodes(const double *nd, const long long *ni, size_t n)
{
  CLCHK(clEnqueueWriteBuffer(g_queue, g_bNodeD, CL_TRUE, 0, n * NF_ND * sizeof(double), nd, 0, NULL, NULL));
  CLCHK(clEnqueueWriteBuffer(g_queue, g_bNodeI, CL_TRUE, 0, n * NF_NI * sizeof(long long), ni, 0, NULL, NULL));
}

long long nfBackendRun(size_t a, size_t b, const unsigned long long *cum, std::vector<NfSurv> &out, bool *overflow)
{
  unsigned long long total = cum[b - a];
  *overflow = false;

  if (total > g_cntCap) {
    if (g_bCnt) clReleaseMemObject(g_bCnt);
    cl_int st;
    g_cntCap = (size_t)(total + (total >> 2) + 1024);
    g_bCnt = clCreateBuffer(g_ctx, CL_MEM_READ_WRITE, g_cntCap * sizeof(unsigned), NULL, &st);  CLALLOC(st);
  }

  unsigned zero = 0;
  CLCHK(clEnqueueWriteBuffer(g_queue, g_bCum, CL_TRUE, 0, (b - a + 1) * sizeof(unsigned long long), cum, 0, NULL, NULL));
  CLCHK(clEnqueueWriteBuffer(g_queue, g_bOutCnt, CL_TRUE, 0, sizeof(unsigned), &zero, 0, NULL, NULL));

  cl_uint nnodes = (cl_uint)(b - a), nodebase = (cl_uint)a, np = g_np, cap = NF_OUT_CAP;
  cl_ulong tot = total;
  int ai = 0;
  CLCHK(clSetKernelArg(g_kern, ai++, sizeof(cl_mem), &g_bNodeD));
  CLCHK(clSetKernelArg(g_kern, ai++, sizeof(cl_mem), &g_bNodeI));
  CLCHK(clSetKernelArg(g_kern, ai++, sizeof(cl_mem), &g_bCum));
  CLCHK(clSetKernelArg(g_kern, ai++, sizeof(cl_uint), &nnodes));
  CLCHK(clSetKernelArg(g_kern, ai++, sizeof(cl_uint), &nodebase));
  CLCHK(clSetKernelArg(g_kern, ai++, sizeof(cl_ulong), &tot));
  CLCHK(clSetKernelArg(g_kern, ai++, sizeof(cl_mem), &g_bPrimes));
  CLCHK(clSetKernelArg(g_kern, ai++, sizeof(cl_uint), &np));
  CLCHK(clSetKernelArg(g_kern, ai++, sizeof(cl_mem), &g_bOut));
  CLCHK(clSetKernelArg(g_kern, ai++, sizeof(cl_mem), &g_bOutCnt));
  CLCHK(clSetKernelArg(g_kern, ai++, sizeof(cl_uint), &cap));
  CLCHK(clSetKernelArg(g_kern, ai++, sizeof(cl_mem), &g_bCnt));

  size_t local = 64;
  size_t global = (size_t)((total + local - 1) / local * local);
  boinc_begin_critical_section();
  cl_event ev;
  CLCHK(clEnqueueNDRangeKernel(g_queue, g_kern, 1, NULL, &global, &local, 0, NULL, &ev));
  waitKernel(ev);

  unsigned outcnt = 0;
  CLCHK(clEnqueueReadBuffer(g_queue, g_bOutCnt, CL_TRUE, 0, sizeof(unsigned), &outcnt, 0, NULL, NULL));
  boinc_end_critical_section();

  if (outcnt > NF_OUT_CAP) { *overflow = true; return 0; }

  size_t base = out.size();
  out.resize(base + outcnt);
  if (outcnt) {
    std::vector<long long> raw((size_t)outcnt * 12);
    CLCHK(clEnqueueReadBuffer(g_queue, g_bOut, CL_TRUE, 0, raw.size() * sizeof(long long), raw.data(), 0, NULL, NULL));
    for (unsigned k = 0; k < outcnt; k++) {
      NfSurv &s = out[base + k];
      s.node = raw[k * 12 + 0];  s.m2 = raw[k * 12 + 1];  s.m1 = raw[k * 12 + 2];
      for (int j = 0; j < 9; j++) s.b[j] = raw[k * 12 + 3 + j];
    }
  }

  std::vector<unsigned> cnt((size_t)total);
  CLCHK(clEnqueueReadBuffer(g_queue, g_bCnt, CL_TRUE, 0, (size_t)total * sizeof(unsigned), cnt.data(), 0, NULL, NULL));
  long long sum = 0;
  for (size_t k = 0; k < (size_t)total; k++) sum += cnt[k];
  return sum;
}
