// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 Alperen Yavuz
// nf_gpu.h -- OpenCL back end for the innermost enumeration loops.
#ifndef NF_GPU_H
#define NF_GPU_H

#include <vector>

#define NF_ND 18
#define NF_NI 35

// Everything the m2/m1 loops need from the outer loops (layout shared with the kernel, nf_kernel.h).
struct NfNode {
  double d[NF_ND];
  long long i[NF_NI];
};

// Pick the GPU (the one BOINC assigned; run by hand: init_data.xml, --device N or the first
// GPU) and build the kernel.  A GPU program (OpenCL, CUDA) returns true or stops the task
// when it cannot use the GPU (nf_gpu_backend.h: nfGpuFail).  It returns false, and the task
// runs on the CPU, only with NF_GPU=0 (a developer switch) or in the MUSA program, which
// BOINC runs as a CPU app version.
bool nfGpuInit(int argc, char **argv, int numP, const int *pSet);
bool nfGpuActive();

// Pending-node queue.
bool nfGpuWantFlush();
bool nfGpuEmpty();
void nfGpuAdd(const NfNode &node);

// Runs all pending nodes, clears the queue and returns the number of
// polynomials of the original enumeration that were inspected.  The survivors of
// the GPU pre-filter are appended to polys (11 values each, in enumeration
// order: 1, bb1..bb10).
long long nfGpuRun(std::vector<long long> &polys);

// bb[0..10] (bb[0] = 1): does the GPU's filter keep this polynomial?  The kernel's own filter compiled for the CPU, with the
// table the device runs with; the checked window at the start of a task compares the device's survivors with it.
bool nfGpuRefSurvives(const long long *bb);

void nfGpuShutdown();

#endif
