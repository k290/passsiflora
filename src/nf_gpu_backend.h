// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 Alperen Yavuz
// nf_gpu_backend.h -- what a GPU back end (OpenCL, CUDA, ...) has to provide.
#ifndef NF_GPU_BACKEND_H
#define NF_GPU_BACKEND_H

#include <vector>

#include "nf_gpu.h"

struct NfSurv { long long node, m2, m1, b[9]; };

// The prime-table format this program builds, and the kernel entry point for it (nf_kernel.h).  Device code built for
// another format has no entry point of this name (CUDA, OpenCL) or reports another format (MUSA): the back end refuses it.
#define NF_TABLE_FORMAT 2
#define NF_KERNEL_NAME "nf_enum_t2"
#define NF_SELFTEST_NAME "nf_selftest_t2"

static const size_t NF_MAX_NODES = 1u << 17;   // nodes per batch
static const unsigned NF_OUT_CAP = 1u << 18;   // survivor records per launch

// Pick the device and build / load the kernel.  primes = the filter-prime table, 8 ints per prime
// (q, 2^16, 2^32, 2^48, 2^64 mod q, (q-1)/2, bits of -1/q mod 2^32, -2^16 mod q).  Returns false when no device of
// this kind is assigned to the task (the app then runs on the CPU); exits on real errors.
bool nfBackendInit(int argc, char **argv, const std::vector<int> &primes);

// Upload the nodes of the batch (n nodes, NF_ND doubles and NF_NI integers each).
void nfBackendSetNodes(const double *nd, const long long *ni, size_t n);

// Run nodes [a, b) of the uploaded batch.  cum[0..b-a] = prefix sums of the nodes' m2 counts.
// Appends the survivors to out and returns the number of polynomials inspected; if the survivor
// buffer overflowed it sets *overflow and appends nothing.
long long nfBackendRun(size_t a, size_t b, const unsigned long long *cum,
                       std::vector<NfSurv> &out, bool *overflow);

void nfBackendShutdown();

// Start-up self-test of the filter on the device, run by every back end at the end of nfBackendInit: the kernel
// nf_selftest_t2 (layout of its output: nf_kernel.h) evaluates redc, coef_mod, chi_reject and survives on the inputs below,
// and the host compares the output with values computed independently of the kernel's code (exact arithmetic modulo q;
// the symbol of the discriminant from the resultant by Euclid with inverses).
struct NfSelfTest {
  unsigned np, nt, nx, npoly;
  std::vector<unsigned> ts;        // np * nt inputs of redc, grouped by prime
  std::vector<long long> xs;       // nx inputs of coef_mod (the same for every prime)
  std::vector<long long> polys;    // npoly * 10 coefficients c1..c10 (leading coefficient 1)
  std::vector<unsigned> expect;    // the kernel's output, work item by work item
};
void nfSelfTestBuild(const std::vector<int> &primes, NfSelfTest &t);
// true if out equals t.expect; otherwise prints the first differences and returns false: the back end then refuses GPU mode
// and the application runs on the CPU.  NF_SELFTEST_FAULT=1 changes one value of out first (test hook).
bool nfSelfTestVerify(const NfSelfTest &t, std::vector<unsigned> &out);

#endif
