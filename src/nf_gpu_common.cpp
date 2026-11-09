// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 Alperen Yavuz
// nf_gpu_common.cpp -- everything about the GPU mode that does not depend on the API used:
// the filter-prime table, the node queue, the batch sizing, the sorting of the survivors.

#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

// NF_GPU_NO_BOINC: the CPU-side tests (pdfilter_Test/test_gpu_filter.cpp) include this file without BOINC's libraries.
#ifndef NF_GPU_NO_BOINC
#include "boinc_api.h"
#include "app_ipc.h"
#endif
#include "nf_gpu.h"
#include "nf_gpu_backend.h"

static const size_t MAX_PRIMES = 24;
static const size_t MIN_PRIMES = 8;
static const unsigned long long CHUNK_ITEMS = 1ull << 20;

static bool g_active = false;
static std::vector<NfNode> g_nodes;
static unsigned long long g_items = 0;
static unsigned long long g_target = 1ull << 18;

// ---- filter primes -----------------------------------------------------------

static long long powmod(long long b, long long e, long long m)
{
  long long r = 1;  b %= m;  if (b < 0) b += m;
  while (e > 0) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
  return r;
}

static bool isPrime(long long n)
{
  if (n < 2) return false;
  for (long long d = 2; d * d <= n; d++) if (n % d == 0) return false;
  return true;
}

// The kernel's row for an odd prime q < 8192: q; 2^16, 2^32, 2^48, 2^64 mod q (the loading weights; 2^16 mod q is also R mod q,
// the Montgomery one); (q-1)/2; qneg = -1/q mod 2^32 (the bits of an unsigned 32-bit value); negr = -R mod q
// (docs/whitepaper.tex, Theorem E).
static void primeRow(long long q, int row[8])
{
  unsigned inv = (unsigned)q;                          // Newton: q*inv = 1 mod 2^k doubles k (q*q = 1 mod 8)
  for (int k = 0; k < 4; k++) inv *= 2u - (unsigned)q * inv;
  unsigned qneg = 0u - inv;
  int qnegbits;  memcpy(&qnegbits, &qneg, 4);
  long long r = powmod(2, 16, q);
  int v[8] = { (int)q, (int)r, (int)powmod(2, 32, q), (int)powmod(2, 48, q),
               (int)powmod(2, 64, q), (int)((q - 1) / 2), qnegbits, (int)(q - r) };
  memcpy(row, v, sizeof v);
}

// q = 1 (mod 4) prime in (256, 8192), every p in S a nonzero square mod q; smallest first (Euler: one step per bit of (q-1)/2).
static std::vector<int> buildPrimeTable(int numP, const int *pSet)
{
  std::vector<int> t;
  for (long long q = 257; q < 8192 && t.size() / 8 < MAX_PRIMES; q++) {
    if (q % 4 != 1 || !isPrime(q)) continue;
    bool ok = true;
    for (int k = 0; k < numP; k++) {
      long long p = pSet[k];
      if (p % q == 0 || powmod(p, (q - 1) / 2, q) != 1) { ok = false; break; }
    }
    if (!ok) continue;
    int row[8];
    primeRow(q, row);
    t.insert(t.end(), row, row + 8);
  }
  return t;
}

// ---- the kernel's filter on the host, and the start-up self-test ---------------------------
// The kernel's source compiled for the CPU (nf_port_cpu.h): the checked window at the start of a task compares the device's
// survivors with its survives(); the self-test compares the device with arithmetic that does not use the kernel's code.

#include <climits>
#include "nf_port_cpu.h"
#include "nf_kernel.h"

static std::vector<int> g_table;            // the table the device runs with

bool nfGpuRefSurvives(const long long *bb)
{
  i64 c[11];
  for (int k = 0; k <= 10; k++) c[k] = bb[k];
  return survives(c, g_table.data(), (nfuint)(g_table.size() / 8)) != 0;
}

// chi(Res(f, f')) for f = x^10 + c[1] x^9 + ... + c[10] over F_q: 1, -1, or 0 when the resultant vanishes; by Euclid with
// inverses (the recursion Res(A, B) = (-1)^(ab) lc(B)^(a - deg R) Res(B, R), R = A mod B), not by the kernel's chain.  As
// q = 1 (mod 4), it is chi(disc f), and the kernel must reject at q exactly when it is -1 (docs/whitepaper.tex, Theorem E).
static int exactChi(const long long *c, long long q)
{
  long long A[11], B[11], R[11];
  for (int k = 0; k <= 10; k++) A[k] = k == 10 ? 1 : ((c[10 - k] % q) + q) % q;
  for (int k = 1; k <= 10; k++) B[k - 1] = k * A[k] % q;
  int a = 10, b = 9;
  while (b >= 0 && B[b] == 0) b--;
  long long r = 1;
  for (;;) {
    if (b < 0) return 0;
    if (b == 0) { r = r * powmod(B[0], a, q) % q; break; }
    for (int k = 0; k <= a; k++) R[k] = A[k];
    long long inv = powmod(B[b], q - 2, q);
    for (int k = a; k >= b; k--) {
      long long f = R[k] * inv % q;
      if (f) for (int j = 0; j <= b; j++) R[k - b + j] = ((R[k - b + j] - f * B[j]) % q + q) % q;
    }
    int rd = b - 1;
    while (rd >= 0 && R[rd] == 0) rd--;
    if (rd < 0) return 0;
    if ((a * b) & 1) r = (q - r) % q;
    r = r * powmod(B[b], a - rd, q) % q;
    for (int k = 0; k <= b; k++) A[k] = B[k];
    for (int k = 0; k <= rd; k++) B[k] = R[k];
    a = b;  b = rd;
  }
  if (r == 0) return 0;
  return powmod(r, (q - 1) / 2, q) == 1 ? 1 : -1;
}

namespace {
struct NfRng {                              // splitmix64 with a fixed seed: the same test on every device
  unsigned long long s;
  unsigned long long next() {
    unsigned long long z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }
  long long range(long long b) { return (long long)(next() % (unsigned long long)(2 * b + 1)) - b; }
};
}

void nfSelfTestBuild(const std::vector<int> &primes, NfSelfTest &t)
{
  NfRng rng{20261030ull};
  t.np = (unsigned)(primes.size() / 8);
  t.nt = 32;  t.nx = 48;  t.npoly = 240;

  // redc: the extremes of [0, 8q^2) and of the 16-bit halves, then pseudo-random inputs
  t.ts.clear();
  for (unsigned p = 0; p < t.np; p++) {
    unsigned long long q = (unsigned)primes[p * 8], lim = 8 * q * q;
    const unsigned long long fixed[] = { 0, 1, 2, q - 1, q, q + 1, 2 * q - 1, 2 * q, (q - 1) * (q - 1), 4 * q * q, 65535, 65536,
                                         lim - 2, lim - 1 };
    for (unsigned long long v : fixed) if (v < lim) t.ts.push_back((unsigned)v);
    while (t.ts.size() < (size_t)(p + 1) * t.nt) t.ts.push_back((unsigned)(rng.next() % lim));
  }

  // coef_mod: the extremes of int64 and of its 16-bit limbs, then pseudo-random inputs (any int64, and small ones)
  t.xs = { 0, 1, -1, 2, -2, 65535, 65536, -65535, -65536, 4294967295LL, 4294967296LL, -4294967296LL, -4294967297LL,
           281474976710655LL, 281474976710656LL, -281474976710656LL, -281474976710657LL, LLONG_MAX, LLONG_MIN, LLONG_MAX - 1,
           LLONG_MIN + 1, 1000000, -1000000 };
  while (t.xs.size() < t.nx) t.xs.push_back(t.xs.size() & 1 ? (long long)rng.next() : rng.range(1000000));

  // polynomials, six families in turn: random coefficients up to 10^6; any int64; depressed with |c| <= 3 (abnormal chains);
  // x^10 + b x^7 + c x^4 + d x; coefficients that are multiples of a table prime; (x - r)^2 g + q h, so that disc = 0 mod q
  t.polys.clear();
  for (unsigned j = 0; j < t.npoly; j++) {
    long long c[11] = { 1 };
    long long qs = primes[8 * (rng.next() % t.np)];
    switch (j % 6) {
      case 0: for (int k = 1; k <= 10; k++) c[k] = rng.range(1000000); break;
      case 1: for (int k = 1; k <= 10; k++) c[k] = (long long)rng.next(); break;
      case 2: c[1] = 0; for (int k = 2; k <= 10; k++) c[k] = rng.range(3); break;
      case 3: for (int k = 1; k <= 10; k++) c[k] = (k == 3 || k == 6 || k == 9) ? rng.range(1000000) : 0; break;
      case 4: for (int k = 1; k <= 10; k++) c[k] = (rng.next() & 1) ? qs * rng.range(1000000) : rng.range(1000000); break;
      default: {
        long long r = rng.range(1000), g[9], f[11] = { 0 };
        g[0] = 1;  for (int k = 1; k <= 8; k++) g[k] = rng.range(1000);          // g monic of degree 8, leading coefficient first
        const long long s[3] = { 1, -2 * r, r * r };
        for (int u = 0; u < 3; u++) for (int v = 0; v <= 8; v++) f[u + v] += s[u] * g[v];
        for (int k = 1; k <= 10; k++) c[k] = f[k] + qs * rng.range(1000);
      }
    }
    t.polys.insert(t.polys.end(), c + 1, c + 11);
  }

  // the expected output
  t.expect.clear();
  for (unsigned p = 0; p < t.np; p++) {            // redc(t) = (t + m q) / 2^16 with m = -t / q mod 2^16, in 64 bits
    unsigned long long q = (unsigned)primes[p * 8], inv = q;
    for (int k = 0; k < 5; k++) inv *= 2 - q * inv; // Newton: q * inv = 1 mod 2^k, k doubling from 3
    for (unsigned i = 0; i < t.nt; i++) {
      unsigned long long tv = t.ts[(size_t)p * t.nt + i], m = (0x10000 - (tv * inv & 0xFFFF)) & 0xFFFF;
      t.expect.push_back((unsigned)((tv + m * q) >> 16));
    }
  }
  for (unsigned p = 0; p < t.np; p++) {            // coef_mod(x) = x mod q in [0, q)
    long long q = primes[p * 8];
    for (unsigned i = 0; i < t.nx; i++) t.expect.push_back((unsigned)(((t.xs[i] % q) + q) % q));
  }
  std::vector<char> keep(t.npoly, 1);
  for (unsigned p = 0; p < t.np; p++)              // chi_reject: 1 exactly when chi(disc f) = -1
    for (unsigned j = 0; j < t.npoly; j++) {
      long long c[11];
      c[0] = 1;
      for (int k = 1; k <= 10; k++) c[k] = t.polys[(size_t)j * 10 + k - 1];
      bool rej = exactChi(c, primes[p * 8]) == -1;
      t.expect.push_back(rej ? 1u : 0u);
      if (rej) keep[j] = 0;
    }
  for (unsigned j = 0; j < t.npoly; j++) t.expect.push_back(keep[j] ? 1u : 0u);   // survives: no prime rejects
}

bool nfSelfTestVerify(const NfSelfTest &t, std::vector<unsigned> &out)
{
  const char *fault = getenv("NF_SELFTEST_FAULT");
  if (fault && *fault == '1' && !out.empty()) out[out.size() / 2] ^= 1u;
  static const char *const part[4] = { "redc", "coef_mod", "chi_reject", "survives" };
  size_t n1 = (size_t)t.np * t.nt, n2 = (size_t)t.np * t.nx, n3 = (size_t)t.np * t.npoly, bad = 0;
  for (size_t g = 0; g < t.expect.size(); g++) {
    unsigned v = g < out.size() ? out[g] : 0u;
    if (v == t.expect[g]) continue;
    if (bad++ < 5)
      fprintf(stderr, "GPU self-test: %s, work item %u: device %u, exact %u\n",
              part[g < n1 ? 0 : g < n1 + n2 ? 1 : g < n1 + n2 + n3 ? 2 : 3], (unsigned)g, v, t.expect[g]);
  }
  if (bad) {
    fprintf(stderr, "GPU self-test FAILED: %u of %u values differ, the device does not compute the filter exactly.\n",
            (unsigned)bad, (unsigned)t.expect.size());
    return false;
  }
  fprintf(stderr, "GPU self-test passed (%u values: redc, coef_mod and the symbol at each of the %u filter primes, %u polynomials).\n",
          (unsigned)t.expect.size(), t.np, t.npoly);
  return true;
}

// ---- device choice and failures -----------------------------------------------------

int nfGpuDeviceArg(int argc, char **argv)
{
  for (int i = 1; i + 1 < argc; i++)
    if (!strcmp(argv[i], "--device") || !strcmp(argv[i], "-device")) return atoi(argv[i + 1]);
  return -1;
}

#ifndef NF_GPU_NO_BOINC
NfGpuChoice nfGpuChoose(int argc, char **argv)
{
  APP_INIT_DATA aid;
  boinc_get_init_data(aid);
  NfGpuChoice c;
  c.type = aid.gpu_type;
  int arg = nfGpuDeviceArg(argc, argv);
  c.byHand = boinc_is_standalone() != 0;
  if (!c.type.empty()) {                       // assigned by the BOINC client, or by an init_data.xml copied in by hand
    c.device = aid.gpu_device_num >= 0 ? aid.gpu_device_num : (arg >= 0 ? arg : 0);
    return c;
  }
  if (!c.byHand && arg < 0)
    nfGpuFail(false, 1, "BOINC started this GPU program without a GPU.  Its app_version in app_info.xml needs a <coproc> "
                        "(use the app_info.xml of the package).");
  c.device = arg >= 0 ? arg : 0;               // run by hand (or a BOINC client older than 7.0.12): --device N, else 0
  return c;
}
#endif

void nfGpuFail(bool temporary, int code, const char *fmt, ...)
{
  char msg[1024];
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(msg, sizeof(msg), fmt, ap);
  va_end(ap);
  fprintf(stderr, "GPU mode: %s\n", msg);
#ifndef NF_GPU_NO_BOINC
  if (temporary && !boinc_is_standalone()) {
    fprintf(stderr, "This is usually temporary: BOINC restarts the task in a minute.\n");
    fflush(stderr);
    boinc_temporary_exit(60, msg);
  }
#else
  (void)temporary;
#endif
  fprintf(stderr, "The task is stopped.\n");
  fflush(stderr);
  exit(code);
}

// ---- init -------------------------------------------------------------------------

bool nfGpuActive() { return g_active; }

bool nfGpuInit(int argc, char **argv, int numP, const int *pSet)
{
  const char *force = getenv("NF_GPU");
  if (force && *force == '0') {                // developer switch: the GPU program's CPU path
    fprintf(stderr, "GPU mode off (NF_GPU=0): running on the CPU.\n");
    return false;
  }

  std::vector<int> primes = buildPrimeTable(numP, pSet);
  size_t np = primes.size() / 8;
  if (np < MIN_PRIMES) {                       // does not happen for |S| <= 3 (the program refuses larger S)
    if (!nfBackendCpuApp)
      nfGpuFail(false, 1, "only %u usable filter primes for this work unit; the GPU kernel needs %u.", (unsigned)np,
                (unsigned)MIN_PRIMES);
    fprintf(stderr, "GPU mode: only %u usable filter primes for this work unit; running on the CPU.\n", (unsigned)np);
    return false;
  }
  g_table = primes;
  if (!nfBackendInit(argc, argv, primes)) return false;   // only a CPU app version (MUSA) returns false

  const char *tg = getenv("NF_GPU_ITEMS");
  if (tg) g_target = strtoull(tg, NULL, 10);
  g_active = true;
  return true;
}

void nfGpuShutdown()
{
  if (!g_active) return;
  nfBackendShutdown();
  g_active = false;
}

// ---- queue ----------------------------------------------------------------------------

bool nfGpuWantFlush() { return g_nodes.size() + 4 >= NF_MAX_NODES || g_items >= g_target; }
bool nfGpuEmpty() { return g_nodes.empty(); }

void nfGpuAdd(const NfNode &node)
{
  long long lo = node.i[0], hi = node.i[1];
  // Very long m2 ranges are cut into consecutive chunks (keeps enumeration order).
  while ((unsigned long long)(hi - lo + 1) > CHUNK_ITEMS) {
    NfNode c = node;
    c.i[0] = lo;  c.i[1] = lo + (long long)CHUNK_ITEMS - 1;
    g_nodes.push_back(c);
    g_items += CHUNK_ITEMS;
    lo += (long long)CHUNK_ITEMS;
  }
  NfNode c = node;
  c.i[0] = lo;  c.i[1] = hi;
  g_nodes.push_back(c);
  g_items += (unsigned long long)(hi - lo + 1);
}

// ---- execution -----------------------------------------------------------------------------

static long long runRange(size_t a, size_t b, std::vector<NfSurv> &out)
{
  std::vector<unsigned long long> cum(b - a + 1);
  cum[0] = 0;
  for (size_t k = a; k < b; k++)
    cum[k - a + 1] = cum[k - a] + (unsigned long long)(g_nodes[k].i[1] - g_nodes[k].i[0] + 1);
  if (cum[b - a] == 0) return 0;

  bool overflow = false;
  long long inspected = nfBackendRun(a, b, cum.data(), out, &overflow);
  if (overflow) {
    if (b - a == 1) { fprintf(stderr, "GPU survivor buffer overflow on a single node.\n"); exit(1); }
    size_t mid = (a + b) / 2;
    return runRange(a, mid, out) + runRange(mid, b, out);
  }
  return inspected;
}

long long nfGpuRun(std::vector<long long> &polys)
{
  size_t n = g_nodes.size();
  if (n == 0) return 0;
  auto t0 = std::chrono::steady_clock::now();

  std::vector<double> nd(n * NF_ND);
  std::vector<long long> ni(n * NF_NI);
  for (size_t k = 0; k < n; k++) {
    memcpy(&nd[k * NF_ND], g_nodes[k].d, sizeof(double) * NF_ND);
    memcpy(&ni[k * NF_NI], g_nodes[k].i, sizeof(long long) * NF_NI);
  }
  nfBackendSetNodes(nd.data(), ni.data(), n);

  std::vector<NfSurv> surv;
  long long inspected = runRange(0, n, surv);

  std::sort(surv.begin(), surv.end(), [](const NfSurv &x, const NfSurv &y) {
    if (x.node != y.node) return x.node < y.node;
    if (x.m2 != y.m2) return x.m2 < y.m2;
    return x.m1 < y.m1;
  });
  for (const NfSurv &s : surv) {
    polys.push_back(1);
    for (int j = 0; j < 9; j++) polys.push_back(s.b[j]);
    polys.push_back(g_nodes[(size_t)s.node].i[34]);
  }

  // Aim for ~100 ms per launch (display drivers kill much longer ones).
  double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  double want = g_target * (100.0 / (ms > 1.0 ? ms : 1.0));
  if (want > 2.0 * g_target) want = 2.0 * g_target;
  if (want < 0.5 * g_target) want = 0.5 * g_target;
  if (want < 65536) want = 65536;
  if (want > 8.0e6) want = 8.0e6;
  if (!getenv("NF_GPU_ITEMS")) g_target = (unsigned long long)want;

  g_nodes.clear();
  g_items = 0;
  return inspected;
}
