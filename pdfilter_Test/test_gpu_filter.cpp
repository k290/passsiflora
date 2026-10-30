// test_gpu_filter.cpp -- the GPU kernel's pre-filter (survives / chi_reject / redc / coef_mod of src/nf_kernel.h, with the prime
// table of src/nf_gpu_common.cpp) compiled for the CPU and compared with an exact oracle (pdfilter_chains.tex, Theorem E): for every
// polynomial and every prime the kernel must reject at q exactly when chi(disc f mod q) = -1.
//   TEST=test_gpu_filter.cpp bash build-linux.sh
//   build/test_gpu_filter COUNT SEED [THREADS] [tiny [TOP]]
// The table is the kernel's for S = {2, 5}: the 24 smallest primes q = 1 mod 4 in (256, 8192) at which 2 and 5 are squares.  "tiny" uses
// instead primeRow for every prime q = 1 mod 4 in [13, 199] and the 8 largest such primes up to TOP (default 8191, the limit of
// Theorem E, set by redc; a larger TOP shows what happens beyond it).  The kernel's sources are included,
// not copied; src/nf_port_cpu.h is the CPU port, in the manner of nf_port_opencl.h and nf_port_cuda.h.  First, for a table below
// 8192, the start-up self-test of the GPU mode (nf_selftest_t2 and nfSelfTestBuild / nfSelfTestVerify) runs on the CPU: its
// expected values, computed without the kernel's code, must agree with the kernel compiled here.
// Families (round robin): 0 random coefficients in [-10^6, 10^6]; 1 any int64; 2 depressed with |a_i| <= 3 (abnormal chains);
// 3 x^10 + b x^7 + c x^4 + d x (b, c, d in [-10^6, 10^6]); 4 coefficients that are multiples of one table prime with probability 1/2;
// 5 (x - r)^2 g(x) + q h(x) for one table prime q (disc = 0 mod q).  Exit status 1 on any mismatch.
#include <cstdint>
#include <cmath>
#include <math.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <thread>
#include <mutex>

#include "../src/nf_gpu_common.cpp"      // with the kernel, compiled for the CPU through src/nf_port_cpu.h
// the back end is not used: buildPrimeTable and primeRow are called directly
bool nfBackendInit(int, char **, const std::vector<int> &) { return false; }
void nfBackendShutdown() {}
void nfBackendSetNodes(const double *, const long long *, size_t) {}
long long nfBackendRun(size_t, size_t, const unsigned long long *, std::vector<NfSurv> &, bool *) { return 0; }

struct Rng { u64 s; u64 next() { s += 0x9E3779B97F4A7C15ull; u64 z = s; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }
             i64 range(i64 b) { return (i64)(next() % (u64)(2 * b + 1)) - b; } };

static i64 pw(i64 b, i64 e, i64 m) { i64 r = 1; b %= m; if (b < 0) b += m; while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; } return r; }

// Res(A, B) over F_q by Euclid with inverses (the recursion of Lemma "Resultants over F"); A, B coefficient vectors, index = degree.
static i64 resq(std::vector<i64> A, std::vector<i64> B, i64 q)
{
  auto deg = [](const std::vector<i64> &P) { int d = (int)P.size() - 1; while (d >= 0 && P[d] == 0) d--; return d; };
  i64 r = 1;
  for (;;) {
    int a = deg(A), b = deg(B);
    if (b < 0) return 0;
    if (b == 0) return r * pw(B[0], a, q) % q;
    std::vector<i64> R(A.begin(), A.begin() + a + 1);
    i64 inv = pw(B[b], q - 2, q);
    for (int k = a; k >= b; k--) {
      i64 f = R[k] * inv % q;
      if (f) for (int j = 0; j <= b; j++) R[k - b + j] = ((R[k - b + j] - f * B[j]) % q + q) % q;
    }
    R.resize(b > 0 ? b : 1);
    int rd = deg(R);
    if (rd < 0) return 0;
    // Res(A, B) = (-1)^(ab) lc(B)^(a - deg R) Res(B, R)
    if ((a * b) & 1) r = (q - r) % q;
    r = r * pw(B[b], a - rd, q) % q;
    A = B; B = R;
  }
}

// chi(disc f mod q) for f = x^10 + c[1] x^9 + ... + c[10] (q = 1 mod 4: chi(disc) = chi(Res(f, f')))
static int chi_oracle(const i64 *c, i64 q)
{
  std::vector<i64> A(11), B(10);
  for (int k = 0; k <= 10; k++) A[k] = k == 10 ? 1 : (((c[10 - k] % q) + q) % q);
  for (int k = 1; k <= 10; k++) B[k - 1] = (k * A[k]) % q;
  i64 r = resq(A, B, q);
  return r == 0 ? 0 : pw(r, (q - 1) / 2, q) == 1 ? 1 : -1;
}

int main(int argc, char **argv)
{
  if (argc < 3) { fprintf(stderr, "usage: test_gpu_filter COUNT SEED [THREADS] [tiny]\n"); return 2; }
  i64 count = atoll(argv[1]); u64 seed = strtoull(argv[2], 0, 10);
  int threads = argc > 3 ? atoi(argv[3]) : 1;
  bool tiny = argc > 4 && !strcmp(argv[4], "tiny");
  std::vector<int> table;
  if (tiny) {
    std::vector<long long> qs;
    for (long long q = 13; q < 200; q += 4) if (isPrime(q)) qs.push_back(q);
    long long top = argc > 5 ? atoll(argv[5]) : 8191;       // the largest q of Theorem E (redc: (t + m q) / 2^16 < 2q needs q < 8192)
    int big = 0; for (long long q = top - (top - 1) % 4; big < 8; q -= 4) if (isPrime(q)) { qs.push_back(q); big++; }
    for (long long q : qs) { int row[8]; primeRow(q, row); table.insert(table.end(), row, row + 8); }
  } else {
    int pSet[2] = {2, 5};
    table = buildPrimeTable(2, pSet);
  }
  int np = (int)table.size() / 8;
  printf("primes %d:", np); for (int p = 0; p < np; p++) printf(" %d", table[8 * p]); printf("\n");
  bool small = true;
  for (int p = 0; p < np; p++) small = small && table[8 * p] < 8192;
  if (small) {                                   // the start-up self-test, with the kernel compiled for the CPU as the device
    NfSelfTest st;
    nfSelfTestBuild(table, st);
    std::vector<unsigned> out(st.expect.size());
    for (size_t g = 0; g < out.size(); g++) {
      nf_cpu_gid = g;
      nf_selftest_t2(st.ts.data(), st.nt, (const i64 *)st.xs.data(), st.nx, (const i64 *)st.polys.data(), st.npoly, table.data(),
                     (nfuint)np, out.data(), (nfuint)out.size());
    }
    bool ok = nfSelfTestVerify(st, out);
    printf("SELFTEST %s: %u values\n", ok ? "passed" : "FAILED", (unsigned)out.size());
    if (!ok) return 1;
  }
  const int NF = 6;
  i64 polys[NF] = {0}, checks[NF] = {0}, rej[NF] = {0}, minus[NF] = {0}, zero[NF] = {0}, falsr[NF] = {0}, missed[NF] = {0}, incons[NF] = {0};
  std::mutex mx;
  auto work = [&](int t) {
    Rng rng{seed * 1000003ull + (u64)t};
    i64 P[NF] = {0}, C[NF] = {0}, Rj[NF] = {0}, M[NF] = {0}, Z[NF] = {0}, Fa[NF] = {0}, Mi[NF] = {0}, In[NF] = {0};
    for (i64 it = t; it < count; it += threads) {
      int fam = (int)(it % NF);
      i64 c[11]; c[0] = 1;
      i64 qs = table[8 * (int)(rng.next() % np)];
      switch (fam) {
        case 0: for (int k = 1; k <= 10; k++) c[k] = rng.range(1000000); break;
        case 1: for (int k = 1; k <= 10; k++) c[k] = (i64)rng.next(); break;
        case 2: c[1] = 0; for (int k = 2; k <= 10; k++) c[k] = rng.range(3); break;
        case 3: for (int k = 1; k <= 10; k++) c[k] = (k == 3 || k == 6 || k == 9) ? rng.range(1000000) : 0; break;
        case 4: for (int k = 1; k <= 10; k++) c[k] = (rng.next() & 1) ? qs * rng.range(1000000) : rng.range(1000000); break;
        default: {                                   // (x - r)^2 g + q h, g monic of degree 8
          i64 r = rng.range(1000), g[9], f[11] = {0};
          g[0] = 1; for (int k = 1; k <= 8; k++) g[k] = rng.range(1000);      // leading first
          i64 s[3] = {1, -2 * r, r * r};
          for (int u = 0; u < 3; u++) for (int v = 0; v <= 8; v++) f[u + v] += s[u] * g[v];
          for (int k = 1; k <= 10; k++) c[k] = f[k] + qs * rng.range(1000);
        }
      }
      P[fam]++;
      bool any = false;
      for (int p = 0; p < np; p++) {
        i64 q = table[8 * p];
        int chi = chi_oracle(c, q);
        bool r = !survives(c, &table[8 * p], 1);
        C[fam]++; Rj[fam] += r; M[fam] += chi == -1; Z[fam] += chi == 0; any |= r;
        if (r && chi != -1) { Fa[fam]++; if (Fa[fam] <= 3) { std::lock_guard<std::mutex> g(mx); fprintf(stderr, "FALSE REJECTION q=%lld chi=%d family %d\n", (long long)q, chi, fam); } }
        if (!r && chi == -1) { Mi[fam]++; if (Mi[fam] <= 3) { std::lock_guard<std::mutex> g(mx); fprintf(stderr, "MISSED REJECTION q=%lld family %d\n", (long long)q, fam); } }
      }
      if ((survives(c, table.data(), (nfuint)np) != 0) == any) In[fam]++;
    }
    std::lock_guard<std::mutex> g(mx);
    for (int f = 0; f < NF; f++) { polys[f] += P[f]; checks[f] += C[f]; rej[f] += Rj[f]; minus[f] += M[f]; zero[f] += Z[f]; falsr[f] += Fa[f]; missed[f] += Mi[f]; incons[f] += In[f]; }
  };
  std::vector<std::thread> th;
  for (int t = 0; t < threads; t++) th.emplace_back(work, t);
  for (auto &x : th) x.join();
  i64 bad = 0;
  for (int f = 0; f < NF; f++) {
    printf("GPUFILTER %s family %d: polys %lld checks %lld rejections %lld chi=-1 %lld chi=0 %lld false %lld missed %lld survives-inconsistent %lld\n",
           tiny ? "tiny" : "table", f, (long long)polys[f], (long long)checks[f], (long long)rej[f], (long long)minus[f], (long long)zero[f],
           (long long)falsr[f], (long long)missed[f], (long long)incons[f]);
    bad += falsr[f] + missed[f] + incons[f];
  }
  return bad ? 1 : 0;
}
