// test_tiny.cpp -- the compiled filter against an exact oracle on complete sets of polynomials at a tiny prime
// (pdfilter_chains.tex, Section "The compiled code at tiny primes").  tiny_checks.sh builds src/pdfilter.f90 with one tiny
// prime in its tables (the search constants QS_START ... changed by sed; nothing else) and runs this program in parallel.
//
// Usage: test_tiny Q KIND MODE START COUNT
//   Q      the prime of the filter build (13, 17, 29, 37, 41: = 1 mod 4 and not dividing n = 10); it is checked against the hook
//   KIND   small | large: the stage the build has (chain_block_s/chain_general_s with mont, or chain_block/chain_general with sredc)
//   MODE   depressed  the monic polynomials with a9 = 0, numbered START .. START+COUNT-1 (base-Q digits = a0 .. a8); every monic
//                     polynomial is a shift f(x+t) of exactly one of them, since Q does not divide 10
//          zerosK     the same with a9 = ... = a(10-K) = 0 (K = 2..9): the first pivot is 0, every chain is abnormal at step 1
//          shifts     COUNT random depressed polynomials (seed START), each with all Q shifts f(x+t), every coefficient
//                     lifted by a random multiple of Q below 2^40 (the residues are what counts; the loading sees large values)
// Every polynomial goes through four runs: the AVX2 and the baseline build, each with the normal routing and with
// NF_GENERAL_ALL set (every lane through chain_general), in calls of varying size (block sizes 1 .. 4096).  pdf_init is called
// before every call, so every call is the first of a routing history (a probe): with the normal routing pass 1 runs
// chain_block on every lane and pass 2 runs chain_general on the lanes that left it.  The flag and the
// stage from the test hook are compared with chi(disc f mod Q), computed by Euclid's algorithm with inverses in F_Q.
// Exit status 1 on any disagreement.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <set>
#include <cstdint>

extern "C" void pdf_init_v3(int, const int*);
extern "C" void pdf_filter_v3(const long long*, int, signed char*);
extern "C" void pdf_init_base(int, const int*);
extern "C" void pdf_filter_base(const long long*, int, signed char*);
extern "C" void pdf_stages_v3(int, int*, long long*, long long*, int*, int);
extern "C" void pdf_stages_base(int, int*, long long*, long long*, int*, int);

typedef long long i64;
static int Q;

static int md(i64 x) { x %= Q; if (x < 0) x += Q; return (int)x; }
static int pw(int b, int e) { i64 r = 1, x = md(b); while (e) { if (e & 1) r = r * x % Q; x = x * x % Q; e >>= 1; } return (int)r; }
static int inv(int a) { return pw(a, Q - 2); }

// chi(disc f) for f = x^10 + ... with f[k] the coefficient of x^k (residues), from Res(f, f') by Euclid with inverses:
// Res(A, B) = +-lc(B)^(deg A - deg R) Res(B, R) for R = A mod B, Res(A, c) = c^deg A; disc f = +-Res(f, f'), and the signs do
// not matter because Q = 1 mod 4.  *key: the degrees of the remainder sequence (bit d) and bit 11 if it ends with 0.
static int oracle(const int *f, unsigned *key)
{
  int A[11], B[11], da = 10, db = 9;
  for (int k = 0; k <= 10; k++) A[k] = f[k];
  for (int k = 1; k <= 10; k++) B[k - 1] = md((i64)k * f[k]);
  i64 res = 1;
  *key = (1u << 10) | (1u << 9);
  for (;;) {
    if (db == 0) { res = res * pw(B[0], da) % Q; break; }
    int il = inv(B[db]);
    for (int i = da; i >= db; i--) {          // A = A mod B
      int t = (int)((i64)A[i] * il % Q);
      if (t) for (int k = 0; k <= db; k++) A[i - db + k] = md(A[i - db + k] - (i64)t * B[k]);
    }
    int dr = db - 1;
    while (dr >= 0 && A[dr] == 0) dr--;
    if (dr < 0) { *key |= 1u << 11; return 0; }
    *key |= 1u << dr;
    res = res * pw(B[db], da - dr) % Q;
    int T[11]; memcpy(T, A, sizeof T); memcpy(A, B, sizeof A); memcpy(B, T, sizeof B);
    da = db; db = dr;
  }
  if (res == 0) return 0;
  return pw((int)res, (Q - 1) / 2) == 1 ? 1 : -1;
}

struct Variant {
  const char *name; void (*init)(int, const int*); void (*filter)(const long long*, int, signed char*);
  void (*stages)(int, int*, long long*, long long*, int*, int); bool general;
  i64 polys = 0, rejected = 0, falserej = 0, missed = 0, wrongstage = 0, cnt[7] = {0};
};

static uint64_t rs;
static uint64_t rnd() { rs += 0x9E3779B97F4A7C15ull; uint64_t z = rs; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }

int main(int argc, char **argv)
{
  if (argc < 6) { fprintf(stderr, "usage: test_tiny Q small|large depressed|zerosK|shifts START COUNT\n"); return 2; }
  Q = atoi(argv[1]);
  bool large = !strcmp(argv[2], "large");
  const char *mode = argv[3];
  i64 start = atoll(argv[4]), count = atoll(argv[5]);
  int zeros = 1;                                   // a9 = 0 always; zerosK: a9 .. a(10-K) = 0
  bool shifts = !strcmp(mode, "shifts");
  if (!strncmp(mode, "zeros", 5)) zeros = atoi(mode + 5);
  int free_ = 10 - zeros;                          // coefficients a0 .. a(free_-1) run over all residues
  Variant V[4] = { {"AVX2 normal ", pdf_init_v3, pdf_filter_v3, pdf_stages_v3, false},
                   {"AVX2 general", pdf_init_v3, pdf_filter_v3, pdf_stages_v3, true},
                   {"base normal ", pdf_init_base, pdf_filter_base, pdf_stages_base, false},
                   {"base general", pdf_init_base, pdf_filter_base, pdf_stages_base, true} };
  const int sizes[] = {4096, 1, 7, 63, 64, 65, 127, 255, 1000, 2048, 4096, 64, 1024};
  const int nsizes = sizeof sizes / sizeof sizes[0];
  int dummy[1] = {0};
  long long counts[7 * 30], primes[30];
  int nst;
  // the build must hold exactly the prime Q at the right stage
  for (Variant &v : V) {
    unsetenv("NF_GENERAL_ALL");
    v.init(0, dummy);
    v.stages(0, nullptr, counts, primes, &nst, 1);
    if (nst != 1 || primes[0] != Q) { fprintf(stderr, "the filter build does not hold the single prime %d (stages %d, first %lld)\n", Q, nst, primes[0]); return 2; }
  }
  rs = (uint64_t)start * 0x2545F4914F6CDD1Dull + 12345;
  std::vector<i64> pol; std::vector<int> chi; std::vector<signed char> flag; std::vector<int> st;
  std::set<unsigned> keys; i64 abnormal = 0, zero = 0, total = 0, minus = 0;
  i64 idx = start, done = 0, todo = shifts ? count * Q : count; int bi = 0;
  i64 nshift = 0; int fbase[11];                   // shifts mode: the current depressed polynomial and the next shift
  while (done < todo) {
    int n = (int)std::min<i64>(sizes[bi++ % nsizes], todo - done);
    pol.assign((size_t)n * 11, 0); chi.assign(n, 0); flag.assign(n, 0); st.assign(n, 0);
    for (int i = 0; i < n; i++) {
      int f[11];
      if (!shifts) {
        i64 x = idx++;
        for (int k = 0; k < 10; k++) {
          if (k < free_) { f[k] = (int)(x % Q); x /= Q; } else f[k] = 0;
        }
        f[10] = 1;
      } else {
        if (nshift == 0) { for (int k = 0; k < 9; k++) fbase[k] = (int)(rnd() % Q); fbase[9] = 0; fbase[10] = 1; }
        // f(x + t) by repeated synthetic division (Taylor shift), t = nshift
        int t = (int)nshift, g[11];
        memcpy(g, fbase, sizeof g);
        for (int i0 = 0; i0 < 10; i0++) for (int k = 9; k >= i0; k--) g[k] = md(g[k] + (i64)t * g[k + 1]);
        memcpy(f, g, sizeof g);
        nshift = (nshift + 1) % Q;
      }
      unsigned key; chi[i] = oracle(f, &key);
      keys.insert(key); total++;
      if (key != 0x7FFu) abnormal++;
      if (chi[i] == 0) zero++;
      if (chi[i] == -1) minus++;
      long long *p = &pol[(size_t)i * 11];
      p[0] = 1;
      for (int col = 1; col <= 10; col++) {
        long long c = f[10 - col] - (Q - 1) / 2 >= 0 ? f[10 - col] - Q : f[10 - col];   // symmetric representative
        if (shifts) c += (long long)Q * (long long)((int64_t)(rnd() % (1ull << 40)) - (1ll << 39));
        p[col] = c;
      }
    }
    for (Variant &v : V) {
      if (v.general) setenv("NF_GENERAL_ALL", "1", 1); else unsetenv("NF_GENERAL_ALL");
      v.init(0, dummy);
      v.filter(pol.data(), n, flag.data());
      v.stages(n, st.data(), counts, primes, &nst, 1);
      for (int k = 0; k < 7; k++) v.cnt[k] += counts[k];
      for (int i = 0; i < n; i++) {
        bool rej = flag[i] == 0, want = chi[i] == -1;
        v.polys++; v.rejected += rej;
        if (rej && !want) { v.falserej++; if (v.falserej <= 3) fprintf(stderr, "%s: false rejection, poly index %lld\n", v.name, (long long)(idx - n + i)); }
        if (!rej && want) { v.missed++; if (v.missed <= 3) fprintf(stderr, "%s: missed rejection, poly index %lld\n", v.name, (long long)(idx - n + i)); }
        if (st[i] != (rej ? 1 : 0)) v.wrongstage++;
      }
    }
    done += n;
  }
  // machine-readable lines (tiny_checks.sh sums them over the chunks)
  bool bad = false;
  const char *kind = large ? "large" : "small";
  printf("ORACLE %d %s %s %lld %lld %lld %lld\n", Q, kind, mode, (long long)total, (long long)minus, (long long)abnormal, (long long)zero);
  for (Variant &v : V) {
    printf("RUN %d %s %s %s %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld\n", Q, kind, mode, v.name[0] == 'A' ? (v.general ? "AVX2-general" : "AVX2-normal") :
           (v.general ? "base-general" : "base-normal"), (long long)v.polys, (long long)v.rejected, (long long)v.falserej, (long long)v.missed,
           (long long)v.wrongstage, (long long)v.cnt[1], (long long)v.cnt[2], (long long)v.cnt[3], (long long)v.cnt[5], (long long)v.cnt[6]);
    if (v.falserej || v.missed || v.wrongstage || v.rejected != minus) bad = true;
  }
  for (unsigned k : keys) printf("SEQ %d %s %s %x\n", Q, kind, mode, k);
  return bad ? 1 : 0;
}
