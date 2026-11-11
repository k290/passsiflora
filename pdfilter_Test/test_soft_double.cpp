// test_soft_double.cpp -- the GPU kernel's software doubles (src/nf_kernel.h with NF_SOFT_DOUBLE, included, not copied)
// against the CPU's IEEE doubles, bit for bit.  A device without double precision (no cl_khr_fp64) runs the enumeration's
// bounds with these functions, so they must give exactly the hardware's results: every operation the kernel uses (add,
// subtract, multiply, divide, square root, conversion from int64, ceil, floor, the x86 rule of double -> int64, <= and >=)
// is compared with the CPU's own operation, NaN against NaN.
//   TEST=test_soft_double.cpp bash build-linux.sh && build/test_soft_double [MILLIONS] [SEED]
// MILLIONS (default 10) random operand pairs per operation, plus every pair of special values (zeros, infinities, NaN,
// subnormal and normal limits, 2^52, 2^63, ...).  The operands are drawn to reach every path: any bit pattern, values near 1,
// subnormals, large and small exponents, integers of every size, short significands (exact results and ties), operands of
// equal exponent (cancellation) and nearly equal operands; conversions from int64 include values above 2^53 that need
// rounding, exact ties and INT64_MIN.  The reference is the compiler's double arithmetic: build without x87 (x86-64 uses
// SSE2) and without fused multiply-add.  Exit status 1 on any difference.
#define NF_SOFT_DOUBLE 1
#include "../src/nf_port_cpu.h"
#include "../src/nf_kernel.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <random>

static double bd(u64 b) { double d; memcpy(&d, &b, 8); return d; }
static u64 db(double d) { u64 b; memcpy(&b, &d, 8); return b; }
static bool same(u64 a, u64 b) { return a == b || (sd_isnan(a) && sd_isnan(b)); }
// the hardware version of the kernel's d2l (nf_kernel.h without NF_SOFT_DOUBLE)
static i64 hw_d2l(double x) { if (!(x > -9223372036854775808.0 && x < 9223372036854775808.0)) return INT64_MIN; return (i64)x; }

static std::mt19937_64 rng;
enum { ADD, SUB, MUL, DIV, SQRT, FROMI, CEIL, FLOOR, D2L, LE, GE, NOPS };
static long long fails[NOPS], runs[NOPS];
static const char *nm[NOPS] = { "add", "sub", "mul", "div", "sqrt", "from_i64", "ceil", "floor", "d2l", "le", "ge" };
static const u64 SIGN_FRAC = 0x800FFFFFFFFFFFFFULL;

static u64 pick()
{
  switch (rng() & 7) {
    case 0: return rng();                                                                  // any bit pattern
    case 1: return (rng() & SIGN_FRAC) | ((u64)(1023 - 32 + rng() % 64) << 52);            // near 1
    case 2: return rng() & SIGN_FRAC;                                                      // subnormal or zero
    case 3: return (rng() & SIGN_FRAC) | ((u64)(1023 + rng() % 120) << 52);               // up to 2^120
    case 4: return db((double)(i64)rng());                                                 // integers of every size
    case 5: return db((double)(int)(rng() % 2001) - 1000.0 + (double)(rng() % 8) / 8.0);  // small values with halves
    case 6: return (rng() & SIGN_FRAC) | ((u64)(rng() % 2046 + 1) << 52);                 // any normal number
    default: return (rng() & 0x8000000000000FFFULL) | ((u64)(1000 + rng() % 100) << 52); // short significands
  }
}

static void report(int op, u64 a, u64 b, long long got, long long want)
{
  if (fails[op]++ < 5)
    printf("%s(%016llx, %016llx): soft %016llx, hardware %016llx\n", nm[op], (unsigned long long)a, (unsigned long long)b,
           (unsigned long long)got, (unsigned long long)want);
}
static void check(int op, u64 got, u64 want, u64 a, u64 b) { runs[op]++; if (!same(got, want)) report(op, a, b, (long long)got, (long long)want); }
static void checki(int op, long long got, long long want, u64 a, u64 b) { runs[op]++; if (got != want) report(op, a, b, got, want); }

static void pair(u64 a, u64 b)
{
  double x = bd(a), y = bd(b);
  check(ADD, sd_add(a, b), db(x + y), a, b);
  check(SUB, sd_sub(a, b), db(x - y), a, b);
  check(MUL, sd_mul(a, b), db(x * y), a, b);
  check(DIV, sd_div(a, b), db(x / y), a, b);
  check(SQRT, sd_sqrt(a), db(sqrt(x)), a, 0);
  check(CEIL, sd_round_int(a, 1), db(ceil(x)), a, 0);
  check(FLOOR, sd_round_int(a, 0), db(floor(x)), a, 0);
  checki(D2L, sd_d2l(a), hw_d2l(x), a, 0);
  checki(LE, sd_le(a, b), x <= y, a, b);
  checki(GE, sd_le(b, a), x >= y, a, b);
}

static i64 pick_int(long long i)
{
  switch (i & 3) {
    case 0: return (i64)rng();
    case 1: { i64 v = (i64)(rng() >> (rng() % 64)); return (rng() & 1) ? -v : v; }
    case 2: return (i64)(((rng() & 0x1FFFFFFFFFFFFFULL) << (rng() % 11)) | ((rng() & 1) ? (1ULL << (rng() % 11)) : 0));  // ties
    default: { static const i64 e[] = { 0, 1, -1, INT64_MIN, INT64_MAX, (1LL << 53) + 1, -(1LL << 53) - 1, (1LL << 62) + (1LL << 9) };
               return e[rng() % 8]; }
  }
}

int main(int argc, char **argv)
{
  long long n = (argc > 1 ? atoll(argv[1]) : 10) * 1000000LL;
  rng.seed(argc > 2 ? strtoull(argv[2], NULL, 10) : 12345);
  static const u64 sp[] = { 0, 0x8000000000000000ULL, 0x7FF0000000000000ULL, 0xFFF0000000000000ULL, 0x7FF8000000000000ULL,
                            0x0000000000000001ULL, 0x000FFFFFFFFFFFFFULL, 0x0010000000000000ULL, 0x7FEFFFFFFFFFFFFFULL,
                            0x3FF0000000000000ULL, 0xBFF0000000000000ULL, 0x3FE0000000000000ULL, 0x4014000000000000ULL,
                            0x43E0000000000000ULL, 0xC3E0000000000000ULL, 0x4340000000000000ULL, 0x433FFFFFFFFFFFFFULL,
                            0x3FEFFFFFFFFFFFFFULL, 0x8010000000000000ULL, 0x800FFFFFFFFFFFFFULL };
  for (u64 a : sp) for (u64 b : sp) pair(a, b);
  for (long long i = 0; i < n; i++) {
    u64 a = pick(), b = pick();
    if ((i & 3) == 0) b = (b & SIGN_FRAC) | (a & 0x7FF0000000000000ULL);  // equal exponents: cancellation
    if ((i & 7) == 1) b = a ^ (rng() & 0x3FULL);                          // nearly equal
    pair(a, b);
    i64 v = pick_int(i);
    checki(FROMI, (long long)sd_from_i64(v), (long long)db((double)v), (u64)v, 0);
  }
  long long total = 0;
  for (int k = 0; k < NOPS; k++) { printf("%-9s %12lld cases  %lld different\n", nm[k], runs[k], fails[k]); total += fails[k]; }
  printf(total ? "FAILED: the software doubles differ from the hardware\n" : "OK: every result equal to the hardware's\n");
  return total != 0;
}
