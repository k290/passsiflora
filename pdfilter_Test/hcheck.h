// hcheck.h -- does hypothesis (H) hold for one polynomial?  Replays every int64 expression of upstream's first two
// sub-resultant iterations (polDiscTest_cpu.cpp at 638ca4e, the UNWRAP_TWO_ITERS branch) in the same order, with GCC's
// overflow builtins, which compute the wrapped result and report whether the exact one fits.  No undefined behaviour.
// Returns 1 if no expression leaves the int64 range (H holds), 0 otherwise.  c = polBuf row, leading coefficient first.
#pragma once
#include <stdint.h>

static inline int nf_h_ok(const long long *c)
{
  typedef long long i64;
  i64 A[11], B[10], R[10], t, u;
  int ov = 0;
#define MUL(r, x, y) (ov |= __builtin_mul_overflow((i64)(x), (i64)(y), &(r)))
#define SUB(r, x, y) (ov |= __builtin_sub_overflow((i64)(x), (i64)(y), &(r)))
  for (int col = 0; col < 11; col++) A[10 - col] = c[col];
  for (int k = 1; k <= 10; k++) MUL(B[k - 1], k, A[k]);                       // B[k-1] = k*A[k]
  int cont2 = 2;
  for (int k = 0; k < 10; k += 2) if ((B[k] & 1) == 1) { cont2 = 1; break; }
  int cont5 = 5;
  for (int k = 0; k < 9; k++) if ((B[k] % 5) != 0) { cont5 = 1; break; }
  int contB = cont2 * cont5;
  // B[k] = B[k] / contB, with the divisor as a constant (a 64-bit idiv per coefficient would cost more than all the rest)
  switch (contB) {
    case 2:  for (int k = 0; k < 10; k++) B[k] = B[k] / 2;  break;
    case 5:  for (int k = 0; k < 10; k++) B[k] = B[k] / 5;  break;
    case 10: for (int k = 0; k < 10; k++) B[k] = B[k] / 10; break;
    default: break;
  }
  MUL(R[0], B[9], A[0]);                                                       // R[0] = B[9]*A[0]
  for (int k = 1; k <= 9; k++) { MUL(t, B[9], A[k]); SUB(R[k], t, B[k - 1]); } // R[k] = B[9]*A[k] - B[k-1]
  for (int k = 0; k <= 8; k++) { MUL(t, B[9], R[k]); MUL(u, R[9], B[k]); SUB(R[k], t, u); }  // R[k] = B[9]*R[k] - R[9]*B[k]
  if (ov) return 0;                                                            // degR below would use wrapped values
  int degR = 0;
  for (int k = 8; k > 0; k--) if (R[k] != 0) { degR = k; break; }
  if (degR == 8) {
    for (int k = 0; k < 10; ++k) A[k] = B[k];
    for (int k = 0; k <= degR; ++k) B[k] = R[k];
    MUL(R[0], B[8], A[0]);
    for (int k = 1; k <= 8; ++k) { MUL(t, B[8], A[k]); MUL(u, A[9], B[k - 1]); SUB(R[k], t, u); }
    for (int k = 0; k <= 7; ++k) { MUL(t, B[8], R[k]); MUL(u, R[8], B[k]); SUB(R[k], t, u); }
    MUL(t, B[9], B[9]);                                                        // Rscale = B[9]*B[9]
  } else if (degR == 7) {
    for (int k = 0; k < 10; ++k) A[k] = B[k];
    for (int k = 0; k <= degR; ++k) B[k] = R[k];
    MUL(R[0], B[7], A[0]); MUL(R[1], B[7], A[1]);
    for (int k = 2; k <= 8; ++k) { MUL(t, B[7], A[k]); MUL(u, A[9], B[k - 2]); SUB(R[k], t, u); }
    MUL(R[0], B[7], R[0]);
    for (int k = 1; k <= 7; ++k) { MUL(t, B[7], R[k]); MUL(u, R[8], B[k - 1]); SUB(R[k], t, u); }
    for (int k = 0; k <= 6; ++k) { MUL(t, B[7], R[k]); MUL(u, R[7], B[k]); SUB(R[k], t, u); }
    MUL(t, B[9], B[9]); MUL(t, t, B[9]);                                       // Rscale = B[9]*B[9]*B[9]
    MUL(t, B[7], B[7]);                                                        // h = (B[7]*B[7])/B[9]
  }
  // degR < 7: upstream goes to GMP with these values; nothing more in int64.  The divisions R[k]/Rscale and
  // (B[7]*B[7])/B[9] cannot overflow (their divisors are positive resp. B[9] = 10/contB > 0).
#undef MUL
#undef SUB
  return !ov;
}

// The same check without data-dependent branches in the content and the division (exact division by 2 is an arithmetic
// shift, by 5 a multiplication by 5^-1 mod 2^64; both give B[k]/contB because contB divides every B[k]).  The two
// remainder degrees still branch, as upstream does.
static inline int nf_h_ok_fast(const long long *c)
{
  typedef long long i64; typedef unsigned long long u64;
  i64 a[11], B[10], R[10], t, u;
  int ov = 0;
#define MUL(r, x, y) (ov |= __builtin_mul_overflow((i64)(x), (i64)(y), &(r)))
#define SUB(r, x, y) (ov |= __builtin_sub_overflow((i64)(x), (i64)(y), &(r)))
  for (int k = 0; k <= 10; k++) a[k] = c[10 - k];                              // a[k]: coefficient of x^k
  for (int k = 1; k <= 10; k++) MUL(B[k - 1], k, a[k]);
  u64 odd = 0, n5 = 0;
  for (int k = 0; k < 10; k += 2) odd |= (u64)B[k];
  for (int k = 0; k < 9; k++) n5 |= (u64)((u64)B[k] * 0xCCCCCCCCCCCCCCCDull + 0x1999999999999999ull > 0x3333333333333332ull);  // 5 does not divide B[k]
  int s2 = !(odd & 1);                                                         // cont2 == 2
  u64 m5 = n5 ? 1ull : 0xCCCCCCCCCCCCCCCDull;                                  // cont5 == 5: multiply by 5^-1
  for (int k = 0; k < 10; k++) B[k] = (i64)((u64)(B[k] >> s2) * m5);
  MUL(R[0], B[9], a[0]);
  for (int k = 1; k <= 9; k++) { MUL(t, B[9], a[k]); SUB(R[k], t, B[k - 1]); }
  for (int k = 0; k <= 8; k++) { MUL(t, B[9], R[k]); MUL(u, R[9], B[k]); SUB(R[k], t, u); }
  if (ov) return 0;
  if (R[8] != 0) {                                                             // degR == 8
    i64 S[10];
    MUL(S[0], R[8], B[0]);
    for (int k = 1; k <= 8; ++k) { MUL(t, R[8], B[k]); MUL(u, B[9], R[k - 1]); SUB(S[k], t, u); }
    for (int k = 0; k <= 7; ++k) { MUL(t, R[8], S[k]); MUL(u, S[8], R[k]); SUB(S[k], t, u); }
  } else if (R[7] != 0) {                                                      // degR == 7
    i64 S[10];
    MUL(S[0], R[7], B[0]); MUL(S[1], R[7], B[1]);
    for (int k = 2; k <= 8; ++k) { MUL(t, R[7], B[k]); MUL(u, B[9], R[k - 2]); SUB(S[k], t, u); }
    MUL(S[0], R[7], S[0]);
    for (int k = 1; k <= 7; ++k) { MUL(t, R[7], S[k]); MUL(u, S[8], R[k - 1]); SUB(S[k], t, u); }
    for (int k = 0; k <= 6; ++k) { MUL(t, R[7], S[k]); MUL(u, S[7], R[k]); SUB(S[k], t, u); }
    MUL(t, R[7], R[7]);
  }
  // B[9] = 10/contB <= 10, so Rscale = B[9]^2 or B[9]^3 cannot overflow.
#undef MUL
#undef SUB
  return !ov;
}
