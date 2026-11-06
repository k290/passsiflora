// chain_trace.cpp -- the two Euclidean chains of the pre-filter, for any degree n and prime q, with a printed trace.
//
//   chain_block    follows subroutine chain_block of src/pdfilter.f90, for one polynomial
//   chain_general  follows subroutine chain_general of src/pdfilter.f90, for one polynomial
//
// It is for reading, not for speed.  Numbers are reduced to the range -(q-1)/2 ... (q-1)/2, as red() does.
//
// Usage:
//   chain_trace show q c_n ... c_0 [block | blockfull | general]
//        traces one monic polynomial, coefficients from the leading one down (c_n = 1); `block` stops after the first failed check,
//        `blockfull` does not, no third word traces both chains.      e.g.   chain_trace show 41 1 0 3 2 8 block      (x^4 + 3x^2 + 2x + 8 mod 41)
//   chain_trace test
//        compares both chains with the exact discriminant (PARI) for random polynomials of degree 4 and 10 and 13 primes; exit status 1 if any differs
//
// Build (this directory):   TEST=chain_trace.cpp bash build-linux.sh      (build-windows.sh on Windows)      -> build/chain_trace
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
// i64 is declared before pari.h: on 64-bit Windows pari.h redefines `long` as `long long`
typedef long long i64;
#include "pari.h"

static i64 Q = 41;
static i64 rd(i64 x) { x %= Q; if (x < 0) x += Q; if (x > (Q - 1) / 2) x -= Q; return x; }   // symmetric residue
static i64 powm(i64 b, i64 e) { i64 r = 1; b = rd(b); while (e > 0) { if (e & 1) r = rd(r * b); b = rd(b * b); e >>= 1; } return r; }

static std::string poly_asc(const std::vector<i64> &c, int deg)       // c[k] = coefficient of x^k
{
  std::string s; bool first = true;
  for (int k = deg; k >= 0; k--) {
    i64 v = c[k]; if (v == 0 && deg != 0) continue;
    if (!first) s += v < 0 ? " - " : " + "; else if (v < 0) s += "-";
    i64 a = v < 0 ? -v : v;
    if (k == 0) s += std::to_string(a);
    else { if (a != 1) s += std::to_string(a); s += "x"; if (k > 1) s += "^" + std::to_string(k); }
    first = false;
  }
  if (first) s = "0";
  return s;
}
static std::string list_desc(const std::vector<i64> &c, int len)
{
  std::string s = "["; for (int i = 0; i < len; i++) { if (i) s += ", "; s += std::to_string(c[i]); } return s + "]";
}

// ---- chain_block: a[k] = coefficient of x^k, a[n] = 1
struct BlockResult { bool ok; i64 base, acc; };
static BlockResult chain_block(int n, const std::vector<i64> &a, bool trace, bool stopAfterFail)
{
  std::vector<i64> A(n + 2, 0), B(n + 2, 0), R(n + 2, 0), R2(n + 2, 0);
  for (int k = 0; k <= n; k++) A[k] = rd(a[k]);
  for (int k = 1; k <= n; k++) B[k - 1] = rd(k * A[k]);
  bool ok = true;
  if (trace) printf("**Start.** A = f = %s and B = A' = %s.\n\n", poly_asc(A, n).c_str(), poly_asc(B, n - 1).c_str());
  for (int m = n - 1; m >= 1; m--) {
    i64 lc = B[m], top = A[m + 1];
    if (trace) printf("#### Step m = %d\n\nA = %s (degree %d), B = %s (degree %d).  `lc = B[%d] = %lld`, `top = A[%d] = %lld`.\n\n", m, poly_asc(A, m + 1).c_str(), m + 1, poly_asc(B, m).c_str(), m, m, lc, m + 1, top);
    std::fill(R.begin(), R.end(), 0);
    R[0] = rd(lc * A[0]);
    if (trace) printf("**Round 1.**  R = lc×A - top×x×B, that is R[0] = lc×A[0] and R[k] = lc×A[k] - top×B[k-1] for k = 1..%d:\n\n| k | lc×A[k] | top×B[k-1] | difference | R[k] mod %lld |\n|---|---|---|---|---|\n| 0 | %lld×%lld = %lld | - | %lld | %lld |\n", m, Q, lc, A[0], lc * A[0], lc * A[0], R[0]);
    for (int k = 1; k <= m; k++) {
      R[k] = rd(lc * A[k] - top * B[k - 1]);
      if (trace) printf("| %d | %lld×%lld = %lld | %lld×%lld = %lld | %lld | %lld |\n", k, lc, A[k], lc * A[k], top, B[k - 1], top * B[k - 1], lc * A[k] - top * B[k - 1], R[k]);
    }
    if (trace) printf("\nSo R = %s.\n\n", poly_asc(R, m).c_str());
    i64 t = R[m];
    R2 = R;
    for (int k = 0; k <= m - 1; k++) R2[k] = rd(lc * R[k] - t * B[k]);
    R2[m] = 0;
    if (trace) {
      printf("**Round 2.**  t = R[%d] = %lld;  R2 = lc×R - t×B, that is R2[k] = lc×R[k] - t×B[k] for k = 0..%d:\n\n| k | lc×R[k] | t×B[k] | difference | R2[k] mod %lld |\n|---|---|---|---|---|\n", m, t, m - 1, Q);
      for (int k = 0; k <= m - 1; k++) printf("| %d | %lld×%lld = %lld | %lld×%lld = %lld | %lld | %lld |\n", k, lc, R[k], lc * R[k], t, B[k], t * B[k], lc * R[k] - t * B[k], R2[k]);
      printf("\nSo R2 = %s, of degree at most %d.\n\n", poly_asc(R2, m - 1).c_str(), m - 1);
    }
    bool failedNow = false;
    if (m >= 2) {
      bool nz = R2[m - 1] != 0;
      failedNow = !nz && ok;
      ok = ok && nz;
      if (trace) printf("**Check** (pdfilter.f90 line 303): the new leading coefficient R2[%d] = %lld is %s, so `ok` %s.\n\n", m - 1, R2[m - 1], nz ? "not zero" : "ZERO", ok ? "stays true" : "is now FALSE");
    }
    if (trace && failedNow && stopAfterFail) { printf("(The loop would carry on with `lc = 0` and meaningless numbers, but `ok` stays false; the result is not used.)\n\n"); BlockResult r; r.ok = false; r.base = 0; r.acc = 0; return r; }
    if (trace && m > 1) printf("Next step: A := B, B := R2.\n\n");
    A = B; B = R2;
  }
  BlockResult r; r.ok = ok; r.base = B[0]; r.acc = powm(B[0], (Q - 1) / 2);
  if (trace) printf("**End.**  The final constant is %lld.  Euler's criterion: %lld^%lld mod %lld = %lld.  `ok` is %s, so `keep = .not. (ok .and. acc == -1)` is **%s**.\n\n", r.base, r.base, (Q - 1) / 2, Q, r.acc, ok ? "true" : "false", (!(ok && r.acc == -1)) ? "true (kept)" : "false (rejected)");
  return r;
}

// ---- chain_general: polynomials are stored with the leading coefficient first
struct GenResult { bool done, resz; i64 mm, acc; };
static GenResult chain_general(int n, const std::vector<i64> &a, bool trace)
{
  std::vector<i64> U(n + 3, 0), V(n + 3, 0), W(n + 3, 0);
  for (int i = 0; i <= n; i++) U[i] = rd(a[n - i]);
  for (int i = 0; i <= n - 1; i++) V[i] = rd((n - i) * a[n - i]);
  int du = n, dv = n - 1;
  i64 mm = 1; bool done = false, resz = false;
  if (trace) printf("**Start.**  U = f = %s and V = f' = %s, written with the leading coefficient first; `mm = 1`.\n\n", list_desc(U, du + 1).c_str(), list_desc(V, dv + 1).c_str());
  for (int it = 1; it <= 40 && !done; it++) {
    i64 lu = U[0], lv = V[0];
    int imax = du;
    if (trace) printf("#### Iteration %d\n\nU = %s (degree %d), V = %s (degree %d), `mm = %lld`.  `lu = %lld`, `lv = %lld`.\n\n", it, list_desc(U, du + 1).c_str(), du, list_desc(V, dv + 1).c_str(), dv, mm, lu, lv);
    if (trace) printf("W[i] = lv×U[i+1] - lu×V[i+1] for i = 0..%d (the leading terms lv×lu - lu×lv cancelled, so W starts one place later):\n\n| i | lv×U[i+1] | lu×V[i+1] | difference | W[i] mod %lld |\n|---|---|---|---|---|\n", imax - 1, Q);
    for (int i = 0; i < imax; i++) {
      W[i] = rd(lv * U[i + 1] - lu * V[i + 1]);
      if (trace) printf("| %d | %lld×%lld = %lld | %lld×%lld = %lld | %lld | %lld |\n", i, lv, U[i + 1], lv * U[i + 1], lu, V[i + 1], lu * V[i + 1], lv * U[i + 1] - lu * V[i + 1], W[i]);
    }
    for (int i = imax; i < n + 3; i++) W[i] = 0;
    int dw = du - 1;
    int shifts = 0;
    std::string before = list_desc(W, imax);
    while (W[0] == 0 && dw >= 0) { for (int i = 0; i < imax; i++) W[i] = W[i + 1]; W[imax] = 0; dw--; shifts++; }
    if (trace) {
      printf("\n * W = %s, nominally of degree %d.\n", before.c_str(), du - 1);
      printf(" * Leading zeros shifted out: %d, so W = %s of degree dw = %d.\n", shifts, dw >= 0 ? list_desc(W, dw + 1).c_str() : "[]", dw);
    }
    if (dw < 0) { resz = true; done = true; if (trace) printf(" * W is the zero polynomial: the resultant is 0 mod q and the polynomial is kept.\n\n"); break; }
    int ev = std::abs(dw - du + dv);
    bool odd = (ev % 2) != 0;
    if (odd) mm = rd(mm * lv);
    if (trace) printf(" * Exponent of lv in this step: ev = |dw - du + dv| = |%d - %d + %d| = %d, %s%s.\n", dw, du, dv, ev, odd ? "odd" : "even", odd ? (", so `mm := mm×lv = " + std::to_string(mm) + "`").c_str() : ", so `mm` is unchanged");
    bool sw = dw < dv;
    if (!sw) { U = W; du = dw; if (trace) printf(" * dw >= dv, so U := W (degree %d); V stays (degree %d).\n", du, dv); }
    else { U = V; V = W; du = dv; dv = dw; if (trace) printf(" * dw < dv, so U and V swap roles: U := old V (degree %d), V := W (degree %d).\n", du, dv); }
    if (dv == 0) {
      if (du % 2) { mm = rd(mm * V[0]); if (trace) printf(" * V is now the constant %lld and deg U = %d is odd, so `mm := mm×V = %lld`; the chain is finished.\n", V[0], du, mm); }
      else if (trace) printf(" * V is now the constant %lld and deg U = %d is even, so `mm` is unchanged; the chain is finished.\n", V[0], du);
      done = true;
    }
    if (trace) printf("\n");
  }
  GenResult r; r.done = done; r.resz = resz; r.mm = mm; r.acc = powm(mm, (Q - 1) / 2);
  if (trace) printf("**End.**  `mm = %lld`.  Euler's criterion: %lld^%lld mod %lld = %lld, so `keep` is **%s**.\n\n", mm, mm, (Q - 1) / 2, Q, r.acc, (resz || r.acc != -1) ? "true (kept)" : "false (rejected)");
  return r;
}

static GEN toPoly(const std::vector<i64> &a, int n)
{
  GEN f = cgetg(n + 3, t_POL); f[1] = evalsigne(1) | evalvarn(0);
  for (int k = 0; k <= n; k++) gel(f, 2 + k) = stoi((long)a[k]);
  return f;
}
static bool is_prime(i64 q) { if (q < 3 || q % 2 == 0) return false; for (i64 d = 3; d * d <= q; d += 2) if (q % d == 0) return false; return true; }

// symbol of the chain result: acc is 1, -1 or 0 after reduction
static i64 sym(i64 acc) { return acc == 1 ? 1 : acc == 0 ? 0 : -1; }

static int run_test()
{
  setrand(utoi(2024));
  const i64 primes[] = {41, 73, 89, 97, 113, 761, 769, 33554249, 13, 17, 29, 37, 53};
  const int degrees[] = {4, 10};
  i64 allWrong = 0;
  for (int di = 0; di < 2; di++) {
    int n = degrees[di];
    for (size_t pi = 0; pi < sizeof(primes) / sizeof(primes[0]); pi++) {
      Q = primes[pi];
      i64 tested = 0, blockOk = 0, blockWrong = 0, abnormal = 0, genWrong = 0;
      for (int rep = 0; rep < 4000; rep++) {
        pari_sp av = avma;
        std::vector<i64> a(n + 1, 0); a[n] = 1;
        int kind = rep % 4;
        for (int k = 0; k < n; k++) a[k] = (i64)random_Fl(41) - 20;
        if (kind >= 1) a[n - 1] = 0;                                 // depressed
        if (kind == 2) a[n - 2] = 0;                                 // the x^(n-2) coefficient zeroed as well (the situation of item 4)
        if (kind == 3) for (int k = 1; k < n; k += 2) a[k] = 0;      // even polynomial
        GEN D = poldisc0(toPoly(a, n), -1);
        if (gequal0(D)) { avma = av; continue; }
        i64 expect = (i64)kronecker(D, utoi((ulong)Q));
        BlockResult b = chain_block(n, a, false, false);
        GenResult g = chain_general(n, a, false);
        tested++;
        if (b.ok) { blockOk++; if (sym(b.acc) != expect) blockWrong++; } else abnormal++;
        bool genGood = g.resz ? (expect == 0) : (sym(g.acc) == expect);
        if (!genGood) genWrong++;
        avma = av;
      }
      allWrong += blockWrong + genWrong;
      printf("degree %2d, q = %-9lld: %4lld polynomials; chain_block valid for %4lld (wrong %lld), not valid for %4lld; chain_general wrong %lld\n", n, Q, tested, blockOk, blockWrong, abnormal, genWrong);
    }
  }
  printf(allWrong ? "MISMATCH: %lld wrong results\n" : "all results agree with the exact discriminant (%lld wrong)\n", allWrong);
  return allWrong ? 1 : 0;
}

int main(int argc, char **argv)
{
  pari_init(20000000, 100);
  std::string mode = argc >= 2 ? argv[1] : "";
  if (mode == "test") return run_test();
  if (mode == "show" && argc >= 5) {
    Q = atoll(argv[2]);
    std::string which = "both";
    std::vector<i64> lead;
    for (int i = 3; i < argc; i++) {
      std::string s = argv[i];
      if (s == "block" || s == "general" || s == "blockfull") which = s; else lead.push_back(atoll(argv[i]));
    }
    int n = (int)lead.size() - 1;
    if (!is_prime(Q) || n < 2 || lead[0] != 1) { fprintf(stderr, "q must be an odd prime and the polynomial monic of degree at least 2 (the first coefficient is 1)\n"); return 2; }
    std::vector<i64> a(lead.size()); for (size_t i = 0; i < lead.size(); i++) a[lead.size() - 1 - i] = lead[i];
    pari_sp av = avma;
    GEN D = poldisc0(toPoly(a, n), -1);
    i64 kr = (i64)kronecker(D, utoi((ulong)Q));
    printf("<!-- f = %s, q = %lld, disc(f) = %s, disc mod q = %lld, Kronecker = %lld -->\n\n", poly_asc(a, n).c_str(), Q, GENtostr(D), (i64)smodis(D, Q), kr);
    if (which == "block" || which == "blockfull" || which == "both") chain_block(n, a, true, which != "blockfull");
    if (which == "general" || which == "both") chain_general(n, a, true);
    avma = av;
    return 0;
  }
  fprintf(stderr, "usage:  chain_trace show q c_n ... c_0 [block | blockfull | general]\n        chain_trace test\n");
  return 2;
}
