# Worked example: the filter's two Euclidean chains

For a student who wants to see what `chain_block` (`src/pdfilter.f90`) and `chain_general` (`src/chain_general.inc`) actually do.  Three small polynomials are
followed step by step: one that the fast chain handles, and two that it cannot, which are sent to the general chain.  Line numbers are those of `src/pdfilter.f90`
and `src/chain_general.inc` (`inc` below).  Every number in the traces was produced by `chain_trace.cpp` in this directory, a program that follows the Fortran step
by step; "How the numbers were checked" at the end says how to run it and how it was tested.

## 1. What the chain is for

The application wants polynomials f whose discriminant, with every factor of 2 and 5 removed, is a perfect square.  The filter rejects most other polynomials
cheaply.  For a prime q with q = 1 (mod 4) and 2 and 5 squares mod q, a polynomial that can pass has disc(f) either 0 or a square mod q.  So if disc(f) is a
**non-square** mod q, f can never pass and is rejected (header, lines 10-12).

What the filter needs is therefore only the quadratic character of disc(f) mod q, the Legendre symbol: +1 for a non-zero square, -1 for a non-square, 0 for zero.
It gets it from Euclid's algorithm applied to f and its derivative f' over the field F_q, because disc(f) is (up to a sign and a leading coefficient) the
resultant of f and f', and Euclid's algorithm computes resultants.

Euclid's algorithm normally divides by leading coefficients, which needs a modular inverse.  The filter avoids that: before each subtraction it multiplies by the
leading coefficient (pseudo-division).  This changes the answer by known factors, but the Legendre symbol does not see perfect squares (the symbol of c^2 is 1), so
only an **odd** power of a leading coefficient would matter.  That is the whole difference between the two chains:

 - **`chain_block`** (fast): for the normal case in which each remainder has degree exactly one lower than the one before.  It then multiplies by squares only, so
   nothing needs to be remembered (header, lines 17-20) and 64 polynomials can be processed by the same instructions.
 - **`chain_general`** (slow): for any sequence of degrees.  It remembers the odd powers in a running product `mm` (the "Bookkeeping" of `chain_general.inc`, lines 14 and 78).

## 2. The set-up of the examples

 - **Degree 4 instead of 10**, to keep the arithmetic readable.  The Fortran is written for degree 10 (the loop runs m = 9 down to 1, arrays have 11 places); for
   degree 4 the loop runs m = 3, 2, 1 and every step is the same.
 - **The prime q = 41**: the smallest prime for which the filter's conditions hold (41 = 1 mod 4, and 2 and 5 are squares mod 41).  The real primes are 761 to 1009
   (single precision) and about 2^25 (double precision).
 - **Numbers are reduced to the range -20 ... 20**, as `red()` does (lines 162-165): 32 becomes -9, 24 becomes -17.
 - **The squares mod 41** (r has symbol +1 exactly when it is in this list): 1 2 4 5 8 9 10 16 18 20 21 23 25 31 32 33 36 37 39 40.  The code does not look them
   up: it computes r^((q-1)/2) mod q by repeated squaring (lines 308-316), which is +1 for a square and -1 for a non-square.
 - **Storage.**  `chain_block` keeps a polynomial as the coefficients of x^0, x^1, x^2, ... (`p(:, k, slot)` is the coefficient of x^k).  `chain_general` keeps the
   leading coefficient first (index 0).  The traces below use the same conventions as the code.

## 3. The fast chain, `chain_block`

A is f and B is its derivative f'.  Each step uses A of degree m+1 and B of degree m, and produces a new B of degree m-1.  In pseudo-code, with the lines of the
Fortran on the right (n = 10 in the code, n = 4 here):

```
A := f ;  B := f'                      B[k-1] = k*f[k]                       lines 278-282
ok := true                                                                   line 283
for m = n-1 down to 1:                                                       line 285
    lc  := B[m]                        leading coefficient of B              line 287
    top := A[m+1]                      leading coefficient of A              line 288
    R   := lc*A - top*x*B              round 1: removes the x^(m+1) term     lines 289-295
    t   := R[m]
    R2  := lc*R - t*B                  round 2: removes the x^m term         lines 296-302
    if m >= 2:  ok := ok and (R2[m-1] != 0)                                  line 303
    A := B ;  B := R2                                                        lines 304-305
base := B[0]                           the final constant                    line 309
acc  := base^((q-1)/2) mod q           Euler's criterion: +1 or -1           lines 310-316
keep := not (ok and acc == -1)                                               line 318
```

Two rounds per step because one round only cancels the leading term of A, and the remainder can still have degree m; the second round cancels that term too.  Both
rounds multiply by `lc`, so a step multiplies the remainder by lc^2, a square, which the Legendre symbol does not see: that is why no bookkeeping is needed.

The line `ok := ok and (R2[m-1] != 0)` is the safety check.  The next step will use R2[m-1], the leading coefficient of the new B, as its `lc`.  If it is zero the
degree of B has dropped by more than one, the pattern "each remainder is one degree lower" is broken, and the squares argument no longer holds.  The loop cannot
branch (64 polynomials share it), so it only records `ok = false`; line 318 then never rejects such a polynomial, and `pdf_filter` sends it to `chain_general`
(see section 5).

## 4. Example 1: f1 = x^4 + 3x^2 + 2x + 8, which the fast chain solves

Exact discriminant: disc(f1) = 80672, which is 25 = 5^2 mod 41, a square, so the symbol is +1 and f1 must survive the prime 41.  The chain should end with a square.

### The trace

**Start.** A = f = x^4 + 3x^2 + 2x + 8 and B = A' = 4x^3 + 6x + 2.

#### Step m = 3

A = x^4 + 3x^2 + 2x + 8 (degree 4), B = 4x^3 + 6x + 2 (degree 3).  `lc = B[3] = 4`, `top = A[4] = 1`.

**Round 1.**  R = lc×A - top×x×B, that is R[0] = lc×A[0] and R[k] = lc×A[k] - top×B[k-1] for k = 1..3:

| k | lc×A[k] | top×B[k-1] | difference | R[k] mod 41 |
|---|---|---|---|---|
| 0 | 4×8 = 32 | - | 32 | -9 |
| 1 | 4×2 = 8 | 1×2 = 2 | 6 | 6 |
| 2 | 4×3 = 12 | 1×6 = 6 | 6 | 6 |
| 3 | 4×0 = 0 | 1×0 = 0 | 0 | 0 |

So R = 6x^2 + 6x - 9.

**Round 2.**  t = R[3] = 0;  R2 = lc×R - t×B, that is R2[k] = lc×R[k] - t×B[k] for k = 0..2:

| k | lc×R[k] | t×B[k] | difference | R2[k] mod 41 |
|---|---|---|---|---|
| 0 | 4×-9 = -36 | 0×2 = 0 | -36 | 5 |
| 1 | 4×6 = 24 | 0×6 = 0 | 24 | -17 |
| 2 | 4×6 = 24 | 0×0 = 0 | 24 | -17 |

So R2 = -17x^2 - 17x + 5, of degree at most 2.

**Check** (pdfilter.f90 line 303): the new leading coefficient R2[2] = -17 is not zero, so `ok` stays true.

Next step: A := B, B := R2.

#### Step m = 2

A = 4x^3 + 6x + 2 (degree 3), B = -17x^2 - 17x + 5 (degree 2).  `lc = B[2] = -17`, `top = A[3] = 4`.

**Round 1.**  R = lc×A - top×x×B, that is R[0] = lc×A[0] and R[k] = lc×A[k] - top×B[k-1] for k = 1..2:

| k | lc×A[k] | top×B[k-1] | difference | R[k] mod 41 |
|---|---|---|---|---|
| 0 | -17×2 = -34 | - | -34 | 7 |
| 1 | -17×6 = -102 | 4×5 = 20 | -122 | 1 |
| 2 | -17×0 = 0 | 4×-17 = -68 | 68 | -14 |

So R = -14x^2 + x + 7.

**Round 2.**  t = R[2] = -14;  R2 = lc×R - t×B, that is R2[k] = lc×R[k] - t×B[k] for k = 0..1:

| k | lc×R[k] | t×B[k] | difference | R2[k] mod 41 |
|---|---|---|---|---|
| 0 | -17×7 = -119 | -14×5 = -70 | -49 | -8 |
| 1 | -17×1 = -17 | -14×-17 = 238 | -255 | -9 |

So R2 = -9x - 8, of degree at most 1.

**Check** (pdfilter.f90 line 303): the new leading coefficient R2[1] = -9 is not zero, so `ok` stays true.

Next step: A := B, B := R2.

#### Step m = 1

A = -17x^2 - 17x + 5 (degree 2), B = -9x - 8 (degree 1).  `lc = B[1] = -9`, `top = A[2] = -17`.

**Round 1.**  R = lc×A - top×x×B, that is R[0] = lc×A[0] and R[k] = lc×A[k] - top×B[k-1] for k = 1..1:

| k | lc×A[k] | top×B[k-1] | difference | R[k] mod 41 |
|---|---|---|---|---|
| 0 | -9×5 = -45 | - | -45 | -4 |
| 1 | -9×-17 = 153 | -17×-8 = 136 | 17 | 17 |

So R = 17x - 4.

**Round 2.**  t = R[1] = 17;  R2 = lc×R - t×B, that is R2[k] = lc×R[k] - t×B[k] for k = 0..0:

| k | lc×R[k] | t×B[k] | difference | R2[k] mod 41 |
|---|---|---|---|---|
| 0 | -9×-4 = 36 | 17×-8 = -136 | 172 | 8 |

So R2 = 8, of degree at most 0.

**End.**  The final constant is 8.  Euler's criterion: 8^20 mod 41 = 1.  `ok` is true, so `keep = .not. (ok .and. acc == -1)` is **true (kept)**.


**What to notice.**

 - The degrees go 4, 3, 2, 1, 0: every B is one degree lower than the one before, which is what `chain_block` assumes.
 - The check passed twice (at m = 3 with R2[2] = -17, at m = 2 with R2[1] = -9), so `ok` is still true at the end.
 - The final constant 8 is not disc(f1) mod 41 (that is 25).  They differ by a square (the `lc^2` factors), so they have the same symbol: 8 = 7^2 mod 41, and
   Euler's criterion gives +1, as the discriminant said.
 - `keep` is true: f1 survives this prime.  The filter would go on with the next prime.

## 5. Example 2: f2 = x^4 + 2x + 8, which sends the filter to `chain_general`

Take f1 and **remove the x^2 term**.  This is what item 4 of `CHANGES.md` does to the polynomials of the test (there the x^8 coefficient of a degree-10
polynomial; the first check below is the degree-4 version of it).  Exact discriminant: disc(f2) = 130640, which is 14 mod 41, not in the list of squares, so
the symbol is -1 and f2 must be **rejected** at q = 41.

### Pass 1: `chain_block`

**Start.** A = f = x^4 + 2x + 8 and B = A' = 4x^3 + 2.

#### Step m = 3

A = x^4 + 2x + 8 (degree 4), B = 4x^3 + 2 (degree 3).  `lc = B[3] = 4`, `top = A[4] = 1`.

**Round 1.**  R = lc×A - top×x×B, that is R[0] = lc×A[0] and R[k] = lc×A[k] - top×B[k-1] for k = 1..3:

| k | lc×A[k] | top×B[k-1] | difference | R[k] mod 41 |
|---|---|---|---|---|
| 0 | 4×8 = 32 | - | 32 | -9 |
| 1 | 4×2 = 8 | 1×2 = 2 | 6 | 6 |
| 2 | 4×0 = 0 | 1×0 = 0 | 0 | 0 |
| 3 | 4×0 = 0 | 1×0 = 0 | 0 | 0 |

So R = 6x - 9.

**Round 2.**  t = R[3] = 0;  R2 = lc×R - t×B, that is R2[k] = lc×R[k] - t×B[k] for k = 0..2:

| k | lc×R[k] | t×B[k] | difference | R2[k] mod 41 |
|---|---|---|---|---|
| 0 | 4×-9 = -36 | 0×2 = 0 | -36 | 5 |
| 1 | 4×6 = 24 | 0×0 = 0 | 24 | -17 |
| 2 | 4×0 = 0 | 0×0 = 0 | 0 | 0 |

So R2 = -17x + 5, of degree at most 2.

**Check** (pdfilter.f90 line 303): the new leading coefficient R2[2] = 0 is ZERO, so `ok` is now FALSE.

(The loop would carry on with `lc = 0` and meaningless numbers, but `ok` stays false; the result is not used.)


**Why the check failed.**  At the first step the value tested, R2[2], works out to 8×a2 - 3×a3^2 (a2, a3 the coefficients of x^2 and x^3 of f; for degree 10 the
same step gives 20×a8 - 9×a9^2).  For f1 it was 8×3 = 24, not zero.  For f2 both a2 and a3 are 0, so it is 0: the new remainder lost its leading terms and has
degree 1 instead of 2, so the degrees go 4, 3, 1 instead of 4, 3, 2.  `chain_block` stops being valid for this polynomial at its very first step, whatever the prime.

**What `pdf_filter` does with it** (prime stage q = 41, lines 426-456):

```
pass 1 (lines 431-444):  run_block(all polynomials, general = .false.)  ->  chain_block        line 433
          f1:  ok = true   and  keep = true       ->  kept, goes on to the next prime        lines 439-441
          f2:  ok = false                         ->  put in the list uidx                   lines 436-438
pass 2 (lines 446-456):  run_block(uidx, general = .true.)               ->  chain_general    line 449
          f2:  decided by chain_general (below)
```

### Pass 2: `chain_general`

It does the same job with a different method.  Differences from `chain_block`, visible in the trace:

 - A polynomial is stored with its leading coefficient first, and U (degree du) and V (degree dv) are reduced **one elementary step at a time**:
   W = lv×U - lu×V, where the leading terms cancel, so W starts one place later (inc lines 57-62).  There is no pairing of two rounds.
 - If W then starts with zeros, they are **shifted out** and the degree falls by one for each (inc lines 63-72).  This is what lets it follow a remainder whose
   degree drops by more than one.  If W becomes the zero polynomial the resultant is 0 and the polynomial is kept (inc lines 73-77).
 - It computes `ev = |dw - du + dv|`, the exponent of `lv` that this step contributed.  If `ev` is odd, `mm` is multiplied by `lv` (inc lines 78-81).  Even exponents
   are squares and are ignored.
 - If the new W has lower degree than V, U and V swap roles (inc lines 82-106); this only changes a sign, and -1 is a square mod 41.
 - When V is a constant c the chain is finished: the resultant is c^(deg U), so `mm` is multiplied by c if deg U is odd (inc lines 107-110).
 - Finally the Euler criterion is applied to `mm` (inc lines 113-121), and the polynomial is rejected only if the chain finished normally and the symbol is -1 (inc line 122).

**Start.**  U = f = [1, 0, 0, 2, 8] and V = f' = [4, 0, 0, 2], written with the leading coefficient first; `mm = 1`.

#### Iteration 1

U = [1, 0, 0, 2, 8] (degree 4), V = [4, 0, 0, 2] (degree 3), `mm = 1`.  `lu = 1`, `lv = 4`.

W[i] = lv×U[i+1] - lu×V[i+1] for i = 0..3 (the leading terms lv×lu - lu×lv cancelled, so W starts one place later):

| i | lv×U[i+1] | lu×V[i+1] | difference | W[i] mod 41 |
|---|---|---|---|---|
| 0 | 4×0 = 0 | 1×0 = 0 | 0 | 0 |
| 1 | 4×0 = 0 | 1×0 = 0 | 0 | 0 |
| 2 | 4×2 = 8 | 1×2 = 2 | 6 | 6 |
| 3 | 4×8 = 32 | 1×0 = 0 | 32 | -9 |

 * W = [0, 0, 6, -9], nominally of degree 3.
 * Leading zeros shifted out: 2, so W = [6, -9] of degree dw = 1.
 * Exponent of lv in this step: ev = |dw - du + dv| = |1 - 4 + 3| = 0, even, so `mm` is unchanged.
 * dw < dv, so U and V swap roles: U := old V (degree 3), V := W (degree 1).

#### Iteration 2

U = [4, 0, 0, 2] (degree 3), V = [6, -9] (degree 1), `mm = 1`.  `lu = 4`, `lv = 6`.

W[i] = lv×U[i+1] - lu×V[i+1] for i = 0..2 (the leading terms lv×lu - lu×lv cancelled, so W starts one place later):

| i | lv×U[i+1] | lu×V[i+1] | difference | W[i] mod 41 |
|---|---|---|---|---|
| 0 | 6×0 = 0 | 4×-9 = -36 | 36 | -5 |
| 1 | 6×0 = 0 | 4×0 = 0 | 0 | 0 |
| 2 | 6×2 = 12 | 4×0 = 0 | 12 | 12 |

 * W = [-5, 0, 12], nominally of degree 2.
 * Leading zeros shifted out: 0, so W = [-5, 0, 12] of degree dw = 2.
 * Exponent of lv in this step: ev = |dw - du + dv| = |2 - 3 + 1| = 0, even, so `mm` is unchanged.
 * dw >= dv, so U := W (degree 2); V stays (degree 1).

#### Iteration 3

U = [-5, 0, 12] (degree 2), V = [6, -9] (degree 1), `mm = 1`.  `lu = -5`, `lv = 6`.

W[i] = lv×U[i+1] - lu×V[i+1] for i = 0..1 (the leading terms lv×lu - lu×lv cancelled, so W starts one place later):

| i | lv×U[i+1] | lu×V[i+1] | difference | W[i] mod 41 |
|---|---|---|---|---|
| 0 | 6×0 = 0 | -5×-9 = 45 | -45 | -4 |
| 1 | 6×12 = 72 | -5×0 = 0 | 72 | -10 |

 * W = [-4, -10], nominally of degree 1.
 * Leading zeros shifted out: 0, so W = [-4, -10] of degree dw = 1.
 * Exponent of lv in this step: ev = |dw - du + dv| = |1 - 2 + 1| = 0, even, so `mm` is unchanged.
 * dw >= dv, so U := W (degree 1); V stays (degree 1).

#### Iteration 4

U = [-4, -10] (degree 1), V = [6, -9] (degree 1), `mm = 1`.  `lu = -4`, `lv = 6`.

W[i] = lv×U[i+1] - lu×V[i+1] for i = 0..0 (the leading terms lv×lu - lu×lv cancelled, so W starts one place later):

| i | lv×U[i+1] | lu×V[i+1] | difference | W[i] mod 41 |
|---|---|---|---|---|
| 0 | 6×-10 = -60 | -4×-9 = 36 | -96 | -14 |

 * W = [-14], nominally of degree 0.
 * Leading zeros shifted out: 0, so W = [-14] of degree dw = 0.
 * Exponent of lv in this step: ev = |dw - du + dv| = |0 - 1 + 1| = 0, even, so `mm` is unchanged.
 * dw < dv, so U and V swap roles: U := old V (degree 1), V := W (degree 0).
 * V is now the constant -14 and deg U = 1 is odd, so `mm := mm×V = -14`; the chain is finished.

**End.**  `mm = -14`.  Euler's criterion: -14^20 mod 41 = -1, so `keep` is **false (rejected)**.


**What to notice.**

 - Iteration 1 produced W = [0, 0, 6, -9]: two leading zeros, shifted out, so the degree dropped from 3 straight to 1 (the sequence 4, 3, 1, 0).  This is the step
   `chain_block` cannot make.
 - The exponent `ev` was even in every iteration, so `mm` stayed 1 until the very end, when V became the constant -14 and deg U = 1 was odd.
 - The result `mm = -14` is, mod 41, **minus the discriminant** (disc = 130640 = 14 mod 41): for x^4 + a1 x + a0 the chain computes 27 a1^4 - 256 a0^3 = -disc exactly.
   Its symbol is the symbol of disc because -1 is a square mod 41.
 - 14 is not a square mod 41, so the Euler criterion gives -1 and f2 is **rejected**, as the discriminant predicted.  The general chain is slower, but it decides.

## 6. Example 3: f3 = x^4 + 8, a shorter case with an odd exponent

This is the shape the test actually generates: an even polynomial g(x^2) with the x^2 coefficient zero.  disc(f3) = 131072 = 36 mod 41, a square, so f3 must be kept.

### Pass 1: `chain_block`

**Start.** A = f = x^4 + 8 and B = A' = 4x^3.

#### Step m = 3

A = x^4 + 8 (degree 4), B = 4x^3 (degree 3).  `lc = B[3] = 4`, `top = A[4] = 1`.

**Round 1.**  R = lc×A - top×x×B, that is R[0] = lc×A[0] and R[k] = lc×A[k] - top×B[k-1] for k = 1..3:

| k | lc×A[k] | top×B[k-1] | difference | R[k] mod 41 |
|---|---|---|---|---|
| 0 | 4×8 = 32 | - | 32 | -9 |
| 1 | 4×0 = 0 | 1×0 = 0 | 0 | 0 |
| 2 | 4×0 = 0 | 1×0 = 0 | 0 | 0 |
| 3 | 4×0 = 0 | 1×0 = 0 | 0 | 0 |

So R = -9.

**Round 2.**  t = R[3] = 0;  R2 = lc×R - t×B, that is R2[k] = lc×R[k] - t×B[k] for k = 0..2:

| k | lc×R[k] | t×B[k] | difference | R2[k] mod 41 |
|---|---|---|---|---|
| 0 | 4×-9 = -36 | 0×0 = 0 | -36 | 5 |
| 1 | 4×0 = 0 | 0×0 = 0 | 0 | 0 |
| 2 | 4×0 = 0 | 0×0 = 0 | 0 | 0 |

So R2 = 5, of degree at most 2.

**Check** (pdfilter.f90 line 303): the new leading coefficient R2[2] = 0 is ZERO, so `ok` is now FALSE.

(The loop would carry on with `lc = 0` and meaningless numbers, but `ok` stays false; the result is not used.)


The check fails at the first step again (a2 = a3 = 0), so f3 goes to pass 2.

### Pass 2: `chain_general`

**Start.**  U = f = [1, 0, 0, 0, 8] and V = f' = [4, 0, 0, 0], written with the leading coefficient first; `mm = 1`.

#### Iteration 1

U = [1, 0, 0, 0, 8] (degree 4), V = [4, 0, 0, 0] (degree 3), `mm = 1`.  `lu = 1`, `lv = 4`.

W[i] = lv×U[i+1] - lu×V[i+1] for i = 0..3 (the leading terms lv×lu - lu×lv cancelled, so W starts one place later):

| i | lv×U[i+1] | lu×V[i+1] | difference | W[i] mod 41 |
|---|---|---|---|---|
| 0 | 4×0 = 0 | 1×0 = 0 | 0 | 0 |
| 1 | 4×0 = 0 | 1×0 = 0 | 0 | 0 |
| 2 | 4×0 = 0 | 1×0 = 0 | 0 | 0 |
| 3 | 4×8 = 32 | 1×0 = 0 | 32 | -9 |

 * W = [0, 0, 0, -9], nominally of degree 3.
 * Leading zeros shifted out: 3, so W = [-9] of degree dw = 0.
 * Exponent of lv in this step: ev = |dw - du + dv| = |0 - 4 + 3| = 1, odd, so `mm := mm×lv = 4`.
 * dw < dv, so U and V swap roles: U := old V (degree 3), V := W (degree 0).
 * V is now the constant -9 and deg U = 3 is odd, so `mm := mm×V = 5`; the chain is finished.

**End.**  `mm = 5`.  Euler's criterion: 5^20 mod 41 = 1, so `keep` is **true (kept)**.


**What to notice.**

 - One iteration is enough: W = [0, 0, 0, -9] has **three** leading zeros, which are all shifted out, so the degree falls from 4 to 0 at once (the sequence 4, 3, 0).
 - Here `ev = |0 - 4 + 3| = 1` is **odd**: this is the case the bookkeeping exists for.  `mm` becomes lv = 4, and the final constant -9 with deg U = 3 odd makes it 4×(-9) = -36 = 5.
 - 5 is a square mod 41, so the symbol is +1 and f3 is kept, as its discriminant (36, a square) says.  Note that `mm` = 5 differs from disc = 36: equal symbols, not equal numbers.

## 7. How this relates to the degree-10 code and to the test

 - **The same first check.**  In degree 10 the first step (m = 9) tests R2[8] = 20×a8 - 9×a9^2.  The polynomials of the test have a9 = 0 (before the shift), so it is
   20×a8: not zero for a random x^8 coefficient from 1 to 99 (all the filter primes are at least 761), zero when the x^8 coefficient is 0.  The shift x -> x+t does not
   change it, because a9 becomes 10t and a8 becomes a8 + 45t^2 and 20(a8 + 45t^2) - 9(10t)^2 = 20×a8.  This is why a quarter of the polynomials get a zero x^8
   coefficient (item 4 of `CHANGES.md`): it sends the polynomial to `chain_general` at every prime.
 - **Seen in the real filter.**  With `NF_PDFDEBUG=1` the first call prints one line per prime (line 458): lanes in, lanes re-run in `chain_general` ("retried"), lanes
   kept.  For the degree-10 polynomial x^10 + c x^8 + 5x^6 + 7x^4 + 11x^2 + 4: with c = 3 every one of the 30 lines says "retried 0" (all in `chain_block`); with c = 0 every
   line says "retried 1" (all sent to `chain_general`).
 - **Not shown here.**  The real code runs 64 polynomials at once, so each of these steps is a loop over 64 lanes; and after one call at a prime where more than half
   of the lanes left `chain_block`, `prefgen` (lines 63, 429, 445) sends all lanes of the next calls straight to `chain_general` for that prime (except every 64th call).

## How the numbers were checked

The traces come from `chain_trace.cpp`, which follows the two chains step by step for any degree and prime.  Build it with the usual scripts and run it
(the first line of its output is a comment giving the exact discriminant, which this document leaves out):

```
TEST=chain_trace.cpp bash build-linux.sh          # or build-windows.sh; makes build/chain_trace
build/chain_trace show 41 1 0 3 2 8 block         # Example 1, pass 1 (the whole trace)
build/chain_trace show 41 1 0 0 2 8 block         # Example 2, pass 1 (stops at the failed check)
build/chain_trace show 41 1 0 0 2 8 general       # Example 2, pass 2
build/chain_trace show 41 1 0 0 0 8 block         # Example 3, pass 1
build/chain_trace show 41 1 0 0 0 8 general       # Example 3, pass 2
build/chain_trace test                            # the comparison with the exact discriminant described next
```

The arguments after `show` are the prime and then the coefficients from the leading one down, so `1 0 3 2 8` is x^4 + 3x^2 + 2x + 8.  Its output was checked on
Windows (PARI 2.8) and on Linux (PARI 2.17.3): the five traces above are identical to the ones in this document.  The program was compared with the exact
discriminant (PARI), which is what `chain_trace test` repeats: for about 4,000 random polynomials of degree 4 and of degree 10 at each of 13 primes, both chains
gave the symbol of the exact discriminant every time, and the polynomials included zeroed x^(n-1) and x^(n-2) coefficients and even polynomials.  It was also
compared with the real Fortran filter, with a throw-away harness that is not in the repository: for 250 random degree-10 polynomials the real filter's "retried"
count at q = 1009 and its reject decision agreed with the program in every case.  To check an example yourself in GP:

```
? f = x^4 + 3*x^2 + 2*x + 8;  poldisc(f)             \\ 80672
? kronecker(poldisc(f), 41)                          \\ 1   (disc is a square mod 41, so f survives)
? f = x^4 + 2*x + 8;  kronecker(poldisc(f), 41)      \\ -1  (a non-square: f is rejected)
? f = x^4 + 8;  kronecker(poldisc(f), 41)            \\ 1
```
