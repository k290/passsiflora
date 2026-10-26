# Changes to Eric's test, and why

> Written for the floating-point filter (commit 327b04f): its line numbers and the fault counts below are of that version.  The integer
> filter has the same algorithms and routing at different lines (`pdfilter_chains.tex`, Appendix A); its fault table, with the oracle test
> `test_pdfilter_oracle.cpp`, is in `README.md` and in Section 13 of `pdfilter_chains.tex`.

Eric's second `test_pdfilter.cpp` (f = g(x^2), g a random quintic) against `test_pdfilter_linux.cpp`; `git diff pdfilter_Test/test_pdfilter_linux.cpp`
shows every line.  Numbers are valid polynomials wrongly rejected by a deliberately broken filter, per 1,000,000 polynomials, AVX2 / baseline build
(`fault_checks.sh`).

## 1. Command line: optional seed, seed printed

**Change.** `if(argc!=2)` became `if(argc<2 || argc>3)`: `N` alone (2 words, counting the program name) or `N seed` (3) is accepted, anything else prints the usage.
`seed` is `argv[2]`, or `time(NULL)` if absent, and is used for both generators (`setrand` for PARI, `srand` for `rand()`).  It is printed as `seed = ...`.
The usage now says `N*10000` polynomials, since the blocks no longer have a fixed size.

**Why.** Eric's version seeded from the clock and did not print it, so a failing run could not be repeated, and two runs started in the same second got the
same polynomials.  With a seed, a failure is reproduced by `test_pdfilter N seed`, and `fault_checks.sh` gives every broken filter the same polynomials, so
the counts can be compared.  The seed must be a positive integer (PARI refuses 0, which is also what text that is not a number turns into).
The same seed gives the same polynomials only with the same C library (glibc and the MinGW runtime have different `rand()`).

## 2. Varied block sizes

**Change.** Each call of the filter gets one of 1, 7, 63, 64, 65, 127, 255, 1000, 4096, 10000 polynomials (`blockSizes`), until `N*10000` have been tested.

**Why.** The filter works in chunks of 64, so a fixed 10000 always ends with a chunk of 16 polynomials in the same positions, and a fault that only affects
the last chunk is seen at one size only.  Measured (x87 arithmetic): 2 / 1 with fixed blocks, 39 / 44 with varied ones.

## 3. Shift x -> x+t in three quarters of the polynomials

**Change.** After the irreducibility test, `f = ZX_translate(f, t)` with 1 <= |t| <= 30.

**Why.** The discriminant does not change, so the polynomial still must pass, but it becomes unlike Eric's: signed, with large coefficients (the real passing
polynomials have 15 to 30 bits, 73% have a negative coefficient; Eric's had at most 12 bits and 16% a negative one), with the odd powers no longer zero, and with
the zeros that the Euclidean chain meets no longer being zero coefficients but appearing by cancellation, which is where a rounding error can turn a zero into
a non-zero.  Measured (rounding mode up): 0 without, 2,490 with.

## 4. Zero x^8 coefficient in a quarter of the polynomials

**Change.** `if( rand()%4 == 0 )  c[4] = 0;`

**Why.** With a zero x^8 coefficient the Euclidean chain of (f, f') has the degrees 10, 9, 6, ... instead of 10, 9, 8, ..., and the filter handles it in its general
chain (`chain_general.inc`) instead of `chain_block`.  25% of the polynomials that pass the discriminant test in the application are like that; Eric's generator
made 0.8% of them, by chance.  Measured (rounding mode up, with the shift): 2,554 without, 37,697 with.  How the filter tells the two paths apart, with line
numbers: "Background for item 4" below.

## 5. Exit status

**Change.** `return totalFailures ? 1 : 0;`

**Why.** The program always exited 0, even after printing failures, so a script could not tell.  (Exit status 1 is also what the usage error gives; the
`Total failures` line tells them apart.)

## 6. Comment

**Change.** `g[0]^5` in the comment above `genPoly()` is `g[0]`.

**Why.** disc(g(x^2)) = -2^10 * g[0] * disc(g)^2.  The fifth power has the same square class, so the construction was right, but the formula was not.

## Background for item 4: `chain_block` or `chain_general` (src/pdfilter.f90)

A step-by-step numerical example, with a polynomial for each chain, is in [WORKED_EXAMPLE.md](WORKED_EXAMPLE.md).

For each of the 30 primes q the filter runs Euclid's algorithm on (f, f') over F_q.  The last remainder, raised to the power (q-1)/2 (Euler's criterion, lines
308-316), gives the Legendre symbol of disc(f); -1 means the polynomial is rejected.  There are two implementations of the Euclidean chain:

 - **`chain_block`** (single precision 192-239, double 266-319): the fast one, for the normal case where each remainder has degree exactly one less than the
   last (10, 9, 8, ..., 0).  Then every leading-coefficient power that appears has an even exponent, so the symbol of the last constant is the answer and
   no inverses are needed (the comment, lines 17-20), and 64 polynomials ("lanes", `BLK = 64`, line 44) go through the same instructions at once.
 - **`chain_general`** (`chain_general.inc`, included at 346-353 for single and 361-368 for double precision): the slow one, for any degree sequence.  The
   leading-coefficient powers no longer pair up, so it multiplies the odd ones into a running product `mm` (the "Bookkeeping" of the comments at inc 14-19 and
   78, explained below) and takes the Euler criterion of `mm` at the end (inc 113-121).
 - **The check** (line 225 single / 303 double): `if (m >= 2) ok = ok .and. (p(:, m-1, ir) /= 0)`.  After every step, is the leading coefficient of the new
   remainder non-zero mod q?  If not, `ok` is false for that lane: `chain_block`'s answer is not valid, and line 238 / 318,
   `keep = .not. (ok .and. acc == -1)`, makes sure it never rejects it.

Who calls which:

 - `run_block` (372-394) runs one block of lanes through one prime and picks the chain by its argument `general`: small primes at 381-385 (`chain_general_s` at 382,
   `chain_block_s` at 384), large primes at 388-392 (`chain_general` at 389, `chain_block` at 391).
 - `pdf_filter` (398-468), for each prime (loop at 426): **pass 1** (431-444) sends all surviving polynomials through `run_block` with `general = usegen`, normally
   false, so through `chain_block` (433).  A lane with `.not. ok(l)` is put in the list `uidx` (436-438); the others are kept or rejected right there (439-442).
   **Pass 2** (446-456) runs exactly the `uidx` lanes again with `general = .true.`, through `chain_general` (449).
 - `prefgen` (declared 63, set 445, used 429): if more than half of the lanes of the last call left `chain_block` at a prime, the next calls skip pass 1 and go
   straight to `chain_general` for that prime, except on every 64th call (`probe`, 425), so a change in the mix is noticed.

What the "bookkeeping" is (the word is in `chain_general.inc` lines 14 and 78, not in `pdfilter.f90`).  To avoid dividing by a leading coefficient mod q, which would
need a modular inverse, every reduction step first multiplies by lc(V): U' = lc(V)*U - lc(U)*x^k*V (inc 8-12).  That multiplies the resultant by lc(V)^e, with
e = deg U' - deg U + deg V (inc 14-16; computed as `ev` at inc 79).  Only the Legendre symbol of the discriminant is wanted, and chi(c^2) = 1, so a factor with an
even exponent can be ignored; only an odd exponent changes the answer.  `chain_block` does two rounds per step with the same lc (lines 287-302), which scales by lc^2,
so there is nothing to remember.  When the degrees drop irregularly an odd exponent has no partner, and `chain_general` multiplies lc(V) into `mm` (set to 1 at
inc 50; updated at inc 80-81 when `ev` is odd).  When V has become a constant c, Res = c^deg(U), so c is multiplied in as well if deg(U) is odd (inc 108-109).
The symbol of `mm` is then the symbol of disc(f).  Swapping U and V only changes a sign, and -1 is a square because the filter primes are 1 mod 4 (inc 16; line 10).

Why x^8 = 0 sends a lane to `chain_general`: the first step (m = 9) works on A = f and B = f'; `lc` is B's leading coefficient (10) and `top` is A's (1)
(lines 287-288).  The x^8 coefficient of the new remainder, which line 303 tests as `p(:, 8, ir)`, is 20*a8 - 9*a9^2 (a8, a9 the coefficients of x^8, x^9 of f).
The polynomials of `genPoly()` have a9 = 0, so it is 20*a8 mod q: zero exactly when a8 is zero.  A random a8 from 1 to 99 is never zero mod the filter primes
(all above 761), so the lane stays in `chain_block`; with a8 = 0 the test fails at the first step for every prime, the next remainder has degree 6 instead of 8
(10, 9, 6, ...), and the lane goes to `chain_general`.  The shift x -> x+t does not change this: for g(x) = f(x+t), a9' = 10t and a8' = a8 + 45t^2, and
20*a8' - 9*a9'^2 = 20*a8 again.

To see it, set `NF_PDFDEBUG=1` (read at line 123); the first call then prints one line per prime (line 458): lanes in, lanes re-run in `chain_general` ("retried"),
lanes kept.  For the single polynomial x^10 + c*x^8 + 5x^6 + 7x^4 + 11x^2 + 4 (a valid one, checked): with c = 3 every one of the 30 lines says "retried 0";
with c = 0 every line says "retried 1".

## Around the test file (not in the diff above)

 - `test_pdfilter.cpp`: the same changes, plus what Windows needs: `i64` and the filter's declarations before `pari.h` (which redefines `long` there), `iabs`
   instead of `labs`, and `isirreducible` for PARI older than 2.9.
 - `build-linux.sh`, `build-windows.sh`, `pari-config.sh`: build the filter and the test with a PARI source tree or the system PARI, into `build/`.
 - `fault_checks.sh`, `rmshim.c`: produce the numbers quoted here, by breaking the filter on purpose.
 - `README.md`: files, how to run, and the results on Linux and Windows.
 - `WORKED_EXAMPLE.md` and `chain_trace.cpp`: a step-by-step model of `chain_block` and `chain_general` (any degree and prime, printed trace) and the worked
   example for students that was generated with it.
