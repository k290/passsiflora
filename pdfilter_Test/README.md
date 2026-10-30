# pdfilter_Test

Tests of the pre-filter `src/pdfilter.f90` (with `src/chain_general.inc`; the integer version).  Every test script compiles `../src/pdfilter.f90`
itself: there is no copy of the filter in this directory, so a test always tests what ships.

* `test_pdfilter_oracle.cpp` compares the filter with an **exact oracle**: for every polynomial and every filter prime it computes the
  Legendre symbol of the discriminant independently, so it sees both a filter that rejects too much and one that rejects too little, and,
  with the filter's test hook, the stage of every rejection and the routing between `chain_block` and `chain_general`.  Its generators are
  built to reach every stage and every path of the filter (random and real polynomials hardly ever get past the first few stages).
* The older tests (`test_pdfilter.cpp` and Eric's versions) generate only polynomials whose discriminant is (a number built from the
  primes 2 and 5) times a perfect square, so the filter **must keep** them; a failure means that the filter wrongly rejected a valid
  polynomial.  They cannot see a filter that rejects too little.

`pdfilter_chains.tex` (and its PDF) proves the algorithms, the integer arithmetic and what the oracle test's generators do, and gives the results.

## Files

| file | what it is |
|---|---|
| `test_pdfilter_oracle.cpp` | **the test to use**: the oracle test (strata, metamorphic and routing checks, coverage report; `-h` lists the options). |
| `capshim.c` | captures the candidates the application hands to the filter, with its verdicts, for `test_pdfilter_oracle -cand` (instructions in the file). |
| `test_pdfilter.cpp` | Eric's second version plus four small changes (below), rewritten so that it builds on Linux and on Windows ("must pass" polynomials only). |
| `test_pdfilter_linux.cpp` | Eric's second version with the same four changes applied as a small patch (+47 / -15 lines; Linux only, PARI 2.9 or later).  `git diff ccb629f:pdfilter_Test/test_pdfilter_eric.cpp HEAD:pdfilter_Test/test_pdfilter_linux.cpp` shows the diff against Eric's file.  It gives the same numbers as `test_pdfilter.cpp` on Linux. |
| `CHANGES.md` | the changes to Eric's test, each with its reason, and how the filter picks `chain_block` or `chain_general`. |
| `WORKED_EXAMPLE.md` | a worked example for students: small polynomials followed step by step through `chain_block` and `chain_general`. |
| `chain_trace.cpp` | the program that produced those traces (`TEST=chain_trace.cpp bash build-linux.sh`, then `build/chain_trace show ...` or `build/chain_trace test`). |
| Eric's original | not in the tree; `git show ccb629f:pdfilter_Test/test_pdfilter_eric.cpp` (it is the file he sent). |
| `test_pdfilter_full.cpp` | an earlier, larger harness (several generators, a corpus of real polynomials, rounding-mode options). |
| `build-linux.sh`, `build-windows.sh`, `pari-config.sh` | build the filter (AVX2, `-mavx2`, and baseline, with the test hook) and the test; `FILTER=dir` builds another filter (e.g. the floating-point one of commit 327b04f), `REF=dir` links a second filter for the oracle test's differential check. |
| `arith_checks.sh`, `arith_checks.f90`, `limit_probe.c` | the compiled reductions against exact models, and a model of the chains beyond the size limits (`pdfilter_chains.tex`, Section 8.8). |
| `tiny_checks.sh`, `test_tiny.cpp` | the compiled filter on **every** polynomial of complete sets at tiny primes (13: all depressed monic degree-10 polynomials; 17-41: structured sets; shifts), both builds, normal routing and every lane through `chain_general`, against an exact Euclid oracle. |
| `poldisc_checks.sh`, `test_poldisc.cpp` | upstream's compiled `polDiscTest_cpu` against the definition (PARI), with a 128-bit replay of its int64 steps, and the filter against both (`UPSTREAM=<GetDecicsSrc>`). |
| `test_gpu_filter.cpp` | the GPU kernel's integer filter (`src/nf_kernel.h`, included, not copied) compiled for the CPU, prime by prime against an exact oracle; first it runs the GPU mode's start-up self-test with the kernel compiled for the CPU in the role of the device. |
| `checks.py` | exact checks of the mathematics (C0-C10), the source-line citations of the document (C11, `--cites`, `--fix-cites`) and its ledger (C12, `--ledger`). |
| `fault_checks.sh`, `rmshim.c` | break the filter on purpose (floating-point faults, prime ranges, seven faults of the integer arithmetic) and count what each test reports. |
| `repro/` | the fault guide of `pdfilter_chains.tex` (floating-point filter of commit 327b04f; `repro.sh` builds that filter). |
| `test_pdfilter` | Eric's old Linux binary.  It is truncated (4 MB of the ~29 MB its ELF header expects) and segfaults; use the build scripts. |

## Run

Linux (gfortran, g++ and PARI; `PARI=/path/to/compiled/pari-tree`, otherwise Eric's tree, otherwise the system PARI, `apt install libpari-dev`):

    TEST=test_pdfilter_oracle.cpp bash build-linux.sh && build/test_pdfilter_oracle 100 12345   # 1,000,000 polynomials, seed 12345
    build/test_pdfilter_oracle 100 12345 -stages -real hits_vec.txt -cand cand.txt              # + table per stage, real results, candidates
    bash build-linux.sh && build/test_pdfilter 100          # 100 blocks of 10000 polynomials; a second argument is the seed
    TEST=test_pdfilter_linux.cpp bash build-linux.sh        # builds the small patch of Eric's file instead (build/test_pdfilter_linux)
    git show ccb629f:pdfilter_Test/test_pdfilter_eric.cpp > /tmp/eric_original.cpp && TEST=/tmp/eric_original.cpp bash build-linux.sh   # Eric's original

Windows (MSYS2 / MinGW-w64, PARI compiled for it):

    PARI=/c/pari-2.17.3 bash build-windows.sh && build/test_pdfilter.exe 100

Which faults does a version of the test see?

    bash fault_checks.sh 100 [test_pdfilter_oracle.cpp | test_pdfilter.cpp | test_pdfilter_linux.cpp | /tmp/eric_original.cpp]

## The four changes to Eric's version

1. **Varied block size.**  The filter works through its input in chunks of 64 polynomials, so calls with exactly 10000 polynomials always end
   with a chunk of 16.  Each call now uses one of 1, 7, 63, 64, 65, 127, 255, 1000, 4096, 10000.
2. **Shift.**  Three quarters of the polynomials are replaced by f(x+t), 1 <= |t| <= 30.  The discriminant does not change, but the
   coefficients become signed and large, the odd powers are no longer zero, and the zeros met by the Euclidean chain no longer are zero
   coefficients but appear by cancellation, which is where rounding errors act.
3. **Zero x^8 coefficient** in a quarter of the polynomials.  The Euclidean chain of (f, f') over F_q normally has the degrees 10, 9, 8, ..., 0;
   with a zero x^8 coefficient it has 10, 9, 6, ..., which the filter handles in its general chain (`chain_general.inc`) instead of
   `chain_block`.  25% of the polynomials that pass the discriminant test in the application behave like that.
4. **Seed and exit status.**  The seed is printed and can be given as a second argument; the exit status is 1 if the filter rejects a polynomial.

Why: Eric's first version took `tschirnhaus()` of `x^10 - 12x^5 - 4`.  Every output has the degenerate chain 10, 9, 7, ... (the seed is a
polynomial in x^5 and `g` has degree 2), so all lanes left the fast chain at every prime and the fast chain was never tested on a valid polynomial.
In the application, 74% of the passing polynomials have the generic chain 10, 9, 8, ..., 25% have 10, 9, 6, ..., 0.4% have 10, 9, 7, ....
Eric's second version (f = g(x^2), g a random quintic) is 99% generic.

## What the tests detect: the integer filter

`bash fault_checks.sh 100 test_pdfilter_oracle.cpp` and `bash fault_checks.sh 100 test_pdfilter.cpp` (Linux, seed 12345, 1,000,000 polynomials per
entry, 3,000,000 for the large primes; AVX2 build / baseline build).  For the oracle test the entry is "Total failures" with the missed rejections in
brackets; its total counts false and missed rejections, wrong stages, transformed copies that differ, and routing mismatches (stages of calls), so it
measures how much the test sees.  For `test_pdfilter.cpp` it is the number of valid polynomials rejected.  The faults are listed at the top of
`fault_checks.sh`.

| deliberate fault of the integer filter | oracle test | `test_pdfilter.cpp` |
|---|---|---|
| none (must be 0) | 0 (0) / 0 (0) | 0 / 0 |
| x87 arithmetic (`-mfpmath=387`) | 0 (0) / 0 (0) | 0 / 0 |
| floating point rounding mode up | 0 (0) / 0 (0) | 0 / 0 |
| large primes up to 2^29 (3,000,000) | 0 (0) / 0 (0) | 0 / 0 |
| small primes up to 2^14 | 0 (0) / 0 (0) | 0 / 0 |
| zero test +-q (`chain_block_s`) | 360 (0) / 546 (0) | 0 / 0 |
| zero test +-q (`chain_general_s`) | 23,283 (0) / 33,647 (0) | 349 / 349 |
| Euler: -R only as `negr` | 57,624 (1,515) / 83,008 (2,235) | 0 / 0 |
| top piece shifted unsigned (the start-up check switches the filter off) | 643,008 (643,008) / 643,008 (643,008) | 0 / 0 |
| loading constants 2^(20i) (the start-up check switches the filter off) | 643,008 (643,008) / 643,008 (643,008) | 0 / 0 |
| `chain_general` parity | 236,849 (4,217) / 281,816 (4,220) | 161,538 / 161,538 |
| last stage skipped | 7,280 (6,619) / 7,651 (6,619) | 0 / 0 |

This is the present filter (24 small primes below 2^13, 6 large; the run without a fault has 643,008 rejections).
The integer filter does no floating-point arithmetic, so the x87 and rounding faults cannot reach it; the large primes up to 2^29 are inside the range
where its arithmetic is proved exact.  The small primes up to 2^14 are above every proved limit, and `mont` itself does fail there (at q = 16369, 2.4 % of
the pairs in (-2q, 2q) give a result outside (-q, q)), but the filter's operands do not come near those pairs: no wrong verdict appeared in either test
or in a directed search, and whether one is possible there is open (see the size limits below).
The table was run again after every change of the filter.  With six small primes below 2^10: after `pdf_init`'s guard and start-up check were added
(2026-10-27) the loading constants failed the start-up check, so the filter switched itself off (every flag 1): the oracle test saw every rejection of
the run missed, the "must pass" test nothing, and the application stays correct, only slower (before: 1,132,612 / 1,383,679 and 751,240 / 751,240).
After the loads were changed to gather whole polynomials and the start-up check to call the loads themselves (2026-10-28), the unsigned top piece
(now in `ldcoef`) failed it too (before: 1,012,353 / 1,249,505 and 412,963 / 412,963).  With 24 small primes (2026-10-29) the faults that act at the
small stages are seen more often; with six they gave: zero test in `chain_block_s` 156 / 130 (routing only), in `chain_general_s` 3,461 / 4,766 and
3 / 3, Euler 186,381 (165,045) / 144,825 (108,090), parity 236,578 / 301,907 and 174,075 / 174,075, last stage 7,983 (7,354) / 8,340 (7,354).
Of the seven faults of the arithmetic, `test_pdfilter.cpp` sees one gross one that reaches the filter (the `chain_general` parity; before the start-up
check covered them, also the top piece and the loading constants), 349 polynomials of a second and nothing of the others; the oracle test sees all seven.
The Euler fault only loses rejections, which with 24 small stages mostly become later ones (47,854 wrong stages, 1,515 missed), the last-stage fault only
shows on polynomials that reach stage 30 (all 6,619 from the oracle test's `target` stratum), and the `chain_block_s` zero test changes few verdicts
(11 wrong stages, all in `crtdeg`) and is seen mainly by the routing check.

With 24 small primes the generated polynomials reach the large stages only through the strata `valid`, `rep` and `target`.  The large-stage code did
not change (`objdump -d`: the same instructions for the chains and loads as the build with six small primes below 2^10), and that build, run through the
oracle test with its primes, sends 788,246 of 1,000,000 polynomials through the large stages, including the CRT strata, with no failure.  Under
`-rm up`, `-rm down` and `-rm zero` the present filter has no failure either.

The oracle test on the floating-point filter of commit 327b04f (`FILTER=<its src>`, polynomials for its primes with `-oldrule`, AVX2 / baseline): x87
66 / 56 false and 296 / 311 missed rejections (61 and 304 in the AVX2 build when `build-linux.sh` still used `-march=x86-64-v3`); rounding up 37,329 /
37,320 false and 7,506 / 7,506 missed; rounding down and toward zero 717,978 missed rejections each, every rejection of the run (the "must pass" tests
report 0 there; under x87 `test_pdfilter.cpp` reports 194 / 194 with the present flags).  Details in Section 13 of `pdfilter_chains.tex`.

### How large may the primes be?

Section 8.8 of `pdfilter_chains.tex` ("Where the size limits come from").  Proved: `mont` (R = 2^16) is exact on operands in (-2q, 2q) for q <= 2^13,
`sredc` (R = 2^32) on differences of products for q <= 2^30 + 23171; for the whole filter, by following the bounds of the values through the chains,
small primes up to 8385 and large primes up to 1,097,010,265.  Computed exactly (complete searches, on the compiled code by `arith_checks.sh` and in
Python by `checks.py` C10): `mont` first fails at q = 8241 (first prime 8263, first prime = 1 mod 4 8269, first prime `pdf_init` could choose for
S = {2,5} 8329), `sredc` on a difference at q = 2^30 + 23279 (first prime 2^30 + 23293, first usable 2^30 + 23377).  Not known: whether the whole-filter
limits are sharp.  `pdf_init` enforces them: a usable prime outside [11, 8385] (small) or [11, 1,097,010,265] (large) is skipped with one
line on stderr (function `chosen`; the shipped search, below 2^13 and 2^25, never meets that case).  Test builds may switch the guard off with
`-DPDF_NO_LIMITS`; `fault_checks.sh` does so for the rows beyond the limits, and its row SG checks that the guard skips the primes above 8385.

The rows S8 ... L31 of `fault_checks.sh` move the primes across these limits (the expected result of each row was written down before it ran).
Oracle test / `test_pdfilter.cpp`, AVX2 / baseline, 1,000,000 polynomials, seed 12345; the last column runs directed inputs found by `limit_probe.c`
(a model of the chains with a hill-climbing search) through the same broken build with the oracle test.

| row | primes | proved? | oracle test | `test_pdfilter.cpp` | directed inputs |
|---|---|---|---|---|---|
| S8 | 241, 89, 41 | yes | 0 / 0 | 0 / 0 | |
| S13 | 8161 ... 8009 (4) | yes | 0 / 0 | 0 / 0 | |
| SC | 8369, 8329 (`mont` has counterexamples) | yes, whole-filter theorem | 0 / 0 | 0 / 0 | |
| SF | 8761 ... 8521 (6) | no | 0 / 0 | 0 / 0 | none found |
| S2 | 16369 ... 14449 | no | 0 / 0 | 0 / 0 | none found |
| SD | 19961 ... 19001 (14) | no | 0 / 0 | 0 / 0 | 1 input: 2,970 / 3,517 |
| SE | 23929 ... 21881 | no | 3,887 / 5,345 | 372 / 372 | 20 inputs: 2,191 / 3,590 |
| S15 | 32609 ... 30809 | no | 831,770 / 1,039,077 | 722,889 / 722,889 | 20 inputs: 2,206 / 3,366 |
| L24 | 16777121 ... 16776481 | yes | 0 / 0 | 0 / 0 | |
| L30 | 1073741689 ... 1073740609 | yes | 0 / 0 | 0 / 0 | |
| LC | 2^30 + 26385 ... 2^30 + 25617 (`sredc` has counterexamples) | yes, whole-filter theorem | 0 / 0 | 0 / 0 | |
| LF | 1097013329 ... 1097012809 | no | 0 / 0 | 0 / 0 | none found |
| L31 | 2147483489 ... 2147482481 | no | 0 / 0 | 0 / 0 | 10 inputs: 3,492 / 3,542 |
| SG | 8761 ... 8521 skipped by the guard; 8369 ... 8009 (7) | yes | 0 / 0 | 0 / 0 | |

Each row takes up to 24 small primes (fewer, in brackets, where its range has fewer usable ones) and 6 large ones.  At SD and L31 the filter gives
wrong verdicts that neither test sees in 1,000,000 polynomials; only the directed inputs show them.  The L31 inputs are lifted so that the 24 small
primes pass them (`limit_probe.c` tries up to 2^28 lifts).  With six small primes per row (until 2026-10-29) the numbers were: SD 0 / 0, 0 / 0,
1 input 2,240 / 3,495; SE 3,475 / 4,460, 252 / 252, 20 inputs 2,143 / 3,602; S15 666,466 / 814,620, 504,958 / 504,958, 20 inputs 2,145 / 3,409;
L31 16 inputs 3,493 / 3,493; every other row 0.

## What the versions of the test detect: the floating-point filter (commit 327b04f)

Valid polynomials wrongly rejected with a deliberately broken filter, 1,000,000 polynomials per entry (3,000,000 for the rare fault);
AVX2 build / baseline build.  `bash fault_checks.sh 100 [test]` reproduced this table with the floating-point filter.  The unbroken filter must give 0 / 0.

| deliberate fault of the filter | Eric's version, Linux | `test_pdfilter.cpp`, Linux | Eric's version, Windows | `test_pdfilter.cpp`, Windows |
|---|---|---|---|---|
| none (must be 0) | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| x87 arithmetic (`-mfpmath=387`) | 3 / 2 | 186 / 194 | 2 / 1 | 165 / 137 |
| floating point rounding mode up | 0 / 0 | 27,936 / 27,936 | 0 / 0 | 27,936 / 27,936 |
| large primes up to 2^29 (rare, 3,000,000) | 14 / 2 | 7 / 1 | 5 / 0 | 10 / 2 |
| small primes up to 2^14 (gross) | 71,857 / 8,074 | 90,384 / 8,242 | 72,474 / 8,139 | 89,967 / 8,486 |

Linux: Ubuntu 26.04 on WSL2, GCC 15.2, PARI 2.17.3 (Ubuntu package), seed 12345 for `test_pdfilter.cpp`, clock seed for Eric's version.
Windows: MSYS2 MinGW-w64 GCC 16.2, PARI 2.8 (Eric's file and `test_pdfilter_linux.cpp` do not compile there as they are: they need `long long` written before `pari.h`, `llabs` and
`isirreducible`, which is how it was run).  The identical rounding-up count on both systems is a coincidence: the failing polynomials differ and other seeds
give different counts.

Reading the table: the unbroken filter passes everything, so the two versions agree there.  Eric's version sees the gross small-prime fault and the rare large-prime
fault, but almost nothing of the x87 fault (a fixed block size of 10000 always ends in the same chunk) and nothing of the rounding-mode fault.  The shift and the
degenerate chains find the rounding-mode fault; the varied block sizes find the x87 fault.

A filter that rejects too little (rounding mode down, `-Ofast`) is not seen by any of them: the test only uses polynomials that must be kept.
