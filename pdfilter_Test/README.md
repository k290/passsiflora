# pdfilter_Test

Tests of the pre-filter `pdfilter.f90` (a copy of `src/pdfilter.f90`).  Every polynomial that the test generates has a discriminant of
the form (a number built from the primes 2 and 5) times a perfect square, so the filter **must keep** it.  A failure means that the
filter wrongly rejected a valid polynomial.

## Files

| file | what it is |
|---|---|
| `test_pdfilter.cpp` | **the test to use**: Eric's second version plus four small changes (below), rewritten so that it builds on Linux and on Windows. |
| `test_pdfilter_linux.cpp` | Eric's second version with the same four changes applied as a small patch (+47 / -15 lines; Linux only, PARI 2.9 or later).  `git diff ccb629f:pdfilter_Test/test_pdfilter_eric.cpp HEAD:pdfilter_Test/test_pdfilter_linux.cpp` shows the diff against Eric's file.  It gives the same numbers as `test_pdfilter.cpp` on Linux. |
| `CHANGES.md` | the changes to Eric's test, each with its reason, and how the filter picks `chain_block` or `chain_general`. |
| `WORKED_EXAMPLE.md` | a worked example for students: small polynomials followed step by step through `chain_block` and `chain_general`. |
| `chain_trace.cpp` | the program that produced those traces (`TEST=chain_trace.cpp bash build-linux.sh`, then `build/chain_trace show ...` or `build/chain_trace test`). |
| Eric's original | not in the tree; `git show ccb629f:pdfilter_Test/test_pdfilter_eric.cpp` (it is the file he sent). |
| `test_pdfilter_full.cpp` | an earlier, larger harness (several generators, a corpus of real polynomials, rounding-mode options). |
| `build-linux.sh`, `build-windows.sh`, `pari-config.sh` | build the filter (AVX2/FMA and baseline) and the test. |
| `fault_checks.sh`, `rmshim.c` | break the filter on purpose and count how many valid polynomials each version of the test sees rejected. |
| `pdfilter.f90` | the filter. |
| `test_pdfilter` | Eric's old Linux binary.  It is truncated (4 MB of the ~29 MB its ELF header expects) and segfaults; use the build scripts. |

## Run

Linux (gfortran, g++ and PARI; `PARI=/path/to/compiled/pari-tree`, otherwise Eric's tree, otherwise the system PARI, `apt install libpari-dev`):

    bash build-linux.sh && build/test_pdfilter 100          # 100 blocks of 10000 polynomials; a second argument is the seed
    TEST=test_pdfilter_linux.cpp bash build-linux.sh        # builds the small patch of Eric's file instead (build/test_pdfilter_linux)
    git show ccb629f:pdfilter_Test/test_pdfilter_eric.cpp > /tmp/eric_original.cpp && TEST=/tmp/eric_original.cpp bash build-linux.sh   # Eric's original

Windows (MSYS2 / MinGW-w64, PARI compiled for it):

    PARI=/c/pari-2.17.3 bash build-windows.sh && build/test_pdfilter.exe 100

Which faults does a version of the test see?

    bash fault_checks.sh 100 [test_pdfilter.cpp | test_pdfilter_linux.cpp | /tmp/eric_original.cpp]

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

## What the versions of the test detect

Valid polynomials wrongly rejected with a deliberately broken filter, 1,000,000 polynomials per entry (3,000,000 for the rare fault);
AVX2 build / baseline build.  `bash fault_checks.sh 100 [test]` reproduces a table.  The unbroken filter must give 0 / 0.

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
