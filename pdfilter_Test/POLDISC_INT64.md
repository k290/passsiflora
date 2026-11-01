# The int64 steps of upstream's exact test

Upstream's `polDiscTest_cpu` runs the first two sub-resultant iterations on f and f′ in int64 machine integers, then switches to GMP.
It never checks that those products and differences fit.  The search regions keep them small: the largest value seen on real
candidates is 2^59.45, against a limit of 2^63.  The author of GetDecics has said that the real work units stay inside.  The
whitepaper names this condition hypothesis (H).

Passiflora relies on it, as the stock program does.  This note says what is in place so that nobody forgets that, why it matters,
and which other options were measured.

## What is in place

* **The source says it** where a reader needs it.  There is a comment above upstream's first int64 iteration in
  `polDiscTest_cpu.cpp`, added through the patch, and a comment at both calls of the pre-filter (`testPolys`, `nfGpuFlushBatch` in
  `TgtMartinet.cpp`).
* **The documents say it.**  The whitepaper now names the assumption in Section 1 ("The one assumption about the data") as well as in
  Section 3, the contract and the ledger.  `docs/TECHNICAL.md` has a short subsection on it.
* **`-fwrapv`** for `polDiscTest_cpu.cpp` in every build, and in `poldisc_checks.sh`.  An overflow is then defined behaviour: it
  wraps modulo 2^64, as it does in practice anyway.  It is a compile flag, not code, and it changed no flag on 319,488 polynomials,
  143,411 of which overflow.
* **`poldisc_checks.sh`** reports every overflow, on any file of candidates (`CAND=`).  Its family `wrap` builds the case below.

## Why it matters

Outside the int64 range upstream's verdict can be wrong in either direction.

* **It rejects a polynomial that passes.**  This happened on 844,265 of 10^6 in the family `passbig`.  The stock program does the
  same, so the output is unchanged, and the field is lost in both programs.
* **It passes a polynomial that does not.**  The filter is exact and may reject it, so Passiflora's counters would differ from the
  stock program's.  The polynomials written out stay the same, because PARI's field-discriminant test rejects such a polynomial.
  Example: upstream sees the constant term a0 only through `B9^2*a0`, with `B9 = 10/content(f')`.  When B9 is even, adding 2^62 to
  a0 changes nothing after the wrap:

      g = x^10 + x^9 - 3x^8 - 2x^7 + x^6 - x^5 - x^3 + 3x^2 + 2x - 2      passes
      g + 2^62                                                          upstream passes it, the filter rejects it

  Family `wrap`, `QUICK=1` run: on all 200 polynomials upstream passes, the definition does not, and the filter rejects (column
  `rej-up-ovf`).

**Evidence that the real work units stay inside:**
* No overflow on 733,855 real candidates and shifts (whitepaper, Table 1).
* With a temporary copy of the application that checked every polynomial (built outside the repository), none of the
  14,846,153,955 polynomials of the 31 test work units overflows, and all 31 gave their truth files.

## Options that were measured

Measured on an Intel Xeon at 2.8 GHz (cloud VM) with GCC 13.3, with programs kept outside the repository.  The exact test takes
6.07 µs per polynomial it runs on; the AVX2 filter about 141 ns per polynomial.

| option | code | cost | effect |
|---|---|---|---|
| **document it, `-fwrapv`** (chosen) | comments, one compile flag | none | (H) stays an assumption, stated where it matters |
| sample at run time: every survivor and every 1024th polynomial; after a find, check every polynomial and leave rejected ones that overflow to the exact test | a 120-line header and 3 calls in the patch | 0.1% of the filter | stock flags from the first find on; misses rare overflows before it; in GPU mode not the kernel's rejections |
| check every polynomial (as above, from the start) | the same | about half the filter's time (one unit: polDisc 3.1 → 4.8 s) | (C1) without (H) on the CPU; the GPU would need the check in the kernel |
| stop the task on an overflow | the same check | as above | turns work units the stock program completes into errors |
| a cheap bound instead of the exact check (every term in absolute value) | small | small | does not work: 6–9 bits above the true values, against 3.5 bits of headroom |
| **fix upstream: the int64 steps with overflow checks, GMP on overflow** | the int64 lines rewritten with overflow builtins, plus a few lines for the GMP start (a few dozen lines in `polDiscTest_cpu.cpp`) | +1% per exact test (6.07 → 6.13 µs, measured as the check plus the test as shipped) | the exact test equals the definition for every input |
| fix upstream: drop the int64 shortcut, all in GMP | removes the int64 code (about 170 lines), adds 7, enlarges 4 arrays | +26% per exact test (7.67 µs) | as above |

**Fixing upstream.**  Both fixes give the definition's verdict on every polynomial; they agreed with PARI on all families,
`passbig` and `wrap` included, where the test as shipped is wrong on 16,869 of 20,000 and on all 300.  In Passiflora the exact test
runs only on the pre-filter's survivors, so either fix costs nothing measurable there.  In the stock program it runs on every
polynomial, so only the checked variant (+1%) is cheap enough.

Done in Passiflora alone, either fix would make its output differ from the stock program's on overflowing polynomials, where the
stock program is wrong.  That is why it is not applied here.  The checked variant is the change to propose upstream: stock and
Passiflora would then change together and stay identical, at about 1% for the stock program.
