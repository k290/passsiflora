# Checking (H) at run time instead of assuming it

**Question.** Hypothesis (H) of the whitepaper (Proposition "What the exact test computes") says that no `int64` expression of the first two
sub-resultant iterations of upstream's `polDiscTest_cpu` overflows.  It is the one named assumption about the *input data*.  Could it be
removed by having the application detect the overflow and stop ("bail") when it happens?

**Short answer.** It can be done, and it would take (H) out of (C1) for the CPU build.  As proposed it is not worth doing:

1. Stopping is the wrong response.  A per-polynomial bypass gives the same guarantee without failing any task: the filter leaves the
   flag at 1 for a polynomial outside (H), and upstream's test decides it, as in the stock program.
2. Any detection that removes (H) has to check **every polynomial the filter rejects**, which is almost every polynomial, not just the
   survivors.  On the CPU that is 41% of the AVX2 filter's time, about +20% of the whole run.  On the GPU the rejected polynomials never
   leave the device, so the check would have to go into the kernel.
3. What it protects is narrow.  Outside (H), Passiflora differs from the stock program only when upstream *wrongly accepts* a polynomial
   that the filter rejects.  This can be constructed (below), but it needs coefficients near 2^62.  The harmful effect of the overflow is
   upstream *wrongly rejecting* polynomials that pass, which loses fields.  Stock and Passiflora lose those alike, and no bail recovers them.

Two cheaper steps give most of the benefit: compile `polDiscTest_cpu.cpp` with `-fwrapv` (no cost, no flag changed), and add a sampled
monitor.  Both are described in the last section.

Numbers are from `hcheck_checks.sh` (`test_hcheck.cpp`, `hcheck.h`), run on an Intel Xeon at 2.8 GHz (cloud VM) with GCC 13.3.  The
ratios matter here, not the absolute times.

## What (H) carries

The filter is proved against the definition of passing and does not need (H).  (H) enters in two places:

* **(C1)**: "flag 0 only where upstream's `polDiscTest_cpu` rejects".  The filter's rejection proves that the polynomial does not pass.
  Under (H), upstream's flag is the definition's, so it rejects too.
* **The upstream half of (C2)**, in the ledger row "for upstream (H)".  A signed overflow is undefined behaviour in C.

Outside (H), upstream's result can be wrong in both directions:

* **Wrongly rejects** a polynomial that passes: 844,265 of 1,000,000 large polynomials built to pass (`passbig`, Table "poldisc").  The
  filter never rejects a polynomial that passes, so Passiflora hands it to upstream, which rejects it as in the stock program.  The output
  is the same, and the field is lost in both programs.
* **Wrongly accepts** a polynomial that fails.  The filter may reject it, and then Passiflora's flag is 0 where the stock program's is 1.
  The polynomials written to the output stay the same, because PARI's field-discriminant test rejects such a polynomial afterwards.  But
  the counters "passing PolDisc test" and "passing irreducibility test" in the output file differ.  The whitepaper says this did not
  happen in its tests.  It can be constructed:

  Upstream sees the constant term a0 only through `B[9]^2*a0`, with `B[9] = 10/cont(f')`.  When `B[9]` is even, adding 2^62 to a0 changes
  that product by a multiple of 2^64.  So with wrap-around upstream computes exactly the same numbers as for the original polynomial.  Take
  a small g that passes, with `B[9]` even, and set f = g ± 2^62 at the constant term:

      g = x^10 + x^9 - 3x^8 - 2x^7 + x^6 - x^5 - x^3 + 3x^2 + 2x - 2
      f = g + 2^62:  upstream passes f (also with -fwrapv), the filter rejects f, (H) fails on f

  `test_hcheck` part 3 finds 100 such pairs in about a million small random polynomials.  In all 100, upstream passes f and the filter
  rejects it.  So in part 4, the application as it is today differs from the stock call on 100 of the 200 polynomials, and the bypass
  differs on 0.  None of the 83 results of the test units can be used for g: they are all ≡ x^10 mod 2, so their `B[9]` is odd.  This also
  shows that the column "rejects passing" of `test_poldisc.cpp` should be split by overflow, as its column "differ" is.  It counts filter
  rejections of polynomials that *upstream* passes, and this family would make `poldisc_checks.sh` fail on a case that (H) excludes.

## What a check has to cover

`hcheck.h` holds `nf_h_ok`, which replays every `int64` expression of upstream's lines 109–135, 154–156 and 201–266 (patched numbering) in
the same order.  It uses GCC's `__builtin_mul_overflow` / `__builtin_sub_overflow`, so the check itself has no undefined behaviour.
`nf_h_ok_fast` is the same check with a branch-free content and division.  On 5.2 million polynomials in 13 families, about 1.1 million of
them outside (H), both agree with the 128-bit replay of `test_poldisc.cpp` on every polynomial (part 1).

For (C1) without (H), the check must hold for every polynomial that gets **flag 0**.  The polynomials passed on to upstream do not matter,
because upstream decides them in both programs.  The polynomials with flag 0 are almost all of them.  So it cannot be limited to the
survivors, which would cost nothing.

## The ways to "bail"

| response when (H) fails | (H) out of (C1) | (H) out of (C2), upstream | output vs stock | the work unit |
|---|---|---|---|---|
| stop the task (bail completely) | yes, for runs that complete | yes, upstream never runs on it | identical when it completes | computation error, resent; the volunteer's quota drops |
| switch the filter off for the rest of the run | yes | no | identical | completes at stock speed (10–25× slower) after the first one |
| **per polynomial: leave the flag at 1** | yes | no | identical (*) | completes, no slowdown |
| nothing (today) | no, (H) assumed | no | identical where (H) holds | completes |

(*) "Identical" means: equal to `polDiscTest_cpu(..., NULL, ...)` as compiled into Passiflora, which is what the glue proposition states.
Whether the stock *binary* does the same thing on an overflow is not defined by its source (undefined behaviour).  In practice both wrap.

Stopping also fires on every overflow where upstream wrongly rejects, where Passiflora would have been identical anyway.  Its one benefit
over the bypass is that upstream's undefined behaviour never runs.  `-fwrapv` (below) gives that benefit for free.

## Cost

Per polynomial, in blocks of 10,240 (the application's buffer), on polynomials of the size of the real candidates (family `w8`, all
inside (H)), best of 5 runs:

| | ns |
|---|---|
| filter, AVX2 build | 156.8 |
| filter, baseline build | 188.1 |
| `nf_h_ok` (follows upstream literally) | 84.5 |
| `nf_h_ok_fast` | 64.4 (41% of the AVX2 filter, 34% of the baseline filter) |
| AVX2 filter + bypass | 232.2 |

The filter is about half of the CPU application's run time.  Whitepaper Table "speed" gives 84 ns of polDisc time per polynomial on
`sf5`.  The README gives 93 s for the 576,322,581 polynomials of `sf5_DS-16x12_Grp90of1600000`, which is 161 ns per polynomial.  So the
bypass would cost roughly **+20% of the CPU run time**.  This is an estimate across two machines and units, not a measurement in the
application.

The check is about 70–85 overflow-checked 64-bit multiplications and 30–50 subtractions.  AVX2 has no 64-bit multiply, so it stays
scalar.  A cheap screen cannot replace it:

* **Triangle-inequality bound** (every term taken in absolute value).  It lies 6.1 bits above the true largest value at the median and
  9.4 bits at the 99th percentile (part 2).  The real data have 3.5 bits of headroom (largest value 2^59.45).  For the `w8` family, 35% of
  polynomials would have a bound of 2^62 or more, although their true values stay below 2^59.
* **Floating-point screen.**  It would compute the signed values in doubles with a running error bound, use AVX2, and fall back to the
  exact check.  It would probably cost a few ns per polynomial, but I have not built or measured it.  It would bring floating-point
  reasoning back into the correctness argument, which the integer filter was written to remove.

The degR = 7 branch, where `R[8] = 0`, for example a9 = a8 = 0, has expressions of weight 18 instead of 13.  This is where the largest values
of the real results come from: 2^48.5 in an `sf5` result with a9 = a8 = 0.  It is also why random polynomials of moderate size (`w4`–`w10`)
overflow now and then.

## The GPU

The kernel rejects a polynomial without showing it to the CPU, so a check on the CPU covers only the GPU's survivors.  Removing (H) for the
GPU builds would need `nf_h_ok` in the kernel for every counted polynomial.  That is about 80 checked 64-bit products, each a
multi-instruction sequence on a GPU, plus a rebuilt CUDA fat binary and MUSA library.  It is not measured here (no GPU).  Without it the
GPU builds keep (H) whatever the CPU does.

## Recommendation

* **Do not bail.**  If (H) is to be removed for the CPU build, use the per-polynomial bypass: `hcheck.h`'s `nf_h_ok_fast` on every flag 0
  after `nfPdfFilter`, in `testPolys` and in `nfGpuFlushBatch`.  Accept about +20% CPU time and leave (H) in place for the GPU builds.  On
  the evidence (no overflow in 733,855 real candidates), that is a high price.
* **Compile `polDiscTest_cpu.cpp` with `-fwrapv`.**  Signed overflow then has a defined result, so upstream's behaviour outside (H) is
  determined by its source.  The upstream half of (C2) then holds without (H).  It changed no flag on 319,488 polynomials, 143,411 of them
  outside (H) (part 5), and the test only runs on survivors, so it costs nothing measurable.  The source is unchanged, only its compile line.
* **Optionally monitor (H) instead of assuming it blindly.**  Run `nf_h_ok_fast` on every survivor (free) and on every 1024th polynomial
  (about 0.06 ns per polynomial).  On the first failure, write one line on stderr and switch that run to the per-polynomial bypass.  A
  search region that drifts outside (H) produces overflows in bulk (`w11`: 1.2%, `w12`: 40%), and a sample finds them quickly.  (H) would
  remain a named assumption, but a monitored one.  In GPU mode the monitor sees the survivors and the start-up window, when the CPU loops
  also run.
* **Whitepaper:** "that did not happen in the tests" can be replaced by the constructed example.  The difference is possible, needs
  coefficients near 2^62, and only changes counters.

## Reproduce

    UPSTREAM=<GetDecicsSrc> bash hcheck_checks.sh 400000 1 -hits hits.txt     # about 1 minute; hits.txt: the 83 results, one per line

(`-hits` is optional.  The file is the same as that of `poldisc_checks.sh`.)
