# How Passiflora works — and why the results are identical

## The stock CPU app spends 98 % of its time in one test

`polDiscTest_cpu` computes the discriminant of each candidate decic with GMP
(numbers grow to ~1000 bits), divides out the primes in the work unit's set S, and
asks whether the rest is a perfect square. Only about one polynomial in 10^7 passes.

## 1. An exact modular pre-filter (`src/pdfilter.f90`)

The value upstream tests differs from Disc(A) only by rational square factors
(the content of B contributes c^10, the last sub-resultant step (B0/h)^(2k)), so the
question is the same for Disc(A). This holds whenever upstream's first two sub-resultant
iterations, which it does in int64 before switching to GMP, do not overflow; when they do,
the number wraps around and upstream's result can be wrong (`docs/whitepaper.tex`,
Proposition "What the exact test computes", measured by `pdfilter_Test/poldisc_checks.sh`;
see "The int64 numbers of the exact test" below). Take a prime q ≡ 1 (mod 4) such that every p in S
is a quadratic residue mod q. If the polynomial passed, Disc(A) would be
± (product of primes in S) × (a square), so its Legendre symbol mod q could only be 0
or 1. **A Legendre symbol of −1 is therefore a proof that the polynomial fails.**

The filter evaluates that symbol for a few primes and clears the flag only on such a
proof. Every other polynomial — including any one where the computation is
inconclusive — goes on to the *unchanged* GMP test (through the `polyMask`
argument upstream already uses for its GPU builds), so the set of polynomials that
pass is exactly the stock set. Each prime rejects about half of the candidates.

Disc(A) mod q comes from a Euclid chain over F_q:

* **Fast path** — generic degree sequence (each remainder one degree lower), two
  pseudo-remainder rounds per step. Every leading-coefficient power that appears has
  an even exponent, so χ(Disc A) = χ(final constant) and no modular inverses are needed.
* **General path** (`chain_general`) — for lanes where a leading coefficient
  vanishes (some families do this structurally, e.g. x^10 + b3·x^7 + …). Polynomials
  are stored left-aligned so one elementary step has the same formula in every SIMD
  lane; the odd powers of lc(V) from Res(V,U′) = lc(V)^(deg U′ − deg U + deg V)·Res(V,U)
  are tracked explicitly.

All arithmetic is integer arithmetic: signed Montgomery multiplication with R = 2^16 in
16-bit lanes for the 30 primes (2^8 < q < 2^13), and a three-piece loading of the 64-bit
coefficients (one reduction with R = 2^32 and 64-bit products). The filter does no floating-point operation, so its result does not depend
on the rounding mode, on x87 or SSE code, or on `-ffast-math`. `pdf_init` chooses a prime
only inside the proved limits of the arithmetic (11 ≤ q ≤ 8385; any other usable prime is skipped with a line on stderr), and at the end checks the
compiled arithmetic on known inputs; if that check fails it uses no prime, so every flag is
1 and the exact test alone decides. The proofs (exactness, ranges, limits, the start-up
check, the masked loop of the general path) are in `docs/whitepaper.tex`.
The CPU build contains an AVX2 copy (compiled with `-mavx2`; used when the CPU has AVX2
and FMA) and a baseline (SSE2) copy and picks one at run time.

### The one assumption: upstream's int64 steps do not overflow

Upstream's test does its first two sub-resultant iterations in int64 and does not check that the
products fit. The search regions keep them below 2^60: the largest value seen on real candidates is
2^59.45, and none of the 14.8 billion polynomials of the 31 test work units overflows. Passiflora relies
on this, as the stock program does, and says so in the source: above upstream's first iteration
(`polDiscTest_cpu.cpp`, through the patch) and where the filter is called (`testPolys`,
`nfGpuFlushBatch`).

If a polynomial did overflow, the value would wrap around (`polDiscTest_cpu.cpp` is compiled with
`-fwrapv`, so this is defined behaviour) and upstream's verdict could be wrong in either direction.
Usually it rejects a polynomial that passes; the stock program does the same, so the output is unchanged
and the field is lost in both. Rarely it passes one that does not; the filter, which is exact, would then
reject it, and the counters of the output would differ from the stock program's (constant terms near
2^62 do this). Before running a work unit with much larger polynomials than today's, check them with
`pdfilter_Test/poldisc_checks.sh` (`CAND=<file>`, which reports every overflow); `pdfilter_Test/POLDISC_INT64.md` compares the ways the
program itself could check.

## 2. Exact shortcuts in the enumeration (patch)

* `sqrt(X) <= B` is replaced by `0 <= X <= T`, where T is the largest double with
  `sqrt(T) <= B` (found once per bound with `nextafter`). `sqrt` is correctly rounded
  and monotone, so both forms take the same branch for every double X, including NaN,
  negative values and −0.
* In the innermost (m1) loop both a4 conditions are convex quadratics in m1. The range
  is narrowed to where they can hold, solving with a margin (1e-9 relative to the
  largest term, ±2 in m1) far larger than the floating-point error of the original
  test; every skipped m1 provably fails the original test. Nothing is narrowed when
  int64 overflow could matter.

The enumeration is compiled exactly like upstream (no FMA contraction), so all
floating-point loop bounds are bit-identical.

## 3. The GPU version (`src/nf_kernel.h`, `src/nf_gpu_*.cpp`)

The stock GPU app offloads only the discriminant test and feeds it from one CPU thread
(88 bytes per polynomial over PCIe, ~50 GB for one work unit), so a modern card idles.

Passiflora's CPU code walks only the outer loops and queues a *node* per fixed
(a1…a3) state — 18 doubles and 35 integers. The OpenCL kernel runs one work item per
(node, m2): the complete m1 loop, the bb4…bb9 bound tests in the original order, and an
integer version of the filter above (unsigned 32-bit Montgomery arithmetic with R = 2^16,
residues kept in [0, 2q)) at the 24 smallest usable primes in (256, 8192); it is exact like the CPU filter
(`docs/whitepaper.tex`, Theorem E). It counts every
polynomial of the original enumeration (the application reports that number) and returns
only the ~1 in 10^5 that survive the filter. The CPU then finishes those exactly as
before: the Fortran filter with its own 30 primes, the GMP test, PARI's irreducibility test and
the field discriminant.

Bit-for-bit agreement with the CPU code relies on:

* doubles evaluated in the same order, `FP_CONTRACT OFF`, IEEE `sqrt` and division;
* integer arithmetic done in `ulong` so it wraps like the CPU's `long long`;
* double → int64 conversion following the x86 rule (NaN or out of range → INT64_MIN),
  which the CPU code can hit in degenerate cases;
* survivors sorted back into enumeration order before the exact tests.

Checkpoints are written between batches and resume at the first node of the next
batch. The batch size adapts to about 100 ms per kernel launch.

## GPU back ends

The kernel is written once (`src/nf_kernel.h`) in a small portable dialect of OpenCL C / CUDA C; a port header per
back end maps the few differences (address spaces, atomics, `as_float`, 64-bit types):

| Back end | Port header | Host code | Notes |
|---|---|---|---|
| OpenCL | `nf_port_opencl.h` | `nf_gpu_opencl.cpp` | kernel compiled at run time, `FP_CONTRACT OFF` |
| MUSA (Moore Threads) | `nf_port_cuda.h` | `nf_gpu_musa.cpp` + `nf_musa_dev.cu` | the device part (kernel + CUDA-runtime-style wrapper) is converted with `musify-text` and built with `mcc` into `libnf_musa.so`; the program `dlopen`s it, so a computer without the vendor runtime or card just runs on the CPU. BOINC cannot see these GPUs: the task is a CPU task that uses the card by itself |
| CUDA | `nf_port_cuda.h` | `nf_gpu_cuda.cpp` | kernel compiled once with `nvcc --fmad=false` into a fat binary (cubin for sm_50 … sm_120 + PTX) that is embedded in the program and loaded through the CUDA *driver* API, which is itself loaded at run time (`nvcuda.dll` / `libcuda.so.1`): nothing CUDA is linked, so the program starts without a CUDA driver and falls back to the CPU |

The CUDA fat binary and the MUSA library are built apart from the program, so the program checks at start-up that they read
the prime table it builds: the table layout has a number (`NF_TABLE_FORMAT`), the kernel entry point carries it (`nf_enum_t2`),
and the MUSA library reports it and checks the table.  Device code built for another layout is never run.

The queueing, batch sizing, filter-prime table and the sorting of the survivors are shared (`nf_gpu_common.cpp`);
a back end implements the five functions of `nf_gpu_backend.h`. `nf_gpu_cuda.cpp` also builds against other
CUDA-driver-API-compatible runtimes (`-DNF_DRV_PREFIX`, `-DNF_DRV_LIB`; HIP, MUSA) given the code object of that
toolchain.

## GPU self-check

At start-up every back end runs the self-test kernel `nf_selftest_t2` once on the device: `redc`, `coef_mod`,
the symbol at each filter prime and `survives` on fixed inputs (the extremes and pseudo-random values; 240
polynomials, among them abnormal chains and discriminants divisible by a table prime), 7,920 values for S = {2,5}.
The host computes what they must be without the kernel's code (exact arithmetic modulo q; the symbol from the
resultant by Euclid with inverses). On any difference, or with device code that has no self-test (an older CUDA
fat binary, OpenCL kernel or `libnf_musa.so`), the task runs on the CPU.

Then, at the start of every GPU task, the CPU loops run next to the GPU for the first 5 million polynomials and
the two counts are compared batch by batch. With `NF_GPUCHECK=1` (the whole task) the CPU also runs the kernel's own
filter, compiled for the CPU, on every polynomial it counts; the polynomials it keeps must be exactly those the GPU
returns, in the same order. A mismatch stops the task with an error (`GPU self-check FAILED`, exit status 6) instead
of sending results from a GPU or driver that miscomputes. The default window compares only the counts: the filter
comparison costs about a microsecond per polynomial on the CPU, and with 24 filter primes few polynomials survive in
the first 5 million (0 to 10 in the test units), so it is the self-test that checks the filter's values.
`NF_SELFCHECK=<n>` changes the number of polynomials (0 = off), `NF_GPUCHECK=1` checks the whole run;
`NF_SELFCHECK_FAULT=1` injects a miscount, `NF_SELFCHECK_FAULT=2` a lost survivor (with `NF_GPUCHECK=1`) and
`NF_SELFTEST_FAULT=1` a wrong self-test value, to test the failure paths.

## Other platforms

* **aarch64 Linux / macOS arm64:** same sources. Compiled with `-ffp-contract=off`: the enumeration must round
  exactly like the x86 build. Where `long double` is only a double (Apple silicon) the m1-range narrowing is
  disabled, because its error margin assumes an extended-precision `long double` (x86: 64-bit mantissa, aarch64
  Linux: 113 bits).
* **macOS:** one universal program: the arm64 build and an x86-64 build like Linux's (AVX2 and baseline filter,
  picked at run time), joined with `lipo`. Both BOINC platform names, `arm64-apple-darwin` and `x86_64-apple-darwin`,
  list it: macOS starts the half that matches the processor.

## Validation

* All 31 work units of upstream's `Test/dataFiles` against their `.truth` outputs
  (polynomials and every counter): CPU on Linux, CPU on Windows (subset), GPU on Windows.
* `NF_VERIFY=1` also runs the exact GMP test on every polynomial and aborts if the filter
  ever disagrees; `NF_GPUCHECK=1` runs the CPU loops next to the GPU ones and aborts if
  the polynomial counts or the filter's survivors differ.
* Debug builds with `-DNF_CHECK_SQRT -DNF_CHECK_PRUNE` evaluate the original enumeration
  conditions next to the shortcuts.
* Kill/resume tests for the CPU and the GPU app.

## Environment switches (testing)

| Variable | Effect |
|---|---|
| `NF_NOFILTER=1` | stock code path (no pre-filter) |
| `NF_VERIFY=1` | also run the exact test on every polynomial, abort on mismatch |
| `NF_FORCE_BASE=1` | use the baseline (SSE2) filter build |
| `NF_GENERAL_ALL=1` | send every lane through the general chain |
| `NF_GPU=0` | never use the GPU |
| `NF_GPUCHECK=1` | GPU mode plus the CPU loops; abort if the counts or the filter's survivors differ |
| `NF_SELFCHECK=N` | CPU loops next to the GPU for the first N polynomials of a task (default 5000000; 0 = off) |
| `NF_SELFCHECK_FAULT=1`, `=2` | test hooks: pretend the GPU miscounted (1) or, with `NF_GPUCHECK=1`, lost a filter survivor (2) |
| `NF_SELFTEST_FAULT=1` | test hook: pretend the device got one value of the start-up self-test wrong |
| `NF_OCL_ANY=1` | GPU mode on the first OpenCL device without BOINC (`NF_OCL_PLATFORM`, `NF_OCL_DEVICE` pick another) |
| `NF_GPU_ITEMS=N` | fixed work items per launch instead of the adaptive size |
| `NF_MAXPOLYS=N` | stop after N polynomials and print a timing split (CPU mode) |
| `NF_PDFDEBUG=1` | print per-prime survivor counts for the first batch |
