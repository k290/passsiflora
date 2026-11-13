# 🌼 Passiflora — an exact pre-filter for GetDecics

Passiflora is a faster build of the NumberFields@home application `GetDecics` ("Get Decic Fields"), packaged as a BOINC
anonymous-platform app.  Its output is the stock program's, byte for byte apart from the "Elapsed Time" line.  The speed comes from a
**pre-filter**: before GetDecics' exact discriminant test, which dominates the stock application's running time, Passiflora discards
every polynomial whose discriminant is provably not of the required form.  The GPU version also runs the enumeration itself on the
graphics card (OpenCL for NVIDIA, AMD and Intel; CUDA; MUSA).

> The passion flower has 5 sepals and 5 petals — ten parts. NumberFields@home searches for number fields of degree ten, built as
> degree-5 extensions of a quadratic field.

[![Latest release](https://img.shields.io/github/v/release/k290/passsiflora?label=latest&color=8e44ad)](https://github.com/k290/passsiflora/releases/latest)
[![License: GPL v2+](https://img.shields.io/badge/license-GPLv2%2B-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)

---

## How it works

GetDecics enumerates monic degree-10 polynomials and keeps those whose discriminant is, up to the primes of a set `S` given by the
work unit (`{2, 5}` in the current searches), a perfect square.  If a polynomial passes, its discriminant is a square or zero modulo
every prime `q ≡ 1 (mod 4)` at which the primes of `S` are squares.  So a Legendre symbol of −1 at one such prime proves that the
polynomial fails.  The pre-filter computes these symbols:

* **CPU:** at 30 primes below 2^13, in 16-bit SIMD lanes, each by one of two Euclidean chains on
  `f` and `f'`: a fast one for the usual case and a general one for every other.  Arithmetic is integer-only (signed Montgomery
  reduction); there is no floating point anywhere in the filter.
* **GPU:** at the 24 smallest suitable primes in (256, 8192), inside the kernel, with unsigned Montgomery arithmetic.
* Every polynomial that is not rejected goes to GetDecics' own exact test, unchanged.

## Why the results can be trusted

The full argument, with every proof, is in [docs/whitepaper.pdf](docs/whitepaper.pdf).  In short:

**Proved.**  The filter rejects a polynomial exactly when its discriminant has symbol −1 at some filter prime, so it never rejects one
that passes (Theorems A–C).  The program as written makes exactly those decisions for every input, independently of rounding mode,
x87 or SSE code and compiler flags such as `-ffast-math` (Theorem D); the same holds for the GPU kernel (Theorem E).  The integer
arithmetic is exact for primes up to 8385; `pdf_init` refuses any prime outside this limit.

**Checked.**
* The compiled filter on **every** polynomial of complete sets at tiny primes: 25.5 billion polynomials, 204 billion verdicts, none wrong.
* An oracle test that computes every symbol independently and reaches every stage and path of the filter; it detects all seven
  deliberately injected arithmetic faults.
* GetDecics' exact test against the mathematical definition (PARI) on 10.7 million polynomials.
* All **31** upstream test work units give their reference outputs in both CPU builds (AVX2 and baseline); on an NVIDIA GPU
  (GTX 1650), all 31 through OpenCL and CUDA on Windows and through CUDA on Linux, with the programs of the release workflow.

**Assumed.**  GetDecics' exact test runs its first two steps in 64-bit integers without overflow checks; Passiflora, like the stock
program, relies on the search regions keeping them in range (no overflow among the 14.8 billion polynomials of the 31 test work units).
The compiler must implement two integer conventions the Fortran standard leaves open (gfortran does); the program tests them at
every start (below).

**At run time.**  The start-up checks guard the build, not the machine.  Integer instructions have no rounding modes, so a binary
computes the same values on every processor of its architecture; what can be wrong is the code its compiler produced.  The filter's
own tests ran on x86-64 with GCC 15.2 and 16.2, but the releases are also built with MinGW-w64 and for ARM
Linux, macOS and FreeBSD.  The release workflow runs each of those programs on the upstream work units (all 31, or two small ones where
the system is emulated); beyond that, the CPU filter's start-up check of its compiled arithmetic is the one test of their code.  The GPU kernel was tested on one card, and under OpenCL each machine's graphics driver compiles it, so every GPU task
first runs a self-test of the kernel's filter against values computed on the host and compares its polynomial count with the CPU over
the first 5 million polynomials.  A failed check reports nothing wrong: the CPU filter switches off and the exact test decides every
polynomial; a GPU task stops with an error (exit status 6).  `NF_GPUCHECK=1` compares the GPU's survivors with the CPU for the whole task.

## How fast?

The integer filter is faster than the floating-point filter it replaced, on the CPU and on the GPU.  The measurements, and where the
time goes, are in the whitepaper, Section 9 "Performance".

## Downloads

**[Latest release](https://github.com/k290/passsiflora/releases/latest).**  Releases from v0.4.0 on contain the integer filter described
above; v0.1.0 to v0.3.0 ([alplix/passiflora](https://github.com/alplix/passiflora/releases)) contain the earlier floating-point filter.
The number after the prefix is the version (`062` = v0.6.2).

| starts with | what it is | for |
|---|---|---|
| `CPU_` | the CPU app | everybody |
| `GPU_` | the OpenCL GPU app, with an `app_config.xml` | NVIDIA, AMD or Intel cards with double-precision OpenCL |
| `CPU-GPU_` | both apps, with an `app_info.xml` that lists both and an `app_config.xml` for the GPU | as `GPU_`; BOINC runs one app on the card and one on the processor |
| `GPU-CUDA_`, `CPU-GPU-CUDA_` | the same with CUDA instead of OpenCL (driver only, no toolkit) | NVIDIA cards, Maxwell to Blackwell |
| `GPU-MUSA_` | the Moore Threads app (a CPU task to BOINC; one task per card also uses the card) | Moore Threads cards, Linux x86-64 |

Systems: Windows x86-64; Linux x86-64 (CPU fully static; GPU needs glibc 2.31+, MUSA 2.34+); Linux aarch64; Android aarch64 (Android 10+; NumberFields'
Android app); FreeBSD 14+ x86-64; Linux armhf; macOS 15+ on Intel and Apple silicon (one
universal program, CPU only).  Install steps are in `INSTALL.txt`, attached to every release; each package has a `README.txt`.

## Running the programs

**Under BOINC (production).**  Install a package as `INSTALL.txt` describes: the program files, `app_info.xml` and, in the GPU
packages, `app_config.xml` go into the project's folder in BOINC's data directory.  Both XML files live on the volunteer's computer:
under anonymous platform the project sends neither, and a package only provides starting versions of them, which the volunteer may
edit.

* `app_info.xml` (required) declares the app versions: which program runs under which BOINC platform names and plan classes, on
  which kind of GPU, and how much CPU BOINC sets aside for each task.  The client reports these versions to the project when it asks
  for work.  After editing it, restart the client.
* `app_config.xml` (optional) tunes how the client schedules the GPU versions it already has: `gpu_usage` 0.25 runs four tasks per
  card, `cpu_usage` is the CPU share of a GPU task and overrides the one in `app_info.xml` (both are 0.5 in the packages).  BOINC
  Manager's *Options → Read config files* applies an edit without a restart.  A volunteer who already has an `app_config.xml` for
  NumberFields should add the package's `<app>` entry to it instead of replacing it.

BOINC then runs the CPU program as a CPU task and a GPU program on the GPU it assigns to the task (it passes the device in
`init_data.xml` and as `--device N`).  A GPU program that cannot use its GPU ends the task with an error whose cause is in
the task's stderr (shown on the project's website); a temporary problem, such as GPU memory in use by another program, makes BOINC
restart the task a minute later.  The Moore Threads program is a CPU task to BOINC; one task at a time also uses each card.

**By hand (development and testing).**  Copy a work unit to `in` in an empty directory and start the program there; it writes `out`,
and its log to `stderr.txt`:

```bash
mkdir run && cd run && cp ../upstream/Test/dataFiles/sf5_DS-14x11_Grp16409of705883.dat in
/path/to/GetDecics_cpu                 # the CPU program
/path/to/GetDecics_gpu                 # OpenCL: the first OpenCL GPU; stderr.txt lists them, --device N picks another
/path/to/GetDecics_cuda --device 1     # CUDA: the second NVIDIA card
```

A GPU program uses the GPU that an `init_data.xml` in the directory names (as the BOINC client writes it; upstream's examples are in
`Test/exampleInitDataXmlFiles`), else `--device N`, else the first GPU; the first lines of `stderr.txt` say which.  Compare `out` with
the unit's `wu_<name>.out.truth`, apart from the "Elapsed Time" line: `build/test-suite-unix.sh` does that for all 31 units (below).
`NF_GPU=0` runs a GPU program on the CPU; the other `NF_*` variables are test switches ([docs/TECHNICAL.md](docs/TECHNICAL.md)).
What Passiflora changes in how the programs are built and run, compared with upstream GetDecics, and why:
[docs/CHANGES-FROM-GETDECICS.md](docs/CHANGES-FROM-GETDECICS.md).

## The repository

Only code written for Passiflora is here; it is applied on top of the upstream GetDecics source (commit `638ca4e`).

| path | what it is |
|---|---|
| `src/` | Passiflora's sources, copied next to upstream's: the CPU pre-filter (`pdfilter.f90`), the GPU kernel (`nf_kernel.h`, `nf_kernel.cu`, `nf_port_*.h`), the GPU host code and back ends (`nf_gpu_common.cpp`, `nf_gpu_opencl.cpp`, `nf_gpu_cuda.cpp`, `nf_gpu_musa.cpp`, `nf_musa_dev.cu`) and the CPU-only stub (`nf_gpu_stub.cpp`). |
| `patches/passiflora.patch` | the changes to four upstream files: `GetDecics.cpp`, `TgtMartinet.cpp`, `TgtMartinet.h` (the calls of the filter and the GPU mode) and `polDiscTest_cpu.cpp` (32-bit `long` support and comments). |
| `build/` | `prepare-src.sh` assembles a source tree (upstream + patch + `src/`); `test-suite-unix.sh` runs work units against their reference outputs; the `build-*.sh` scripts and `make-release.py` are the release builds, with the paths of the release machine. |
| `docs/` | `whitepaper.tex`/`.pdf` (algorithms, proofs, evidence), `TECHNICAL.md` (a shorter technical overview), `CHANGES-FROM-GETDECICS.md` (what differs from upstream and why), `BUILDING.md` (how to build), `INSTALL.txt` (for users). |
| `bench/` | micro-benchmarks of the filter (`bench_filter.f90`, `timing_probe.f90`). |
| `pdfilter_Test/` | the test suite of the pre-filter and of GetDecics' exact test; see its [README](pdfilter_Test/README.md). |

**Build.**  Check out upstream at `638ca4e` (step 1 of [docs/BUILDING.md](docs/BUILDING.md)), assemble a source tree, and compile it
as BUILDING.md describes:

```bash
UPSTREAM=$PWD/upstream OUT=$PWD/build_src bash passiflora/build/prepare-src.sh
```

**Check a build against the 31 reference work units** (one line per unit, then a summary; the program needs an absolute path;
`NF_FORCE_BASE=1` tests the baseline build):

```bash
bash passiflora/build/test-suite-unix.sh $PWD/GetDecics_cpu $PWD/upstream/Test/dataFiles 8
```

**Run the verification** (Linux or WSL, in `pdfilter_Test`, with gfortran, g++ and PARI; each script builds what it needs from `src/`):

```bash
TEST=test_pdfilter_oracle.cpp bash build-linux.sh && build/test_pdfilter_oracle 100 12345   # oracle test, 10^6 polynomials
bash tiny_checks.sh                     # complete sets at tiny primes (QUICK=1: a few minutes)
bash arith_checks.sh                    # the compiled reductions and the size limits (QUICK=1: about a minute)
bash fault_checks.sh 100 test_pdfilter_oracle.cpp                     # deliberate faults and moved primes
UPSTREAM=<build_src>/GetDecicsSrc bash poldisc_checks.sh                 # GetDecics' exact test against the definition
TEST=test_gpu_filter.cpp bash build-linux.sh && build/test_gpu_filter 6000000 11 8   # the GPU kernel's filter on the CPU
python3 checks.py                       # the mathematics, with independent exact models (about a minute)
```

## License and credits

Passiflora's own code is **GPL-2.0-or-later**.  Passiflora was created by Alperen Yavuz; who wrote which part is in
[NOTICE.md](NOTICE.md).  The original application, GetDecics, is by Eric Driver
([drivere/get-decics-numberfields](https://github.com/drivere/get-decics-numberfields)); its code is not part of this repository.  The
binaries link PARI/GP (GPL-2.0-or-later), GMP and BOINC's libraries; see [THIRD_PARTY.md](THIRD_PARTY.md) and [NOTICE.md](NOTICE.md).
The source that corresponds to the binaries is this repository plus the upstream repository and PARI 2.8 named there.
