# How Passiflora differs from GetDecics

Passiflora is applied on top of Eric Driver's GetDecics ([drivere/get-decics-numberfields](https://github.com/drivere/get-decics-numberfields),
commit `638ca4e`). This document lists what Passiflora changes in how the programs compute, are built, are packaged and run, why
each change was made, and how it fits BOINC's conventions. Each entry names the release that introduced it. The mathematics and its proofs are in [whitepaper.pdf](whitepaper.pdf); [TECHNICAL.md](TECHNICAL.md)
is the technical overview.

BOINC's conventions here come from BOINC's documentation (the wiki pages AppCoprocessor, OpenCL-Applications,
API-for-native-apps, GpuSched, Client-configuration and Anonymous-platform) and from its source: the client (`client/app_start.cpp`)
gives each GPU task its device in `init_data.xml` and as `--device N`, and the project scheduler (`sched/sched_types.cpp`,
`sched/sched_send.cpp`) ignores coprocessor types it does not know.

## At a glance

| Change | Since | Upstream | Fits BOINC? |
|---|---|---|---|
| CPU pre-filter before the exact test | v0.1 (exact integer filter: v0.4) | exact test on every polynomial | yes: same output, shorter tasks |
| Exact shortcuts in the enumeration loops | v0.1 | `sqrt` comparisons, full m1 range | yes |
| The whole search on the GPU | v0.1 | GPU runs only the discriminant test | yes |
| GPU self-checks that stop a miscomputing GPU | v0.3, v0.4 | none | yes: an error, not a wrong result |
| One CPU program choosing AVX2 or baseline code at run time | v0.1 | one generic build | yes |
| Static programs; GPU drivers loaded at run time | v0.1, v0.2 | partly static | yes |
| One OpenCL program for NVIDIA, AMD and Intel | v0.1 | one program per vendor | yes (BOINC supports it) |
| Kernel inside the program, batch size adapts | v0.1 | kernel file and per-card table shipped next to it | yes |
| CUDA through the driver API, kernel for every generation embedded | v0.2 | CUDA runtime, compiled for listed architectures | yes |
| More platforms (ARM, macOS, Android, FreeBSD) | v0.3, v0.5.1, v0.5.2 | Windows, Linux (ARM64 in beta), Intel Mac | yes: what anonymous platform is for |
| Release workflow: tag → every package, tested | v0.4 | built by hand | — |
| Packages with CPU and GPU app versions in one `app_info.xml` | v0.1 | served by the project | yes |
| 4 GPU tasks per card | v0.1 | the project's setting | yes |
| A GPU task's CPU share: half a core instead of a whole one | v0.6.0 | 0.1 CPU in upstream's example `init_data.xml` | yes, measured |
| The CPU sleeps while the GPU computes | v0.6.0 (upstream did this for CUDA) | CUDA: blocking sync | yes |
| Device choice by BOINC's rules, also when run by hand | v0.6.0 | BOINC's rules | yes |
| A GPU that cannot be used stops the task (temporary problems: retry) | v0.6.0 (v0.1–v0.5.2 ran on the CPU) | stops the task | yes; the CPU fallback did not |
| MUSA: a CPU app version that uses one card per task | v0.3, one task per card: v0.6.0 | — | as far as BOINC allows |

## Upstream in brief

NumberFields serves GetDecics as one program per resource: a CPU program (version 4.00) and GPU programs (version 4.02) for CUDA
(plan class `cuda30`) and for OpenCL on NVIDIA, AMD and Intel (`opencl_nvidia`, `opencl_amd`, `opencl_intel_gpu`), each with
its vendor compiled in. The CPU program enumerates candidate polynomials and runs an exact GMP discriminant test on each; the test
takes about 98 % of the time and rejects all but about one polynomial in 10^7. The GPU programs keep the enumeration on the CPU and
send batches of polynomials to the GPU, which runs the discriminant test. The OpenCL kernel (`pdtKernel.cl` and its headers) and a
per-card tuning table (`gpuLookupTable.txt`) are files of the app version. A GPU program takes its device from BOINC
(`init_data.xml`, or `--device N` from old clients); without one it stops with an error. To run one by hand, upstream's test
instructions copy an example `init_data.xml` into the directory (`Test/exampleInitDataXmlFiles`). GPU work runs in BOINC critical
sections; every error ends the task.

## 1. What the programs compute

**The pre-filter** (v0.1; exact integer filter v0.4; one prime size v0.5). Before the exact test, `src/pdfilter.f90` computes the
Legendre symbol of the polynomial's discriminant at 30 primes. A symbol of −1 proves that the polynomial fails, so the filter
rejects it; every other polynomial goes to the unchanged exact test. *Why:* the exact test dominates the run time and fails almost
always; a cheap proof of failure avoids most of it, and since the filter only rejects with a proof, the output is the stock
program's. v0.1 to v0.3 computed the symbols in floating point; v0.4 replaced that with integer arithmetic, which needs no
assumption about rounding modes or compiler flags and comes with proofs (whitepaper). v0.5 dropped a second, slower stage of larger
primes that only the floating-point filter had needed.

**Exact shortcuts in the enumeration** (v0.1). `sqrt(X) <= B` becomes `0 <= X <= T` with a precomputed T, and the innermost (m1)
loop skips the part of its range where the a4 conditions provably fail. *Why:* speed; both are exact (TECHNICAL.md, section 2).

**The whole search on the GPU** (v0.1). The kernel runs the innermost loops and an integer version of the filter; the CPU walks
the outer loops and finishes the roughly one polynomial in 10^5 that survives, with the same exact tests as the CPU program.
*Why:* upstream's GPU programs moved 88 bytes per polynomial over PCIe (about 50 GB per work unit) from one enumerating CPU
thread, so a modern card waited most of the time. The output is unchanged: the kernel follows the CPU's floating-point order and
integer wrap-around, and the survivors are sorted back into enumeration order.

**GPU self-checks** (cross-check at the start of every task: v0.3; start-up self-test: v0.4). The self-test compares the kernel's
filter on fixed inputs with values computed independently on the host; during the first 5 million polynomials of every task the
CPU runs the same loops and the counts must agree. *Why:* each machine's driver compiles the OpenCL kernel, so the code that runs
was never tested by the project. *BOINC:* a GPU that computes wrong values must produce an error, not a result: both checks stop the
task with exit status 6 (the self-test does so from v0.6.0 on; before, a failed self-test ran the task on the CPU).

## 2. How the programs are built

**One CPU program per platform that picks AVX2 or baseline code at run time** (v0.1). `pdfilter.f90` is compiled twice on x86-64;
the program uses the AVX2 copy when the processor has AVX2 and FMA. *Why:* one file for anonymous-platform users, who would
otherwise have to choose; it cannot run the wrong code. A project would use CPU plan classes for the same effect.

**Static programs** (v0.1, v0.2). The Linux CPU program is fully static; the GPU programs are built on Ubuntu 20.04 (glibc 2.31) with
everything static except the C library and the driver, which is loaded at run time. *Why:* the programs run on any distribution without installing
libraries. From v0.2 to v0.6.0 the GPU programs also loaded `libquadmath.so.0` by mistake (`-lquadmath` without `-static` picks the
shared library), and on systems without gfortran's runtime their tasks stopped with exit code 127 before starting. Since v0.6.1
it is static too, the build scripts refuse any shared library other than the C library and the OpenCL loader, and the release
workflow starts every Linux x86-64 program in a bare Debian 13 container. Up to v0.6.1 the GPU programs were built on Ubuntu
22.04 and needed glibc 2.34: on Ubuntu 20.04, Debian 11 and other systems with an older C library their tasks stopped at start
("version \`GLIBC_2.34' not found"). After v0.6.1 they are built on Ubuntu 20.04 with the same compiler (GCC 11, from Ubuntu's
toolchain PPA) and need glibc 2.31; the build scripts refuse a program that needs a newer one, and the release workflow also
starts each program in a bare Ubuntu 20.04 container (MUSA excepted: Moore Threads' SDK needs Ubuntu 22.04).

**One OpenCL program for NVIDIA, AMD and Intel** (v0.1). The vendor comes from BOINC's `gpu_type`; `app_info.xml` lists the same
program under the three OpenCL plan classes. *Why:* one file to build and test. BOINC supports this explicitly
(`boinc_get_opencl_ids` with type 0 takes the vendor from `init_data.xml`).

**The kernel inside the program; the batch size adapts** (v0.1). The OpenCL source is compiled into the program, and the number of
work items per launch adapts to about 100 ms. *Why:* no kernel or tuning files to keep in step with the program, no table of cards
to maintain, and launches short enough for display drivers' watchdogs.

**CUDA through the driver API** (v0.2). The CUDA program loads the driver (`nvcuda.dll`, `libcuda.so.1`) at run time and carries a
fat binary of the kernel (cubins from sm_50 to sm_120 plus PTX). *Why:* one program for every NVIDIA generation, new ones included
through PTX; nothing of CUDA is linked. Upstream's nvcc build links the CUDA runtime statically, so it too needs only the driver.

**More platforms** (v0.3: Linux ARM64, macOS on Apple silicon, RISC-V, PowerPC, IBM Z, armhf, Android, FreeBSD; v0.5.1: one
universal macOS program for Intel and Apple silicon, built for macOS 11; v0.5.2: Android built by the release workflow with a gfortran
for Android built from GNU's sources; after v0.6.0: macOS built for macOS 15, the oldest macOS the release workflow tests, and no more RISC-V, PowerPC or IBM Z
packages). *Why:* anonymous platform exists for platforms a project does not serve.

**The release workflow** (v0.4). Pushing a version tag builds every package on GitHub's runners, runs the CPU programs on the 31
upstream work units and makes a draft release.

## 3. Packaging and running under BOINC

**Anonymous-platform packages** (v0.1, CUDA v0.2). `CPU_`, `GPU_`, `CPU-GPU_`, `GPU-CUDA_` and `CPU-GPU-CUDA_`; `app_info.xml` uses the
project's app name, version numbers (400 for the CPU, 402 for the GPU) and plan-class names. *Why:* under anonymous platform BOINC
runs only the app versions listed, so a package with only the GPU program would leave the CPU idle for the project; the combined
packages let BOINC give each task to the CPU or the GPU, as it does with the project's own versions.

**4 tasks per card** (v0.1, `app_config.xml`: `gpu_usage` 0.25). *Why:* a task's CPU phases (the outer loops, the exact tests of the
survivors) leave the card idle; with four tasks the others fill those gaps. Measured on a GTX 1650: four tasks finished a work unit
each in 83 s, one task alone in 29.7 s, 43 % more work per hour.

**A GPU task's CPU share** (v0.6.0; was a whole core since v0.1). `app_config.xml` sets `cpu_usage` to 0.5 and `app_info.xml`
`avg_ncpus` to 0.5 for the GPU versions. *Why:* BOINC sets this much CPU aside for each GPU task and runs fewer CPU tasks
accordingly; BOINC also overcommits the CPU on purpose to keep GPUs busy (GpuSched), so the number should be what a task uses.
Measured on a GTX 1650 with an i5-12500H, once the task sleeps while it waits (next entry): 0.3 to 0.6 of a core, depending on
the work unit and on how many tasks share the card. A faster card waits less for its CPU work, so a task uses more there;
`app_config.xml` explains how to read the right value from a task's CPU time and run time.

**The CPU sleeps while the GPU computes** (v0.6.0). The CUDA program asks the driver to block (`CU_CTX_SCHED_BLOCKING_SYNC`), as
upstream's CUDA program did (`cudaDeviceScheduleBlockingSync`), and the MUSA library does the same. NVIDIA's OpenCL has no such
mode, so the OpenCL program polls the kernel's event every millisecond (on Windows with a high-resolution timer). *Why:* both
NVIDIA drivers by default spin a CPU core for as long as a kernel runs, so every GPU task took a whole core: 0.91 to 0.96 of a core
in the measurements, against 0.3 to 0.6 now, at the same speed (one task: 27.5 s of CPU time for 30.1 s before, 14.1 s for 29.7 s
now with OpenCL; four tasks: 83 s instead of 86 s).

**Which GPU a task uses** (v0.6.0; before: BOINC's assignment only, or `NF_OCL_ANY=1`/`NF_CUDA_ANY=1`). Under BOINC the GPU is the
one the client assigned (`init_data.xml`: `gpu_type`, `gpu_device_num`). Run by hand it is the one an `init_data.xml` in the
directory names (upstream's way of testing), else the one `--device N` picks, else the first: the OpenCL program lists the OpenCL
GPUs of all platforms on stderr, numbered for `--device`. *Why:* these are BOINC's own conventions, so the program behaves like
upstream's and like BOINC's samples, and running it by hand needs no environment variables.

**When the GPU cannot be used** (v0.6.0; v0.1 to v0.5.2 ran the task on the CPU). The GPU programs now use the GPU or stop the task with
a message that names the cause: no GPU assigned (the program was listed without a `<coproc>`), a GPU type the program does not
use, no driver, no such device, no double precision, a kernel the driver cannot compile, or a failed self-test (exit status 6).
Problems that pass by themselves (out of GPU or host memory, the device in use by another program) end the task with
`boinc_temporary_exit`, so BOINC restarts it a minute later from its checkpoint, as BOINC recommends for GPU memory allocation
failures. Run by hand, every problem ends the program with an error. `NF_GPU=0` still runs the GPU program on the CPU, as a
developer switch. *Why the fallback had to go:* BOINC gives a GPU task a share of a GPU and a fraction of a CPU core. A task that
quietly ran on the CPU left its GPU share idle, ran about ten times longer (v0.2's figures: 338 s for the 31 test units on an RTX
5070 Ti, 3,155 s on one core), skewed BOINC's run-time estimates for the GPU version and so its work fetch, and hid a broken
driver or card behind slow successful tasks. It also did the CPU version's job outside BOINC's CPU accounting: CPU work belongs to
the CPU app version, which the combined packages provide. (The fallbacks came in v0.1: no GPU assigned, too few filter primes; v0.2:
no CUDA driver; v0.4: failed self-test. The too-few-primes case cannot happen for sets S of up to three primes, which is all the
program accepts.)

**MUSA** (v0.3; one task per card, `--device` and blocking sync: v0.6.0). BOINC does not detect Moore Threads cards. Its documented
way for such a card is a `<coproc>` in the client's `cc_config.xml` with the same type in `app_info.xml`, but a project's scheduler
ignores coprocessor types it does not know: it counts such an app version as a CPU version and skips the client's requests for
that coprocessor's work, so the card would never get tasks. The MUSA program therefore stays what it was: to BOINC a CPU app version,
which also uses the card. What changes is that only one task at a time uses each card (an exclusive lock on a file next to the
program, released when the process ends); the other tasks run on the CPU, which is what BOINC set aside a core for. Without a
usable card, or when its self-test fails, every task runs on the CPU.

**Critical sections** (v0.1, as upstream). Every kernel launch and the reads of its results run in a BOINC critical section, so the
client does not kill a task in the middle of a GPU call.

## 4. Running by hand

Copy a work unit to `in` in an empty directory and start the program there; it writes `out` and its log to `stderr.txt`. A GPU
program uses the GPU an `init_data.xml` in the directory names (as BOINC writes it, or upstream's examples in
`Test/exampleInitDataXmlFiles`), else `--device N`, else the first GPU, and says which on its first lines. `NF_GPU=0` runs a GPU program
on the CPU. The other environment variables are test switches for developers (TECHNICAL.md); none is needed to run the programs.
