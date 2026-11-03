# Building Passiflora

Passiflora is a set of files that you add to the upstream source, plus a small patch.
Nothing of upstream is copied into this repository.

The quickest route is `build/prepare-src.sh`: it copies the five upstream files that are needed, applies
`patches/passiflora.patch` and adds `src/` (set `UPSTREAM` to your checkout). The steps below do the same by hand.

## 1. Get the sources

```bash
git clone --filter=blob:none --no-checkout https://github.com/drivere/get-decics-numberfields upstream
cd upstream
git checkout 638ca4e                        # the commit Passiflora was built from
git sparse-checkout init --cone
git sparse-checkout set GetDecicsSrc Test Pari/pari-2.8-1711-ge5c317c
git checkout
patch -p1 < ../passiflora/patches/passiflora.patch   # changes TgtMartinet.cpp/.h, GetDecics.cpp and polDiscTest_cpu.cpp
cp ../passiflora/src/* GetDecicsSrc/
```

`Test/dataFiles` contains the 31 work units and their `.truth` outputs used for the tests.

## 2. Libraries

* **PARI 2.8-1711** (static): `cd Pari/pari-2.8-1711-ge5c317c && ./Configure --static --with-gmp
  --without-readline && cd Olinux-x86_64 && make lib-sta`. On a recent distribution set
  `PERL5LIB=<pari>/src/desc` first (the old build scripts expect `.` in Perl's path) and have
  `bison` installed.
* **GMP**, **BOINC's application libraries** (`libboinc.a`, `libboinc_api.a`, and
  `libboinc_opencl.a` for the GPU build) and, for the GPU build, an OpenCL library and the
  Khronos OpenCL headers.
* **gfortran** for `pdfilter.f90`.  Checked: GCC 15.2 (Ubuntu, WSL) and GCC 16.2 (MSYS2 MinGW-w64) with
  `pdfilter_Test/arith_checks.sh` (the compiled reductions against exact models, all three flag sets); GCC 15.2 also with
  `pdfilter_Test/tiny_checks.sh` (the whole filter on complete sets of polynomials at tiny primes; Linux only).  The filter relies on
  gfortran reducing an out-of-range `int(x, kind)` modulo 2^k and on `shifta` being arithmetic (`docs/whitepaper.tex`,
  Section 4, "What the language must do"); with another compiler or version run those two scripts first.  At run time
  `pdf_init` checks the compiled arithmetic on known inputs and switches the filter off (every polynomial goes to the exact
  test) if it disagrees.

## 3. Compile (Linux)

Inside `GetDecicsSrc` (with an empty `config.h` next to the sources):

```bash
# the Fortran filter, twice: AVX2 and baseline; the C++ side picks one at run time (the AVX2 one when the CPU reports
# AVX2 and FMA).  The AVX2 object is compiled with -mavx2 only (not -march=x86-64-v3, which would also allow BMI1/2,
# LZCNT, MOVBE and F16C), so the run-time check covers every instruction set extension it may use.
gfortran -O3 -march=x86-64 -mavx2 -mno-popcnt -cpp -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" \
         "-DPDF_FILTER_NAME='pdf_filter_v3'" -c pdfilter.f90 -o pdfilter_v3.o
gfortran -O3 -march=x86-64 -mtune=generic -cpp -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" \
         "-DPDF_FILTER_NAME='pdf_filter_base'" -c pdfilter.f90 -o pdfilter_base.o

INC="-I. -I$BOINC_INCLUDE -I$PARI/src/headers -I$PARI/Olinux-x86_64"
# upstream's exact test with -fwrapv: an int64 overflow in it wraps by definition (docs/TECHNICAL.md)
g++ -O3 -fwrapv $INC -DAPP_VERSION_CPU_STD -c polDiscTest_cpu.cpp
for f in TgtMartinet GetDecics; do g++ -O3 $INC -DAPP_VERSION_CPU_STD -c $f.cpp; done
g++ -O2 $INC -c nf_gpu_stub.cpp                      # CPU-only build (no OpenCL)

g++ -O3 -static -o GetDecics_cpu GetDecics.o TgtMartinet.o polDiscTest_cpu.o nf_gpu_stub.o \
    pdfilter_v3.o pdfilter_base.o libboinc_api.a libboinc.a -L$PARI/Olinux-x86_64 -l:libpari.a \
    -l:libgfortran.a -lquadmath -l:libgmp.a -lm -lpthread
```

Do **not** add `-march`, `-ffast-math` or FMA contraction to the C++ files: the enumeration
must evaluate its floating-point bounds exactly as upstream does.

CUDA kernel: `build/build-cuda-kernel.sh` compiles `src/nf_kernel.cu` with `nvcc --fmad=false` (CUDA 12.x, with a host
gcc <= 13; Ubuntu 22.04's gcc 11 works) and writes `src/nf_cuda_fatbin.h`; compile `nf_gpu_common.cpp` and
`nf_gpu_cuda.cpp` instead of the stub and link `-ldl` on Linux (nothing from CUDA is linked).  `nf_cuda_fatbin.h` and the
MUSA library (`build/build-musa.sh`) are built from `src/nf_kernel.h` and `src/nf_port_cuda.h`: rebuild them whenever
either changes (the OpenCL kernel is embedded as source by every GPU build and needs nothing extra).  The program checks at
start-up that the device code reads the prime table it builds (`NF_TABLE_FORMAT`, the kernel entry point `nf_enum_t2`): a CUDA
program with a fat binary left over from the floating-point kernel stops with a message, and a MUSA program refuses an older
`libnf_musa.so` and runs on the CPU.

OpenCL GPU build: embed the kernel (`nf_enum_cl.h` is `R"NFCL(` + `nf_port_opencl.h` + `nf_kernel.h` + `)NFCL"`), compile
`nf_gpu_common.cpp` and `nf_gpu_opencl.cpp` instead of the stub, and link `libboinc_opencl.a` and `-lOpenCL` (everything else
static). To keep the binary usable on older distributions build in an Ubuntu 22.04 environment
(glibc 2.34).

## 3b. Other systems

* **Linux aarch64:** `build/build-linux-aarch64.sh`, run natively on an Ubuntu 22.04 arm64 machine or VM.
* **macOS arm64:** `build/build-macos-arm64.sh` (Xcode tools, Homebrew gcc + gmp, BOINC's libraries, PARI built in place).
* **MUSA:** `build/build-musa.sh` inside the `mthreads/musa` container (Ubuntu 22.04); produces the program,
  `libnf_musa.so` and the vendor runtime.
* The package builder is `build/make-release.py`; tests: `build/test-suite-unix.sh` (any Unix).

## 4. Windows

Cross-compile with MinGW-w64. PARI's `Configure` must be run on a path the Windows side can
execute (it runs small test programs) with two wrappers: one that drops the `.exe` suffix from
the compiler's output and one that strips CRLF from the test programs' output (otherwise PARI
decides `sizeof(long) = 4`). The GPU build needs an import library for `OpenCL.dll`
(`dlltool -d OpenCL.def -l libOpenCL.a`).

`build/` contains the scripts used for the released binaries. They contain paths from the
author's machine and are meant as a reference, not as a turnkey build.

## 5. Test

Run each work unit of `Test/dataFiles/*.dat` (copy it to `in` in an empty directory, run the
binary, compare `out` with the `wu_*.out.truth` file ignoring the "Elapsed Time" line):
`build/test-wu.sh` and `build/test-suite.sh` do that for Linux, `build/test-gpu-windows.ps1`
and `build/test-gpu-suite-windows.ps1` for the GPU on Windows.
