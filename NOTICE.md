# Notice

**Passiflora** is an optimized build of the NumberFields@home application `GetDecics`.

* The original application — the search algorithm, its bounds, the BOINC wrapper,
  the polynomial discriminant test and the CPU/GPU test programs — was written by
  **Eric Driver** and is published at
  <https://github.com/drivere/get-decics-numberfields>. That repository does not
  state a license. **Its code is not part of this repository.** `patches/passiflora.patch`
  contains only the lines Passiflora adds to or removes from four of its files
  (`TgtMartinet.cpp`, `TgtMartinet.h`, `GetDecics.cpp`, `polDiscTest_cpu.cpp`); to build Passiflora you
  fetch the upstream source yourself (see `docs/BUILDING.md`).
* Passiflora was written by **Alperen Yavuz**: the files in `src/` and `bench/`, the
  build scripts in `build/`, the GPU version and the release packages.
* **Matthew** made the pre-filter exact and verified it:
  * the integer-only pre-filter (`src/pdfilter.f90`, `src/chain_general.inc`), which
    replaced the floating-point one: Montgomery arithmetic, a guard that keeps its
    primes inside the proved size limits, a start-up check of the compiled
    arithmetic, vectorised loads and the split into 24 small and 6 large primes;
  * the integer filter of the GPU kernel (`src/nf_kernel.h`, `src/nf_gpu_common.cpp`,
    `src/nf_port_cpu.h` and the back ends): a table-format number checked at start-up,
    the wider prime table, a start-up self-test on the device and the survivor
    comparison of `NF_GPUCHECK=1`;
  * in `patches/passiflora.patch`, the refusal of work units with more than three
    primes in S and the comments on the 64-bit steps of the exact test; in the build
    scripts, the compile flags `-mavx2` (instead of `-march=x86-64-v3`) and `-fwrapv`;
  * the proofs in `docs/whitepaper.tex` and the test suite in `pdfilter_Test/`.
* Passiflora's own code, by both, is licensed GPL-2.0-or-later (see `LICENSE`).
* Passiflora is not affiliated with or endorsed by NumberFields@home or Arizona
  State University. The released binaries contain upstream code compiled together
  with Passiflora's, linked statically with PARI/GP; they are therefore distributed
  under the terms of the GNU GPL, version 2 or later.

If you are the author of the upstream code and want something changed, please open
an issue.
