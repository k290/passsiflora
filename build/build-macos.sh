#!/bin/bash
# macOS CPU build for the architecture it runs on: arm64 (Apple silicon) or x86_64 (Intel).  CPU only: Apple GPUs have no
# double precision.  The release builds each architecture natively on its own GitHub runner, tests it there and joins the
# two programs with lipo into one universal program (.github/workflows/release.yml, job macos-universal).
# Needs: Xcode command line tools, Homebrew gcc (gfortran) and gmp, BOINC's application libraries
# (libboinc.a, libboinc_api.a + headers; built from source, see docs/BUILDING.md), PARI 2.8-1711 built in $W/pari
# (./Configure --static --with-gmp=<Homebrew prefix> --without-readline; make lib-sta) and the source tree in $W/src
# (output of prepare-src.sh).  Build BOINC and PARI with MACOSX_DEPLOYMENT_TARGET=11.0, the program's minimum macOS.
# Floating-point contraction is off so the enumeration rounds like the x86 build.
set -e
ARCH=$(uname -m)
if [ "$ARCH" = arm64 ]; then HB=/opt/homebrew; else HB=/usr/local; fi   # Homebrew's prefix on each architecture
export PATH=$HB/bin:$PATH
W=${W:-$HOME/passiflora-work}
BOINC=${BOINC:-$HOME/juniper-work/boinc-install}
P=$W/pari
OBJ=$(basename "$(ls -d $P/O* | head -1)")
MIN=-mmacosx-version-min=11.0
cd "$W/src"
: > config.h
INC="-I. -I$BOINC/include/boinc -I$BOINC/include -I$P/src/headers -I$P/$OBJ -I$HB/include"
rm -f *.o *.mod GetDecics_*
if [ "$ARCH" = arm64 ]; then
  gfortran -O3 -mcpu=apple-m1 -ffp-contract=off -cpp $MIN -c pdfilter.f90 -o pdfilter.o
  PDF=pdfilter.o
else
  # x86-64: an AVX2 and a baseline filter, picked at run time as on Linux and Windows (NF_FORCE_BASE=1 forces the baseline)
  gfortran -O3 -march=x86-64 -mavx2 -mno-popcnt -cpp $MIN -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" "-DPDF_FILTER_NAME='pdf_filter_v3'" -c pdfilter.f90 -o pdfilter_v3.o
  gfortran -O3 -march=x86-64 -mtune=generic -cpp $MIN -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" "-DPDF_FILTER_NAME='pdf_filter_base'" -c pdfilter.f90 -o pdfilter_base.o
  PDF="pdfilter_v3.o pdfilter_base.o"
fi
# upstream's exact test with -fwrapv: an int64 overflow in it wraps by definition (docs/TECHNICAL.md)
clang++ -O3 -fwrapv -ffp-contract=off -std=c++11 -w $MIN $INC -DAPP_VERSION_CPU_STD -c -o polDiscTest_cpu.o polDiscTest_cpu.cpp
for f in TgtMartinet GetDecics; do
  clang++ -O3 -ffp-contract=off -std=c++11 -w $MIN $INC -DAPP_VERSION_CPU_STD -c -o $f.o $f.cpp
done
clang++ -O2 -std=c++11 -w $MIN $INC -c nf_gpu_stub.cpp
GF=$(gfortran -print-file-name=libgfortran.a)
QM=$(gfortran -print-file-name=libquadmath.a)
GCC=$(gfortran -print-libgcc-file-name)
clang++ -O3 $MIN -o GetDecics_cpu GetDecics.o TgtMartinet.o polDiscTest_cpu.o nf_gpu_stub.o $PDF \
  $BOINC/lib/libboinc_api.a $BOINC/lib/libboinc.a $P/$OBJ/libpari.a $HB/lib/libgmp.a \
  "$GF" "$QM" "$GCC" -lm -lc++
strip GetDecics_cpu
ls -la GetDecics_cpu
lipo -archs GetDecics_cpu
otool -L GetDecics_cpu
echo PK_MACOS_OK
