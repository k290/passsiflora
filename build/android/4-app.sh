#!/bin/bash
# Android step 4: the CPU program.  pdfilter.f90 is compiled with the gfortran of step 0, everything else with the NDK's
# clang, linked against libgfortran.a and libgcc.a of step 0 and the TLS padding of android_tls_pad.c.  API level 24; the
# native TLS of the Fortran runtime needs Android 10 or newer at run time.
# Env: NDK, GT (bin directory of step 0's toolchain), BOINC (libraries in BOINC/lib, headers in BOINC/include/boinc:
# build-boinc-libs.sh with HOST and PREFIX), W (work directory of steps 1-3), SRC (prepare-src.sh's GetDecicsSrc).
# Result: $W/src/GetDecics_cpu.
set -e
NDK=${NDK:-/root/android/android-ndk-r27c}
GT=${GT:-/root/android/toolchain-aarch64-linux-android/bin}
BOINC=${BOINC:-/root/android/boinc-android}
W=${W:-/root/nfa}
SRC=${SRC:-/root/nf/build_src/GetDecicsSrc}
TC=$NDK/toolchains/llvm/prebuilt/linux-x86_64/bin
P=$W/pari
G=$W/gmp-inst
rm -rf "$W/src" && mkdir -p "$W/src" && cp "$SRC"/* "$W/src/" && cd "$W/src"
: > config.h
CXX=$TC/aarch64-linux-android24-clang++
CC=$TC/aarch64-linux-android24-clang
INC="-I. -I$BOINC/include/boinc -I$P/src/headers -I$P/Olinux-aarch64 -I$G/include"
$GT/aarch64-linux-android-gfortran -O3 -ffp-contract=off -cpp -c pdfilter.f90 -o pdfilter.o
# upstream's exact test with -fwrapv: an int64 overflow in it wraps by definition (docs/TECHNICAL.md)
$CXX -O3 -fwrapv -ffp-contract=off -std=c++11 -w $INC -DAPP_VERSION_CPU_STD -c -o polDiscTest_cpu.o polDiscTest_cpu.cpp
for f in TgtMartinet GetDecics; do
  $CXX -O3 -ffp-contract=off -std=c++11 -w $INC -DAPP_VERSION_CPU_STD -c -o $f.o $f.cpp
done
$CXX -O2 -std=c++11 -w $INC -c nf_gpu_stub.cpp
$CC -O2 -fno-emulated-tls -c "$(cd "$(dirname "$0")" && pwd)/android_tls_pad.c" -o tlspad.o
GF=$($GT/aarch64-linux-android-gfortran -print-file-name=libgfortran.a)
GC=$($GT/aarch64-linux-android-gfortran -print-libgcc-file-name)
$CXX -O3 -static-libstdc++ -o GetDecics_cpu GetDecics.o TgtMartinet.o polDiscTest_cpu.o nf_gpu_stub.o pdfilter.o tlspad.o \
  "$BOINC/lib/libboinc_api.a" "$BOINC/lib/libboinc.a" "$P/Olinux-aarch64/libpari.a" "$G/lib/libgmp.a" "$GF" "$GC" -lm -ldl
$TC/llvm-strip GetDecics_cpu
ls -la GetDecics_cpu
$TC/llvm-readelf -l GetDecics_cpu | grep -E "TLS|interp"
echo ANDROID_APP_OK
