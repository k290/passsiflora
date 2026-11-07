#!/bin/bash
# Android step 1: GMP 6.3.0 (static) with the NDK's clang.  Result: $W/gmp-inst.
# Env: NDK (the Android NDK), W (work directory); the defaults are the original author's machine.
set -e
NDK=${NDK:-/root/android/android-ndk-r27c}
W=${W:-/root/nfa}
TC=$NDK/toolchains/llvm/prebuilt/linux-x86_64/bin
mkdir -p "$W" && cd "$W"
[ -f gmp-6.3.0.tar.xz ] || wget -q https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
rm -rf gmp-6.3.0 gmp-inst && tar xf gmp-6.3.0.tar.xz && cd gmp-6.3.0
export CC=$TC/aarch64-linux-android24-clang AR=$TC/llvm-ar RANLIB=$TC/llvm-ranlib STRIP=$TC/llvm-strip
./configure --host=aarch64-linux-android --prefix="$W/gmp-inst" --disable-shared --enable-static --with-pic > "$W/gmp_conf.log" 2>&1
make -j"$(nproc)" > "$W/gmp_make.log" 2>&1
make install > "$W/gmp_inst.log" 2>&1
ls -la "$W/gmp-inst/lib"
echo GMP_ANDROID_OK
