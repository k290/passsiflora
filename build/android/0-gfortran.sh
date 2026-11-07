#!/bin/bash
# Android step 0: a gfortran for Android.  The NDK has no Fortran compiler, so GNU binutils and GCC (C and Fortran only) are
# built from GNU's sources as a cross compiler for aarch64-linux-android, against the NDK's sysroot.  It compiles pdfilter.f90
# and provides libgfortran.a and libgcc.a; everything else is compiled and linked with the NDK's clang (4-app.sh).
#   NDK=<android-ndk-r27c> PREFIX=<install dir> bash 0-gfortran.sh        (WORK: build directory, JOBS: parallel jobs)
# Result: $PREFIX/bin/aarch64-linux-android-gfortran.  GCC 15.2 is one of the two versions the filter's own tests ran on.
set -e
GCC_V=15.2.0
BINUTILS_V=2.45
API=24                                   # the NDK's start files and libraries of this API level (as 4-app.sh's clang)
T=aarch64-linux-android
NDK=${NDK:?set NDK to the Android NDK}
PREFIX=${PREFIX:?set PREFIX to the install directory}
WORK=${WORK:-$PWD/gfortran-build}
JOBS=${JOBS:-$(nproc)}
export PATH=$PREFIX/bin:$PATH
mkdir -p "$WORK" "$PREFIX"
cd "$WORK"
[ -f gcc-$GCC_V.tar.xz ] || wget -q https://ftp.gnu.org/gnu/gcc/gcc-$GCC_V/gcc-$GCC_V.tar.xz
[ -f binutils-$BINUTILS_V.tar.xz ] || wget -q https://ftp.gnu.org/gnu/binutils/binutils-$BINUTILS_V.tar.xz

# The sysroot: the NDK's headers, with the architecture's headers and the API level's start files and stub libraries in
# usr/include and usr/lib, where GCC looks for them.
N=$NDK/toolchains/llvm/prebuilt/linux-x86_64/sysroot/usr
S=$PREFIX/sysroot
rm -rf "$S" && mkdir -p "$S/usr/lib"
cp -r "$N/include" "$S/usr/" && cp -r "$N/include/$T/." "$S/usr/include/"
cp -r "$N/lib/$T/$API/." "$S/usr/lib/"

# binutils: the assembler and archiver.  The linker is the NDK's lld instead: binutils' ld cannot read the NDK's libraries
# (their debug sections are zstd-compressed), and lld knows Android.
rm -rf binutils-$BINUTILS_V b-binutils && tar xf binutils-$BINUTILS_V.tar.xz && mkdir b-binutils && cd b-binutils
../binutils-$BINUTILS_V/configure --target=$T --prefix="$PREFIX" --with-sysroot="$S" --disable-nls --disable-werror \
  --disable-gdb --disable-gprofng --disable-sim > ../binutils.log 2>&1
make -j"$JOBS" >> ../binutils.log 2>&1 && make install >> ../binutils.log 2>&1
for f in "$PREFIX/$T/bin/ld" "$PREFIX/$T/bin/ld.bfd" "$PREFIX/bin/$T-ld" "$PREFIX/bin/$T-ld.bfd"; do
  ln -sf "$NDK/toolchains/llvm/prebuilt/linux-x86_64/bin/ld.lld" "$f"
done
cd ..

rm -rf gcc-$GCC_V b-gcc && tar xf gcc-$GCC_V.tar.xz
( cd gcc-$GCC_V && ./contrib/download_prerequisites > ../prereq.log 2>&1 )   # GMP, MPFR, MPC and ISL for the compiler itself
mkdir b-gcc && cd b-gcc
# Only what the Fortran runtime needs; static libraries only (the program links them statically), position-independent
# code as Android requires.  The NDK's headers are written for clang: __ANDROID_VENDOR__ turns their availability
# attributes off (the API level is fixed by the libraries linked), BIONIC_IOCTL_NO_SIGNEDNESS_OVERLOAD a clang-only overload
# of ioctl, and the nullability keywords are defined away, so GCC can compile the runtime against them.
../gcc-$GCC_V/configure --target=$T --prefix="$PREFIX" --with-sysroot="$S" --enable-languages=c,fortran --enable-default-pie \
  --disable-lto --disable-bootstrap --disable-multilib --disable-nls --disable-shared --disable-werror --disable-libssp \
  --disable-libgomp --disable-libquadmath --disable-libsanitizer --disable-libatomic --disable-libitm --disable-libvtv \
  --disable-libstdcxx \
  CFLAGS_FOR_TARGET="-O2 -D__ANDROID_API__=$API -D__ANDROID_VENDOR__ -D_Nonnull= -D_Nullable= -D_Null_unspecified= -DBIONIC_IOCTL_NO_SIGNEDNESS_OVERLOAD" > ../gcc.log 2>&1
make -j"$JOBS" all-gcc >> ../gcc.log 2>&1 || { grep -n -m20 "error:" ../gcc.log; exit 1; }
# GCC knows Android only for 32-bit ARM.  In the specs that the build's compiler reads (and that are installed with it): bionic's
# start files instead of glibc's crt1.o, crti.o and crtn.o, and the 64-bit dynamic linker /system/bin/linker64.
awk -v s='%{shared:crtbegin_so.o%s;static:crtbegin_static.o%s;:crtbegin_dynamic.o%s}' \
    -v e='%{shared:crtend_so.o%s;:crtend_android.o%s}' \
    '{print} /^\*startfile:$/{getline; print s} /^\*endfile:$/{getline; print e}' gcc/specs \
  | sed 's#/system/bin/linker\([^6]\)#/system/bin/linker64\1#g' > gcc/specs.android
mv gcc/specs.android gcc/specs
make -j"$JOBS" >> ../gcc.log 2>&1 || { grep -n -m20 "error:" ../gcc.log; exit 1; }
make install >> ../gcc.log 2>&1
cp gcc/specs "$PREFIX/lib/gcc/$T/$GCC_V/specs"
"$PREFIX/bin/$T-gfortran" --version | head -1
ls -la "$("$PREFIX/bin/$T-gfortran" -print-file-name=libgfortran.a)"
echo GFORTRAN_ANDROID_OK
