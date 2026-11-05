#!/bin/bash
# Builds the BOINC application libraries (libboinc_api.a, libboinc.a) from source, for architectures whose
# distribution packages ship empty archives (Ubuntu 22.04 riscv64), or for macOS.  Result: $HOME/pk/boinclibs/*.a, or with
# PREFIX=dir the libraries in dir/lib and the headers in dir/include/boinc.
# MINGW=1: for Windows instead, with BOINC's lib/Makefile.mingw (MinGW-w64 cross compilers on Linux, or MSYS2 MINGW64).
# Result: $B/lib/libboinc.a and $B/lib/libboinc_api.a, B the BOINC source (the B of build-windows-gpu.sh); AR, RANLIB as needed.
set -e
if [ -n "$MINGW" ]; then
  B=${B:-/root/boincsrc/boinc}
  rm -rf "$B" && git clone -q --depth 1 --branch client_release/7.20/7.20.7 https://github.com/BOINC/boinc.git "$B"
  cd "$B"
  # version.h is otherwise made by configure, which does not run for MinGW
  sed -e 's/@BOINC_MAJOR_VERSION@/7/; s/@BOINC_MINOR_VERSION@/20/; s/@BOINC_RELEASE@/7/; s/@WRAPPER_RELEASE@/0/; s/@VBOXWRAPPER_RELEASE@/0/' \
      -e 's/@BOINC_VERSION_STRING@/7.20.7/; s/@PACKAGE@/boinc/; s/@PACKAGE_BUGREPORT@//; s/@PACKAGE_NAME@/BOINC/' \
      -e 's/@PACKAGE_STRING@/BOINC 7.20.7/; s/@PACKAGE_TARNAME@/boinc/; s/@PACKAGE_VERSION@/7.20.7/' version.h.in > version.h
  sh ./generate_svn_version.sh
  make -C lib -f Makefile.mingw -j4 MINGW=x86_64-w64-mingw32 AR=${AR:-x86_64-w64-mingw32-ar} RANLIB=${RANLIB:-x86_64-w64-mingw32-ranlib} \
       libboinc.a libboinc_api.a
  ls -la lib/libboinc.a lib/libboinc_api.a
  echo BOINCLIBS_MINGW_OK
  exit 0
fi
if command -v apt-get >/dev/null; then
  apt-get install -y -qq git autoconf automake libtool m4 pkg-config make g++ libssl-dev libcurl4-openssl-dev >/dev/null
fi
cd /tmp
rm -rf boinc && git clone -q --depth 1 --branch client_release/7.20/7.20.7 https://github.com/BOINC/boinc.git || git clone -q --depth 1 https://github.com/BOINC/boinc.git
cd boinc
./_autosetup >/dev/null 2>&1
./configure --enable-libraries --disable-server --disable-client --disable-manager --disable-shared --enable-static --disable-fcgi >/dev/null 2>&1
make -j3 >/dev/null 2>&1 || true
if [ -n "$PREFIX" ]; then   # PREFIX=dir: libraries in dir/lib, headers in dir/include/boinc (the BOINC of build-macos-arm64.sh)
  mkdir -p "$PREFIX/lib" "$PREFIX/include/boinc"
  cp api/.libs/libboinc_api.a lib/.libs/libboinc.a "$PREFIX/lib/" 2>/dev/null || cp api/libboinc_api.a lib/libboinc.a "$PREFIX/lib/"
  cp lib/*.h api/*.h version.h config.h "$PREFIX/include/boinc/"
  ls -la "$PREFIX/lib"
  echo BOINCLIBS_OK
  exit 0
fi
mkdir -p $HOME/pk/boinclibs
cp api/.libs/libboinc_api.a lib/.libs/libboinc.a $HOME/pk/boinclibs/ 2>/dev/null || cp api/libboinc_api.a lib/libboinc.a $HOME/pk/boinclibs/
ls -la $HOME/pk/boinclibs
echo BOINCLIBS_OK
