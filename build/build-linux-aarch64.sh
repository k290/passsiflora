#!/bin/bash
# Linux aarch64 build, run natively on an Ubuntu 20.04 arm64 machine, VM or container with GCC 11 (the release workflow's
# setup: the OpenCL program then needs glibc 2.31+; on a newer system set MAX_GLIBC, see below).
# Needs: build-essential gfortran libgmp-dev libboinc-app-dev bison perl python3 (+ ocl-icd-opencl-dev opencl-headers for OpenCL)
# Expects $HOME/pk/armpack.tar.gz made from the output of prepare-src.sh (src/), PARI 2.8-1711 (pari/) and Test/.
# Floating-point contraction is switched off: the enumeration must round exactly like the x86 build.
set -e
W=$HOME/pk
cd "$W"
[ -d src ] || tar xzf armpack.tar.gz
if [ ! -f pari/Olinux-arm/libpari.a ] && ! ls pari/O*/libpari.a >/dev/null 2>&1; then
  ( cd pari && ./Configure --static --with-gmp --without-readline > "$W/pari_conf.log" 2>&1 )
  export PERL5LIB="$W/pari/src/desc"
  ( cd pari/O* && make -j4 lib-sta > "$W/pari_make.log" 2>&1 )
fi
OBJ=$(basename "$(ls -d $W/pari/O* | head -1)")
P=$W/pari
cd "$W/src"
: > config.h
python3 - <<'PY'
s = open('nf_port_opencl.h').read() + open('nf_kernel.h').read()
nl = chr(10)
open('nf_enum_cl.h', 'w').write('R"NFCL(' + nl + s + ')NFCL"' + nl)
PY
INC="-I. -I/usr/include/boinc -I$P/src/headers -I$P/$OBJ"
BL=/usr/lib/aarch64-linux-gnu
rm -f *.o *.mod GetDecics_*
gfortran -O3 -ffp-contract=off -cpp -c pdfilter.f90 -o pdfilter.o
# upstream's exact test with -fwrapv: an int64 overflow in it wraps by definition (docs/TECHNICAL.md)
g++ -O3 -fwrapv -ffp-contract=off -w $INC -DAPP_VERSION_CPU_STD -c -o polDiscTest_cpu.o polDiscTest_cpu.cpp
for f in TgtMartinet GetDecics; do
  g++ -O3 -ffp-contract=off -w $INC -DAPP_VERSION_CPU_STD -c -o $f.o $f.cpp
done
g++ -O2 -std=c++11 -w $INC -c nf_gpu_stub.cpp
g++ -O2 -std=c++11 -w $INC -c nf_gpu_common.cpp
g++ -O2 -std=c++11 -w $INC -c nf_gpu_opencl.cpp
# CPU: fully static (the pthread symbols are forced in: static libgfortran otherwise crashes at exit with glibc 2.35 and newer)
U=""
for s in pthread_once pthread_getspecific pthread_setspecific pthread_create pthread_join pthread_equal pthread_self pthread_detach pthread_cancel pthread_mutex_lock pthread_mutex_trylock pthread_mutex_unlock pthread_mutex_init pthread_mutex_destroy pthread_cond_init pthread_cond_broadcast pthread_cond_signal pthread_cond_wait pthread_cond_destroy pthread_mutexattr_init pthread_mutexattr_settype pthread_mutexattr_destroy pthread_key_create pthread_key_delete; do U="$U -Wl,-u,$s"; done
g++ -O3 -static $U -o GetDecics_cpu GetDecics.o TgtMartinet.o polDiscTest_cpu.o nf_gpu_stub.o pdfilter.o \
  $BL/libboinc_api.a $BL/libboinc.a -L$P/$OBJ -l:libpari.a -l:libgfortran.a -l:libgmp.a -lm -lpthread
# OpenCL: loaded from the system, everything else static
g++ -O3 -o GetDecics_gpu GetDecics.o TgtMartinet.o polDiscTest_cpu.o nf_gpu_common.o nf_gpu_opencl.o pdfilter.o \
  $BL/libboinc_opencl.a $BL/libboinc_api.a $BL/libboinc.a -L$P/$OBJ -l:libpari.a \
  -static-libgcc -static-libstdc++ -l:libgfortran.a -l:libgmp.a -lOpenCL -lm -lpthread
strip GetDecics_cpu GetDecics_gpu
ls -la GetDecics_cpu GetDecics_gpu
# The OpenCL program needs at least the glibc it was linked against: 2.31 when built on Ubuntu 20.04 (build-linux-chroot-inner.sh).
MAX_GLIBC=${MAX_GLIBC:-2.31}
g=$(objdump -T GetDecics_gpu | grep -o 'GLIBC_[0-9][0-9.]*' | sed 's/GLIBC_//' | sort -uV | tail -1)
echo "GetDecics_gpu needs glibc $g"
if [ "$(printf '%s\n%s\n' "$g" "$MAX_GLIBC" | sort -V | tail -1)" != "$MAX_GLIBC" ]; then
  echo "GetDecics_gpu needs glibc $g, newer than $MAX_GLIBC: build it on an older system"; exit 1
fi
# Only the C library and the OpenCL loader may be shared: any other library stops every task with exit code 127 on systems
# without it (build-linux-chroot-inner.sh).
libs=$(readelf -d GetDecics_gpu | sed -n 's/.*Shared library: \[\(.*\)\]/\1/p' | tr '\n' ' ')
echo "GetDecics_gpu loads: $libs"
extra=$(echo $libs | tr ' ' '\n' | grep -v -x -E 'libc\.so\.6|libm\.so\.6|libpthread\.so\.0|libdl\.so\.2|ld-linux-aarch64\.so\.1|libOpenCL\.so\.1' || true)
if [ -n "$extra" ]; then echo "GetDecics_gpu needs libraries a volunteer may not have: $extra"; exit 1; fi
echo PK_AARCH64_OK
