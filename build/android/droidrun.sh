#!/bin/bash
# Runs an Android arm64 program on Linux, in Termux's Docker image (termux/termux-docker), whose /system holds Android's C
# library (bionic) and dynamic linker.  The program runs in the current directory:
#   droidrun.sh <program> [arguments]
#   NF_ANDROID_PROG=<program> droidrun.sh [arguments]     (as the program of build/test-suite-unix.sh, which takes one path)
# Natively on an arm64 machine; on x86-64 through qemu-user registered with binfmt (Ubuntu: qemu-user-static).
# 2-pari-configure.sh gives it to PARI's Configure as RUNTEST, which runs its test programs with it.
if [ -n "$NF_ANDROID_PROG" ]; then P=$NF_ANDROID_PROG; else P=$1; shift; fi
P=$(realpath "$P")
exec docker run --rm --platform linux/arm64 --user "$(id -u):$(id -g)" --security-opt seccomp=unconfined \
  -v "$PWD:$PWD" -w "$PWD" -v "$P:/droidrun/prog:ro" --entrypoint /droidrun/prog termux/termux-docker:aarch64 "$@"
