#!/bin/bash
# Runs an Android arm64 program on Linux, in Termux's Docker image (termux/termux-docker), whose /system holds Android's C
# library (bionic) and dynamic linker.  The program runs in the current directory:
#   droidrun.sh <program> [arguments]
#   NF_ANDROID_PROG=<program> droidrun.sh [arguments]     (as the program of build/test-suite-unix.sh, which takes one path)
# Natively on an arm64 machine; on x86-64 through qemu-user registered with binfmt (e.g. docker/setup-qemu-action).
# 2-pari-configure.sh gives it to PARI's Configure as RUNTEST, which runs its test programs with it.
# The program runs as root in the container: /system lies in Termux's app directory, which (as on Android) only its owner can
# enter.  The files it writes in the current directory therefore belong to root.
if [ -n "$NF_ANDROID_PROG" ]; then P=$NF_ANDROID_PROG; else P=$1; shift; fi
P=$(realpath "$P")
exec docker run --rm --platform linux/arm64 --security-opt seccomp=unconfined \
  -v "$PWD:$PWD" -w "$PWD" -v "$P:/droidrun/prog:ro" --entrypoint /droidrun/prog termux/termux-docker:aarch64 "$@"
