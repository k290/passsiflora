#!/bin/bash
# Android step 2: configures PARI 2.8 for Android with the NDK's clang.  PARI's Configure runs small test programs; RUNTEST
# runs them as Android programs (droidrun.sh: in Termux's Docker image, through qemu-user on x86-64).
# Env: NDK, W (work directory, with gmp-inst from step 1), PARISRC (upstream's Pari/pari-2.8-1711-ge5c317c).
set -e
NDK=${NDK:-/root/android/android-ndk-r27c}
W=${W:-/root/nfa}
PARISRC=${PARISRC:-/root/nf/upstream/Pari/pari-2.8-1711-ge5c317c}
TC=$NDK/toolchains/llvm/prebuilt/linux-x86_64/bin
rm -rf "$W/pari" && cp -r "$PARISRC" "$W/pari"
cd "$W/pari" && rm -rf O*                    # upstream commits a prebuilt x86-64 PARI
export CC=$TC/aarch64-linux-android24-clang AR=$TC/llvm-ar RANLIB=$TC/llvm-ranlib
RUNTEST="$(cd "$(dirname "$0")" && pwd)/droidrun.sh" ./Configure --static --host=aarch64-linux --with-gmp="$W/gmp-inst" \
  --without-readline > "$W/pari_conf.log" 2>&1 || true
tail -15 "$W/pari_conf.log"
ls -d O*
