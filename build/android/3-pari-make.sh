#!/bin/bash
# Android step 3: builds libpari.a (env W as in step 2).
set -e
W=${W:-/root/nfa}
cd "$W/pari/Olinux-aarch64"
grep -E "sizeof_long|doubleformat|kernlvl|^CC=|gmp|GMP" pari.cfg | head
grep -E "kernel|Kernel|gmp" "$W/pari_conf.log" | head -5
export PERL5LIB=$W/pari/src/desc:$W/pari/src/desc/PARI
make -j"$(nproc)" lib-sta > "$W/pari_make.log" 2>&1 || { tail -20 "$W/pari_make.log"; exit 1; }
tail -3 "$W/pari_make.log"
ls -la libpari.a
echo PARI_ANDROID_DONE
