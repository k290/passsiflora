# Builds PARI 2.8-1711 (bundled with upstream) as a static library for MinGW-w64: cross-compiled on Linux (Configure runs its
# test programs through WSL's Windows interop) or natively in MSYS2 MINGW64.
# Environment (the defaults are the author's machine):
#   PARI_SRC  upstream's Pari/pari-2.8-1711-ge5c317c     WP     where PARI is built (the WP of build-windows-gpu.sh)
#   GMP       GMP for MinGW (include/, lib/)             TOOLS  where the helper scripts mingwcc and winrun are written
#   AR, RANLIB   LOG (- = no log file)
# mingwcc: Configure wants the -o target without the .exe suffix, and PARI 2.8 needs gnu11 and -fcommon with newer compilers.
PARI_SRC=${PARI_SRC:-/root/nf/upstream/Pari/pari-2.8-1711-ge5c317c}
WP=${WP:-/mnt/c/Users/Alp/AppData/Local/Temp/nfw/wpari}
GMP=${GMP:-/root/nf/wgmp}
TOOLS=${TOOLS:-/root/nf}
AR=${AR:-x86_64-w64-mingw32-ar}
RANLIB=${RANLIB:-x86_64-w64-mingw32-ranlib}
LOG=${LOG:-/root/nf/winpari4.log}
[ "$LOG" = "-" ] || exec > "$LOG" 2>&1
mkdir -p "$TOOLS"
cat > "$TOOLS/winrun" <<'W'
#!/bin/bash
# Run a Windows test exe from PARI's Configure and strip CRLF from its output.
"$@" | tr -d '\r'; exit ${PIPESTATUS[0]}
W
cat > "$TOOLS/mingwcc" <<'W'
#!/bin/bash
out=""; prev=""
for a in "$@"; do if [ "$prev" = "-o" ]; then out="$a"; fi; prev="$a"; done
x86_64-w64-mingw32-gcc -std=gnu11 -fcommon -Wno-implicit-function-declaration -Wno-int-conversion -Wno-incompatible-pointer-types "$@"
rc=$?
if [ -n "$out" ] && [ -f "$out.exe" ] && [ ! -f "$out" ]; then mv "$out.exe" "$out"; fi
exit $rc
W
chmod +x "$TOOLS/winrun" "$TOOLS/mingwcc"
rm -rf "$WP" && cp -r "$PARI_SRC" "$WP"
cd "$WP" && rm -rf Olinux-x86_64
RUNTEST="$TOOLS/winrun" CC="$TOOLS/mingwcc" RANLIB=$RANLIB AR=$AR ./Configure --static --host=x86_64-mingw --with-gmp="$GMP" --without-readline 2>&1 | grep -E 'sizeof|kernel|Kernel|gmp|GMP|Abort'
O=$(ls -d O* | head -1); echo "objdir=$O"; grep -E "sizeof_long|doubleformat|kernlvl" $O/pari.cfg
export PERL5LIB="$WP/src/desc:$WP/src/desc/PARI"
cd $O && make -j4 lib-sta 2>&1 | grep -E 'error|Error' | head
ls -la libpari.a && echo WINPARI_OK
