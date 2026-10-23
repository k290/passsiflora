# Sourced by build-linux.sh, build-windows.sh and fault_checks.sh: finds PARI and sets the tools and the link flags.
#   PARI  a PARI source tree that has been compiled (it contains src/headers and, after the build, Olinux-x86_64 or Omingw-x86_64).
#         Linux:    default Eric's tree, and if that does not exist the PARI installed on the system (apt install libpari-dev).
#         Windows:  required (MSYS2 / MinGW-w64), e.g.  PARI=/c/pari-2.17.3 bash build-windows.sh
# Sets  CXX CC FC EXE INC LINK LIBS.
CXX=${CXX:-g++}
CC=${CC:-gcc}
FC=${FC:-gfortran}
case "$(uname -s)" in
  MINGW*|MSYS*|CYGWIN*)
    [ -d /c/msys64/mingw64/bin ] && export PATH=/c/msys64/mingw64/bin:$PATH
    EXE=.exe
    if [ ! -d "$PARI/src/headers" ]; then
      echo "Set PARI to a PARI source tree compiled for MinGW-w64 (it must contain src/headers and Omingw-x86_64)."
      exit 1
    fi
    INC="-I$PARI/src/headers -I$PARI/Omingw-x86_64"
    LINK="-static -static-libgcc -static-libstdc++"
    LIBS="-L$PARI/Omingw-x86_64 -l:libpari.a -l:libgmp.a -l:libgfortran.a -l:libquadmath.a -lm -lws2_32 -lpsapi"
    ;;
  *)
    EXE=
    PARI=${PARI:-/home/eric/BOINC_Project/Pari/pari-2.17.3}
    if [ -d "$PARI/src/headers" ]; then
      INC="-I$PARI/src/headers -I$PARI/Olinux-x86_64"
      LINK="-static -static-libgcc -static-libstdc++"
      LIBS="-L$PARI/Olinux-x86_64 -l:libpari.a -l:libgfortran.a -lquadmath -lm -l:libgmp.a"
    else
      echo "No PARI tree at $PARI: using the PARI installed on the system."
      INC=""
      for d in /usr/include/pari /usr/include/x86_64-linux-gnu/pari /usr/local/include/pari; do   # Debian/Ubuntu use the second
        if [ -f "$d/pari.h" ]; then INC="-I$d"; break; fi
      done
      [ -n "$INC" ] || { echo "pari.h not found: install libpari-dev, or set PARI to a compiled PARI source tree."; exit 1; }
      LINK=""
      LIBS="-lpari -lgmp -lm -lgfortran"
    fi
    ;;
esac
