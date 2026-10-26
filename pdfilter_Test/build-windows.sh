# Builds the pre-filter and the test program (Windows, MSYS2 / MinGW-w64 gfortran and g++).  Run it in this directory:
#   PARI=/c/pari-2.17.3 bash build-windows.sh
#   PARI          a PARI source tree compiled for MinGW-w64 (it must contain src/headers and Omingw-x86_64); required
#   FFLAGS_EXTRA  extra gfortran flags for the filter, for experiments, e.g.  FFLAGS_EXTRA=-mfpmath=387
#   OUT           directory for the object files and the program (default: build)
#   FILTER, REF   as for build-linux.sh
# Then run   build/test_pdfilter.exe 100   (100 blocks of 10000 polynomials; a second argument is the seed).
# It is the same as build-linux.sh; pari-config.sh chooses the Windows settings.
cd "$(dirname "$0")"
. ./pari-config.sh
OUT=${OUT:-build}
mkdir -p "$OUT"

SRC=${FILTER:+$FILTER/}pdfilter.f90
$FC -O3 -march=x86-64-v3 $FFLAGS_EXTRA -I../src -cpp -J "$OUT" -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" "-DPDF_FILTER_NAME='pdf_filter_v3'" "-DPDF_STAGES_NAME='pdf_stages_v3'" -c "$SRC" -o "$OUT/pdfilter_v3.o" || exit 1
$FC -O3 -march=x86-64 -mtune=generic $FFLAGS_EXTRA -I../src -cpp -J "$OUT" -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" "-DPDF_FILTER_NAME='pdf_filter_base'" "-DPDF_STAGES_NAME='pdf_stages_base'" -c "$SRC" -o "$OUT/pdfilter_base.o" || exit 1

REFOBJ=
if [ -n "$REF" ]; then
  $FC -O3 -march=x86-64 -mtune=generic -cpp -J "$OUT" -DPDF_MOD=pdfilter_ref "-DPDF_INIT_NAME='pdf_init_ref'" "-DPDF_FILTER_NAME='pdf_filter_ref'" -c "$REF/pdfilter.f90" -o "$OUT/pdfilter_ref.o" || exit 1
  REFOBJ="$OUT/pdfilter_ref.o"
fi

TEST=${TEST:-test_pdfilter.cpp}      # the test source; TEST=test_pdfilter_linux.cpp builds Eric's file with the small patch (Linux only)
PROG=$(basename "${TEST%.cpp}")
$CXX -O3 -w $INC -c -o "$OUT/$PROG.o" "$TEST" || exit 1
$CXX -O3 $LINK -o "$OUT/$PROG$EXE" "$OUT/$PROG.o" "$OUT/pdfilter_v3.o" "$OUT/pdfilter_base.o" $REFOBJ $LIBS || exit 1
echo "built $OUT/$PROG$EXE"
