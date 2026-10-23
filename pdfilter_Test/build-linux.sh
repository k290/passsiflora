# Builds the pre-filter and the test program (Linux).  Run it in this directory:   bash build-linux.sh
#   PARI          see pari-config.sh (a compiled PARI source tree; default Eric's, else the PARI installed on the system)
#   FFLAGS_EXTRA  extra gfortran flags for the filter, for experiments, e.g.  FFLAGS_EXTRA=-mfpmath=387 bash build-linux.sh
#                 (x87 arithmetic, as in a 32-bit x86 build; the test should then report failures)
#   OUT           directory for the object files and the program (default: build)
# Then run   build/test_pdfilter 100   (100 blocks of 10000 polynomials; a second argument is the seed).
cd "$(dirname "$0")"
. ./pari-config.sh
OUT=${OUT:-build}
mkdir -p "$OUT"

# Build the fortran pdfilter.  chain_general.inc, which pdfilter.f90 includes, is not in this directory: -I../src finds it there
# (a copy next to pdfilter.f90 would be used first).
$FC -O3 -march=x86-64-v3 $FFLAGS_EXTRA -I../src -cpp -J "$OUT" -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" "-DPDF_FILTER_NAME='pdf_filter_v3'" -c pdfilter.f90 -o "$OUT/pdfilter_v3.o" || exit 1
$FC -O3 -march=x86-64 -mtune=generic $FFLAGS_EXTRA -I../src -cpp -J "$OUT" -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" "-DPDF_FILTER_NAME='pdf_filter_base'" -c pdfilter.f90 -o "$OUT/pdfilter_base.o" || exit 1

# Build the test program
TEST=${TEST:-test_pdfilter.cpp}      # the test source; TEST=test_pdfilter_eric.cpp builds Eric's version
PROG=${TEST%.cpp}
$CXX -O3 -w $INC -c -o "$OUT/$PROG.o" "$TEST" || exit 1
$CXX -O3 $LINK -o "$OUT/$PROG$EXE" "$OUT/$PROG.o" "$OUT/pdfilter_v3.o" "$OUT/pdfilter_base.o" $LIBS || exit 1
echo "built $OUT/$PROG$EXE"
