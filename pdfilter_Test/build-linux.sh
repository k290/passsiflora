# Builds the pre-filter and the test program (Linux).  Run it in this directory:   bash build-linux.sh
#   PARI          see pari-config.sh (a compiled PARI source tree; default Eric's, else the PARI installed on the system)
#   FFLAGS_EXTRA  extra gfortran flags for the filter, for experiments, e.g.  FFLAGS_EXTRA=-mfpmath=387 bash build-linux.sh
#                 (x87 arithmetic, as in a 32-bit x86 build; the test should then report failures)
#   OUT           directory for the object files and the program (default: build)
#   FILTER        a directory whose pdfilter.f90 (with the chain_general.inc next to it) is built instead of ../src/pdfilter.f90,
#                 e.g. the filter of another commit:  git archive main src | tar -x -C /tmp/main && FILTER=/tmp/main/src bash build-linux.sh
#   REF           a directory with a second pdfilter.f90 (and its chain_general.inc), linked as pdf_init_ref / pdf_filter_ref for
#                 the differential check of test_pdfilter_oracle.cpp (other tests ignore it)
#   EXTRA_CXX     further C++ sources compiled and linked with the test, LIBS_EXTRA further libraries; test_poldisc.cpp needs upstream's
#                 polDiscTest_cpu.cpp and GMP:  EXTRA_CXX=<GetDecicsSrc>/polDiscTest_cpu.cpp LIBS_EXTRA=-lgmp TEST=test_poldisc.cpp bash build-linux.sh
#   EXTRA_CXXFLAGS further flags for the EXTRA_CXX sources (poldisc_checks.sh: -fwrapv, as the application compiles polDiscTest_cpu.cpp)
# Then run   build/test_pdfilter 100   (100 blocks of 10000 polynomials; a second argument is the seed).
cd "$(dirname "$0")"
. ./pari-config.sh
OUT=${OUT:-build}
mkdir -p "$OUT"

# Build the fortran pdfilter from ../src (the tests never compile a copy).  chain_general.inc is found next to it, or through -I../src
# for a FILTER without one.  PDF_STAGES_NAME compiles the test hook that test_pdfilter_oracle.cpp reads (the
# application's build does not define it; a filter without the hook is tested on its verdicts only).
SRC=${FILTER:-../src}/pdfilter.f90      # the shipped source: src/ is the only copy of the filter
$FC -O3 -march=x86-64 -mavx2 -mno-popcnt $FFLAGS_EXTRA -I../src -cpp -J "$OUT" -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" "-DPDF_FILTER_NAME='pdf_filter_v3'" "-DPDF_STAGES_NAME='pdf_stages_v3'" -c "$SRC" -o "$OUT/pdfilter_v3.o" || exit 1
$FC -O3 -march=x86-64 -mtune=generic $FFLAGS_EXTRA -I../src -cpp -J "$OUT" -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" "-DPDF_FILTER_NAME='pdf_filter_base'" "-DPDF_STAGES_NAME='pdf_stages_base'" -c "$SRC" -o "$OUT/pdfilter_base.o" || exit 1

REFOBJ=
if [ -n "$REF" ]; then
  $FC -O3 -march=x86-64 -mtune=generic -cpp -J "$OUT" -DPDF_MOD=pdfilter_ref "-DPDF_INIT_NAME='pdf_init_ref'" "-DPDF_FILTER_NAME='pdf_filter_ref'" -c "$REF/pdfilter.f90" -o "$OUT/pdfilter_ref.o" || exit 1
  REFOBJ="$OUT/pdfilter_ref.o"
fi

# Build the test program
TEST=${TEST:-test_pdfilter.cpp}      # the test source; TEST=test_pdfilter_linux.cpp builds Eric's file with the small patch (Linux only)
PROG=$(basename "${TEST%.cpp}")
$CXX -O3 -w $INC -c -o "$OUT/$PROG.o" "$TEST" || exit 1
EXTRAOBJ=
for x in $EXTRA_CXX; do
  o="$OUT/extra_$(basename "${x%.cpp}").o"
  $CXX -O3 -w $EXTRA_CXXFLAGS -c -o "$o" "$x" || exit 1
  EXTRAOBJ="$EXTRAOBJ $o"
done
$CXX -O3 $LINK -o "$OUT/$PROG$EXE" "$OUT/$PROG.o" "$OUT/pdfilter_v3.o" "$OUT/pdfilter_base.o" $REFOBJ $EXTRAOBJ $LIBS $LIBS_EXTRA || exit 1
echo "built $OUT/$PROG$EXE"
