# Which faults of the pre-filter can the test detect?   Linux, or Windows with MSYS2 / MinGW-w64.
#   bash fault_checks.sh [N] [test source]        default: N = 100 (1,000,000 polynomials), test_pdfilter.cpp
#   e.g.  bash fault_checks.sh 100 test_pdfilter_eric.cpp     to see what Eric's version of the test detects
# The test is built once and linked against pdfilter.f90 as it is and against deliberately broken versions of it.  It is run
# with the AVX2/FMA build and with the baseline build of the filter, and the number of valid polynomials that the broken
# filter wrongly REJECTED is printed ("none" is the unbroken filter and must be 0).  A test that finds 0 for a fault cannot see it.
#   x87           the filter compiled with x87 arithmetic (-mfpmath=387), as a 32-bit x86 build would be
#   round up      the floating point rounding mode set to "up" while the filter runs (rmshim.c)
#   large primes  the large filter primes taken up to 2^29 instead of below 2^25: rare (the run is 3 times longer)
#   small primes  the small filter primes taken up to 2^14 instead of below 2^10: a gross error
# A filter that rejects too little (rounding down, -Ofast) cannot be seen by a test that only uses polynomials that must be kept.
# Environment: PARI and OUT as for build-linux.sh / build-windows.sh (OUT default: fault_build), SEED (default 12345).
cd "$(dirname "$0")"
. ./pari-config.sh
N=${1:-100}
TEST=${2:-test_pdfilter.cpp}
OUT=${OUT:-fault_build}
SEED=${SEED:-12345}
mkdir -p "$OUT"
F3="-O3 -march=x86-64-v3 -I../src -cpp"
FB="-O3 -march=x86-64 -mtune=generic -I../src -cpp"

filter() {   # name source.f90 extra-flags entry-suffix  ->  $OUT/<name>_v3.o  $OUT/<name>_base.o
  $FC $F3 $3 -J "$OUT" -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" "-DPDF_FILTER_NAME='pdf_filter_v3$4'" -c "$2" -o "$OUT/$1_v3.o" || exit 1
  $FC $FB $3 -J "$OUT" -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" "-DPDF_FILTER_NAME='pdf_filter_base$4'" -c "$2" -o "$OUT/$1_base.o" || exit 1
}
mutate() {   # out.f90 sed-expression count : checks that every expression matched
  sed -e "$2" pdfilter.f90 > "$1"
  [ "$(diff pdfilter.f90 "$1" | grep -c '^>')" = "$3" ] || { echo "the sed expression did not change $1 as expected"; exit 1; }
}
mutate "$OUT/pdfilter_L2.f90" 's/q = 2_i8\*\*25 - 1/q = 2_i8**29 - 1/; s/q > 2_i8\*\*24/q > 2_i8**28/' 2
mutate "$OUT/pdfilter_S2.f90" 's/q = 2_i8\*\*10 - 1/q = 2_i8**14 - 1/; s/q > 256)/q > 8192)/' 2
filter none pdfilter.f90
filter x87 pdfilter.f90 -mfpmath=387
filter up pdfilter.f90 "" _in
filter L2 "$OUT/pdfilter_L2.f90"
filter S2 "$OUT/pdfilter_S2.f90"
$CC -O2 -c rmshim.c -o "$OUT/rmshim.o" || exit 1
$CXX -O3 -w $INC -c -o "$OUT/test.o" "$TEST" || exit 1
for f in none x87 L2 S2; do
  $CXX -O3 $LINK -o "$OUT/test_$f$EXE" "$OUT/test.o" "$OUT/${f}_v3.o" "$OUT/${f}_base.o" $LIBS || exit 1
done
$CXX -O3 $LINK -o "$OUT/test_up$EXE" "$OUT/test.o" "$OUT/up_v3.o" "$OUT/up_base.o" "$OUT/rmshim.o" $LIBS || exit 1

# Tests that take a seed get the same one for every run (paired comparison); the others seed from the clock.
if grep -q "seed" "$TEST"; then ARGS="$SEED"; else ARGS=""; fi
run() {   # fault build blocks
  local fe=""; [ "$2" = base ] && fe=1
  if [ -n "$ARGS" ]; then NF_FORCE_BASE=$fe "$OUT/test_$1$EXE" "$3" $ARGS > "$OUT/res_$1_$2.txt" 2>&1
  else NF_FORCE_BASE=$fe "$OUT/test_$1$EXE" "$3" > "$OUT/res_$1_$2.txt" 2>&1; fi
}
launch() { run "$@" & [ -z "$ARGS" ] && sleep 1.2; }
for f in none x87 up S2; do launch $f avx2 $N; launch $f base $N; done
launch L2 avx2 $((3*N)); launch L2 base $((3*N))
wait
failures() { grep "Total failures" "$OUT/res_$1_$2.txt" | sed 's/.*= //'; }
tested()   { grep "Total polynomials tested" "$OUT/res_$1_$2.txt" | sed 's/.*= //'; }
echo "test: $TEST   seed: ${ARGS:-clock}   $(grep -h 'pre-filter' "$OUT/res_none_avx2.txt" | head -1)"
printf "%-26s %12s %12s %12s\n" "fault" "polynomials" "AVX2 build" "baseline"
printf "%-26s %12s %12s %12s\n" "none (must be 0)"        "$(tested none avx2)" "$(failures none avx2)" "$(failures none base)"
printf "%-26s %12s %12s %12s\n" "x87 arithmetic"          "$(tested x87 avx2)"  "$(failures x87 avx2)"  "$(failures x87 base)"
printf "%-26s %12s %12s %12s\n" "rounding mode up"        "$(tested up avx2)"   "$(failures up avx2)"   "$(failures up base)"
printf "%-26s %12s %12s %12s\n" "large primes (rare)"     "$(tested L2 avx2)"   "$(failures L2 avx2)"   "$(failures L2 base)"
printf "%-26s %12s %12s %12s\n" "small primes (gross)"    "$(tested S2 avx2)"   "$(failures S2 avx2)"   "$(failures S2 base)"
