# Which faults of the pre-filter can the test detect?   Linux, or Windows with MSYS2 / MinGW-w64.
#   bash fault_checks.sh [N] [test source]        default: N = 100 (1,000,000 polynomials), test_pdfilter.cpp
#   e.g.  bash fault_checks.sh 100 test_pdfilter_oracle.cpp    (the oracle test: both kinds of error, every row below)
#         bash fault_checks.sh 100 test_pdfilter_linux.cpp     (Eric's file with the small patch; Linux only)
#   Eric's original file is in the history:  git show ccb629f:pdfilter_Test/test_pdfilter_eric.cpp > /tmp/eric_original.cpp
#         and then  bash fault_checks.sh 100 /tmp/eric_original.cpp  shows what it detects.
# The test is built once and linked against pdfilter.f90 as it is and against deliberately broken versions of it.  It is run
# with the AVX2/FMA build and with the baseline build of the filter, and the number of failures it reports is printed ("none" is
# the unbroken filter and must be 0).  A test that finds 0 for a fault cannot see it.  For a test that prints "Missed
# rejections" (test_pdfilter_oracle.cpp) the number in brackets is the part of the failures that are polynomials the broken
# filter wrongly KEPT; the tests that use only valid polynomials count the ones it wrongly REJECTED.
# Floating-point faults (the integer filter does not use floating point, so they must give 0 there):
#   x87             the filter compiled with x87 arithmetic (-mfpmath=387), as a 32-bit x86 build would be
#   round up        the floating point rounding mode set to "up" while the filter runs (rmshim.c)
# Prime ranges:
#   large primes    the large filter primes taken up to 2^29 instead of below 2^25 (the run is 3 times longer).  For the
#                   floating-point filter beyond its exact range; for the integer filter inside it (sredc needs q < 2^30), so
#                   a test that uses the filter's own primes (the oracle test reads them from the filter) must report 0
#   small primes    the small filter primes taken up to 2^14 instead of below 2^10.  A gross error for the floating-point
#                   filter.  For the integer filter outside the proven range (mont needs q < 2^13): its results can leave
#                   (-q, q), but so rarely that no wrong verdict was seen in 1,000,000 polynomials (pdfilter_chains.tex, 13.7)
# Faults of the integer arithmetic (one line each; the oracle test reports which strata found them):
#   zero test +-q   chain_block_s takes a pivot of +-q (0 mod q) for nonzero
#   general +-q     chain_general_s takes a value of +-q for nonzero (degree drops missed)
#   Euler -R        chain_block and chain_block_s recognise -1 only as negr, not as negr - q
#   unsigned top    the top 21-bit piece of a coefficient is shifted logically (wrong for negative coefficients)
#   load 2^20       the loading constants are 2^(20i)*R instead of 2^(21i)*R (wrong for coefficients of 2^21 and more and all negative ones)
#   general parity  chain_general multiplies the even powers of lc(V) into its symbol instead of the odd ones
#   last stage      the filter skips its last stage (only polynomials that reach stage 30 can show it)
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
  $FC $F3 $3 -J "$OUT" -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" "-DPDF_FILTER_NAME='pdf_filter_v3$4'" "-DPDF_STAGES_NAME='pdf_stages_v3'" -c "$2" -o "$OUT/$1_v3.o" || exit 1
  $FC $FB $3 -J "$OUT" -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" "-DPDF_FILTER_NAME='pdf_filter_base$4'" "-DPDF_STAGES_NAME='pdf_stages_base'" -c "$2" -o "$OUT/$1_base.o" || exit 1
}
mutate() {   # out.f90 sed-expression count [file] : checks that every expression matched
  local src=${4:-pdfilter.f90}
  sed -e "$2" "$src" > "$1"
  [ "$(diff "$src" "$1" | grep -c '^>')" = "$3" ] || { echo "the sed expression did not change $1 as expected"; exit 1; }
}
mutate "$OUT/pdfilter_L2.f90" 's/q = 2_i8\*\*25 - 1/q = 2_i8**29 - 1/; s/q > 2_i8\*\*24/q > 2_i8**28/' 2
mutate "$OUT/pdfilter_S2.f90" 's/q = 2_i8\*\*10 - 1/q = 2_i8**14 - 1/; s/q > 256)/q > 8192)/' 2
mutate "$OUT/pdfilter_Z1.f90" 's/ \.and\. abs(p(:, m-1, ir)) \/= q)/)/' 1
mutate "$OUT/pdfilter_Z2.f90" 's/^#define GZERO(x) (x == 0 \.or\. abs(x) == q)/#define GZERO(x) (x == 0)/' 1
mutate "$OUT/pdfilter_E1.f90" 's/(acc == negr \.or\. acc == negr - q)/(acc == negr)/' 2
mutate "$OUT/pdfilter_U1.f90" 's/p2A(l, k) = int(shifta(x, 42), i4)/p2A(l, k) = int(shiftr(x, 42), i4)/' 1
mutate "$OUT/pdfilter_D1.f90" 's/r = modulo(r\*2_i8\*\*21, q)/r = modulo(r*2_i8**20, q)/' 2
mutate "$OUT/pdfilter_K1.f90" 's/^    do j = 1, ns + nq$/    do j = 1, ns + nq - 1/' 1
# chain_general.inc is included from the directory of the including file first: a broken copy next to a copy of the filter
mkdir -p "$OUT/G1"
cp pdfilter.f90 "$OUT/G1/pdfilter.f90"
mutate "$OUT/G1/chain_general.inc" 's/act > 0 \.and\. iand(ev, int(1, GK)) \/= 0)/act > 0 .and. iand(ev, int(1, GK)) == 0)/' 1 ../src/chain_general.inc
filter none pdfilter.f90
filter x87 pdfilter.f90 -mfpmath=387
filter up pdfilter.f90 "" _in
filter L2 "$OUT/pdfilter_L2.f90"
filter S2 "$OUT/pdfilter_S2.f90"
INT="Z1 Z2 E1 U1 D1 K1"
for f in $INT; do filter $f "$OUT/pdfilter_$f.f90"; done
filter G1 "$OUT/G1/pdfilter.f90"
INT="$INT G1"
$CC -O2 -c rmshim.c -o "$OUT/rmshim.o" || exit 1
$CXX -O3 -w $INC -c -o "$OUT/test.o" "$TEST" || exit 1
for f in none x87 L2 S2 $INT; do
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
for f in none x87 up S2 $INT; do launch $f avx2 $N; launch $f base $N; done
launch L2 avx2 $((3*N)); launch L2 base $((3*N))
wait
failures() { grep "Total failures" "$OUT/res_$1_$2.txt" | sed 's/.*= //; s/ .*//'; }
missed()   { grep "Missed rejections" "$OUT/res_$1_$2.txt" | sed 's/.*= //; s/ .*//'; }
tested()   { grep "Total polynomials tested" "$OUT/res_$1_$2.txt" | sed 's/.*= //'; }
cell()     { local m; m=$(missed $1 $2); echo "$(failures $1 $2)${m:+ ($m)}"; }
row()      { printf "%-26s %12s %16s %16s\n" "$1" "$(tested $2 avx2)" "$(cell $2 avx2)" "$(cell $2 base)"; }
echo "test: $TEST   seed: ${ARGS:-clock}   $(grep -h 'pre-filter' "$OUT/res_none_avx2.txt" | head -1)"
printf "%-26s %12s %16s %16s\n" "fault" "polynomials" "AVX2 build" "baseline"
row "none (must be 0)"       none
row "x87 arithmetic"         x87
row "rounding mode up"       up
row "large primes to 2^29"   L2
row "small primes to 2^14"   S2
row "zero test +-q"          Z1
row "general +-q"            Z2
row "Euler -R"               E1
row "unsigned top"           U1
row "load 2^20"              D1
row "general parity"         G1
row "last stage"             K1
