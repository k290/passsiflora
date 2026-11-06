# Which faults of the pre-filter can the test detect?   Linux, or Windows with MSYS2 / MinGW-w64.
#   bash fault_checks.sh [N] [test source]        default: N = 100 (1,000,000 polynomials), test_pdfilter.cpp
#   e.g.  bash fault_checks.sh 100 test_pdfilter_oracle.cpp    (the oracle test: both kinds of error, every row below)
#         bash fault_checks.sh 100 test_pdfilter_linux.cpp     (Eric's file with the small patch; Linux only)
#         ROWS="S13 SC SD" bash fault_checks.sh 100 test_pdfilter_oracle.cpp      (only those rows, and "none")
#   Eric's original file is in the history:  git show ccb629f:pdfilter_Test/test_pdfilter_eric.cpp > /tmp/eric_original.cpp
#         and then  bash fault_checks.sh 100 /tmp/eric_original.cpp  shows what it detects.
# The test is built once and linked against ../src/pdfilter.f90 as it is and against deliberately broken versions of it.  It is run
# with the AVX2/FMA build and with the baseline build of the filter, and the number of failures it reports is printed ("none" is
# the unbroken filter and must be 0).  A test that finds 0 for a fault cannot see it.  For a test that prints "Missed
# rejections" (test_pdfilter_oracle.cpp) the number in brackets is the part of the failures that are polynomials the broken
# filter wrongly KEPT; the tests that use only valid polynomials count the ones it wrongly REJECTED.
# Floating-point faults (the integer filter does not use floating point, so they must give 0 there):
#   x87             the filter compiled with x87 arithmetic (-mfpmath=387), as a 32-bit x86 build would be
#   up              the floating point rounding mode set to "up" while the filter runs (rmshim.c)
# Prime ranges (docs/whitepaper.tex: Section 4 "Exact integer arithmetic" has the proofs, its Section 4.4 "The size limits,
# and the guard" the limits; the expected result of each row was written down before the row ran; "directed" = a polynomial
# found by limit_probe.c that must make the row fail):
#   name  primes                                   proof (whitepaper)                             expected
#   S8    primes below 2^8 (241, 89, 41)             the lower bounds of pdf_init only choose     0
#   S13   primes in [8000, 8192)                     Lemma mont (q <= 8192)                       0
#   SC    8369, 8329: mont has counterexamples       Theorem "limits of the whole filter" (q <= 8385)  0
#   SF    primes in (8385, 8761]                     none (above the proved bound)                0 (no event in simulation)
#   S2    primes up to 2^14 (the 30 from 16369 down) none                                         0 (no event in simulation)
#   SD    primes in (19000, 19961]                   none                                         0 or a few; directed: failures
#   SE    primes in (21000, 23981]                   none                                         failures (seen in simulation)
#   S15   primes in (30000, 32767]                   none                                         many failures
#   SG    primes in (8000, 8761], guard on           the guard skips the 6 above 8385 (6 lines on 0
#                                                    stderr); the filter takes 8369 ... 8009
# The rows beyond the proved limits (NOLIM below) are built with -DPDF_NO_LIMITS, which switches off pdf_init's guard (the
# primes outside the limits would otherwise be skipped); every other row keeps the guard, and SG shows that it works.
# Faults of the integer arithmetic (one line each; the oracle test reports which strata found them):
#   Z1  zero test +-q   chain_block takes a pivot of +-q (0 mod q) for nonzero
#   Z2  general +-q     chain_general takes a value of +-q for nonzero (degree drops missed)
#   E1  Euler -R        chain_block recognises -1 only as negr, not as negr - q
#   U1  unsigned top    the top 21-bit piece of a coefficient is shifted logically (wrong for negative coefficients)
#   D1  load 2^20       the loading constants are 2^(20i)*R instead of 2^(21i)*R (wrong for coefficients of 2^21 and more and all negative ones)
#   G1  general parity  chain_general multiplies the even powers of lc(V) into its symbol instead of the odd ones
#   K1  last stage      the filter skips its last stage (only polynomials that reach stage 30 can show it)
# Environment: PARI and OUT as for build-linux.sh / build-windows.sh (OUT default: fault_build), SEED (default 12345), ROWS.
cd "$(dirname "$0")"
. ./pari-config.sh
N=${1:-100}
TEST=${2:-test_pdfilter.cpp}
OUT=${OUT:-fault_build}
SEED=${SEED:-12345}
ALL="none x87 up S2 Z1 Z2 E1 U1 D1 G1 K1 S8 S13 SC SF SD SE S15 SG"
NOLIM="S2 SF SD SE S15"
ROWS="none ${ROWS:-$ALL}"
ROWS=$(echo $ROWS | tr ' ' '\n' | awk '!seen[$0]++' | tr '\n' ' ')
mkdir -p "$OUT"
F3="-O3 -march=x86-64 -mavx2 -mno-popcnt -I../src -cpp"
FB="-O3 -march=x86-64 -mtune=generic -I../src -cpp"

filter() {   # name source.f90 extra-flags entry-suffix  ->  $OUT/<name>_v3.o  $OUT/<name>_base.o
  $FC $F3 $3 -J "$OUT" -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" "-DPDF_FILTER_NAME='pdf_filter_v3$4'" "-DPDF_STAGES_NAME='pdf_stages_v3'" -c "$2" -o "$OUT/$1_v3.o" || exit 1
  $FC $FB $3 -J "$OUT" -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" "-DPDF_FILTER_NAME='pdf_filter_base$4'" "-DPDF_STAGES_NAME='pdf_stages_base'" -c "$2" -o "$OUT/$1_base.o" || exit 1
}
mutate() {   # out.f90 sed-expression count [file] : checks that every expression matched
  local src=${4:-../src/pdfilter.f90}
  sed -e "$2" "$src" > "$1"
  [ "$(diff "$src" "$1" | grep -c '^>')" = "$3" ] || { echo "the sed expression did not change $1 as expected"; exit 1; }
}
small() { echo "s/QSTART = 2_i8\*\*13 - 1, QSTOP = 256$/QSTART = $1_i8, QSTOP = $2/"; }   # start, lower bound (exclusive)
for f in $ROWS; do
  case $f in
    none) filter none ../src/pdfilter.f90 ;;
    x87)  filter x87 ../src/pdfilter.f90 -mfpmath=387 ;;
    up)   filter up ../src/pdfilter.f90 "" _in ;;
    S2)   mutate "$OUT/pdfilter_S2.f90" "$(small 16383 8192)" 1 ;;
    S8)   mutate "$OUT/pdfilter_S8.f90" "$(small 255 2)" 1 ;;
    S13)  mutate "$OUT/pdfilter_S13.f90" "$(small 8191 7999)" 1 ;;
    SC)   mutate "$OUT/pdfilter_SC.f90" "$(small 8385 8328)" 1 ;;
    SF)   mutate "$OUT/pdfilter_SF.f90" "$(small 8761 8385)" 1 ;;
    SG)   mutate "$OUT/pdfilter_SG.f90" "$(small 8761 8000)" 1 ;;
    SD)   mutate "$OUT/pdfilter_SD.f90" "$(small 19961 19000)" 1 ;;
    SE)   mutate "$OUT/pdfilter_SE.f90" "$(small 23981 21000)" 1 ;;
    S15)  mutate "$OUT/pdfilter_S15.f90" "$(small 32767 30000)" 1 ;;
    Z1)   mutate "$OUT/pdfilter_Z1.f90" 's/ \.and\. abs(p(:, m-1, ir)) \/= q)/)/' 1 ;;
    Z2)   mutate "$OUT/pdfilter_Z2.f90" 's/(p(:, 0, iw) == 0 \.or\. abs(p(:, 0, iw)) == q)/(p(:, 0, iw) == 0)/' 1 ;;
    E1)   mutate "$OUT/pdfilter_E1.f90" 's/(ok \.and\. (acc == negr \.or\. acc == negr - q))/(ok .and. (acc == negr))/' 1 ;;
    U1)   mutate "$OUT/pdfilter_U1.f90" 's/int(shifta(c, 42), i4), &$/int(shiftr(c, 42), i4), \&/' 1 ;;
    D1)   mutate "$OUT/pdfilter_D1.f90" 's/r = modulo(r\*2_i8\*\*21, q)/r = modulo(r*2_i8**20, q)/' 1 ;;
    K1)   # the stage loop of pdf_filter (arith_ok has a loop with the same line): the line after "probe = ..."
          mutate "$OUT/pdfilter_K1.f90" '/probe = (ncall == 1)/{n;s/^    do j = 1, ns$/    do j = 1, ns - 1/;}' 1 ;;
    G1)   mutate "$OUT/pdfilter_G1.f90" 's/act > 0 \.and\. iand(ev, int(1, i2)) \/= 0)/act > 0 .and. iand(ev, int(1, i2)) == 0)/' 1 ;;
    *)    echo "unknown row $f"; exit 1 ;;
  esac
  case $f in
    none|x87|up) ;;
    *) case " $NOLIM " in *" $f "*) filter $f "$OUT/pdfilter_$f.f90" -DPDF_NO_LIMITS ;; *) filter $f "$OUT/pdfilter_$f.f90" ;; esac ;;
  esac
done
$CC -O2 -c rmshim.c -o "$OUT/rmshim.o" || exit 1
$CXX -O3 -w $INC -c -o "$OUT/test.o" "$TEST" || exit 1
for f in $ROWS; do
  if [ $f = up ]; then $CXX -O3 $LINK -o "$OUT/test_up$EXE" "$OUT/test.o" "$OUT/up_v3.o" "$OUT/up_base.o" "$OUT/rmshim.o" $LIBS || exit 1
  else $CXX -O3 $LINK -o "$OUT/test_$f$EXE" "$OUT/test.o" "$OUT/${f}_v3.o" "$OUT/${f}_base.o" $LIBS || exit 1; fi
done

# Tests that take a seed get the same one for every run (paired comparison); the others seed from the clock.
if grep -q "seed" "$TEST"; then ARGS="$SEED"; else ARGS=""; fi
run() {   # fault build blocks
  local fe=""; [ "$2" = base ] && fe=1
  if [ -n "$ARGS" ]; then NF_FORCE_BASE=$fe "$OUT/test_$1$EXE" "$3" $ARGS > "$OUT/res_$1_$2.txt" 2>&1
  else NF_FORCE_BASE=$fe "$OUT/test_$1$EXE" "$3" > "$OUT/res_$1_$2.txt" 2>&1; fi
}
launch() { run "$@" & [ -z "$ARGS" ] && sleep 1.2; }
for f in $ROWS; do
  launch $f avx2 $N; launch $f base $N
done
wait
failures() { grep "Total failures" "$OUT/res_$1_$2.txt" | sed 's/.*= //; s/ .*//'; }
missed()   { grep "Missed rejections" "$OUT/res_$1_$2.txt" | sed 's/.*= //; s/ .*//'; }
tested()   { grep "Total polynomials tested" "$OUT/res_$1_$2.txt" | sed 's/.*= //'; }
cell()     { local m; m=$(missed $1 $2); echo "$(failures $1 $2)${m:+ ($m)}"; }
row()      { case " $ROWS " in *" $2 "*) printf "%-26s %12s %16s %16s\n" "$1" "$(tested $2 avx2)" "$(cell $2 avx2)" "$(cell $2 base)" ;; esac; }
echo "test: $TEST   seed: ${ARGS:-clock}   $(grep -h 'pre-filter' "$OUT/res_none_avx2.txt" | head -1)"
printf "%-26s %12s %16s %16s\n" "fault" "polynomials" "AVX2 build" "baseline"
row "none (must be 0)"       none
row "x87 arithmetic"         x87
row "rounding mode up"       up
row "primes to 2^14"         S2
row "zero test +-q"          Z1
row "general +-q"            Z2
row "Euler -R"               E1
row "unsigned top"           U1
row "load 2^20"              D1
row "general parity"         G1
row "last stage"             K1
row "S8  primes < 2^8"       S8
row "S13 primes [8000,8192)" S13
row "SC  primes 8329, 8369"  SC
row "SF  primes (8385,8761]" SF
row "SD  (19000,19961]"      SD
row "SE  (21000,23981]"      SE
row "S15 (30000,32767]"      S15
row "SG  (8000,8761], guard on" SG
