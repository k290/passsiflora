# The compiled filter on complete sets of polynomials at tiny primes (docs/whitepaper.tex, Section 8, "The
# compiled filter on complete sets at tiny primes").  For each tiny prime q (= 1 mod 4, not dividing 10) a filter is built from
# ../src/pdfilter.f90 with build-linux.sh whose prime table holds only q.  The only change is the two search constants QSTART
# and QSTOP (sed; the guard QMIN = 11 and every other line are as shipped).  test_tiny.cpp then runs every polynomial of a set
# through the AVX2 and the baseline build, each with the normal routing and with every lane through chain_general, and
# compares flag and stage with chi(disc f mod q).
#   bash tiny_checks.sh              all sets (25.5 billion polynomials, 4 runs each: 2.7 hours with JOBS=12 of 16 logical CPUs)
#   QUICK=1 bash tiny_checks.sh      a few minutes
# Environment: OUT (default tiny_build), JOBS (default: the number of CPUs), PARI as for build-linux.sh.  Exit status 1 on a failure.
cd "$(dirname "$0")"
OUT=${OUT:-tiny_build}
JOBS=${JOBS:-$(nproc 2>/dev/null || echo 4)}
mkdir -p "$OUT"
SRC=../src/pdfilter.f90
build() {   # q
  local d="$OUT/f_$1"; mkdir -p "$d"
  e="s/QSTART = 2_i8\*\*13 - 1, QSTOP = 256$/QSTART = $1_i8, QSTOP = $(($1 - 1))/"
  sed -e "$e" "$SRC" > "$d/pdfilter.f90"
  [ "$(diff "$SRC" "$d/pdfilter.f90" | grep -c '^>')" = 1 ] || { echo "the sed expression did not change the search constants"; exit 1; }
  FILTER="$d" OUT="$OUT/b_$1" TEST=test_tiny.cpp bash build-linux.sh > "$OUT/build_$1.log" 2>&1 || { cat "$OUT/build_$1.log"; exit 1; }
}
# jobs: "q mode start count", split into chunks
JOBLIST="$OUT/jobs.txt"; : > "$JOBLIST"
chunks() {   # q mode total nchunks
  local q=$1 mode=$2 total=$3 nc=$4 s=0 c
  for c in $(seq 1 $nc); do e=$(( total * c / nc )); echo "$q $mode $s $(( e - s ))" >> "$JOBLIST"; s=$e; done
}
if [ -n "$QUICK" ]; then
  PR="13 17"
  chunks 13 depressed 20000000 16
  chunks 17 zeros2 4000000 8
  for q in $PR; do chunks $q shifts 20000 2; done
else
  PR="13 17 29 37 41"
  chunks 13 depressed 10604499373 64                  # 13^9: every depressed polynomial mod 13
  chunks 17 zeros2 6975757441 48                      # 17^8: a9 = a8 = 0
  chunks 29 zeros4 594823321 8                        # 29^6
  chunks 37 zeros4 2565726409 24                      # 37^6
  chunks 41 zeros4 4750104241 32                      # 41^6
  for q in $PR; do chunks $q shifts 200000 4; done     # 200,000 random polynomials with all q shifts, lifted coefficients
fi
for q in $PR; do build $q; done
echo "built; $(wc -l < "$JOBLIST") jobs on $JOBS CPUs"
run_job() {   # q mode start count
  local f="$OUT/res_$1_$2_$3.txt"
  if "$OUT/b_$1/test_tiny" "$@" > "$f" 2>&1; then echo ok > "$f.status"; else echo FAIL > "$f.status"; fi
}
export -f run_job; export OUT
rm -f "$OUT"/res_*.txt "$OUT"/res_*.status
start=$SECONDS
xargs -P "$JOBS" -L 1 bash -c 'run_job "$@"' _ < "$JOBLIST"
# summary per (q, mode): sums over the chunks
nfail=$(cat "$OUT"/res_*.status | grep -c FAIL)
cat "$OUT"/res_*.txt | awk '
  $1 == "ORACLE" { k = $2 " " $3; if (!(k in T)) K[++n] = k; T[k] += $4; M[k] += $5; A[k] += $6; Z[k] += $7 }
  $1 == "SEQ"    { k = $2 " " $3; if (!((k, $4) in S)) { S[k, $4] = 1; NS[k]++ } }
  $1 == "RUN"    { r = $2 " " $3 " " $4; for (i = 5; i <= 14; i++) R[r, i] += $i }
  END {
    split("AVX2-normal AVX2-general base-normal base-general", vv, " ")
    for (j = 1; j <= n; j++) {
      k = K[j]
      printf "q=%s: %.0f polynomials, chi=-1 %.0f, abnormal chains %.0f, disc=0 %.0f, degree sequences %d\n", k, T[k], M[k], A[k], Z[k], NS[k]
      for (v = 1; v <= 4; v++) { r = k " " vv[v]
        printf "    %-13s polys %.0f rejected %.0f false %.0f missed %.0f wrong-stage %.0f | chain_block %.0f chain_general %.0f retried %.0f lockstep-it %.0f masked-it %.0f\n",
               vv[v], R[r,5], R[r,6], R[r,7], R[r,8], R[r,9], R[r,10], R[r,11], R[r,12], R[r,13], R[r,14] }
    }
  }'
grep -h "false rejection\|missed rejection" "$OUT"/res_*.txt | head -5
echo "tiny_checks: $(wc -l < "$JOBLIST") jobs, $nfail failed, $((SECONDS - start)) s"
[ "$nfail" = 0 ]
