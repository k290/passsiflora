# Checks the compiled integer arithmetic of the pre-filter (mont, sredc, ldres, inv32, signedmod of ../src/pdfilter.f90) against
# exact models, and the prime-size limits of pdfilter_chains.tex, Section 8.  Linux, or Windows with MSYS2 / MinGW-w64.
#   bash arith_checks.sh            all checks (about 10 minutes on 16 logical CPUs; bounds-large alone is about 5 CPU minutes)
#   QUICK=1 bash arith_checks.sh    a reduced set (about a minute)
# Environment: FC (default gfortran), OUT (default arith_build), JOBS (default: the number of CPUs).
# The function texts are not retyped: they are cut out of ../src/pdfilter.f90 (two independent ways, which must agree) into
# $OUT/arith_funcs.inc and included in arith_checks.f90, so that the compiler inlines and vectorises them as in the filter.
# They are compiled with the two flag sets of fault_checks.sh (AVX2 and baseline) and with -O0.  Exit status 1 if any check fails.
cd "$(dirname "$0")"
FC=${FC:-gfortran}
OUT=${OUT:-arith_build}
JOBS=${JOBS:-$(nproc 2>/dev/null || echo 4)}
SRC=../src/pdfilter.f90
mkdir -p "$OUT"

# ---- extraction
cut_out() {   # first-line last-line : the lines of $SRC from first to last, which must each occur exactly once
  local n1 n2
  [ "$(grep -cxF "$1" "$SRC")" = 1 ] && [ "$(grep -cxF "$2" "$SRC")" -ge 1 ] || { echo "extraction: '$1' not found exactly once in $SRC"; exit 1; }
  n1=$(grep -nxF "$1" "$SRC" | cut -d: -f1)
  n2=$(awk -v n1="$n1" -v e="$2" 'NR > n1 && $0 == e { print NR; exit }' "$SRC")
  sed -n "${n1},${n2}p" "$SRC" > "$OUT/cut_a.txt"
  awk -v s="$1" -v e="$2" '$0 == s { p = 1 } p { print } p && $0 == e { exit }' "$SRC" > "$OUT/cut_b.txt"
  cmp -s "$OUT/cut_a.txt" "$OUT/cut_b.txt" || { echo "extraction: the two methods differ for '$1'"; exit 1; }
  cat "$OUT/cut_a.txt"
}
KINDS='  integer, parameter :: i2 = c_int16_t, i4 = c_int32_t, i8 = c_int64_t'
[ "$(grep -cxF "$KINDS" "$SRC")" = 1 ] || { echo "extraction: the kind parameters of $SRC changed"; exit 1; }
echo "$KINDS" > "$OUT/arith_kinds.inc"
{
  cut_out '  pure integer(i8) function inv32(q)' '  end function inv32'
  cut_out '  pure integer(i8) function signedmod(x, b)' '  end function signedmod'
  cut_out '  elemental integer(i2) function mont(a, b, q, qinv)' '  end function mont'
  cut_out '  elemental integer(i4) function sredc(x, q, qinv)' '  end function sredc'
  cut_out '  elemental integer(i4) function ldres(p0, p1, p2, d0, d1, d2, q, qinv)' '  end function ldres'
} > "$OUT/arith_funcs.inc" || exit 1
echo "extracted from $SRC: $(grep -c 'function' "$OUT/arith_funcs.inc") function lines, $(wc -l < "$OUT/arith_funcs.inc") lines in all"

# ---- builds: the flag sets of fault_checks.sh, and -O0
F3="-O3 -march=x86-64-v3"
FB="-O3 -march=x86-64 -mtune=generic"
for v in avx2 base O0; do
  case $v in avx2) F=$F3 ;; base) F=$FB ;; O0) F=-O0 ;; esac
  mkdir -p "$OUT/$v"
  $FC $F -cpp -I"$OUT" -J "$OUT/$v" arith_checks.f90 -o "$OUT/arith_checks_$v" || exit 1
done

# ---- jobs: "variant part arguments", run in parallel
JOBLIST="$OUT/jobs.txt"; : > "$JOBLIST"
if [ -n "$QUICK" ]; then
  for v in avx2 base O0; do
    echo "$v mont-pairs 3 1100 prime" >> "$JOBLIST"
    echo "$v mont-pairs 8237 8243 odd" >> "$JOBLIST"
    echo "$v inv32" >> "$JOBLIST"
  done
  echo "avx2 mont-scan 8193 9001" >> "$JOBLIST"
  echo "avx2 bounds-small 3 9001" >> "$JOBLIST"
else
  # all primes below 2^13 (16 ranges of about equal work, q^2 each) and all odd q in [8193, 8300], every pair
  B=$(awk 'BEGIN { for (k = 1; k <= 16; k++) printf "%d ", 8192 * (k / 16) ^ (1 / 3) }')
  for v in avx2 base; do
    lo=3; for b in $B; do echo "$v mont-pairs $lo $b prime" >> "$JOBLIST"; lo=$((b + 1)); done
    for lo in 8193 8221 8249 8277; do echo "$v mont-pairs $lo $((lo + 26)) odd" >> "$JOBLIST"; done
    echo "$v mont-count 16369" >> "$JOBLIST"
    echo "$v mont-scan 8193 16383" >> "$JOBLIST"
    echo "$v mont-extremes 8193 16383" >> "$JOBLIST"
    echo "$v sredc" >> "$JOBLIST"
    echo "$v ldres" >> "$JOBLIST"
    echo "$v inv32" >> "$JOBLIST"
  done
  echo "O0 mont-pairs 3 1100 prime" >> "$JOBLIST"
  echo "O0 mont-pairs 8237 8243 odd" >> "$JOBLIST"
  echo "O0 sredc" >> "$JOBLIST"
  echo "O0 ldres" >> "$JOBLIST"
  echo "O0 inv32" >> "$JOBLIST"
  echo "avx2 bounds-small 3 32767" >> "$JOBLIST"
  lo=3; for k in $(seq 1 16); do hi=$(( (1100000001 * k / 16) | 1 )); echo "avx2 bounds-large $lo $hi" >> "$JOBLIST"; lo=$((hi + 2)); done
fi
run_job() {   # variant part args...
  local v=$1; shift
  local tag="$v $*"; local f="$OUT/res_$(echo "$tag" | tr ' ' '_').txt"
  if "$OUT/arith_checks_$v" "$@" > "$f" 2>&1; then echo "ok    $tag" >> "$f.status"; else echo "FAIL  $tag" >> "$f.status"; fi
}
export -f run_job; export OUT
rm -f "$OUT"/res_*.txt "$OUT"/res_*.status
start=$SECONDS
xargs -P "$JOBS" -L 1 bash -c 'run_job "$@"' _ < "$JOBLIST"

# ---- report
nfail=$(cat "$OUT"/res_*.status | grep -c '^FAIL')
for f in "$OUT"/res_*.txt; do
  sed "s/^/[$(head -1 "$f.status" | awk '{ $1 = ""; print substr($0, 2) }')] /" "$f" | grep -v '^\[.*\] mont q='
done | sort
echo "arith_checks: $(wc -l < "$JOBLIST") jobs, $nfail failed, $((SECONDS - start)) s"
[ "$nfail" = 0 ]
