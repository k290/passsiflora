# Upstream's exact test against Definition "Passing" (docs/whitepaper.tex, Proposition "What the exact test computes"): builds
# test_poldisc.cpp with build-linux.sh, linked with upstream's polDiscTest_cpu.cpp and GMP, and runs every family of polynomials of
# test_poldisc.cpp in parallel chunks.  Linux.
#   UPSTREAM=<GetDecicsSrc> CAND=<candidates> HITS=<hits> bash poldisc_checks.sh      about 10.7 million polynomials
#   QUICK=1 ...                                                                       about 0.3 million
#   S=2,11 ...                                                                        another set S (one to three primes; the cand and
#                                                                                     hits files, polynomials of S = {2,5}, are not used)
# polDiscTest_cpu.cpp is compiled with -fwrapv, as the application compiles it.
# UPSTREAM is upstream's source directory, with the patch applied (the patch changes only polDiscTest_cpu.cpp's conversions for a
# 32-bit long, which do not matter on 64-bit Linux), CAND a file of
# captured candidates and HITS a file of results, one polynomial per line, leading coefficient first (Appendix "Reproducing the
# tables").  Environment: OUT (default poldisc_build), JOBS (default: the number of CPUs), PARI as for build-linux.sh.
# Exit status 1 if upstream and the definition disagree on a polynomial where upstream's int64 steps do not overflow, or the filter
# rejects a polynomial that the definition passes, or one that upstream passes without an overflow.  rej-up-ovf counts the filter's
# rejections of polynomials that upstream passes only because its int64 steps overflow (family wrap); they do not fail the run.
cd "$(dirname "$0")"
OUT=${OUT:-poldisc_build}
JOBS=${JOBS:-$(nproc 2>/dev/null || echo 4)}
[ -f "$UPSTREAM/polDiscTest_cpu.cpp" ] || { echo "set UPSTREAM to upstream's GetDecicsSrc directory"; exit 1; }
mkdir -p "$OUT"
OUT="$OUT" EXTRA_CXX="$UPSTREAM/polDiscTest_cpu.cpp" EXTRA_CXXFLAGS=-fwrapv LIBS_EXTRA=-lgmp TEST=test_poldisc.cpp bash build-linux.sh > "$OUT/build.log" 2>&1 || { cat "$OUT/build.log"; exit 1; }
# jobs: "family count seed [file options]"
JOBLIST="$OUT/jobs.txt"; : > "$JOBLIST"
chunks() {   # family total nchunks
  local c; for c in $(seq 1 $3); do echo "$1 $(( $2 * c / $3 - $2 * (c - 1) / $3 )) $c${S:+ -S $S}" >> "$JOBLIST"; done
}
if [ -n "$QUICK" ]; then
  for f in rand3 rand100 rand1e6 pass passbig zero; do chunks $f 50000 2; done
  chunks wrap 200 1
else
  chunks rand3 3000000 16
  for f in rand100 pass; do chunks $f 2000000 8; done
  for f in rand1e6 passbig zero; do chunks $f 1000000 8; done
  chunks wrap 2000 4
fi
[ -z "$S" ] && [ -n "$CAND" ] && for c in $(seq 1 16); do echo "cand 0 $c -cand $CAND -part $c/16" >> "$JOBLIST"; done
[ -z "$S" ] && [ -n "$HITS" ] && echo "hits 5 1 -hits $HITS" >> "$JOBLIST"
run_job() {
  local f="$OUT/res_$(echo "$1_$2_$3" | tr ' /' '__').txt"
  if "$OUT/test_poldisc" "$@" > "$f" 2>&1; then echo ok > "$f.status"; else echo FAIL > "$f.status"; fi
}
export -f run_job; export OUT
rm -f "$OUT"/res_*.txt "$OUT"/res_*.status
start=$SECONDS
xargs -P "$JOBS" -L 1 bash -c 'run_job "$@"' _ < "$JOBLIST"
nfail=$(cat "$OUT"/res_*.status | grep -c FAIL)
cat "$OUT"/res_*.txt | awk '
  $1 == "POLDISC" { if (!($2 in N)) { N[$2] = 1; F[++n] = $2 }; for (i = 3; i <= 13; i++) V[$2, i] += $i }
  END {
    printf "%-8s %11s %9s %9s %10s %9s %9s %11s %9s %9s %6s %10s\n", "family", "polys", "upstream", "definit.", "int64-ovf", "diff-ok", "diff-ovf", "filter-rej", "rej-def", "rej-up", "v3/base", "rej-up-ovf"
    for (j = 1; j <= n; j++) { f = F[j]; printf "%-8s", f; for (i = 3; i <= 13; i++) { printf " %*.0f", (i == 3 ? 11 : i == 9 ? 11 : i == 6 ? 10 : i == 12 ? 6 : i == 13 ? 10 : 9), V[f, i]; T[i] += V[f, i] }; printf "\n" }
    printf "%-8s", "all"; for (i = 3; i <= 13; i++) printf " %*.0f", (i == 3 ? 11 : i == 9 ? 11 : i == 6 ? 10 : i == 12 ? 6 : i == 13 ? 10 : 9), T[i]; printf "\n"
  }'
grep -h "DISAGREE\|leaves the int64" "$OUT"/res_*.txt | head -5
echo "poldisc_checks${S:+ (S = $S)}: $(wc -l < "$JOBLIST") jobs, $nfail failed, $((SECONDS - start)) s"
[ "$nfail" = 0 ]
