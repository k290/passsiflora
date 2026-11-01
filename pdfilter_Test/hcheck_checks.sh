# Could hypothesis (H) be checked at run time instead of assumed?  Builds test_hcheck.cpp with upstream's polDiscTest_cpu.cpp (twice:
# as shipped, and with -fwrapv, renamed polDiscTest_wrapv), GMP and both builds of ../src/pdfilter.f90 (the application's flags), and
# runs it.  No PARI.  Linux.  H_CHECK.md explains the parts and gives the results.
#   UPSTREAM=<GetDecicsSrc> bash hcheck_checks.sh [COUNT [SEED]] [-hits FILE]      COUNT polynomials per family (default 400000),
#                                                                                   about 3 minutes; -hits as in test_hcheck.cpp
# Environment: OUT (default build/hcheck), CXX (g++), FC (gfortran).  Exit status 1 if the check disagrees with test_poldisc.cpp's
# 128-bit replay, the bypass changes a flag of the stock call, or -fwrapv changes a flag.
cd "$(dirname "$0")"
CXX=${CXX:-g++}
FC=${FC:-gfortran}
OUT=${OUT:-build/hcheck}
[ -f "$UPSTREAM/polDiscTest_cpu.cpp" ] || { echo "set UPSTREAM to upstream's GetDecicsSrc directory"; exit 1; }
mkdir -p "$OUT"
$FC -O3 -march=x86-64 -mavx2 -mno-popcnt -I../src -cpp -J "$OUT" -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" "-DPDF_FILTER_NAME='pdf_filter_v3'" -c ../src/pdfilter.f90 -o "$OUT/pdfilter_v3.o" || exit 1
$FC -O3 -march=x86-64 -mtune=generic -I../src -cpp -J "$OUT" -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" "-DPDF_FILTER_NAME='pdf_filter_base'" -c ../src/pdfilter.f90 -o "$OUT/pdfilter_base.o" || exit 1
$CXX -O3 -w -c "$UPSTREAM/polDiscTest_cpu.cpp" -o "$OUT/polDiscTest_cpu.o" || exit 1
$CXX -O3 -w -fwrapv -c "$UPSTREAM/polDiscTest_cpu.cpp" -o "$OUT/polDiscTest_wrapv.o" || exit 1
objcopy --redefine-sym polDiscTest_cpu=polDiscTest_wrapv "$OUT/polDiscTest_wrapv.o" || exit 1
$CXX -O3 -march=x86-64 -mtune=generic -Wall -c test_hcheck.cpp -o "$OUT/test_hcheck.o" || exit 1
$CXX -O3 -o "$OUT/test_hcheck" "$OUT/test_hcheck.o" "$OUT/pdfilter_v3.o" "$OUT/pdfilter_base.o" "$OUT/polDiscTest_cpu.o" "$OUT/polDiscTest_wrapv.o" -lgfortran -lgmp || exit 1
grep -m1 "model name" /proc/cpuinfo 2>/dev/null; $CXX --version | head -1
"$OUT/test_hcheck" "$@"
