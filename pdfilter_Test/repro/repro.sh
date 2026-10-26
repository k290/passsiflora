#!/bin/bash
# repro.sh -- runs every step of the fault guide (the section "When the floating-point assumptions fail" of pdfilter_chains.tex)
# on Linux (Ubuntu, or WSL), and compares what each step prints with the verbatim output in the guide (expected/<part>.txt).
#
#   R=<this repository> W=<an empty work directory> bash repro.sh <part>...       exit status 1 if any output differs
#
# parts, in the order of the guide:
#   sources   upstream sparse checkout at 638ca4e (only if $W/upstream is absent) + the repository's build/prepare-src.sh,
#             with the floating-point filter of commit 327b04f (the filter the guide is about)
#   check     the source tree is as expected (commit, 31 work units, the shipped PARI untouched)
#   build     recipe step 3: good filter, application objects, feeder, programs, stub OpenCL library
#   step1     red_demo and signtest (the directed-rounding lemma at the 30 filter primes, and the x87 lemma)
#   step2     the repository's test script on two work units; R1-R8 through the good pipeline
#   step3     the x87 filter (rebuild, relink); R1-R8; batch sizes; NF_VERIFY
#   inst      instrumented temporary copy of the filter (build) and the lines it inserts
#   step3i    inside the x87 filter: R6
#   step4     the stub OpenCL library: S1-S8, their unshifted forms, T1-T3, the filter alone, inside S1
#   step5     round down: junk polynomials
#   step6     complete work units (good, x87, library modes; about 15 minutes) and the printf effect
#   step7     the filter alone over the 83 real results and their 996 shifted copies
#   appendix  the filter alone over 20,000 generated valid and 5,000 junk polynomials (about 10 minutes)
#   all       sources check build step1 step2 step3 inst step3i step4 step5 step6 step7 appendix
# What is compared: the transcript ("$ command" followed by its output) with expected/<part>.txt, after replacing the paths of
# $W and $R by those names and eliding what varies between runs: the numbers of the application's NF_TIMING line (if one is
# shown at all), the wall=<s> of build/test-suite-unix.sh and the "took <s> s" of step 6.  $W/transcripts/<part>.raw keeps the
# transcript before eliding.  REPRO_ACCEPT=1 records the new transcripts as the expected output instead of comparing.
# The last line is a summary: how many parts match the guide, and which differ.
set -u
R=${R:-$(cd "$(dirname "$0")/../.." && pwd)}
W=${W:?set W to a work directory}
export R W
RP=$R/pdfilter_Test/repro
EXP=$RP/expected
mkdir -p "$W" "$W/transcripts"
cd "$W" || exit 2
TR=$W/transcripts
status=0; n_ok=0; n_rec=0; differ=""

# the variables and helpers that the guide defines (recipe step 3 and guide step 2); the transcripts define them again
P=$W/upstream/Pari/pari-2.8-1711-ge5c317c; PO=$P/Olinux-x86_64; S=$W/build_src/GetDecicsSrc
INC="-I$S -I$RP/stubs -I$P/src/headers -I$PO"
LIBS="-L$PO -l:libpari.a -lgfortran -lquadmath -l:libgmp.a -lm -lpthread -ldl"
APP="obj/TgtMartinet.o obj/polDiscTest_cpu.o obj/nf_gpu_stub.o"
DF=$W/upstream/Test/dataFiles
WUS="sf5_DS-14x11_Grp33043of705883 sf5_DS-14x11_Grp10973of705883 sf4_DS-13x271-25_Grp98of200000 sf3_DS-16x270_Grp74of3932160 sf3_DS-15x271_Grp47of6553600 sf3_DS-16x271-1_Grp99of2000000 sf3_DS-16x271-25_Grp79of2000000"
export P PO S RP INC LIBS APP DF WUS
run()  { "$@" 2>&1 | grep -v '^NF_TIMING'; }
each() { i=0; while read -r p; do i=$((i+1)); echo "$p" > one.txt; echo "$1$i: $("$2" one.txt out.txt 2>/dev/null), result lines: $(wc -l < out.txt)"; done < "$3"; }
stats() { d=$(mktemp -d); cp $DF/$1.dat $d/in; ( cd $d && $2 > so.txt 2>/dev/null; echo "== $1"; grep StatVec so.txt; diff <(grep -v 'Elapsed Time' out) <(grep -v 'Elapsed Time' $DF/wu_$1.out.truth) ); rm -rf $d; }
ocl()  { d=$(mktemp -d); cp $DF/$1.dat $d/in; ( cd $d && s=$SECONDS; PROXY_MODE=$2 OCL_LIB=$W/libOpenCL.so.1 $W/cpu_ocl_good > so.txt 2>/dev/null; echo "== PROXY_MODE=$2: $(grep StatVec so.txt), took $((SECONDS-s)) s"; diff <(grep -v 'Elapsed Time' out) <(grep -v 'Elapsed Time' $DF/wu_$1.out.truth) ); rm -rf $d; }
export -f run each stats ocl

elide() {
  sed -e "s|$W|\$W|g" -e "s|$R|\$R|g" -e 's/wall=[0-9]* s/wall=... s/' -e 's/took [0-9]* s/took ... s/' \
      -e 's/^NF_TIMING total=.*/NF_TIMING total=... s polDisc=... s other=... s/'
}
# transcript <part> : runs the command lines read from stdin (one command per line, run in this order in one shell) and
# records "$ command" and the output; then compares the record with expected/<part>.txt
transcript() {
  local part=$1 line
  local -a cmds=()
  while IFS= read -r line; do [ -n "$line" ] && cmds+=("$line"); done
  { cd "$W"; for line in "${cmds[@]}"; do printf '$ %s\n' "$line"; eval "$line" < /dev/null 2>&1; done; } | tee "$TR/$part.raw" | elide > "$TR/$part.txt"
  if [ -n "${REPRO_ACCEPT:-}" ]; then mkdir -p "$EXP"; cp "$TR/$part.txt" "$EXP/$part.txt"; echo "$part: transcript recorded as the expected output"; n_rec=$((n_rec+1))
  elif [ ! -f "$EXP/$part.txt" ]; then echo "$part: no expected output (transcript in $TR/$part.txt)"; status=1; differ="$differ $part"
  elif diff "$EXP/$part.txt" "$TR/$part.txt" > "$TR/$part.diff"; then echo "$part: output as in the guide"; n_ok=$((n_ok+1))
  else echo "$part: OUTPUT DIFFERS from the guide (see $TR/$part.diff)"; head -20 "$TR/$part.diff"; status=1; differ="$differ $part"; fi
}
summary() {
  if [ -n "${REPRO_ACCEPT:-}" ]; then echo "summary: $n_rec part(s) recorded as the expected output"
  else echo "summary: $n_ok part(s) as in the guide, $(echo $differ | wc -w) differ${differ:+:$differ}"; fi
}

part_sources() {
  if [ -d "$W/upstream/.git" ]; then echo "sources: $W/upstream exists, not cloned again"; else
    git clone --filter=blob:none --no-checkout https://github.com/drivere/get-decics-numberfields upstream &&
    ( cd upstream && git sparse-checkout init --cone &&
      git sparse-checkout set GetDecicsSrc Test Pari/pari-2.8-1711-ge5c317c && git checkout 638ca4e ) || { status=1; return; }
  fi
  # a checkout interrupted by the network leaves files missing: check before going on
  [ -f $W/upstream/Test/dataFiles/sf5_DS-14x11_Grp33043of705883.dat ] && [ -f $W/upstream/Pari/pari-2.8-1711-ge5c317c/Olinux-x86_64/libpari.a ] ||
    { echo "sources: the checkout is incomplete (run 'git checkout -f 638ca4e' in $W/upstream again)"; status=1; return; }
  UPSTREAM=$W/upstream PASSIFLORA=$R OUT=$W/build_src bash $R/build/prepare-src.sh || status=1
  # The guide is about the floating-point filter of commit 327b04f.  The filter of later commits is integer-only, and the faults of
  # the guide do not touch it; so the filter sources are taken from that commit, and the transcripts stay those of the guide.
  for f in pdfilter.f90 chain_general.inc; do git -C "$R" show 327b04f:src/$f > $W/build_src/GetDecicsSrc/$f || status=1; done
}

part_check() { transcript check <<'EOF'
git -C upstream log --oneline -1
git -C upstream status --short | wc -l
ls upstream/Test/dataFiles | grep -c '\.dat$'
ls upstream/Pari/pari-2.8-1711-ge5c317c/Olinux-x86_64 | grep -E '^(libpari\.a|paricfg\.h)$'
EOF
}

part_build() { transcript build <<'EOF'
P=$W/upstream/Pari/pari-2.8-1711-ge5c317c; PO=$P/Olinux-x86_64; S=$W/build_src/GetDecicsSrc; RP=$R/pdfilter_Test/repro
: > $S/config.h
INC="-I$S -I$RP/stubs -I$P/src/headers -I$PO"
mkdir -p good obj
gfortran -O3 -march=x86-64-v3 -cpp -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" "-DPDF_FILTER_NAME='pdf_filter_v3'" -I$S -J good -c $S/pdfilter.f90 -o good/pdfilter_v3.o
gfortran -O3 -march=x86-64 -mtune=generic -cpp -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" "-DPDF_FILTER_NAME='pdf_filter_base'" -I$S -J good -c $S/pdfilter.f90 -o good/pdfilter_base.o
for f in polDiscTest_cpu TgtMartinet; do g++ -O3 -w $INC -DAPP_VERSION_CPU_STD -c $S/$f.cpp -o obj/$f.o; done
g++ -O2 -w $INC -c $S/nf_gpu_stub.cpp -o obj/nf_gpu_stub.o
for f in app_feed main_ocl stubs/main; do g++ -O3 -w $INC -DAPP_VERSION_CPU_STD -c $RP/$f.cpp -o obj/$(basename $f).o; done
LIBS="-L$PO -l:libpari.a -lgfortran -lquadmath -l:libgmp.a -lm -lpthread -ldl"; APP="obj/TgtMartinet.o obj/polDiscTest_cpu.o obj/nf_gpu_stub.o"
g++ -O3 -o app_feed_good obj/app_feed.o $APP good/pdfilter_v3.o good/pdfilter_base.o $LIBS
g++ -O3 -o cpu_good obj/main.o $APP good/pdfilter_v3.o good/pdfilter_base.o $LIBS
g++ -O3 -o cpu_ocl_good obj/main_ocl.o $APP good/pdfilter_v3.o good/pdfilter_base.o $LIBS
g++ -O2 -o pdf_check_good $RP/pdf_check.cpp good/pdfilter_v3.o good/pdfilter_base.o -lgfortran -ldl
gcc -shared -fPIC -O2 -o libOpenCL.so.1 $RP/stub_opencl.c
cp $RP/R.txt $RP/S.txt $RP/S0.txt $RP/T.txt $RP/junk8.txt .
ls app_feed_good cpu_good cpu_ocl_good pdf_check_good libOpenCL.so.1
EOF
}

part_step1() { transcript step1 <<'EOF'
gcc -O0 -o red_demo $RP/red_demo.c -lm && ./red_demo
gcc -O0 -o signtest $RP/signtest.c -lm && ./signtest
EOF
}

part_step2() { transcript step2 <<'EOF'
DF=$W/upstream/Test/dataFiles
bash $R/build/test-suite-unix.sh $W/cpu_good $DF 1 sf5_DS-14x11_Grp33043of705883 sf5_DS-14x11_Grp10973of705883
run() { "$@" 2>&1 | grep -v '^NF_TIMING'; }
each() { i=0; while read -r p; do i=$((i+1)); echo "$p" > one.txt; echo "$1$i: $("$2" one.txt out.txt 2>/dev/null), result lines: $(wc -l < out.txt)"; done < "$3"; }
each R ./app_feed_good R.txt
sed -n 6p R.txt > R6.txt; run ./app_feed_good R6.txt out.txt; cat out.txt
EOF
}

part_step3() { transcript step3 <<'EOF'
mkdir -p x87
gfortran -O3 -march=x86-64-v3 -mfpmath=387 -cpp -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" "-DPDF_FILTER_NAME='pdf_filter_v3'" -I$S -J x87 -c $S/pdfilter.f90 -o x87/pdfilter_v3.o
gfortran -O3 -march=x86-64 -mtune=generic -mfpmath=387 -cpp -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" "-DPDF_FILTER_NAME='pdf_filter_base'" -I$S -J x87 -c $S/pdfilter.f90 -o x87/pdfilter_base.o
g++ -O3 -o app_feed_x87 obj/app_feed.o $APP x87/pdfilter_v3.o x87/pdfilter_base.o $LIBS
g++ -O3 -o cpu_x87 obj/main.o $APP x87/pdfilter_v3.o x87/pdfilter_base.o $LIBS
g++ -O2 -o pdf_check_x87 $RP/pdf_check.cpp x87/pdfilter_v3.o x87/pdfilter_base.o -lgfortran -ldl
each R ./app_feed_x87 R.txt
run ./app_feed_x87 R.txt out.txt; echo "result lines: $(wc -l < out.txt)"
FEED_ALONE=1 run ./app_feed_x87 R.txt out.txt; echo "result lines: $(wc -l < out.txt)"
NF_VERIFY=1 ./app_feed_x87 R6.txt out.txt 2>&1 | grep -v '^NF_TIMING'; echo "exit status ${PIPESTATUS[0]}"
for n in 1 2 3 4 5 6 7 8; do head -n $n R.txt > b.txt; echo "R1-R$n in one batch: $(./app_feed_x87 b.txt out.txt 2>/dev/null), baseline build: $(NF_FORCE_BASE=1 ./app_feed_x87 b.txt out.txt 2>/dev/null)"; done
for o in good/pdfilter_v3.o x87/pdfilter_v3.o; do echo "== x87 instructions per function in $o"; objdump -d --no-show-raw-insn $o | awk '/^[0-9a-f]+ </{f=$2} /\tf[a-z]+ /{n[f]++} END{for(k in n) print n[k], k}' | sort -k2; done
EOF
}

# instrumented temporary copy of the filter: what lane 1 sees, written to standard error
part_inst() {
  rm -rf $W/inst && mkdir -p $W/inst/good $W/inst/x87 || { status=1; return; }
  cp $S/pdfilter.f90 $S/chain_general.inc $W/inst/
  ( cd $W/inst
    # the Euler loops (EU): first squaring at which base(1) is an integer, acc(1) then, final acc(1)
    EUW="write(0, '(a,i3,i4,2f22.6)') 'EU  j=', j, first_, acc_at_, real(acc(1), 8)"
    EU1='s|^    integer :: m, k, l, ia, ib, ir, tmp$|&\n    integer :: it_, first_\n    real(8) :: acc_at_|'
    EU2='s|^    e = q\(s\)\{0,1\}exp(j)$|&\n    it_ = 0; first_ = -1; acc_at_ = 0|'
    EU3='s|^      if (iand(e, 1_i8) == 1) acc = red\(s\)\{0,1\}(acc\*base, q, qi)$|      if (first_ < 0 .and. base(1) == aint(base(1))) then\n        first_ = it_; acc_at_ = acc(1)\n      end if\n&|'
    EU4='s|^      e = shiftr(e, 1)$|&; it_ = it_ + 1|'
    sed -i -e "s|^\( *\)a(:, 10) = 1.0_sp|&\n\1if (c == 1) write(0, '(a,i3,10f9.4)') 'LDS j=', j, a(1, 0:9)|" \
           -e "s|^\( *\)a(:, 10) = 1.0_dp|&\n\1if (c == 1) write(0, '(a,i3,10f14.4)') 'LDD j=', j, a(1, 0:9)|" \
           -e "s|^\( *\)keep = .not. (ok .and. acc == -1.0_sp)|&\n\1$EUW\n\1write(0, '(a,i3,f18.6,f12.6,l2)') 'CBS j=', j, p(1, 0, ib), acc(1), ok(1)|" \
           -e "s|^\( *\)keep = .not. (ok .and. acc == -1.0_dp)|&\n\1$EUW\n\1write(0, '(a,i3,f18.6,f12.6,l2)') 'CBD j=', j, p(1, 0, ib), acc(1), ok(1)|" \
           -e "/^  subroutine chain_block/,/^  end subroutine chain_block/{" -e "$EU1" -e "$EU2" -e "$EU3" -e "$EU4" -e "}" pdfilter.f90
    sed -i -e "s|^\( *\)mm = merge(tmp, mm, act > 0 .and. ev - 2\*aint(real(0.5, GK)\*ev) /= 0)|&\n\1write(0, '(a,i3,a,i3,3i4,2f14.4)') 'CGIT j=', j, ' it=', it, int(du(1)), int(dv(1)), int(dw(1)), real(mm(1), 8), real(lv(1), 8)|" \
           -e "s|^\( *\)keep = .not. (done > 0 .and. resz == 0 .and. acc == -1)|&\n\1$EUW\n\1write(0, '(a,i3,2f18.6,2f4.0)') 'CG  j=', j, real(mm(1), 8), real(acc(1), 8), real(done(1), 8), real(resz(1), 8)|" \
           -e 's|^    real(GK) :: du(BLK), dv(BLK), dw(BLK), mm(BLK), lu(BLK), lv(BLK), ev(BLK)$|&\n    integer :: it_, first_\n    real(8) :: acc_at_|' \
           -e 's|^    e = GEXP(j)$|&\n    it_ = 0; first_ = -1; acc_at_ = 0|' \
           -e 's|^      if (iand(e, 1_i8) == 1) acc = GRED(acc\*base, q, qi)$|      if (first_ < 0 .and. base(1) == aint(base(1))) then\n        first_ = it_; acc_at_ = acc(1)\n      end if\n&|' \
           -e "$EU4" chain_general.inc
    [ "$(diff $S/pdfilter.f90 pdfilter.f90 | grep -c '^>')" = 20 ] && [ "$(diff $S/chain_general.inc chain_general.inc | grep -c '^>')" = 10 ] || { echo "inst: the sed expressions did not apply"; exit 1; }
    for v in good x87; do
      X=""; [ $v = x87 ] && X=-mfpmath=387
      gfortran -O3 -march=x86-64-v3 $X -ffree-line-length-none -cpp -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" "-DPDF_FILTER_NAME='pdf_filter_v3'" -I. -J $v -c pdfilter.f90 -o $v/pdfilter_v3.o &&
      gfortran -O3 -march=x86-64 -mtune=generic $X -ffree-line-length-none -cpp -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" "-DPDF_FILTER_NAME='pdf_filter_base'" -I. -J $v -c pdfilter.f90 -o $v/pdfilter_base.o &&
      ( cd $W && g++ -O3 -o inst/app_feed_$v obj/app_feed.o $APP inst/$v/pdfilter_v3.o inst/$v/pdfilter_base.o $LIBS ) || exit 1
    done ) || { status=1; return; }
  transcript inst <<'EOF'
diff $S/pdfilter.f90 inst/pdfilter.f90
diff $S/chain_general.inc inst/chain_general.inc
EOF
}

part_step3i() { transcript step3i <<'EOF'
NF_PDFDEBUG=1 inst/app_feed_good R6.txt out.txt 2>&1 | grep -E '^(LDS j=  1|LDD j= 19|CBD j= 19|pdf prime 2[45])'
NF_PDFDEBUG=1 inst/app_feed_x87 R6.txt out.txt 2>&1 | grep -E '^(LDS j=  1|LDD j= 19|EU  j= 19|CBD j= 1[789]|pdf prime 2[45])'
NF_FORCE_BASE=1 NF_PDFDEBUG=1 inst/app_feed_x87 R6.txt out.txt 2>&1 | grep -E '^(CBD j= 19|pdf prime 2[45])'
echo 'q=1009; printf("%.4f\n", apply(r->r-q/2048*round(2048*r/q), [429,223,-413,-268,-100])); q=33552361; printf("%.4f\n", apply(r->r-q/2048*round(2048*r/q), [625000,1390625,662500]))' | gp -q
rm -f polys.txt; gp -q $RP/genpolys.gp < /dev/null 2> /dev/null; head -n 2000 polys.txt > p2000.txt
FEED_ALONE=1 run ./app_feed_good p2000.txt out.txt
FEED_ALONE=1 run ./app_feed_x87 p2000.txt out.txt
FEED_ALONE=1 inst/app_feed_x87 p2000.txt out.txt 2> trace.txt; awk -f $RP/classify.awk trace.txt
EOF
}

part_step4() { transcript step4 <<'EOF'
for m in none up down zero; do echo "== PROXY_MODE=$m"; PROXY_MODE=$m OCL_LIB=./libOpenCL.so.1 run ./app_feed_good S.txt out.txt; echo "result lines: $(wc -l < out.txt)"; done
PROXY_MODE=up OCL_LIB=./libOpenCL.so.1 run ./app_feed_good S0.txt out.txt; echo "result lines: $(wc -l < out.txt)"
PROXY_MODE=up OCL_LIB=./libOpenCL.so.1 each S ./app_feed_good S.txt
PROXY_MODE=up OCL_LIB=./libOpenCL.so.1 NF_VERIFY=1 ./app_feed_good S.txt out.txt 2>&1 | grep -v '^NF_TIMING'; echo "exit status ${PIPESTATUS[0]}"
for m in none up; do echo "== PROXY_MODE=$m"; PROXY_MODE=$m OCL_LIB=./libOpenCL.so.1 run ./app_feed_good T.txt out.txt; echo "result lines: $(wc -l < out.txt)"; done
./pdf_check_good alone < T.txt 2>&1 > /dev/null
./pdf_check_good alone up < T.txt 2>&1 > /dev/null
PROXY_MODE=up ./pdf_check_good alone ocl < T.txt 2>&1 > /dev/null
PROXY_MODE=up ./pdf_check_good alone ocl-late < T.txt 2>&1 > /dev/null
sed -n 1p S.txt > S1.txt
inst/app_feed_good S1.txt out.txt 2>&1 | grep '^CGIT j=  1' | head -n 3
PROXY_MODE=up OCL_LIB=./libOpenCL.so.1 inst/app_feed_good S1.txt out.txt 2>&1 | grep '^CGIT j=  1' | head -n 3
PROXY_MODE=up OCL_LIB=./libOpenCL.so.1 NF_PDFDEBUG=1 inst/app_feed_good S1.txt out.txt 2>&1 | grep '^pdf prime'
EOF
}

part_step5() { transcript step5 <<'EOF'
NF_PDFDEBUG=1 run ./app_feed_good junk8.txt out.txt
PROXY_MODE=down OCL_LIB=./libOpenCL.so.1 NF_PDFDEBUG=1 run ./app_feed_good junk8.txt out.txt | grep -v 'kept *8$'
PROXY_MODE=down OCL_LIB=./libOpenCL.so.1 NF_PDFDEBUG=1 run ./app_feed_good junk8.txt out.txt | grep -c 'kept *8$'
PROXY_MODE=zero OCL_LIB=./libOpenCL.so.1 NF_PDFDEBUG=1 run ./app_feed_good junk8.txt out.txt | grep -c 'kept *8$'
EOF
}

part_step6() { transcript step6 <<'EOF'
bash $R/build/test-suite-unix.sh $W/cpu_good $DF 7 $WUS | sort
bash $R/build/test-suite-unix.sh $W/cpu_x87 $DF 7 $WUS | sort
stats() { d=$(mktemp -d); cp $DF/$1.dat $d/in; ( cd $d && $2 > so.txt 2>/dev/null; echo "== $1"; grep StatVec so.txt; diff <(grep -v 'Elapsed Time' out) <(grep -v 'Elapsed Time' $DF/wu_$1.out.truth) ); rm -rf $d; }
for wu in sf5_DS-14x11_Grp33043of705883 sf5_DS-14x11_Grp10973of705883 sf4_DS-13x271-25_Grp98of200000 sf3_DS-16x270_Grp74of3932160; do stats $wu $W/cpu_x87 > stats_$wu.txt & done; wait; cat stats_*.txt
for m in none up; do PROXY_MODE=$m OCL_LIB=$W/libOpenCL.so.1 bash $R/build/test-suite-unix.sh $W/cpu_ocl_good $DF 1 sf5_DS-14x11_Grp33043of705883 | head -n 1; done
gcc -o printf_rounding $RP/printf_rounding.c -lm && ./printf_rounding
ocl() { d=$(mktemp -d); cp $DF/$1.dat $d/in; ( cd $d && s=$SECONDS; PROXY_MODE=$2 OCL_LIB=$W/libOpenCL.so.1 $W/cpu_ocl_good > so.txt 2>/dev/null; echo "== PROXY_MODE=$2: $(grep StatVec so.txt), took $((SECONDS-s)) s"; diff <(grep -v 'Elapsed Time' out) <(grep -v 'Elapsed Time' $DF/wu_$1.out.truth) ); rm -rf $d; }
for m in none up down zero; do ocl sf5_DS-14x11_Grp33043of705883 $m > ocl_$m.txt & done; wait; cat ocl_none.txt ocl_up.txt ocl_down.txt ocl_zero.txt
EOF
}

part_step7() { transcript step7 <<'EOF'
grep -h '^x' $DF/wu_*.out.truth > hits.txt; rm -f hits_vec.txt hits_shift_vec.txt; gp -q $RP/hits.gp < /dev/null 2> /dev/null
./pdf_check_good alone < hits_vec.txt 2>&1 > /dev/null
PROXY_MODE=up ./pdf_check_good alone ocl < hits_vec.txt 2>&1 > /dev/null
for b in 0 1; do bad=""; i=0; while read -r p; do i=$((i+1)); r=$(echo "$p" | NF_FORCE_BASE=$b ./pdf_check_x87 alone 2>/dev/null | cut -d' ' -f2); [ "$r" = 0 ] && bad="$bad $i"; done < hits_vec.txt; echo "NF_FORCE_BASE=$b: $(echo $bad | wc -w) rejected:$bad"; done
./pdf_check_good alone < hits_shift_vec.txt 2>&1 > /dev/null
./pdf_check_good alone up < hits_shift_vec.txt 2>&1 > /dev/null
PROXY_MODE=up ./pdf_check_good alone ocl < hits_shift_vec.txt 2>&1 > /dev/null
EOF
}

# fresh process per polynomial over a file; prints how many were rejected
fresh() { n=0; while read -r p; do [ "$(echo "$p" | "$@" 2>/dev/null | cut -d' ' -f2)" = 0 ] && n=$((n+1)); done; echo $n; }
export -f fresh
part_appendix() {
  if [ ! -f $W/polys.txt ]; then ( cd $W && gp -q $RP/genpolys.gp < /dev/null ) || { status=1; return; }; fi
  if [ ! -f $W/junk.txt ];  then ( cd $W && gp -q $RP/genjunk.gp < /dev/null ) || { status=1; return; }; fi
  transcript appendix <<'EOF'
wc -l < polys.txt; wc -l < junk.txt
fresh() { n=0; while read -r p; do [ "$(echo "$p" | "$@" 2>/dev/null | cut -d' ' -f2)" = 0 ] && n=$((n+1)); done; echo $n; }
echo "good:          $(fresh ./pdf_check_good alone < polys.txt)"
echo "x87 AVX2:      $(NF_FORCE_BASE=0 fresh ./pdf_check_x87 alone < polys.txt)"
echo "x87 baseline:  $(NF_FORCE_BASE=1 fresh ./pdf_check_x87 alone < polys.txt)"
echo "up, shim:      $(fresh ./pdf_check_good alone up < polys.txt)"
echo "up, lib first: $(PROXY_MODE=up fresh ./pdf_check_good alone ocl < polys.txt)"
echo "down, shim:    $(fresh ./pdf_check_good alone down < polys.txt)"
echo "zero, shim:    $(fresh ./pdf_check_good alone zero < polys.txt)"
echo "junk, good:    $(./pdf_check_good alone < junk.txt 2>&1 > /dev/null)"
echo "junk, down:    $(PROXY_MODE=down ./pdf_check_good alone ocl < junk.txt 2>&1 > /dev/null)"
echo "junk, zero:    $(PROXY_MODE=zero ./pdf_check_good alone ocl < junk.txt 2>&1 > /dev/null)"
EOF
}

[ $# -eq 0 ] && { sed -n '2,24p' "$0"; exit 2; }
for p in "$@"; do
  [ "$p" = all ] && set -- sources check build step1 step2 step3 inst step3i step4 step5 step6 step7 appendix
done
for p in "$@"; do
  case $p in
    sources|check|build|step1|step2|step3|inst|step3i|step4|step5|step6|step7|appendix) part_$p ;;
    *) echo "unknown part: $p"; status=2 ;;
  esac
  [ $p = sources ] && [ $status -ne 0 ] && { echo "sources failed: stopping"; summary; exit 1; }
  [ $p = build ] && [ ! -x $W/app_feed_good ] && { echo "build failed: stopping"; summary; exit 1; }
done
summary
exit $status
