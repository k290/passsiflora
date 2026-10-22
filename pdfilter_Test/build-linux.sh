# Builds the pre-filter and the test program (Linux).  Run it in this directory:   bash build-linux.sh
#   PARI          where PARI 2.17.3 was built (default: Eric's tree), e.g.   PARI=/opt/pari-2.17.3 bash build-linux.sh
#   FFLAGS_EXTRA  extra gfortran flags for the filter, for experiments, e.g.  FFLAGS_EXTRA=-mfpmath=387 bash build-linux.sh
#                 (x87 arithmetic, as in a 32-bit x86 build; the test should then report failures)
PARI=${PARI:-/home/eric/BOINC_Project/Pari/pari-2.17.3}
INC="-I$PARI/src/headers -I$PARI/Olinux-x86_64"

# Build the fortran pdfilter.  chain_general.inc, which pdfilter.f90 includes, is not in this directory: -I../src finds it there
# (a copy next to pdfilter.f90 would be used first).
gfortran -O3 -march=x86-64-v3 $FFLAGS_EXTRA -I../src -cpp -DPDF_MOD=pdfilter_v3 "-DPDF_INIT_NAME='pdf_init_v3'" "-DPDF_FILTER_NAME='pdf_filter_v3'" -c pdfilter.f90 -o pdfilter_v3.o
gfortran -O3 -march=x86-64 -mtune=generic $FFLAGS_EXTRA -I../src -cpp -DPDF_MOD=pdfilter_base "-DPDF_INIT_NAME='pdf_init_base'" "-DPDF_FILTER_NAME='pdf_filter_base'" -c pdfilter.f90 -o pdfilter_base.o

# Build the test program
g++ -O3 -w $INC -c -o test_pdfilter.o test_pdfilter.cpp
g++ -O3 -static -static-libgcc -static-libstdc++ -o test_pdfilter test_pdfilter.o pdfilter_v3.o pdfilter_base.o -L$PARI/Olinux-x86_64 -l:libpari.a -l:libgfortran.a -lquadmath -lm -l:libgmp.a
