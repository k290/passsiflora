/* Android on arm64 needs the program's thread-local storage (TLS) segment aligned to at least 64 bytes: bionic keeps its
   own thread-local slots in the first 64 bytes after the thread pointer, and Android's linker refuses a program whose
   segment is aligned less.  The Fortran runtime from 0-gfortran.sh has thread-local variables aligned to 8; this one
   variable, compiled with -fno-emulated-tls (4-app.sh), raises the segment's alignment to 64.  Native TLS in a program
   needs Android 10 or newer. */
__attribute__((aligned(64), used)) __thread char nf_android_tls_pad;
