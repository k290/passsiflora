// This program is designed to robustly test the super fast pdfilter.  In the GetDecics
// program, so few polynomials pass the polynomial discriminant test, that we are not
// perceptive to small numbers of failures (if they are present).  To get a statistically
// significant number of polynomials for testing, we start with a representative polynomial
// that we know should be retained by the filter.  We apply a random Tschirnhausen
// transformation to the polynomial which creates a new polynomial whose discriminant is
// the same modulo a square.  We do this inside a loop, a fixed number of times (specified
// by the user).  Every one of these transformed polynomials should also pass the pdfilter.
// From the Euclidean chain perspective the polynomials will appear to be completely
// random.
//
// CHANGED after measuring which deliberately injected faults of the filter this test detects (fault_checks.sh):
// the blocks passed to the filter have varying sizes (blockSizes); genPoly() shifts x -> x+t in three quarters of
// the polynomials and sets the x^8 coefficient to 0 in a quarter of them; the seed is printed and can be given.


#include <iostream>
#include <math.h>
#include "pari.h"


using namespace std;

#define BLOCK_SIZE 10000

// The sizes of the blocks passed to the filter: one is chosen at random for each call.  The filter works through
// its input in chunks of 64 polynomials, so a block size that is always 10000 always ends with the same chunk of 16.
static const long blockSizes[] = {1, 7, 63, 64, 65, 127, 255, 1000, 4096, 10000};

#if defined(__x86_64__) || defined(_M_X64)
  #define NF_PDF_DISPATCH 1
  extern "C" void pdf_init_v3(int, const int*);
  extern "C" void pdf_filter_v3(const long long*, int, signed char*);
  extern "C" void pdf_init_base(int, const int*);
  extern "C" void pdf_filter_base(const long long*, int, signed char*);
#else
  extern "C" void pdf_init(int, const int*);
  extern "C" void pdf_filter(const long long*, int, signed char*);
#endif

GEN genPoly();

static void (*nfPdfFilter)(const long long*, int, signed char*) = 0;


// Picks the instruction-set variant of the Fortran filter and initialises it (once).
static void nfPdfSetup(int numPrimes, int *pSet)
{
  if(nfPdfFilter) return;
  void (*pdfInit)(int, const int*);
  #ifdef NF_PDF_DISPATCH
    const char *fb = getenv("NF_FORCE_BASE");
    __builtin_cpu_init();
    if(!(fb && *fb=='1') && __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma")) {
      pdfInit = pdf_init_v3;
      nfPdfFilter = pdf_filter_v3;
      fprintf(stderr, "polDisc pre-filter: AVX2/FMA build.\n");
      }
    else {
      pdfInit = pdf_init_base;
      nfPdfFilter = pdf_filter_base;
      fprintf(stderr, "polDisc pre-filter: baseline build.\n");
      }
  #else
    pdfInit = pdf_init;
    nfPdfFilter = pdf_filter;
  #endif
  pdfInit(numPrimes, pSet);
}




int main(int argc, char **argv)
  {

  if(argc<2 || argc>3)  {
    cout << "Usage: " << argv[0] << " N [seed]\n";
    cout << "where N*" << BLOCK_SIZE << " is the number of polynomials to test.\n";
    exit(1);  }

  long N = atoi(argv[1]);
  unsigned long seed = (argc==3) ? strtoul(argv[2],NULL,10) : (unsigned long)time(NULL);

  pari_init(100000000,10000);  // 100 MB

  // Seed the ransom number generator (both Pari and C)
  setrand(utoi(seed));
  srand(seed);
  cout << "seed = " << seed << "\n";
  
  // Setup the pdf filter.  Hardcode the set S = {2,5}
  int numPrimes = 2;
  int pSet[2] = {2, 5};
  nfPdfSetup(numPrimes, pSet);
  
  // Define the test polynomial:  f = x^10 - 12*x^5 - 4
  //GEN testPoly = mkpoln(11,gen_1,gen_0,gen_0,gen_0,gen_0,stoi(-12),gen_0,gen_0,gen_0,gen_0,stoi(-4));
  //cout << "testPoly = " << GENtostr(testPoly) << "\n";

  // The block of polynomials and the corresponding flag.
  long long polyBlock[11*BLOCK_SIZE];
  char polyFlag[11*BLOCK_SIZE];
  long totalPolysTested = 0;
  long totalFailures = 0;

  long ltop = avma;
  
  long remaining = N*BLOCK_SIZE;   // the number of polynomials still to test

  while(remaining > 0) {

    // Reset the Pari stack for each block.
    avma = ltop;

    // Choose the size of this block at random (but not more than remains to be tested).
    long blockSize = blockSizes[random_Fl(sizeof(blockSizes)/sizeof(blockSizes[0]))];
    if(blockSize > remaining)  blockSize = remaining;

    // Generate a block of random polynomials
    for(long b=0; b<blockSize; ++b) {
      //GEN f = tschirnhaus(testPoly);
      GEN f = genPoly();

      // This is for debug purposes only
      if(0) {
        cout << "f = " << GENtostr(f) << "\n";
	GEN pd = discsr(f);
        cout << "disc(f) = " << GENtostr(pd) << "\n";
        cout << "factored disc = " << GENtostr(factor(pd)) << "\n";
        }

      // I am disabling this. It didn't seem to be needed for tschirnhaus, and definitely not for genPoly().
      // Make sure the new coefficients don't exceed 64 bits.
      // Pari polynomials have 2 code words, so f[2] is the constant and f[12] is the x^10 coeff.
      if(0) {
        bool badPoly = false;
        for(int k=2; k<=12; ++k) {
          if (lgefint((GEN)f[k]) > 3)  { badPoly = true; break; }
          }
        if(badPoly) {
          b = b - 1;
          cout << "Discarding a bad test polynomial:\n";
          cout << "    f = " << GENtostr(f) << "\n";
	  continue;
          }
        }

      // Add the polynomial to the block.  Note that C coeffs ordering is reversed relative to Pari.
      for(int k=0; k<=10; ++k)  polyBlock[b*11+k] = itos((GEN)f[12-k]);
      }

    // Call the pd filter on the whole block
    nfPdfFilter(polyBlock, blockSize, (signed char*)polyFlag);

    // Every polynomial should be good.  Check for failures.  Print failures to screen.
    long accum = 0;  for(long k=0; k<blockSize; ++k)  accum = accum + polyFlag[k];
    if( accum != blockSize ) {
      totalFailures = totalFailures + (blockSize - accum);
      cout << "pdfilter failed.  Failed polynomials:\n";
      for(long b=0; b<blockSize; ++b) {
        if(polyFlag[b]==0) {
	  long g[11];
          for(int k=0; k<=10; ++k)  g[k] = polyBlock[b*11+k];
	  cout << "   x^10 ";
	  for(int k=1; k<=9; ++k) { cout << (g[k]<0 ? "-" : "+") << " " << labs(g[k]) << "*x^" << 10-k << " ";}
          cout << (g[10]<0 ? "-" : "+") << " " << labs(g[10]) << "\n";
          }
        }
      }

    // Increment the test count.
    totalPolysTested = totalPolysTested + blockSize;
    remaining = remaining - blockSize;
  
    } // End loop on blocks

  cout << "Total polynomials tested = " << totalPolysTested << "\n";
  cout << "Total failures = " << totalFailures << "\n";

  return totalFailures ? 1 : 0;
  }



// Generate an irreducible degree 10 polynomial with square discriminant.
// If g(x) is an irreducible quintic, then set f = g(x^2).
// disc(f) = (-1) * 2^10 * (disc(g))^2 * g[0]. So we set g[0] to be a square.
// (The comment had g[0]^5, which has the same square class.)
GEN genPoly()
  {
  GEN f;

  // Keep generating random polynomials until we get an irreducible one.
  do {
    long c[5];
    for(int i=0; i<=4; ++i)  c[i] = rand()%100;

    // A quarter of the polynomials get a zero x^8 coefficient: the Euclidean chain then has the degrees 10, 9, 6, ...
    // as in 25% of the polynomials that pass the discriminant test in the application.
    if( rand()%4 == 0 )  c[4] = 0;

    // Force the constant coeff to be a square mod q for any q we are interested in
    // So it can be an integral square or from the set {-1,2,5}.
    // We give each of {-1,2,5} a 15% chance.
    if( c[0]>=55 && c[0]<70) c[0] = -1;
    else if( c[0]>=70 && c[0]<85) c[0] = 2;
    else if( c[0]>=85) c[0] = 5;
    else c[0] = c[0]*c[0];

    // Generate the polynomial
    f = mkpoln(11,gen_1,gen_0,stoi(c[4]),gen_0,stoi(c[3]),gen_0,stoi(c[2]),gen_0,stoi(c[1]),gen_0,stoi(c[0]));

  } while ( !polisirreducible(f) );

  // Replace x by x+t in three quarters of the polynomials.  The discriminant, and with it the square class, does not
  // change.  The coefficients become signed and large, the odd powers are no longer zero, and the zeros of the
  // Euclidean chain appear by cancellation instead of being zero coefficients.  |t| <= 30 keeps them below 2^63.
  if( rand()%4 != 0 ) {
    long t = 1 + rand()%30;
    if( rand()%2 )  t = -t;
    f = ZX_translate(f, stoi(t));
    }

  return(f);
  }
