/* glibc's printf rounds the decimal digits it prints with the CURRENT floating-point rounding mode.
   This is why the application's output file differs from the truth file in the last digit of one percentage when the
   stub OpenCL library has put the thread into round-up mode:  93/109*100 = 85.3211009174...  prints as 85.3212.
     gcc -o printf_rounding printf_rounding.c -lm && ./printf_rounding
   Output on glibc (Ubuntu 26.04):
     nearest  85.3211  (value 85.321100917431)
     up       85.3212  (value 85.321100917432)
     down     85.3211  (value 85.321100917431)
     zero     85.3211  (value 85.321100917431)                                                                    */
#include <stdio.h>
#include <fenv.h>

int main(void)
{
  volatile double x = 9300.0 / 109.0;      /* computed once, in round-to-nearest */
  const int mode[4] = { FE_TONEAREST, FE_UPWARD, FE_DOWNWARD, FE_TOWARDZERO };
  const char *name[4] = { "nearest", "up", "down", "zero" };
  for (int i = 0; i < 4; i++) {
    fesetround(mode[i]);
    printf("%-8s %.4f  (value %.12f)\n", name[i], x, x);
  }
  return 0;
}
