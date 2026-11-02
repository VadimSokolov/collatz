/* Bin counts of fractional parts {t*theta}, theta = log2 3, for 0 <= t < R,
   in 2^16 bins, using exact 128 bit fixed point accumulation.
   Usage: br_sum R1 R2 ... (increasing). Prints bounds on B_R = sum 2^{-{t theta}}
   and the count of t with {t theta} > 1 - 2^-30. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
typedef unsigned __int128 u128;
#define NB (1<<16)
int main(int argc,char**argv){
  u128 th = ((u128)HI<<64) | (u128)LO;
  static uint64_t cnt[NB];
  u128 x = 0; uint64_t t = 0; uint64_t near = 0;
  for(int a=1;a<argc;a++){
    uint64_t R = strtoull(argv[a],0,10);
    for(; t<R; t++){
      uint64_t top = (uint64_t)(x>>64);
      cnt[top>>48]++;
      if(top > 0xFFFFFFFC00000000ULL) near++;   /* {t theta} > 1 - 2^-30 */
      x += th;
    }
    long double lo=0, hi=0;
    for(int b=0;b<NB;b++){ if(!cnt[b]) continue;
      hi += cnt[b]*exp2l(-(long double)b/NB);
      lo += cnt[b]*exp2l(-(long double)(b+1)/NB); }
    long double asym = (long double)R/(2*logl(2));
    printf("R=%llu  B_lower=%.6Lf  B_upper=%.6Lf  R/(2ln2)=%.6Lf  ratio_up=%.10Lf  near=%llu\n",
      (unsigned long long)R, lo, hi, asym, hi/R, (unsigned long long)near);
    fflush(stdout);
  }
  return 0;
}
