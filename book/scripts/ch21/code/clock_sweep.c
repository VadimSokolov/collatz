// Sweep of the clock tau(n) for n <= N, with 128 bit costs (c(v) exceeds 64 bits for the largest values).
// Usage: ./clock_sweep N out.bin. Writes four doubles per n: tau(n), max under T, steps of T, odd steps.
#include <stdio.h>
#include <stdlib.h>
typedef unsigned __int128 u128;
static u128 cost(unsigned long long v){ u128 k=v/2; return v%2 ? 5*k*k+9*k+4 : 2*k*(k+1); }
int main(int argc,char**argv){
  long N=atol(argv[1]); FILE*f=fopen(argv[2],"wb");
  for(long n=1;n<=N;n++){
    unsigned long long x=n,mx=n; u128 a=0; long st=0,od=0;
    while(x!=1){ a+=cost(x); if(x%2){od++;x=(3*x+1)/2;} else x/=2; st++; if(x>mx)mx=x; }
    double r[4]={(double)a,(double)mx,(double)st,(double)od}; fwrite(r,8,4,f);
  }
  fclose(f); return 0;
}
