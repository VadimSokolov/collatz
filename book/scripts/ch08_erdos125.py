"""Erdos 125: exact density of A+B in [0,N], A = base 3 digits in {0,1}, B = base 4 digits in {0,1}.
Usage: python3 -I ch08_erdos125.py EXP   (computes up to N = 4^EXP).  Saves a thinned curve for the figure."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import numpy as np, sys
E=int(sys.argv[1]) if len(sys.argv)>1 else 11
N=4**E
def digitset(base,limit):
    s=np.array([0],dtype=np.int64); p=1
    while p<=limit:
        s=np.concatenate([s,s+p]); p*=base
    return np.sort(s[s<=limit])
A=digitset(3,N); B=digitset(4,N)
mark=np.zeros(N+1,dtype=bool)
for b in B:
    a=A[A<=N-b]; mark[a+b]=True
cnt=np.cumsum(mark,dtype=np.int64)
dens=lambda n: cnt[n]/(n+1)
print("N = 4^%d = %d, |A| = %d, |B| = %d"%(E,N,len(A),len(B)))
for k in range(7,E+1): print("4^%d  %.4f"%(k,dens(4**k)))
m=10
while 3**m<=N: print("3^%d  %.4f"%(m,dens(3**m))); m+=1
lo=10**4
n=np.arange(lo,N+1); d=cnt[lo:]/(n+1.0); i=int(np.argmin(d))
print("min density over [1e4, N]: %.6f at N = %d"%(d[i],lo+i))
# thinned curve on a log grid, keeping local minima just below powers of 3
grid=np.unique(np.round(np.logspace(3,np.log10(N),4000)).astype(np.int64))
extra=[3**m-1 for m in range(7,40) if 3**m-1<=N]
grid=np.unique(np.concatenate([grid,extra]))
np.save(BOOK + "/scripts/ch08_erdos125_curve.npy",np.vstack([grid,cnt[grid]/(grid+1.0)]))
