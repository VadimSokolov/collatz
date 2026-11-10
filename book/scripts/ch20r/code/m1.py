from mpmath import mp, mpf, log, floor, ceil
import math
mp.prec=200
th=log(3)/log(2)
N=200000
# delta(t), with high precision frac then to float
import numpy as np
d=np.empty(N+1); 
for t in range(1,N+1):
    x=t*th; d[t]=float(ceil(x)-x)
d[0]=1.0  # dummy
def B(r,dl):
    dd=d[1:r]
    return 1.0+np.where(dd>dl,2.0**(dd-1),1.0).sum()
R=np.zeros(N+1)
# O(N^2) too slow for 2e5 -> restrict
N2=20000
rungs=[];m=2
for r in range(1,N2+1):
    if d[r]<m: m=d[r]; rungs.append(r)
print('rungs',rungs)
best=-1;rec=[]
for r in range(1,N2+1):
    R[r]=B(r,d[r])/(3*math.log(2)*d[r])
    if R[r]>best: best=R[r];rec.append(r)
print('records',rec)
print('records subset of rungs',set(rec)<=set(rungs), 'rungs not records',[x for x in rungs if x not in rec])
# closest calls: non-rungs with R(r)/max_{r'<r} R(r') largest
run=0;worst=[]
for r in range(1,N2+1):
    if r not in set(rungs) and run>0: worst.append((R[r]/run,r))
    run=max(run,R[r])
worst.sort(reverse=True);print(worst[:8])
# test candidate inequalities for non-rung r: R(r)<=max(R(ri),R(u))
import bisect
bad=0;bad2=0
for r in range(2,N2+1):
    i=bisect.bisect_right(rungs,r)-1; ri=rungs[i]
    if ri==r: continue
    u=r-ri
    if R[r]>max(R[ri],R[u])*(1+1e-12): bad+=1
    if R[r]>R[ri]*(1+1e-12): bad2+=1
print('mediant fails',bad,'R(r)>R(ri) count',bad2)
