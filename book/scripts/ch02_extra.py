import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Extra numbers and the Syracuse figure for chapter 2
import math
from collections import Counter, defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys as _s; _s.path.insert(0, BOOK + '/scripts'); from bookstyle import *
FIG=BOOK + "/figures/"
def T(n): return n//2 if n%2==0 else (3*n+1)//2
def v2(n):
    c=0
    while n%2==0: n//=2; c+=1
    return c
def parvec(n,k):
    v=[]
    for _ in range(k): v.append(n%2); n=T(n)
    return v
# (6) k=3 table
for r in range(8): print("k=3",r,parvec(r,3))
for r in range(16): print("k=4",r,parvec(r,4))
# (3) inverse formula n = -c_k 3^{-m} mod 2^k
def ck(v):
    c=0
    for i,b in enumerate(v):
        if b: c=3*c+2**i
    return c
for k in range(1,15):
    for r in range(2**k):
        v=parvec(r,k); m=sum(v)
        assert (-ck(v)*pow(3,-m,2**k))%2**k==r
print("inverse formula ok k<15")
# c_70(27)
v=parvec(27,70); m=sum(v); print("m",m,"c70(27)",ck(v),"2^70-27*3^41",2**70-27*3**41, ck(v)==2**70-27*3**41)
# Syracuse orbit of 27: the climb happens in the first 32 strides (ex:ch2:27)
y=[27];ek=[]
while y[-1]!=1:
    t=3*y[-1]+1; k=v2(t); ek.append(k); y.append(t>>k)
print("27: strides",len(ek),"sum k",sum(ek),"first 32 strides: sum k",sum(ek[:32]),"value",y[32],"largest odd value",max(y))
# (2) mean log(S/n) over odd n < 2^20 and the stats of k
N=2**20
ks=[];incs=[]
for n in range(1,N,2):
    t=3*n+1; k=v2(t); ks.append(k); incs.append(math.log((t>>k)/n))
ks=np.array(ks); incs=np.array(incs)
print("odd n<2^20: mean k",ks.mean(),"mean log S/n",incs.mean(),"var",incs.var(), "frac ascents",(incs>0).mean())
print("log(3/4)",math.log(3/4),"2 log2^2",2*math.log(2)**2,"(log3)^2/4",math.log(3)**2/4)
print("odd n<2^20: arithmetic mean of S(n)/n", np.exp(incs).mean())
# (4) P(tau>k)
th=math.log(2)/math.log(3); D=th*math.log(2*th)+(1-th)*math.log(2*(1-th))
st={0:1};P=[]
for j in range(1,201):
    new=defaultdict(int)
    for mm,c in st.items():
        for b in (0,1):
            m2=mm+b
            if 3**m2>=2**j: new[m2]+=c
    st=new; P.append(sum(st.values())/2**j)
for k in [20,40,60,80,100,150,200]:
    print("k",k,"P(tau>k)",P[k-1],"exp(-Dk)",math.exp(-D*k),"ratio",P[k-1]/math.exp(-D*k))
# (5) tau values
vals=sorted(set(math.floor(m*math.log2(3))+1 for m in range(0,40)))
nz=[j for j in range(1,61) if (P[j-2] if j>1 else 1.0)-P[j-1]>0]
print("values floor(m lt)+1 <=60:",[x for x in vals if x<=60]); print("k with new stopped residues <=60:",nz)
# The Syracuse figure ch2_syracuse.pdf is drawn by ch02_figures.py in the shared book style.
