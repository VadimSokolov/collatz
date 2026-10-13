"""Records of f(r) = r/(3 Lambda_r), Lambda_r = (ceil(r theta) - r theta) log 2,
computed in exact integer arithmetic with theta truncated to 10^-120.
Compares records with the one sided ladder q_{k-1} + j q_k (k even, 1<=j<=a_{k+1})."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
import mpmath as mp, sys, json
mp.mp.dps = 200
TH = mp.log(3)/mp.log(2); L2 = mp.log(2)
D = 10**120
T = int(mp.floor(TH*D))          # theta ~ T/D, error < 1e-120
def cf(x,n):
    a=[];y=mp.mpf(x)
    for _ in range(n):
        ai=int(mp.floor(y));a.append(ai);y=1/(y-ai)
    return a
A=cf(TH,60)
q=[0,1]; p=[1,A[0]]
for ai in A[1:45]:
    q.append(ai*q[-1]+q[-2]); p.append(ai*p[-1]+p[-2])
Q=lambda k:q[k+1]
# upper ladder: k even
upper=[]; both=set()
for k in range(0,40):
    for j in range(0,A[k+1]+1):
        v=Q(k-1)+j*Q(k)
        if v>=1: both.add(v)
        if k%2==0 and j>=1: upper.append(v)
upper=sorted(set(upper))
R = int(sys.argv[1]) if len(sys.argv)>1 else 3*10**6
best_r, best_d = None, None   # record of r/d with d = ceil(rT/D)*D - rT (integer)
recs=[]; dmin=None; dmins=[]
for r in range(1,R+1):
    x=r*T; d = (-x) % D       # = D*(ceil(r theta) - r theta), exact for truncated theta
    if dmin is None or d<dmin: dmin=d; dmins.append(r)
    if best_r is None or r*best_d > best_r*d:
        best_r,best_d=r,d; recs.append(r)
U=set(upper)
print("R =",R)
print("number of strict records of f:",len(recs))
print("records:",recs)
print("records not on upper ladder:",[r for r in recs if r not in U])
print("strict one sided minima of d:",dmins)
print("one sided minima equal upper ladder up to R:", dmins==[u for u in upper if u<=R])
print("upper rungs <= R that are not records of f:",[u for u in upper if u<=R and u not in set(recs)])
conv=set(q[1:])
print("records that are not convergent denominators:",[r for r in recs if r not in conv][:20])
json.dump({"recs":recs},open(BOOK + "/scripts/data/ch08_records.json","w"))
