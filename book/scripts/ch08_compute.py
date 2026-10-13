"""Chapter 8: continued fraction of log_2 3 at 300 digits (first 60 partial
quotients confirmed at 150) and the convergents with their sides."""
import mpmath as mp
mp.mp.dps = 300
TH = mp.log(3)/mp.log(2)
L2 = mp.log(2)

def cf(x, n):
    a=[]; y=mp.mpf(x)
    for _ in range(n):
        ai=int(mp.floor(y)); a.append(ai); y=1/(y-ai)
    return a
A = cf(TH, 60)
print("CF 60 terms:", A)
# check stability at lower precision
mp.mp.dps=150; A2=cf(mp.log(3)/mp.log(2),60); mp.mp.dps=300
print("agree with 150 digits for first", next((i for i in range(60) if A[i]!=A2[i]),60), "terms")
p=[1,A[0]]; q=[0,1]   # p_{-1},p_0
for ai in A[1:40]:
    p.append(ai*p[-1]+p[-2]); q.append(ai*q[-1]+q[-2])
# index k -> p[k+1]
print("k a_k p_k q_k side delta_k")
for k in range(0,30):
    pk,qk=p[k+1],q[k+1]; d=qk*TH-pk
    print(k, A[k], pk, qk, "below" if d>0 else "above", mp.nstr(abs(d),6))
