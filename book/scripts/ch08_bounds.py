"""Block bounds B_k (standard indexing q_0 = 1), sharp bound over the upper ladder,
ties at 150 digits, negative cycles of T_0, trivial cycle check, Eliahou check."""
import mpmath as mp
mp.mp.dps = 150
TH = mp.log(3)/mp.log(2); L2=mp.log(2); L3=mp.log(3)
V = mp.mpf(2)**68
def cf(x,n):
    a=[];y=mp.mpf(x)
    for _ in range(n):
        ai=int(mp.floor(y));a.append(ai);y=1/(y-ai)
    return a
A=cf(TH,50)
q=[0,1];p=[1,A[0]]
for ai in A[1:40]:
    q.append(ai*q[-1]+q[-2]); p.append(ai*p[-1]+p[-2])
Q=lambda k:q[k+1]
def nn(r):
    f=r*TH-mp.floor(r*TH); return min(f,1-f)
def lam(r):
    s=int(mp.ceil(r*TH)); return s, s*L2-r*L3
def f(r):
    s,l=lam(r); return s, l, r/(3*l)
print("Block bounds, standard indexing:")
for k in range(0,25):
    B=Q(k+1)/(3*nn(Q(k))*L2)
    print(k, Q(k), Q(k+1), mp.nstr(B,4), "<V" if B<V else ">V")
print("2^68 =", mp.nstr(V,8))
print("f at convergent denominators:")
for k in range(0,25):
    s,l,M=f(Q(k)); print(k,Q(k),s,mp.nstr(l,5),mp.nstr(M,5), "d=",mp.nstr(mp.ceil(Q(k)*TH)-Q(k)*TH,6))
# upper ladder
upper=[]
for k in range(0,36,2):
    for j in range(1,A[k+1]+1):
        upper.append(Q(k-1)+j*Q(k))
upper=sorted(set(upper))
lim=2*10**11
up=[u for u in upper if u<=lim]
print("upper rungs <= 2e11:",len(up))
both=set()
for k in range(0,36):
    for j in range(0,A[k+1]+1):
        v=Q(k-1)+j*Q(k)
        if 1<=v<=lim: both.add(v)
print("all rungs (both parities) <= 2e11:",len(both))
run=mp.mpf(0); prev=None
for u in up:
    s,l,M=f(u)
    if M>V:
        print("first upper rung exceeding 2^68: r =",u,"s =",s,"bound",mp.nstr(M,6)); print("running max before:",mp.nstr(run,6),"at",prev); break
    run=M; prev=u
print("record values along ladder up to failure (r, s, f):")
for u in up:
    s,l,M=f(u)
    print(u,s,mp.nstr(M,5))
    if M>V: break
# ties
for (m,qq) in [(3,15601),(3,190537),(5,190537)]:
    _,l1=lam(qq); _,lm=lam(m*qq); print("tie",m,qq,mp.nstr(lm-m*l1,5))
# double precision artefact
import math
th=math.log(3)/math.log(2)
def fd(r):
    return r/(3*(math.ceil(r*th)-r*th)*math.log(2))
for (m,qq) in [(3,15601),(3,190537),(5,190537)]:
    print("double:",m*qq, fd(m*qq)-fd(qq))
# Eliahou check: verification 2^40
V40=mp.mpf(2)**40
run=0
for u in up:
    s,l,M=f(u)
    if M>V40: print("V=2^40 first failing upper rung r=",u,"s=",s); break
for k in range(0,25):
    B=Q(k+1)/(3*nn(Q(k))*L2)
    if B>V40: print("V=2^40 block: r >= q_k =",Q(k), "s there =", int(mp.ceil(Q(k)*TH))); break
