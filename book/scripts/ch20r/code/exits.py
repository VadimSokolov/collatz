import numpy as np, math, sys, cmath
def full_table(nmax,Bmax,pad):
    Bmin=-pad; W=Bmax-Bmin+1; u=[1+0j]*W; tab=[np.array(u)]
    for n in range(1,nmax+1):
        q=3**n; inv2=(q+1)//2; x=pow(inv2,pad,q); v=[0j]*W
        for i in range(W):
            fr=x/q; v[i]=cmath.exp(2j*math.pi*fr)*u[i]; x=(x<<1)%q
        c=0j; nl=[0j]*W
        for i in range(1,W):
            c=(c+v[i-1])*0.5; nl[i]=c
        u=nl; tab.append(np.array(u))
    return tab,Bmin
n=int(sys.argv[1]); k=int(sys.argv[2]); pad=int(2*n+14*math.sqrt(n)+60)
tab,Bmin=full_table(n,k,pad)
U=lambda L,B: tab[L][B-Bmin]
amp={k:1+0j}; contrib=np.zeros(n+1,dtype=complex)   # contrib[L]: exit while stepping from level L to L-1
for L in range(n,0,-1):
    q=3**L; inv2=(q+1)//2; new={}
    for C,am in amp.items():
        w=0.5
        for a in range(1,C+70):
            Cp=C-a
            if Cp>=0:
                ph=cmath.exp(2j*math.pi*(pow(2,Cp,q)/q)); new[Cp]=new.get(Cp,0)+am*w*ph
            else:
                if Cp<Bmin: break
                ph=cmath.exp(2j*math.pi*(pow(inv2,-Cp,q)/q)); contrib[L]+=am*w*ph*U(L-1,Cp)
            w*=0.5
            if w<1e-40: break
    amp=new
never=sum(amp.values())
tot=contrib.sum()+never
print('n',n,'k',k,'total',abs(tot),'direct',abs(U(n,k)),'never exit share',abs(never)/abs(tot))
cum=0j
print('exit at step with l=L-1 levels remaining: |partial sum up to l| / |total|')
for l in range(0,n):
    cum+=contrib[l+1]
    if l in (0,1,2,3,5,8,12,16,20,30,40,60,80,120,n-1): print(' l<=',l, round(abs(cum+never)/abs(tot),4), ' |contrib_l|/|tot|',round(abs(contrib[l+1])/abs(tot),4))
