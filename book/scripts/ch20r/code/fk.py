import cmath, math, sys, numpy as np, pickle
def run(nmax, Bmax=None, pad=None):
    """u_n(B)=phi_n(2^B)=E exp(2 pi i 2^B Syrac_n/3^n), exact Feynman-Kac recursion
       u_n(B)=sum_{a>=1} 2^{-a} e(2^{B-a}/3^n) u_{n-1}(B-a), u_0=1."""
    if Bmax is None: Bmax=int(1.75*nmax)+5
    if pad is None: pad=int(2*nmax+14*math.sqrt(nmax)+60)
    Bmin=-pad; W=Bmax-Bmin+1
    u=np.ones(W,dtype=complex); out={}
    twopi=2*math.pi
    for n in range(1,nmax+1):
        q=3**n; inv2=(q+1)//2
        # v(C)=e(2^{C}/3^n) u_{n-1}(C), C index
        ph=np.empty(W,dtype=complex)
        # 2^C mod q for C=Bmin..Bmax : start at 2^Bmin = inv2^pad
        x=pow(inv2,pad,q)
        fr=[0.0]*W
        for i in range(W):
            fr[i]=x/q if n<=300 else float(x*(1<<60)//q)/(1<<60)
            x=(x<<1)%q
        ph=np.exp(1j*twopi*np.array(fr))
        v=ph*u
        # u_n(B)=sum_{a>=1}2^{-a} v(B-a): c(B)=(c(B-1)+v(B-1))/2, c(Bmin)=0 (truncation)
        new=np.zeros(W,dtype=complex); c=0j
        vl=v.tolist(); nl=[0j]*W
        for i in range(1,W):
            c=(c+vl[i-1])*0.5; nl[i]=c
        u=np.array(nl)
        a=np.abs(u[-Bmin:])   # B>=0
        k=int(a.argmax()); out[n]=(float(a[k]),k,a[:].copy() if n in keep else None)
    return out
keep=set()
if __name__=="__main__":
    nmax=int(sys.argv[1]); keep=set(int(x) for x in sys.argv[2:])
    out=run(nmax)
    pickle.dump({n:(v[0],v[1]) for n,v in out.items()},open(f'fk_{nmax}.pkl','wb'))
    prev=None
    for n in sorted(out):
        m,k,_=out[n]
        if n<=20 or n%10==0: print(n,'max',f'{m:.6g}','k',k,'k/n',round(k/n,3),'d^2max',round(m*n*n,3),'ratio',round(m/prev,4) if prev else '')
        prev=m
