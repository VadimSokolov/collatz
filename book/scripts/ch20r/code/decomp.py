import numpy as np, math, sys, pickle
def run(nmax, mode, Bmax=None):
    """mode 'full': all paths; 'real': only paths whose exponent stays >=0 (all phases are real numbers 2^B/3^l);
       'barrier': exponent >=0 and phase 2^B/3^l < 1/2 at every step."""
    if Bmax is None: Bmax=int(1.75*nmax)+5
    pad=int(2*nmax+14*math.sqrt(nmax)+60) if mode=='full' else 0
    Bmin=-pad; W=Bmax-Bmin+1
    u=[1+0j]*W; out={}
    for n in range(1,nmax+1):
        q=3**n; inv2=(q+1)//2; x=pow(inv2,pad,q); v=[0j]*W
        half=q//2
        for i in range(W):
            ok=True
            if mode=='barrier' and 2*x>=q: ok=False     # here C>=0 so x=2^C mod q; need 2^C<q/2 as integers
            if mode=='barrier' and (i> q.bit_length()+2): ok=False   # 2^C >= q certainly
            if ok:
                fr=x/q if n<=300 else float((x<<60)//q)/(1<<60)
                v[i]=complex(math.cos(2*math.pi*fr),math.sin(2*math.pi*fr))*u[i]
            x=(x<<1)%q
        c=0j; nl=[0j]*W
        for i in range(1,W):
            c=(c+v[i-1])*0.5; nl[i]=c
        u=nl
        out[n]=np.array(u[-Bmin:])
    return out
if __name__=="__main__":
    nmax=int(sys.argv[1])
    full=run(nmax,'full'); real=run(nmax,'real'); bar=run(nmax,'barrier')
    pickle.dump({n:(full[n],real[n],bar[n]) for n in full if n%10==0 or n<=20},open('decomp.pkl','wb'))
    for n in [10,18,30,50,80,120,160,200,250,300]:
        if n>nmax: break
        a=np.abs(full[n]); k=int(a.argmax())
        print(n,'k',k,'|u|',f'{a[k]:.4e}','|real|/|u|',round(abs(real[n][k])/a[k],4),'|barrier|/|u|',round(abs(bar[n][k])/a[k],4),'|u-real|/|u|',round(abs(full[n][k]-real[n][k])/a[k],4),'|u-bar|/|u|',round(abs(full[n][k]-bar[n][k])/a[k],4))
