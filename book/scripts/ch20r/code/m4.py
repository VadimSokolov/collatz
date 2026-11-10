import numpy as np, math, sys
from scipy.optimize import minimize_scalar
th=math.log2(3)
def setup(r):
    t=np.arange(r+1,dtype=np.float64)
    fl=np.floor(t*th); x=t*th-fl
    m=2.0**(-x[:r]); g=(fl[1:]-fl[:-1]).astype(int)   # g_t for t=0..r-1
    return m,g
def logZ(r,lams,J=64):
    """log of sum over admissible words (s=ceil(r th)) of exp(lam*Phi), vectorised over lams"""
    m,g=setup(r); L=np.asarray(lams)[:,None]
    v=np.zeros((len(lams),J+1)); v[:,0]=1.0; logs=np.zeros(len(lams))
    pw=2.0**(-np.arange(J+1))[None,:]
    for t in range(r):
        w=v*np.exp(L*m[t]*pw)           # weight for index t at depth j
        if t==r-1:
            tot=w.sum(1); return logs+np.log(tot)   # last step forced
        # j' <= j+g-1 : new[j']=sum_{j>=j'-g+1} w[j]
        suf=np.cumsum(w[:,::-1],1)[:,::-1]
        if g[t]==1: v=suf
        else:
            v=np.empty_like(w); v[:,0]=suf[:,0]; v[:,1:]=suf[:,:-1]; v[:,-1]+=suf[:,-1]  # clamp depth at J: upper bound
        s=v.max(1); v/=s[:,None]; logs+=np.log(s)
def Br(r): return setup(r)[0].sum()
if __name__=="__main__":
    cs=[0.5,0.6,0.7,0.75,0.751,0.8,0.85,0.9,0.95,0.99]
    lam=np.concatenate([[0],np.geomspace(0.01,60,260)])
    for r in [306,2301,190537]:
        B=Br(r); lz=logZ(r,lam)
        row=[]
        for c in cs:
            vals=(lz-lam*c*B)/(r*math.log(10)); i=vals.argmin(); row.append((round(vals[i],4),round(float(lam[i]),2)))
        print(r,'B/r',round(B/r,5),'total',round(lz[0]/(r*math.log(10)),4)); print('   ',row); sys.stdout.flush()
