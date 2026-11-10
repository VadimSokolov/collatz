import numpy as np, math, sys
from scipy.optimize import brentq
beta=math.log(1.5)/math.log(2); L=1.0; alpha=1/L
def stat(mu,c,s,W=None,dx=0.1):
    """stationary mean square errors of GD and SMM (phase averaged) by quadrature on a grid; 1-d f=mu x^2/2+(L-mu)log cosh(x-c)"""
    g=lambda x: mu*x+(L-mu)*np.tanh(x-c)
    xs=brentq(lambda x: mu*x+(L-mu)*math.tanh(x-c),-1e4,1e4)
    if W is None: W=14*s*math.sqrt(alpha/mu)+40
    x=np.arange(xs-W,xs+W+dx,dx); n=len(x)
    def kernel(mapped,sd):
        # K[i,j]=P(next in cell i | current x_j): next = mapped_j + sd*z
        K=np.exp(-0.5*((x[:,None]-mapped[None,:])/sd)**2); K/=K.sum(0,keepdims=True); return K
    G1=x-alpha*g(x); G2=x-(1+beta)*alpha*g(x)
    Kgd=kernel(G1,alpha*s)
    def stationary(K):
        p=np.ones(n)/n
        for _ in range(200000):
            q=K@p; q/=q.sum()
            if np.abs(q-p).sum()<1e-13: break
            p=q
        return q
    pg=stationary(Kgd); mgd=(pg*(x-xs)**2).sum()
    K2=kernel(G2,(1+beta)*alpha*s)       # X -> u
    Ke=Kgd@K2                            # u -> X' uses G1 and sd alpha*s : same kernel as GD step
    pe=stationary(Ke); me=(pe*(x-xs)**2).sum()
    po=Kgd@pe; po/=po.sum(); mo=(po*(x-xs)**2).sum()
    return mgd,0.5*(me+mo),xs,n,(pg*(x-xs)).sum(),(pe*(x-xs)).sum()
if __name__=="__main__":
    for mu,s,dx in [(0.05,2.0,0.1),(0.01,3.0,0.15),(0.005,3.0,0.2),(0.005,3.0,0.3)]:
        mg,ms,xs,n,bg,be=stat(mu,1.5,s,dx=dx)
        print(f'mu {mu} sigma {s} grid {n} dx {dx}: GD mse {mg:.4f} SMM mse {ms:.4f} ratio {ms/mg:.5f}  bias GD {bg:.3f} SMM even {be:.3f}',flush=True)
