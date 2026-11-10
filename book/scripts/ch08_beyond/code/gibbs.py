import numpy as np, math
from scipy.optimize import brentq
beta=math.log(1.5)/math.log(2); V=beta+2/(2+beta); L=1.0
def m(T,F,xs,lo,hi,n=400001):
    x=np.linspace(lo,hi,n); w=np.exp(-(F(x)-F(xs))/T); return (w*(x-xs)**2).sum()/w.sum()
def logcosh(z): return np.abs(z)+np.log1p(np.exp(-2*np.abs(z)))-math.log(2)
print('V',V,'V^2',V*V)
for mu in (0.05,0.005,0.0005):
    c=1.5; F=lambda x: mu*x*x/2+(L-mu)*logcosh(x-c)
    xs=brentq(lambda x: mu*x+(L-mu)*math.tanh(x-c),-1e5,1e5)
    best=(0,0); row=[]
    for T in [0.00125,0.045,0.5,2.0,4.5,12.5,50,200,1000]:
        W=60*math.sqrt(V*T/mu)+50
        r=m(V*T,F,xs,xs-W,xs+W)/m(T,F,xs,xs-W,xs+W); row.append((T,round(r,4))); best=max(best,(r,T))
    print('mu',mu,row); print('   max ratio',best)
# asymmetric: curvature L on the left, kink then mu on the right
for mu in (0.005,0.0005):
    F=lambda x: np.where(x<0, L*x*x/2, mu*x*x/2+(L-mu)*(logcosh(x)+0*x))  # right side: |x|-like then mu
    best=(0,0)
    for T in np.geomspace(0.01,2000,60):
        W=60*math.sqrt(V*T/mu)+50
        r=m(V*T,F,0.0,-W,W)/m(T,F,0.0,-W,W); best=max(best,(r,T))
    print('asymmetric mu',mu,'max ratio',best)
