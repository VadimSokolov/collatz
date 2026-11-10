import numpy as np, math
from scipy.optimize import brentq
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from local53b import run, Vb
L=1.0
for mu in [0.01,0.005,0.002]:
    c=1.5
    g=lambda x,mu=mu,c=c: mu*x+(L-mu)*np.tanh(x-c)
    xs=brentq(lambda x: mu*x+(L-mu)*math.tanh(x-c),-500,500)
    for s in [1.5,2.0,3.0]:
        T=int(40/mu)+20000; r=[]
        for seed in (0,1,2):
            mg,ms=run(g,xs,1/L,s,nch=3000,T=T,burn=T//3,seed=seed); r.append(ms/mg)
        print(f'mu {mu} sigma {s}: ratios '+' '.join(f'{x:.4f}' for x in r), flush=True)
print('V(beta)',Vb,'5/3',5/3)
