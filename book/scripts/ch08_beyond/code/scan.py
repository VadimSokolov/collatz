import numpy as np, math
from scipy.optimize import brentq
from local53b import run, Vb
L=1.0; best=(0,)
for mu in [0.02,0.05]:
    for c in [1.5,4.0]:
        g=lambda x,mu=mu,c=c: mu*x+(L-mu)*np.tanh(x-c)
        xs=brentq(lambda x: mu*x+(L-mu)*math.tanh(x-c),-50,50)
        for s in [1.0,2.0,4.0,8.0]:
            r=[]
            for seed in (0,1):
                mg,ms=run(g,xs,1/L,s,nch=3000,T=30000,burn=12000,seed=seed); r.append(ms/mg)
            print(f'mu {mu} c {c} sigma {s}: ratio {r[0]:.4f} {r[1]:.4f}'); best=max(best,(max(r),mu,c,s))
print('max',best,'V(beta)',Vb,'5/3',5/3)
