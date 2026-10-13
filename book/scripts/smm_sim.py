import numpy as np
rng=np.random.default_rng(1); beta=np.log(1.5)/np.log(2)
a=0.02; T=2_000_000; e=0.0; ep=0.0; acc=[]
xi=rng.standard_normal(T)
for t in range(T):
    if t%2==1:
        y=e+beta*(e-ep); new=(1-a)*y-xi[t]
    else: new=(1-a)*e-xi[t]
    ep,e=e,new
    if t>1000: acc.append(e*e)
v=np.mean(acc); 
th=( (2+2*beta+beta**2)/(2*(2+beta)) )/a
print("sim",v,"lead theory",th, "ratio",v/th)
