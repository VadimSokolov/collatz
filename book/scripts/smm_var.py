import numpy as np, sympy as sp
beta=np.log(1.5)/np.log(2)
def mats(a,b):
    D=np.array([[1-a,0],[1,0]]); N=np.array([[(1-a)*(1+b),-(1-a)*b],[1,0]])
    H=np.array([[1-a+b,-b],[1,0]]); return D,N,H
def stat(As,Q):
    # periodic stationary covariance, iterate
    P=np.zeros((2,2))
    for _ in range(200000):
        P0=P.copy()
        for A in As: P=A@P@A.T+Q
        if np.abs(P-P0).max()<1e-15: break
    # average var over phases
    vs=[];Pc=P.copy()
    for A in As:
        vs.append(Pc[0,0]); Pc=A@Pc@A.T+Q
    return np.mean(vs)
Q=np.array([[1.,0],[0,0]])  # alpha^2 sigma^2 =1
print("a, GD, NAG585, HB585, SMM585, NAG9, SMM9 ; ratios to GD; rates per step")
for a in [0.001,0.01,0.04,0.1,0.25,0.5,1.0]:
    D,N,H=mats(a,beta); D9,N9,H9=mats(a,0.9)
    gd=stat([D],Q); nag=stat([N],Q); hb=stat([H],Q); smm=stat([D,N],Q); nag9=stat([N9],Q); smm9=stat([D9,N9],Q)
    r=lambda M: max(abs(np.linalg.eigvals(M)))
    print(f"{a:6} {gd:10.4g} nag {nag/gd:6.3f} hb {hb/gd:6.3f} smm {smm/gd:6.3f} nag9 {nag9/gd:6.3f} smm9 {smm9/gd:6.3f} | rate gd {1-a:.4f} nag {r(N):.4f} smm {np.sqrt(r(N@D)):.4f}")
# symbolic small-a limit for SMM
a,b=sp.symbols('a b',positive=True)
D=sp.Matrix([[1-a,0],[1,0]]); N=sp.Matrix([[(1-a)*(1+b),-(1-a)*b],[1,0]])
p11,p12,p22=sp.symbols('p11 p12 p22'); P=sp.Matrix([[p11,p12],[p12,p22]]); Qs=sp.Matrix([[1,0],[0,0]])
P1=D*P*D.T+Qs; P2=N*P1*N.T+Qs
sol=sp.solve([sp.Eq(P2[0,0],p11),sp.Eq(P2[0,1],p12),sp.Eq(P2[1,1],p22)],[p11,p12,p22],dict=True)[0]
v_even=sp.simplify(sol[p11]); v_odd=sp.simplify(P1[0,0].subs(sol))
avg=sp.simplify((v_even+v_odd)/2)
print("SMM avg var exact:", sp.factor(avg))
lead=sp.limit(avg*a,a,0); print("SMM leading a*Var as a->0:", sp.simplify(lead))
# NAG every step
Pn=N*P*N.T+Qs
soln=sp.solve([sp.Eq(Pn[0,0],p11),sp.Eq(Pn[0,1],p12),sp.Eq(Pn[1,1],p22)],[p11,p12,p22],dict=True)[0]
vn=sp.factor(soln[p11]); print("NAG var exact:", vn); print("NAG lead:", sp.simplify(sp.limit(vn*a,a,0)))
print("SMM lead at beta:", float(lead.subs(b,beta)), "NAG lead:", float(sp.limit(vn*a,a,0).subs(b,beta)), "GD lead 0.5")
