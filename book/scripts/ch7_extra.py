import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Extra checks for Chapter 7: closed forms, bounds, slow mode cyclic step formula, logistic predictions.
import numpy as np, sympy as sp
a,b=sp.symbols('a b',positive=True)
rho=(1-a)*(1-(1+b)*a)
ve=((1-a)**2*(1+b)**2+1)/(1-rho**2)
vo=(1-a)**2*ve+1
avg=(ve+vo)/2
Q=(2+2*b+b**2)-(1+b)*(2+3*b)*a+(1+b)*(1+3*b)*a**2-b*(1+b)*a**3
den=a*(2+b-(1+b)*a)*(2-(2+b)*a+(1+b)*a**2)
print("avg == Q/den:", sp.simplify(avg-Q/den)==0)
print("1-rho^2 factor:", sp.factor(1-rho**2))
print("ve*a limit", sp.simplify(sp.limit(a*ve,a,0)), " vo*a limit", sp.simplify(sp.limit(a*vo,a,0)))
print("V = b + 2/(2+b):", sp.simplify((2+2*b+b**2)/(2+b)-(b+2/(2+b)))==0)
print("2(2+2b+b^2)=(2+b)^2+b^2:", sp.expand(2*(2+2*b+b**2)-(2+b)**2-b**2)==0)
# NAG and HB closed forms
def ar2(a1,a2): return (1-a2)/((1+a2)*((1-a2)**2-a1**2))
vh=sp.factor(ar2(1+b-a,-b)); vn=sp.factor(ar2((1-a)*(1+b),-(1-a)*b))
print("HB var:", vh); print("HB target check:", sp.simplify(vh-(1+b)/((1-b)*a*(2*(1+b)-a)))==0)
print("NAG var:", vn)
print("NAG slow lim a*v:", sp.limit(a*vn,a,0))
gd=1/(a*(2-a))
# check ratio SMM/GD <= V(b) on grid, and monotone in a
V=lambda bb:(2+2*bb+bb**2)/(2+bb)
f=sp.lambdify((a,b),avg/gd,'numpy')
A=np.linspace(1e-4,1,2000); worst=-1
for bb in np.linspace(0,0.999,400):
    r=f(A,bb); worst=max(worst,np.max(r-V(bb))); 
print("max over grid of ratio - V(beta):", worst)
diff=sp.factor(sp.together(V(b)*gd-avg)); print("V*gd - avg factored:", diff)
# heuristic
for bb in [np.log(1.5)/np.log(2),0.9]: print("heuristic half 1/(1-b):",0.5/(1-bb)," exact V:",V(bb))
# exact at a=0.02, beta0
b0=np.log(1.5)/np.log(2)
print("exact avg a=0.02:", float(avg.subs({a:0.02,b:b0})), "lead:", V(b0)/(2*0.02))
# logistic predictions at a = 1/67.4: iteration ratio GD/method = log r_GD / log r_method
aa=1/67.4
def rate_nag(aa,bb):
    return max(abs(np.roots([1,-(1-aa)*(1+bb),(1-aa)*bb])))
def rate_hb(aa,bb):
    return max(abs(np.roots([1,-(1+bb-aa),bb])))
r_smm=np.sqrt(abs((1-aa)*(1-(1+b0)*aa)))
lg=np.log(1-aa)
print("pred speedup SMM", lg/np.log(r_smm))
for bb in [0.5,b0,2/3,0.9]:
    print(f"b={bb:.3f} NAG {lg/np.log(rate_nag(aa,bb)):.3f} HB {lg/np.log(rate_hb(aa,bb)):.3f}")
# HB at alpha=1/L, but logistic table: HB fixed beta uses alpha = 1/L presumably
# proof check for uniform bound Var_SMM <= V(beta) Var_GD on 0<a<=1
P=(b+2)*a**3+(b**2-2*b-7)*a**2+(-2*b**2+b+10)*a+b**2-6
num=sp.numer(diff); print("numerator matches -b(b+1)P:", sp.expand(num + b*(b+1)*P)==0 or sp.expand(num - b*(b+1)*P)==0, sp.expand(num + b*(b+1)*P)==0)
qd=-(5+b-b**2)*a**2+(10+b-2*b**2)*a-(6-b**2)
print("disc:", sp.expand((10+b-2*b**2)**2-4*(5+b-b**2)*(6-b**2)))
print("P - qd = (b+2)a^2(a-1):", sp.expand(P-qd-(b+2)*a**2*(a-1))==0)
print("den of diff:", sp.denom(diff))
print("speedups over GD (iterations GD / iterations method), slow mode a=1/67.4:")
print("SMM %.3f"%(np.log(r_smm)/lg))
for bb in [0.5,b0,2/3,0.9]:
    print("b=%.3f NAG %.3f HB %.3f  underdamped NAG? %s"%(bb,np.log(rate_nag(aa,bb))/lg,np.log(rate_hb(aa,bb))/lg,(1-aa)*(1+bb)**2<4*bb))
# lower bounds on every eigendirection: V_SMM - V_GD (prop:ch7:uniform) and V_NAG - V_GD (reconciliation)
lb=b*(1+b)*(1-a)**2*((1-a)**2+1)/(a*(2-a)*(2+b-(1+b)*a)*(1+rho))
print("V_SMM - V_GD identity:", sp.simplify(avg-gd-lb)==0)
nagd=2*b*(1+b)*(1-a)**2/(a*(2-a)*(1-b+a*b)*(2*(1+b)-a*(1+2*b)))
print("V_NAG - V_GD identity:", sp.simplify(vn-gd-nagd)==0)
# observed speedups from the logistic experiment of ch7_experiments.py
import json
lgc=json.load(open(BOOK + '/scripts/ch7_experiments.json'))['logistic']['counts']
for k in ["HB","NAG"]: print(k,"observed",[round(lgc["GD"][0]/lgc[k+" "+nm][0],3) for nm in ["0.5","0.585","0.667","0.9"]])
print("SMM observed",round(lgc["GD"][0]/lgc["SMM"][0],3))
