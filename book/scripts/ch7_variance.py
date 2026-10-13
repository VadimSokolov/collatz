import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
# Stationary variance of SMM, NAG and HB under additive gradient noise; amplification figure,
# simulation check, candidate table, and predicted noise floor ratios for the kappa = 25 experiment.
import numpy as np, sympy as sp
beta0=np.log(1.5)/np.log(2)
# exact formulas in units of alpha^2 sigma^2
def v_gd(a): return 1/(a*(2-a))
def v_ar2(a1,a2): return (1-a2)/((1+a2)*((1-a2)**2-a1**2))
def v_hb(a,b): return v_ar2(1+b-a,-b)
def v_nag(a,b): return v_ar2((1-a)*(1+b),-(1-a)*b)
def v_smm(a,b):
    rho=(1-a)*(1-(1+b)*a); ve=((1-a)**2*(1+b)**2+1)/(1-rho**2); vo=(1-a)**2*ve+1
    return ve,vo,(ve+vo)/2
# symbolic check of the closed form against the polynomial given in VERIFIED
a,b=sp.symbols('a b',positive=True)
rho=(1-a)*(1-(1+b)*a); ve=((1-a)**2*(1+b)**2+1)/(1-rho**2); avg=((1+(1-a)**2)*ve+1)/2
poly=(a**3*b**2+a**3*b-3*a**2*b**2-4*a**2*b-a**2+3*a*b**2+5*a*b+2*a-b**2-2*b-2)/(a*(a*b+a-b-2)*(a**2*b+a**2-a*b-2*a+2))
print("closed form equals polynomial form:", sp.simplify(avg-poly)==0)
print("slow limit:", sp.simplify(sp.limit(a*avg,a,0)))
# AR(2) Yule Walker check symbolic
a1,a2,s=sp.symbols('a1 a2 s')
g0,g1=sp.symbols('g0 g1'); sol=sp.solve([sp.Eq(g0,a1*g1+a2*(a1*g0+a2*g1)+s**2*0+ (a1*g1+a2*g0)*0), ],[g0])
# direct: g1=a1 g0/(1-a2); g0 = a1 g1 + a2 g2 + s^2, g2 = a1 g1 + a2 g0
G0,G1,G2=sp.symbols('G0 G1 G2')
S=sp.solve([sp.Eq(G1,a1*G0+a2*G1),sp.Eq(G2,a1*G1+a2*G0),sp.Eq(G0,a1*G1+a2*G2+s**2)],[G0,G1,G2],dict=True)[0]
print("AR2 var:", sp.factor(S[G0]), " check:", sp.simplify(S[G0]-(1-a2)*s**2/((1+a2)*((1-a2)**2-a1**2)))==0)
# table of ratios
for bb in [beta0,0.9]:
    for aa in [0.001,0.01,0.04,0.1,0.25,0.5,1.0]:
        print(f"b={bb:.3f} a={aa}: NAG {v_nag(aa,bb)/v_gd(aa):.3f} HB {v_hb(aa,bb)/v_gd(aa):.3f} SMM {v_smm(aa,bb)[2]/v_gd(aa):.3f}" if aa<1 or True else "")
# candidate table: V, m, V/m
print("beta, V_smm, m_smm, V/m, V_every")
for bb in [0.0,0.5,beta0,2/3,0.9,0.99]:
    V=(2+2*bb+bb**2)/(2+bb); m=(2+bb)/2
    print(f"{bb:.3f} {V:.4f} {m:.4f} {V/m:.4f} {1/(1-bb):.3f}")
# simulation check at a=0.01, many chains
rng=np.random.default_rng(2026)
def sim(method,b,a,R=4000,T=6000,burn=2000):
    e=np.zeros(R); ep=np.zeros(R); acc=0.0; n=0
    for t in range(T):
        xi=rng.standard_normal(R)
        if method=="gd": new=(1-a)*e-xi
        elif method=="hb": new=(1-a+b)*e-b*ep-xi
        elif method=="nag": y=e+b*(e-ep); new=(1-a)*y-xi
        elif method=="smm":
            if t%2==1: y=e+b*(e-ep); new=(1-a)*y-xi
            else: new=(1-a)*e-xi
        ep,e=e,new
        if t>=burn: acc+=np.mean(e*e); n+=1
    return acc/n
aa=0.01; gdv=v_gd(aa)
sims={}
for bb in [0.25,0.5,beta0,0.75,0.9]:
    rs=sim("smm",bb,aa)/gdv; rn=sim("nag",bb,aa,T=8000,burn=3000)/gdv
    sims[bb]=(rs,rn)
    print(f"sim a=0.01 b={bb:.3f}: SMM {rs:.3f} (exact {v_smm(aa,bb)[2]/gdv:.3f}, slow {(2+2*bb+bb**2)/(2+bb):.3f}); NAG {rn:.3f} (exact {v_nag(aa,bb)/gdv:.3f}, slow {1/(1-bb):.3f})")
# figure (shared book style)
import sys; sys.path.insert(0, BOOK + '/scripts')
from bookstyle import *
B=np.linspace(0,0.97,300)
fig,ax=plt.subplots(figsize=(WIDTH,2.9))
ax.fill_between(B,1,(2+2*B+B**2)/(2+B),color=PAL['coral'],alpha=0.10,lw=0)
ax.plot(B,1/(1-B),color=PAL['violet'],lw=1.8,ls="--",label=r"every step momentum, $1/(1-\beta)$")
ax.plot(B,(2+2*B+B**2)/(2+B),color=PAL['coral'],lw=2.2,label=r"SMM, $(2+2\beta+\beta^2)/(2+\beta)$")
ax.plot(B,1+B**2/(2+B)**2,color=PAL['teal'],lw=1.3,ls="-.",label=r"SMM at matched speed, $1+\beta^2/(2+\beta)^2$")
ax.axhline(5/3,color=PAL['coral'],lw=0.8,ls=":")
ax.text(0.015,5/3*1.07,r"$5/3$",color=PAL['coral'],fontsize=8)
ks=sorted(sims); ax.plot(ks,[sims[k][0] for k in ks],"o",color=PAL['coral'],ms=5,mfc="white",mew=1.2,label=r"simulated SMM, $a=0.01$")
ax.plot(ks,[sims[k][1] for k in ks],"s",color=PAL['violet'],ms=5,mfc="white",mew=1.2,label=r"simulated NAG, $a=0.01$")
ax.axvline(beta0,color=PAL['grey'],lw=0.8); ax.text(beta0+0.012,22,r"$\log_2 3-1$",color=PAL['ink'],fontsize=8)
ax.set_yscale("log"); ax.set_xlabel(r"momentum coefficient $\beta$"); ax.set_ylabel("variance relative to GD")
ax.set_xlim(0,1); ax.set_ylim(0.9,40); ax.legend(loc="upper left",fontsize=7.2)
fig.savefig(BOOK + "/figures/ch7_amplification.pdf")
# kappa=25 experiment predictions (stationary objective gap ratio), eigenvalues geomspace(1,25,25), alpha = 1/25
lam=np.geomspace(1,25,25); al=1/25; A_=al*lam
def gap(vfun): return 0.5*np.sum(lam*al**2*np.array([vfun(x) for x in A_]))
g_gd=gap(v_gd)
print("kappa25 predicted gap / sigma^2: GD", g_gd)
for name,f in [("SMM",lambda x: v_smm(x,beta0)[2]),("NAG585",lambda x: v_nag(x,beta0)),("NAG05",lambda x: v_nag(x,0.5)),("NAG667",lambda x: v_nag(x,2/3)),("NAG09",lambda x: v_nag(x,0.9)),("HB585",lambda x: v_hb(x,beta0)),("HB09",lambda x: v_hb(x,0.9))]:
    print(name, "ratio to GD", gap(f)/g_gd)
k=67.4; bt=(np.sqrt(k)-1)/(np.sqrt(k)+1); print("tuned NAG beta at kappa 67.4:",bt,"1/(1-b)",1/(1-bt))
