import sympy as sp
a,b,x=sp.symbols('a b x')
D=sp.Matrix([[1-a,0],[1,0]]); N=sp.Matrix([[(1-a)*(1+b),-(1-a)*b],[1,0]])
M=N*D; cp=sp.factor(M.charpoly(x).as_expr()); print("charpoly ND:",cp)
# slow eigenvalue series in a
lam=sp.symbols('lam'); 
ev=sp.solve(sp.Eq(M.charpoly(x).as_expr(),0),x)
for e in ev:
    s=sp.series(e,a,0,3).removeO(); print(sp.simplify(s))
