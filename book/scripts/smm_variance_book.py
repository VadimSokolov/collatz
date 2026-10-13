# Exact stationary variance of Switching MM on one eigendirection of a
# quadratic, a = alpha * lambda, unit noise (alpha^2 sigma^2 = 1).
import sympy as sp
a, b = sp.symbols('a b', positive=True)
rho = (1 - a) * (1 - (1 + b) * a)            # two step contraction
q = (1 - a)**2 * (1 + b)**2 + 1               # two step noise variance
v_even = q / (1 - rho**2)
# odd phase: one descent step after an even time
v_odd = (1 - a)**2 * v_even + 1
v_avg = sp.simplify((v_even + v_odd) / 2)
slow = sp.simplify(sp.limit(a * v_avg, a, 0))  # GD gives 1/2
print("slow mode a*Var:", slow)
V = sp.simplify(2 * slow)
print("amplification vs GD:", V, " at beta = log2(3) - 1:",
      float(V.subs(b, sp.log(3, 2) - 1)), " limit beta -> 1:", V.subs(b, 1))
