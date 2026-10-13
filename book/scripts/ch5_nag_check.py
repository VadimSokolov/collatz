import numpy as np
def rad(beta, a):
    # NAG on scalar mode a = alpha*lambda: z^2 - (1-a)(1+beta) z + (1-a) beta
    return max(abs(np.roots([1, -(1-a)*(1+beta), (1-a)*beta])))
def hb(beta, a):
    return max(abs(np.roots([1, -(1+beta-a), beta])))
for kappa in (4, 14.58, 100):
    A = np.linspace(1/kappa, 1, 400)
    bs = np.linspace(0, 0.99, 2000)
    worst = [max(rad(b, a) for a in A) for b in bs]
    i = int(np.argmin(worst))
    print(kappa, "best beta", bs[i], "worst rad", worst[i], "formula", (np.sqrt(kappa)-1)/(np.sqrt(kappa)+1), 1-1/np.sqrt(kappa))
