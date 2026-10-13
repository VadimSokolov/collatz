"""Chapter 3: exact mean and variance of the passage time of the Bernoulli model,
for the comparisons with the map in tab:ch3:clt, fig:ch3:paths and the census
of sec:ch3:missing.

tau is the first n with log x_n = u + sum Delta_i <= 0, Delta = log(3/2) or
-log 2 with probability 1/2 each.  For a given u = log x_0 the law of tau is
computed exactly by dynamic programming over the number of up steps
(passage_law of ch3_o4_model.py, truncated when the surviving mass is below
1e-17).  For a range [a, b) of starting values, u is distributed as log n with
n uniform on the range (density proportional to e^u on [log a, log b)), and the
law of tau is the mixture over NPTS midpoints in u.  Printed: E log n, the
leading terms gamma_0 E log n and (sigma^2/|delta|^3) E log n, the exact mean
and variance of tau, the excess of the mean over gamma_0 E log n, and the ratio
of the variance to its leading term."""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import numpy as np, sys
sys.path.insert(0, BOOK + '/scripts')
from ch3_o4_model import passage_law

d = 0.5*np.log(4/3)            # |delta^*|
s2 = 0.25*np.log(3)**2          # sigma^2
g0 = 1/d                        # gamma_0 = 6.952
cv = s2/d**3                    # 101.387


def range_moments(a, b, npts=400):
    """Mean and variance of tau when x_0 = n is uniform on [a, b)."""
    us = np.log(a) + (np.arange(npts) + 0.5)/npts*np.log(b/a)
    w = np.exp(us); w /= w.sum()
    m1 = m2 = 0.0
    for wi, u in zip(w, us):
        p = passage_law(u); t = np.arange(len(p)); p = p/p.sum()
        m1 += wi*(t*p).sum(); m2 += wi*(t*t*p).sum()
    return (w*us).sum(), m1, m2 - m1*m1


def report(label, a, b, npts):
    eu, m, v = range_moments(a, b, npts)
    print('%-22s E log n %8.4f | gamma0 E log n %8.3f  mean %8.3f  excess %6.3f | '
          'lead var %8.1f  var %8.1f  ratio %.4f' % (label, eu, g0*eu, m, m - g0*eu, cv*eu, v, v/(cv*eu)), flush=True)


if __name__ == '__main__':
    # the five ranges of tab:ch3:clt
    for j in range(3, 8):
        report('[10^%d, 2*10^%d)' % (j, j), 10.0**j, 2*10.0**j, 400)
    # dyadic blocks of the census, tab:ch3:missing
    for k in (10, 15, 20, 25, 29):
        report('[2^%d, 2^%d)' % (k, k + 1), 2.0**k, 2.0**(k + 1), 400)
    # larger blocks: the excess tends to the stationary overshoot constant
    for k in (50, 100, 200):
        report('[2^%d, 2^%d)' % (k, k + 1), 2.0**k, 2.0**(k + 1), 100)
    # fig:ch3:paths: model paths started at x_0 = 10^12 exactly
    u = 12*np.log(10.0)
    p = passage_law(u); t = np.arange(len(p)); p = p/p.sum()
    m = (t*p).sum(); v = (t*t*p).sum() - m*m
    print('x_0 = 10^12: gamma0 log x_0 %.2f  mean %.2f | leading sd %.2f  sd %.2f'
          % (g0*u, m, np.sqrt(cv*u), np.sqrt(v)), flush=True)
