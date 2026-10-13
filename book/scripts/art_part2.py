"""Part II opener: stationary noise clouds of three first order methods.

Quadratic f(x) = x^T A x / 2, eigenvalues 1 and kappa = 20, eigenvectors
rotated by 72 degrees; step alpha = 1/L = 1/20; the gradient is observed with
additive isotropic Gaussian noise (sigma = 1).  Methods, all from x* = 0:
  gradient descent : x+ = x - alpha g(x)
  Nesterov         : y = x + beta (x - x_prev), x+ = y - alpha g(y), every step
  Switching MM     : Nesterov step on odd t (injection), gradient step on even t
with beta = log2 3 - 1 = 0.585 for both momentum methods.
Dots: 2500 iterates of each method sampled from a long run after burn in.
Curves: exact two standard deviation ellipses of the (phase averaged)
stationary covariance, from discrete Lyapunov equations.
Background: filled contours of log10 f.  The number under each panel is the
stationary variance of the slow mode relative to gradient descent.
"""
from art_common import *
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Ellipse
from scipy.linalg import solve_discrete_lyapunov as dlyap

rng = np.random.default_rng(20261010)
kappa, ang = 20.0, np.deg2rad(72)
R = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
lam = np.array([1.0, kappa])
A = R @ np.diag(lam) @ R.T
alpha = 1 / kappa
beta = np.log2(3) - 1
sigma = 1.0


def mode_var(a, method):
    """Exact stationary variance (phase averaged) of one eigenmode e_t, a = alpha*lambda."""
    b0 = np.array([[-alpha * sigma], [0.0]])
    Mn = np.array([[(1 - a) * (1 + beta), -(1 - a) * beta], [1, 0]])
    Me = np.array([[1 - a, 0], [1, 0]])
    if method == 'gd':
        return alpha**2 * sigma**2 / (1 - (1 - a) ** 2)
    if method == 'nag':
        return dlyap(Mn, b0 @ b0.T)[0, 0]
    # SMM: even step Me then odd step Mn
    Phi = Mn @ Me
    Q = Mn @ b0 @ b0.T @ Mn.T + b0 @ b0.T
    Pe = dlyap(Phi, Q)
    Po = Me @ Pe @ Me.T + b0 @ b0.T
    return 0.5 * (Pe[0, 0] + Po[0, 0])


# check against the book's table (beta = 0.585, a = 0.01: NAG 2.372, SMM 1.355)
gd01 = mode_var(0.01, 'gd')
print('check a=0.01: NAG/GD', mode_var(0.01, 'nag') / gd01, ' SMM/GD', mode_var(0.01, 'smm') / gd01)

cov = {}
for m in ('gd', 'nag', 'smm'):
    v = np.array([mode_var(alpha * l, m) for l in lam])
    cov[m] = R @ np.diag(v) @ R.T
    print(m, 'slow mode var', v[0], 'fast', v[1], 'ratio to GD slow', v[0] / mode_var(alpha, 'gd'))


def simulate(method, T=400000, burn=2000, every=157):
    x = np.zeros(2); xp = np.zeros(2); out = []
    xi = rng.normal(0, sigma, (T, 2))
    for t in range(T):
        if method == 'gd' or (method == 'smm' and t % 2 == 0):
            xn = x - alpha * (A @ x + xi[t])
        else:
            y = x + beta * (x - xp)
            xn = y - alpha * (A @ y + xi[t])
        xp, x = x, xn
        if t > burn and t % every == 0:
            out.append(x)
    return np.array(out)


pts = {m: simulate(m) for m in ('nag', 'smm', 'gd')}
for m in pts:
    print(m, 'n', len(pts[m]), 'empirical cov trace / exact', np.trace(np.cov(pts[m].T)) / np.trace(cov[m]))

LAND = LinearSegmentedColormap.from_list(
    'land', ['#07061A', '#14123A', '#1B1F5E', '#24407A', '#1F6F86', '#2E9C93'])
W, H = 4.5, 3.5
fig = plt.figure(figsize=(W, H), facecolor=NIGHT)
hy = 0.82
gap = 0.012
pw = (1 - 2 * gap) / 3
order = ('gd', 'smm', 'nag')
style = {'nag': (PAL['coral'], 'Nesterov'), 'smm': ('#7FE0FF', 'Switching MM'),
         'gd': (PAL['gold'], 'gradient descent')}
for i, m in enumerate(order):
    ax = fig.add_axes([i * (pw + gap), 0, pw, 1]); ax.set_facecolor(NIGHT); ax.axis('off')
    hx = hy * (pw * W) / H
    ax.set_xlim(-hx, hx); ax.set_ylim(-hy, hy)
    gx, gy = np.meshgrid(np.linspace(-1.03 * hx, 1.03 * hx, 300), np.linspace(-1.02 * hy, 1.02 * hy, 700))
    F = 0.5 * (A[0, 0] * gx**2 + 2 * A[0, 1] * gx * gy + A[1, 1] * gy**2)
    lv = np.log10(F + 1e-5)
    levels = np.linspace(-3.0, 1.3, 22)
    cf = ax.contourf(gx, gy, lv, levels=levels, cmap=LAND.reversed(), extend='both')
    cf.set_rasterized(True)
    ax.contour(gx, gy, lv, levels=levels, colors='#9FD8E8', linewidths=0.3, alpha=0.3)
    c, lab = style[m]
    P = pts[m][:2500]
    sc = ax.scatter(P[:, 0], P[:, 1], s=2.2, color=c, alpha=0.6, lw=0)
    sc.set_rasterized(True)
    w, V = np.linalg.eigh(cov[m])
    angd = np.degrees(np.arctan2(V[1, 1], V[0, 1]))
    for lw_, ec in ((1.4, c), (0.35, 'white')):
        ax.add_patch(Ellipse((0, 0), 4 * np.sqrt(w[1]), 4 * np.sqrt(w[0]), angle=angd,
                             fill=False, ec=ec, lw=lw_))
    ax.plot(0, 0, marker='*', ms=7, mfc='white', mec=NIGHT, mew=0.5, zorder=6)
    ax.text(0, -hy * 0.9, lab, ha='center', va='center', color='white', fontsize=8.5,
            fontstyle='italic')
    ax.text(0, -hy * 0.9 - 0.075, r'$\times\,%.2f$' % (mode_var(alpha, m) / mode_var(alpha, 'gd')),
            ha='center', va='center', color=c, fontsize=7.5)
save(fig, 'art_part2', facecolor=NIGHT, bbox_inches=None)
