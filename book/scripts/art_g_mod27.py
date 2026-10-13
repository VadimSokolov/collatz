"""Gallery: the model chain on the 18 units modulo 27.

For odd y with y = r (mod 27) and k = v2(3y+1), S(y) = (3y+1)/2^k is congruent
to (3r+1) 2^{-k} modulo 27, and k is Geometric(1/2) independently of r (the
one step matrix P_3 of chap:divergence, which is exact).  Hence
P(r -> s) = sum over k >= 1 with (3r+1) 2^{-k} = s (mod 27) of 2^{-k}
          = 2^{-k0} / (1 - 2^{-18}),  k0 in 1..18,
since 2 has order 18 modulo 27.  Units are placed round the circle in the
order 2^0, 2^1, ..., 2^17 (2 is a primitive root mod 27), counterclockwise, so
multiplying by 2^{-1} is one step clockwise: the out edges of r fan clockwise
from the node 3r+1.  Edges with k = 1, 2, 3 are drawn (mass 7/8 per
row; loops as rings outside the node), width and opacity proportional to probability; node area proportional
to the stationary law pi_3.
"""
from art_common import *
from art_coral import EMBER
from matplotlib.patches import FancyArrowPatch, Circle

M = 27
units = [pow(2, j, M) for j in range(18)]
pos = {u: j for j, u in enumerate(units)}
inv2 = pow(2, -1, M)
P = np.zeros((18, 18))   # P[i, j] = prob from units[i] to units[j]
for i, r in enumerate(units):
    for k in range(1, 200):
        s = (3 * r + 1) * pow(inv2, k, M) % M
        P[i, pos[s]] += 2.0 ** -k
assert np.allclose(P.sum(1), 1)
w, V = np.linalg.eig(P.T)
pi = np.real(V[:, np.argmin(abs(w - 1))]); pi /= pi.sum()
ev = np.sort(abs(w))[::-1]
print('eigenvalue moduli', np.round(ev[:4], 6))
print('pi_3 x 63*9? :', {u: round(pi[pos[u]], 5) for u in sorted(units)})
# check projection to mod 9 against pi_2 = (8,16,11,4,2,22)/63
pi2 = {a: sum(pi[pos[u]] for u in units if u % 9 == a) for a in (1, 2, 4, 5, 7, 8)}
print('projection mod 9 x63:', {a: round(63 * v, 4) for a, v in pi2.items()})
from fractions import Fraction
# exact pi_3 via rational arithmetic: power iteration 3 steps gives exact rank one
Pf = [[Fraction(0)] * 18 for _ in range(18)]
den = Fraction(1) - Fraction(1, 2**18)
for i, r in enumerate(units):
    for k in range(1, 19):
        s = (3 * r + 1) * pow(inv2, k, M) % M
        Pf[i][pos[s]] += Fraction(1, 2**k) / den
row = [Fraction(1, 18)] * 18
for _ in range(3):
    row = [sum(row[i] * Pf[i][j] for i in range(18)) for j in range(18)]
print('exact pi_3:', {units[j]: str(row[j]) for j in range(18)})
print('max pi', max(row), 'at', units[int(np.argmax([float(x) for x in row]))],
      'min pi', min(row), 'at', units[int(np.argmin([float(x) for x in row]))])

fig = plt.figure(figsize=(4.5, 4.5), facecolor=NIGHT)
ax = fig.add_axes([0, 0, 1, 1]); ax.set_facecolor(NIGHT); ax.axis('off'); ax.set_aspect('equal')
ax.set_xlim(-1.36, 1.36); ax.set_ylim(-1.36, 1.36)
th = np.pi / 2 + 2 * np.pi * np.arange(18) / 18
X, Y = np.cos(th), np.sin(th)
srccol = {a: EMBER(c) for a, c in zip((1, 2, 4, 5, 7, 8), (0.28, 0.42, 0.55, 0.68, 0.8, 0.92))}
edges = []
for i, r in enumerate(units):
    for k in range(1, 4):
        s = (3 * r + 1) * pow(inv2, k, M) % M
        edges.append((2.0 ** -k, i, pos[s], r))
edges.sort()
for p, i, j, r in edges:
    if i == j:
        rl = 0.075
        ax.add_patch(Circle((X[i] * (1.0 + 0.03 + 0.22 * np.sqrt(float(row[i])) + rl * 0.6),
                             Y[i] * (1.0 + 0.03 + 0.22 * np.sqrt(float(row[i])) + rl * 0.6)), rl,
                            fill=False, ec=srccol[r % 9], lw=0.3 + 3.2 * p, alpha=0.25 + 0.7 * p, zorder=3))
        continue
    a = FancyArrowPatch((X[i], Y[i]), (X[j], Y[j]), connectionstyle='arc3,rad=0.28',
                        arrowstyle='-|>,head_length=%.2f,head_width=%.2f' % (2.2 + 6 * p, 1.2 + 3 * p),
                        lw=0.3 + 3.2 * p, color=srccol[r % 9], alpha=0.25 + 0.7 * p,
                        shrinkA=7, shrinkB=9, zorder=2)
    ax.add_patch(a)
for i, r in enumerate(units):
    for k in range(1, 6):
        if (3 * r + 1) * pow(inv2, k, M) % M == r:
            print('self loop at', r, 'k', k)
size = np.array([float(row[j]) for j in range(18)])
for j, u in enumerate(units):
    rad = 0.03 + 0.22 * np.sqrt(size[j])
    ax.add_patch(Circle((X[j], Y[j]), rad, facecolor=srccol[u % 9], edgecolor='#FFF3C4', lw=0.6, zorder=4))
    lr = 1.0 + rad + 0.08 + (0.17 if any(e[1] == j and e[2] == j for e in edges) else 0)
    ax.text(lr * X[j], lr * Y[j], f'{u}', ha='center', va='center', fontsize=8, color='#FFF3C4')
ax.text(0, 0.06, r'units modulo $27$', ha='center', va='center', fontsize=9, color='#E6E6EE', fontstyle='italic')
ax.text(0, -0.07, r'$r\;\mapsto\;(3r+1)\,2^{-k}$', ha='center', va='center', fontsize=8.5, color='#C9CCD8')
save(fig, 'art_g_mod27', facecolor=NIGHT, bbox_inches=None)
