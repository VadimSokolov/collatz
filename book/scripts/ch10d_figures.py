# Road map figure for the synthesis chapter: two lanes (cycles, divergence) from a shared kernel
# to the conjecture, with the two mathematical gates marked.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bookstyle import *
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
fig, ax = plt.subplots(figsize=(WIDTH, 5.0))
ax.set_xlim(-0.1, 10.5); ax.set_ylim(0, 12.6); ax.axis('off')

def box(x, y, w, h, text, col, fs=5.8, alpha=0.16, ls='-'):
    p = FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0.02,rounding_size=0.18',
                       fc=col, ec=col, lw=1.0, ls=ls)
    p.set_facecolor(mpl.colors.to_rgba(col, alpha))
    ax.add_patch(p)
    ax.text(x, y, text, ha='center', va='center', fontsize=fs, color=PAL['ink'], linespacing=1.25)

def arrow(x0, y0, x1, y1, col=PAL['ink'], ls='-', lw=0.9):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=7,
                                 color=col, lw=lw, ls=ls, shrinkA=0, shrinkB=0))

T, C, D, G, K = PAL['teal'], PAL['coral'], PAL['violet'], PAL['gold'], PAL['ink']
# kernel
box(4.85, 11.85, 9.6, 1.25, 'Stage 0: the kernel\nLean statements, exact census to $2^{30}$,\n'
    'stopping count recursion, perfect sampler', K, alpha=0.08)
ax.text(2.2, 10.72, 'cycle lane', ha='right', fontsize=7.4, color=T, style='italic')
ax.text(7.2, 10.72, 'divergence lane', ha='right', fontsize=7.4, color=C, style='italic')
# cycle lane
xc, xd, w, h = 2.35, 7.35, 4.6, 1.12
cyc = [('C1  verification past $2^{71.88}$\n$r\\geq137\\,528\\,045\\,312$', 9.6),
       ('C2  Lean replay of the 44 rung\ninteger certificate', 8.15),
       ('C3  $m$ cycles beyond $m=91$\nBaker bounds and lattice reduction', 6.7),
       ('C4  divisibility data: residues of\n$c_k(w)$ modulo factors of $2^k-3^m$', 5.25)]
for t, y in cyc: box(xc, y, w, h, t, T)
div = [('D1  hard instances: exceptional\nclasses at depth $k$, slow trees', 9.6),
       ('D2  a weighted automaton\nreading binary and ternary digits', 8.15),
       ('D3  SAT search for matrix\ninterpretations, sampler in the loop', 6.7),
       ('D4  independent check of the\ncertificate in Lean', 5.25)]
for t, y in div: box(xd, y, w, h, t, C)
for x in (xc, xd):
    arrow(x, 11.2, x, 10.16)
    for y0, y1 in ((9.6, 8.15), (8.15, 6.7), (6.7, 5.25)):
        arrow(x, y0 - h / 2, x, y1 + h / 2)
# loop arrow D3 -> D1 (counterexamples feed back)
ax.add_patch(FancyArrowPatch((xd + w / 2, 6.7), (xd + w / 2, 9.6), connectionstyle='arc3,rad=0.32',
                             arrowstyle='-|>', mutation_scale=7, color=C, lw=0.8, ls='--'))
ax.text(10.45, 8.15, 'counter\nexamples', fontsize=5.6, color=C, ha='center', va='center', rotation=90)
# gates
box(xc, 3.55, w, 1.15, 'Gate A, mathematics: an arithmetic\nmargin of order $m/3$ in the window,\nor no word with $2^k-3^m\\mid c_k(w)$', G, alpha=0.22)
box(xd, 3.55, w, 1.15, 'Gate B, search: a termination\ncertificate for the mixed base\nrewriting system', G, alpha=0.22)
arrow(xc, 5.25 - h / 2, xc, 3.55 + 0.575); arrow(xd, 5.25 - h / 2, xd, 3.55 + 0.575)
# outcomes
box(xc, 1.75, w, 0.95, 'no nontrivial cycle', T, alpha=0.10)
box(5.6, 0.45, 5.4, 0.8, 'the conjecture, machine checked', K, alpha=0.10)
arrow(xc, 3.55 - 0.575, xc, 1.75 + 0.475)
arrow(xd, 3.55 - 0.575, 6.3, 0.45 + 0.4)
ax.text(8.55, 1.95, 'settles both halves', fontsize=6.0, color=C, ha='center', style='italic')
arrow(xc + 0.6, 1.75 - 0.475, 3.6, 0.45 + 0.4, ls=':')
ax.text(1.25, 0.62, 'divergence half\nstill open', fontsize=5.6, color=T, ha='center', style='italic')
fig.savefig(os.path.join(FIG, 'ch10d_roadmap.pdf')); plt.close(fig)
print('saved')
