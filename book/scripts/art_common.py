"""Shared helpers for the art_* figures (cover, part openers, gallery).

Imports the book style, switches the serif font to TeX Gyre Pagella (a Palatino
clone matching the mathpazo text face), and provides fast Collatz utilities.
"""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import sys, os
sys.path.insert(0, BOOK + '/scripts')
from bookstyle import *          # noqa: F401,F403
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numba

FIG = BOOK + '/figures'
GYRE = '/usr/share/texmf/fonts/opentype/public/tex-gyre/'
for f in ('texgyrepagella-regular.otf', 'texgyrepagella-italic.otf'):
    p = GYRE + f
    if os.path.exists(p):
        font_manager.fontManager.addfont(p)
mpl.rcParams.update({'font.serif': ['TeX Gyre Pagella', 'DejaVu Serif'],
                     'mathtext.fontset': 'custom',
                     'mathtext.rm': 'TeX Gyre Pagella',
                     'mathtext.it': 'TeX Gyre Pagella:italic',
                     'mathtext.cal': 'TeX Gyre Pagella:italic',
                     'pdf.fonttype': 42})

INK = PAL['ink']
NIGHT = '#0B0A1A'      # deep background for dark plates


@numba.njit(cache=True)
def total_stopping(N):
    """sigma_inf(n) (number of 3n+1 and n/2 steps to reach 1) and the orbit
    maximum M(n), for 1 <= n <= N."""
    st = np.zeros(N + 1, np.int32)
    mx = np.zeros(N + 1, np.int64)
    for n in range(2, N + 1):
        x = n; s = 0; m = n
        while x >= n and x != 1:
            if x & 1:
                x = 3 * x + 1
            else:
                x >>= 1
            if x > m:
                m = x
            s += 1
        if x < n:              # fell below n, reuse the stored values
            st[n] = s + st[x]
            mx[n] = max(m, mx[x])
        else:
            st[n] = s
            mx[n] = m
    mx[1] = 4  # orbit 1 -> 4 -> 2 -> 1
    return st, mx


def collatz_orbit(n):
    out = [n]
    while n != 1:
        n = 3 * n + 1 if n & 1 else n // 2
        out.append(n)
    return out


def save(fig, name, png=True, dpi=400, **kw):  # the plates are drawn 4.5in wide and printed wider, so 400 dpi keeps them above 300 ppi
    fig.savefig(f'{FIG}/{name}.pdf', dpi=dpi, **kw)
    if png:
        fig.savefig(f'{FIG}/{name}.png', dpi=110, **kw)
    plt.close(fig)
