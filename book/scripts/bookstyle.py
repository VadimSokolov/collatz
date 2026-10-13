"""Shared figure style for the book. Usage:
    import sys; sys.path.insert(0, BOOK + '/scripts'); from bookstyle import *
"""
import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
import matplotlib as mpl
import matplotlib.pyplot as plt

PAL = dict(teal='#0B6E79', coral='#E8573F', gold='#F2A541', violet='#6C4AB6',
           sky='#3FA7D6', green='#59A96A', plum='#B5407A', ink='#22223B',
           sand='#F4EDE4', grey='#8D99AE')
CYCLE = [PAL[k] for k in ('teal', 'coral', 'violet', 'gold', 'sky', 'green', 'plum')]
CMAP_SEQ = 'viridis'      # sequential
CMAP_HOT = 'magma'        # dramatic sequential
CMAP_DIV = 'RdYlBu_r'     # diverging
CMAP_CYC = 'twilight'     # cyclic (angles, residues)

mpl.rcParams.update({
    'figure.dpi': 150, 'savefig.dpi': 400, 'savefig.bbox': 'tight',   # figures are drawn WIDTH wide and printed up to 1.32 times wider, so rasterized parts (colour bars, dense scatters) need 400 dpi for 300 ppi
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'mathtext.fontset': 'dejavuserif',
    'font.size': 9, 'axes.titlesize': 9, 'axes.labelsize': 9, 'legend.fontsize': 8,
    'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'axes.prop_cycle': mpl.cycler(color=CYCLE),
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': PAL['ink'], 'axes.labelcolor': PAL['ink'],
    'xtick.color': PAL['ink'], 'ytick.color': PAL['ink'],
    'axes.grid': True, 'grid.color': '#E3E3E3', 'grid.linewidth': 0.5,
    'legend.frameon': False, 'lines.linewidth': 1.4,
})
WIDTH = 4.5   # inches, full text width
