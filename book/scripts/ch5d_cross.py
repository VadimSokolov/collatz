import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ch5d_pep import pep, smm, cyclic_gd, similar_triangles
out = []
for N in (24, 28, 32, 36, 40):
    g = pep(cyclic_gd([1.0] * N)); s = pep(smm(N))
    Hs, As = similar_triangles(N, [t % 2 == 0 for t in range(N)]); w = pep(Hs)
    Hf, Af = similar_triangles(N, [True] * N); f = pep(Hf)
    out.append((N, g, s, w, f, 1 / (2 * As), 1 / (2 * Af))); print(out[-1], flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'ch5d_cross.json'), 'w'))
