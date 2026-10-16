"""The exact 3 adic offset law for the chapter on positive density from the backward tree (adapted from N. Polson's
offset_law_fourier.py, offset_law_fourier_profile.py and offset_law_level18.py). Letters a >= 1 with P(a) = 2^-a, offsets
off(a v) = 2^-a (1 + 3 off(v)); the law mu_d of the offset modulo 3^d is computed from mu_{d-1} exactly, with the letters
truncated at a <= 60 (mass below 2^-60) and renormalized. For each level d <= DMAX (the first argument, 18 by default): the
largest Fourier coefficient over all frequencies prime to 3 (a full FFT) and the frequency attaining it (as an integer in
[1, 3^d/2]), whether that integer is a power of 2, the mean distances of rho_d = 3^d mu_d from the lifts of rho_{d-1} and
of rho_{ceil(d/2)}, max rho_d, and for d >= 14 the profile |mu_d^(2^k)| over every k with 2^k < 3^d. About 10 minutes and
under 30 GB to d = 18. At d = 19 (3^19 = 1.16e9 entries) the build fits in 30 GB but numpy's FFT ended in a segmentation
fault on Hopper (job 2516663). Writes data/ch9b_offset_law.json after every level."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data'); os.makedirs(DATA, exist_ok=True)
AMAX = 60; DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 18
OUT = os.path.join(DATA, 'ch9b_offset_law.json')

def mean_gap(mu, m):
    """Mean over Z/3^d of |rho_d - rho_m o pr|, where rho_m is the density of mu modulo 3^m; computed in blocks."""
    M = len(mu); V = mu.reshape(-1, 3**m); rm = 3**m * V.sum(0); tot = 0.0; B = max(1, 2**24 // 3**m)
    for r in range(0, V.shape[0], B):
        tot += float(np.abs(M * V[r:r + B] - rm).sum())
    return tot / M

res = []; mu = np.array([1.0])
for d in range(1, DMAX + 1):
    M = 3**d; inv2 = pow(2, -1, M)
    z = (1 + 3 * np.arange(3**(d - 1), dtype=np.int64)) % M
    new = np.zeros(M); w = 1.0
    for a in range(1, AMAX + 1):
        z *= inv2; z %= M; w *= 0.5
        new += np.bincount(z, weights=w * mu, minlength=M)
    del z
    new /= new.sum(); mu = new; del new
    F = np.abs(np.fft.rfft(mu)); F[::3] = -1.0              # frequencies divisible by 3 (and 0) left out
    k = int(np.argmax(F))
    entry = dict(d=d, fmax=float(F[k]), argxi=k, pow2=bool(k & (k - 1) == 0), log2xi=float(np.log2(k)),
                 L1=mean_gap(mu, d - 1), L1half=mean_gap(mu, (d + 1) // 2) if d > 1 else 0.0, rhomax=float(M * mu.max()))
    if d >= 14:
        entry['profile'] = {str(j): float(F[min(2**j, M - 2**j)]) for j in range(1, 2 * d) if 2**j < M}
    del F
    res.append(entry)
    print({q: v for q, v in entry.items() if q != 'profile'}, flush=True)
    json.dump(res, open(OUT, 'w'))
