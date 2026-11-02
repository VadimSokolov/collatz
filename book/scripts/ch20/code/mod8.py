"""Residues modulo 8 of the odd iterates along long Syracuse orbits, and the
Proposition 20.9 margin: kbar >= f3*1 + f1*2 + f5*3 + f7*1 must stay above log2 3.
Classes: y = 3,7 mod 8 give k = 1; y = 1 mod 8 gives k = 2; y = 5 mod 8 gives k >= 3."""
import random, math, sys
TH = math.log2(3)
def run(y, label, burn=0):
    cnt = {1:0,3:0,5:0,7:0}; n = 0; ksum = 0; worst_dev = 0.0; worst_lb = 9.0
    while y != 1:
        r = y & 7; cnt[r] += 1; n += 1
        z = 3*y + 1; k = (z & -z).bit_length() - 1; ksum += k; y = z >> k
        if n >= 1000 and n % 97 == 0:
            f = {c: cnt[c]/n for c in cnt}
            dev = max(abs(f[c]-0.25) for c in f)
            lb = f[3] + f[7] + 2*f[1] + 3*f[5]
            worst_dev = max(worst_dev, dev); worst_lb = min(worst_lb, lb)
    f = {c: cnt[c]/n for c in cnt}
    print(f"{label:>28s}: strides={n:>8d} kbar={ksum/n:.5f} freq(1,3,5,7 mod 8)="
          f"({f[1]:.4f},{f[3]:.4f},{f[5]:.4f},{f[7]:.4f}) max|f-1/4| after 1000={worst_dev:.4f} "
          f"min lower bound={worst_lb:.4f}")
random.seed(20)
run(3**3155, "3^3155")
for bits in (20000, 50000, 100000):
    y = random.getrandbits(bits) | 1 | (1 << (bits-1))
    run(y, f"random {bits} bit odd")
y = (1 << 3000) - 1
run(y, "2^3000 - 1 (long ascent)")
