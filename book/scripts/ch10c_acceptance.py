import os
BOOK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the book folder
os.makedirs(BOOK + '/scripts/data', exist_ok=True)  # intermediate results that later scripts load
# Direct test of coalescence beyond the horizon: for every integer 2 <= n < 2^B with sigma(n) > m and
# every bit position i < B, flip bit i and record the fraction of these proposals that stay in the
# level set {sigma > m} (exact enumeration, not sampling; flips to 0 or 1 are not proposals).  Compare
# with the acceptance rate Z(m) of rejection sampling and with e^{DB} Z(m).  Also: the hazard of the tau ladder.
import numpy as np, math, json
B = 26; M = 1 << B
tau = np.load(BOOK + '/scripts/data/ch3_tau_small.npy')
sig = tau[2:].astype(np.int64); cnt = np.bincount(sig, minlength=602)
Z = (cnt.sum() - np.cumsum(cnt))[:601] / cnt.sum()
L2, L3 = math.log(2), math.log(3); q = L2 / L3
D = L2 - (-q * math.log(q) - (1 - q) * math.log(1 - q))
ms = list(range(20, 460, 20)); acc = []; acc_hi = []; acc_lo = []
s16 = tau.astype(np.int16)                       # s16[n] = sigma(n), indexed by n itself
for m in ms:
    pool = np.nonzero(s16[2:] > m)[0].astype(np.int64) + 2
    if len(pool) == 0: break
    cnt = {'lo': [0, 0], 'hi': [0, 0]}            # [proposals, accepted]
    for i in range(B):
        ys = pool ^ (1 << i); ys = ys[ys >= 2]
        g = 'lo' if i < B // 2 else 'hi'
        cnt[g][0] += len(ys); cnt[g][1] += int(np.count_nonzero(s16[ys] > m))
    stay = (cnt['lo'][1] + cnt['hi'][1]) / (cnt['lo'][0] + cnt['hi'][0])
    s_lo = cnt['lo'][1] / cnt['lo'][0]; s_hi = cnt['hi'][1] / cnt['hi'][0]
    acc.append(float(stay)); acc_hi.append(float(s_hi)); acc_lo.append(float(s_lo))
    print(f"m={m:4d}  |level set|={len(pool):9d}  Z={Z[m]:.3e}  accept={stay:.5f}  high bits {s_hi:.6f}  low bits {s_lo:.5f}  accepted flips {cnt['lo'][1] + cnt['hi'][1]}  accept/Z={stay/Z[m]:.2f}  e^(DB)={math.exp(D*B):.2f}")
ms = ms[:len(acc)]
# hazard of the tau ladder from the exact Z_tau
p = json.load(open(BOOK + '/scripts/data/ch10c_perfect.json')); Zt = np.array(p['Z'])
haz = 1 - Zt[1:] / Zt[:-1]; rungs = np.nonzero(haz > 0)[0] + 1
print("rungs up to 400:", len(rungs), " mean hazard on rungs (k>100):", haz[rungs[rungs > 100] - 1].mean(),
      " prediction 1-exp(-D log2 3)=", 1 - math.exp(-D * L3 / L2))
json.dump({'ms': ms, 'acc': acc, 'acc_hi': acc_hi, 'acc_lo': acc_lo, 'Z': [float(Z[m]) for m in ms], 'D': D, 'B': B,
           'rungs': rungs.tolist(), 'hazard': haz[rungs - 1].tolist()},
          open(BOOK + '/scripts/data/ch10c_acceptance.json', 'w'))
