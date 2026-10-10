"""Importance weights correct the bias of sequential stub matching.

Four hubs of degree 12 among 44 vertices of degree 2. Statistic: number of hub-hub edges.
Compare the exact (rejection) sampler, the raw sequential sampler, and the self-normalized
importance-weighted estimate, with delta-method standard errors and the effective sample size.

    PYTHONPATH=python python3 examples/importance_weights.py
"""
import numpy as np

import hypergraph_sampler as hs

deg = [12] * 4 + [2] * 44
hub = lambda e: int(((e[:, 0] < 4) & (e[:, 1] < 4)).sum())
R = 20000

exact = np.array([hub(hs.config_degrees(deg, 2, seed=s)) for s in range(R)])
f, lw = np.zeros(R), np.zeros(R)
for s in range(R):
    e, lw[s] = hs.sequential(deg, 2, seed=10**6 + s)
    f[s] = hub(e)
ok = np.isfinite(lw)
w = np.zeros(R)
w[ok] = np.exp(lw[ok] - lw[ok].max())          # stuck runs keep weight 0
est = (w * f).sum() / w.sum()
se = np.sqrt((w**2 * (f - est) ** 2).sum()) / w.sum()
ess = w.sum() ** 2 / (w**2).sum()

print(f"exact sampler       {exact.mean():.3f} +- {exact.std() / np.sqrt(R):.3f}")
print(f"raw sequential      {f[ok].mean():.3f} +- {f[ok].std() / np.sqrt(ok.sum()):.3f}   (stuck {100 * (~ok).mean():.0f}%)")
print(f"weighted sequential {est:.3f} +- {se:.3f}   (effective sample size {ess:.0f} of {R})")
