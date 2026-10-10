"""Acceptance of the configuration model: mean trials vs the asymptotic 1/P(simple).

P(simple) -> exp(-(k-1)(d-1)/2) for k >= 3 (repeated vertices only), and
exp(-(d-1)/2 - (d-1)^2/4) for graphs, where repeated edges also stay Poisson.

    PYTHONPATH=python python3 examples/acceptance.py
"""
import numpy as np

import hypergraph_sampler as hs

reps = 40


def predicted(k, d):
    lam = (k - 1) * (d - 1) / 2 + ((d - 1) ** 2 / 4 if k == 2 else 0.0)
    return np.exp(lam)


print(" k  d      n   mean trials   asymptotic 1/P(simple)")
for k, d in [(2, 3), (2, 5), (3, 2), (3, 4), (3, 6), (3, 8)]:
    for n in [60, 600, 3000]:
        if (n * d) % k:
            continue
        t = [hs.config_regular(n, k, d, seed=s, return_trials=True)[1] for s in range(reps)]
        print(f"{k:2d} {d:2d} {n:6d}   {np.mean(t):11.1f}   {predicted(k, d):11.1f}")
