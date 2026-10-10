"""Statistical tests against exhaustive enumeration on small cases.

Each exactness test enumerates every simple (hyper)graph with the given degrees, draws many
samples, and applies a chi-square goodness-of-fit test against the uniform distribution at a
strict threshold (p > 1e-4), so a correct sampler fails with probability about 1e-4.
"""
import itertools
from collections import Counter

import numpy as np
import pytest
from scipy.stats import binom, chi2

import hypergraph_sampler as hs


def enumerate_simple(n, k, deg):
    """All simple k-uniform hypergraphs on n vertices with degree sequence deg."""
    deg = list(deg)
    m = sum(deg) // k
    ksets = list(itertools.combinations(range(n), k))
    out = []
    for combo in itertools.combinations(ksets, m):
        c = Counter(v for e in combo for v in e)
        if all(c[v] == deg[v] for v in range(n)):
            out.append(frozenset(combo))
    return out


def key(edges):
    return frozenset(tuple(int(x) for x in row) for row in edges)


def chi2_uniform(counts, support):
    N = sum(counts.values())
    assert set(counts) <= set(support), "sampler produced something outside the support"
    exp = N / len(support)
    stat = sum((counts.get(h, 0) - exp) ** 2 / exp for h in support)
    return chi2.sf(stat, len(support) - 1)


CASES = [  # (n, k, degree sequence)
    (6, 2, [3, 3, 2, 2, 2, 2]),   # irregular graphs
    (6, 2, [3] * 6),               # cubic graphs on 6 vertices
    (6, 3, [2] * 6),               # 2-regular 3-graphs on 6 vertices
    (6, 3, [3] * 6),               # 3-regular 3-graphs on 6 vertices (330)
    (7, 2, [2] * 7),               # 2-regular graphs on 7 vertices (465)
]


@pytest.mark.parametrize("n,k,deg", CASES)
def test_config_is_exactly_uniform(n, k, deg):
    support = enumerate_simple(n, k, deg)
    N = 400 * len(support) if len(support) < 200 else 60 * len(support)
    counts = Counter(key(hs.config_degrees(deg, k, seed=s)) for s in range(N))
    assert set(counts) == set(support)
    assert chi2_uniform(counts, support) > 1e-4


def test_config_regular_matches_degrees_and_is_simple():
    H, trials = hs.config_regular(500, 3, 12, seed=1, return_trials=True)
    assert H.shape == (2000, 3) and hs.is_simple(H, 500)
    assert (hs.degrees(H, 500) == 12).all() and trials >= 1
    assert (hs.config_regular(500, 3, 12, seed=1) == H).all()          # reproducible
    G = hs.config_regular(1000, 2, 4, seed=2)
    assert hs.is_simple(G, 1000) and (hs.degrees(G, 1000) == 4).all()


def test_saturation_raises():
    with pytest.raises(hs.SaturationError):
        hs.config_regular(8, 3, 12, seed=1, max_trials=10**4)


@pytest.mark.parametrize("n,k,deg", CASES)
def test_sequential_weights_correct_bias(n, k, deg):
    """Raw sequential output is biased; importance-weighted frequencies are uniform."""
    support = enumerate_simple(n, k, deg)
    index = {h: i for i, h in enumerate(support)}
    N = 400 * len(support) if len(support) < 200 else 60 * len(support)
    lab, w = np.full(N, -1), np.zeros(N)
    for s in range(N):
        e, lw = hs.sequential(deg, k, seed=s)
        if np.isfinite(lw):                      # stuck runs stay in with weight 0
            lab[s] = index[key(e)]               # KeyError would mean output outside the support
            w[s] = np.exp(lw)
    W = w.sum()
    S = len(support)
    for h in range(S):
        ind = (lab == h).astype(float)
        p = (w * ind).sum() / W
        se = np.sqrt((w**2 * (ind - p) ** 2).sum()) / W       # delta method, self-normalized IS
        assert abs(p - 1 / S) < 4.5 * se + 1e-12, (h, p, 1 / S, se)


def test_sequential_raw_biased_weighted_unbiased_on_hub_graphs():
    """Four degree-12 hubs among degree-2 vertices: the raw sequential sampler over-produces
    hub-hub edges; the importance-weighted estimate agrees with the exact sampler."""
    deg = [12] * 4 + [2] * 44
    hub = lambda e: int(((e[:, 0] < 4) & (e[:, 1] < 4)).sum())
    exact = np.array([hub(hs.config_degrees(deg, 2, seed=s)) for s in range(6000)])
    f, lw = [], []
    for s in range(20000):
        e, l = hs.sequential(deg, 2, seed=10**6 + s)
        f.append(hub(e)); lw.append(l)
    f, lw = np.array(f, float), np.array(lw)
    ok = np.isfinite(lw)
    w = np.zeros_like(lw); w[ok] = np.exp(lw[ok] - lw[ok].max())
    mu, se_mu = exact.mean(), exact.std() / np.sqrt(len(exact))
    raw, se_raw = f[ok].mean(), f[ok].std() / np.sqrt(ok.sum())
    est = (w * f).sum() / w.sum()
    se = np.sqrt((w**2 * (f - est) ** 2).sum()) / w.sum()
    assert abs(raw - mu) > 8 * np.hypot(se_raw, se_mu)          # raw output is biased
    assert abs(est - mu) < 4.5 * np.hypot(se, se_mu)            # weights remove the bias


def test_switch_chain_mixes_to_uniform_on_graphs():
    deg = [3, 3, 2, 2, 2, 2]
    support = enumerate_simple(6, 2, deg)
    start = np.array(sorted(next(iter(support))), dtype=np.int32)
    counts = Counter()
    for s in range(400 * len(support)):
        e, _ = hs.switch_chain(start, 6, 200, seed=s)
        counts[key(e)] += 1
    assert set(counts) == set(support)
    assert chi2_uniform(counts, support) > 1e-4


def test_switch_chain_preserves_degrees_and_simplicity():
    H = hs.config_regular(300, 3, 6, seed=3)
    E, acc = hs.switch_chain(H, 300, 50000, seed=4)
    assert acc > 0 and hs.is_simple(E, 300) and (hs.degrees(E, 300) == 6).all()
    assert not (E == H).all()


@pytest.mark.parametrize("n,k,p", [(6, 2, 0.3), (7, 3, 0.2), (12, 2, 1.0), (10, 3, 0.0)])
def test_gnp_marginals_and_counts(n, k, p):
    N = hs.binomial(n, k)
    R = 20000 if 0 < p < 1 else 50
    samples = [hs.gnp(n, p, k, seed=s) for s in range(R)]
    sizes = np.array([len(x) for x in samples])
    if p in (0.0, 1.0):
        assert (sizes == round(p * N)).all()
        return
    # edge count ~ Binomial(N, p)
    obs = np.bincount(sizes, minlength=N + 1)
    exp = R * binom.pmf(np.arange(N + 1), N, p)
    keep = exp > 5
    stat = ((obs[keep] - exp[keep]) ** 2 / exp[keep]).sum() + (obs[~keep].sum() - exp[~keep].sum()) ** 2 / max(exp[~keep].sum(), 1e-9)
    assert chi2.sf(stat, keep.sum()) > 1e-4
    # every k-set equally likely, and pairs of k-sets independent
    inc = Counter(tuple(r) for x in samples for r in x)
    freq = np.array([inc[c] for c in itertools.combinations(range(n), k)]) / R
    assert np.abs(freq - p).max() < 5 * np.sqrt(p * (1 - p) / R)
    a, b = tuple(range(k)), tuple(range(1, k + 1))
    both = sum(1 for x in samples if a in map(tuple, x) and b in map(tuple, x)) / R
    assert abs(both - p * p) < 5 * np.sqrt(p * p * (1 - p * p) / R)


def test_gnp_reproducible_and_buffer_retry():
    a = hs.gnp(200, 0.05, 2, seed=11)
    b = hs.gnp(200, 0.05, 2, seed=11)
    assert (a == b).all() and hs.is_simple(a, 200)


@pytest.mark.parametrize("n,k,m", [(5, 2, 3), (5, 2, 8), (6, 3, 4), (6, 3, 18)])
def test_gnm_uniform_over_subsets(n, k, m):
    ksets = list(itertools.combinations(range(n), k))
    support = [frozenset(c) for c in itertools.combinations(ksets, m)]
    R = 200 * len(support) if len(support) < 400 else 40 * len(support)
    counts = Counter(key(hs.gnm(n, m, k, seed=s)) for s in range(R))
    assert chi2_uniform(counts, support) > 1e-4


def test_is_simple_detects_violations():
    assert hs.is_simple(np.array([[0, 1], [1, 2]]), 3)
    assert not hs.is_simple(np.array([[0, 1], [1, 0]]), 3)       # repeated edge (unsorted)
    assert not hs.is_simple(np.array([[0, 0, 1]]), 3)            # repeated vertex
    with pytest.raises(hs.HypersampleError):
        hs.is_simple(np.array([[0, 5]]), 3)                       # vertex out of range


def test_binomial():
    assert hs.binomial(64, 32) == 1832624140942590534
    with pytest.raises(OverflowError):
        hs.binomial(70, 35)
