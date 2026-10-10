"""Random graphs and uniform hypergraphs (Python binding of the hypersample C core).

A k-uniform hypergraph on vertices 0..n-1 is an (m, k) int32 array, one sorted row per edge;
k = 2 gives ordinary graphs. Every sampler takes a seed (an int below 2**64; None draws one
from the OS) and is reproducible in it.

    from hypergraph_sampler import config_regular
    H = config_regular(n=1000, k=3, d=12, seed=1)     # exactly uniform simple 12-regular 3-graph
"""
from __future__ import annotations

import ctypes
import os

import numpy as np

from ._lib import load

__all__ = ["config_regular", "config_degrees", "sequential", "switch_chain", "gnp", "gnm",
           "is_simple", "degrees", "binomial", "version", "HypersampleError", "SaturationError"]

_lib = load()


class HypersampleError(RuntimeError):
    def __init__(self, status: int, where: str):
        self.status = status
        super().__init__(f"{where}: {_lib.hs_status_string(status).decode()} (status {status})")


class SaturationError(HypersampleError):
    """Rejection exhausted max_trials (the parameters are too close to saturation)."""


def version() -> str:
    return _lib.hs_version_string().decode()


def _seed(seed):
    if seed is None:
        return int.from_bytes(os.urandom(8), "little")
    seed = int(seed)
    if not 0 <= seed < 2**64:
        raise ValueError("seed must be an integer in [0, 2**64)")
    return seed


def _p32(a):
    return a.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))


def _check(status, where):
    if status == 1:
        raise SaturationError(status, where)
    if status != 0:
        raise HypersampleError(status, where)


def _deg_array(deg):
    deg = np.ascontiguousarray(deg, dtype=np.int32)
    if deg.ndim != 1 or len(deg) == 0:
        raise ValueError("degree sequence must be a nonempty 1-d array")
    return deg


def config_degrees(deg, k: int = 2, seed=None, max_trials: int = 10**9, *, return_trials: bool = False):
    """Exactly uniform simple k-uniform hypergraph with degree sequence `deg`.

    Configuration model conditioned on simplicity, by rejection with early abort. For
    d-regular hypergraphs P(simple) -> exp(-(k-1)(d-1)/2) when k >= 3, and
    exp(-(d-1)/2 - (d-1)**2/4) for graphs (k = 2); the expected number of trials is the
    inverse. Raises SaturationError after max_trials.
    """
    deg = _deg_array(deg)
    if deg.sum() % k:
        raise ValueError("sum of degrees must be divisible by k")
    m = int(deg.sum()) // k
    out = np.empty((m, k), dtype=np.int32)
    trials = ctypes.c_int64(0)
    _check(_lib.hs_config_degrees(len(deg), k, _p32(deg), _seed(seed), int(max_trials), _p32(out),
                                  ctypes.byref(trials)), "config_degrees")
    return (out, trials.value) if return_trials else out


def config_regular(n: int, k: int = 2, d: int = 3, seed=None, max_trials: int = 10**9, *, return_trials: bool = False):
    """Exactly uniform simple d-regular k-uniform hypergraph on n vertices (n*d divisible by k)."""
    if (n * d) % k:
        raise ValueError("n*d must be divisible by k")
    out = np.empty((n * d // k, k), dtype=np.int32)
    trials = ctypes.c_int64(0)
    _check(_lib.hs_config_regular(n, k, d, _seed(seed), int(max_trials), _p32(out), ctypes.byref(trials)),
           "config_regular")
    return (out, trials.value) if return_trials else out


def sequential(deg, k: int = 2, seed=None):
    """Sequential stub matching with an exact importance weight (k <= 4).

    Returns (edges, log_w). Each edge is chosen uniformly among the currently valid stub
    k-sets; w = prod_t V_t makes  sum(w f) / sum(w)  a consistent estimate of the uniform
    average of f. A stuck run returns (partial edges, -inf): keep it with weight 0.
    """
    deg = _deg_array(deg)
    if deg.sum() % k:
        raise ValueError("sum of degrees must be divisible by k")
    m = int(deg.sum()) // k
    out = np.zeros((m, k), dtype=np.int32)
    lw = ctypes.c_double(0.0)
    s = _lib.hs_sequential_degrees(len(deg), k, _p32(deg), _seed(seed), _p32(out), ctypes.byref(lw))
    if s not in (0, 4):
        _check(s, "sequential")
    return out, lw.value


def switch_chain(edges, n: int, steps: int, seed=None):
    """Run `steps` steps of the degree-preserving switch chain from `edges` (not modified).

    Returns (edges, accepted). Uniform is stationary; irreducible for graphs (k = 2).
    """
    e = np.array(edges, dtype=np.int32, order="C", copy=True)
    if e.ndim != 2:
        raise ValueError("edges must be an (m, k) array")
    acc = ctypes.c_int64(0)
    _check(_lib.hs_switch_chain(n, e.shape[1], e.shape[0], _p32(e), int(steps), _seed(seed), ctypes.byref(acc)),
           "switch_chain")
    return e, acc.value


def gnp(n: int, p: float, k: int = 2, seed=None):
    """Erdos-Renyi G(n, p) on k-sets: each k-subset independently with probability p."""
    seed = _seed(seed)
    cap = int(min(binomial(n, k), max(16, 1.2 * p * binomial(n, k) + 10 * (p * binomial(n, k)) ** 0.5 + 16)))
    m = ctypes.c_int64(0)
    out = np.empty((cap, k), dtype=np.int32)
    s = _lib.hs_gnp(n, k, float(p), seed, _p32(out), cap, ctypes.byref(m))
    if s == 5:                                           # buffer too small: same seed, exact size
        out = np.empty((m.value, k), dtype=np.int32)
        s = _lib.hs_gnp(n, k, float(p), seed, _p32(out), m.value, ctypes.byref(m))
    _check(s, "gnp")
    return out[: m.value].copy()


def gnm(n: int, m: int, k: int = 2, seed=None):
    """Uniform G(n, m) on k-sets: m distinct k-subsets, uniformly over all such choices."""
    out = np.empty((m, k), dtype=np.int32)
    _check(_lib.hs_gnm(n, k, int(m), _seed(seed), _p32(out)), "gnm")
    return out


def is_simple(edges, n: int) -> bool:
    e = np.ascontiguousarray(edges, dtype=np.int32)
    r = _lib.hs_is_simple(n, e.shape[1], e.shape[0], _p32(e))
    if r < 0:
        _check(-r, "is_simple")
    return bool(r)


def degrees(edges, n: int):
    e = np.ascontiguousarray(edges, dtype=np.int32)
    out = np.empty(n, dtype=np.int32)
    _check(_lib.hs_degrees(n, e.shape[1], e.shape[0], _p32(e), _p32(out)), "degrees")
    return out


def binomial(n: int, k: int) -> int:
    c = _lib.hs_binomial(n, k)
    if c < 0:
        raise OverflowError("C(n, k) >= 2**63")
    return c
