# hypergraph-sampler

Random graphs and uniform hypergraphs: exact and approximate samplers in a small C core with a
stable C ABI, and thin bindings for Python and MATLAB (or GNU Octave).

The C library (`hypersample`) holds all the algorithms. The bindings only convert arrays and
call it, so a new language (R, Go, Julia, ...) needs only a binding, and the core can be
retargeted to other hardware independently of who calls it.

## What it samples

A k-uniform hypergraph on vertices 0..n−1 (1..n in MATLAB) is an m × k integer array, one
sorted edge per row; k = 2 gives ordinary graphs. *Simple* means no edge repeats a vertex and
no edge appears twice.

| function | distribution | exact? |
|---|---|---|
| `config_degrees`, `config_regular` | simple k-uniform hypergraph with a given degree sequence (e.g. d-regular) | exactly uniform |
| `gnp` | Erdős–Rényi G(n, p) on k-sets: each k-subset independently with probability p | exact |
| `gnm` | G(n, m) on k-sets: m distinct k-subsets, uniform | exactly uniform |
| `sequential` | sequential stub matching with an exact importance weight | biased, weight corrects it |
| `switch_chain` | degree-preserving switch Markov chain | uniform is stationary |

Every sampler takes an explicit 64-bit seed and is reproducible in it.

## Quick start

Python (needs numpy and a C compiler; the library is compiled on first import, or run `make`):

```python
import hypergraph_sampler as hs                       # with PYTHONPATH=python, or pip install -e .
H = hs.config_regular(n=1000, k=3, d=12, seed=1)      # 4000 x 3 array, exactly uniform
G = hs.config_degrees([3, 3, 2, 2, 2, 2], k=2, seed=2)
R = hs.gnp(500, 0.01, seed=3)                         # ordinary random graph
E, log_w = hs.sequential([12] * 4 + [2] * 44, k=2, seed=4)
```

MATLAB or Octave:

```matlab
addpath matlab; build_hypersample            % once: compiles hs_mex from ../src
H = hs_config_regular(1000, 3, 12, 1);       % 4000-by-3 int32, vertices 1..1000
test_hypersample; demo_hypersample           % checks, and a classroom demo
```

## Methods

**Configuration model with rejection (exact).** Give vertex v deg(v) labelled stubs and split
the S = Σ deg stubs uniformly into blocks of k. Every simple hypergraph comes from exactly
Π_v deg(v)! configurations, so a configuration *conditioned on being simple* is uniform. The
sampler draws blocks one at a time by a partial Fisher–Yates shuffle and abandons a trial at
the first block that repeats a vertex or an edge; the configuration is already non-simple, so
the accepted law is unchanged. For d-regular k-uniform hypergraphs the number of blocks with a
repeated vertex is asymptotically Poisson with mean (k−1)(d−1)/2. Repeated edges matter only
for graphs: for k = 2 their number is asymptotically Poisson with mean (d−1)²/4, while for
k ≥ 3 they are O(1/n). Hence

  k = 2:  P(simple) → exp(−(d−1)/2 − (d−1)²/4)   (e⁻² for cubic graphs; Bender–Canfield),
  k ≥ 3:  P(simple) → exp(−(k−1)(d−1)/2)          (e⁻¹¹ for 12-regular 3-graphs).

The expected work is about 1/P(simple) times the length of a failed trial: linear in n
with a constant set by k and d. Close to saturation (small n for the degree) the acceptance
collapses; the sampler then stops after `max_trials` rather than change the distribution.

**Sequential stub matching with importance weights.** Build the edges one at a time, each
uniform among the V_t currently *valid* stub k-sets (k free stubs at distinct vertices whose
vertex set is not yet an edge). A run has probability Π_t 1/V_t and each simple hypergraph has
the same number of runs, so w = Π_t V_t is an exact importance weight for the uniform law:

  Σ_i w_i f(H_i) / Σ_i w_i → E_uniform[f].

V_t = e_k(r) − Σ_{e ∈ E} Π_{v∈e} r_v, where r is the vector of free stubs and e_k the
elementary symmetric polynomial; it is updated in O(k + d_max) per stub from power sums by
Newton's identities, in exact 128-bit integer arithmetic (k ≤ 4). A run can get stuck
(V_t = 0); it then has weight 0 and must stay in the estimator. The raw output is biased —
`tests/test_samplers.py` shows a case where it over-produces hub-to-hub edges by more than 40 standard
errors and the weighted estimate agrees with the exact sampler.

**Switch chain.** Pick an ordered pair of edges and one vertex of each uniformly; swap the two
vertices; accept if both edges stay simple and new. The proposal is symmetric, so this
Metropolis chain has the uniform law on simple hypergraphs with the given degrees as a
stationary distribution. For graphs it is the double-edge-swap chain, which is irreducible;
for k > 2 irreducibility is not guaranteed in general.

**G(n, p) by geometric skipping.** The gaps between successive chosen k-subsets in a fixed order
are geometric, ⌊log U / log(1−p)⌋, so the cost is proportional to the number of edges, not to
C(n, k) (Batagelj–Brandes for graphs). k-subsets are indexed by the combinatorial number system
and unranked by binary search.

**G(n, m).** Uniform k-subsets by Floyd's algorithm, rejecting repeats (fewer than two draws per
edge while m ≤ C(n, k)/2); for denser cases the complement is sampled instead.

## Layout

```
include/hypersample.h   the C ABI (the contract every binding uses)
src/                    C core: config.c, sequential.c, switch.c, erdos_renyi.c, common.c
python/hypergraph_sampler/   ctypes binding
matlab/                 MEX gateway (hs_mex.c), hs_*.m wrappers, build, test and demo scripts
tests/                  smoke.c (C) and test_samplers.py (exactness against full enumeration)
examples/               short Python scripts for teaching
```

## Tests

`make test` builds the library, runs the C smoke test and the Python tests. The exactness tests
enumerate every simple (hyper)graph for small degree sequences and apply chi-square
goodness-of-fit tests to the configuration sampler, G(n, m), the importance-weighted sequential
sampler and the switch chain. In MATLAB or Octave, `test_hypersample` checks the binding.

## License

MIT, see [LICENSE](LICENSE).
