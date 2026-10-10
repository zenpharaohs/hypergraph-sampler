/* hypersample: random graphs and uniform hypergraphs, C core with a stable C ABI.
 *
 * Conventions shared by every entry point (bindings rely on them):
 *   - A k-uniform hypergraph on vertices 0..n-1 with m edges is an int32 array of
 *     m*k entries, row-major: edge t is edges[t*k .. t*k+k-1], sorted ascending.
 *     k = 2 is an ordinary graph.
 *   - "Simple" means: the k vertices of an edge are distinct, and no edge repeats.
 *   - Every random routine takes an explicit 64-bit seed and is deterministic in it
 *     (same seed, same arguments, same platform word size -> same output).
 *   - The caller owns all memory; outputs are written into caller-provided buffers.
 *   - Return value is an hs_status code (0 = success).
 *
 * ABI stability: functions are only added, never changed, within one
 * HS_ABI_VERSION. A binding should check hs_abi_version() at load time.
 */
#ifndef HYPERSAMPLE_H
#define HYPERSAMPLE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) && defined(HS_BUILD_DLL)
#define HS_API __declspec(dllexport)
#else
#define HS_API
#endif

#define HS_ABI_VERSION 1
#define HS_MAX_K 16

enum hs_status {
    HS_OK = 0,
    HS_EXHAUSTED = 1,   /* rejection sampler used max_trials without a simple configuration */
    HS_EINVAL = 2,      /* invalid arguments (see each function) */
    HS_ENOMEM = 3,      /* allocation failed */
    HS_ESTUCK = 4,      /* sequential sampler reached a state with no valid edge */
    HS_EBUFFER = 5,     /* output buffer too small; *m_out holds the required edge count */
    HS_EOVERFLOW = 6    /* a count needed by the method does not fit in 64 bits */
};

HS_API int hs_abi_version(void);
HS_API const char *hs_version_string(void);
HS_API const char *hs_status_string(int status);

/* ---------------------------------------------------------------- exact samplers */

/* Exactly uniform simple k-uniform hypergraph with the given degree sequence
 * (configuration model conditioned on simplicity, by rejection with early abort).
 * deg[0..n-1] >= 0, sum(deg) divisible by k; m = sum(deg)/k edges are written.
 * Expected trials ~ 1/P(simple); for d-regular, P(simple) -> exp(-(k-1)(d-1)/2) when k >= 3
 * and exp(-(d-1)/2 - (d-1)^2/4) for graphs (k = 2), as n -> infinity.
 * *trials_out (may be NULL) receives the number of configurations tried. */
HS_API int hs_config_degrees(int32_t n, int32_t k, const int32_t *deg, uint64_t seed,
                             int64_t max_trials, int32_t *edges, int64_t *trials_out);

/* Same, every degree equal to d (m = n*d/k). */
HS_API int hs_config_regular(int32_t n, int32_t k, int32_t d, uint64_t seed,
                             int64_t max_trials, int32_t *edges, int64_t *trials_out);

/* Erdos-Renyi G(n, p) generalized to k-sets: each of the C(n,k) k-subsets is an edge
 * independently with probability p. Geometric skipping over the k-subset index space
 * (Batagelj-Brandes for k = 2), unranked by the combinatorial number system; expected
 * time O(k^2 log n * (1 + p C(n,k))). Edges are written in increasing colex order.
 * If more than `capacity` edges occur, writing stops, counting continues, *m_out is the
 * full count and HS_EBUFFER is returned: call again with the same seed and a buffer of
 * that size to get the identical sample. */
HS_API int hs_gnp(int32_t n, int32_t k, double p, uint64_t seed, int32_t *edges,
                  int64_t capacity, int64_t *m_out);

/* Uniform G(n, m) on k-sets: m distinct k-subsets, uniform over all C(C(n,k), m) choices.
 * Requires m <= C(n,k). Edges in arbitrary order (each sorted). */
HS_API int hs_gnm(int32_t n, int32_t k, int64_t m, uint64_t seed, int32_t *edges);

/* ---------------------------------------------------------------- approximate samplers */

/* Sequential stub matching with importance weights (k <= 4, sum(deg) < 2^31).
 * Edge t is chosen uniformly among the V_t valid stub k-sets of the remaining stubs
 * (distinct vertices, not an existing edge). A run that produces hypergraph H has
 * probability prod_t 1/V_t, and every simple H has the same number of runs, so
 *     w = prod_t V_t   (written as log w to *log_weight)
 * is an importance weight for the uniform distribution: E[w f] / E[w] = E_unif[f].
 * If some V_t = 0 the run is stuck: HS_ESTUCK, *log_weight = -INFINITY (weight 0),
 * edges partially written. Stuck runs must be kept in the estimator with weight 0. */
HS_API int hs_sequential_degrees(int32_t n, int32_t k, const int32_t *deg, uint64_t seed,
                                 int32_t *edges, double *log_weight);

/* Degree-preserving switch chain (Metropolis, symmetric proposal), in place.
 * One step: choose an ordered pair of distinct edges (e, f), u in e and v in f uniformly;
 * propose e' = e - u + v, f' = f - v + u; accept if both are simple and new.
 * The uniform distribution on simple hypergraphs with this degree sequence is stationary;
 * for graphs (k = 2) the chain is irreducible (connected by double-edge swaps); for k > 2
 * irreducibility is not guaranteed in general. Edges stay sorted. *accepted may be NULL. */
HS_API int hs_switch_chain(int32_t n, int32_t k, int64_t m, int32_t *edges, int64_t steps,
                           uint64_t seed, int64_t *accepted);

/* ---------------------------------------------------------------- utilities */

/* 1 if simple, 0 if not, negative -HS_EINVAL on bad arguments. Edges need not be sorted. */
HS_API int hs_is_simple(int32_t n, int32_t k, int64_t m, const int32_t *edges);

/* deg_out[0..n-1] = vertex degrees. */
HS_API int hs_degrees(int32_t n, int32_t k, int64_t m, const int32_t *edges, int32_t *deg_out);

/* C(n, k) if it is below 2^63, else -1. */
HS_API int64_t hs_binomial(int32_t n, int32_t k);

#ifdef __cplusplus
}
#endif
#endif
