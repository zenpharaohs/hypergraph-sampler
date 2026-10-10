/* Degree-preserving switch chain on simple k-uniform hypergraphs.
 *
 * Proposal: an ordered pair of distinct edge indices (a, b) uniformly, then u uniform in
 * edge a and v uniform in edge b; the move swaps u and v between the two edges. It is
 * proposed with probability 1/(m(m-1)k^2) and its reverse (same pair, v in a', u in b')
 * with the same probability, so accepting exactly the moves that keep the hypergraph simple
 * (Metropolis with a symmetric proposal) leaves the uniform distribution invariant.
 * Rejected: v already in a, u already in b (repeated vertex), or a new edge that already
 * exists (repeated edge). For k = 2 this is the double-edge swap chain.
 */
#include "hs_internal.h"

static int pos_in(const int32_t *e, int k, int32_t x) {
    for (int i = 0; i < k; i++) if (e[i] == x) return i;
    return -1;
}

int hs_switch_chain(int32_t n, int32_t k, int64_t m, int32_t *edges, int64_t steps,
                    uint64_t seed, int64_t *accepted) {
    if (accepted) *accepted = 0;
    if (n < 1 || k < 1 || k > HS_MAX_K || m < 0 || !edges || steps < 0) return HS_EINVAL;
    if (m > 0x7fffffffLL) return HS_EOVERFLOW;
    if (hs_is_simple(n, k, m, edges) != 1) return HS_EINVAL;
    if (m < 2 || steps == 0) return HS_OK;
    for (int64_t t = 0; t < m; t++) hs_sort_small(edges + t * k, k);
    hs_kset set;
    if (!hs_kset_init(&set, m, edges, k)) return HS_ENOMEM;
    for (int64_t t = 0; t < m; t++) { int f; int64_t at = hs_kset_find(&set, edges + t * k, &f); hs_kset_put(&set, at, t); }
    hs_rng g; hs_rng_seed(&g, seed);
    int32_t na[HS_MAX_K], nb[HS_MAX_K];
    int64_t acc = 0;
    for (int64_t s = 0; s < steps; s++) {
        int64_t a = hs_below32(&g, (uint32_t)m), b = hs_below32(&g, (uint32_t)(m - 1));
        if (b >= a) b++;
        int32_t *ea = edges + a * k, *eb = edges + b * k;
        int iu = (int)hs_below32(&g, (uint32_t)k), iv = (int)hs_below32(&g, (uint32_t)k);
        int32_t u = ea[iu], v = eb[iv];
        if (u == v || pos_in(ea, k, v) >= 0 || pos_in(eb, k, u) >= 0) continue;
        memcpy(na, ea, sizeof(int32_t) * (size_t)k); na[iu] = v; hs_sort_small(na, k);
        memcpy(nb, eb, sizeof(int32_t) * (size_t)k); nb[iv] = u; hs_sort_small(nb, k);
        int f1, f2;
        hs_kset_find(&set, na, &f1);
        if (f1) continue;
        hs_kset_find(&set, nb, &f2);
        if (f2) continue;
        if (memcmp(na, nb, sizeof(int32_t) * (size_t)k) == 0) continue;   /* cannot happen for u != v, kept as a guard */
        /* remove a and b from the set, rewrite them, reinsert */
        int f; int64_t at;
        at = hs_kset_find(&set, ea, &f); hs_kset_remove_at(&set, at);
        at = hs_kset_find(&set, eb, &f); hs_kset_remove_at(&set, at);
        memcpy(ea, na, sizeof(int32_t) * (size_t)k);
        memcpy(eb, nb, sizeof(int32_t) * (size_t)k);
        at = hs_kset_find(&set, ea, &f); set.slot[at] = a;
        at = hs_kset_find(&set, eb, &f); set.slot[at] = b;
        acc++;
    }
    if (accepted) *accepted = acc;
    hs_kset_free(&set);
    return HS_OK;
}
