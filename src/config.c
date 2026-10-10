/* Configuration model conditioned on simplicity, by rejection with early abort.
 *
 * A configuration is a uniformly random partition of the S = sum(deg) labelled stubs into
 * m = S/k blocks of size k. Each simple hypergraph H with this degree sequence arises from
 * exactly prod_v deg(v)! configurations (relabel the stubs at each vertex), so a
 * configuration conditioned on being simple is uniform over simple hypergraphs.
 *
 * Blocks are drawn one at a time by a partial Fisher-Yates shuffle of the stub array, and a
 * trial is abandoned at the first block that repeats a vertex or duplicates an earlier block:
 * the configuration is already non-simple, so abandoning it does not change the accepted
 * distribution. Fisher-Yates from ANY arrangement of the stub multiset is uniform, so the
 * stub array is not reset between trials.
 */
#include "hs_internal.h"

int hs_config_degrees(int32_t n, int32_t k, const int32_t *deg, uint64_t seed,
                      int64_t max_trials, int32_t *edges, int64_t *trials_out) {
    if (trials_out) *trials_out = 0;
    if (n < 1 || k < 1 || k > HS_MAX_K || !deg || !edges || max_trials < 1) return HS_EINVAL;
    int64_t S = 0;
    for (int32_t v = 0; v < n; v++) { if (deg[v] < 0) return HS_EINVAL; S += deg[v]; }
    if (S % k) return HS_EINVAL;
    if (S > 0x7fffffffLL) return HS_EOVERFLOW;
    int64_t m = S / k;
    if (m == 0) return HS_OK;
    int32_t *st = (int32_t *)malloc(sizeof(int32_t) * (size_t)S);
    hs_kset set;
    if (!st) return HS_ENOMEM;
    if (!hs_kset_init(&set, m, edges, k)) { free(st); return HS_ENOMEM; }
    for (int32_t v = 0, p = 0; v < n; v++) for (int32_t j = 0; j < deg[v]; j++) st[p++] = v;
    hs_rng r; hs_rng_seed(&r, seed);
    int status = HS_EXHAUSTED;
    int64_t trials = 0;
    while (trials < max_trials) {
        trials++;
        int64_t left = S; int ok = 1;
        for (int64_t t = 0; t < m && ok; t++) {
            int32_t *e = edges + t * k;
            for (int j = 0; j < k; j++) {                       /* draw k stubs without replacement */
                uint32_t q = hs_below32(&r, (uint32_t)left);
                int32_t v = st[q]; left--; st[q] = st[left]; st[left] = v; e[j] = v;
            }
            if (!hs_distinct_small(e, k)) { ok = 0; break; }
            hs_sort_small(e, k);
            int found; int64_t at = hs_kset_find(&set, e, &found);
            if (found) { ok = 0; break; }
            hs_kset_put(&set, at, t);
        }
        hs_kset_clear(&set);
        if (ok) { status = HS_OK; break; }
    }
    if (trials_out) *trials_out = trials;
    free(st); hs_kset_free(&set);
    return status;
}

int hs_config_regular(int32_t n, int32_t k, int32_t d, uint64_t seed,
                      int64_t max_trials, int32_t *edges, int64_t *trials_out) {
    if (n < 1 || d < 0 || ((int64_t)n * d) % (k > 0 ? k : 1)) { if (trials_out) *trials_out = 0; return HS_EINVAL; }
    int32_t *deg = (int32_t *)malloc(sizeof(int32_t) * (size_t)n);
    if (!deg) return HS_ENOMEM;
    for (int32_t v = 0; v < n; v++) deg[v] = d;
    int s = hs_config_degrees(n, k, deg, seed, max_trials, edges, trials_out);
    free(deg);
    return s;
}
