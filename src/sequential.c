/* Sequential stub matching with exact importance weights.
 *
 * State: r[v] = unused stubs at vertex v, E = edges chosen so far. A valid stub k-set is k
 * unused stubs at k distinct vertices whose vertex set is not already in E. Their number is
 *
 *     V = e_k(r) - sum_{e in E} prod_{v in e} r[v],
 *
 * e_k the k-th elementary symmetric polynomial of r. Each step picks one valid stub k-set
 * uniformly (rejection from uniform stub k-sets, which is exact), so a run has probability
 * prod_t 1/V_t. Each simple hypergraph is produced by the same number of runs
 * (m! * prod_v deg(v)!), hence w = prod_t V_t is an importance weight for the uniform law.
 *
 * e_k is computed from the power sums p_j = sum_v r[v]^j by Newton's identities; p_j and the
 * correction term are updated in O(k + deg) per stub, so a run costs O(m k (k + d_max)) plus
 * the rejection draws. Exact integer arithmetic (__int128) needs sum(deg)^k < 2^126: k <= 4
 * and sum(deg) < 2^31.
 */
#include <math.h>
#include "hs_internal.h"

typedef __int128 i128;

typedef struct {
    int k;
    i128 p[5];            /* power sums p_1..p_k */
    i128 corr;            /* sum over chosen edges of prod r */
} seq_state;

static i128 ipow(i128 x, int j) { i128 y = 1; while (j--) y *= x; return y; }

static i128 elem_sym(const seq_state *s) {
    i128 e[5]; e[0] = 1;
    for (int j = 1; j <= s->k; j++) {
        i128 acc = 0;
        for (int i = 1; i <= j; i++) acc += ((i & 1) ? 1 : -1) * e[j - i] * s->p[i];
        e[j] = acc / j;                                   /* exact: Newton's identities over Z */
    }
    return e[s->k];
}

int hs_sequential_degrees(int32_t n, int32_t k, const int32_t *deg, uint64_t seed,
                          int32_t *edges, double *log_weight) {
    if (log_weight) *log_weight = 0.0;
    if (n < 1 || k < 2 || k > 4 || !deg || !edges) return HS_EINVAL;
    int64_t S = 0; int32_t dmax = 0;
    for (int32_t v = 0; v < n; v++) { if (deg[v] < 0) return HS_EINVAL; S += deg[v]; if (deg[v] > dmax) dmax = deg[v]; }
    if (S % k) return HS_EINVAL;
    if (S > 0x7fffffffLL) return HS_EOVERFLOW;
    int64_t m = S / k;
    if (m == 0) return HS_OK;

    int32_t *r = (int32_t *)malloc(sizeof(int32_t) * (size_t)n);
    int32_t *st = (int32_t *)malloc(sizeof(int32_t) * (size_t)S);
    int64_t *inc = (int64_t *)malloc(sizeof(int64_t) * (size_t)S);    /* incident edge lists */
    int64_t *off = (int64_t *)malloc(sizeof(int64_t) * (size_t)(n + 1));
    int32_t *cnt = (int32_t *)calloc((size_t)n, sizeof(int32_t));
    hs_kset set;
    int status = HS_OK;
    if (!r || !st || !inc || !off || !cnt) { status = HS_ENOMEM; goto done_nofree_set; }
    if (!hs_kset_init(&set, m, edges, k)) { status = HS_ENOMEM; goto done_nofree_set; }

    off[0] = 0;
    for (int32_t v = 0; v < n; v++) { r[v] = deg[v]; off[v + 1] = off[v] + deg[v]; }
    for (int32_t v = 0, q = 0; v < n; v++) for (int32_t j = 0; j < deg[v]; j++) st[q++] = v;

    seq_state s; s.k = k; s.corr = 0;
    for (int j = 1; j <= k; j++) { s.p[j] = 0; for (int32_t v = 0; v < n; v++) s.p[j] += ipow(deg[v], j); }

    hs_rng g; hs_rng_seed(&g, seed);
    double lw = 0.0;
    int64_t left = S;
    int32_t cand[HS_MAX_K];
    for (int64_t t = 0; t < m; t++) {
        i128 V = elem_sym(&s) - s.corr;
        if (V <= 0) { status = HS_ESTUCK; lw = -INFINITY; break; }
        lw += log((double)V);
        for (;;) {                                         /* uniform valid stub k-set by rejection */
            int64_t L = left;
            for (int j = 0; j < k; j++) {
                uint32_t q = hs_below32(&g, (uint32_t)L);
                int32_t v = st[q]; L--; st[q] = st[L]; st[L] = v; cand[j] = v;
            }
            if (!hs_distinct_small(cand, k)) continue;     /* the stubs stay in the pool (only permuted) */
            hs_sort_small(cand, k);
            memcpy(edges + t * k, cand, sizeof(int32_t) * (size_t)k);
            int found; int64_t at = hs_kset_find(&set, edges + t * k, &found);
            if (found) continue;
            hs_kset_put(&set, at, t);
            break;
        }
        left -= k;                                         /* accepted stubs are st[left-k .. left-1] */
        for (int j = 0; j < k; j++) {                      /* r[v] -= 1, update p_j and corr */
            int32_t v = edges[t * k + j];
            i128 a = r[v], b = a - 1;
            for (int i = 1; i <= k; i++) s.p[i] -= ipow(a, i) - ipow(b, i);
            for (int32_t q = 0; q < cnt[v]; q++) {        /* existing edges at v lose one factor */
                const int32_t *e = edges + inc[off[v] + q] * k;
                i128 prod = 1;
                for (int i = 0; i < k; i++) if (e[i] != v) prod *= r[e[i]];
                s.corr -= prod;
            }
            r[v]--;
        }
        i128 prod = 1;
        for (int j = 0; j < k; j++) prod *= r[edges[t * k + j]];
        s.corr += prod;
        for (int j = 0; j < k; j++) { int32_t v = edges[t * k + j]; inc[off[v] + cnt[v]++] = t; }
    }
    if (log_weight) *log_weight = lw;
    hs_kset_free(&set);
done_nofree_set:
    free(r); free(st); free(inc); free(off); free(cnt);
    return status;
}
