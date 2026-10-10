/* Erdos-Renyi samplers on k-sets: G(n, p) by geometric skipping, G(n, m) uniform.
 *
 * k-subsets of {0..n-1} are indexed 0..C(n,k)-1 in colex order by the combinatorial number
 * system: the set c_1 < ... < c_k has index sum_j C(c_j, j).
 */
#include <math.h>
#include "hs_internal.h"

/* C(n, k) saturated at UINT64_MAX. */
static uint64_t binom_sat(int64_t n, int64_t k) {
    if (k < 0 || n < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    unsigned __int128 c = 1;
    for (int64_t i = 1; i <= k; i++) {
        c = c * (unsigned __int128)(n - k + i) / (unsigned __int128)i;   /* exact at every step */
        if (c > (unsigned __int128)UINT64_MAX) return UINT64_MAX;
    }
    return (uint64_t)c;
}

int64_t hs_binomial(int32_t n, int32_t k) {
    uint64_t c = binom_sat(n, k);
    return c >= ((uint64_t)1 << 63) ? -1 : (int64_t)c;
}

/* index -> sorted k-set (colex). */
static void unrank(uint64_t idx, int32_t n, int32_t k, int32_t *out) {
    int32_t hi = n - 1;
    for (int32_t j = k; j >= 1; j--) {
        int32_t lo = j - 1, h = hi;                       /* largest c in [j-1, hi] with C(c, j) <= idx */
        while (lo < h) {
            int32_t mid = lo + (h - lo + 1) / 2;
            if (binom_sat(mid, j) <= idx) lo = mid; else h = mid - 1;
        }
        out[j - 1] = lo;
        idx -= binom_sat(lo, j);
        hi = lo - 1;
    }
}

int hs_gnp(int32_t n, int32_t k, double p, uint64_t seed, int32_t *edges,
           int64_t capacity, int64_t *m_out) {
    if (m_out) *m_out = 0;
    if (n < 1 || k < 1 || k > HS_MAX_K || k > n || !(p >= 0.0 && p <= 1.0) || capacity < 0 || (capacity > 0 && !edges))
        return HS_EINVAL;
    int64_t N = hs_binomial(n, k);
    if (N < 0) return HS_EOVERFLOW;
    hs_rng g; hs_rng_seed(&g, seed);
    int64_t count = 0;
    int32_t tmp[HS_MAX_K];
    if (p > 0.0) {
        double lq = log1p(-p);                             /* -inf when p = 1: every skip is 0 */
        int64_t i = -1;
        for (;;) {
            double skip = (p >= 1.0) ? 0.0 : floor(log(hs_unit_open0(&g)) / lq);
            if (skip >= (double)(N - 1 - i)) break;       /* next index would be >= N */
            i += 1 + (int64_t)skip;
            int32_t *dst = (count < capacity) ? edges + count * k : tmp;
            unrank((uint64_t)i, n, k, dst);
            count++;
        }
    }
    if (m_out) *m_out = count;
    return count > capacity ? HS_EBUFFER : HS_OK;
}

/* Uniform random k-subset of {0..n-1} by Floyd's algorithm, sorted. */
static void random_kset(hs_rng *g, int32_t n, int32_t k, int32_t *out) {
    int c = 0;
    for (int32_t j = n - k; j < n; j++) {
        int32_t t = (int32_t)hs_below32(g, (uint32_t)(j + 1));
        int dup = 0;
        for (int i = 0; i < c; i++) if (out[i] == t) { dup = 1; break; }
        out[c++] = dup ? j : t;
    }
    hs_sort_small(out, k);
}

int hs_gnm(int32_t n, int32_t k, int64_t m, uint64_t seed, int32_t *edges) {
    if (n < 1 || k < 1 || k > HS_MAX_K || k > n || m < 0 || (m > 0 && !edges)) return HS_EINVAL;
    int64_t N = hs_binomial(n, k);
    if (N < 0) return HS_EOVERFLOW;
    if (m > N) return HS_EINVAL;
    if (m == 0) return HS_OK;
    hs_rng g; hs_rng_seed(&g, seed);
    if (2 * m <= N) {                                      /* rejection of repeats: < 2 draws per edge on average */
        hs_kset set;
        if (!hs_kset_init(&set, m, edges, k)) return HS_ENOMEM;
        for (int64_t t = 0; t < m; ) {
            random_kset(&g, n, k, edges + t * k);
            int f; int64_t at = hs_kset_find(&set, edges + t * k, &f);
            if (f) continue;
            hs_kset_put(&set, at, t); t++;
        }
        hs_kset_free(&set);
        return HS_OK;
    }
    /* dense case: choose the N - m k-sets to leave out, then list the rest in colex order */
    int64_t out_n = N - m;
    int32_t *gone = (int32_t *)malloc(sizeof(int32_t) * (size_t)((out_n > 0 ? out_n : 1) * k));
    hs_kset set;
    if (!gone) return HS_ENOMEM;
    if (!hs_kset_init(&set, out_n, gone, k)) { free(gone); return HS_ENOMEM; }
    for (int64_t t = 0; t < out_n; ) {
        random_kset(&g, n, k, gone + t * k);
        int f; int64_t at = hs_kset_find(&set, gone + t * k, &f);
        if (f) continue;
        hs_kset_put(&set, at, t); t++;
    }
    int64_t w = 0;
    for (int64_t i = 0; i < N; i++) {
        int32_t *e = edges + w * k;
        unrank((uint64_t)i, n, k, e);
        int f; hs_kset_find(&set, e, &f);
        if (!f) w++;
    }
    hs_kset_free(&set); free(gone);
    return HS_OK;
}
