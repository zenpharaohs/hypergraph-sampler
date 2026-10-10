/* Version, status strings and checks. */
#include "hs_internal.h"

#define HS_VERSION_STRING "0.1.0"

int hs_abi_version(void) { return HS_ABI_VERSION; }
const char *hs_version_string(void) { return HS_VERSION_STRING; }

const char *hs_status_string(int status) {
    switch (status) {
    case HS_OK: return "ok";
    case HS_EXHAUSTED: return "no simple configuration within max_trials (too close to saturation for rejection?)";
    case HS_EINVAL: return "invalid arguments";
    case HS_ENOMEM: return "out of memory";
    case HS_ESTUCK: return "sequential sampler stuck: no valid edge remains (weight 0)";
    case HS_EBUFFER: return "output buffer too small; retry with the reported edge count";
    case HS_EOVERFLOW: return "a required count does not fit in 64 bits";
    default: return "unknown status";
    }
}

int hs_degrees(int32_t n, int32_t k, int64_t m, const int32_t *edges, int32_t *deg_out) {
    if (n < 1 || k < 1 || m < 0 || (m > 0 && !edges) || !deg_out) return HS_EINVAL;
    memset(deg_out, 0, sizeof(int32_t) * (size_t)n);
    for (int64_t i = 0; i < m * k; i++) {
        if (edges[i] < 0 || edges[i] >= n) return HS_EINVAL;
        deg_out[edges[i]]++;
    }
    return HS_OK;
}

int hs_is_simple(int32_t n, int32_t k, int64_t m, const int32_t *edges) {
    if (n < 1 || k < 1 || k > HS_MAX_K || m < 0 || (m > 0 && !edges)) return -HS_EINVAL;
    if (m == 0) return 1;
    int32_t *copy = (int32_t *)malloc(sizeof(int32_t) * (size_t)(m * k));
    hs_kset set;
    if (!copy) return -HS_ENOMEM;
    memcpy(copy, edges, sizeof(int32_t) * (size_t)(m * k));
    if (!hs_kset_init(&set, m, copy, k)) { free(copy); return -HS_ENOMEM; }
    int simple = 1;
    for (int64_t t = 0; t < m && simple; t++) {
        int32_t *e = copy + t * k;
        for (int j = 0; j < k; j++) if (e[j] < 0 || e[j] >= n) { simple = -HS_EINVAL; break; }
        if (simple != 1) break;
        if (!hs_distinct_small(e, k)) { simple = 0; break; }
        hs_sort_small(e, k);
        int f; int64_t at = hs_kset_find(&set, e, &f);
        if (f) { simple = 0; break; }
        hs_kset_put(&set, at, t);
    }
    hs_kset_free(&set); free(copy);
    return simple;
}
