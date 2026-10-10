/* C-level smoke test: every entry point runs and returns well-formed output. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "hypersample.h"

#define CHECK(c) do { if (!(c)) { fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); return 1; } } while (0)

int main(void) {
    CHECK(hs_abi_version() == HS_ABI_VERSION);
    int32_t n = 200, k = 3, d = 12, m = n * d / k;
    int32_t *e = malloc(sizeof(int32_t) * (size_t)m * (size_t)k), deg[200];
    int64_t trials;
    CHECK(hs_config_regular(n, k, d, 42, 100000000, e, &trials) == HS_OK);
    CHECK(hs_is_simple(n, k, m, e) == 1);
    CHECK(hs_degrees(n, k, m, e, deg) == HS_OK);
    for (int v = 0; v < n; v++) CHECK(deg[v] == d);
    int64_t acc;
    CHECK(hs_switch_chain(n, k, m, e, 10000, 7, &acc) == HS_OK && acc > 0);
    CHECK(hs_is_simple(n, k, m, e) == 1);
    CHECK(hs_degrees(n, k, m, e, deg) == HS_OK);
    for (int v = 0; v < n; v++) CHECK(deg[v] == d);
    double lw;
    int s = hs_sequential_degrees(n, k, deg, 3, e, &lw);
    CHECK(s == HS_OK || s == HS_ESTUCK);
    if (s == HS_OK) CHECK(hs_is_simple(n, k, m, e) == 1 && isfinite(lw));
    int64_t cnt;
    int32_t g[3 * 2000];
    CHECK(hs_gnp(30, 3, 0.1, 5, g, 2000, &cnt) == HS_OK && hs_is_simple(30, 3, cnt, g) == 1);
    CHECK(hs_gnm(30, 2, 100, 9, g) == HS_OK && hs_is_simple(30, 2, 100, g) == 1);
    CHECK(hs_binomial(64, 32) == 1832624140942590534LL);
    CHECK(hs_binomial(70, 35) == -1);
    CHECK(hs_config_regular(8, 3, 12, 1, 1000, e, &trials) == HS_EXHAUSTED);
    printf("smoke ok (config trials %lld, switch accepted %lld)\n", (long long)trials, (long long)acc);
    free(e);
    return 0;
}
