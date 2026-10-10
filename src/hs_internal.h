/* Internal helpers shared by the hypersample sources: RNG and a hash set of k-sets. */
#ifndef HS_INTERNAL_H
#define HS_INTERNAL_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "hypersample.h"

/* ---- RNG: xoshiro256** seeded by splitmix64 (Blackman and Vigna) ---- */
typedef struct { uint64_t s[4]; } hs_rng;

static inline uint64_t hs_splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline void hs_rng_seed(hs_rng *r, uint64_t seed) {
    uint64_t x = seed;
    for (int i = 0; i < 4; i++) r->s[i] = hs_splitmix64(&x);
    if (!(r->s[0] | r->s[1] | r->s[2] | r->s[3])) r->s[0] = 1;
}

static inline uint64_t hs_rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static inline uint64_t hs_next(hs_rng *r) {
    uint64_t *s = r->s, out = hs_rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = hs_rotl(s[3], 45);
    return out;
}

/* Uniform integer in [0, bound), bound >= 1, exact (Lemire's method with rejection). */
static inline uint32_t hs_below32(hs_rng *r, uint32_t bound) {
    uint64_t x = (hs_next(r) >> 32) * (uint64_t)bound;
    uint32_t low = (uint32_t)x;
    if (low < bound) {
        uint32_t t = (uint32_t)(-bound) % bound;
        while (low < t) { x = (hs_next(r) >> 32) * (uint64_t)bound; low = (uint32_t)x; }
    }
    return (uint32_t)(x >> 32);
}

/* Uniform double in (0, 1]. */
static inline double hs_unit_open0(hs_rng *r) {
    return ((double)(hs_next(r) >> 11) + 1.0) * (1.0 / 9007199254740992.0);
}

/* ---- small sort for one edge ---- */
static inline void hs_sort_small(int32_t *a, int k) {
    for (int i = 1; i < k; i++) {
        int32_t v = a[i]; int j = i - 1;
        while (j >= 0 && a[j] > v) { a[j + 1] = a[j]; j--; }
        a[j + 1] = v;
    }
}

/* ---- hash set of sorted k-sets, stored by index into an edge array ---- */
typedef struct {
    int64_t *slot;          /* -1 empty, else edge index */
    int64_t cap;            /* power of two */
    int64_t *used;          /* occupied slots, for O(size) clearing */
    int64_t nused;
    const int32_t *edges;   /* the edge array the indices refer to */
    int k;
} hs_kset;

static inline uint64_t hs_hash_edge(const int32_t *e, int k) {
    uint64_t h = 0x84222325CBF29CE4ULL;
    for (int i = 0; i < k; i++) { h ^= (uint64_t)(uint32_t)e[i]; h *= 0x100000001B3ULL; h ^= h >> 29; }
    return h * 0x9E3779B97F4A7C15ULL;
}

static inline int hs_kset_init(hs_kset *s, int64_t max_items, const int32_t *edges, int k) {
    int64_t cap = 16;
    while (cap < 2 * max_items + 16) cap <<= 1;
    s->slot = (int64_t *)malloc(sizeof(int64_t) * (size_t)cap);
    s->used = (int64_t *)malloc(sizeof(int64_t) * (size_t)(max_items + 1));
    if (!s->slot || !s->used) { free(s->slot); free(s->used); return 0; }
    for (int64_t i = 0; i < cap; i++) s->slot[i] = -1;
    s->cap = cap; s->nused = 0; s->edges = edges; s->k = k;
    return 1;
}

static inline void hs_kset_free(hs_kset *s) { free(s->slot); free(s->used); s->slot = s->used = NULL; }

static inline void hs_kset_clear(hs_kset *s) {
    for (int64_t i = 0; i < s->nused; i++) s->slot[s->used[i]] = -1;
    s->nused = 0;
}

/* Slot holding edge e (sorted), or the empty slot where it would go. *found set accordingly. */
static inline int64_t hs_kset_find(const hs_kset *s, const int32_t *e, int *found) {
    int64_t i = (int64_t)(hs_hash_edge(e, s->k) >> 1) & (s->cap - 1);
    while (s->slot[i] >= 0) {
        if (memcmp(s->edges + s->slot[i] * s->k, e, sizeof(int32_t) * (size_t)s->k) == 0) { *found = 1; return i; }
        i = (i + 1) & (s->cap - 1);
    }
    *found = 0; return i;
}

static inline void hs_kset_put(hs_kset *s, int64_t at, int64_t edge_index) {
    s->slot[at] = edge_index; s->used[s->nused++] = at;
}

/* Remove the entry at slot `at` (linear probing deletion with back-shift). */
static inline void hs_kset_remove_at(hs_kset *s, int64_t at) {
    int64_t i = at, j = at, mask = s->cap - 1;
    s->slot[i] = -1;
    for (;;) {
        j = (j + 1) & mask;
        if (s->slot[j] < 0) return;
        int64_t h = (int64_t)(hs_hash_edge(s->edges + s->slot[j] * s->k, s->k) >> 1) & mask;
        /* move slot j back to i if its home h is not cyclically in (i, j] */
        int in_range = (i <= j) ? (h > i && h <= j) : (h > i || h <= j);
        if (!in_range) { s->slot[i] = s->slot[j]; s->slot[j] = -1; i = j; }
    }
}

static inline int hs_distinct_small(const int32_t *a, int k) {
    for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) if (a[i] == a[j]) return 0;
    return 1;
}

#endif
