/* MATLAB MEX gateway for the hypersample C core (one entry point, string commands).
 *
 * Vertices are 1-based on the MATLAB side and 0-based in C; this file converts.
 * Edge matrices are m-by-k (one edge per row), int32 on output; input edges may be
 * double or int32. Seeds may be uint64 or a double holding an integer < 2^53.
 * Build with build_hypersample.m. Use the hs_*.m wrappers rather than calling this directly.
 */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include "hypersample.h"

static void fail(int status) {
    const char *id = status == HS_EXHAUSTED ? "hypersample:saturation" : "hypersample:status";
    mexErrMsgIdAndTxt(id, "hypersample: %s (status %d)", hs_status_string(status), status);
}

static double scalar(const mxArray *a, const char *what) {
    if (!mxIsNumeric(a) || mxGetNumberOfElements(a) != 1) mexErrMsgIdAndTxt("hypersample:args", "%s must be a numeric scalar", what);
    return mxGetScalar(a);
}

static int32_t iscalar(const mxArray *a, const char *what) {
    double x = scalar(a, what);
    if (x != floor(x) || x < -2147483648.0 || x > 2147483647.0) mexErrMsgIdAndTxt("hypersample:args", "%s must be an int32 integer", what);
    return (int32_t)x;
}

static uint64_t seedarg(const mxArray *a) {
    if (mxIsUint64(a) && mxGetNumberOfElements(a) == 1) return *(const uint64_t *)mxGetData(a);
    double x = scalar(a, "seed");
    if (x < 0 || x != floor(x) || x > 9007199254740992.0) mexErrMsgIdAndTxt("hypersample:args", "seed must be a nonnegative integer (use uint64 above 2^53)");
    return (uint64_t)x;
}

/* numeric vector -> new int32 C array */
static int32_t *ivector(const mxArray *a, mwSize *len, const char *what) {
    if (!mxIsNumeric(a)) mexErrMsgIdAndTxt("hypersample:args", "%s must be numeric", what);
    mwSize L = mxGetNumberOfElements(a);
    int32_t *out = (int32_t *)mxMalloc(sizeof(int32_t) * (L ? L : 1));
    if (mxIsInt32(a)) memcpy(out, mxGetData(a), sizeof(int32_t) * L);
    else if (mxIsDouble(a)) { const double *p = mxGetPr(a); for (mwSize i = 0; i < L; i++) out[i] = (int32_t)p[i]; }
    else mexErrMsgIdAndTxt("hypersample:args", "%s must be double or int32", what);
    *len = L;
    return out;
}

/* m-by-k MATLAB edge matrix (1-based) -> row-major 0-based C array */
static int32_t *edges_in(const mxArray *a, int64_t *m, int32_t *k) {
    if (mxGetNumberOfDimensions(a) != 2) mexErrMsgIdAndTxt("hypersample:args", "edges must be an m-by-k matrix");
    mwSize M = mxGetM(a), K = mxGetN(a), L;
    int32_t *col = ivector(a, &L, "edges");
    int32_t *e = (int32_t *)mxMalloc(sizeof(int32_t) * (M * K ? M * K : 1));
    for (mwSize i = 0; i < M; i++) for (mwSize j = 0; j < K; j++) e[i * K + j] = col[i + j * M] - 1;
    mxFree(col);
    *m = (int64_t)M; *k = (int32_t)K;
    return e;
}

/* row-major 0-based C array -> m-by-k int32 MATLAB matrix (1-based) */
static mxArray *edges_out(const int32_t *e, int64_t m, int32_t k) {
    mxArray *a = mxCreateNumericMatrix((mwSize)m, (mwSize)k, mxINT32_CLASS, mxREAL);
    int32_t *p = (int32_t *)mxGetData(a);
    for (int64_t i = 0; i < m; i++) for (int32_t j = 0; j < k; j++) p[i + j * m] = e[i * k + j] + 1;
    return a;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    char cmd[64];
    if (nrhs < 1 || !mxIsChar(prhs[0]) || mxGetString(prhs[0], cmd, sizeof cmd))
        mexErrMsgIdAndTxt("hypersample:args", "first argument must be a command string");
    if (hs_abi_version() != HS_ABI_VERSION) mexErrMsgIdAndTxt("hypersample:abi", "ABI mismatch: rebuild with build_hypersample");
#define NEED(c) do { if (nrhs != (c)) mexErrMsgIdAndTxt("hypersample:args", "%s: expected %d inputs", cmd, (c) - 1); } while (0)

    if (!strcmp(cmd, "version")) {
        plhs[0] = mxCreateString(hs_version_string());
    } else if (!strcmp(cmd, "config_degrees") || !strcmp(cmd, "sequential")) {
        int seq = !strcmp(cmd, "sequential");
        NEED(seq ? 4 : 5);
        mwSize n; int32_t *deg = ivector(prhs[1], &n, "deg");
        int32_t k = iscalar(prhs[2], "k");
        uint64_t seed = seedarg(prhs[3]);
        int64_t S = 0; for (mwSize v = 0; v < n; v++) S += deg[v];
        if (k < 1 || S % k) mexErrMsgIdAndTxt("hypersample:args", "sum(deg) must be divisible by k");
        int64_t m = S / k;
        int32_t *e = (int32_t *)mxCalloc((size_t)(m * k ? m * k : 1), sizeof(int32_t));
        int s;
        if (seq) {
            double lw = 0;
            s = hs_sequential_degrees((int32_t)n, k, deg, seed, e, &lw);
            if (s != HS_OK && s != HS_ESTUCK) fail(s);
            plhs[0] = edges_out(e, m, k);
            if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(lw);
        } else {
            int64_t trials = 0;
            s = hs_config_degrees((int32_t)n, k, deg, seed, (int64_t)scalar(prhs[4], "maxTrials"), e, &trials);
            if (s) fail(s);
            plhs[0] = edges_out(e, m, k);
            if (nlhs > 1) plhs[1] = mxCreateDoubleScalar((double)trials);
        }
        mxFree(e); mxFree(deg);
    } else if (!strcmp(cmd, "config_regular")) {
        NEED(6);
        int32_t n = iscalar(prhs[1], "n"), k = iscalar(prhs[2], "k"), d = iscalar(prhs[3], "d");
        if (k < 1 || ((int64_t)n * d) % k) mexErrMsgIdAndTxt("hypersample:args", "n*d must be divisible by k");
        int64_t m = (int64_t)n * d / k, trials = 0;
        int32_t *e = (int32_t *)mxCalloc((size_t)(m * k ? m * k : 1), sizeof(int32_t));
        int s = hs_config_regular(n, k, d, seedarg(prhs[4]), (int64_t)scalar(prhs[5], "maxTrials"), e, &trials);
        if (s) fail(s);
        plhs[0] = edges_out(e, m, k);
        if (nlhs > 1) plhs[1] = mxCreateDoubleScalar((double)trials);
        mxFree(e);
    } else if (!strcmp(cmd, "switch_chain")) {
        NEED(5);
        int64_t m; int32_t k; int32_t *e = edges_in(prhs[1], &m, &k);
        int64_t acc = 0;
        int s = hs_switch_chain(iscalar(prhs[2], "n"), k, m, e, (int64_t)scalar(prhs[3], "steps"), seedarg(prhs[4]), &acc);
        if (s) fail(s);
        plhs[0] = edges_out(e, m, k);
        if (nlhs > 1) plhs[1] = mxCreateDoubleScalar((double)acc);
        mxFree(e);
    } else if (!strcmp(cmd, "gnp")) {
        NEED(5);
        int32_t n = iscalar(prhs[1], "n"), k = iscalar(prhs[3], "k");
        double p = scalar(prhs[2], "p");
        uint64_t seed = seedarg(prhs[4]);
        int64_t N = hs_binomial(n, k);
        if (N < 0) fail(HS_EOVERFLOW);
        double mean = p * (double)N;
        int64_t cap = (int64_t)fmin((double)N, 1.2 * mean + 10.0 * sqrt(mean) + 16.0), m = 0;
        int32_t *e = (int32_t *)mxMalloc(sizeof(int32_t) * (size_t)((cap ? cap : 1) * k));
        int s = hs_gnp(n, k, p, seed, e, cap, &m);
        if (s == HS_EBUFFER) {
            mxFree(e); cap = m;
            e = (int32_t *)mxMalloc(sizeof(int32_t) * (size_t)((cap ? cap : 1) * k));
            s = hs_gnp(n, k, p, seed, e, cap, &m);
        }
        if (s) fail(s);
        plhs[0] = edges_out(e, m, k);
        mxFree(e);
    } else if (!strcmp(cmd, "gnm")) {
        NEED(5);
        int32_t n = iscalar(prhs[1], "n"), k = iscalar(prhs[3], "k");
        int64_t m = (int64_t)scalar(prhs[2], "m");
        int32_t *e = (int32_t *)mxMalloc(sizeof(int32_t) * (size_t)((m ? m : 1) * k));
        int s = hs_gnm(n, k, m, seedarg(prhs[4]), e);
        if (s) fail(s);
        plhs[0] = edges_out(e, m, k);
        mxFree(e);
    } else if (!strcmp(cmd, "is_simple") || !strcmp(cmd, "degrees")) {
        NEED(3);
        int64_t m; int32_t k; int32_t *e = edges_in(prhs[1], &m, &k);
        int32_t n = iscalar(prhs[2], "n");
        if (!strcmp(cmd, "is_simple")) {
            int r = hs_is_simple(n, k, m, e);
            if (r < 0) fail(-r);
            plhs[0] = mxCreateLogicalScalar(r == 1);
        } else {
            plhs[0] = mxCreateNumericMatrix((mwSize)n, 1, mxINT32_CLASS, mxREAL);
            int s = hs_degrees(n, k, m, e, (int32_t *)mxGetData(plhs[0]));
            if (s) fail(s);
        }
        mxFree(e);
    } else {
        mexErrMsgIdAndTxt("hypersample:args", "unknown command '%s'", cmd);
    }
}
