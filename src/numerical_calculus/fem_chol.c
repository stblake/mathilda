/* Mathilda -- sparse Cholesky factorisation for finite element matrices.
 * See fem_chol.h. */
#include "fem_chol.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct FemChol {
    int64_t n;
    int64_t* perm;       /* new -> old */
    int64_t* Lp; int64_t* Li; double* Lx;   /* L by columns, diagonal first */
    double* work;
};

/* ------------------------------------------------------------------ */
/*  Nested dissection ordering                                         */
/* ------------------------------------------------------------------ */

typedef struct {
    const int64_t* rp; const int64_t* ci; const double* xy;
    int64_t* side;       /* per node: stamp-coded side of the current split */
    int64_t stamp;
    int64_t* out; int64_t nout;
    double* key;         /* scratch for the median */
} Nd;

static void swap64(int64_t* a, int64_t* b) { int64_t t = *a; *a = *b; *b = t; }

/* Partially order idx[lo..hi] by coordinate axis so idx[k] is the k-th
 * smallest (quickselect, median-of-three pivot). */
static void nd_select(const double* xy, int axis, int64_t* idx, int64_t lo, int64_t hi, int64_t k) {
    while (hi > lo) {
        int64_t mid = lo + (hi - lo) / 2;
        double a = xy[2 * idx[lo] + axis], b = xy[2 * idx[mid] + axis], c = xy[2 * idx[hi] + axis];
        int64_t pi = (a < b) ? ((b < c) ? mid : (a < c ? hi : lo)) : ((a < c) ? lo : (b < c ? hi : mid));
        double pv = xy[2 * idx[pi] + axis];
        int64_t i = lo, j = hi;
        while (i <= j) {
            while (xy[2 * idx[i] + axis] < pv) i++;
            while (xy[2 * idx[j] + axis] > pv) j--;
            if (i <= j) { swap64(&idx[i], &idx[j]); i++; j--; }
        }
        if (k <= j) hi = j;
        else if (k >= i) lo = i;
        else return;
    }
}

static void nd_order(Nd* D, int64_t* S, int64_t m) {
    if (m <= 64) {
        for (int64_t i = 0; i < m; i++) D->out[D->nout++] = S[i];
        return;
    }
    double x0 = INFINITY, x1 = -INFINITY, y0 = INFINITY, y1 = -INFINITY;
    for (int64_t i = 0; i < m; i++) {
        double x = D->xy[2 * S[i]], y = D->xy[2 * S[i] + 1];
        if (x < x0) x0 = x;
        if (x > x1) x1 = x;
        if (y < y0) y0 = y;
        if (y > y1) y1 = y;
    }
    int axis = (x1 - x0) >= (y1 - y0) ? 0 : 1;
    int64_t half = m / 2;
    nd_select(D->xy, axis, S, 0, m - 1, half);
    /* S[0..half-1] left, S[half..m-1] right */
    int64_t sl = ++D->stamp, sr = ++D->stamp;
    for (int64_t i = 0; i < half; i++) D->side[S[i]] = sl;
    for (int64_t i = half; i < m; i++) D->side[S[i]] = sr;
    int64_t* L = malloc(sizeof(int64_t) * (size_t)half);
    int64_t* Sep = malloc(sizeof(int64_t) * (size_t)half);
    int64_t nl = 0, ns = 0;
    for (int64_t i = 0; i < half; i++) {
        int64_t v = S[i];
        bool sep = false;
        for (int64_t p = D->rp[v]; p < D->rp[v + 1] && !sep; p++)
            if (D->side[D->ci[p]] == sr) sep = true;
        if (sep) Sep[ns++] = v; else L[nl++] = v;
    }
    int64_t nr = m - half;
    int64_t* R = malloc(sizeof(int64_t) * (size_t)nr);
    memcpy(R, S + half, sizeof(int64_t) * (size_t)nr);
    if (nl == 0 || nr == 0) {          /* degenerate split: stop dissecting */
        for (int64_t i = 0; i < m; i++) D->out[D->nout++] = S[i];
    } else {
        nd_order(D, L, nl);
        nd_order(D, R, nr);
        for (int64_t i = 0; i < ns; i++) D->out[D->nout++] = Sep[i];
    }
    free(L); free(R); free(Sep);
}

/* ------------------------------------------------------------------ */
/*  Factorisation                                                      */
/* ------------------------------------------------------------------ */

FemChol* fem_chol_factor(int64_t n, const int64_t* rp, const int64_t* ci, const double* x,
                         const double* xy, const char** why) {
    *why = NULL;
    if (n <= 0) { *why = "empty matrix"; return NULL; }
    FemChol* F = calloc(1, sizeof *F);
    F->n = n;
    F->perm = malloc(sizeof(int64_t) * (size_t)n);
    /* ordering */
    {
        Nd D; memset(&D, 0, sizeof D);
        D.rp = rp; D.ci = ci; D.xy = xy;
        D.side = calloc((size_t)n, sizeof(int64_t));
        D.out = F->perm;
        int64_t* S = malloc(sizeof(int64_t) * (size_t)n);
        for (int64_t i = 0; i < n; i++) S[i] = i;
        nd_order(&D, S, n);
        free(S); free(D.side);
    }
    int64_t* pinv = malloc(sizeof(int64_t) * (size_t)n);
    for (int64_t k = 0; k < n; k++) pinv[F->perm[k]] = k;
    /* C = upper triangle of P A P^T by columns */
    int64_t* Cp = calloc((size_t)n + 1, sizeof(int64_t));
    for (int64_t r = 0; r < n; r++)
        for (int64_t p = rp[r]; p < rp[r + 1]; p++) {
            int64_t i = pinv[r], j = pinv[ci[p]];
            if (i <= j) Cp[j + 1]++;
        }
    for (int64_t j = 0; j < n; j++) Cp[j + 1] += Cp[j];
    int64_t* Ci = malloc(sizeof(int64_t) * (size_t)(Cp[n] ? Cp[n] : 1));
    double* Cx = malloc(sizeof(double) * (size_t)(Cp[n] ? Cp[n] : 1));
    int64_t* fill = malloc(sizeof(int64_t) * (size_t)n);
    memcpy(fill, Cp, sizeof(int64_t) * (size_t)n);
    for (int64_t r = 0; r < n; r++)
        for (int64_t p = rp[r]; p < rp[r + 1]; p++) {
            int64_t i = pinv[r], j = pinv[ci[p]];
            if (i <= j) { int64_t q = fill[j]++; Ci[q] = i; Cx[q] = x[p]; }
        }
    free(pinv);
    /* elimination tree */
    int64_t* parent = malloc(sizeof(int64_t) * (size_t)n);
    int64_t* anc = fill;           /* reuse */
    for (int64_t k = 0; k < n; k++) {
        parent[k] = -1; anc[k] = -1;
        for (int64_t p = Cp[k]; p < Cp[k + 1]; p++) {
            int64_t i = Ci[p];
            while (i != -1 && i < k) {
                int64_t inext = anc[i];
                anc[i] = k;
                if (inext == -1) parent[i] = k;
                i = inext;
            }
        }
    }
    /* column counts by the row patterns */
    int64_t* mark = malloc(sizeof(int64_t) * (size_t)n);
    int64_t* s = malloc(sizeof(int64_t) * (size_t)n);
    int64_t* cnt = calloc((size_t)n, sizeof(int64_t));
    for (int64_t k = 0; k < n; k++) mark[k] = -1;
    for (int64_t k = 0; k < n; k++) {
        cnt[k]++;
        mark[k] = k;
        for (int64_t p = Cp[k]; p < Cp[k + 1]; p++) {
            int64_t i = Ci[p];
            if (i > k) continue;
            for (; mark[i] != k; i = parent[i]) { mark[i] = k; cnt[i]++; }
        }
    }
    F->Lp = malloc(sizeof(int64_t) * (size_t)(n + 1));
    F->Lp[0] = 0;
    for (int64_t k = 0; k < n; k++) F->Lp[k + 1] = F->Lp[k] + cnt[k];
    int64_t lnz = F->Lp[n];
    if (lnz > INT64_C(400000000)) {
        *why = "the Cholesky factor would be too large";
        free(Cp); free(Ci); free(Cx); free(parent); free(mark); free(s); free(cnt); free(fill);
        fem_chol_free(F);
        return NULL;
    }
    F->Li = malloc(sizeof(int64_t) * (size_t)lnz);
    F->Lx = malloc(sizeof(double) * (size_t)lnz);
    int64_t* c = cnt;               /* next free slot per column */
    for (int64_t k = 0; k < n; k++) c[k] = F->Lp[k];
    double* w = calloc((size_t)n, sizeof(double));
    for (int64_t k = 0; k < n; k++) mark[k] = -1;
    bool ok = true;
    for (int64_t k = 0; k < n && ok; k++) {
        /* pattern of row k of L, topologically ordered, in s[top..n-1] */
        int64_t top = n;
        mark[k] = k;
        for (int64_t p = Cp[k]; p < Cp[k + 1]; p++) {
            int64_t i = Ci[p];
            if (i > k) continue;
            w[i] += Cx[p];
            int64_t len = 0;
            for (; mark[i] != k; i = parent[i]) { s[len++] = i; mark[i] = k; }
            while (len > 0) s[--top] = s[--len];
        }
        double d = w[k];
        w[k] = 0.0;
        for (; top < n; top++) {
            int64_t i = s[top];
            double lki = w[i] / F->Lx[F->Lp[i]];
            w[i] = 0.0;
            for (int64_t p = F->Lp[i] + 1; p < c[i]; p++) w[F->Li[p]] -= F->Lx[p] * lki;
            d -= lki * lki;
            int64_t p = c[i]++;
            F->Li[p] = k; F->Lx[p] = lki;
        }
        if (!(d > 0.0) || !isfinite(d)) { ok = false; break; }
        int64_t p = c[k]++;
        F->Li[p] = k; F->Lx[p] = sqrt(d);
    }
    free(w); free(Cp); free(Ci); free(Cx); free(parent); free(mark); free(s); free(cnt); free(fill);
    if (!ok) { *why = "the matrix is not positive definite"; fem_chol_free(F); return NULL; }
    F->work = malloc(sizeof(double) * (size_t)n);
    return F;
}

void fem_chol_solve(const FemChol* F, const double* b, double* x) {
    int64_t n = F->n;
    double* y = F->work;
    for (int64_t k = 0; k < n; k++) y[k] = b[F->perm[k]];
    for (int64_t j = 0; j < n; j++) {
        double yj = y[j] / F->Lx[F->Lp[j]];
        y[j] = yj;
        for (int64_t p = F->Lp[j] + 1; p < F->Lp[j + 1]; p++) y[F->Li[p]] -= F->Lx[p] * yj;
    }
    for (int64_t j = n - 1; j >= 0; j--) {
        double t = y[j];
        for (int64_t p = F->Lp[j] + 1; p < F->Lp[j + 1]; p++) t -= F->Lx[p] * y[F->Li[p]];
        y[j] = t / F->Lx[F->Lp[j]];
    }
    for (int64_t k = 0; k < n; k++) x[F->perm[k]] = y[k];
}

int64_t fem_chol_nnz(const FemChol* F) { return F ? F->Lp[F->n] : 0; }

void fem_chol_free(FemChol* F) {
    if (!F) return;
    free(F->perm); free(F->Lp); free(F->Li); free(F->Lx); free(F->work);
    free(F);
}
