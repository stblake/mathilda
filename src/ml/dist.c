/* dist.c -- distribution objects, RandomVariate and PDF. See dist.h. */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "numarray.h"
#include "random.h"
#include "mlutil.h"
#include "encode.h"
#include "pca.h"       /* ml_column_mean */
#include "gmm.h"       /* ml_gmm_fit / ml_gmm_param_count -- extracted in iteration 9 for this */
#include "dist.h"
#include "numeric.h"   /* expr_is_numeric_quantity, numericalize_machine_unplanned */
#include "eval.h"
#include "opform.h"   /* PDF[dist][x] */
#include "common.h"   /* common_scan_inexact */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- Box-Muller ---------------------------------------------------------- */

static bool  bm_have = false;
static double bm_spare = 0.0;

void ml_dist_reset_cache(void) { bm_have = false; bm_spare = 0.0; }

double ml_normal_deviate(void) {
    if (bm_have) { bm_have = false; return bm_spare; }
    /* The polar (Marsaglia) form rather than the trigonometric one: it needs no sin or
     * cos, and rejecting the corners of the square is cheaper than two transcendental
     * calls. The loop also excludes s == 0, where the log would be -infinity. */
    double u, v, s;
    do {
        u = 2.0 * random_uniform_01() - 1.0;
        v = 2.0 * random_uniform_01() - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    double f = sqrt(-2.0 * log(s) / s);
    bm_spare = v * f;                 /* keep the second deviate; see dist.h */
    bm_have = true;
    return u * f;
}

/* ---- Distribution objects ------------------------------------------------ */

/* Recognised distributions. A distribution object is an ordinary expression whose
 * head names the family and whose arguments are its parameters -- specified, not
 * fitted, so it prints in full. */
typedef enum { ML_D_NONE, ML_D_NORMAL, ML_D_UNIFORM, ML_D_MULTINORMAL,
               ML_D_MIXTURE, ML_D_KDE } MlDistKind;

typedef struct {
    MlDistKind kind;
    double a, b;          /* Normal: mu, sigma. Uniform: lo, hi. */
    size_t dim;           /* Multinormal / Mixture */
    const double* mu;     /* borrowed. Mixture: k means, row-major */
    const double* chol;   /* borrowed, lower factor(s) */
    size_t k;             /* Mixture: component count. KDE: sample size */
    const double* w;      /* borrowed. Mixture: k weights */
    const double* logdet; /* borrowed. Mixture: k log-determinants */
    const double* band;   /* borrowed. KDE: per-dimension bandwidths */
} MlDist;

/* Read a LearnedDistribution["SmoothKernel", {bandwidths, sample rows...}, dim, n].
 *
 * The sample IS the model, as it is for a nearest-neighbour predictor: nothing is fitted
 * except the bandwidth, which is why a KDE costs nothing to build and everything to
 * evaluate. */
static bool ml_read_kde(Expr* e, MlDist* d, double** owned) {
    Expr* pay  = e->data.function.args[1];
    Expr* dimx = e->data.function.args[2];
    Expr* nx   = e->data.function.args[3];
    if (!pay || pay->type != EXPR_FUNCTION) return false;
    if (!dimx || dimx->type != EXPR_INTEGER || dimx->data.integer <= 0) return false;
    if (!nx || nx->type != EXPR_INTEGER || nx->data.integer <= 0) return false;
    size_t dim = (size_t)dimx->data.integer, n = (size_t)nx->data.integer;
    if (pay->data.function.arg_count != n + 1) return false;

    double* buf = malloc(sizeof(double) * (dim + n * dim));
    if (!buf) return false;
    double* band = buf, *smp = buf + dim;
    double im = 0.0; bool ok = true;
    Expr* br = pay->data.function.args[0];
    if (!br || br->type != EXPR_FUNCTION || br->data.function.arg_count != dim) ok = false;
    for (size_t a = 0; ok && a < dim; a++) {
        ok = na_read_scalar(br->data.function.args[a], &band[a], &im) && im == 0.0;
        if (ok && !(band[a] > 0.0)) ok = false;      /* a zero bandwidth has no density */
    }
    for (size_t i = 0; ok && i < n; i++) {
        Expr* r = pay->data.function.args[i + 1];
        if (!r || r->type != EXPR_FUNCTION || r->data.function.arg_count != dim) {
            ok = false; break;
        }
        for (size_t a = 0; ok && a < dim; a++)
            ok = na_read_scalar(r->data.function.args[a], &smp[i * dim + a], &im)
              && im == 0.0;
    }
    if (!ok) { free(buf); return false; }
    d->kind = ML_D_KDE; d->dim = dim; d->k = n; d->band = band; d->mu = smp;
    *owned = buf;
    return true;
}

/* Read a LearnedDistribution["GaussianMixture", {weights, (mean, covRows...) x k}, dim, k].
 *
 * The Cholesky factors and log-determinants are DERIVED, so they are recomputed here
 * rather than stored: keeping them in the payload would double its size and create two
 * places for the same fact to live. That is the same choice the Multinormal path makes.
 *
 * A component whose covariance will not factorise makes the whole density undefined, so
 * this declines. ml_gmm_fit cannot produce one -- its ridge guarantees positive
 * definiteness -- but a hand-written LearnedDistribution could. */
static bool ml_read_mixture(Expr* e, MlDist* d, double** owned) {
    Expr* pay  = e->data.function.args[1];
    Expr* dimx = e->data.function.args[2];
    Expr* kx   = e->data.function.args[3];
    if (!pay || pay->type != EXPR_FUNCTION) return false;
    if (!dimx || dimx->type != EXPR_INTEGER || dimx->data.integer <= 0) return false;
    if (!kx || kx->type != EXPR_INTEGER || kx->data.integer <= 0) return false;
    size_t dim = (size_t)dimx->data.integer, k = (size_t)kx->data.integer;
    if (pay->data.function.arg_count != 1 + k * (1 + dim)) return false;

    /* One allocation: weights | means | chol | logdet. */
    double* buf = malloc(sizeof(double) * (k + k * dim + k * dim * dim + k));
    double* cov = malloc(sizeof(double) * dim * dim);
    if (!buf || !cov) { free(buf); free(cov); return false; }
    double* w = buf, *mu = buf + k, *ch = mu + k * dim, *ld = ch + k * dim * dim;

    double im = 0.0; bool ok = true;
    Expr* wr = pay->data.function.args[0];
    if (!wr || wr->type != EXPR_FUNCTION || wr->data.function.arg_count != k) ok = false;
    for (size_t j = 0; ok && j < k; j++)
        ok = na_read_scalar(wr->data.function.args[j], &w[j], &im) && im == 0.0;

    for (size_t j = 0; ok && j < k; j++) {
        size_t base = 1 + j * (1 + dim);
        Expr* mr = pay->data.function.args[base];
        if (!mr || mr->type != EXPR_FUNCTION || mr->data.function.arg_count != dim) {
            ok = false; break;
        }
        for (size_t a = 0; ok && a < dim; a++)
            ok = na_read_scalar(mr->data.function.args[a], &mu[j * dim + a], &im)
              && im == 0.0;
        for (size_t a = 0; ok && a < dim; a++) {
            Expr* cr = pay->data.function.args[base + 1 + a];
            if (!cr || cr->type != EXPR_FUNCTION || cr->data.function.arg_count != dim) {
                ok = false; break;
            }
            for (size_t b = 0; ok && b < dim; b++)
                ok = na_read_scalar(cr->data.function.args[b], &cov[a * dim + b], &im)
                  && im == 0.0;
        }
        if (ok) ok = ml_chol(cov, ch + j * dim * dim, dim);
        if (ok) {
            double t = 0.0;
            for (size_t a = 0; a < dim; a++) t += log(ch[j * dim * dim + a * dim + a]);
            ld[j] = 2.0 * t;
        }
    }
    free(cov);
    if (!ok) { free(buf); return false; }
    d->kind = ML_D_MIXTURE; d->dim = dim; d->k = k;
    d->w = w; d->mu = mu; d->chol = ch; d->logdet = ld;
    *owned = buf;
    return true;
}

/* Read a LearnedDistribution["Multinormal", {mean, covRows...}, dim, 0] into d.
 *
 * `owned` receives one allocation holding the mean followed by the Cholesky factor, so
 * the caller frees exactly one thing. Returns false if the covariance will not
 * factorise -- which is information rather than an error, since it means the fitted
 * covariance is singular (fewer points than dimensions, or collinear ones) and no
 * density exists. */
static bool ml_read_learned(Expr* e, MlDist* d, double** owned) {
    size_t argc = e->data.function.arg_count;
    if (argc != 4) return false;
    Expr* mname = e->data.function.args[0];
    Expr* pay   = e->data.function.args[1];
    Expr* dimx  = e->data.function.args[2];
    if (!mname || mname->type != EXPR_STRING) return false;
    if (strcmp(mname->data.string, "GaussianMixture") == 0)
        return ml_read_mixture(e, d, owned);
    if (strcmp(mname->data.string, "SmoothKernel") == 0)
        return ml_read_kde(e, d, owned);
    if (strcmp(mname->data.string, "Multinormal") != 0) return false;
    if (!pay || pay->type != EXPR_FUNCTION) return false;
    if (!dimx || dimx->type != EXPR_INTEGER || dimx->data.integer <= 0) return false;
    size_t dim = (size_t)dimx->data.integer;
    if (pay->data.function.arg_count != dim + 1) return false;

    double* buf = malloc(sizeof(double) * (dim + dim * dim + dim * dim));
    if (!buf) return false;
    double* mu = buf, *cov = buf + dim, *chol = buf + dim + dim * dim;
    double im = 0.0; bool ok = true;
    Expr* mrow = pay->data.function.args[0];
    if (!mrow || mrow->type != EXPR_FUNCTION || mrow->data.function.arg_count != dim)
        ok = false;
    for (size_t j = 0; ok && j < dim; j++)
        ok = na_read_scalar(mrow->data.function.args[j], &mu[j], &im) && im == 0.0;
    for (size_t a = 0; ok && a < dim; a++) {
        Expr* cr = pay->data.function.args[a + 1];
        if (!cr || cr->type != EXPR_FUNCTION || cr->data.function.arg_count != dim) {
            ok = false; break;
        }
        for (size_t b = 0; ok && b < dim; b++)
            ok = na_read_scalar(cr->data.function.args[b], &cov[a * dim + b], &im)
              && im == 0.0;
    }
    if (ok) ok = ml_chol(cov, chol, dim);
    if (!ok) { free(buf); return false; }
    d->kind = ML_D_MULTINORMAL; d->dim = dim; d->mu = mu; d->chol = chol;
    *owned = buf;
    return true;
}

/* A real distribution parameter: a number, or an exact numeric constant (Pi,
 * Sqrt[2]) at its machine value -- RandomVariate[UniformDistribution[{0, Pi}]]
 * is a draw from [0, 3.14159...], not an unevaluated call. */
static bool ml_read_param(Expr* e, double* v) {
    double im = 0.0;
    if (na_read_scalar(e, v, &im)) return im == 0.0;
    if (!expr_is_numeric_quantity(e)) return false;
    Expr* nv = numericalize_machine_unplanned(e);
    bool ok = nv && nv->type == EXPR_REAL;
    if (ok) *v = nv->data.real;
    if (nv) expr_free(nv);
    return ok;
}

static bool ml_read_dist(Expr* e, MlDist* d, double** owned) {
    *owned = NULL;
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    const char* hn = h->data.symbol.name;
    size_t argc = e->data.function.arg_count;

    if (hn == SYM_LearnedDistribution) return ml_read_learned(e, d, owned);

    if (hn == SYM_NormalDistribution) {
        d->kind = ML_D_NORMAL; d->dim = 1;
        if (argc == 0) { d->a = 0.0; d->b = 1.0; return true; }   /* the standard normal */
        if (argc != 2) return false;
        if (!ml_read_param(e->data.function.args[0], &d->a)) return false;
        if (!ml_read_param(e->data.function.args[1], &d->b)) return false;
        /* A non-positive standard deviation is not a distribution. Declining beats
         * returning NaNs that propagate silently through a whole sample. */
        return d->b > 0.0;
    }
    if (hn == SYM_UniformDistribution) {
        d->kind = ML_D_UNIFORM; d->dim = 1;
        if (argc == 0) { d->a = 0.0; d->b = 1.0; return true; }
        if (argc != 1) return false;
        Expr* r = e->data.function.args[0];
        if (!r || r->type != EXPR_FUNCTION || r->data.function.arg_count != 2) return false;
        if (!ml_read_param(r->data.function.args[0], &d->a)) return false;
        if (!ml_read_param(r->data.function.args[1], &d->b)) return false;
        return d->b > d->a;
    }
    return false;
}

/* ---- PDF ---------------------------------------------------------------- */

static bool ml_pdf_at(const MlDist* d, double x, double* out) {
    switch (d->kind) {
        case ML_D_NORMAL: {
            double z = (x - d->a) / d->b;
            *out = exp(-0.5 * z * z) / (d->b * sqrt(2.0 * M_PI));
            return true;
        }
        case ML_D_UNIFORM:
            /* Zero strictly outside, 1/(b-a) inside. The endpoints are included, which
             * is Wolfram's convention for a continuous uniform. */
            *out = (x < d->a || x > d->b) ? 0.0 : 1.0 / (d->b - d->a);
            return true;
        default: return false;
    }
}

/* Multinormal log-density at a dim-vector. Kept in log space and exponentiated once,
 * for the reason gmm.c's E-step is: a product of dim Gaussian factors underflows for a
 * point a few standard deviations out, and the log form does not. */
static bool ml_multinormal_pdf(const MlDist* d, const double* x, double* out) {
    double* y = malloc(sizeof(double) * d->dim);
    if (!y) return false;
    double q = ml_mahalanobis(d->chol, d->mu, x, d->dim, y);
    free(y);
    double logdet = 0.0;
    for (size_t a = 0; a < d->dim; a++) logdet += log(d->chol[a * d->dim + a]);
    logdet *= 2.0;
    *out = exp(-0.5 * ((double)d->dim * log(2.0 * M_PI) + logdet + q));
    return true;
}

/* Mixture density: log-sum-exp over the components, exponentiated once. The log form is
 * not optional here -- each component's density is a product of dim factors, so a point
 * a few standard deviations from every mean underflows in linear space, and the answer
 * would be a flat zero exactly where the tail behaviour matters. */
static bool ml_mixture_pdf(const MlDist* d, const double* x, double* out) {
    double* y = malloc(sizeof(double) * d->dim);
    double* lp = malloc(sizeof(double) * d->k);
    if (!y || !lp) { free(y); free(lp); return false; }
    double best = -INFINITY;
    for (size_t j = 0; j < d->k; j++) {
        if (!(d->w[j] > 0.0)) { lp[j] = -INFINITY; continue; }
        double q = ml_mahalanobis(d->chol + j * d->dim * d->dim,
                                 d->mu + j * d->dim, x, d->dim, y);
        lp[j] = log(d->w[j])
              - 0.5 * ((double)d->dim * log(2.0 * M_PI) + d->logdet[j] + q);
        if (lp[j] > best) best = lp[j];
    }
    double sacc = 0.0;
    if (best > -INFINITY)
        for (size_t j = 0; j < d->k; j++)
            if (lp[j] > -INFINITY) sacc += exp(lp[j] - best);
    free(y); free(lp);
    *out = (sacc > 0.0) ? exp(best + log(sacc)) : 0.0;
    return true;
}

/* KDE density: the mean of product-Gaussian kernels centred on the sample points.
 *
 * A PRODUCT (diagonal) kernel with a per-dimension bandwidth rather than a full-covariance
 * one. That is the standard choice, and the honest one here: a full-covariance kernel would
 * need a bandwidth MATRIX, and estimating one from the same sample it smooths is a
 * different and much harder problem than the normal-reference rule below solves.
 *
 * Summed in log space per kernel and combined by log-sum-exp, for the reason the mixture
 * is: with dim factors per kernel, a point several bandwidths from every sample point
 * underflows in linear space, and the tail would read as a flat zero. */
static bool ml_kde_pdf(const MlDist* d, const double* x, double* out) {
    size_t n = d->k, dim = d->dim;
    double lognorm = -log((double)n);
    for (size_t a = 0; a < dim; a++)
        lognorm -= log(d->band[a] * sqrt(2.0 * M_PI));
    double best = -INFINITY;
    double* lp = malloc(sizeof(double) * n);
    if (!lp) return false;
    for (size_t i = 0; i < n; i++) {
        double q = 0.0;
        for (size_t a = 0; a < dim; a++) {
            double z = (x[a] - d->mu[i * dim + a]) / d->band[a];
            q += z * z;
        }
        lp[i] = -0.5 * q;
        if (lp[i] > best) best = lp[i];
    }
    double sacc = 0.0;
    if (best > -INFINITY)
        for (size_t i = 0; i < n; i++) sacc += exp(lp[i] - best);
    free(lp);
    *out = (sacc > 0.0) ? exp(lognorm + best + log(sacc)) : 0.0;
    return true;
}

/* PDF of a ContingencyTable: look the outcome up in the stored table.
 *
 * Separate from the MlDist path on purpose. MlDist holds doubles, and the argument here is
 * a nominal outcome -- a string, a symbol, a pair of them -- so there is nothing to read
 * into a double buffer. Reusing the numeric reader would mean either declining the input
 * this method exists for, or coercing it into numbers it is not.
 *
 * Returns true and sets *out when `dist` IS a ContingencyTable (including the miss case,
 * where the probability is an exact 0); false means "not this method, try the others". */
static bool ct_pdf(Expr* dist, Expr* xe, double* out) {
    if (!dist || dist->type != EXPR_FUNCTION || dist->data.function.arg_count != 4)
        return false;
    Expr* h = dist->data.function.head;
    if (!h || h->type != EXPR_SYMBOL
        || h->data.symbol.name != SYM_LearnedDistribution) return false;
    Expr* m = dist->data.function.args[0];
    if (!m || m->type != EXPR_STRING
        || strcmp(m->data.string, "ContingencyTable") != 0) return false;

    Expr* pay = dist->data.function.args[1];
    if (!pay || pay->type != EXPR_FUNCTION || pay->data.function.arg_count != 2) return false;
    Expr* outs = pay->data.function.args[0];
    Expr* prs  = pay->data.function.args[1];
    if (!outs || outs->type != EXPR_FUNCTION || !prs || prs->type != EXPR_FUNCTION)
        return false;
    if (outs->data.function.arg_count != prs->data.function.arg_count) return false;

    for (size_t k = 0; k < outs->data.function.arg_count; k++) {
        if (expr_eq(outs->data.function.args[k], xe)) {
            double v = 0.0, im = 0.0;
            if (!na_read_scalar(prs->data.function.args[k], &v, &im) || im != 0.0)
                return false;
            *out = v;
            return true;
        }
    }
    /* Never observed. Exactly 0, for the reason the fit records: smoothing would need to
     * know how many outcomes were possible but unseen, which for arbitrary expressions is
     * unknowable. */
    *out = 0.0;
    return true;
}

static Expr* pdf_numeric(Expr* res) {
    /* Nominal first: a ContingencyTable's argument is an outcome, not a number. */
    {
        double cp = 0.0;
        if (ct_pdf(res->data.function.args[0], res->data.function.args[1], &cp))
            return expr_new_real(cp);
    }
    MlDist d; double* owned = NULL;
    if (!ml_read_dist(res->data.function.args[0], &d, &owned)) return NULL;
    Expr* xe = res->data.function.args[1];
    double x = 0.0, im = 0.0;

    if (d.kind == ML_D_MULTINORMAL || d.kind == ML_D_MIXTURE
        || d.kind == ML_D_KDE) {
        /* The argument is a POINT, so a list here is one observation rather than many.
         * That is the opposite reading from the scalar case below, where a list is many
         * points -- and it has to be, because a multinormal's argument is itself a
         * list. Getting this backwards would silently treat each coordinate as a
         * separate observation. */
        size_t n, dm; double* px = NULL; bool vec = false;
        Expr* out = NULL;
        bool (*pdf)(const MlDist*, const double*, double*) =
            (d.kind == ML_D_MIXTURE) ? ml_mixture_pdf :
            (d.kind == ML_D_KDE)     ? ml_kde_pdf : ml_multinormal_pdf;
        /* A univariate learned distribution takes a scalar point, as the
         * specified ones do: PDF[m, 3.5] is PDF[m, {3.5}]. */
        if (d.dim == 1 && na_read_scalar(xe, &x, &im) && im == 0.0) {
            double p = 0.0;
            if (pdf(&d, &x, &p)) out = expr_new_real(p);
            free(owned);
            return out;
        }
        if (ml_read_data(xe, &n, &dm, &px, &vec)) {
            double p = 0.0;
            if (vec && n == d.dim && pdf(&d, px, &p))
                out = expr_new_real(p);
            else if (!vec && dm == d.dim) {         /* a matrix: many points */
                double* o = malloc(sizeof(double) * n);
                bool ok = o != NULL;
                for (size_t i = 0; ok && i < n; i++)
                    ok = pdf(&d, px + i * dm, &o[i]);
                if (ok) out = ml_list_of_reals(o, n);
                free(o);
            }
            free(px);
        }
        free(owned);
        return out;
    }

    /* A list of points threads, giving a list of densities -- the shape a caller
     * plotting a density actually wants. */
    if (xe && xe->type == EXPR_FUNCTION && xe->data.function.head->type == EXPR_SYMBOL
        && xe->data.function.head->data.symbol.name == SYM_List) {
        size_t n = xe->data.function.arg_count;
        double* o = malloc(sizeof(double) * (n ? n : 1));
        if (!o) { free(owned); return NULL; }
        bool ok = true;
        for (size_t i = 0; i < n && ok; i++)
            ok = na_read_scalar(xe->data.function.args[i], &x, &im) && im == 0.0
              && ml_pdf_at(&d, x, &o[i]);
        Expr* out = ok ? ml_list_of_reals(o, n) : NULL;
        free(o); free(owned);
        return out;
    }
    if (!na_read_scalar(xe, &x, &im) || im != 0.0) { free(owned); return NULL; }
    double p = 0.0;
    bool ok = ml_pdf_at(&d, x, &p);
    free(owned);
    return ok ? expr_new_real(p) : NULL;
}

/* True iff e holds an approximate number anywhere. */
static bool ml_has_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (ml_has_inexact(e->data.function.args[i])) return true;
    return false;
}

/* PDF[dist, x]. Exact or symbolic input to a specified distribution gets the
 * closed form (dist_formula.c) -- PDF[NormalDistribution[], 0] is 1/Sqrt[2 Pi],
 * and PDF[NormalDistribution[m, s], x] the density formula; approximate input
 * takes the machine path, and falls back to the closed form when that path
 * cannot read it (a symbolic coordinate, a Pi parameter). */
static Expr* builtin_pdf(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* dist = res->data.function.args[0];
    Expr* xe = res->data.function.args[1];
    bool exact = !ml_has_inexact(dist) && !ml_has_inexact(xe);
    /* Arbitrary-precision input keeps its precision: the closed form,
     * evaluated at the input precision, not the machine density. */
    bool bignum = !exact && common_scan_inexact(res).min_bits > 53;
    if (exact || bignum) {
        Expr* cf = ml_dist_closed_form(ML_PROP_PDF, dist, xe);
        if (cf) return cf;
    }
    Expr* r = pdf_numeric(res);
    if (r || exact) return r;
    return ml_dist_closed_form(ML_PROP_PDF, dist, xe);
}

/* CDF[dist, x]: the closed form for the specified distributions; a learned
 * one-dimensional SmoothKernel gives the mean of the kernel CDFs. */
static Expr* kde_cdf(Expr* dist, Expr* xe);
static Expr* builtin_cdf(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* r = ml_dist_closed_form(ML_PROP_CDF, res->data.function.args[0],
                                  res->data.function.args[1]);
    if (r) return r;
    return kde_cdf(res->data.function.args[0], res->data.function.args[1]);
}

/* ---- Univariate SmoothKernel properties ---------------------------------- *
 * The KDE is the uniform mixture of N(x_i, h^2) over the n sample points, so
 *   Mean     = mean(x_i)
 *   Variance = (1/n) sum (x_i - mean)^2 + h^2      (the mixture's variance)
 *   CDF[x]   = (1/n) sum Phi((x - x_i)/h)
 *   Quantile = the root of CDF[x] = q, by bisection (CDF is strictly monotone).
 * Machine numbers: the distribution itself is machine-numeric. */

/* The 1-D KDE behind dist, or false. On success *owned must be freed. */
static bool kde1(Expr* dist, MlDist* d, double** owned) {
    if (!ml_read_dist(dist, d, owned)) return false;
    if (d->kind == ML_D_KDE && d->dim == 1 && d->k > 0) return true;
    free(*owned);
    *owned = NULL;
    return false;
}

static double kde1_cdf_at(const MlDist* d, double x) {
    double s = 0.0, h = d->band[0];
    for (size_t i = 0; i < d->k; i++)
        s += 0.5 * erfc(-(x - d->mu[i]) / (h * sqrt(2.0)));
    return s / (double)d->k;
}

static Expr* kde_cdf(Expr* dist, Expr* xe) {
    MlDist d; double* owned = NULL;
    double x, im = 0.0;
    if (!na_read_scalar(xe, &x, &im) || im != 0.0) return NULL;
    if (!kde1(dist, &d, &owned)) return NULL;
    Expr* out = expr_new_real(kde1_cdf_at(&d, x));
    free(owned);
    return out;
}

/* Mean / Variance / Median / Quantile of a 1-D SmoothKernel, or NULL. */
static Expr* kde_property(MlDistProp prop, Expr* dist, const Expr* x) {
    MlDist d; double* owned = NULL;
    double q = 0.5, im = 0.0;
    if (prop == ML_PROP_QUANTILE &&
        (!x || !na_read_scalar((Expr*)x, &q, &im) || im != 0.0 || !(q > 0.0 && q < 1.0)))
        return NULL;
    if (!kde1(dist, &d, &owned)) return NULL;
    double mean = 0.0, var = 0.0, h = d.band[0];
    for (size_t i = 0; i < d.k; i++) mean += d.mu[i];
    mean /= (double)d.k;
    for (size_t i = 0; i < d.k; i++) var += (d.mu[i] - mean) * (d.mu[i] - mean);
    var = var / (double)d.k + h * h;
    Expr* out = NULL;
    if (prop == ML_PROP_MEAN) out = expr_new_real(mean);
    else if (prop == ML_PROP_VARIANCE) out = expr_new_real(var);
    else if (prop == ML_PROP_MEDIAN || prop == ML_PROP_QUANTILE) {
        double lo = d.mu[0], hi = d.mu[0];
        for (size_t i = 1; i < d.k; i++) {
            if (d.mu[i] < lo) lo = d.mu[i];
            if (d.mu[i] > hi) hi = d.mu[i];
        }
        lo -= 40.0 * h; hi += 40.0 * h;         /* CDF is 0 / 1 to rounding here */
        for (int it = 0; it < 200 && hi - lo > 4.0 * DBL_EPSILON * (fabs(lo) + fabs(hi)); it++) {
            double mid = 0.5 * (lo + hi);
            if (kde1_cdf_at(&d, mid) < q) lo = mid; else hi = mid;
        }
        out = expr_new_real(0.5 * (lo + hi));
    }
    free(owned);
    return out;
}

/* InverseCDF[dist, q]: for these continuous univariate distributions the
 * quantile function, so Quantile[dist, q]. */
static Expr* builtin_inverse_cdf(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    return ml_dist_property(ML_PROP_QUANTILE, res->data.function.args[0],
                            res->data.function.args[1]);
}

/* RandomVariate[UniformDistribution[{{a1, b1}, ...}]] (n draws when n given):
 * a point drawn uniformly from the box, one uniform per coordinate.  NULL
 * when dist is not a numeric box. */
static Expr* box_draw_one(const double* lo, const double* w, size_t d) {
    Expr** el = malloc(sizeof(Expr*) * d);
    if (!el) return NULL;
    for (size_t i = 0; i < d; i++) el[i] = expr_new_real(lo[i] + w[i] * random_uniform_01());
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), el, d);
    free(el);
    return out;
}
static Expr* box_randomvariate(const Expr* dist, bool one, size_t n) {
    if (!dist || dist->type != EXPR_FUNCTION || dist->data.function.arg_count != 1 ||
        dist->data.function.head->type != EXPR_SYMBOL ||
        dist->data.function.head->data.symbol.name != SYM_UniformDistribution) return NULL;
    const Expr* r = dist->data.function.args[0];
    if (r->type != EXPR_FUNCTION || r->data.function.head->type != EXPR_SYMBOL ||
        r->data.function.head->data.symbol.name != SYM_List || r->data.function.arg_count == 0)
        return NULL;
    size_t d = r->data.function.arg_count;
    double* lo = malloc(sizeof(double) * d);
    double* w = malloc(sizeof(double) * d);
    bool ok = lo && w;
    for (size_t i = 0; ok && i < d; i++) {
        const Expr* p = r->data.function.args[i];
        double a = 0, b = 0, im = 0;
        ok = p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL &&
             p->data.function.head->data.symbol.name == SYM_List && p->data.function.arg_count == 2 &&
             na_read_scalar(p->data.function.args[0], &a, &im) && im == 0.0 &&
             na_read_scalar(p->data.function.args[1], &b, &im) && im == 0.0 && a < b;
        if (ok) { lo[i] = a; w[i] = b - a; }
    }
    Expr* out = NULL;
    if (ok) {
        if (one) out = box_draw_one(lo, w, d);
        else {
            Expr** el = malloc(sizeof(Expr*) * (n ? n : 1));
            if (el) {
                for (size_t i = 0; i < n; i++) el[i] = box_draw_one(lo, w, d);
                out = expr_new_function(expr_new_symbol(SYM_List), el, n);
                free(el);
            }
        }
    }
    free(lo); free(w);
    return out;
}

Expr* ml_dist_property(MlDistProp prop, Expr* dist, const Expr* x) {
    Expr* r = ml_dist_closed_form(prop, dist, x);
    if (r) return r;
    if (prop == ML_PROP_CDF) return x ? kde_cdf(dist, (Expr*)x) : NULL;
    if (prop == ML_PROP_PDF) return NULL;
    return kde_property(prop, dist, x);
}

/* ---- RandomVariate ------------------------------------------------------ */

/* A uniform index in [0, n): floor(u n), with the u -> 1 rounding edge folded back. */
static size_t ml_draw_index(size_t n) {
    size_t i = (size_t)(random_uniform_01() * (double)n);
    return i < n ? i : n - 1;
}

/* mu + L z for z a standard normal dim-vector and L the lower Cholesky factor of the
 * covariance: Cov(L z) = L L' = Sigma. */
static void ml_draw_gauss(const double* mu, const double* l, size_t dim, double* z,
                          double* out) {
    for (size_t a = 0; a < dim; a++) z[a] = ml_normal_deviate();
    for (size_t a = 0; a < dim; a++) {
        double v = mu[a];
        for (size_t b = 0; b <= a; b++) v += l[a * dim + b] * z[b];
        out[a] = v;
    }
}

/* One draw of d->dim values into out. `z` is dim scratch doubles.
 *
 * Every kind ml_read_dist accepts is sampled here; there is no default value. The learned
 * kinds used to fall through to a constant 0.0, which is not a sample of anything.
 *   Multinormal   mu + L z.
 *   Mixture       pick component j with probability w_j (inverse CDF of the weights), then
 *                 mu_j + L_j z -- the generative definition of a mixture.
 *   SmoothKernel  pick a sample point uniformly, add an independent N(0, h_a) per dimension:
 *                 the KDE is exactly the uniform mixture of those product kernels. */
static void ml_draw(const MlDist* d, double* z, double* out) {
    switch (d->kind) {
        case ML_D_NORMAL:  out[0] = d->a + d->b * ml_normal_deviate(); return;
        case ML_D_UNIFORM: out[0] = d->a + (d->b - d->a) * random_uniform_01(); return;
        case ML_D_MULTINORMAL:
            ml_draw_gauss(d->mu, d->chol, d->dim, z, out);
            return;
        case ML_D_MIXTURE: {
            double tot = 0.0;
            for (size_t j = 0; j < d->k; j++) if (d->w[j] > 0.0) tot += d->w[j];
            double u = random_uniform_01() * tot, acc = 0.0;
            size_t pick = d->k - 1;
            for (size_t j = 0; j < d->k; j++) {
                if (!(d->w[j] > 0.0)) continue;
                acc += d->w[j];
                pick = j;
                if (u < acc) break;
            }
            ml_draw_gauss(d->mu + pick * d->dim, d->chol + pick * d->dim * d->dim, d->dim,
                          z, out);
            return;
        }
        case ML_D_KDE: {
            size_t i = ml_draw_index(d->k);
            for (size_t a = 0; a < d->dim; a++)
                out[a] = d->mu[i * d->dim + a] + d->band[a] * ml_normal_deviate();
            return;
        }
        default:
            return;     /* unreachable: ml_read_dist only accepts the kinds above */
    }
}

/* RandomVariate[dist] / RandomVariate[dist, n]. A univariate distribution gives a real (a list
 * of n reals); a dim > 1 one gives a dim-vector (an n x dim matrix). */
/* The outcome and probability lists of a ContingencyTable LearnedDistribution
 * (borrowed), or false when dist is not one. */
static bool ct_table(Expr* dist, Expr** outs, Expr** prs) {
    if (!dist || dist->type != EXPR_FUNCTION || dist->data.function.arg_count != 4) return false;
    Expr* h = dist->data.function.head;
    if (!h || h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_LearnedDistribution) return false;
    Expr* m = dist->data.function.args[0];
    if (!m || m->type != EXPR_STRING || strcmp(m->data.string, "ContingencyTable") != 0) return false;
    Expr* pay = dist->data.function.args[1];
    if (!pay || pay->type != EXPR_FUNCTION || pay->data.function.arg_count != 2) return false;
    *outs = pay->data.function.args[0];
    *prs = pay->data.function.args[1];
    return (*outs)->type == EXPR_FUNCTION && (*prs)->type == EXPR_FUNCTION &&
           (*outs)->data.function.arg_count == (*prs)->data.function.arg_count &&
           (*outs)->data.function.arg_count > 0;
}

/* One outcome drawn from a ContingencyTable: invert the cumulative of its
 * empirical probabilities at one uniform draw (so SeedRandom reproduces it). */
static Expr* ct_draw(Expr* outs, const double* cum, size_t k) {
    double u = random_uniform_01() * cum[k - 1];
    size_t j = 0;
    while (j + 1 < k && u >= cum[j]) j++;
    return expr_copy(outs->data.function.args[j]);
}

/* RandomVariate of a nominal LearnedDistribution: one outcome, or a List of
 * n.  NULL when the table's probabilities do not read. */
static Expr* ct_randomvariate(Expr* outs, Expr* prs, bool one, size_t n) {
    size_t k = outs->data.function.arg_count;
    double* cum = malloc(sizeof(double) * k);
    if (!cum) return NULL;
    double acc = 0.0;
    for (size_t j = 0; j < k; j++) {
        double v = 0.0, im = 0.0;
        if (!na_read_scalar(prs->data.function.args[j], &v, &im) || im != 0.0 || v < 0.0) {
            free(cum); return NULL;
        }
        acc += v; cum[j] = acc;
    }
    Expr* out = NULL;
    if (acc > 0.0) {
        if (one) out = ct_draw(outs, cum, k);
        else {
            Expr** el = malloc(sizeof(Expr*) * (n ? n : 1));
            if (el) {
                for (size_t i = 0; i < n; i++) el[i] = ct_draw(outs, cum, k);
                out = expr_new_function(expr_new_symbol(SYM_List), el, n);
                free(el);
            }
        }
    }
    free(cum);
    return out;
}

static Expr* builtin_randomvariate(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    {
        /* a nominal distribution draws outcomes, not numbers */
        Expr *outs, *prs;
        if (ct_table(res->data.function.args[0], &outs, &prs)) {
            if (argc == 1) return ct_randomvariate(outs, prs, true, 0);
            Expr* ne = res->data.function.args[1];
            if (ne->type != EXPR_INTEGER || ne->data.integer < 0 ||
                ne->data.integer > INT64_C(100000000)) return NULL;
            return ct_randomvariate(outs, prs, false, (size_t)ne->data.integer);
        }
        /* the uniform distribution on a box draws points */
        if (argc == 1) {
            Expr* b = box_randomvariate(res->data.function.args[0], true, 0);
            if (b) return b;
        } else if (res->data.function.args[1]->type == EXPR_INTEGER &&
                   res->data.function.args[1]->data.integer >= 0 &&
                   res->data.function.args[1]->data.integer <= INT64_C(100000000)) {
            Expr* b = box_randomvariate(res->data.function.args[0], false,
                                        (size_t)res->data.function.args[1]->data.integer);
            if (b) return b;
        }
    }
    size_t n = 1;
    /* RandomVariate[dist, {n1, n2, ...}]: an n1 x n2 x ... array of draws. */
    size_t shape[16];
    size_t rank = 0;
    if (argc == 2) {
        Expr* ne = res->data.function.args[1];
        if (ne && ne->type == EXPR_FUNCTION && ne->data.function.head->type == EXPR_SYMBOL
            && ne->data.function.head->data.symbol.name == SYM_List) {
            rank = ne->data.function.arg_count;
            if (rank == 0 || rank > 16) return NULL;
            for (size_t i = 0; i < rank; i++) {
                Expr* di = ne->data.function.args[i];
                if (di->type != EXPR_INTEGER || di->data.integer < 0) return NULL;
                shape[i] = (size_t)di->data.integer;
                if (shape[i] && n > (size_t)INT64_C(100000000) / shape[i]) return NULL;
                n *= shape[i];
            }
        } else {
            if (!ne || ne->type != EXPR_INTEGER || ne->data.integer < 0) return NULL;
            n = (size_t)ne->data.integer;
        }
    }
    MlDist d; double* owned = NULL;
    if (!ml_read_dist(res->data.function.args[0], &d, &owned)) return NULL;
    size_t dim = d.dim ? d.dim : 1;

    double* z = malloc(sizeof(double) * dim);
    double* o = malloc(sizeof(double) * (n ? n : 1) * dim);
    Expr* out = NULL;
    if (z && o) {
        for (size_t i = 0; i < n; i++) ml_draw(&d, z, o + i * dim);
        if (argc == 1)
            out = (dim == 1) ? expr_new_real(o[0]) : ml_list_of_reals(o, dim);
        else
            out = (dim == 1) ? ml_list_of_reals(o, n) : ml_list_matrix(o, n, dim);
        if (out && rank > 1) {
            /* Fold the flat list of n draws into the requested shape, the last
             * dimension innermost. */
            for (size_t r = rank; r-- > 1; ) {
                size_t w = shape[r], groups = 1;
                for (size_t q = 0; q < r; q++) groups *= shape[q];
                Expr** gs = malloc(sizeof(Expr*) * (groups ? groups : 1));
                for (size_t g = 0; g < groups; g++) {
                    Expr** el = malloc(sizeof(Expr*) * (w ? w : 1));
                    for (size_t j = 0; j < w; j++)
                        el[j] = expr_copy(out->data.function.args[g * w + j]);
                    gs[g] = expr_new_function(expr_new_symbol(SYM_List), el, w);
                    free(el);
                }
                Expr* nxt = expr_new_function(expr_new_symbol(SYM_List), gs, groups);
                free(gs);
                expr_free(out);
                out = nxt;
            }
        }
    }
    free(z); free(o); free(owned);
    return out;
}

/* ---- LearnDistribution ---------------------------------------------------- */

/* LearnDistribution[data] fits a distribution and returns a LearnedDistribution.
 *
 * A learned distribution IS a fitted model, so it uses the trained-model
 * representation from src/ml/predict.h and prints elided. That is the deliberate
 * opposite of a SPECIFIED distribution like NormalDistribution[mu, sigma], which prints
 * in full because its parameters are what the user wrote rather than what was derived.
 * The two must not be unified, and a test pins each convention.
 *
 * Only "Multinormal" here; GaussianMixture is the next piece and will reuse
 * ml_gmm_fit / ml_gmm_logpdf, which were extracted for exactly that.
 */
/* The variance floor for a standalone mixture fit: the SQUARED MEDIAN NEAREST-NEIGHBOUR
 * DISTANCE BETWEEN DISTINCT POINTS.
 *
 * This is load-bearing, not defensive. A Gaussian mixture's likelihood is unbounded
 * above -- a component collapsing onto a single point drives its density, and hence the
 * likelihood, to infinity -- so with a floor set merely "small" the BIC search buys
 * arbitrarily many near-singular spikes. That was measured in the clustering path before
 * its floor existed: six components for eight points.
 *
 * The median nearest-neighbour distance is the standalone analogue of the median
 * spanning-tree edge weight that fc_gmm_ndim uses, and it says the same honest thing:
 * structure finer than the spacing between samples is not resolvable from this data. The
 * MEDIAN rather than the mean, so one tight pair cannot drag the floor to nearly zero
 * and reopen the same hole.
 *
 * DISTINCT points, because a repeated value is not a sample spacing of zero. Rounded or
 * otherwise discretised data -- iris petal lengths to 0.1 cm -- repeats most values, and
 * counting a duplicate's zero distance to its twin put the median at 0 once more than
 * half the points had one; the floor then fell to the 1e-300 fallback, and the BIC search
 * bought a component per stack of duplicates (nine for 45 petal lengths, several of them
 * a single point at weight 1/45). The spacing between distinct values is the data's real
 * resolution. That is also exactly the tie problem fc_scale_ndim meets, which it answers
 * by falling back to the mean spanning-tree edge.
 *
 * A second guard makes the floor scale-aware as well: it is at least a small fraction of
 * the average per-coordinate variance, so a single anomalously tight pair of distinct
 * points among otherwise coarse data cannot reopen the hole either. With fewer than two
 * distinct points there is no spread at all, and the floor is left to the caller.
 *
 * *n_distinct receives the number of distinct points, which also bounds k: a component
 * needs a location of its own. */
static double ml_nn_floor(const double* x, size_t n, size_t dim, size_t* n_distinct) {
    *n_distinct = n;
    if (n < 2) return 0.0;
    size_t* rep = malloc(sizeof(size_t) * n);           /* indices of distinct points */
    double* nn = malloc(sizeof(double) * n);
    if (!rep || !nn) { free(rep); free(nn); return 0.0; }
    size_t m = 0;
    for (size_t i = 0; i < n; i++) {                     /* O(n^2), as the search below */
        bool dup = false;
        for (size_t r = 0; r < m && !dup; r++)
            dup = memcmp(x + i * dim, x + rep[r] * dim, sizeof(double) * dim) == 0
                  || ml_sqdist(x + i * dim, x + rep[r] * dim, dim) == 0.0;
        if (!dup) rep[m++] = i;
    }
    *n_distinct = m;
    double med = 0.0;
    if (m >= 2) {
        for (size_t a = 0; a < m; a++) {
            double best = -1.0;
            for (size_t b = 0; b < m; b++) {
                if (b == a) continue;
                double d2 = ml_sqdist(x + rep[a] * dim, x + rep[b] * dim, dim);
                if (best < 0.0 || d2 < best) best = d2;
            }
            nn[a] = best;
        }
        for (size_t a = 1; a < m; a++) {                /* insertion sort, n is bounded */
            double v = nn[a]; size_t b = a;
            while (b > 0 && nn[b - 1] > v) { nn[b] = nn[b - 1]; b--; }
            nn[b] = v;
        }
        med = (m % 2) ? nn[m / 2] : 0.5 * (nn[m / 2 - 1] + nn[m / 2]);

        /* Scale guard: 1e-4 of the mean per-coordinate variance, i.e. a component
         * standard deviation of at least 1% of the data's. */
        double var = 0.0;
        for (size_t a = 0; a < dim; a++) {
            double mu = 0.0, ss = 0.0;
            for (size_t i = 0; i < n; i++) mu += x[i * dim + a];
            mu /= (double)n;
            for (size_t i = 0; i < n; i++) {
                double d = x[i * dim + a] - mu;
                ss += d * d;
            }
            var += ss / (double)n;
        }
        var /= (double)dim;
        if (med < 1e-4 * var) med = 1e-4 * var;
    }
    free(rep); free(nn);
    return med;                                          /* already a SQUARED distance */
}

/* True when no leaf of data (a List, or a List of Lists) is a number: the
 * outcomes are nominal (strings, symbols), which the default method learns as
 * a ContingencyTable instead of declining in the numeric reader. */
static bool ld_all_nominal(const Expr* data) {
    if (!data || data->type != EXPR_FUNCTION || data->data.function.head->type != EXPR_SYMBOL
        || data->data.function.head->data.symbol.name != SYM_List
        || data->data.function.arg_count == 0) return false;
    for (size_t i = 0; i < data->data.function.arg_count; i++) {
        const Expr* r = data->data.function.args[i];
        bool row = r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL
                   && r->data.function.head->data.symbol.name == SYM_List;
        size_t m = row ? r->data.function.arg_count : 1;
        for (size_t j = 0; j < m; j++) {
            const Expr* x = row ? r->data.function.args[j] : r;
            if (x->type != EXPR_STRING &&
                !(x->type == EXPR_SYMBOL && !expr_is_numeric_quantity(x))) return false;
        }
    }
    return true;
}

static Expr* builtin_learn_distribution(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    bool mixture = false;
    bool contingency = false;
    bool automatic = true;    /* no Method, or Method -> Automatic */
    if (argc == 2) {
        Expr* o = res->data.function.args[1];
        if (!o || o->type != EXPR_FUNCTION || o->data.function.arg_count != 2) return NULL;
        Expr* h = o->data.function.head;
        if (!h || h->type != EXPR_SYMBOL) return NULL;
        const char* hn = h->data.symbol.name;
        if (hn != SYM_Rule && hn != SYM_RuleDelayed) return NULL;
        Expr* lhs = o->data.function.args[0];
        Expr* rhs = o->data.function.args[1];
        if (!lhs || lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != SYM_Method)
            return NULL;
        if (rhs && rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) {
            /* the default method, as with no option */
        } else {
            if (!rhs || rhs->type != EXPR_STRING) return NULL;
            const char* mm = rhs->data.string;
            automatic = false;
            if (strcmp(mm, "Multinormal") == 0)          mixture = false;
            else if (strcmp(mm, "GaussianMixture") == 0) mixture = true;
            else if (strcmp(mm, "ContingencyTable") == 0) contingency = true;
            else return NULL;
        }
    }
    /* The default method detects nominal data. */
    if (automatic && ld_all_nominal(res->data.function.args[0])) contingency = true;

    if (contingency) {
        /* NOMINAL data: a probability per distinct outcome, which is what Wolfram's
         * "ContingencyTable" method is -- and in one dimension is a categorical
         * distribution.
         *
         * This branch comes BEFORE ml_read_data because ml_read_data reads NUMBERS, and
         * the whole point here is outcomes that are not numbers. "red", "blue", a symbol
         * or a pair {"a", "x"} are the intended input, and forcing them through a numeric
         * reader would decline exactly the data this method exists for.
         *
         * The label vocabulary does all the work: ml_labels_build already collects distinct
         * expressions in first-appearance order using expr_eq, which is precisely a
         * cross-tabulation over arbitrary outcomes. This is the third consumer of that
         * module and it needed no change to serve it.
         *
         * Counts are divided by n, so the probabilities are EMPIRICAL frequencies with no
         * smoothing. That is a deliberate choice, and the consequence is stated rather than
         * hidden: an outcome absent from the training data gets probability exactly 0, not a
         * small positive number. Add-one smoothing would need a claim about the size of the
         * outcome space, which for arbitrary expressions is unknowable -- there is no way to
         * know how many nominal values were possible but unseen. A hard 0 for "never
         * observed" is the honest answer, and PDF returns it as an exact 0. */
        Expr* data = res->data.function.args[0];
        if (!data || data->type != EXPR_FUNCTION) return NULL;
        Expr* dh = data->data.function.head;
        if (!dh || dh->type != EXPR_SYMBOL || dh->data.symbol.name != SYM_List) return NULL;
        size_t cn = data->data.function.arg_count;
        if (cn == 0) return NULL;

        /* The outcome arity: a list-valued outcome is that many nominal variables, anything
         * else is one. Ragged outcomes are a different table per row, which is not a
         * distribution, so they decline. */
        Expr* f0 = data->data.function.args[0];
        bool listy = f0 && f0->type == EXPR_FUNCTION && f0->data.function.head
                     && f0->data.function.head->type == EXPR_SYMBOL
                     && f0->data.function.head->data.symbol.name == SYM_List;
        size_t arity = listy ? f0->data.function.arg_count : 1;
        for (size_t i = 0; i < cn; i++) {
            Expr* r = data->data.function.args[i];
            if (!r) return NULL;
            bool rl = r->type == EXPR_FUNCTION && r->data.function.head
                      && r->data.function.head->type == EXPR_SYMBOL
                      && r->data.function.head->data.symbol.name == SYM_List;
            if (rl != listy) return NULL;
            if (rl && r->data.function.arg_count != arity) return NULL;
        }

        MlLabels vocab = {0, NULL};
        size_t* idx = malloc(sizeof(size_t) * cn);
        if (!idx) return NULL;
        if (!ml_labels_build(data->data.function.args, cn, &vocab, idx)) {
            free(idx); return NULL;
        }
        double* prob = calloc(vocab.count, sizeof(double));
        Expr* outc = NULL;
        if (prob) {
            for (size_t i = 0; i < cn; i++) prob[idx[i]] += 1.0;
            for (size_t k = 0; k < vocab.count; k++) prob[k] /= (double)cn;
            Expr** rows = malloc(sizeof(Expr*) * 2);
            if (rows) {
                rows[0] = ml_labels_to_list(&vocab);
                rows[1] = ml_list_of_reals(prob, vocab.count);
                Expr* pay = (rows[0] && rows[1])
                    ? expr_new_function(expr_new_symbol(SYM_List), rows, 2) : NULL;
                if (!pay) { expr_free(rows[0]); expr_free(rows[1]); }
                free(rows);
                if (pay) {
                    Expr* a4[4];
                    a4[0] = expr_new_string("ContingencyTable");
                    a4[1] = pay;
                    a4[2] = expr_new_integer((int64_t)arity);
                    a4[3] = expr_new_integer((int64_t)cn);
                    if (a4[0] && a4[1] && a4[2] && a4[3])
                        outc = expr_new_function(expr_new_symbol("LearnedDistribution"), a4, 4);
                    else { expr_free(a4[0]); expr_free(a4[1]);
                           expr_free(a4[2]); expr_free(a4[3]); }
                }
            }
        }
        free(prob); free(idx); ml_labels_free(&vocab);
        return outc;
    }

    size_t n, dim; double* x = NULL; bool vec = false;
    if (!ml_read_data(res->data.function.args[0], &n, &dim, &x, &vec)) return NULL;
    /* A flat list is n observations of ONE variable, which is a perfectly good
     * univariate normal -- so unlike PrincipalComponents this does not decline it. */
    if (n < 2) { free(x); return NULL; }        /* one point has no dispersion to fit */

    if (mixture) {
        /* BIC over k = 1..k_max, exactly as the clustering path does. k_max is bounded by
         * the PARAMETER COUNT rather than the point count: a full covariance costs
         * dim*(dim+1)/2 numbers per component, so a component needs more points than
         * dimensions before it means anything. */
        size_t n_distinct = n;
        double vfloor = ml_nn_floor(x, n, dim, &n_distinct);
        if (!(vfloor > 0.0)) vfloor = 1e-300;           /* all points identical: k = 1 */
        size_t kmax = n / (dim + 1);
        if (kmax > n_distinct) kmax = n_distinct;       /* as FindClusters caps it */
        if (kmax > 10) kmax = 10;
        if (kmax < 1) kmax = 1;

        MlGmm* best = NULL; double best_bic = INFINITY;
        for (size_t kk = 1; kk <= kmax; kk++) {
            MlGmm* g = ml_gmm_fit(x, n, dim, kk, vfloor, NULL);
            if (!g) continue;
            double bic = ml_gmm_param_count(kk, dim) * log((double)n) - 2.0 * g->loglik;
            if (bic < best_bic) {
                best_bic = bic;
                if (best) ml_gmm_free(best);
                best = g;
            } else {
                ml_gmm_free(g);
            }
        }
        free(x);
        if (!best) return NULL;

        size_t k = best->k;
        Expr** rows = malloc(sizeof(Expr*) * (1 + k * (1 + dim)));
        Expr* outm = NULL;
        if (rows) {
            rows[0] = ml_list_of_reals(best->w, k);
            for (size_t j = 0; j < k; j++) {
                size_t base = 1 + j * (1 + dim);
                rows[base] = ml_list_of_reals(best->mu + j * dim, dim);
                for (size_t a = 0; a < dim; a++)
                    rows[base + 1 + a] =
                        ml_list_of_reals(best->cov + j * dim * dim + a * dim, dim);
            }
            Expr* pay = expr_new_function(expr_new_symbol(SYM_List), rows,
                                          1 + k * (1 + dim));
            free(rows);
            if (pay) {
                Expr* a4[4];
                a4[0] = expr_new_string("GaussianMixture");
                a4[1] = pay;
                a4[2] = expr_new_integer((int64_t)dim);
                a4[3] = expr_new_integer((int64_t)k);
                if (a4[0] && a4[1] && a4[2] && a4[3])
                    outm = expr_new_function(expr_new_symbol("LearnedDistribution"), a4, 4);
                else { expr_free(a4[0]); expr_free(a4[1]); expr_free(a4[2]); expr_free(a4[3]); }
            }
        }
        ml_gmm_free(best);
        return outm;
    }

    double* mean = malloc(sizeof(double) * dim);
    double* cov  = malloc(sizeof(double) * dim * dim);
    double* chk  = malloc(sizeof(double) * dim * dim);
    if (!mean || !cov || !chk) { free(mean); free(cov); free(chk); free(x); return NULL; }
    ml_column_mean(x, n, dim, mean);
    /* Sample covariance with the n-1 divisor, matching ml_column_sd and Variance, so a
     * one-variable fit agrees with StandardDeviation squared. */
    double den = (double)(n - 1);
    for (size_t a = 0; a < dim; a++)
        for (size_t b = 0; b <= a; b++) {
            double sacc = 0.0;
            for (size_t i = 0; i < n; i++)
                sacc += (x[i * dim + a] - mean[a]) * (x[i * dim + b] - mean[b]);
            cov[a * dim + b] = cov[b * dim + a] = sacc / den;
        }

    /* Refuse a singular covariance rather than returning an object with no density.
     * Collinear features or fewer observations than dimensions land here, and a
     * pseudo-density would be a fiction. */
    Expr* out = NULL;
    /* A sample covariance has rank at most n - 1, so n <= dim is singular whatever the
     * arithmetic says. Beyond that, an exact Cholesky test (pivot > 0) is not a rank test: for
     * two points in the plane the second pivot is 0.5 - 0.5 evaluated in floating point, a
     * roundoff residue of 1e-16 that "passed" and produced a density spike of 2e7 at the
     * mean. A pivot is accepted only when it carries a non-negligible share of its own
     * variance, L_aa^2 > 64 eps Sigma_aa -- far below any genuine conditional variance. */
    bool spd = n > dim && ml_chol(cov, chk, dim);
    for (size_t a = 0; spd && a < dim; a++)
        if (!(chk[a * dim + a] * chk[a * dim + a] > 64.0 * DBL_EPSILON * cov[a * dim + a]))
            spd = false;
    if (spd) {
        Expr** rows = malloc(sizeof(Expr*) * (dim + 1));
        if (rows) {
            rows[0] = ml_list_of_reals(mean, dim);
            for (size_t a = 0; a < dim; a++)
                rows[a + 1] = ml_list_of_reals(cov + a * dim, dim);
            Expr* pay = expr_new_function(expr_new_symbol(SYM_List), rows, dim + 1);
            free(rows);
            if (pay) {
                Expr* a4[4];
                a4[0] = expr_new_string("Multinormal");
                a4[1] = pay;
                a4[2] = expr_new_integer((int64_t)dim);
                a4[3] = expr_new_integer(0);
                if (a4[0] && a4[1] && a4[2] && a4[3])
                    out = expr_new_function(expr_new_symbol("LearnedDistribution"), a4, 4);
                else { expr_free(a4[0]); expr_free(a4[1]); expr_free(a4[2]); expr_free(a4[3]); }
            }
        }
    }
    free(mean); free(cov); free(chk); free(x);
    return out;
}

/* SmoothKernelDistribution[data] -- a kernel density estimate.
 *
 * Bandwidth: the multivariate NORMAL-REFERENCE rule,
 *   h_a = sigma_a * (4 / ((dim + 2) n))^(1/(dim + 4))
 * which in one dimension is exactly Silverman's 1.06 sigma n^(-1/5) --
 * (4/3)^(1/5) = 1.0592. Naming it as the multivariate rule rather than as Silverman's is
 * the honest description, since it is applied per dimension for any dim.
 *
 * It is a NORMAL-reference rule, so it is the right default only insofar as the data is
 * not wildly non-normal; on strongly multimodal data it oversmooths, which is a known
 * property of the rule rather than a defect here. A user-supplied bandwidth is the
 * remedy, and is accepted as a second argument.
 */
static Expr* builtin_smooth_kernel(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    size_t n, dim; double* x = NULL; bool vec = false;
    if (!ml_read_data(res->data.function.args[0], &n, &dim, &x, &vec)) return NULL;
    if (n < 2) { free(x); return NULL; }    /* one point has no scale to estimate */

    double* band = malloc(sizeof(double) * dim);
    double* mean = malloc(sizeof(double) * dim);
    if (!band || !mean) { free(band); free(mean); free(x); return NULL; }

    bool ok = true;
    if (argc == 2) {
        /* An explicit bandwidth: one number for every dimension, or one per dimension. */
        double h = 0.0, im = 0.0;
        Expr* be = res->data.function.args[1];
        if (na_read_scalar(be, &h, &im) && im == 0.0 && h > 0.0) {
            for (size_t a = 0; a < dim; a++) band[a] = h;
        } else {
            size_t bn, bd; double* bb = NULL; bool bvec = false;
            if (ml_read_data(be, &bn, &bd, &bb, &bvec) && bvec && bn == dim) {
                for (size_t a = 0; a < dim; a++) {
                    band[a] = bb[a];
                    if (!(band[a] > 0.0)) ok = false;
                }
            } else ok = false;
            free(bb);
        }
    } else {
        ml_column_mean(x, n, dim, mean);
        ml_column_sd(x, n, dim, mean, band);          /* n-1 divisor, matching Variance */
        double expo = 1.0 / ((double)dim + 4.0);
        double factor = pow(4.0 / (((double)dim + 2.0) * (double)n), expo);
        for (size_t a = 0; a < dim; a++) {
            band[a] *= factor;
            /* A constant column has zero spread, so the normal-reference rule gives a
             * zero bandwidth and there is no density. Declining beats dividing by zero
             * and returning infinities. */
            if (!(band[a] > 0.0)) ok = false;
        }
    }

    Expr* out = NULL;
    if (ok) {
        Expr** rows = malloc(sizeof(Expr*) * (n + 1));
        if (rows) {
            rows[0] = ml_list_of_reals(band, dim);
            for (size_t i = 0; i < n; i++) rows[i + 1] = ml_list_of_reals(x + i * dim, dim);
            Expr* pay = expr_new_function(expr_new_symbol(SYM_List), rows, n + 1);
            free(rows);
            if (pay) {
                Expr* a4[4];
                a4[0] = expr_new_string("SmoothKernel");
                a4[1] = pay;
                a4[2] = expr_new_integer((int64_t)dim);
                a4[3] = expr_new_integer((int64_t)n);
                if (a4[0] && a4[1] && a4[2] && a4[3])
                    out = expr_new_function(expr_new_symbol("LearnedDistribution"), a4, 4);
                else { expr_free(a4[0]); expr_free(a4[1]); expr_free(a4[2]); expr_free(a4[3]); }
            }
        }
    }
    free(band); free(mean); free(x);
    return out;
}

void ml_dist_init(void) {
    symtab_add_builtin("SmoothKernelDistribution", builtin_smooth_kernel);
    symtab_get_def("SmoothKernelDistribution")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("SmoothKernelDistribution",
        "SmoothKernelDistribution[data] gives a kernel density estimate as a "
        "LearnedDistribution, usable with PDF. The kernel is a product Gaussian with a "
        "per-dimension bandwidth from the multivariate normal-reference rule, which in "
        "one dimension is Silverman's 1.06 sigma n^(-1/5). "
        "SmoothKernelDistribution[data, h] sets the bandwidth explicitly, as one number "
        "or one per dimension. Being a normal-reference rule the default oversmooths "
        "strongly multimodal data -- a known property of the rule, and the reason the "
        "explicit form exists. A constant column has no scale and returns unevaluated.");

    symtab_add_builtin("LearnDistribution", builtin_learn_distribution);
    symtab_get_def("LearnDistribution")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LearnDistribution",
        "LearnDistribution[data] fits a distribution to data and returns a "
        "LearnedDistribution, usable with PDF. Method -> \"Multinormal\" is the default; "
        "Method -> \"GaussianMixture\" fits a mixture, choosing the component count by "
        "BIC. Multinormal fits a mean vector and a sample "
        "covariance (n-1 divisor, matching Variance). Rows are observations and columns "
        "are variables; a flat list is n observations of one variable. A singular "
        "covariance -- collinear columns, or fewer observations than dimensions -- "
        "returns unevaluated, because no density exists rather than because of an "
        "error. Method -> \"ContingencyTable\" is for NOMINAL data instead of numeric: it "
        "stores a probability per distinct outcome, which in one dimension is a "
        "categorical distribution. Outcomes may be any expressions -- strings, symbols, or "
        "equal-length lists of them -- compared structurally, and are kept in "
        "first-appearance order. Probabilities are empirical frequencies with no "
        "smoothing, so PDF of an outcome never observed is exactly 0; smoothing would "
        "require knowing how many outcomes were possible but unseen, which for arbitrary "
        "expressions is unknowable. Ragged outcomes decline.");

    symtab_get_def("LearnedDistribution")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LearnedDistribution",
        "LearnedDistribution[method, parameters, dimension, extra] is the fitted "
        "distribution LearnDistribution returns. Use it with PDF. Unlike a SPECIFIED "
        "distribution such as NormalDistribution[mu, sigma], it prints elided, because "
        "its parameters are derived rather than user-supplied.");

    symtab_add_builtin("RandomVariate", builtin_randomvariate);
    symtab_get_def("RandomVariate")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("RandomVariate",
        "RandomVariate[dist] draws one value from dist; RandomVariate[dist, n] draws a "
        "list of n, and RandomVariate[dist, {n1, n2, ...}] an n1 x n2 x ... array. "
        "Parameters may be exact numeric constants (Pi, Sqrt[2]). Supports NormalDistribution[mu, sigma] and "
        "UniformDistribution[{lo, hi}], each also usable with no arguments for the "
        "standard case, and every LearnedDistribution LearnDistribution or "
        "SmoothKernelDistribution returns (Multinormal, GaussianMixture, SmoothKernel); a "
        "multivariate one gives a vector per draw. Draws come from the same stream as RandomReal, so SeedRandom "
        "makes them reproducible. A non-positive standard deviation, or an inverted "
        "range, returns unevaluated rather than producing NaNs.");

    symtab_add_builtin("PDF", builtin_pdf);
    symtab_get_def("PDF")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PDF",
        "PDF[dist, x] gives the probability density of dist at x, and threads over a "
        "list of x. Supports NormalDistribution and UniformDistribution (on an interval "
        "{a, b} or a box {{a1, b1}, ...}) -- exactly for exact input, at the input "
        "precision for arbitrary-precision input, and as a formula for symbolic input -- "
        "and every LearnedDistribution; a univariate learned distribution takes a scalar "
        "point. PDF[dist][x] is PDF[dist, x]. A parameter set that is no distribution "
        "(NormalDistribution[m, s] with s <= 0, decreasing bounds) gives a message.");
    /* PDF[dist][x] and CDF[dist][x]: the point is the second argument */
    opform_register("PDF", 1, 1, 1);
    opform_register("CDF", 1, 1, 1);

    symtab_add_builtin("InverseCDF", builtin_inverse_cdf);
    symtab_get_def("InverseCDF")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("InverseCDF",
        "InverseCDF[dist, q] gives the inverse of the CDF of dist at the probability q: "
        "for the continuous univariate distributions here, Quantile[dist, q] "
        "(NormalDistribution, UniformDistribution and a univariate SmoothKernel "
        "LearnedDistribution).");

    symtab_add_builtin("CDF", builtin_cdf);
    symtab_get_def("CDF")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CDF",
        "CDF[dist, x] gives the cumulative distribution function of dist at x, and threads "
        "over a list of x. NormalDistribution[m, s] gives Erfc[(m - x)/(Sqrt[2] s)]/2 and "
        "UniformDistribution[{a, b}] the piecewise-linear ramp, exactly for exact input; a "
        "univariate SmoothKernelDistribution gives the mean of its kernel CDFs.");

    symtab_get_def("NormalDistribution")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("NormalDistribution",
        "NormalDistribution[mu, sigma] represents a normal distribution; "
        "NormalDistribution[] is the standard normal. Unlike a fitted model it prints "
        "its parameters in full, because they are what the user specified rather than "
        "an implementation detail.");

    symtab_get_def("UniformDistribution")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("UniformDistribution",
        "UniformDistribution[{lo, hi}] represents a continuous uniform distribution; "
        "UniformDistribution[] is uniform on {0, 1}.");
}
