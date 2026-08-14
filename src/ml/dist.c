/* dist.c -- distribution objects, RandomVariate and PDF. See dist.h. */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "numarray.h"
#include "random.h"
#include "mlutil.h"
#include "pca.h"       /* ml_column_mean */
#include "gmm.h"       /* ml_gmm_fit / ml_gmm_param_count -- extracted in iteration 9 for this */
#include "dist.h"

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

static bool ml_read_dist(Expr* e, MlDist* d, double** owned) {
    *owned = NULL;
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    const char* hn = h->data.symbol.name;
    size_t argc = e->data.function.arg_count;
    double im = 0.0;

    if (strcmp(hn, "LearnedDistribution") == 0) return ml_read_learned(e, d, owned);

    if (strcmp(hn, "NormalDistribution") == 0) {
        d->kind = ML_D_NORMAL; d->dim = 1;
        if (argc == 0) { d->a = 0.0; d->b = 1.0; return true; }   /* the standard normal */
        if (argc != 2) return false;
        if (!na_read_scalar(e->data.function.args[0], &d->a, &im) || im != 0.0) return false;
        if (!na_read_scalar(e->data.function.args[1], &d->b, &im) || im != 0.0) return false;
        /* A non-positive standard deviation is not a distribution. Declining beats
         * returning NaNs that propagate silently through a whole sample. */
        return d->b > 0.0;
    }
    if (strcmp(hn, "UniformDistribution") == 0) {
        d->kind = ML_D_UNIFORM; d->dim = 1;
        if (argc == 0) { d->a = 0.0; d->b = 1.0; return true; }
        if (argc != 1) return false;
        Expr* r = e->data.function.args[0];
        if (!r || r->type != EXPR_FUNCTION || r->data.function.arg_count != 2) return false;
        if (!na_read_scalar(r->data.function.args[0], &d->a, &im) || im != 0.0) return false;
        if (!na_read_scalar(r->data.function.args[1], &d->b, &im) || im != 0.0) return false;
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

static Expr* builtin_pdf(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
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
        if (ml_read_data(xe, &n, &dm, &px, &vec)) {
            double p = 0.0;
            bool (*pdf)(const MlDist*, const double*, double*) =
                (d.kind == ML_D_MIXTURE) ? ml_mixture_pdf :
                (d.kind == ML_D_KDE)     ? ml_kde_pdf : ml_multinormal_pdf;
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

/* ---- RandomVariate ------------------------------------------------------ */

static double ml_draw(const MlDist* d) {
    switch (d->kind) {
        case ML_D_NORMAL:  return d->a + d->b * ml_normal_deviate();
        case ML_D_UNIFORM: return d->a + (d->b - d->a) * random_uniform_01();
        default:           return 0.0;
    }
}

static Expr* builtin_randomvariate(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    MlDist d; double* owned = NULL;
    if (!ml_read_dist(res->data.function.args[0], &d, &owned)) return NULL;

    if (argc == 1) {
        double v = ml_draw(&d);
        free(owned);
        return expr_new_real(v);
    }
    Expr* ne = res->data.function.args[1];
    if (!ne || ne->type != EXPR_INTEGER || ne->data.integer < 0) { free(owned); return NULL; }
    size_t n = (size_t)ne->data.integer;
    double* o = malloc(sizeof(double) * (n ? n : 1));
    if (!o) { free(owned); return NULL; }
    for (size_t i = 0; i < n; i++) o[i] = ml_draw(&d);
    Expr* out = ml_list_of_reals(o, n);
    free(o); free(owned);
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
 * DISTANCE.
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
 * and reopen the same hole. */
static double ml_nn_floor(const double* x, size_t n, size_t dim) {
    if (n < 2) return 0.0;
    double* nn = malloc(sizeof(double) * n);
    if (!nn) return 0.0;
    for (size_t i = 0; i < n; i++) {
        double best = -1.0;
        for (size_t j = 0; j < n; j++) {
            if (j == i) continue;
            double d2 = ml_sqdist(x + i * dim, x + j * dim, dim);
            if (best < 0.0 || d2 < best) best = d2;
        }
        nn[i] = (best > 0.0) ? best : 0.0;
    }
    for (size_t a = 1; a < n; a++) {                    /* insertion sort, n is bounded */
        double v = nn[a]; size_t b = a;
        while (b > 0 && nn[b - 1] > v) { nn[b] = nn[b - 1]; b--; }
        nn[b] = v;
    }
    double med = (n % 2) ? nn[n / 2] : 0.5 * (nn[n / 2 - 1] + nn[n / 2]);
    free(nn);
    return med;                                          /* already a SQUARED distance */
}

static Expr* builtin_learn_distribution(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    bool mixture = false;
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
        if (!rhs || rhs->type != EXPR_STRING) return NULL;
        const char* mm = rhs->data.string;
        if (strcmp(mm, "Multinormal") == 0)          mixture = false;
        else if (strcmp(mm, "GaussianMixture") == 0) mixture = true;
        else return NULL;
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
        double vfloor = ml_nn_floor(x, n, dim);
        if (!(vfloor > 0.0)) vfloor = 1e-300;
        size_t kmax = n / (dim + 1);
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
    if (ml_chol(cov, chk, dim)) {
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
        "error.");

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
        "list of n. Supports NormalDistribution[mu, sigma] and "
        "UniformDistribution[{lo, hi}], each also usable with no arguments for the "
        "standard case. Draws come from the same stream as RandomReal, so SeedRandom "
        "makes them reproducible. A non-positive standard deviation, or an inverted "
        "range, returns unevaluated rather than producing NaNs.");

    symtab_add_builtin("PDF", builtin_pdf);
    symtab_get_def("PDF")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PDF",
        "PDF[dist, x] gives the probability density of dist at x, and threads over a "
        "list of x. Supports NormalDistribution and UniformDistribution.");

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
