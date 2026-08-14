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
typedef enum { ML_D_NONE, ML_D_NORMAL, ML_D_UNIFORM, ML_D_MULTINORMAL } MlDistKind;

typedef struct {
    MlDistKind kind;
    double a, b;          /* Normal: mu, sigma. Uniform: lo, hi. */
    size_t dim;           /* Multinormal */
    const double* mu;     /* borrowed */
    const double* chol;   /* borrowed, lower factor */
} MlDist;

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

static Expr* builtin_pdf(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    MlDist d; double* owned = NULL;
    if (!ml_read_dist(res->data.function.args[0], &d, &owned)) return NULL;
    Expr* xe = res->data.function.args[1];
    double x = 0.0, im = 0.0;

    if (d.kind == ML_D_MULTINORMAL) {
        /* The argument is a POINT, so a list here is one observation rather than many.
         * That is the opposite reading from the scalar case below, where a list is many
         * points -- and it has to be, because a multinormal's argument is itself a
         * list. Getting this backwards would silently treat each coordinate as a
         * separate observation. */
        size_t n, dm; double* px = NULL; bool vec = false;
        Expr* out = NULL;
        if (ml_read_data(xe, &n, &dm, &px, &vec)) {
            double p = 0.0;
            if (vec && n == d.dim && ml_multinormal_pdf(&d, px, &p))
                out = expr_new_real(p);
            else if (!vec && dm == d.dim) {         /* a matrix: many points */
                double* o = malloc(sizeof(double) * n);
                bool ok = o != NULL;
                for (size_t i = 0; ok && i < n; i++)
                    ok = ml_multinormal_pdf(&d, px + i * dm, &o[i]);
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
static Expr* builtin_learn_distribution(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
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
        if (strcmp(rhs->data.string, "Multinormal") != 0) return NULL;
    }

    size_t n, dim; double* x = NULL; bool vec = false;
    if (!ml_read_data(res->data.function.args[0], &n, &dim, &x, &vec)) return NULL;
    /* A flat list is n observations of ONE variable, which is a perfectly good
     * univariate normal -- so unlike PrincipalComponents this does not decline it. */
    if (n < 2) { free(x); return NULL; }        /* one point has no dispersion to fit */

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

void ml_dist_init(void) {
    symtab_add_builtin("LearnDistribution", builtin_learn_distribution);
    symtab_get_def("LearnDistribution")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LearnDistribution",
        "LearnDistribution[data] fits a distribution to data and returns a "
        "LearnedDistribution, usable with PDF. Method -> \"Multinormal\" is the only "
        "method implemented and is the default: it fits a mean vector and a sample "
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
