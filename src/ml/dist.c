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

static bool ml_read_dist(Expr* e, MlDist* d, double** owned) {
    *owned = NULL;
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    const char* hn = h->data.symbol.name;
    size_t argc = e->data.function.arg_count;
    double im = 0.0;

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

static Expr* builtin_pdf(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    MlDist d; double* owned = NULL;
    if (!ml_read_dist(res->data.function.args[0], &d, &owned)) return NULL;
    Expr* xe = res->data.function.args[1];
    double x = 0.0, im = 0.0;

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

void ml_dist_init(void) {
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
