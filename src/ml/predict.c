/* predict.c -- Predict, LinearModelFit, and applying a fitted model.
 *
 * See predict.h for the trained-model representation and why it is a plain
 * EXPR_FUNCTION rather than a new node type.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "numarray.h"
#include "mlutil.h"
#include "pca.h"        /* ml_pca, ml_column_mean -- shared with the reducer */
#include "classify.h"
#include "predict.h"
#include "eval.h"       /* evaluate: basis-function values for LinearModelFit[data, funs, vars] */

/* Solve the normal equations (A'A) c = A'y where A is [1 | x], by Gaussian
 * elimination with partial pivoting.
 *
 * NOT LAPACK's dgels, and the reason is not laziness. dgels is the numerically better
 * route -- a QR of A avoids squaring the condition number -- but A'A here is
 * (dim + 1) x (dim + 1), a handful of rows for any realistic feature count, and the
 * normal equations are what make the SINGULAR case detectable: a zero pivot says
 * "these features are collinear" directly, where dgels' rank handling would have to be
 * interrogated. Since declining on a singular fit is a deliberate behaviour here (see
 * predict.h), the more legible failure mode wins at this size. A feature count large
 * enough for the conditioning to matter would be a reason to switch, and the call site
 * is one function.
 */
bool ml_ols(const double* x, const double* y, size_t n, size_t dim, double* coef) {
    if (!x || !y || !coef || n == 0 || dim == 0) return false;
    if (n < dim + 1) return false;            /* underdetermined: no unique answer */

    /* CENTRE FIRST. With an intercept the slopes are the least-squares solution of the
     * centred problem (x - xbar) b = y - ybar, and the intercept is ybar - xbar . b. Forming
     * the normal equations on the raw columns instead put sum(1) = n next to sum(x^2) ~ n x^2
     * in one matrix: for x ~ 1e8 that is a condition number ~1e17, and an exactly collinear
     * data set came back with slope 1.68 instead of 2. Centring removes the offset, which is
     * what dominated the conditioning, and costs one pass. */
    size_t p = dim;
    double* xbar = calloc(dim, sizeof(double));
    double* a = calloc(p * (p + 1), sizeof(double));   /* augmented [Xc'Xc | Xc'yc] */
    if (!xbar || !a) { free(xbar); free(a); return false; }
    double ybar = 0.0;
    for (size_t i = 0; i < n; i++) {
        ybar += y[i];
        for (size_t j = 0; j < dim; j++) xbar[j] += x[i * dim + j];
    }
    ybar /= (double)n;
    for (size_t j = 0; j < dim; j++) xbar[j] /= (double)n;

    for (size_t i = 0; i < n; i++) {
        double yc = y[i] - ybar;
        for (size_t r = 0; r < p; r++) {
            double ar = x[i * dim + r] - xbar[r];
            for (size_t c = 0; c < p; c++)
                a[r * (p + 1) + c] += ar * (x[i * dim + c] - xbar[c]);
            a[r * (p + 1) + p] += ar * yc;
        }
    }

    for (size_t c = 0; c < p; c++) {
        size_t piv = c; double best = fabs(a[c * (p + 1) + c]);
        for (size_t r = c + 1; r < p; r++) {
            double v = fabs(a[r * (p + 1) + c]);
            if (v > best) { best = v; piv = r; }
        }
        /* A zero pivot is a collinear feature set. Scaled against the largest
         * diagonal so the test is not a bare absolute tolerance on unscaled data. */
        if (!(best > 1e-12)) { free(a); free(xbar); return false; }
        if (piv != c)
            for (size_t k = 0; k <= p; k++) {
                double t = a[c * (p + 1) + k];
                a[c * (p + 1) + k] = a[piv * (p + 1) + k];
                a[piv * (p + 1) + k] = t;
            }
        for (size_t r = 0; r < p; r++) {
            if (r == c) continue;
            double f = a[r * (p + 1) + c] / a[c * (p + 1) + c];
            if (f == 0.0) continue;
            for (size_t k = c; k <= p; k++)
                a[r * (p + 1) + k] -= f * a[c * (p + 1) + k];
        }
    }
    double icpt = ybar;
    for (size_t r = 0; r < p; r++) {
        coef[r + 1] = a[r * (p + 1) + p] / a[r * (p + 1) + r];
        icpt -= coef[r + 1] * xbar[r];
    }
    coef[0] = icpt;
    free(a); free(xbar);
    return true;
}

/* ------------------------------------------------------------------------- */

/* Split `data` into an n x dim feature matrix and an n response vector.
 *
 * Wolfram's Predict takes a list of rules, {features -> value, ...}. Also accepted
 * here is a plain matrix whose LAST column is the response, which is the shape data
 * usually arrives in and saves the caller a Thread. Both are read, and which one was
 * given is not remembered -- the fitted model is the same either way. */
static bool ml_read_training(Expr* data, size_t* n_out, size_t* dim_out,
                             double** x_out, double** y_out) {
    if (!data || data->type != EXPR_FUNCTION) return false;
    Expr* h = data->data.function.head;
    if (!h || h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_List) return false;
    size_t n = data->data.function.arg_count;
    if (n == 0) return false;

    Expr* first = data->data.function.args[0];
    /* arg_count == 2: the width is read off the first rule's lhs
     * (Predict[{Rule[]}] read past an empty node, as Classify did). */
    bool rules = first && first->type == EXPR_FUNCTION
              && first->data.function.arg_count == 2
              && first->data.function.head->type == EXPR_SYMBOL
              && (first->data.function.head->data.symbol.name == SYM_Rule
                  || first->data.function.head->data.symbol.name == SYM_RuleDelayed);

    if (rules) {
        /* Feature width comes from the first row; a ragged row then fails to load and
         * the whole call declines, rather than being padded or truncated. */
        Expr* lhs0 = first->data.function.args[0];
        size_t dim = 1;
        if (lhs0 && lhs0->type == EXPR_FUNCTION
            && lhs0->data.function.head->type == EXPR_SYMBOL
            && lhs0->data.function.head->data.symbol.name == SYM_List)
            dim = lhs0->data.function.arg_count;
        if (dim == 0) return false;

        double* x = malloc(sizeof(double) * n * dim);
        double* y = malloc(sizeof(double) * n);
        if (!x || !y) { free(x); free(y); return false; }
        for (size_t i = 0; i < n; i++) {
            Expr* r = data->data.function.args[i];
            if (!r || r->type != EXPR_FUNCTION || r->data.function.arg_count != 2
                || r->data.function.head->type != EXPR_SYMBOL
                || (r->data.function.head->data.symbol.name != SYM_Rule
                    && r->data.function.head->data.symbol.name != SYM_RuleDelayed)) {
                free(x); free(y); return false;
            }
            Expr* lhs = r->data.function.args[0];
            Expr* rhs = r->data.function.args[1];
            double im = 0.0;
            if (dim == 1 && !(lhs && lhs->type == EXPR_FUNCTION)) {
                if (!na_read_scalar(lhs, &x[i * dim], &im) || im != 0.0) {
                    free(x); free(y); return false;
                }
            } else {
                if (!lhs || lhs->type != EXPR_FUNCTION
                    || lhs->data.function.arg_count != dim) { free(x); free(y); return false; }
                for (size_t j = 0; j < dim; j++)
                    if (!na_read_scalar(lhs->data.function.args[j], &x[i * dim + j], &im)
                        || im != 0.0) { free(x); free(y); return false; }
            }
            if (!na_read_scalar(rhs, &y[i], &im) || im != 0.0) {
                free(x); free(y); return false;
            }
        }
        *n_out = n; *dim_out = dim; *x_out = x; *y_out = y;
        return true;
    }

    /* Matrix form: last column is the response. */
    int rr = 0, cc = 0; double* buf = NULL;
    if (!na_load_matrix(data, 0, 0, &rr, &cc, &buf)) return false;
    if (rr <= 0 || cc < 2) { free(buf); return false; }
    size_t rows = (size_t)rr, cols = (size_t)cc, dim = cols - 1;
    double* x = malloc(sizeof(double) * rows * dim);
    double* y = malloc(sizeof(double) * rows);
    if (!x || !y) { free(x); free(y); free(buf); return false; }
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < dim; j++) x[i * dim + j] = buf[i * cols + j];
        y[i] = buf[i * cols + dim];
    }
    free(buf);
    *n_out = rows; *dim_out = dim; *x_out = x; *y_out = y;
    return true;
}

/* Assemble a fitted model. `params` is already-built and OWNED by this call.
 *
 * The parameter slot holds whatever the method needs and nothing more: a coefficient
 * List for LinearRegression, the training matrix for NearestNeighbors. That the two
 * shapes differ is not a wart -- it is why the representation is positional with a
 * method tag rather than a fixed schema, and ml_model_apply dispatches on the tag. */
static Expr* ml_make_model(const char* method, Expr* params, size_t dim, size_t extra) {
    Expr* args[4];
    args[0] = expr_new_string(method);
    args[1] = params;
    args[2] = expr_new_integer((int64_t)dim);
    args[3] = expr_new_integer((int64_t)extra);
    if (!args[0] || !args[1] || !args[2] || !args[3]) {
        expr_free(args[0]); expr_free(args[1]); expr_free(args[2]); expr_free(args[3]);
        return NULL;
    }
    return expr_new_function(expr_new_symbol("PredictorFunction"), args, 4);
}

static Expr* ml_make_predictor(const char* method, const double* coef, size_t dim) {
    Expr* c = ml_list_of_reals(coef, dim + 1);
    if (!c) return NULL;
    return ml_make_model(method, c, dim, 0);
}

static Expr* builtin_predict(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    bool knn = false;

    size_t kopt = 0;                  /* 0 = use the default */
    if (argc == 2) {
        /* Method -> m, or Method -> {m, subopt -> value} for NearestNeighbors'
         * NeighborsNumber. Anything unrecognised declines rather than silently
         * linear-regressing, so an unsupported method or a mistyped sub-option is
         * visible instead of quietly ignored. */
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
        /* The value is either a method string or {method, subopt -> value, ...}. */
        Expr* mstr = rhs;
        if (rhs && rhs->type == EXPR_FUNCTION
            && rhs->data.function.head->type == EXPR_SYMBOL
            && rhs->data.function.head->data.symbol.name == SYM_List
            && rhs->data.function.arg_count >= 1) {
            mstr = rhs->data.function.args[0];
            for (size_t i = 1; i < rhs->data.function.arg_count; i++) {
                Expr* so = rhs->data.function.args[i];
                if (!so || so->type != EXPR_FUNCTION || so->data.function.arg_count != 2
                    || so->data.function.head->type != EXPR_SYMBOL
                    || (so->data.function.head->data.symbol.name != SYM_Rule
                        && so->data.function.head->data.symbol.name != SYM_RuleDelayed))
                    return NULL;
                Expr* sk = so->data.function.args[0];
                Expr* sv = so->data.function.args[1];
                if (!sk || sk->type != EXPR_STRING) return NULL;
                if (strcmp(sk->data.string, "NeighborsNumber") != 0) return NULL;
                if (!sv || sv->type != EXPR_INTEGER || sv->data.integer <= 0) return NULL;
                kopt = (size_t)sv->data.integer;
            }
        }
        if (!mstr || mstr->type != EXPR_STRING) return NULL;
        const char* mm = mstr->data.string;
        if (strcmp(mm, "LinearRegression") == 0)      knn = false;
        else if (strcmp(mm, "NearestNeighbors") == 0) knn = true;
        else return NULL;
        /* A neighbour count is meaningless for a regression, so it is refused rather
         * than accepted and ignored. */
        if (kopt && !knn) return NULL;
    }

    size_t n, dim; double* x = NULL; double* y = NULL;
    if (!ml_read_training(res->data.function.args[0], &n, &dim, &x, &y)) return NULL;

    Expr* out = NULL;
    if (knn) {
        /* A nearest-neighbour predictor has no fitted parameters -- the training set IS
         * the model, which is why this method costs nothing to fit and everything to
         * apply. Stored as one matrix with the response appended to each row, so the
         * features and their answers cannot drift apart.
         *
         * k defaults to 3: enough to average away a single noisy response, few enough
         * to stay local. Clamped to n, since asking for more neighbours than there are
         * points is a request the data cannot fill. */
        size_t k = kopt ? kopt : 3;
        if (k > n) k = n;
        double* joint = malloc(sizeof(double) * n * (dim + 1));
        if (joint) {
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < dim; j++) joint[i * (dim + 1) + j] = x[i * dim + j];
                joint[i * (dim + 1) + dim] = y[i];
            }
            Expr* tm = ml_list_matrix(joint, n, dim + 1);
            free(joint);
            if (tm) out = ml_make_model("NearestNeighbors", tm, dim, k);
        }
    } else {
        double* coef = malloc(sizeof(double) * (dim + 1));
        if (coef && ml_ols(x, y, n, dim, coef))
            out = ml_make_predictor("LinearRegression", coef, dim);
        free(coef);
    }
    free(x); free(y);
    return out;
}

/* LinearModelFit[data] returns the fitted model too, but Wolfram's version is a
 * FittedModel whose properties are regression diagnostics. Only the shared subset is
 * offered: the same coefficients Predict finds, reached through a name a reader of
 * regression code will look for. Its own properties (RSquared, standard errors,
 * ANOVA) are a separate piece of work and are NOT silently approximated. */
/* ---- LinearModelFit[data, funs, vars]: a fit in the span of basis functions ----
 *
 * The model is y = c0 + sum_k c_k f_k(vars). The constant basis function is always present (the
 * reference's IncludeConstantBasis default), so the fit is exactly ml_ols applied to the VALUES of
 * the non-constant basis functions at the data points: those values are the features, and the
 * centring there handles an offset abscissa here as well as in LinearModelFit[data]. A basis
 * element free of every variable is the constant itself and is not a second column (it would be
 * collinear with the intercept).
 *
 * Result: FittedModel["Linear", {bestFit, params, basis, vars}, nvars, nbasis], with
 * bestFit = params . basis. model[v1, ...] substitutes into bestFit; "BestFit",
 * "BestFitParameters", "BasisFunctions" and "Function" are properties. Regression diagnostics
 * are NOT offered (see the spec page). */

static bool lmf_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/* Does `e` contain any of the variable symbols (interned names compare by pointer)? */
static bool lmf_mentions(const Expr* e, Expr* const* vars, size_t nv) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        for (size_t j = 0; j < nv; j++)
            if (e->data.symbol.name == vars[j]->data.symbol.name) return true;
        return false;
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (lmf_mentions(e->data.function.head, vars, nv)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (lmf_mentions(e->data.function.args[i], vars, nv)) return true;
    return false;
}

/* N[{f1, ...} /. {x1 -> p1, ...}] read into out[0 .. nf-1]. */
static bool lmf_basis_values(Expr* const* fs, size_t nf, Expr* const* vars, size_t nv,
                             const double* pt, double* out) {
    Expr** fl = malloc(sizeof(Expr*) * nf);
    Expr** rl = malloc(sizeof(Expr*) * nv);
    if (!fl || !rl) { free(fl); free(rl); return false; }
    for (size_t k = 0; k < nf; k++) fl[k] = expr_copy(fs[k]);
    for (size_t j = 0; j < nv; j++) {
        Expr* ra[2];
        ra[0] = expr_copy(vars[j]);
        ra[1] = expr_new_real(pt[j]);
        rl[j] = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
    }
    Expr* ra2[2];
    ra2[0] = expr_new_function(expr_new_symbol(SYM_List), fl, nf);
    ra2[1] = expr_new_function(expr_new_symbol(SYM_List), rl, nv);
    free(fl); free(rl);
    Expr* na[1];
    na[0] = expr_new_function(expr_new_symbol(SYM_ReplaceAll), ra2, 2);
    Expr* call = expr_new_function(expr_new_symbol(SYM_N), na, 1);
    Expr* v = evaluate(call);
    expr_free(call);
    bool ok = v && lmf_is_list(v) && v->data.function.arg_count == nf;
    double im = 0.0;
    for (size_t k = 0; ok && k < nf; k++)
        ok = na_read_scalar(v->data.function.args[k], &out[k], &im) && im == 0.0
          && isfinite(out[k]);
    expr_free(v);
    return ok;
}

static Expr* lmf_basis_fit(Expr* res) {
    Expr* data = res->data.function.args[0];
    Expr* fe   = res->data.function.args[1];
    Expr* ve   = res->data.function.args[2];
    if (!lmf_is_list(data) || data->data.function.arg_count == 0) return NULL;

    /* Variables: a symbol or a list of symbols. */
    size_t nv = lmf_is_list(ve) ? ve->data.function.arg_count : 1;
    Expr* const* vars = lmf_is_list(ve) ? ve->data.function.args : &res->data.function.args[2];
    if (nv == 0) return NULL;
    for (size_t j = 0; j < nv; j++)
        if (!vars[j] || vars[j]->type != EXPR_SYMBOL) return NULL;

    /* Non-constant basis functions, in order. */
    size_t nf_all = lmf_is_list(fe) ? fe->data.function.arg_count : 1;
    Expr* const* fall = lmf_is_list(fe) ? fe->data.function.args : &res->data.function.args[1];
    Expr** fs = malloc(sizeof(Expr*) * (nf_all ? nf_all : 1));
    if (!fs) return NULL;
    size_t nf = 0;
    for (size_t k = 0; k < nf_all; k++)
        if (lmf_mentions(fall[k], vars, nv)) fs[nf++] = fall[k];
    if (nf == 0) { free(fs); return NULL; }

    /* Data: {y1, y2, ...} with one variable taking 1, 2, ...; or rows {x1, ..., xnv, y}. */
    size_t n = data->data.function.arg_count;
    double* pts  = malloc(sizeof(double) * n * nv);
    double* y    = malloc(sizeof(double) * n);
    double* fx   = malloc(sizeof(double) * n * nf);
    double* coef = malloc(sizeof(double) * (nf + 1));
    Expr* out = NULL;
    bool ok = pts && y && fx && coef;
    bool rows = ok && lmf_is_list(data->data.function.args[0]);
    if (ok && !rows && nv != 1) ok = false;
    double im = 0.0;
    for (size_t i = 0; ok && i < n; i++) {
        Expr* r = data->data.function.args[i];
        if (rows) {
            if (!lmf_is_list(r) || r->data.function.arg_count != nv + 1) { ok = false; break; }
            for (size_t j = 0; ok && j <= nv; j++) {
                double v = 0.0;
                ok = na_read_scalar(r->data.function.args[j], &v, &im) && im == 0.0;
                if (j < nv) pts[i * nv + j] = v; else y[i] = v;
            }
        } else {
            pts[i] = (double)(i + 1);
            ok = na_read_scalar(r, &y[i], &im) && im == 0.0;
        }
    }
    /* The basis functions are arbitrary user expressions, so their values come from the
     * evaluator, one point at a time. */
    for (size_t i = 0; ok && i < n; i++)
        ok = lmf_basis_values(fs, nf, vars, nv, pts + i * nv, fx + i * nf);

    if (ok && ml_ols(fx, y, n, nf, coef)) {
        /* basis = {1, f1, ...}; bestFit = c0 + c1 f1 + ... */
        Expr** bl = malloc(sizeof(Expr*) * (nf + 1));
        Expr** tl = malloc(sizeof(Expr*) * (nf + 1));
        Expr** vl = malloc(sizeof(Expr*) * nv);
        if (bl && tl && vl) {
            bl[0] = expr_new_integer(1);
            tl[0] = expr_new_real(coef[0]);
            for (size_t k = 0; k < nf; k++) {
                Expr* ta[2];
                bl[k + 1] = expr_copy(fs[k]);
                ta[0] = expr_new_real(coef[k + 1]);
                ta[1] = expr_copy(fs[k]);
                tl[k + 1] = expr_new_function(expr_new_symbol(SYM_Times), ta, 2);
            }
            for (size_t j = 0; j < nv; j++) vl[j] = expr_copy(vars[j]);
            Expr* pay[4];
            pay[0] = expr_new_function(expr_new_symbol(SYM_Plus), tl, nf + 1);
            pay[1] = ml_list_of_reals(coef, nf + 1);
            pay[2] = expr_new_function(expr_new_symbol(SYM_List), bl, nf + 1);
            pay[3] = expr_new_function(expr_new_symbol(SYM_List), vl, nv);
            Expr* a4[4];
            a4[0] = expr_new_string("Linear");
            a4[1] = expr_new_function(expr_new_symbol(SYM_List), pay, 4);
            a4[2] = expr_new_integer((int64_t)nv);
            a4[3] = expr_new_integer((int64_t)(nf + 1));
            out = expr_new_function(expr_new_symbol(SYM_FittedModel), a4, 4);
        }
        free(bl); free(tl); free(vl);
    }
    free(fs); free(pts); free(y); free(fx); free(coef);
    return out;
}

/* FittedModel[...][v1, ...] substitutes into the best fit; FittedModel[...]["prop"] reads a
 * property. Anything else stays unevaluated. */
static Expr* lmf_apply(Expr* head, Expr** args, size_t argc) {
    if (head->data.function.arg_count != 4) return NULL;
    Expr* pay = head->data.function.args[1];
    if (!lmf_is_list(pay) || pay->data.function.arg_count != 4) return NULL;
    Expr* best = pay->data.function.args[0];
    Expr* vl   = pay->data.function.args[3];
    if (!lmf_is_list(vl)) return NULL;
    if (argc == 1 && args[0] && args[0]->type == EXPR_STRING) {
        const char* q = args[0]->data.string;
        if (strcmp(q, "BestFit") == 0) return expr_copy(best);
        if (strcmp(q, "BestFitParameters") == 0) return expr_copy(pay->data.function.args[1]);
        if (strcmp(q, "BasisFunctions") == 0) return expr_copy(pay->data.function.args[2]);
        if (strcmp(q, "Function") == 0) {
            Expr* fa[2];
            fa[0] = expr_copy(vl);
            fa[1] = expr_copy(best);
            return expr_new_function(expr_new_symbol(SYM_Function), fa, 2);
        }
        return NULL;
    }
    if (argc == 0 || argc != vl->data.function.arg_count) return NULL;
    Expr** rl = malloc(sizeof(Expr*) * argc);
    if (!rl) return NULL;
    for (size_t j = 0; j < argc; j++) {
        Expr* ra[2];
        ra[0] = expr_copy(vl->data.function.args[j]);
        ra[1] = expr_copy(args[j]);
        rl[j] = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
    }
    Expr* ra2[2];
    ra2[0] = expr_copy(best);
    ra2[1] = expr_new_function(expr_new_symbol(SYM_List), rl, argc);
    free(rl);
    return expr_new_function(expr_new_symbol(SYM_ReplaceAll), ra2, 2);
}

static Expr* builtin_linear_model_fit(Expr* res) {
    if (res->data.function.arg_count == 3) return lmf_basis_fit(res);
    if (res->data.function.arg_count != 1) return NULL;
    size_t n, dim; double* x = NULL; double* y = NULL;
    if (!ml_read_training(res->data.function.args[0], &n, &dim, &x, &y)) return NULL;
    double* coef = malloc(sizeof(double) * (dim + 1));
    if (!coef) { free(x); free(y); return NULL; }
    Expr* out = NULL;
    if (ml_ols(x, y, n, dim, coef)) out = ml_make_predictor("LinearRegression", coef, dim);
    free(coef); free(x); free(y);
    return out;
}

/* ------------------------------------------------------------------------- */

static Expr* ml_reducer_apply(Expr* head, Expr** args, size_t argc);

/* ------------------------------------------------------------------------- */
/* DimensionReducerFunction -- the second model kind on the shared representation */
/* ------------------------------------------------------------------------- */

/* A reducer is the same positional, method-tagged object a predictor is, and that is
 * the payoff for having designed the representation once: this needed no new node type,
 * no new evaluation concept, and no change to eval.c beyond the probe above recognising
 * one more head.
 *
 * Payload: a List whose FIRST row is the training column means and whose remaining
 * `target` rows are the loadings, one component per row. Those two things together are
 * exactly what it takes to project a point that was not in the training set -- which is
 * the whole difference between a reducer and the reduced data.
 */
static Expr* ml_make_reducer(const char* method, const double* mean,
                             const double* evec, size_t dim, size_t target) {
    Expr** rows = malloc(sizeof(Expr*) * (target + 1));
    if (!rows) return NULL;
    rows[0] = ml_list_of_reals(mean, dim);
    for (size_t t = 0; t < target; t++)
        rows[t + 1] = ml_list_of_reals(evec + t * dim, dim);
    Expr* payload = expr_new_function(expr_new_symbol(SYM_List), rows, target + 1);
    free(rows);
    if (!payload) return NULL;
    Expr* a[4];
    a[0] = expr_new_string(method);
    a[1] = payload;
    a[2] = expr_new_integer((int64_t)dim);
    a[3] = expr_new_integer((int64_t)target);
    if (!a[0] || !a[1] || !a[2] || !a[3]) {
        expr_free(a[0]); expr_free(a[1]); expr_free(a[2]); expr_free(a[3]); return NULL;
    }
    return expr_new_function(expr_new_symbol("DimensionReducerFunction"), a, 4);
}

/* DimensionReduction[data, k]. Wolfram's split is kept: DimensionReduce returns the
 * REDUCED DATA, DimensionReduction returns a REUSABLE REDUCER. */
static Expr* builtin_dimension_reduction(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    Expr* kexpr = res->data.function.args[1];
    if (!kexpr || kexpr->type != EXPR_INTEGER || kexpr->data.integer <= 0) return NULL;
    size_t target = (size_t)kexpr->data.integer;

    size_t n, dim; double* x = NULL; bool vec = false;
    if (!ml_read_data(res->data.function.args[0], &n, &dim, &x, &vec)) return NULL;
    if (vec || target > dim) { free(x); return NULL; }

    double* mean = malloc(sizeof(double) * dim);
    double* evec = malloc(sizeof(double) * dim * dim);
    if (!mean || !evec) { free(mean); free(evec); free(x); return NULL; }
    ml_column_mean(x, n, dim, mean);
    Expr* out = NULL;
    if (ml_pca(x, n, dim, false, NULL, NULL, evec))
        out = ml_make_reducer("PrincipalComponentsAnalysis", mean, evec, dim, target);
    free(mean); free(evec); free(x);
    return out;
}

/* Project one row: subtract the TRAINING means, then contract with each loading.
 *
 * Using the training means rather than the incoming batch's is the entire point of a
 * reusable reducer, and getting it wrong would still look plausible: centring a single
 * new point against itself gives all zeros, and on data that happens to sit near the
 * origin the two agree. The test for this therefore uses training data with a
 * deliberately off-centre mean. */
static void ml_project_row(const double* row, const double* mean, const double* load,
                           size_t dim, size_t target, double* out) {
    for (size_t t = 0; t < target; t++) {
        double s = 0.0;
        for (size_t j = 0; j < dim; j++) s += (row[j] - mean[j]) * load[t * dim + j];
        out[t] = s;
    }
}

static Expr* ml_reducer_apply(Expr* head, Expr** args, size_t argc) {
    if (head->data.function.arg_count != 4 || argc != 1) return NULL;
    Expr* mname = head->data.function.args[0];
    Expr* pay   = head->data.function.args[1];
    Expr* dimx  = head->data.function.args[2];
    Expr* tgtx  = head->data.function.args[3];
    if (!mname || mname->type != EXPR_STRING) return NULL;
    if (!pay || pay->type != EXPR_FUNCTION) return NULL;
    if (!dimx || dimx->type != EXPR_INTEGER || !tgtx || tgtx->type != EXPR_INTEGER)
        return NULL;
    size_t dim = (size_t)dimx->data.integer, target = (size_t)tgtx->data.integer;
    if (pay->data.function.arg_count != target + 1) return NULL;

    if (args[0] && args[0]->type == EXPR_STRING) {
        const char* q = args[0]->data.string;
        if (strcmp(q, "Method") == 0) return expr_copy(mname);
        if (strcmp(q, "FeatureCount") == 0) return expr_copy(dimx);
        if (strcmp(q, "ReducedDimension") == 0) return expr_copy(tgtx);
        return NULL;
    }

    double* mean = malloc(sizeof(double) * dim);
    double* load = malloc(sizeof(double) * target * dim);
    if (!mean || !load) { free(mean); free(load); return NULL; }
    double im = 0.0; bool ok = true;
    Expr* mrow = pay->data.function.args[0];
    if (!mrow || mrow->type != EXPR_FUNCTION || mrow->data.function.arg_count != dim)
        ok = false;
    for (size_t j = 0; ok && j < dim; j++)
        ok = na_read_scalar(mrow->data.function.args[j], &mean[j], &im) && im == 0.0;
    for (size_t t = 0; ok && t < target; t++) {
        Expr* lr = pay->data.function.args[t + 1];
        if (!lr || lr->type != EXPR_FUNCTION || lr->data.function.arg_count != dim) {
            ok = false; break;
        }
        for (size_t j = 0; ok && j < dim; j++)
            ok = na_read_scalar(lr->data.function.args[j], &load[t * dim + j], &im)
              && im == 0.0;
    }
    if (!ok) { free(mean); free(load); return NULL; }

    /* One point or a whole matrix of them. A reducer applied to a batch is the common
     * case, so threading it here is better than making every caller Map. */
    size_t rn, rdim; double* rx = NULL; bool rvec = false;
    Expr* out = NULL;
    if (ml_read_data(args[0], &rn, &rdim, &rx, &rvec)) {
        if (rvec && rn == dim) {
            double* o = malloc(sizeof(double) * target);
            if (o) { ml_project_row(rx, mean, load, dim, target, o);
                     out = ml_list_of_reals(o, target); }
            free(o);
        } else if (!rvec && rdim == dim) {
            double* o = malloc(sizeof(double) * rn * target);
            if (o) {
                for (size_t i = 0; i < rn; i++)
                    ml_project_row(rx + i * dim, mean, load, dim, target, o + i * target);
                out = ml_list_matrix(o, rn, target);
            }
            free(o);
        }
        free(rx);
    }
    free(mean); free(load);
    return out;
}

bool ml_model_apply_probe(Expr* head) {
    if (!head || head->type != EXPR_FUNCTION) return false;
    Expr* hh = head->data.function.head;
    if (!hh || hh->type != EXPR_SYMBOL) return false;
    const char* n = hh->data.symbol.name;
    return n == SYM_PredictorFunction
        || n == SYM_FittedModel
        || n == SYM_DimensionReducerFunction
        || n == SYM_ClassifierFunction;
}

Expr* ml_model_apply(Expr* head, Expr** args, size_t argc) {
    if (!head || head->type != EXPR_FUNCTION) return NULL;
    Expr* hh = head->data.function.head;
    if (!hh || hh->type != EXPR_SYMBOL) return NULL;
    if (hh->data.symbol.name == SYM_DimensionReducerFunction)
        return ml_reducer_apply(head, args, argc);
    if (hh->data.symbol.name == SYM_ClassifierFunction)
        return ml_classifier_apply(head, args, argc);
    if (hh->data.symbol.name == SYM_FittedModel)
        return lmf_apply(head, args, argc);
    if (hh->data.symbol.name != SYM_PredictorFunction) return NULL;
    if (head->data.function.arg_count != 4 || argc != 1) return NULL;

    Expr* mname = head->data.function.args[0];
    Expr* cl    = head->data.function.args[1];
    Expr* dimx  = head->data.function.args[2];
    Expr* extrax = head->data.function.args[3];
    if (!mname || mname->type != EXPR_STRING) return NULL;
    if (!cl || cl->type != EXPR_FUNCTION) return NULL;
    if (!dimx || dimx->type != EXPR_INTEGER) return NULL;
    if (!extrax || extrax->type != EXPR_INTEGER) return NULL;
    size_t dim = (size_t)dimx->data.integer;
    bool knn = strcmp(mname->data.string, "NearestNeighbors") == 0;
    if (!knn && cl->data.function.arg_count != dim + 1) return NULL;

    /* Named property access. This is why an Association payload was unnecessary: the
     * place a user asks for a property is the application, and that is here. */
    if (args[0] && args[0]->type == EXPR_STRING) {
        const char* q = args[0]->data.string;
        if (strcmp(q, "Method") == 0) return expr_copy(mname);
        if (strcmp(q, "Coefficients") == 0) return expr_copy(cl);
        if (strcmp(q, "FeatureCount") == 0) return expr_copy(dimx);
        if (knn && strcmp(q, "NeighborCount") == 0) return expr_copy(extrax);
        if (knn && strcmp(q, "TrainingData") == 0)  return expr_copy(cl);
        return NULL;                       /* unknown property: leave unevaluated */
    }

    /* Numeric application. A single-feature model accepts a bare scalar as well as a
     * one-element list, since writing p[3.] for a one-variable regression is the
     * natural thing and rejecting it would be pedantry. */
    double* xin = malloc(sizeof(double) * (dim ? dim : 1));
    if (!xin) return NULL;
    double im = 0.0;
    Expr* a0 = args[0];
    bool ok;
    if (a0 && a0->type == EXPR_FUNCTION && a0->data.function.head->type == EXPR_SYMBOL
        && a0->data.function.head->data.symbol.name == SYM_List) {
        ok = a0->data.function.arg_count == dim;
        for (size_t j = 0; ok && j < dim; j++)
            ok = na_read_scalar(a0->data.function.args[j], &xin[j], &im) && im == 0.0;
    } else {
        ok = (dim == 1) && na_read_scalar(a0, &xin[0], &im) && im == 0.0;
    }
    if (!ok) { free(xin); return NULL; }    /* wrong shape: unevaluated, not guessed */

    if (knn) {
        /* Mean response of the k nearest training rows.
         *
         * A partial selection rather than a sort: only the k smallest are needed, and k
         * is small, so an insertion into a k-sized list beats ordering all n. Ties keep
         * the incumbent, which is the earlier training row, so the answer does not
         * depend on scan order. */
        size_t n = cl->data.function.arg_count;
        size_t k = (size_t)extrax->data.integer;
        if (k == 0 || k > n) k = n;
        double* bd = malloc(sizeof(double) * k);
        double* bv = malloc(sizeof(double) * k);
        if (!bd || !bv) { free(bd); free(bv); free(xin); return NULL; }
        size_t cnt = 0;
        double* row = malloc(sizeof(double) * (dim + 1));
        if (!row) { free(bd); free(bv); free(xin); return NULL; }
        for (size_t i = 0; i < n; i++) {
            Expr* r = cl->data.function.args[i];
            if (!r || r->type != EXPR_FUNCTION
                || r->data.function.arg_count != dim + 1) {
                free(row); free(bd); free(bv); free(xin); return NULL;
            }
            bool rok = true;
            for (size_t j = 0; j <= dim && rok; j++)
                rok = na_read_scalar(r->data.function.args[j], &row[j], &im) && im == 0.0;
            if (!rok) { free(row); free(bd); free(bv); free(xin); return NULL; }
            double d2 = ml_sqdist(xin, row, dim);
            if (cnt < k) {
                size_t q2 = cnt++;
                while (q2 > 0 && bd[q2 - 1] > d2) { bd[q2] = bd[q2 - 1]; bv[q2] = bv[q2 - 1]; q2--; }
                bd[q2] = d2; bv[q2] = row[dim];
            } else if (d2 < bd[k - 1]) {
                size_t q2 = k - 1;
                while (q2 > 0 && bd[q2 - 1] > d2) { bd[q2] = bd[q2 - 1]; bv[q2] = bv[q2 - 1]; q2--; }
                bd[q2] = d2; bv[q2] = row[dim];
            }
        }
        free(row);
        double s2 = 0.0;
        for (size_t i = 0; i < cnt; i++) s2 += bv[i];
        double mean = (cnt > 0) ? s2 / (double)cnt : 0.0;
        free(bd); free(bv); free(xin);
        return (cnt > 0) ? expr_new_real(mean) : NULL;
    }

    double acc = 0.0, c = 0.0;
    if (!na_read_scalar(cl->data.function.args[0], &c, &im)) { free(xin); return NULL; }
    acc = c;
    for (size_t j = 0; j < dim; j++) {
        if (!na_read_scalar(cl->data.function.args[j + 1], &c, &im)) {
            free(xin); return NULL;
        }
        acc += c * xin[j];
    }
    free(xin);
    return expr_new_real(acc);
}

void ml_predict_init(void) {
    symtab_add_builtin("DimensionReduction", builtin_dimension_reduction);
    symtab_get_def("DimensionReduction")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DimensionReduction",
        "DimensionReduction[data, k] returns a DimensionReducerFunction projecting into "
        "k dimensions, applicable to data it was NOT trained on -- the difference from "
        "DimensionReduce[data, k], which returns the reduced training data. New rows are "
        "centred on the TRAINING column means, which is what makes projections "
        "comparable across batches. Accepts one point or a matrix of points, and answers "
        "\"Method\", \"FeatureCount\" and \"ReducedDimension\".");

    symtab_get_def("DimensionReducerFunction")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DimensionReducerFunction",
        "DimensionReducerFunction[method, {means, loadings...}, featureCount, "
        "reducedDimension] is the reusable reducer DimensionReduction returns. Apply it "
        "to a feature vector, or to a matrix of them, to project into the reduced "
        "space.");

    symtab_add_builtin("Predict", builtin_predict);
    symtab_get_def("Predict")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Predict",
        "Predict[data] fits a predictor to data and returns a PredictorFunction, "
        "which can be stored and applied to new inputs. Data is either a list of "
        "rules {features -> value, ...} or a matrix whose last column is the "
        "response. Method -> \"LinearRegression\" is the only method implemented and "
        "is the default; any other declines rather than silently linear-regressing. "
        "The returned object also answers \"Method\", \"Coefficients\" and "
        "\"FeatureCount\". A collinear feature set has no unique fit and returns "
        "unevaluated.");

    symtab_add_builtin("LinearModelFit", builtin_linear_model_fit);
    symtab_get_def("LinearModelFit")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LinearModelFit",
        "LinearModelFit[data] fits a linear model with an intercept and returns a "
        "PredictorFunction carrying its coefficients. LinearModelFit[data, {f1, ...}, vars] "
        "fits c0 + c1 f1 + ... by least squares (data {y1, ...} at x = 1, 2, ..., or rows "
        "{x1, ..., y}) and returns a FittedModel: model[x] evaluates the fit, and "
        "\"BestFit\", \"BestFitParameters\", \"BasisFunctions\", \"Function\" are its "
        "properties. Wolfram's version returns a "
        "FittedModel with regression diagnostics (RSquared, standard errors, ANOVA); "
        "those are not implemented and are not approximated. The coefficients are the "
        "same ones Predict finds.");

    symtab_get_def("FittedModel")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FittedModel",
        "FittedModel[\"Linear\", {bestFit, params, basis, vars}, nvars, nbasis] is the model "
        "LinearModelFit[data, funs, vars] returns. model[x, ...] evaluates the fit; "
        "\"BestFit\", \"BestFitParameters\", \"BasisFunctions\" and \"Function\" are its "
        "properties.");

    symtab_get_def("PredictorFunction")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PredictorFunction",
        "PredictorFunction[method, coefficients, featureCount] is the fitted object "
        "Predict returns. Apply it to a feature vector to get a prediction, or to "
        "\"Method\", \"Coefficients\" or \"FeatureCount\" to read it. A one-feature "
        "model also accepts a bare scalar.");
}
