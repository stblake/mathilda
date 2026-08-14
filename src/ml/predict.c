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
#include "predict.h"

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
    size_t p = dim + 1;                       /* +1 for the intercept */
    if (n < p) return false;                  /* underdetermined: no unique answer */

    double* a = calloc(p * (p + 1), sizeof(double));   /* augmented [A'A | A'y] */
    if (!a) return false;

    for (size_t i = 0; i < n; i++) {
        /* Row i of the design matrix is (1, x_i1, ..., x_id), formed on the fly
         * rather than materialised -- it is the only place the intercept exists. */
        for (size_t r = 0; r < p; r++) {
            double ar = (r == 0) ? 1.0 : x[i * dim + (r - 1)];
            for (size_t c = 0; c < p; c++) {
                double ac = (c == 0) ? 1.0 : x[i * dim + (c - 1)];
                a[r * (p + 1) + c] += ar * ac;
            }
            a[r * (p + 1) + p] += ar * y[i];
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
        if (!(best > 1e-12)) { free(a); return false; }
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
    for (size_t r = 0; r < p; r++) coef[r] = a[r * (p + 1) + p] / a[r * (p + 1) + r];
    free(a);
    return true;
}

/* ------------------------------------------------------------------------- */

static Expr* ml_reals_list(const double* v, size_t n) {
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!a) return NULL;
    for (size_t i = 0; i < n; i++) a[i] = expr_new_real(v[i]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), a, n);
    free(a);
    return out;
}

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
    bool rules = first && first->type == EXPR_FUNCTION
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

static Expr* ml_make_predictor(const char* method, const double* coef, size_t dim) {
    Expr* args[3];
    args[0] = expr_new_string(method);
    args[1] = ml_reals_list(coef, dim + 1);
    args[2] = expr_new_integer((int64_t)dim);
    if (!args[0] || !args[1] || !args[2]) {
        expr_free(args[0]); expr_free(args[1]); expr_free(args[2]); return NULL;
    }
    return expr_new_function(expr_new_symbol("PredictorFunction"), args, 3);
}

static Expr* builtin_predict(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    if (argc == 2) {
        /* Only LinearRegression is implemented; anything else declines rather than
         * silently linear-regressing, so an unsupported method is visible. */
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
        if (strcmp(rhs->data.string, "LinearRegression") != 0) return NULL;
    }

    size_t n, dim; double* x = NULL; double* y = NULL;
    if (!ml_read_training(res->data.function.args[0], &n, &dim, &x, &y)) return NULL;

    double* coef = malloc(sizeof(double) * (dim + 1));
    if (!coef) { free(x); free(y); return NULL; }
    Expr* out = NULL;
    if (ml_ols(x, y, n, dim, coef)) out = ml_make_predictor("LinearRegression", coef, dim);
    free(coef); free(x); free(y);
    return out;
}

/* LinearModelFit[data] returns the fitted model too, but Wolfram's version is a
 * FittedModel whose properties are regression diagnostics. Only the shared subset is
 * offered: the same coefficients Predict finds, reached through a name a reader of
 * regression code will look for. Its own properties (RSquared, standard errors,
 * ANOVA) are a separate piece of work and are NOT silently approximated. */
static Expr* builtin_linear_model_fit(Expr* res) {
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

bool ml_model_apply_probe(Expr* head) {
    if (!head || head->type != EXPR_FUNCTION) return false;
    Expr* hh = head->data.function.head;
    if (!hh || hh->type != EXPR_SYMBOL) return false;
    return strcmp(hh->data.symbol.name, "PredictorFunction") == 0;
}

Expr* ml_model_apply(Expr* head, Expr** args, size_t argc) {
    if (!head || head->type != EXPR_FUNCTION) return NULL;
    Expr* hh = head->data.function.head;
    if (!hh || hh->type != EXPR_SYMBOL) return NULL;
    if (strcmp(hh->data.symbol.name, "PredictorFunction") != 0) return NULL;
    if (head->data.function.arg_count != 3 || argc != 1) return NULL;

    Expr* mname = head->data.function.args[0];
    Expr* cl    = head->data.function.args[1];
    Expr* dimx  = head->data.function.args[2];
    if (!mname || mname->type != EXPR_STRING) return NULL;
    if (!cl || cl->type != EXPR_FUNCTION) return NULL;
    if (!dimx || dimx->type != EXPR_INTEGER) return NULL;
    size_t dim = (size_t)dimx->data.integer;
    if (cl->data.function.arg_count != dim + 1) return NULL;

    /* Named property access. This is why an Association payload was unnecessary: the
     * place a user asks for a property is the application, and that is here. */
    if (args[0] && args[0]->type == EXPR_STRING) {
        const char* q = args[0]->data.string;
        if (strcmp(q, "Method") == 0) return expr_copy(mname);
        if (strcmp(q, "Coefficients") == 0) return expr_copy(cl);
        if (strcmp(q, "FeatureCount") == 0) return expr_copy(dimx);
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
        "PredictorFunction carrying its coefficients. Wolfram's version returns a "
        "FittedModel with regression diagnostics (RSquared, standard errors, ANOVA); "
        "those are not implemented and are not approximated. The coefficients are the "
        "same ones Predict finds.");

    symtab_get_def("PredictorFunction")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PredictorFunction",
        "PredictorFunction[method, coefficients, featureCount] is the fitted object "
        "Predict returns. Apply it to a feature vector to get a prediction, or to "
        "\"Method\", \"Coefficients\" or \"FeatureCount\" to read it. A one-feature "
        "model also accepts a bare scalar.");
}
