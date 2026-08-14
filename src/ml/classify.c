/* classify.c -- Classify and ClassifierFunction.
 *
 * A ClassifierFunction is the FOURTH head on the model representation designed in
 * src/ml/predict.h, after PredictorFunction, DimensionReducerFunction and
 * LearnedDistribution. It needed no change to that design: the payload shape varies by
 * method, which is exactly what a positional method-tagged representation is for.
 *
 * What IS new is the label vocabulary (src/ml/encode.h). Every earlier family took numeric
 * responses; a class is an arbitrary expression, so the vocabulary is the bridge between
 * "the user's classes" and "indices an algorithm can count with". It lives in its own
 * module because a ContingencyTable and a categorical FEATURE encoder will both need it.
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
#include "encode.h"
#include "classify.h"

/* Read {features -> class, ...} into a feature matrix and a label vocabulary.
 *
 * Only the rule form is accepted, and that is not laziness: a matrix with the class in its
 * last column cannot work here, because a class is not necessarily a number and
 * na_load_matrix would refuse the whole thing. Predict can take a matrix precisely because
 * its response IS numeric. */
static bool ml_read_labelled(Expr* data, size_t* n_out, size_t* dim_out,
                             double** x_out, size_t** y_out, MlLabels* labels) {
    if (!data || data->type != EXPR_FUNCTION) return false;
    Expr* h = data->data.function.head;
    if (!h || h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_List) return false;
    size_t n = data->data.function.arg_count;
    if (n == 0) return false;

    Expr* first = data->data.function.args[0];
    if (!first || first->type != EXPR_FUNCTION
        || first->data.function.head->type != EXPR_SYMBOL
        || (first->data.function.head->data.symbol.name != SYM_Rule
            && first->data.function.head->data.symbol.name != SYM_RuleDelayed))
        return false;

    Expr* lhs0 = first->data.function.args[0];
    size_t dim = 1;
    if (lhs0 && lhs0->type == EXPR_FUNCTION
        && lhs0->data.function.head->type == EXPR_SYMBOL
        && lhs0->data.function.head->data.symbol.name == SYM_List)
        dim = lhs0->data.function.arg_count;
    if (dim == 0) return false;

    double* x = malloc(sizeof(double) * n * dim);
    Expr** cls = malloc(sizeof(Expr*) * n);
    size_t* y = malloc(sizeof(size_t) * n);
    if (!x || !cls || !y) { free(x); free(cls); free(y); return false; }

    double im = 0.0; bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        Expr* r = data->data.function.args[i];
        if (!r || r->type != EXPR_FUNCTION || r->data.function.arg_count != 2
            || r->data.function.head->type != EXPR_SYMBOL
            || (r->data.function.head->data.symbol.name != SYM_Rule
                && r->data.function.head->data.symbol.name != SYM_RuleDelayed)) {
            ok = false; break;
        }
        Expr* lhs = r->data.function.args[0];
        if (dim == 1 && !(lhs && lhs->type == EXPR_FUNCTION)) {
            ok = na_read_scalar(lhs, &x[i * dim], &im) && im == 0.0;
        } else {
            if (!lhs || lhs->type != EXPR_FUNCTION
                || lhs->data.function.arg_count != dim) { ok = false; break; }
            for (size_t j = 0; j < dim && ok; j++)
                ok = na_read_scalar(lhs->data.function.args[j], &x[i * dim + j], &im)
                  && im == 0.0;
        }
        cls[i] = r->data.function.args[1];      /* borrowed; copied by the vocabulary */
    }
    if (ok) ok = ml_labels_build(cls, n, labels, y);
    free(cls);
    if (!ok) { free(x); free(y); return false; }
    *n_out = n; *dim_out = dim; *x_out = x; *y_out = y;
    return true;
}

static Expr* builtin_classify(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    size_t kopt = 0;
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
        Expr* mstr = rhs;
        if (rhs && rhs->type == EXPR_FUNCTION
            && rhs->data.function.head->type == EXPR_SYMBOL
            && rhs->data.function.head->data.symbol.name == SYM_List
            && rhs->data.function.arg_count >= 1) {
            mstr = rhs->data.function.args[0];
            for (size_t i = 1; i < rhs->data.function.arg_count; i++) {
                Expr* so = rhs->data.function.args[i];
                if (!so || so->type != EXPR_FUNCTION || so->data.function.arg_count != 2)
                    return NULL;
                Expr* sk = so->data.function.args[0];
                Expr* sv = so->data.function.args[1];
                if (!sk || sk->type != EXPR_STRING
                    || strcmp(sk->data.string, "NeighborsNumber") != 0) return NULL;
                if (!sv || sv->type != EXPR_INTEGER || sv->data.integer <= 0) return NULL;
                kopt = (size_t)sv->data.integer;
            }
        }
        if (!mstr || mstr->type != EXPR_STRING) return NULL;
        if (strcmp(mstr->data.string, "NearestNeighbors") != 0) return NULL;
    }

    size_t n, dim; double* x = NULL; size_t* y = NULL;
    MlLabels labels = {0, NULL};
    if (!ml_read_labelled(res->data.function.args[0], &n, &dim, &x, &y, &labels))
        return NULL;

    /* k defaults to 1 rather than 3, unlike the k-NN PREDICTOR. A regression averages, so
     * a little smoothing helps; a classifier votes, and at k = 1 it reproduces the
     * training labels exactly -- which is both the standard default for a nearest-neighbour
     * classifier and an exactly checkable property. */
    size_t k = kopt ? kopt : 1;
    if (k > n) k = n;

    /* Payload: the vocabulary, then one row per example holding its features followed by
     * its CLASS INDEX. Indices rather than the labels themselves, so the numeric part of
     * the model stays numeric and the vocabulary is the single place a class is named. */
    double* joint = malloc(sizeof(double) * n * (dim + 1));
    Expr* out = NULL;
    if (joint) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < dim; j++) joint[i * (dim + 1) + j] = x[i * dim + j];
            joint[i * (dim + 1) + dim] = (double)y[i];
        }
        Expr** rows = malloc(sizeof(Expr*) * (n + 1));
        if (rows) {
            rows[0] = ml_labels_to_list(&labels);
            for (size_t i = 0; i < n; i++)
                rows[i + 1] = ml_list_of_reals(joint + i * (dim + 1), dim + 1);
            Expr* pay = expr_new_function(expr_new_symbol(SYM_List), rows, n + 1);
            free(rows);
            if (pay) {
                Expr* a4[4];
                a4[0] = expr_new_string("NearestNeighbors");
                a4[1] = pay;
                a4[2] = expr_new_integer((int64_t)dim);
                a4[3] = expr_new_integer((int64_t)k);
                if (a4[0] && a4[1] && a4[2] && a4[3])
                    out = expr_new_function(expr_new_symbol("ClassifierFunction"), a4, 4);
                else { expr_free(a4[0]); expr_free(a4[1]); expr_free(a4[2]); expr_free(a4[3]); }
            }
        }
        free(joint);
    }
    ml_labels_free(&labels);
    free(x); free(y);
    return out;
}

/* Apply a fitted classifier. Handles ClassifierFunction[...][x] and
 * ClassifierFunction[...][x, "Probabilities"], plus the property queries. */
Expr* ml_classifier_apply(Expr* head, Expr** args, size_t argc) {
    if (head->data.function.arg_count != 4 || argc < 1 || argc > 2) return NULL;
    Expr* mname = head->data.function.args[0];
    Expr* pay   = head->data.function.args[1];
    Expr* dimx  = head->data.function.args[2];
    Expr* kx    = head->data.function.args[3];
    if (!mname || mname->type != EXPR_STRING) return NULL;
    if (!pay || pay->type != EXPR_FUNCTION || pay->data.function.arg_count < 2) return NULL;
    if (!dimx || dimx->type != EXPR_INTEGER || !kx || kx->type != EXPR_INTEGER) return NULL;
    size_t dim = (size_t)dimx->data.integer, k = (size_t)kx->data.integer;
    size_t n = pay->data.function.arg_count - 1;

    MlLabels labels = {0, NULL};
    if (!ml_labels_from_list(pay->data.function.args[0], &labels)) return NULL;

    if (argc == 1 && args[0] && args[0]->type == EXPR_STRING) {
        const char* q = args[0]->data.string;
        Expr* r = NULL;
        if (strcmp(q, "Method") == 0)            r = expr_copy(mname);
        else if (strcmp(q, "FeatureCount") == 0) r = expr_copy(dimx);
        else if (strcmp(q, "NeighborCount") == 0) r = expr_copy(kx);
        else if (strcmp(q, "Classes") == 0)      r = ml_labels_to_list(&labels);
        ml_labels_free(&labels);
        return r;                                /* NULL for an unknown property */
    }

    bool want_probs = false;
    if (argc == 2) {
        if (!args[1] || args[1]->type != EXPR_STRING
            || strcmp(args[1]->data.string, "Probabilities") != 0) {
            ml_labels_free(&labels); return NULL;
        }
        want_probs = true;
    }

    /* Read the query point. A one-feature classifier accepts a bare scalar. */
    double* xin = malloc(sizeof(double) * (dim ? dim : 1));
    if (!xin) { ml_labels_free(&labels); return NULL; }
    double im = 0.0; bool ok;
    Expr* a0 = args[0];
    if (a0 && a0->type == EXPR_FUNCTION && a0->data.function.head->type == EXPR_SYMBOL
        && a0->data.function.head->data.symbol.name == SYM_List) {
        ok = a0->data.function.arg_count == dim;
        for (size_t j = 0; ok && j < dim; j++)
            ok = na_read_scalar(a0->data.function.args[j], &xin[j], &im) && im == 0.0;
    } else {
        ok = (dim == 1) && na_read_scalar(a0, &xin[0], &im) && im == 0.0;
    }
    if (!ok) { free(xin); ml_labels_free(&labels); return NULL; }

    /* k nearest by squared Euclidean distance, kept as a partial selection. Ties keep the
     * earlier training row, so the vote does not depend on scan order. */
    if (k == 0 || k > n) k = n;
    double* bd = malloc(sizeof(double) * k);
    size_t* bc = malloc(sizeof(size_t) * k);
    double* row = malloc(sizeof(double) * (dim + 1));
    if (!bd || !bc || !row) {
        free(bd); free(bc); free(row); free(xin); ml_labels_free(&labels); return NULL;
    }
    size_t cnt = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* r = pay->data.function.args[i + 1];
        if (!r || r->type != EXPR_FUNCTION || r->data.function.arg_count != dim + 1) {
            ok = false; break;
        }
        for (size_t j = 0; j <= dim && ok; j++)
            ok = na_read_scalar(r->data.function.args[j], &row[j], &im) && im == 0.0;
        if (!ok) break;
        double d2 = ml_sqdist(xin, row, dim);
        size_t ci = (size_t)row[dim];
        if (cnt < k) {
            size_t q = cnt++;
            while (q > 0 && bd[q - 1] > d2) { bd[q] = bd[q - 1]; bc[q] = bc[q - 1]; q--; }
            bd[q] = d2; bc[q] = ci;
        } else if (d2 < bd[k - 1]) {
            size_t q = k - 1;
            while (q > 0 && bd[q - 1] > d2) { bd[q] = bd[q - 1]; bc[q] = bc[q - 1]; q--; }
            bd[q] = d2; bc[q] = ci;
        }
    }
    free(row); free(xin);
    if (!ok || cnt == 0) { free(bd); free(bc); ml_labels_free(&labels); return NULL; }

    size_t* votes = calloc(labels.count, sizeof(size_t));
    if (!votes) { free(bd); free(bc); ml_labels_free(&labels); return NULL; }
    for (size_t i = 0; i < cnt; i++)
        if (bc[i] < labels.count) votes[bc[i]]++;
    free(bd); free(bc);

    Expr* out = NULL;
    if (want_probs) {
        /* One Rule per class, class -> share of the k votes. These sum to 1 by
         * construction, which is asserted as an absolute check: a classifier whose
         * probabilities did not sum to 1 would still look plausible per class. */
        Expr** rules = malloc(sizeof(Expr*) * labels.count);
        if (rules) {
            for (size_t c = 0; c < labels.count; c++) {
                Expr* pr[2];
                pr[0] = expr_copy(labels.label[c]);
                pr[1] = expr_new_real((double)votes[c] / (double)cnt);
                rules[c] = expr_new_function(expr_new_symbol(SYM_Rule), pr, 2);
            }
            out = expr_new_function(expr_new_symbol(SYM_List), rules, labels.count);
            free(rules);
        }
    } else {
        /* Majority vote; the lowest class index wins a tie, which is first appearance in
         * the training data and therefore deterministic. */
        size_t best = 0;
        for (size_t c = 1; c < labels.count; c++) if (votes[c] > votes[best]) best = c;
        out = expr_copy(labels.label[best]);
    }
    free(votes);
    ml_labels_free(&labels);
    return out;
}

void ml_classify_init(void) {
    symtab_add_builtin("Classify", builtin_classify);
    symtab_get_def("Classify")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Classify",
        "Classify[data] trains a classifier and returns a ClassifierFunction. Data is a "
        "list of rules {features -> class, ...}; a class may be any expression -- a "
        "string, a symbol, a number -- and the distinct classes are numbered by first "
        "appearance. Method -> \"NearestNeighbors\" is the only method implemented and is "
        "the default, with NeighborsNumber defaulting to 1: a classifier votes rather "
        "than averages, so at k = 1 it reproduces its training labels exactly. Apply the "
        "result to a feature vector for a class, or with \"Probabilities\" for the vote "
        "shares. It also answers \"Classes\", \"Method\", \"FeatureCount\" and "
        "\"NeighborCount\".");

    symtab_get_def("ClassifierFunction")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ClassifierFunction",
        "ClassifierFunction[method, parameters, featureCount, k] is the fitted classifier "
        "Classify returns. Apply it to a feature vector to get a class, or to a feature "
        "vector and \"Probabilities\" to get one rule per class.");
}
