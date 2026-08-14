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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "numarray.h"
#include "mlutil.h"
#include "encode.h"
#include "tree.h"
#include "pca.h"      /* ml_column_mean, ml_column_sd */
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

/* One binary logistic fit by iteratively reweighted least squares -- Newton's method on the
 * log-likelihood, each step a WEIGHTED least squares with weights p(1-p) and working response
 * eta + (y - p)/(p(1-p)).
 *
 * `pos` is the class index treated as the positive outcome, and it is the whole reason this is
 * a function rather than the loop it used to be: one-vs-rest runs the identical iteration K
 * times, changing only which class counts as 1. Extracted at the second real consumer, not in
 * anticipation of one.
 *
 * A SMALL RIDGE, AND WHY IT IS NOT COSMETIC. On linearly separable data the likelihood is
 * UNBOUNDED: driving the coefficients to infinity sends every fitted probability to 0 or 1 and
 * the likelihood to its supremum, so plain Newton diverges and never converges. Same shape of
 * problem as a mixture's unbounded likelihood, and it gets the same honest treatment rather
 * than a silent hang -- a ridge on the non-intercept coefficients makes the penalised objective
 * strictly concave, so the fit is finite and unique even on separable data, with the iteration
 * count capped as a backstop. The intercept is left unpenalised, which is standard: shrinking it
 * would bias the predicted base rate.
 *
 * One-vs-rest makes separability the COMMON case rather than a corner one. With K well-spaced
 * classes every binary sub-problem -- class k against all the others -- is separable by
 * construction, so the ridge is what each of the K fits relies on to come back finite at all.
 * It carries more weight here than it ever did in the two-class fit.
 *
 * Returns false only when a step's normal-equation matrix is singular. `beta` (dim + 1 entries,
 * intercept first) is filled on success. */
static bool logit_irls(const double* x, const size_t* y, size_t n, size_t dim,
                       size_t pos, double* beta) {
    size_t p1 = dim + 1;
    double* eta = malloc(sizeof(double) * n);
    double* aug = malloc(sizeof(double) * p1 * (p1 + 1));
    if (!eta || !aug) { free(eta); free(aug); return false; }
    for (size_t r = 0; r < p1; r++) beta[r] = 0.0;
    bool singular = false;
    const double ridge = 1e-6;
    bool converged = false;
    for (int it = 0; it < 100 && !converged; it++) {
        for (size_t a = 0; a < p1 * (p1 + 1); a++) aug[a] = 0.0;
        for (size_t i = 0; i < n; i++) {
            double e = beta[0];
            for (size_t j = 0; j < dim; j++) e += beta[j + 1] * x[i * dim + j];
            eta[i] = e;
            double pi = 1.0 / (1.0 + exp(-e));
            /* Clamp the weight away from zero: a saturated point contributes
             * p(1-p) ~ 0, and dividing by it in the working response is exactly
             * where a separable fit blows up. */
            double w = pi * (1.0 - pi);
            if (w < 1e-10) w = 1e-10;
            double z = e + ((double)(y[i] == pos ? 1.0 : 0.0) - pi) / w;
            for (size_t r = 0; r < p1; r++) {
                double ar = (r == 0) ? 1.0 : x[i * dim + (r - 1)];
                for (size_t c = 0; c < p1; c++) {
                    double ac = (c == 0) ? 1.0 : x[i * dim + (c - 1)];
                    aug[r * (p1 + 1) + c] += w * ar * ac;
                }
                aug[r * (p1 + 1) + p1] += w * ar * z;
            }
        }
        for (size_t r = 1; r < p1; r++) aug[r * (p1 + 1) + r] += ridge;

        bool ok2 = true;
        for (size_t c = 0; c < p1 && ok2; c++) {
            size_t piv = c; double best = fabs(aug[c * (p1 + 1) + c]);
            for (size_t r = c + 1; r < p1; r++) {
                double v = fabs(aug[r * (p1 + 1) + c]);
                if (v > best) { best = v; piv = r; }
            }
            if (!(best > 1e-14)) { ok2 = false; break; }
            if (piv != c)
                for (size_t kk = 0; kk <= p1; kk++) {
                    double t = aug[c * (p1 + 1) + kk];
                    aug[c * (p1 + 1) + kk] = aug[piv * (p1 + 1) + kk];
                    aug[piv * (p1 + 1) + kk] = t;
                }
            for (size_t r = 0; r < p1; r++) {
                if (r == c) continue;
                double f = aug[r * (p1 + 1) + c] / aug[c * (p1 + 1) + c];
                if (f == 0.0) continue;
                for (size_t kk = c; kk <= p1; kk++)
                    aug[r * (p1 + 1) + kk] -= f * aug[c * (p1 + 1) + kk];
            }
        }
        if (!ok2) { singular = true; break; }
        double delta = 0.0;
        for (size_t r = 0; r < p1; r++) {
            double nb = aug[r * (p1 + 1) + p1] / aug[r * (p1 + 1) + r];
            double d = fabs(nb - beta[r]);
            if (d > delta) delta = d;
            beta[r] = nb;
        }
        if (delta < 1e-10) converged = true;
    }
    free(eta); free(aug);
    return !singular;
}

static Expr* builtin_classify(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    size_t kopt = 0;
    bool bayes = false, logit = false, tree = false;
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
        const char* mm = mstr->data.string;
        if (strcmp(mm, "NearestNeighbors") == 0)         { bayes = false; logit = false; }
        else if (strcmp(mm, "NaiveBayes") == 0)           { bayes = true;  logit = false; }
        else if (strcmp(mm, "LogisticRegression") == 0)   { bayes = false; logit = true; }
        else if (strcmp(mm, "DecisionTree") == 0)          { bayes = false; logit = false; tree = true; }
        else return NULL;
        /* A neighbour count means nothing to a Bayes classifier, so it is refused rather
         * than accepted and ignored -- silently swallowing it would hide a real mistake. */
        if (kopt && (bayes || logit || tree)) return NULL;
    }

    size_t n, dim; double* x = NULL; size_t* y = NULL;
    MlLabels labels = {0, NULL};
    if (!ml_read_labelled(res->data.function.args[0], &n, &dim, &x, &y, &labels))
        return NULL;

    if (logit) {
        /* Two-class logistic regression by iteratively reweighted least squares, which is
         * Newton's method on the log-likelihood: each step solves a WEIGHTED least squares
         * with weights p(1-p) and working response eta + (y - p)/(p(1-p)).
         *
         * TWO CLASSES ONLY. Multi-class would be one-vs-rest or a softmax, and rather than
         * guess which, this declines -- the caller can see that it declined, where a
         * silently-binarised three-class problem would just be wrong.
         *
         * A SMALL RIDGE, AND WHY IT IS NOT COSMETIC. On linearly separable data the
         * likelihood is UNBOUNDED: pushing the coefficients to infinity drives every fitted
         * probability to 0 or 1 and the likelihood to its supremum, so plain Newton diverges
         * and never converges. This is the same shape of problem as a mixture's unbounded
         * likelihood, and it gets the same honest treatment rather than a silent hang: a
         * ridge on the non-intercept coefficients makes the penalised objective strictly
         * concave, so the fit is finite and unique even when the data is separable, and the
         * iteration count is capped as a backstop. The intercept is left unpenalised, which
         * is standard -- shrinking it would bias the predicted base rate. */
        /* TWO CLASSES stay exactly as they were: one fit, one coefficient vector, the same
         * payload shape and the same numbers. At K = 2 one-vs-rest would fit two models that
         * are reflections of each other, so the single fit is not a shortcut but the right
         * answer -- and its coefficients are pinned by tests.
         *
         * MORE THAN TWO is one-vs-rest: K binary fits, class k against everything else, and
         * the predicted class is the arg-max of the K fitted probabilities. Chosen over a
         * softmax for two reasons. It reuses this exact iteration K times instead of needing
         * a different one. And a softmax's parameters are identified only up to an additive
         * constant per feature, so the stored coefficients would not be unique -- meaning no
         * test could pin them, which is precisely the property that caught real bugs in the
         * other four families. */
        if (labels.count < 2) {
            ml_labels_free(&labels); free(x); free(y); return NULL;
        }
        size_t p1 = dim + 1;
        size_t nfit = (labels.count == 2) ? 1 : labels.count;
        double* beta = calloc(nfit * p1, sizeof(double));
        Expr* outl = NULL;
        if (beta) {
            bool okfit = true;
            for (size_t k = 0; k < nfit && okfit; k++) {
                /* Two classes keep the original sign convention, where the positive outcome
                 * is class index 1; one-vs-rest makes it class k. */
                size_t pos = (nfit == 1) ? 1 : k;
                okfit = logit_irls(x, y, n, dim, pos, beta + k * p1);
            }
            Expr* coef = NULL;
            if (okfit) {
                if (nfit == 1) {
                    coef = ml_list_of_reals(beta, p1);
                } else {
                    Expr** vs = malloc(sizeof(Expr*) * nfit);
                    if (vs) {
                        bool allv = true;
                        for (size_t k = 0; k < nfit; k++) {
                            vs[k] = ml_list_of_reals(beta + k * p1, p1);
                            if (!vs[k]) allv = false;
                        }
                        if (allv) {
                            coef = expr_new_function(expr_new_symbol(SYM_List), vs, nfit);
                        } else {
                            for (size_t k = 0; k < nfit; k++) expr_free(vs[k]);
                        }
                        free(vs);
                    }
                }
            }
            if (coef) {
                Expr** rows = malloc(sizeof(Expr*) * 2);
                if (rows) {
                    rows[0] = ml_labels_to_list(&labels);
                    rows[1] = coef;
                    Expr* pay = expr_new_function(expr_new_symbol(SYM_List), rows, 2);
                    free(rows);
                    if (pay) {
                        Expr* a4[4];
                        a4[0] = expr_new_string("LogisticRegression");
                        a4[1] = pay;
                        a4[2] = expr_new_integer((int64_t)dim);
                        a4[3] = expr_new_integer(0);
                        if (a4[0] && a4[1] && a4[2] && a4[3]) {
                            outl = expr_new_function(expr_new_symbol("ClassifierFunction"), a4, 4);
                        } else {
                            expr_free(a4[0]); expr_free(a4[1]);
                            expr_free(a4[2]); expr_free(a4[3]);
                        }
                    }
                } else {
                    expr_free(coef);
                }
            }
        }
        free(beta);
        ml_labels_free(&labels); free(x); free(y);
        return outl;
    }

    if (tree) {
        /* A CART classification tree. The kernel is src/ml/tree.c; this is the Expr side.
         *
         * DEPTH 32 AND MIN-SPLIT 2 MEAN THE TREE GROWS UNTIL EVERY LEAF IS PURE OR
         * UNSPLITTABLE, which is a deliberate default rather than an oversight. It makes
         * "reproduces every training label" an exact property to assert instead of an accuracy
         * figure to hope for -- the role k = 1 plays for the nearest-neighbour classifier. It
         * also overfits, which is the honest trade and is stated in the docs: pruning needs a
         * validation split or a complexity parameter, and inventing either silently would be
         * worse than growing the tree that was asked for. Depth 32 bounds the work without
         * being reachable on data that leaf purity would not stop first.
         *
         * The payload is the vocabulary plus TWO matrices with one row per node: the split
         * (feature, threshold, left, right) and the class counts. Counts at every node rather
         * than only at leaves is what makes "Probabilities" fall out of the same array as the
         * class, with no separate leaf table to keep in step. */
        MlTree* tt = ml_tree_fit(x, y, n, dim, labels.count, 32, 2);
        Expr* outt = NULL;
        if (tt) {
            Expr** nrows = malloc(sizeof(Expr*) * tt->count);
            Expr** drows = malloc(sizeof(Expr*) * tt->count);
            bool ok2 = nrows && drows;
            if (ok2) {
                for (size_t i = 0; i < tt->count; i++) { nrows[i] = NULL; drows[i] = NULL; }
                for (size_t i = 0; i < tt->count && ok2; i++) {
                    Expr* four[4];
                    four[0] = expr_new_integer(tt->feature[i]);
                    four[1] = expr_new_real(tt->thresh[i]);
                    four[2] = expr_new_integer((int64_t)tt->left[i]);
                    four[3] = expr_new_integer((int64_t)tt->right[i]);
                    if (four[0] && four[1] && four[2] && four[3])
                        nrows[i] = expr_new_function(expr_new_symbol(SYM_List), four, 4);
                    else { expr_free(four[0]); expr_free(four[1]);
                           expr_free(four[2]); expr_free(four[3]); }
                    drows[i] = ml_list_of_reals(tt->dist + i * tt->k, tt->k);
                    if (!nrows[i] || !drows[i]) ok2 = false;
                }
            }
            Expr* nodes = NULL; Expr* dists = NULL;
            if (ok2) {
                nodes = expr_new_function(expr_new_symbol(SYM_List), nrows, tt->count);
                dists = expr_new_function(expr_new_symbol(SYM_List), drows, tt->count);
                if (!nodes || !dists) ok2 = false;
            }
            if (!ok2) {
                for (size_t i = 0; i < tt->count; i++) {
                    if (nrows) expr_free(nrows[i]);
                    if (drows) expr_free(drows[i]);
                }
                expr_free(nodes); expr_free(dists);
                nodes = NULL; dists = NULL;
            }
            free(nrows); free(drows);
            if (nodes && dists) {
                Expr** rows = malloc(sizeof(Expr*) * 3);
                if (rows) {
                    rows[0] = ml_labels_to_list(&labels);
                    rows[1] = nodes;
                    rows[2] = dists;
                    Expr* pay = rows[0]
                        ? expr_new_function(expr_new_symbol(SYM_List), rows, 3) : NULL;
                    if (!pay) { expr_free(rows[0]); expr_free(nodes); expr_free(dists); }
                    free(rows);
                    if (pay) {
                        Expr* a4[4];
                        a4[0] = expr_new_string("DecisionTree");
                        a4[1] = pay;
                        a4[2] = expr_new_integer((int64_t)dim);
                        a4[3] = expr_new_integer((int64_t)tt->count);
                        if (a4[0] && a4[1] && a4[2] && a4[3])
                            outt = expr_new_function(expr_new_symbol("ClassifierFunction"), a4, 4);
                        else { expr_free(a4[0]); expr_free(a4[1]);
                               expr_free(a4[2]); expr_free(a4[3]); }
                    }
                } else { expr_free(nodes); expr_free(dists); }
            }
            ml_tree_free(tt);
        }
        ml_labels_free(&labels); free(x); free(y);
        return outt;
    }

    if (bayes) {
        /* Gaussian naive Bayes: per class, a mean and a per-feature variance, plus the
         * class prior. "Naive" is the independence assumption -- the joint density is the
         * PRODUCT of per-feature densities, i.e. a diagonal covariance -- which is why
         * this needs no Cholesky and why it works with far fewer points per class than the
         * full-covariance Multinormal does.
         *
         * THE VARIANCE FLOOR. A class whose feature j takes one value everywhere has zero
         * variance there and therefore infinite density at that value, which would make it
         * win every comparison involving that feature. Same class of problem as the
         * mixture's unbounded likelihood, so the same kind of fix -- but expressed as a
         * fraction of the feature's OVERALL variance across all classes rather than as a
         * fixed epsilon, so it is scale-invariant: a feature measured in millimetres and
         * the same feature in kilometres get proportionate floors, where a fixed epsilon
         * would be enormous for one and negligible for the other. */
        size_t C = labels.count;
        double* gmean = malloc(sizeof(double) * dim);
        double* gsd   = malloc(sizeof(double) * dim);
        double* prior = calloc(C, sizeof(double));
        double* cmean = calloc(C * dim, sizeof(double));
        double* cvar  = calloc(C * dim, sizeof(double));
        size_t* cnt   = calloc(C, sizeof(size_t));
        Expr* outb = NULL;
        if (gmean && gsd && prior && cmean && cvar && cnt) {
            ml_column_mean(x, n, dim, gmean);
            ml_column_sd(x, n, dim, gmean, gsd);

            for (size_t i = 0; i < n; i++) {
                cnt[y[i]]++;
                for (size_t j = 0; j < dim; j++) cmean[y[i] * dim + j] += x[i * dim + j];
            }
            for (size_t c = 0; c < C; c++) {
                prior[c] = (double)cnt[c] / (double)n;
                if (cnt[c] > 0)
                    for (size_t j = 0; j < dim; j++) cmean[c * dim + j] /= (double)cnt[c];
            }
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < dim; j++) {
                    double d = x[i * dim + j] - cmean[y[i] * dim + j];
                    cvar[y[i] * dim + j] += d * d;
                }
            for (size_t c = 0; c < C; c++)
                for (size_t j = 0; j < dim; j++) {
                    /* ML (n) divisor per class: with one point in a class the unbiased
                     * form would divide by zero, and the floor below is what makes the
                     * ML form safe. */
                    cvar[c * dim + j] = (cnt[c] > 0) ? cvar[c * dim + j] / (double)cnt[c]
                                                     : 0.0;
                    double floor_j = 1e-6 * gsd[j] * gsd[j];
                    if (!(floor_j > 0.0)) floor_j = 1e-12;   /* a constant feature */
                    if (!(cvar[c * dim + j] > floor_j)) cvar[c * dim + j] = floor_j;
                }

            /* Payload: vocabulary, priors, then a mean row and a variance row per class. */
            Expr** rows = malloc(sizeof(Expr*) * (2 + 2 * C));
            if (rows) {
                rows[0] = ml_labels_to_list(&labels);
                rows[1] = ml_list_of_reals(prior, C);
                for (size_t c = 0; c < C; c++) {
                    rows[2 + 2 * c]     = ml_list_of_reals(cmean + c * dim, dim);
                    rows[2 + 2 * c + 1] = ml_list_of_reals(cvar + c * dim, dim);
                }
                Expr* pay = expr_new_function(expr_new_symbol(SYM_List), rows, 2 + 2 * C);
                free(rows);
                if (pay) {
                    Expr* a4[4];
                    a4[0] = expr_new_string("NaiveBayes");
                    a4[1] = pay;
                    a4[2] = expr_new_integer((int64_t)dim);
                    a4[3] = expr_new_integer(0);
                    if (a4[0] && a4[1] && a4[2] && a4[3])
                        outb = expr_new_function(expr_new_symbol("ClassifierFunction"),
                                                 a4, 4);
                    else { expr_free(a4[0]); expr_free(a4[1]); expr_free(a4[2]); expr_free(a4[3]); }
                }
            }
        }
        free(gmean); free(gsd); free(prior); free(cmean); free(cvar); free(cnt);
        ml_labels_free(&labels);
        free(x); free(y);
        return outb;
    }

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
    size_t n = pay->data.function.arg_count - 1;   /* k-NN: one row per example */

    MlLabels labels = {0, NULL};
    if (!ml_labels_from_list(pay->data.function.args[0], &labels)) return NULL;

    if (argc == 1 && args[0] && args[0]->type == EXPR_STRING) {
        const char* q = args[0]->data.string;
        Expr* r = NULL;
        if (strcmp(q, "Method") == 0)            r = expr_copy(mname);
        else if (strcmp(q, "FeatureCount") == 0) r = expr_copy(dimx);
        else if (strcmp(q, "NeighborCount") == 0
                 && strcmp(mname->data.string, "NearestNeighbors") == 0) r = expr_copy(kx);
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

    if (strcmp(mname->data.string, "LogisticRegression") == 0) {
        /* TWO CLASSES: p = logistic(intercept + coef . x) is the probability of the SECOND
         * class in the vocabulary, the first being its complement, so the pair sums to 1 by
         * construction. That is why the tests assert the stronger properties instead --
         * exactly 0.5 on the fitted boundary, and monotone along the coefficient direction.
         *
         * MORE THAN TWO: K one-vs-rest fits, hence K independent sigmoids with nothing tying
         * them together. The class is the arg-max; the reported probabilities are those
         * sigmoids normalised to sum to 1. That normalisation is a presentation convention
         * and not a likelihood -- but being monotone it cannot move the arg-max, so the CLASS
         * is the trustworthy part and the tests lean on it hardest.
         *
         * The two shapes are told apart by the payload itself: a flat list of dim + 1 reals is
         * the two-class fit, a list of K such lists is one-vs-rest. No extra tag is stored,
         * because the shape already answers the question without ambiguity. */
        if (labels.count < 2 || pay->data.function.arg_count != 2) {
            free(xin); ml_labels_free(&labels); return NULL;
        }
        Expr* br = pay->data.function.args[1];
        if (!br || br->type != EXPR_FUNCTION) {
            free(xin); ml_labels_free(&labels); return NULL;
        }
        bool multi = labels.count != 2
                     && br->data.function.arg_count == labels.count
                     && br->data.function.args[0]
                     && br->data.function.args[0]->type == EXPR_FUNCTION;
        size_t nfit = multi ? labels.count : 1;
        if (!multi && br->data.function.arg_count != dim + 1) {
            free(xin); ml_labels_free(&labels); return NULL;
        }
        double* pv = malloc(sizeof(double) * nfit);
        if (!pv) { free(xin); ml_labels_free(&labels); return NULL; }
        for (size_t k = 0; k < nfit && ok; k++) {
            Expr* row = multi ? br->data.function.args[k] : br;
            if (!row || row->type != EXPR_FUNCTION
                || row->data.function.arg_count != dim + 1) { ok = false; break; }
            double e = 0.0, cf = 0.0;
            ok = na_read_scalar(row->data.function.args[0], &cf, &im) && im == 0.0;
            e = cf;
            for (size_t j2 = 0; j2 < dim && ok; j2++) {
                ok = na_read_scalar(row->data.function.args[j2 + 1], &cf, &im) && im == 0.0;
                if (ok) e += cf * xin[j2];
            }
            if (ok) pv[k] = 1.0 / (1.0 + exp(-e));
        }
        free(xin);
        if (!ok) { free(pv); ml_labels_free(&labels); return NULL; }
        Expr* outl = NULL;
        if (!multi) {
            double pr = pv[0];
            if (want_probs) {
                Expr** rules = malloc(sizeof(Expr*) * 2);
                if (rules) {
                    Expr* r0[2]; r0[0] = expr_copy(labels.label[0]); r0[1] = expr_new_real(1.0 - pr);
                    Expr* r1[2]; r1[0] = expr_copy(labels.label[1]); r1[1] = expr_new_real(pr);
                    rules[0] = expr_new_function(expr_new_symbol(SYM_Rule), r0, 2);
                    rules[1] = expr_new_function(expr_new_symbol(SYM_Rule), r1, 2);
                    outl = expr_new_function(expr_new_symbol(SYM_List), rules, 2);
                    free(rules);
                }
            } else {
                /* A probability of exactly 0.5 lies on the boundary; it goes to the FIRST
                 * class, matching the lowest-index tie-break the other methods use. */
                outl = expr_copy(labels.label[pr > 0.5 ? 1 : 0]);
            }
        } else if (want_probs) {
            /* Normalise the K sigmoids. Were every one of them to underflow to zero the
             * normaliser would be zero too, so that case DECLINES rather than dividing: a
             * uniform answer invented out of no information is worse than no answer. */
            double tot = 0.0;
            for (size_t k = 0; k < nfit; k++) tot += pv[k];
            if (!(tot > 0.0)) { free(pv); ml_labels_free(&labels); return NULL; }
            Expr** rules = malloc(sizeof(Expr*) * nfit);
            if (rules) {
                bool allr = true;
                for (size_t k = 0; k < nfit; k++) {
                    Expr* rr[2];
                    rr[0] = expr_copy(labels.label[k]);
                    rr[1] = expr_new_real(pv[k] / tot);
                    rules[k] = expr_new_function(expr_new_symbol(SYM_Rule), rr, 2);
                    if (!rules[k]) allr = false;
                }
                if (allr) {
                    outl = expr_new_function(expr_new_symbol(SYM_List), rules, nfit);
                } else {
                    for (size_t k = 0; k < nfit; k++) expr_free(rules[k]);
                }
                free(rules);
            }
        } else {
            /* Arg-max, ties to the lowest index -- the same tie-break as everywhere else
             * here, which is why this is a strict > rather than >=. */
            size_t best = 0;
            for (size_t k = 1; k < nfit; k++) if (pv[k] > pv[best]) best = k;
            outl = expr_copy(labels.label[best]);
        }
        free(pv);
        ml_labels_free(&labels);
        return outl;
    }

    if (strcmp(mname->data.string, "DecisionTree") == 0) {
        /* Walk from the root, then read the reached node's class histogram: arg-max for the
         * class, normalised for "Probabilities". One array answers both, which is the point of
         * storing counts at every node.
         *
         * The walk is over the payload directly rather than by rebuilding an MlTree. Routing a
         * single point touches at most depth nodes, so reconstructing the whole tree per
         * prediction would cost more than the prediction -- and the payload is already exactly
         * the structure the walk needs.
         *
         * The step count is bounded by the node count. A ClassifierFunction can be typed out
         * by hand, so left/right are untrusted input: a cycle would otherwise spin here
         * forever, and this runs on every prediction. */
        if (labels.count == 0 || pay->data.function.arg_count != 3) {
            free(xin); ml_labels_free(&labels); return NULL;
        }
        Expr* nodes = pay->data.function.args[1];
        Expr* dists = pay->data.function.args[2];
        if (!nodes || nodes->type != EXPR_FUNCTION || !dists || dists->type != EXPR_FUNCTION
            || nodes->data.function.arg_count == 0
            || nodes->data.function.arg_count != dists->data.function.arg_count) {
            free(xin); ml_labels_free(&labels); return NULL;
        }
        size_t nn = nodes->data.function.arg_count;

        size_t at = 0;
        for (size_t step = 0; step <= nn && ok; step++) {
            Expr* nd = nodes->data.function.args[at];
            if (!nd || nd->type != EXPR_FUNCTION || nd->data.function.arg_count != 4) {
                ok = false; break;
            }
            double fv = 0.0, th = 0.0, lv = 0.0, rv = 0.0;
            ok = na_read_scalar(nd->data.function.args[0], &fv, &im) && im == 0.0
              && na_read_scalar(nd->data.function.args[1], &th, &im) && im == 0.0
              && na_read_scalar(nd->data.function.args[2], &lv, &im) && im == 0.0
              && na_read_scalar(nd->data.function.args[3], &rv, &im) && im == 0.0;
            if (!ok) break;
            if (fv < 0.0) break;                        /* a leaf: this is the answer */
            size_t f = (size_t)fv;
            if (f >= dim) { ok = false; break; }
            size_t nxt = (xin[f] <= th) ? (size_t)lv : (size_t)rv;
            if (nxt >= nn || nxt == at) break;          /* malformed: stop where we are */
            at = nxt;
        }
        free(xin);
        if (!ok) { ml_labels_free(&labels); return NULL; }

        Expr* row = dists->data.function.args[at];
        if (!row || row->type != EXPR_FUNCTION
            || row->data.function.arg_count != labels.count) {
            ml_labels_free(&labels); return NULL;
        }
        double* h = malloc(sizeof(double) * labels.count);
        if (!h) { ml_labels_free(&labels); return NULL; }
        double tot = 0.0;
        for (size_t c = 0; c < labels.count && ok; c++) {
            ok = na_read_scalar(row->data.function.args[c], &h[c], &im) && im == 0.0;
            if (ok) tot += h[c];
        }
        if (!ok) { free(h); ml_labels_free(&labels); return NULL; }

        Expr* outt = NULL;
        if (want_probs) {
            /* An all-zero histogram would mean a node holding no training points, which
             * ml_tree_fit never produces. Declining beats dividing by zero, and beats
             * inventing a uniform answer from no evidence. */
            if (!(tot > 0.0)) { free(h); ml_labels_free(&labels); return NULL; }
            Expr** rules = malloc(sizeof(Expr*) * labels.count);
            if (rules) {
                bool allr = true;
                for (size_t c = 0; c < labels.count; c++) {
                    Expr* rr[2];
                    rr[0] = expr_copy(labels.label[c]);
                    rr[1] = expr_new_real(h[c] / tot);
                    rules[c] = expr_new_function(expr_new_symbol(SYM_Rule), rr, 2);
                    if (!rules[c]) allr = false;
                }
                if (allr) outt = expr_new_function(expr_new_symbol(SYM_List), rules, labels.count);
                else for (size_t c = 0; c < labels.count; c++) expr_free(rules[c]);
                free(rules);
            }
        } else {
            /* Arg-max, ties to the lowest class index -- the same tie-break the other three
             * methods use, hence the strict >. */
            size_t best = 0;
            for (size_t c = 1; c < labels.count; c++) if (h[c] > h[best]) best = c;
            outt = expr_copy(labels.label[best]);
        }
        free(h);
        ml_labels_free(&labels);
        return outt;
    }

    if (strcmp(mname->data.string, "NaiveBayes") == 0) {
        /* log posterior = log prior + sum over features of the Gaussian log density.
         *
         * Log space throughout, then one softmax. A product of dim densities underflows in
         * linear space for a point several standard deviations out in every feature -- the
         * same reason the mixture and the KDE work in logs -- and here it would silently
         * make every class equally (in)probable exactly where the answer is most obvious. */
        size_t C = labels.count;
        if (pay->data.function.arg_count != 2 + 2 * C) {
            free(xin); ml_labels_free(&labels); return NULL;
        }
        double* lp = malloc(sizeof(double) * C);
        if (!lp) { free(xin); ml_labels_free(&labels); return NULL; }
        Expr* prow = pay->data.function.args[1];
        ok = prow && prow->type == EXPR_FUNCTION && prow->data.function.arg_count == C;
        double best = -INFINITY;
        for (size_t c = 0; ok && c < C; c++) {
            double pr = 0.0;
            ok = na_read_scalar(prow->data.function.args[c], &pr, &im) && im == 0.0;
            if (!ok) break;
            /* A class with zero prior cannot win; log(0) would be -inf and poison the
             * softmax, so it is excluded rather than computed. */
            if (!(pr > 0.0)) { lp[c] = -INFINITY; continue; }
            double acc = log(pr);
            Expr* mr = pay->data.function.args[2 + 2 * c];
            Expr* vr = pay->data.function.args[2 + 2 * c + 1];
            if (!mr || mr->type != EXPR_FUNCTION || mr->data.function.arg_count != dim
                || !vr || vr->type != EXPR_FUNCTION || vr->data.function.arg_count != dim) {
                ok = false; break;
            }
            for (size_t j = 0; j < dim && ok; j++) {
                double mu = 0.0, va = 0.0;
                ok = na_read_scalar(mr->data.function.args[j], &mu, &im) && im == 0.0
                  && na_read_scalar(vr->data.function.args[j], &va, &im) && im == 0.0;
                if (!ok || !(va > 0.0)) { ok = false; break; }
                double z = xin[j] - mu;
                acc += -0.5 * (log(2.0 * M_PI * va) + z * z / va);
            }
            if (!ok) break;
            lp[c] = acc;
            if (acc > best) best = acc;
        }
        free(xin);
        if (!ok || !(best > -INFINITY)) {
            free(lp); ml_labels_free(&labels); return NULL;
        }
        Expr* outb = NULL;
        if (want_probs) {
            double sacc = 0.0;
            for (size_t c = 0; c < C; c++)
                if (lp[c] > -INFINITY) sacc += exp(lp[c] - best);
            Expr** rules = malloc(sizeof(Expr*) * C);
            if (rules) {
                for (size_t c = 0; c < C; c++) {
                    Expr* pr2[2];
                    pr2[0] = expr_copy(labels.label[c]);
                    pr2[1] = expr_new_real((lp[c] > -INFINITY)
                                           ? exp(lp[c] - best) / sacc : 0.0);
                    rules[c] = expr_new_function(expr_new_symbol(SYM_Rule), pr2, 2);
                }
                outb = expr_new_function(expr_new_symbol(SYM_List), rules, C);
                free(rules);
            }
        } else {
            size_t bc2 = 0;
            for (size_t c = 1; c < C; c++) if (lp[c] > lp[bc2]) bc2 = c;
            outb = expr_copy(labels.label[bc2]);
        }
        free(lp);
        ml_labels_free(&labels);
        return outb;
    }

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
        "\"NeighborCount\". Method -> \"NaiveBayes\" fits a Gaussian per class with a "
        "diagonal covariance; Method -> \"LogisticRegression\" fits a logistic model by "
        "iteratively reweighted least squares with a small ridge on the non-intercept "
        "coefficients -- the ridge is load-bearing, because on linearly separable data the "
        "unpenalised likelihood is unbounded and the coefficients would diverge. Two "
        "classes give a single fit; more than two are fitted one-vs-rest, one binary model "
        "per class, and the class is the arg-max of the fitted probabilities. Those "
        "probabilities are normalised to sum to 1, which is a convention rather than a "
        "likelihood -- being monotone it cannot change the arg-max, so the class is the "
        "better-founded of the two answers. A single class declines: it is not a "
        "classification problem.");

    symtab_get_def("ClassifierFunction")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ClassifierFunction",
        "ClassifierFunction[method, parameters, featureCount, k] is the fitted classifier "
        "Classify returns. Apply it to a feature vector to get a class, or to a feature "
        "vector and \"Probabilities\" to get one rule per class.");
}
