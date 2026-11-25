/* dist_formula.c -- closed forms of the specified distributions.
 *
 * NormalDistribution[mu, sigma] and UniformDistribution[{a, b}] are specified by
 * their parameters, so their properties are formulas in those parameters, exact
 * whenever the inputs are:
 *
 *                 NormalDistribution[m, s]            UniformDistribution[{a, b}]
 *   PDF[d, x]     E^(-(x - m)^2/(2 s^2))/(Sqrt[2 Pi] s)  Piecewise[{{1/(b - a), a <= x <= b}}, 0]
 *   CDF[d, x]     Erfc[(m - x)/(Sqrt[2] s)]/2          Piecewise[{{(x - a)/(b - a), a <= x <= b},
 *                                                                  {1, x > b}}, 0]
 *   Quantile[d,q] m - Sqrt[2] s InverseErfc[2 q]       a + (b - a) q
 *   Mean          m                                    (a + b)/2
 *   Median        m                                    (a + b)/2
 *   Variance      s^2                                  (b - a)^2/12
 *
 * StandardDeviation is Sqrt[Variance[d]] (stats), so it follows. Quantile is
 * answered only for a numeric q in [0, 1]. A numeric parameter set that is not
 * a distribution (s <= 0, b <= a) is declined, as the sampler does.
 *
 * The machine-number PDF of dist.c stays the fast path for approximate input;
 * this file answers exact and symbolic input, where a rounded density would be
 * the wrong answer (PDF[NormalDistribution[], 0] is 1/Sqrt[2 Pi], not 0.398942).
 *
 * Memory (SPEC section 4): arguments are borrowed; results are fresh.
 */
#include <stdlib.h>
#include <stdbool.h>

#include "expr.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "common.h"
#include "dist.h"

static Expr* mk1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), &a, 1);
}
static Expr* mk2(const char* h, Expr* a, Expr* b) {
    Expr* v[2] = { a, b };
    return expr_new_function(expr_new_symbol(h), v, 2);
}
static Expr* mk3(const char* h, Expr* a, Expr* b, Expr* c) {
    Expr* v[3] = { a, b, c };
    return expr_new_function(expr_new_symbol(h), v, 3);
}
static Expr* rat(int64_t p, int64_t q) { return mk2(SYM_Rational, expr_new_integer(p), expr_new_integer(q)); }
static Expr* neg(Expr* a) { return mk2(SYM_Times, expr_new_integer(-1), a); }
static Expr* sub(Expr* a, Expr* b) { return mk2(SYM_Plus, a, neg(b)); }
static Expr* inv(Expr* a) { return mk2(SYM_Power, a, expr_new_integer(-1)); }

/* True iff evaluating test gives True (consumes test). */
static bool holds(Expr* test) {
    Expr* r = eval_and_free(test);
    bool t = r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True;
    expr_free(r);
    return t;
}

typedef enum { FD_NONE, FD_NORMAL, FD_UNIFORM } FdKind;

/* The distribution's kind and parameters (borrowed into p[0..1]); FD_NONE when
 * dist is neither family or a numeric parameter set is not a distribution. */
static FdKind fd_read(const Expr* d, Expr** p0, Expr** p1) {
    if (!d || d->type != EXPR_FUNCTION || d->data.function.head->type != EXPR_SYMBOL)
        return FD_NONE;
    const char* h = d->data.function.head->data.symbol.name;
    size_t argc = d->data.function.arg_count;
    if (h == SYM_NormalDistribution) {
        if (argc == 0) { *p0 = expr_new_integer(0); *p1 = expr_new_integer(1); return FD_NORMAL; }
        if (argc != 2) return FD_NONE;
        Expr* s = d->data.function.args[1];
        if (expr_is_numeric_quantity(s) && !holds(mk1(SYM_Positive, expr_copy(s))))
            return FD_NONE;
        *p0 = expr_copy(d->data.function.args[0]);
        *p1 = expr_copy(s);
        return FD_NORMAL;
    }
    if (h == SYM_UniformDistribution) {
        if (argc == 0) { *p0 = expr_new_integer(0); *p1 = expr_new_integer(1); return FD_UNIFORM; }
        if (argc != 1) return FD_NONE;
        const Expr* r = d->data.function.args[0];
        if (!head_is(r, SYM_List) || r->data.function.arg_count != 2) return FD_NONE;
        Expr* a = r->data.function.args[0];
        Expr* b = r->data.function.args[1];
        if (expr_is_numeric_quantity(a) && expr_is_numeric_quantity(b) &&
            !holds(mk2(SYM_Less, expr_copy(a), expr_copy(b))))
            return FD_NONE;
        *p0 = expr_copy(a);
        *p1 = expr_copy(b);
        return FD_UNIFORM;
    }
    return FD_NONE;
}

/* The unevaluated formula for one scalar point x (consumed). */
static Expr* fd_point(MlDistProp prop, FdKind k, const Expr* p0, const Expr* p1, Expr* x) {
    Expr* m = expr_copy((Expr*)p0);
    Expr* s = expr_copy((Expr*)p1);
    if (k == FD_NORMAL) {
        if (prop == ML_PROP_PDF) {
            /* E^(-(x - m)^2/(2 s^2)) (2 Pi)^(-1/2) / s */
            Expr* ex = mk3(SYM_Times, rat(-1, 2),
                           mk2(SYM_Power, sub(x, m), expr_new_integer(2)),
                           mk2(SYM_Power, expr_copy(s), expr_new_integer(-2)));
            return mk3(SYM_Times, mk2(SYM_Power, expr_new_symbol(SYM_E), ex),
                       mk2(SYM_Power, mk2(SYM_Times, expr_new_integer(2), expr_new_symbol(SYM_Pi)),
                           rat(-1, 2)),
                       inv(s));
        }
        if (prop == ML_PROP_CDF) {
            /* Erfc[(m - x)/(Sqrt[2] s)]/2 */
            Expr* arg = mk3(SYM_Times, sub(m, x),
                            mk2(SYM_Power, expr_new_integer(2), rat(-1, 2)), inv(s));
            return mk2(SYM_Times, rat(1, 2), mk1(SYM_Erfc, arg));
        }
        /* Quantile: m - Sqrt[2] s InverseErfc[2 q] */
        Expr* t = mk3(SYM_Times, mk2(SYM_Power, expr_new_integer(2), rat(1, 2)), s,
                      mk1(SYM_InverseErfc, mk2(SYM_Times, expr_new_integer(2), x)));
        return sub(m, t);
    }
    /* Uniform on [a, b] = [m, s] */
    Expr* a = m;
    Expr* b = s;
    if (prop == ML_PROP_QUANTILE)
        return mk2(SYM_Plus, a, mk2(SYM_Times, sub(b, expr_copy(a)), x));
    Expr* inside = mk2(SYM_List,
        prop == ML_PROP_PDF ? inv(sub(expr_copy(b), expr_copy(a)))
                            : mk2(SYM_Times, sub(expr_copy(x), expr_copy(a)),
                                  inv(sub(expr_copy(b), expr_copy(a)))),
        mk3(SYM_LessEqual, expr_copy(a), expr_copy(x), expr_copy(b)));
    Expr* cases;
    if (prop == ML_PROP_PDF) {
        cases = mk1(SYM_List, inside);
    } else {
        Expr* above = mk2(SYM_List, expr_new_integer(1), mk2(SYM_Greater, expr_copy(x), expr_copy(b)));
        cases = mk2(SYM_List, inside, above);
    }
    expr_free(a); expr_free(b); expr_free(x);
    return mk2(SYM_Piecewise, cases, expr_new_integer(0));
}

/* A q Quantile can be read for: a number in [0, 1]. */
static bool fd_prob_ok(const Expr* q) {
    if (!expr_is_numeric_quantity(q)) return false;
    return holds(mk3(SYM_LessEqual, expr_new_integer(0), expr_copy((Expr*)q), expr_new_integer(1)));
}

Expr* ml_dist_closed_form(MlDistProp prop, const Expr* dist, const Expr* x) {
    Expr *p0 = NULL, *p1 = NULL;
    FdKind k = fd_read(dist, &p0, &p1);
    if (k == FD_NONE) return NULL;
    Expr* out = NULL;
    switch (prop) {
        case ML_PROP_MEAN:
        case ML_PROP_MEDIAN:
            out = (k == FD_NORMAL) ? expr_copy(p0)
                : mk2(SYM_Times, rat(1, 2), mk2(SYM_Plus, expr_copy(p0), expr_copy(p1)));
            break;
        case ML_PROP_VARIANCE:
            out = (k == FD_NORMAL) ? mk2(SYM_Power, expr_copy(p1), expr_new_integer(2))
                : mk2(SYM_Times, rat(1, 12),
                      mk2(SYM_Power, sub(expr_copy(p1), expr_copy(p0)), expr_new_integer(2)));
            break;
        case ML_PROP_PDF:
        case ML_PROP_CDF:
        case ML_PROP_QUANTILE: {
            if (!x) break;
            /* A list of points threads (the univariate families only). */
            if (head_is(x, SYM_List)) {
                size_t n = x->data.function.arg_count;
                Expr** el = malloc(sizeof(Expr*) * (n ? n : 1));
                bool ok = el != NULL;
                for (size_t i = 0; ok && i < n; i++) {
                    const Expr* xi = x->data.function.args[i];
                    if (prop == ML_PROP_QUANTILE && !fd_prob_ok(xi)) {
                        for (size_t j = 0; j < i; j++) expr_free(el[j]);
                        ok = false;
                        break;
                    }
                    el[i] = fd_point(prop, k, p0, p1, expr_copy((Expr*)xi));
                }
                if (ok) out = expr_new_function(expr_new_symbol(SYM_List), el, n);
                free(el);
                break;
            }
            if (prop == ML_PROP_QUANTILE && !fd_prob_ok(x)) break;
            out = fd_point(prop, k, p0, p1, expr_copy((Expr*)x));
            break;
        }
    }
    expr_free(p0);
    expr_free(p1);
    return out ? eval_and_free(out) : NULL;
}
