/* integrate_line.c
 *
 * Complex line (contour) integration -- see integrate_line.h for the overview.
 *
 * A contour integral Integrate[f, {x, z0, z1, ..., zn}] is the sum of the
 * straight-segment integrals z_k -> z_{k+1}.  Each segment is handled by
 * `integrate_line_segment`, which:
 *
 *   1. F = Integrate[f, x]                 -- antiderivative in the friendly
 *      variable x (rational/real coefficients).  Bail if it cannot be found.
 *
 *   2. Detect on-path singularities by parametrising the segment as
 *      gamma(t) = a + t (b - a), t in [0, 1], and locating the real roots
 *      t* in (0, 1) of Denominator[Together[f(gamma(t))]].  Such a root is a
 *      singularity strictly on the segment: the contour integral is divergent
 *      (Integrate::idiv), left unevaluated.
 *
 *   3. Evaluate the continuous change of F along the segment.  The primary
 *      method is direct: F(b) - F(a), with endpoint values taken by the Limit
 *      engine as real one-sided limits in t when substitution is singular (so a
 *      complex-ray approach is reduced to a real one-sided limit).  When the
 *      segment crosses a branch cut of F the direct value can be off by the
 *      branch jump; for rational integrands whose antiderivative is a sum of
 *      logarithms/inverse-tangents of AFFINE arguments, the branch-correct value
 *      is recovered by combining each log into a single principal Log of a ratio
 *      Log[(u(b))/(u(a))] -- exact because a straight segment subtends an angle
 *      < pi at any point off it.
 *
 *   4. Cross-check the symbolic value against a complex quadrature of
 *      f(gamma(t)) gamma'(t) over t in [0, 1].  Accept a symbolic candidate only
 *      when it agrees; an uncorrectable branch crossing leaves the integral
 *      unevaluated rather than return a wrong branch.
 */

#include "integrate_line.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "numeric.h"
#include "sym_names.h"
#include "print.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* -------------------------------------------------------------------------
 * Small expression-construction / evaluation helpers (mirroring the
 * Newton-Leibniz driver -- the useful ones there are file-static).
 * ---------------------------------------------------------------------- */

static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }

static Expr* mk_fn1(const char* head, Expr* a) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* mk_fn3(const char* head, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b, c }, 3);
}

/* Evaluate `call`, free the call expression, return the result. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* True iff `e` is the compound `name[...]` (by head name, not interned ptr). */
static bool head_name_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol, name) == 0;
}

/* True iff any subexpression of `e` is a call with head `name`. */
static bool contains_head(const Expr* e, const char* name) {
    if (!e) return false;
    if (head_name_is(e, name)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_head(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_head(e->data.function.args[i], name)) return true;
    return false;
}

/* True iff the symbol `x` occurs anywhere in `e` (by interned pointer). */
static bool contains_symbol(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol == x->data.symbol;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol(e->data.function.args[i], x)) return true;
    return false;
}

/* True iff a symbol with textual name `name` occurs anywhere in `e`. */
static bool has_symbol_name(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol, name) == 0;
    if (e->type != EXPR_FUNCTION) return false;
    if (has_symbol_name(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_symbol_name(e->data.function.args[i], name)) return true;
    return false;
}

/* -------------------------------------------------------------------------
 * Numeric extraction.
 * ---------------------------------------------------------------------- */

/* Extract a machine double from an already-numeric real leaf. */
static bool line_real_double(const Expr* e, double* out) {
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;   return true;
        case EXPR_REAL:    *out = e->data.real;              return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }
    int64_t n, d;
    if (is_rational((Expr*)e, &n, &d) && d != 0) {
        *out = (double)n / (double)d;
        return true;
    }
    return false;
}

/* Numericalize `e` to a machine complex (re, im).  Returns false if the result
 * is not a finite complex number (symbolic, infinite, or indeterminate). */
static bool line_numeric_complex(Expr* e, double* re, double* im) {
    if (!e) return false;
    if (is_infinity_sym(e) || is_neg_infinity_form(e) ||
        is_complex_infinity_sym(e) || is_indeterminate_sym(e))
        return false;

    Expr* n = numericalize(e, numeric_machine_spec());
    if (!n) return false;

    bool ok = false;
    Expr* rp = NULL; Expr* ip = NULL;
    if (is_complex(n, &rp, &ip)) {
        double r, i;
        if (line_real_double(rp, &r) && line_real_double(ip, &i)) {
            *re = r; *im = i; ok = true;
        }
    } else if (line_real_double(n, re)) {
        *im = 0.0; ok = true;
    }
    expr_free(n);
    if (ok && (!isfinite(*re) || !isfinite(*im))) ok = false;
    return ok;
}

/* -------------------------------------------------------------------------
 * Singular-sentinel classification of boundary values.
 * ---------------------------------------------------------------------- */

static bool line_is_infinite(Expr* e) {
    return is_infinity_sym(e) || is_neg_infinity_form(e) ||
           is_complex_infinity_sym(e);
}

/* True iff any subexpression of `e` is a divergence sentinel. */
static bool line_contains_singular(const Expr* e) {
    if (!e) return false;
    if (line_is_infinite((Expr*)e) || is_indeterminate_sym((Expr*)e)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (line_contains_singular(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (line_contains_singular(e->data.function.args[i])) return true;
    return false;
}

/* A boundary value is usable directly only if it is a finite closed form. */
static bool line_is_finite_closed(Expr* v, const Expr* x) {
    if (!v) return false;
    if (line_contains_singular(v)) return false;
    if (head_name_is(v, "Interval")) return false;
    if (contains_head(v, "Limit") || contains_head(v, "Integrate")) return false;
    if (contains_symbol(v, x)) return false;
    return true;
}

/* -------------------------------------------------------------------------
 * Path parametrisation and substitution.
 * ---------------------------------------------------------------------- */

/* ReplaceAll[e, x -> val] then evaluate.  Borrows all; returns a new Expr. */
static Expr* subst_eval(Expr* e, Expr* x, Expr* val) {
    return eval_take(mk_fn2("ReplaceAll", expr_copy(e),
                            mk_fn2("Rule", expr_copy(x), expr_copy(val))));
}

/* gamma(t) = a + t (b - a), with the parameter symbol named `tname`. */
static Expr* line_make_gamma(Expr* a, Expr* b, const char* tname) {
    Expr* bma = mk_fn2("Plus", expr_copy(b),
                       mk_fn2("Times", expr_new_integer(-1), expr_copy(a)));
    Expr* g = mk_fn2("Plus", expr_copy(a),
                     mk_fn2("Times", mk_sym(tname), bma));
    return eval_take(g);
}

/* A parameter-symbol name guaranteed not to occur in f, a, or b. */
static char* line_fresh_param(Expr* f, Expr* a, Expr* b) {
    char buf[48];
    for (int n = 0; n < 10000; n++) {
        snprintf(buf, sizeof(buf), "$lineparam%d$", n);
        if (!has_symbol_name(f, buf) && !has_symbol_name(a, buf) &&
            !has_symbol_name(b, buf)) {
            size_t len = strlen(buf);
            char* out = malloc(len + 1);
            memcpy(out, buf, len + 1);
            return out;
        }
    }
    const char* fb = "$lineparam$";
    size_t len = strlen(fb);
    char* out = malloc(len + 1);
    memcpy(out, fb, len + 1);
    return out;
}

typedef enum { LINE_DIR_ABOVE, LINE_DIR_BELOW } LineDir;

/* Evaluate the antiderivative F at a segment endpoint `endpoint` (= gamma(tv),
 * tv in {0, 1}).  Direct substitution when it yields a finite closed form;
 * otherwise a real one-sided limit of F(gamma(t)) as t -> tv along the segment
 * (Direction "FromAbove" at t=0, "FromBelow" at t=1). */
static Expr* line_eval_endpoint(Expr* F, Expr* x, Expr* gamma, const char* tname,
                                Expr* endpoint, double tv, LineDir dir) {
    arith_warnings_mute_push();
    Expr* sub = subst_eval(F, x, endpoint);
    Expr* out;
    if (line_is_finite_closed(sub, x)) {
        out = sub;
    } else {
        if (sub) expr_free(sub);
        Expr* G = subst_eval(F, x, gamma);        /* F as a function of t */
        Expr* t = mk_sym(tname);
        Expr* rule = mk_fn2("Rule", t, expr_new_real(tv));
        Expr* d = mk_fn2("Rule", mk_sym("Direction"),
                         expr_new_string(dir == LINE_DIR_ABOVE ? "FromAbove"
                                                               : "FromBelow"));
        out = eval_take(mk_fn3("Limit", G, rule, d)); /* G consumed */
    }
    arith_warnings_mute_pop();
    return out;
}

/* -------------------------------------------------------------------------
 * On-path singularity detection.
 * ---------------------------------------------------------------------- */

/* The RHS of the (single) `Rule[t, val]` inside a Solve solution element. */
static Expr* line_rule_rhs_for(Expr* el, const Expr* t) {
    if (head_name_is(el, "Rule") && el->data.function.arg_count == 2) {
        Expr* lhs = el->data.function.args[0];
        if (lhs->type == EXPR_SYMBOL && lhs->data.symbol == t->data.symbol)
            return el->data.function.args[1];
        return NULL;
    }
    if (head_name_is(el, "List")) {
        for (size_t i = 0; i < el->data.function.arg_count; i++) {
            Expr* rr = line_rule_rhs_for(el->data.function.args[i], t);
            if (rr) return rr;
        }
    }
    return NULL;
}

/* Numeric grid scan of |Denominator[f(gamma(t))]| across (0, 1) for a genuine
 * interior zero (both real and imaginary parts vanish).  Returns
 *    1 = an interior zero found,   0 = none found,   -1 = cannot scan. */
static int line_denominator_grid_scan(Expr* den, const char* tname) {
    Expr* t = mk_sym(tname);
    const int N = 512;
    int result = 0;
    double prev_re = 0.0, prev_im = 0.0;
    bool have_prev = false;
    for (int i = 0; i <= N && result == 0; i++) {
        double tv = (double)i / (double)N;
        /* Substitute t -> tv and numericalize. */
        Expr* sub = eval_take(mk_fn2("ReplaceAll", expr_copy(den),
                                     mk_fn2("Rule", expr_copy(t),
                                            expr_new_real(tv))));
        double re, im;
        bool ok = line_numeric_complex(sub, &re, &im);
        if (sub) expr_free(sub);
        if (!ok) { result = -1; break; }
        bool interior = (i != 0 && i != N);
        double mod = hypot(re, im);
        if (interior && mod < 1e-9) { result = 1; break; }
        if (have_prev && (prev_re * re < 0.0) && (prev_im * im < 0.0)) {
            result = 1; break;   /* both components changed sign together */
        }
        prev_re = re; prev_im = im; have_prev = true;
    }
    expr_free(t);
    return result;
}

/* Collect the on-segment singular points of f (as z-values gamma(t*)) for the
 * straight segment a -> b.  Appends owned Expr* z-points to *pts / *n (growing
 * *cap).  Returns:
 *    1 = at least one interior singularity found,
 *    0 = none,
 *   -1 = detection undecidable (symbolic root, or unscannable). */
static int line_segment_singularities(Expr* f, Expr* x, Expr* a, Expr* b,
                                      Expr*** pts, size_t* n, size_t* cap) {
    if (expr_eq(a, b)) return 0;                       /* degenerate segment */

    char* tname = line_fresh_param(f, a, b);
    Expr* gamma = line_make_gamma(a, b, tname);
    Expr* t = mk_sym(tname);

    Expr* fg  = subst_eval(f, x, gamma);               /* f(gamma(t)) */
    Expr* tog = fg ? eval_take(mk_fn1("Together", fg)) : NULL;  /* fg consumed */
    Expr* den = tog ? eval_take(mk_fn1("Denominator", tog)) : NULL;

    int status = 0;
    if (!den || !contains_symbol(den, t)) {
        status = 0;                                    /* no t in denominator */
    } else {
        Expr* eq   = mk_fn2("Equal", expr_copy(den), expr_new_integer(0));
        Expr* sols = eval_take(mk_fn3("Solve", eq, expr_copy(t), mk_sym("Reals")));
        bool undecidable = false;
        bool found = false;
        if (sols && head_name_is(sols, "List")) {
            for (size_t i = 0; i < sols->data.function.arg_count; i++) {
                Expr* val = line_rule_rhs_for(sols->data.function.args[i], t);
                if (!val) { undecidable = true; continue; }
                double re, im;
                if (!line_numeric_complex(val, &re, &im)) { undecidable = true; continue; }
                if (fabs(im) > 1e-12 * (1.0 + fabs(re))) continue; /* non-real t */
                /* Interior of the open segment: 0 < t* < 1 (endpoints handled
                 * separately as improper-integral limits). */
                if (re <= 1e-12 || re >= 1.0 - 1e-12) continue;
                Expr* z = subst_eval(gamma, t, val);   /* gamma(t*) */
                if (!z) continue;
                if (*n == *cap) {
                    size_t nc = *cap ? *cap * 2 : 4;
                    *pts = realloc(*pts, nc * sizeof(**pts));
                    *cap = nc;
                }
                (*pts)[(*n)++] = z;
                found = true;
            }
        } else {
            undecidable = true;
        }
        if (sols) expr_free(sols);

        if (found)            status = 1;
        else if (undecidable) status = line_denominator_grid_scan(den, tname);
        else                  status = 0;
    }

    if (den) expr_free(den);
    expr_free(gamma);
    expr_free(t);
    free(tname);
    return status;
}

/* -------------------------------------------------------------------------
 * Branch-correct evaluation for rational integrands (affine Log / ArcTan).
 * ---------------------------------------------------------------------- */

/* A head that introduces a branch cut when applied to something containing x. */
static bool is_branch_head(const char* h) {
    static const char* const H[] = {
        "Log", "Sqrt", "ArcTan", "ArcCot", "ArcTanh", "ArcCoth",
        "ArcSin", "ArcCos", "ArcSec", "ArcCsc", "ArcSinh", "ArcCosh"
    };
    for (size_t i = 0; i < sizeof H / sizeof *H; i++)
        if (strcmp(h, H[i]) == 0) return true;
    return false;
}

/* True iff `e` contains a multivalued (branch-cut) subexpression involving x:
 * a branch head over x, or a Power with a non-integer exponent whose base holds
 * x.  Entire heads (Exp, Sin, Cos, polynomials) are single-valued and pass. */
static bool has_branch_over_x(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (is_branch_head(h) && contains_symbol(e, x)) return true;
        if (strcmp(h, "Power") == 0 && e->data.function.arg_count == 2) {
            Expr* base = e->data.function.args[0];
            Expr* exp  = e->data.function.args[1];
            if (contains_symbol(base, x)) {
                /* Non-integer exponent -> branch (Sqrt, cube roots, ...). */
                if (exp->type != EXPR_INTEGER && exp->type != EXPR_BIGINT)
                    return true;
            }
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_branch_over_x(e->data.function.args[i], x)) return true;
    return false;
}

/* True iff u is affine in x: u actually contains x and D[u, x] is free of x. */
static bool is_affine_in(Expr* u, Expr* x) {
    if (!contains_symbol(u, x)) return false;
    Expr* d = eval_take(mk_fn2("D", expr_copy(u), expr_copy(x)));
    bool ok = d && !contains_symbol(d, x);
    if (d) expr_free(d);
    return ok;
}

/* Principal Log of the ratio: Log[(u/.x->b)/(u/.x->a)] -- the exact continuous
 * change of Log[u] along the straight segment a->b for AFFINE u. */
static Expr* line_log_ratio(Expr* u, Expr* x, Expr* a, Expr* b) {
    Expr* ub = subst_eval(u, x, b);
    Expr* ua = subst_eval(u, x, a);
    if (!ub || !ua) { if (ub) expr_free(ub); if (ua) expr_free(ua); return NULL; }
    Expr* ratio = mk_fn2("Times", ub, mk_fn2("Power", ua, expr_new_integer(-1)));
    return eval_take(mk_fn1("Log", ratio));
}

/* Continuous change of a single term of F along a->b, or NULL if the term is
 * not branch-safely reducible (the caller then falls back to direct + check). */
static Expr* line_term_delta(Expr* term, Expr* x, Expr* a, Expr* b) {
    if (!contains_symbol(term, x)) return expr_new_integer(0);

    /* Single-valued (rational / entire) term: the plain difference is exact. */
    if (!has_branch_over_x(term, x)) {
        Expr* tb = subst_eval(term, x, b);
        Expr* ta = subst_eval(term, x, a);
        if (!tb || !ta) { if (tb) expr_free(tb); if (ta) expr_free(ta); return NULL; }
        return eval_take(mk_fn2("Plus", tb,
                                mk_fn2("Times", expr_new_integer(-1), ta)));
    }

    /* Otherwise the term must be  coeff * H[u]  with coeff free of x, H one of
     * Log / ArcTan / ArcTanh, and u affine in x. */
    Expr** factors; size_t nf;
    Expr* single[1];
    if (head_name_is(term, "Times")) {
        factors = term->data.function.args;
        nf = term->data.function.arg_count;
    } else {
        single[0] = term; factors = single; nf = 1;
    }

    Expr* coeff = expr_new_integer(1);
    Expr* core = NULL;                      /* the one x-dependent factor */
    bool bad = false;
    for (size_t i = 0; i < nf && !bad; i++) {
        Expr* fac = factors[i];
        if (!contains_symbol(fac, x)) {
            coeff = eval_take(mk_fn2("Times", coeff, expr_copy(fac)));
        } else if (!core) {
            core = fac;
        } else {
            bad = true;                     /* more than one x-dependent factor */
        }
    }
    if (bad || !core) { expr_free(coeff); return NULL; }

    /* core must be H[u] with H a supported branch head and u affine. */
    if (core->type != EXPR_FUNCTION || core->data.function.arg_count != 1 ||
        core->data.function.head->type != EXPR_SYMBOL) {
        expr_free(coeff); return NULL;
    }
    const char* h = core->data.function.head->data.symbol;
    Expr* u = core->data.function.args[0];
    if (!is_affine_in(u, x)) { expr_free(coeff); return NULL; }

    Expr* delta = NULL;
    if (strcmp(h, "Log") == 0) {
        delta = line_log_ratio(u, x, a, b);
    } else if (strcmp(h, "ArcTan") == 0) {
        /* ArcTan[u] = (1/(2I)) ( Log[1+I u] - Log[1-I u] ), 1 +/- I u affine.
         * delta = (-I/2)( Log[(1+I u_b)/(1+I u_a)] - Log[(1-I u_b)/(1-I u_a)] ). */
        Expr* iu = mk_fn2("Times", mk_sym(SYM_I), expr_copy(u));
        Expr* up = eval_take(mk_fn2("Plus", expr_new_integer(1), expr_copy(iu)));
        Expr* um = eval_take(mk_fn2("Plus", expr_new_integer(1),
                                    mk_fn2("Times", expr_new_integer(-1), iu)));
        Expr* lp = up ? line_log_ratio(up, x, a, b) : NULL;
        Expr* lm = um ? line_log_ratio(um, x, a, b) : NULL;
        if (up) expr_free(up);
        if (um) expr_free(um);
        if (lp && lm) {
            Expr* diff = mk_fn2("Plus", lp,
                                mk_fn2("Times", expr_new_integer(-1), lm));
            Expr* half_i = mk_fn2("Times", mk_fn2("Power", expr_new_integer(2),
                                                  expr_new_integer(-1)),
                                  mk_sym(SYM_I));   /* I/2 */
            /* (-I/2) = -1 * (I/2) */
            delta = eval_take(mk_fn2("Times",
                        mk_fn2("Times", expr_new_integer(-1), half_i), diff));
        } else {
            if (lp) expr_free(lp);
            if (lm) expr_free(lm);
        }
    } else if (strcmp(h, "ArcTanh") == 0) {
        /* ArcTanh[u] = (1/2)( Log[1+u] - Log[1-u] ), 1 +/- u affine.
         * delta = (1/2)( Log[(1+u_b)/(1+u_a)] - Log[(1-u_b)/(1-u_a)] ). */
        Expr* up = eval_take(mk_fn2("Plus", expr_new_integer(1), expr_copy(u)));
        Expr* um = eval_take(mk_fn2("Plus", expr_new_integer(1),
                                    mk_fn2("Times", expr_new_integer(-1),
                                           expr_copy(u))));
        Expr* lp = up ? line_log_ratio(up, x, a, b) : NULL;
        Expr* lm = um ? line_log_ratio(um, x, a, b) : NULL;
        if (up) expr_free(up);
        if (um) expr_free(um);
        if (lp && lm) {
            Expr* diff = mk_fn2("Plus", lp,
                                mk_fn2("Times", expr_new_integer(-1), lm));
            delta = eval_take(mk_fn2("Times",
                        mk_fn2("Power", expr_new_integer(2), expr_new_integer(-1)),
                        diff));
        } else {
            if (lp) expr_free(lp);
            if (lm) expr_free(lm);
        }
    }

    if (!delta) { expr_free(coeff); return NULL; }
    return eval_take(mk_fn2("Times", coeff, delta));
}

/* Branch-correct continuous change of F along the straight segment a->b, or
 * NULL if F is not reducible to affine Log / ArcTan / rational parts. */
static Expr* line_branch_combine(Expr* F, Expr* x, Expr* a, Expr* b) {
    Expr** terms; size_t nt;
    Expr* single[1];
    if (head_name_is(F, "Plus")) {
        terms = F->data.function.args;
        nt = F->data.function.arg_count;
    } else {
        single[0] = F; terms = single; nt = 1;
    }

    Expr* total = expr_new_integer(0);
    for (size_t i = 0; i < nt; i++) {
        Expr* d = line_term_delta(terms[i], x, a, b);
        if (!d) { expr_free(total); return NULL; }
        total = eval_take(mk_fn2("Plus", total, d));
        if (!total) return NULL;
    }
    return total;
}

/* -------------------------------------------------------------------------
 * Numeric cross-check.
 * ---------------------------------------------------------------------- */

/* Numeric value of Integrate[f, {x, a, b}] along the segment, by complex
 * quadrature of f(gamma(t)) (b - a) over t in [0, 1] (NIntegrate). */
static bool line_numeric_value(Expr* f, Expr* x, Expr* a, Expr* b,
                               double* re, double* im) {
    char* tname = line_fresh_param(f, a, b);
    Expr* gamma = line_make_gamma(a, b, tname);
    Expr* t = mk_sym(tname);

    Expr* fg  = subst_eval(f, x, gamma);
    Expr* bma = eval_take(mk_fn2("Plus", expr_copy(b),
                                 mk_fn2("Times", expr_new_integer(-1),
                                        expr_copy(a))));
    Expr* integrand = NULL;
    if (fg && bma) {
        integrand = eval_take(mk_fn2("Times", fg, bma));
    } else {
        if (fg)  expr_free(fg);
        if (bma) expr_free(bma);
    }

    bool ok = false;
    if (integrand) {
        arith_warnings_mute_push();
        Expr* spec = mk_fn3("List", expr_copy(t), expr_new_integer(0),
                            expr_new_integer(1));
        Expr* ni = eval_take(mk_fn2("NIntegrate", integrand, spec)); /* consumed */
        arith_warnings_mute_pop();
        if (ni) {
            ok = line_numeric_complex(ni, re, im);
            expr_free(ni);
        }
    }

    expr_free(gamma);
    expr_free(t);
    free(tname);
    return ok;
}

/* True iff the symbolic value `V` matches the numeric segment value (rn, in)
 * to a loose relative tolerance -- large branch jumps (O(2 Pi)) are rejected
 * while quadrature noise passes. */
static bool line_value_matches(Expr* V, double rn, double in) {
    double rv, iv;
    if (!line_numeric_complex(V, &rv, &iv)) return false;
    double err = hypot(rv - rn, iv - in);
    double mag = hypot(rv, iv);
    return err <= 1e-3 * (1.0 + mag);
}

/* -------------------------------------------------------------------------
 * Core single-segment driver.
 * ---------------------------------------------------------------------- */

static void line_emit_idiv(Expr* f, Expr* a, Expr* b) {
    char* fs = expr_to_string(f);
    char* as = expr_to_string(a);
    char* bs = expr_to_string(b);
    fprintf(stderr,
            "Integrate::idiv: Integral of %s does not converge on the contour "
            "{%s, %s}.\n", fs ? fs : "?", as ? as : "?", bs ? bs : "?");
    free(fs); free(as); free(bs);
}

/* Contour integral of f along the single straight segment a -> b. */
static Expr* integrate_line_segment(Expr* f, Expr* x, Expr* a, Expr* b,
                                    const char* method) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (expr_eq(a, b)) return expr_new_integer(0);     /* zero-length segment */

    /* 1. Antiderivative in x. */
    Expr* F;
    if (method) {
        F = eval_take(mk_fn3("Integrate", expr_copy(f), expr_copy(x),
                             mk_fn2("Rule", mk_sym("Method"),
                                    expr_new_string(method))));
    } else {
        F = eval_take(mk_fn2("Integrate", expr_copy(f), expr_copy(x)));
    }
    if (!F) return NULL;
    if (contains_head(F, "Integrate")) { expr_free(F); return NULL; }

    /* 2. On-path singularity: a genuine interior singularity makes the contour
     * integral divergent (undefined). */
    Expr** spts = NULL; size_t sn = 0, scap = 0;
    int scan = line_segment_singularities(f, x, a, b, &spts, &sn, &scap);
    for (size_t i = 0; i < sn; i++) expr_free(spts[i]);
    free(spts);
    if (scan == 1) {
        line_emit_idiv(f, a, b);
        expr_free(F);
        return NULL;
    }

    /* 3a. Numeric segment value for the cross-check. */
    double rn = 0.0, in = 0.0;
    bool have_num = line_numeric_value(f, x, a, b, &rn, &in);

    /* 3b. Direct continuous change F(b) - F(a) (endpoint limits along the path
     * when substitution is singular). */
    char* tname = line_fresh_param(f, a, b);
    Expr* gamma = line_make_gamma(a, b, tname);
    Expr* Fb = line_eval_endpoint(F, x, gamma, tname, b, 1.0, LINE_DIR_BELOW);
    Expr* Fa = line_eval_endpoint(F, x, gamma, tname, a, 0.0, LINE_DIR_ABOVE);
    expr_free(gamma);
    free(tname);

    Expr* direct = NULL;
    if (Fb && Fa) {
        arith_warnings_mute_push();
        direct = eval_take(mk_fn2("Plus", Fb,
                                  mk_fn2("Times", expr_new_integer(-1), Fa)));
        arith_warnings_mute_pop();
    } else {
        if (Fb) expr_free(Fb);
        if (Fa) expr_free(Fa);
    }

    /* Divergence surfaced by an endpoint blow-up. */
    if (direct && (line_is_infinite(direct) || is_indeterminate_sym(direct))) {
        line_emit_idiv(f, a, b);
        expr_free(direct);
        expr_free(F);
        return NULL;
    }

    /* Accept the direct value when it is a clean closed form that the numeric
     * check confirms (or, if quadrature is unavailable, only when singularity
     * detection was clean). */
    if (direct && line_is_finite_closed(direct, x)) {
        if (have_num) {
            if (line_value_matches(direct, rn, in)) { expr_free(F); return direct; }
        } else if (scan == 0) {
            expr_free(F);
            return direct;
        }
    }
    if (direct) expr_free(direct);

    /* 3c. Branch-correct fallback for rational integrands (affine Log/ArcTan). */
    Expr* branch = line_branch_combine(F, x, a, b);
    expr_free(F);
    if (branch) {
        if (line_is_infinite(branch) || is_indeterminate_sym(branch)) {
            line_emit_idiv(f, a, b);
            expr_free(branch);
            return NULL;
        }
        if (line_is_finite_closed(branch, x) && have_num &&
            line_value_matches(branch, rn, in)) {
            return branch;
        }
        expr_free(branch);
    }

    return NULL;   /* unverifiable / unevaluable -> leave unevaluated */
}

/* -------------------------------------------------------------------------
 * Contour driver + spec parsing.
 * ---------------------------------------------------------------------- */

Expr* integrate_line_contour(Expr* f, Expr* x, Expr** pts, size_t npts,
                             const char* method) {
    if (!f || !x || !pts || npts < 2 || x->type != EXPR_SYMBOL) return NULL;

    Expr* total = expr_new_integer(0);
    for (size_t k = 0; k + 1 < npts; k++) {
        Expr* seg = integrate_line_segment(f, x, pts[k], pts[k + 1], method);
        if (!seg) { expr_free(total); return NULL; }
        total = eval_take(mk_fn2("Plus", total, seg));
        if (!total) return NULL;
    }
    return total;
}

bool integrate_line_is_contour_spec(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == SYM_List &&
           e->data.function.arg_count >= 3 &&
           e->data.function.args[0]->type == EXPR_SYMBOL;
}

/* True iff an endpoint is a non-real complex number (a Complex[..] node -- the
 * evaluator collapses Complex[a, 0] back to a, so any surviving Complex head
 * carries a non-zero imaginary part) or a symbolic expression containing the
 * imaginary unit I. */
static bool endpoint_is_nonreal(Expr* e) {
    Expr* re = NULL; Expr* im = NULL;
    if (is_complex(e, &re, &im)) return true;
    if (has_symbol_name(e, "I")) return true;
    return false;
}

bool integrate_line_spec_is_complex(const Expr* spec) {
    if (!integrate_line_is_contour_spec(spec)) return false;
    size_t argc = spec->data.function.arg_count;
    if (argc > 3) return true;                          /* polyline contour */
    for (size_t i = 1; i < argc; i++)
        if (endpoint_is_nonreal(spec->data.function.args[i])) return true;
    return false;
}

Expr* integrate_line_from_spec(Expr* f, Expr* spec, const char* method) {
    if (!integrate_line_is_contour_spec(spec)) return NULL;
    size_t argc = spec->data.function.arg_count;
    Expr* x = spec->data.function.args[0];
    size_t npts = argc - 1;
    Expr** pts = spec->data.function.args + 1;          /* borrowed */
    return integrate_line_contour(f, x, pts, npts, method);
}

/* -------------------------------------------------------------------------
 * Builtins.
 * ---------------------------------------------------------------------- */

Expr* builtin_integrate_line(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* f = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!integrate_line_is_contour_spec(spec)) return NULL;
    return integrate_line_from_spec(f, spec, NULL);
}

Expr* builtin_integrate_path_singular_points(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* f = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!integrate_line_is_contour_spec(spec)) return NULL;

    Expr* x = spec->data.function.args[0];
    size_t argc = spec->data.function.arg_count;

    Expr** pts = NULL; size_t n = 0, cap = 0;
    for (size_t k = 1; k + 1 < argc; k++) {
        line_segment_singularities(f, x, spec->data.function.args[k],
                                   spec->data.function.args[k + 1],
                                   &pts, &n, &cap);
    }
    /* Sort canonically and dedup. */
    for (size_t i = 0; i + 1 < n; i++)
        for (size_t j = i + 1; j < n; j++)
            if (expr_compare(pts[i], pts[j]) > 0) {
                Expr* tmp = pts[i]; pts[i] = pts[j]; pts[j] = tmp;
            }
    size_t m = 0;
    for (size_t i = 0; i < n; i++) {
        if (m > 0 && expr_eq(pts[m - 1], pts[i])) { expr_free(pts[i]); continue; }
        pts[m++] = pts[i];
    }

    Expr* out = expr_new_function(mk_sym("List"), pts, m);
    free(pts);
    return out;
}

void integrate_line_init(void) {
    symtab_add_builtin("Integrate`LineIntegral", builtin_integrate_line);
    symtab_get_def("Integrate`LineIntegral")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`LineIntegral",
        "Integrate`LineIntegral[f, {x, z0, z1, ..., zn}] evaluates the contour "
        "integral of f along the straight segments z0->z1->...->zn in the "
        "complex plane.  Each segment is parametrised by a real parameter, so "
        "singularities on the path are detected exactly and endpoint limits are "
        "taken along the contour; branch cuts of the antiderivative are handled "
        "so closed-contour residues are exact.  Returns unevaluated on a "
        "divergent contour or an unknown antiderivative.");

    symtab_add_builtin("Integrate`PathSingularPoints",
                       builtin_integrate_path_singular_points);
    symtab_get_def("Integrate`PathSingularPoints")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`PathSingularPoints",
        "Integrate`PathSingularPoints[f, {x, z0, ..., zn}] returns the sorted "
        "list of singular points of f lying strictly on the interior of the "
        "straight-line contour z0->z1->...->zn.");
}
