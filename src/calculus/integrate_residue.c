/* integrate_residue.c
 *
 * Definite integration by the residue theorem.  See integrate_residue.h for the
 * families handled and the surface forms that reach this module.
 *
 * Design.  Each family is a narrow, conjunctive recognizer that (i) checks the
 * bounds, (ii) reduces the integrand to a rational function whose relevant poles
 * it finds with Solve, (iii) classifies each pole against the family's contour
 * (upper half-plane / unit disk / real axis) by the numeric sign of Im or the
 * numeric value of Abs, (iv) sums the enclosed residues with the exported
 * residue_compute() primitive, and (v) closes the sum to a scalar (RootReduce
 * for the algebraic families A/C; Re/Im of the Jordan integral for the Fourier
 * family B).  The residue theorem makes each family's value correct BY
 * CONSTRUCTION once its structural gates hold, so there is NO numeric quadrature
 * (NIntegrate) crosscheck.  The only post-hoc gate is a self-consistency check
 * that the closed form actually reduced to a scalar (no surviving x / Root) and,
 * for the real-valued families, that its imaginary part vanishes (a residual Im
 * betrays a mis-classified pole).  Any failure returns NULL and the
 * Newton-Leibniz path takes over -- no existing integral is weakened.
 *
 * Real-axis poles.  A pure rational integrand (Family A) with a real pole is
 * genuinely singular there, so the ordinary improper integral does not converge
 * and the family returns NULL (only a principal value would exist, which plain
 * Integrate does not compute).  In the Fourier family (B) a real pole of the
 * rational part is admitted only when it is a REMOVABLE singularity of f -- a
 * simple pole whose kernel supplies a matching zero (e.g. Sin[x]/x at 0) -- and
 * then contributes a half residue (Pi I Res), the indented-contour value, giving
 * Integrate[Sin[x]/x] = Pi.  A genuine axis pole returns NULL.
 */

#include "integrate_residue.h"
#include "residue.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "sym_names.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* M_PI is POSIX, not C99; glibc hides it under -std=c99 (macOS exposes it). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Fresh variable names, context-qualified to stay clear of user symbols. */
#define RES_W  "IntegrateResidue`$w"   /* Family C: z = Exp[I x] on |z|=1        */

/* Classification tolerance for Im/Abs sign decisions (relative to magnitude). */
static const double RES_TOL = 1e-8;

/* -------------------------------------------------------------------------
 * Small construction / evaluation helpers.
 * ---------------------------------------------------------------------- */

static Expr* mk_sym(const char* s)  { return expr_new_symbol(s); }
static Expr* mk_int(int64_t n)      { return expr_new_integer(n); }

static Expr* mk_fn1(const char* h, Expr* a) {
    return expr_new_function(mk_sym(h), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(h), (Expr*[]){ a, b }, 2);
}

/* Evaluate `call`, free it, return the (owned) result. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* head[a] evaluated (a consumed). */
static Expr* ev1(const char* h, Expr* a) { return eval_take(mk_fn1(h, a)); }

/* True iff `e` is the compound `name[...]` by head name. */
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

/* True iff the symbol `x` occurs anywhere in `e`. */
static bool contains_symbol(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol == x->data.symbol;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol(e->data.function.args[i], x)) return true;
    return false;
}

/* -------------------------------------------------------------------------
 * Numeric extraction.
 * ---------------------------------------------------------------------- */

/* Machine double from a concrete real leaf. */
static bool res_real_double(const Expr* e, double* out) {
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;     return true;
        case EXPR_REAL:    *out = e->data.real;                return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint);   return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }
    int64_t n, d;
    if (is_rational((Expr*)e, &n, &d) && d != 0) { *out = (double)n / (double)d; return true; }
    return false;
}

/* Numeric real & imaginary parts of e (via N, then Re/Im).  Returns false unless
 * both parts are finite reals. */
static bool res_reim(Expr* e, double* re, double* im) {
    Expr* n = ev1("N", expr_copy(e));
    if (!n) return false;
    if (is_infinity_sym(n) || is_neg_infinity_form(n) ||
        is_complex_infinity_sym(n) || is_indeterminate_sym(n)) { expr_free(n); return false; }
    Expr* rr = ev1("Re", expr_copy(n));
    Expr* ii = ev1("Im", n);   /* n consumed */
    bool ok = rr && ii && res_real_double(rr, re) && res_real_double(ii, im);
    if (rr) expr_free(rr);
    if (ii) expr_free(ii);
    return ok;
}

/* Bound classification: 0 = symbolic/complex, 1 = finite (value in *v), 2 = +Inf,
 * 3 = -Inf. */
static int res_bound(Expr* e, double* v) {
    if (is_infinity_sym(e))      return 2;
    if (is_neg_infinity_form(e)) return 3;
    double re, im;
    if (!res_reim(e, &re, &im))  return 0;
    if (fabs(im) > RES_TOL * (1.0 + fabs(re))) return 0;   /* genuinely complex */
    *v = re;
    return 1;
}

/* True iff (a, b) is exactly (-Infinity, +Infinity). */
static bool is_neg_pos_infinity(Expr* a, Expr* b) {
    double v;
    return res_bound(a, &v) == 3 && res_bound(b, &v) == 2;
}

/* True iff (a, b) is (0, +Infinity). */
static bool is_zero_pos_infinity(Expr* a, Expr* b) {
    double av, bv;
    return res_bound(a, &av) == 1 && fabs(av) < RES_TOL && res_bound(b, &bv) == 2;
}

/* True iff (a, b) is a full trig period: (0, 2Pi) or (-Pi, Pi). */
static bool is_full_period(Expr* a, Expr* b) {
    double av, bv;
    if (res_bound(a, &av) != 1 || res_bound(b, &bv) != 1) return false;
    double two_pi = 2.0 * M_PI;
    if (fabs(av) < 1e-9 && fabs(bv - two_pi) < 1e-9) return true;          /* (0, 2Pi)  */
    if (fabs(av + M_PI) < 1e-9 && fabs(bv - M_PI) < 1e-9) return true;     /* (-Pi, Pi) */
    return false;
}

/* -------------------------------------------------------------------------
 * Rational structure: numerator / denominator polynomials and their degrees.
 * ---------------------------------------------------------------------- */

/* Degree of polynomial p in x via Length[CoefficientList[p, x]] - 1, or -1 if p
 * is not a polynomial in x. */
static int res_degree(Expr* p, Expr* x) {
    Expr* cl = eval_take(mk_fn2("CoefficientList", expr_copy(p), expr_copy(x)));
    if (!cl) return -1;
    int deg = -1;
    if (head_name_is(cl, "List")) deg = (int)cl->data.function.arg_count - 1;
    expr_free(cl);
    return deg;
}

/* PolynomialQ[p, x]. */
static bool res_polyq(Expr* p, Expr* x) {
    Expr* v = eval_take(mk_fn2("PolynomialQ", expr_copy(p), expr_copy(x)));
    bool ok = v && v->type == EXPR_SYMBOL && v->data.symbol == SYM_True;
    if (v) expr_free(v);
    return ok;
}

/* Together[expr] -> (Numerator, Denominator).  On success both *num, *den are
 * owned; returns false (nothing owned) if either cannot be formed. */
static bool res_num_den(Expr* expr, Expr** num, Expr** den) {
    *num = *den = NULL;
    Expr* tog = ev1("Together", expr_copy(expr));
    if (!tog) return false;
    Expr* n = ev1("Numerator", expr_copy(tog));
    Expr* d = ev1("Denominator", tog);   /* tog consumed */
    if (!n || !d) { if (n) expr_free(n); if (d) expr_free(d); return false; }
    *num = n; *den = d;
    return true;
}

/* -------------------------------------------------------------------------
 * Pole finding + classification.
 * ---------------------------------------------------------------------- */

typedef struct { Expr** v; size_t n, cap; } ExprVec;
static void ev_init(ExprVec* s) { s->v = NULL; s->n = s->cap = 0; }
static void ev_push(ExprVec* s, Expr* e) {
    if (s->n == s->cap) { s->cap = s->cap ? s->cap * 2 : 8;
                          s->v = realloc(s->v, s->cap * sizeof(*s->v)); }
    s->v[s->n++] = e;
}
static void ev_free(ExprVec* s) {
    for (size_t i = 0; i < s->n; i++) expr_free(s->v[i]);
    free(s->v); ev_init(s);
}

/* The RHS r of a Solve solution element `{x -> r}` (or bare `x -> r`); borrowed,
 * NULL if the shape does not match. */
static Expr* solve_rhs(Expr* el, const Expr* x) {
    if (head_name_is(el, "Rule") && el->data.function.arg_count == 2) {
        Expr* lhs = el->data.function.args[0];
        if (lhs->type == EXPR_SYMBOL && lhs->data.symbol == x->data.symbol)
            return el->data.function.args[1];
        return NULL;
    }
    if (head_name_is(el, "List"))
        for (size_t i = 0; i < el->data.function.arg_count; i++) {
            Expr* r = solve_rhs(el->data.function.args[i], x);
            if (r) return r;
        }
    return NULL;
}

/* Solve[Q == 0, x] and collect the DISTINCT root expressions into `roots`
 * (owned).  Solve reports a multiple root with multiplicity (e.g. (1+x^2)^2
 * yields {-I,-I,I,I}); we keep one representative per numeric location, since
 * residue_compute() already accounts for the full pole order at that point --
 * summing the same pole twice would double-count its residue.  Returns false if
 * Solve did not produce a plain List of solutions (symbolic / conditional). */
static bool solve_roots(Expr* Q, Expr* x, ExprVec* roots) {
    ev_init(roots);
    Expr* eq   = mk_fn2("Equal", expr_copy(Q), mk_int(0));
    Expr* sols = eval_take(mk_fn2("Solve", eq, expr_copy(x)));
    if (!sols) return false;
    if (!head_name_is(sols, "List")) { expr_free(sols); return false; }
    double* kre = NULL; double* kim = NULL; size_t nk = 0, cap = 0;
    for (size_t i = 0; i < sols->data.function.arg_count; i++) {
        Expr* r = solve_rhs(sols->data.function.args[i], x);
        if (!r) { expr_free(sols); ev_free(roots); free(kre); free(kim); return false; }
        double re, im;
        bool numeric = res_reim(r, &re, &im);
        if (numeric) {
            bool dup = false;
            double mag = 1.0 + fabs(re) + fabs(im);
            for (size_t j = 0; j < nk; j++)
                if (fabs(kre[j] - re) <= RES_TOL * mag && fabs(kim[j] - im) <= RES_TOL * mag)
                    { dup = true; break; }
            if (dup) continue;
            if (nk == cap) { cap = cap ? cap * 2 : 8;
                             kre = realloc(kre, cap * sizeof(double));
                             kim = realloc(kim, cap * sizeof(double)); }
            kre[nk] = re; kim[nk] = im; nk++;
        }
        ev_push(roots, expr_copy(r));
    }
    expr_free(sols);
    free(kre); free(kim);
    return true;
}

/* Is the real root r a SIMPLE root of the polynomial Q (i.e. Q'(r) != 0)?  Used
 * to gate the principal-value half residue: a higher-order pole on the axis is
 * not integrable this way. */
static bool res_simple_root(Expr* Q, Expr* x, Expr* r) {
    Expr* dQ = eval_take(mk_fn2("D", expr_copy(Q), expr_copy(x)));
    if (!dQ) return false;
    Expr* at = eval_take(mk_fn2("ReplaceAll", dQ, mk_fn2("Rule", expr_copy(x), expr_copy(r))));
    if (!at) return false;
    double re, im;
    bool ok = res_reim(at, &re, &im);
    expr_free(at);
    return ok && (fabs(re) + fabs(im)) > 1e-7;
}

/* -------------------------------------------------------------------------
 * Residue summation + closure.
 * ---------------------------------------------------------------------- */

/* S = Sum_i weight_i * Res[T, {v, poles[i]}], evaluated.  weight 2 for a fully
 * enclosed pole (2 Pi I), 1 for a half residue on the axis (Pi I).  Returns the
 * evaluated sum (owned) or NULL if any residue could not be computed. */
static Expr* sum_residues(Expr* T, Expr* v, Expr** poles, int* weight, size_t n) {
    Expr* total = mk_int(0);
    for (size_t i = 0; i < n; i++) {
        Expr* r = residue_compute(T, v, poles[i]);
        if (!r) { expr_free(total); return NULL; }
        Expr* term = (weight[i] == 2) ? mk_fn2("Times", mk_int(2), r) : r;
        total = eval_take(mk_fn2("Plus", total, term));
        if (!total) return NULL;
    }
    return total;
}

/* Algebraic closure (families A, C): value = Pi * RootReduce[I * S]. */
static Expr* close_algebraic(Expr* S) {
    Expr* alg = ev1("RootReduce", mk_fn2("Times", mk_sym(SYM_I), expr_copy(S)));
    if (!alg) return NULL;
    return ev1("Simplify", mk_fn2("Times", mk_sym(SYM_Pi), alg));
}

/* Consistency gates on the closed form.  The residue theorem makes each family's
 * value correct BY CONSTRUCTION once its structural gates hold, so no numeric
 * quadrature (NIntegrate) crosscheck is used.  These two checks only reject a
 * result that did not actually close -- a surviving `x`/`Root`, or (for the real
 * families) a non-vanishing imaginary part that betrays a mis-classified pole.
 *
 * True iff `value` is a finite numeric scalar free of x and Root/RootSum;
 * (out) its numeric real and imaginary parts. */
static bool res_is_finite_scalar(Expr* value, Expr* x, double* re, double* im) {
    if (!value) return false;
    if (contains_symbol(value, x)) return false;
    if (contains_head(value, "Root") || contains_head(value, "RootSum")) return false;
    return res_reim(value, re, im);
}

/* True iff `value` is additionally REAL (imaginary part negligible) -- the
 * correctness invariant for the real-valued families A / C and the Cos/Sin
 * Fourier answers.  A surviving imaginary part means a pole was placed in the
 * wrong half-plane, so the value is rejected. */
static bool res_is_real_scalar(Expr* value, Expr* x, double* out) {
    double re, im;
    if (!res_is_finite_scalar(value, x, &re, &im)) return false;
    if (fabs(im) > 1e-6 * (1.0 + fabs(re))) return false;
    *out = re;
    return true;
}

/* True iff the original kernel expression (Cos[a x] / Sin[a x] / Exp[I a x])
 * VANISHES at the real point r, decided symbolically by PossibleZeroQ.  A simple
 * real-axis pole of the rational part R contributes an integrable (removable)
 * singularity of f = R.K only when the kernel supplies a matching zero there
 * (e.g. Sin[x]/x at 0); otherwise f genuinely diverges on the axis and the
 * ordinary integral does not converge (only a principal value would). */
static bool res_kernel_vanishes_at(Expr* kernel, Expr* x, Expr* r) {
    Expr* at = eval_take(mk_fn2("ReplaceAll", expr_copy(kernel),
                                mk_fn2("Rule", expr_copy(x), expr_copy(r))));
    if (!at) return false;
    Expr* pz = ev1("PossibleZeroQ", at);   /* at consumed */
    bool ok = pz && pz->type == EXPR_SYMBOL && pz->data.symbol == SYM_True;
    if (pz) expr_free(pz);
    return ok;
}

/* -------------------------------------------------------------------------
 * Family A -- rational integrands on (-Inf, Inf).
 * ---------------------------------------------------------------------- */

static Expr* residue_family_rational(Expr* f, Expr* x, Expr* a, Expr* b) {
    if (!is_neg_pos_infinity(a, b)) return NULL;

    Expr* P; Expr* Q;
    if (!res_num_den(f, &P, &Q)) return NULL;
    if (!contains_symbol(Q, x) || !res_polyq(P, x) || !res_polyq(Q, x)) {
        expr_free(P); expr_free(Q); return NULL;
    }
    int dP = res_degree(P, x), dQ = res_degree(Q, x);
    expr_free(P);
    if (dP < 0 || dQ < dP + 2) { expr_free(Q); return NULL; }   /* convergence */

    ExprVec roots; ExprVec keep; ev_init(&keep);
    if (!solve_roots(Q, x, &roots)) { expr_free(Q); return NULL; }

    /* Sum the upper-half-plane residues.  A pole ON the real axis is a genuine
     * pole of the rational integrand (Together has already cancelled any common
     * factor), so the ordinary improper integral does not converge -- leave it
     * for Newton-Leibniz rather than silently returning a principal value. */
    bool bad = false;
    for (size_t i = 0; i < roots.n && !bad; i++) {
        double re, im;
        if (!res_reim(roots.v[i], &re, &im)) { bad = true; break; }   /* symbolic root */
        double mag = 1.0 + fabs(re) + fabs(im);
        if (im > RES_TOL * mag)          ev_push(&keep, expr_copy(roots.v[i])); /* UHP */
        else if (im < -RES_TOL * mag)    { /* lower half-plane: ignore */ }
        else                             { bad = true; break; }        /* real pole */
    }
    ev_free(&roots);
    if (bad || keep.n == 0) { expr_free(Q); ev_free(&keep); return NULL; }
    expr_free(Q);

    Expr** poles = malloc(keep.n * sizeof(*poles));
    int* weight  = malloc(keep.n * sizeof(*weight));
    for (size_t i = 0; i < keep.n; i++) { poles[i] = keep.v[i]; weight[i] = 2; }

    Expr* S = sum_residues(f, x, poles, weight, keep.n);
    free(poles); free(weight);
    if (!S) { ev_free(&keep); return NULL; }

    Expr* value = close_algebraic(S);
    expr_free(S);
    ev_free(&keep);

    double vv;
    if (!value || !res_is_real_scalar(value, x, &vv)) { if (value) expr_free(value); return NULL; }
    return value;
}

/* -------------------------------------------------------------------------
 * Family B -- Fourier / Jordan integrands on (-Inf, Inf).
 * ---------------------------------------------------------------------- */

typedef enum { KERN_NONE, KERN_COS, KERN_SIN, KERN_EXP } KernelKind;

/* If `e` is a trig/exp kernel Cos[a x] / Sin[a x] / Exp[I a x] with a a nonzero
 * real constant, set *kind and *a (the real frequency) and return true. */
static bool match_kernel(Expr* e, Expr* x, KernelKind* kind, double* a) {
    if (!head_name_is(e, "Cos") && !head_name_is(e, "Sin") && !head_name_is(e, "Exp"))
        return false;
    if (e->data.function.arg_count != 1) return false;
    Expr* u = e->data.function.args[0];
    if (!contains_symbol(u, x)) return false;

    /* u must be linear in x with zero constant term. */
    Expr* cl = eval_take(mk_fn2("CoefficientList", expr_copy(u), expr_copy(x)));
    if (!cl || !head_name_is(cl, "List") || cl->data.function.arg_count != 2) {
        if (cl) expr_free(cl); return false;
    }
    Expr* c0 = cl->data.function.args[0];   /* constant term (must be 0)   */
    Expr* c1 = expr_copy(cl->data.function.args[1]);   /* slope             */
    double c0r, c0i, c1r, c1i;
    bool okc0 = res_reim(c0, &c0r, &c0i);
    bool okc1 = res_reim(c1, &c1r, &c1i);
    expr_free(cl); expr_free(c1);
    if (!okc0 || !okc1) return false;
    if (fabs(c0r) > 1e-9 || fabs(c0i) > 1e-9) return false;   /* nonzero offset */

    if (head_name_is(e, "Exp")) {
        /* slope must be purely imaginary I a: a = Im(slope), Re ~ 0. */
        if (fabs(c1r) > 1e-9 || fabs(c1i) < 1e-9) return false;
        *kind = KERN_EXP; *a = c1i;
    } else {
        /* Cos/Sin: slope must be real a. */
        if (fabs(c1i) > 1e-9 || fabs(c1r) < 1e-9) return false;
        *kind = head_name_is(e, "Cos") ? KERN_COS : KERN_SIN;
        *a = c1r;
    }
    return true;
}

static Expr* residue_family_fourier(Expr* f, Expr* x, Expr* a, Expr* b) {
    if (!is_neg_pos_infinity(a, b)) return NULL;

    /* Split f into a product of factors and find exactly one trig/exp kernel. */
    Expr** factors; size_t nf;
    Expr* single[1];
    if (head_name_is(f, "Times")) {
        factors = f->data.function.args; nf = f->data.function.arg_count;
    } else { single[0] = f; factors = single; nf = 1; }

    KernelKind kind = KERN_NONE; double freq = 0.0; int kidx = -1;
    for (size_t i = 0; i < nf; i++) {
        KernelKind kk; double kf;
        if (match_kernel(factors[i], x, &kk, &kf)) {
            if (kidx >= 0) return NULL;   /* more than one kernel: not this family */
            kind = kk; freq = kf; kidx = (int)i;
        }
    }
    if (kidx < 0) return NULL;

    /* R = product of the remaining factors. */
    Expr* R;
    if (nf == 1) {
        R = mk_int(1);   /* f was a bare kernel: R = 1 (nonconvergent; handled below) */
    } else {
        Expr** ra = malloc((nf - 1) * sizeof(*ra));
        size_t k = 0;
        for (size_t i = 0; i < nf; i++) if ((int)i != kidx) ra[k++] = expr_copy(factors[i]);
        R = eval_take(expr_new_function(mk_sym("Times"), ra, nf - 1));
        free(ra);
    }
    if (!R) return NULL;

    /* R must be rational in x with deg-drop >= 1 (Jordan: R -> 0). */
    Expr* P; Expr* Q;
    if (!res_num_den(R, &P, &Q) || !contains_symbol(Q, x) ||
        !res_polyq(P, x) || !res_polyq(Q, x)) {
        expr_free(R); if (P) expr_free(P); if (Q) expr_free(Q); return NULL;
    }
    int dP = res_degree(P, x), dQ = res_degree(Q, x);
    expr_free(P);
    if (dP < 0 || dQ < dP + 1) { expr_free(R); expr_free(Q); return NULL; }

    /* g = R * Exp[expo] (the Jordan integrand), where expo is the EXACT exponent
     * built from the kernel's own argument u (u = a x): Cos/Sin -> I u, Exp -> u.
     * Using the exact u (not the numeric freq) keeps residues exact. */
    Expr* kernel = factors[kidx];   /* borrowed: the original Cos/Sin/Exp factor */
    Expr* u = kernel->data.function.args[0];
    Expr* expo = (kind == KERN_EXP) ? expr_copy(u)
                                    : mk_fn2("Times", mk_sym(SYM_I), expr_copy(u));
    Expr* g = eval_take(mk_fn2("Times", expr_copy(R), mk_fn1("Exp", expo)));
    expr_free(R);
    if (!g) { expr_free(Q); return NULL; }

    bool up = freq > 0.0;   /* close upper (a>0) or lower (a<0) half-plane */

    ExprVec roots; ExprVec keep; ExprVec realp;
    ev_init(&keep); ev_init(&realp);
    if (!solve_roots(Q, x, &roots)) { expr_free(Q); expr_free(g); return NULL; }

    bool bad = false;
    for (size_t i = 0; i < roots.n && !bad; i++) {
        double re, im;
        if (!res_reim(roots.v[i], &re, &im)) { bad = true; break; }
        double mag = 1.0 + fabs(re) + fabs(im);
        bool enclosed = up ? (im > RES_TOL * mag) : (im < -RES_TOL * mag);
        bool other    = up ? (im < -RES_TOL * mag) : (im > RES_TOL * mag);
        if (enclosed) ev_push(&keep, expr_copy(roots.v[i]));
        else if (other) { /* opposite half-plane: not enclosed */ }
        else {   /* real axis: valid only as a removable singularity of f */
            /* A real pole of R gives an integrable f = R.K only when it is simple
             * AND the kernel vanishes there (e.g. Sin[x]/x at 0); the half residue
             * then equals the indented-contour contribution.  Otherwise f has a
             * genuine axis pole and the ordinary integral diverges -> bail.  The
             * a<0 (lower-plane) closure is not handled for on-axis poles. */
            if (!up ||
                !res_simple_root(Q, x, roots.v[i]) ||
                !res_kernel_vanishes_at(kernel, x, roots.v[i])) { bad = true; break; }
            ev_push(&realp, expr_copy(roots.v[i]));
        }
    }
    ev_free(&roots);
    expr_free(Q);
    if (bad || keep.n + realp.n == 0) {
        expr_free(g); ev_free(&keep); ev_free(&realp); return NULL;
    }

    size_t nk = keep.n + realp.n;
    Expr** poles = malloc(nk * sizeof(*poles));
    int* weight = malloc(nk * sizeof(*weight));
    size_t kk = 0;
    for (size_t i = 0; i < keep.n; i++)  { poles[kk] = keep.v[i];  weight[kk] = 2; kk++; }
    for (size_t i = 0; i < realp.n; i++) { poles[kk] = realp.v[i]; weight[kk] = 1; kk++; }

    Expr* S = sum_residues(g, x, poles, weight, nk);
    free(poles); free(weight);
    expr_free(g);
    ev_free(&keep); ev_free(&realp);
    if (!S) return NULL;

    /* J = Pi I S, times the closure orientation (upper +, lower -). */
    Expr* sign = up ? mk_int(1) : mk_int(-1);
    Expr* J = ev1("Simplify",
                  ev1("Together",
                      mk_fn2("Times", mk_fn2("Times", mk_sym(SYM_Pi), mk_sym(SYM_I)),
                             mk_fn2("Times", sign, expr_copy(S)))));
    expr_free(S);
    if (!J) return NULL;

    /* Kernel-specific answer: Exp -> J, Cos -> Re[J], Sin -> Im[J], via the exact
     * Conjugate identities (symbolic Re/Im do not reduce Pi/E etc.). */
    Expr* value;
    if (kind == KERN_EXP) {
        value = expr_copy(J);
    } else if (kind == KERN_COS) {
        value = ev1("Simplify",
                    mk_fn2("Times", mk_fn2("Power", mk_int(2), mk_int(-1)),
                           mk_fn2("Plus", expr_copy(J), mk_fn1("Conjugate", expr_copy(J)))));
    } else { /* KERN_SIN */
        Expr* diff = mk_fn2("Plus", expr_copy(J),
                            mk_fn2("Times", mk_int(-1), mk_fn1("Conjugate", expr_copy(J))));
        value = ev1("Simplify", mk_fn2("Times", diff,
                        mk_fn2("Power", mk_fn2("Times", mk_int(2), mk_sym(SYM_I)), mk_int(-1))));
    }
    expr_free(J);

    /* Cos/Sin integrals are real by construction; the Exp integral may be
     * genuinely complex, so only require a finite scalar there. */
    double re, im, vv;
    bool ok = (kind == KERN_EXP) ? (value && res_is_finite_scalar(value, x, &re, &im))
                                 : (value && res_is_real_scalar(value, x, &vv));
    if (!ok) { if (value) expr_free(value); return NULL; }
    return value;
}

/* -------------------------------------------------------------------------
 * Family C -- rational-in-{Sin,Cos} integrands over a full period.
 * ---------------------------------------------------------------------- */

static Expr* residue_family_trig(Expr* f, Expr* x, Expr* a, Expr* b) {
    if (!is_full_period(a, b)) return NULL;

    Expr* z = mk_sym(RES_W);

    /* TrigExpand, then z = Exp[I x]: Cos[x] -> (z+1/z)/2, Sin[x] -> (z-1/z)/(2 I). */
    Expr* te = ev1("TrigExpand", expr_copy(f));
    if (!te) { expr_free(z); return NULL; }

    Expr* cos_rule = mk_fn2("Rule", mk_fn1("Cos", expr_copy(x)),
        mk_fn2("Times", mk_fn2("Plus", expr_copy(z), mk_fn2("Power", expr_copy(z), mk_int(-1))),
               mk_fn2("Power", mk_int(2), mk_int(-1))));
    Expr* sin_rule = mk_fn2("Rule", mk_fn1("Sin", expr_copy(x)),
        mk_fn2("Times", mk_fn2("Plus", expr_copy(z),
                    mk_fn2("Times", mk_int(-1), mk_fn2("Power", expr_copy(z), mk_int(-1)))),
               mk_fn2("Power", mk_fn2("Times", mk_int(2), mk_sym(SYM_I)), mk_int(-1))));
    Expr* rules = mk_fn2("List", cos_rule, sin_rule);
    Expr* subst = eval_take(mk_fn2("ReplaceAll", te, rules));   /* te consumed */
    if (!subst) { expr_free(z); return NULL; }

    /* Contour integrand F(z) = subst / (I z). */
    Expr* Fz = ev1("Together",
                   mk_fn2("Times", subst,
                          mk_fn2("Power", mk_fn2("Times", mk_sym(SYM_I), expr_copy(z)), mk_int(-1))));
    if (!Fz) { expr_free(z); return NULL; }

    /* Must be rational in z (leftover trig heads over x -> reject). */
    if (contains_symbol(Fz, x)) { expr_free(z); expr_free(Fz); return NULL; }
    Expr* P; Expr* Q;
    if (!res_num_den(Fz, &P, &Q) || !res_polyq(P, z) || !res_polyq(Q, z)) {
        expr_free(z); expr_free(Fz); if (P) expr_free(P); if (Q) expr_free(Q); return NULL;
    }
    expr_free(P);

    ExprVec roots; ExprVec keep; ev_init(&keep);
    if (!solve_roots(Q, z, &roots)) { expr_free(z); expr_free(Fz); expr_free(Q); return NULL; }

    bool bad = false;
    for (size_t i = 0; i < roots.n && !bad; i++) {
        double re, im;
        if (!res_reim(roots.v[i], &re, &im)) { bad = true; break; }
        double mag = sqrt(re * re + im * im);
        if (mag < 1.0 - RES_TOL)      ev_push(&keep, expr_copy(roots.v[i]));  /* inside disk */
        else if (mag > 1.0 + RES_TOL) { /* outside: skip */ }
        else { bad = true; break; }   /* on the unit circle: real-axis pole -> divergent */
    }
    ev_free(&roots);
    expr_free(Q);
    if (bad || keep.n == 0) { expr_free(z); expr_free(Fz); ev_free(&keep); return NULL; }

    Expr** poles = malloc(keep.n * sizeof(*poles));
    int* weight  = malloc(keep.n * sizeof(*weight));
    for (size_t i = 0; i < keep.n; i++) { poles[i] = keep.v[i]; weight[i] = 2; }

    Expr* S = sum_residues(Fz, z, poles, weight, keep.n);
    free(poles); free(weight);
    expr_free(Fz);
    if (!S) { expr_free(z); ev_free(&keep); return NULL; }

    Expr* value = close_algebraic(S);
    expr_free(S); expr_free(z);
    ev_free(&keep);

    double vv;
    if (!value || !res_is_real_scalar(value, x, &vv)) { if (value) expr_free(value); return NULL; }
    return value;
}

/* -------------------------------------------------------------------------
 * Half-line [0, Inf) via even symmetry.
 * ---------------------------------------------------------------------- */

/* True iff f(x) == f(-x) (even), decided by PossibleZeroQ. */
static bool is_even_in(Expr* f, Expr* x) {
    Expr* fm = eval_take(mk_fn2("ReplaceAll", expr_copy(f),
                                mk_fn2("Rule", expr_copy(x), mk_fn2("Times", mk_int(-1), expr_copy(x)))));
    if (!fm) return false;
    Expr* diff = mk_fn2("Plus", expr_copy(f), mk_fn2("Times", mk_int(-1), fm));
    Expr* pz = ev1("PossibleZeroQ", diff);
    bool ok = pz && pz->type == EXPR_SYMBOL && pz->data.symbol == SYM_True;
    if (pz) expr_free(pz);
    return ok;
}

static Expr* residue_half_line(Expr* f, Expr* x, Expr* a, Expr* b) {
    if (!is_zero_pos_infinity(a, b)) return NULL;
    if (!is_even_in(f, x)) return NULL;

    /* -Infinity as DirectedInfinity[-1] via -1 * Infinity. */
    Expr* neg = eval_take(mk_fn2("Times", mk_int(-1), mk_sym(SYM_Infinity)));
    Expr* pos = mk_sym(SYM_Infinity);
    Expr* full = integrate_residue_try(f, x, neg, pos);   /* recurse: full line */
    expr_free(neg); expr_free(pos);
    if (!full) return NULL;

    /* full is adopted into the Times node and freed by evaluation -- do not
     * free it again here.  It has already passed the full-line family's own
     * consistency gate, so halving an even integrand needs no further check
     * beyond confirming the halved form is still a real scalar. */
    Expr* half = ev1("Simplify", mk_fn2("Times", mk_fn2("Power", mk_int(2), mk_int(-1)), full));
    double vv;
    if (!half || !res_is_real_scalar(half, x, &vv)) { if (half) expr_free(half); return NULL; }
    return half;
}

/* -------------------------------------------------------------------------
 * Master entry + builtin.
 * ---------------------------------------------------------------------- */

Expr* integrate_residue_try(Expr* f, Expr* x, Expr* a, Expr* b) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;

    /* Full period (0,2Pi)/(-Pi,Pi): the unit-circle trig-substitution family. */
    if (is_full_period(a, b)) return residue_family_trig(f, x, a, b);

    /* Whole line: Fourier/Jordan if a trig/exp kernel is present, else rational. */
    if (is_neg_pos_infinity(a, b)) {
        Expr* r = residue_family_fourier(f, x, a, b);
        if (r) return r;
        return residue_family_rational(f, x, a, b);
    }

    /* Half-line [0,Inf) for even integrands. */
    if (is_zero_pos_infinity(a, b)) return residue_half_line(f, x, a, b);

    return NULL;
}

Expr* builtin_integrate_contour_residue(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!head_name_is(spec, "List") || spec->data.function.arg_count != 3) return NULL;
    Expr* x = spec->data.function.args[0];
    if (x->type != EXPR_SYMBOL) return NULL;
    Expr* a = spec->data.function.args[1];
    Expr* b = spec->data.function.args[2];
    return integrate_residue_try(f, x, a, b);
}

void integrate_residue_init(void) {
    symtab_add_builtin("Integrate`ContourResidue", builtin_integrate_contour_residue);
    symtab_get_def("Integrate`ContourResidue")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`ContourResidue",
        "Integrate`ContourResidue[f, {x, a, b}] evaluates a definite real "
        "integral by the residue theorem, for rational integrands on "
        "(-Infinity, Infinity), Fourier/Jordan integrands R(x) {Cos|Sin|Exp[I .]}"
        "(a x), rational-in-{Sin,Cos} integrands over a full period (0,2Pi) or "
        "(-Pi,Pi), and their principal-value (simple real-axis pole) and even "
        "half-line [0,Infinity) variants.  Sums enclosed residues and verifies "
        "the closed form against NIntegrate; returns unevaluated when no family "
        "applies, the residue sum does not close to a scalar, or the crosscheck "
        "fails.");
}
