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
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
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
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == x->data.symbol.name;
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

/* -------------------------------------------------------------------------
 * Symbolic-parameter mode: sign-consistent numeric instantiation.
 *
 * When Integrate is called with an `Assumptions -> ...` option that constrains
 * the integrand's free parameters (e.g. `a > 0`, `0 < a < 1`, `n > 1`), a family
 * cannot classify a parameter-dependent pole or kernel frequency by a direct
 * numeric evaluation -- N[I a] stays symbolic.  Instead we pick ONE generic
 * value for each parameter consistent with the assumptions (`g_inst`, a List of
 * Rule[param, value]) and classify at that point, while the residue arithmetic
 * stays fully symbolic.  A single generic representative determines the
 * half-plane / sign whenever the assumptions pin it on a connected region (a
 * pole that stayed off the contour cannot have crossed it), so the resulting
 * symbolic closed form is correct by construction -- there is no numeric
 * quadrature crosscheck.  Convergence/applicability gates that could vary across
 * the region are instead checked against the assumption-guaranteed bounds (see
 * param_interval), so an under-constrained problem is refused, not guessed.
 *
 * g_inst is NULL outside symbolic-parameter mode, so every helper below behaves
 * exactly as before when no assumptions are supplied.
 * ---------------------------------------------------------------------- */
static Expr* g_inst = NULL;   /* List[Rule[param, Real], ...] or NULL */
static bool  g_all_pos = false;  /* true iff every instantiated parameter is > 0 */
/* The assumption-guaranteed interval of each parameter (mirrors g_inst), so a
 * family can verify a convergence/applicability gate holds over the WHOLE
 * assumed region -- not merely at the single instantiation point (which would
 * wrongly accept e.g. n > 0 for a gate that needs n > 1). */
typedef struct { const char* sym; double lo, hi; } ParamBound;
static ParamBound g_bounds[16];
static size_t     g_nbounds = 0;

/* e with the parameter instantiation applied and evaluated (owned), or NULL if
 * no instantiation is active. */
static Expr* apply_inst(Expr* e) {
    if (!g_inst) return NULL;
    return eval_take(mk_fn2("ReplaceAll", expr_copy(e), expr_copy(g_inst)));
}

/* Numeric real & imaginary parts of a CONCRETE numeric e (via N, then Re/Im). */
static bool res_reim_direct(Expr* e, double* re, double* im) {
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

/* Numeric real & imaginary parts of e.  Tries a direct evaluation first; when
 * that leaves a symbolic residue and a parameter instantiation is active, it
 * retries under the instantiation, so a parameter-dependent quantity (a pole
 * `I a`, a kernel slope `k`) classifies at the generic point.  Returns false
 * unless both parts come out as finite reals. */
static bool res_reim(Expr* e, double* re, double* im) {
    if (res_reim_direct(e, re, im)) return true;
    if (g_inst) {
        Expr* s = apply_inst(e);
        if (s) { bool ok = res_reim_direct(s, re, im); expr_free(s); return ok; }
    }
    return false;
}

/* PowerExpand e when all parameters are known positive -- exposes the imaginary
 * unit inside radical pole locations (Sqrt[-4 a^2] -> 2 I a) so downstream
 * residue arithmetic and conjugation stay explicit.  Consumes e, returns owned;
 * a no-op outside all-positive symbolic-parameter mode. */
static Expr* res_powerclean(Expr* e) {
    if (g_inst && g_all_pos) return ev1("PowerExpand", e);
    return e;
}

/* Complex conjugate of J.  In symbolic-parameter mode the parameters are real
 * and the imaginaries are the explicit unit I (after res_powerclean), so
 * conjugation is J with I -> -I; the symbolic Conjugate head would not reduce.
 * Outside that mode, the ordinary Conjugate head (numeric params resolve it). */
static Expr* res_conjugate(Expr* J) {
    if (g_inst)
        return eval_take(mk_fn2("ReplaceAll", expr_copy(J),
                                mk_fn2("Rule", mk_sym(SYM_I),
                                       mk_fn2("Times", mk_int(-1), mk_sym(SYM_I)))));
    return mk_fn1("Conjugate", expr_copy(J));
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
    bool ok = v && v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_True;
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
        if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == x->data.symbol.name)
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
    bool ok = pz && pz->type == EXPR_SYMBOL && pz->data.symbol.name == SYM_True;
    if (pz) expr_free(pz);
    return ok;
}

/* -------------------------------------------------------------------------
 * Family A -- rational integrands on (-Inf, Inf).
 * ---------------------------------------------------------------------- */

static Expr* residue_family_rational(Expr* f, Expr* x, Expr* a, Expr* b,
                                     bool* diverges) {
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
     * factor) and, since the degree gate above guarantees convergence at
     * infinity, the ONLY obstruction: the ordinary improper integral does not
     * converge.  Flag that conclusively (the dispatcher emits Integrate::idiv)
     * -- distinct from an undecidable symbolic root, which stays a plain bail. */
    bool real_pole = false, undecidable = false;
    for (size_t i = 0; i < roots.n; i++) {
        double re, im;
        if (!res_reim(roots.v[i], &re, &im)) { undecidable = true; break; } /* symbolic root */
        double mag = 1.0 + fabs(re) + fabs(im);
        if (im > RES_TOL * mag)          ev_push(&keep, expr_copy(roots.v[i])); /* UHP */
        else if (im < -RES_TOL * mag)    { /* lower half-plane: ignore */ }
        else                             { real_pole = true; break; }   /* real-axis pole */
    }
    ev_free(&roots);
    if (real_pole) {
        if (diverges) *diverges = true;
        expr_free(Q); ev_free(&keep); return NULL;
    }
    if (undecidable || keep.n == 0) { expr_free(Q); ev_free(&keep); return NULL; }
    expr_free(Q);

    /* Owned pole copies, PowerExpand-cleaned under all-positive symbolic-
     * parameter mode so a radical location becomes an explicit I a and the
     * residue closes to a clean rational form (Pi/a rather than a Sqrt[-4a^2]
     * surface). */
    Expr** poles = malloc(keep.n * sizeof(*poles));
    int* weight  = malloc(keep.n * sizeof(*weight));
    for (size_t i = 0; i < keep.n; i++) { poles[i] = res_powerclean(expr_copy(keep.v[i])); weight[i] = 2; }

    Expr* S = sum_residues(f, x, poles, weight, keep.n);
    for (size_t i = 0; i < keep.n; i++) expr_free(poles[i]);
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

/* The "phase argument" u of a Fourier kernel: Cos[u] / Sin[u] / Exp[u] carry it
 * as args[0], while the evaluator normalises Exp[I a x] to the Power form
 * Power[E, u] (u the exponent, args[1]).  Returns a borrowed pointer to u, or
 * NULL if `e` is not one of these shapes.  Power is accepted ONLY with base E,
 * so ordinary powers (x^2, denominators (1+x^2)^-1) are rejected. */
static Expr* kernel_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (head_name_is(e, "Cos") || head_name_is(e, "Sin") || head_name_is(e, "Exp")) {
        if (e->data.function.arg_count != 1) return NULL;
        return e->data.function.args[0];
    }
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == SYM_E)
            return e->data.function.args[1];
    }
    return NULL;
}

/* If `e` is a trig/exp kernel Cos[a x] / Sin[a x] / Exp[I a x] (in either the
 * Exp[.] or the normalised Power[E, .] spelling) with a a nonzero real constant,
 * set *kind and *a (the real frequency) and return true. */
static bool match_kernel(Expr* e, Expr* x, KernelKind* kind, double* a) {
    Expr* u = kernel_arg(e);
    if (!u) return false;
    /* Exp-form kernels: literal Exp[.] or the normalised Power[E, .]. */
    bool is_exp = head_name_is(e, "Exp") || head_name_is(e, "Power");
    if (!contains_symbol(u, x)) return false;

    /* u must be linear in x with zero constant term. */
    Expr* cl = eval_take(mk_fn2("CoefficientList", expr_copy(u), expr_copy(x)));
    if (!cl || !head_name_is(cl, "List") || cl->data.function.arg_count != 2) {
        if (cl) { expr_free(cl); } return false;
    }
    Expr* c0 = cl->data.function.args[0];   /* constant term (must be 0)   */
    Expr* c1 = expr_copy(cl->data.function.args[1]);   /* slope             */
    double c0r, c0i, c1r, c1i;
    bool okc0 = res_reim(c0, &c0r, &c0i);
    bool okc1 = res_reim(c1, &c1r, &c1i);
    expr_free(cl); expr_free(c1);
    if (!okc0 || !okc1) return false;
    if (fabs(c0r) > 1e-9 || fabs(c0i) > 1e-9) return false;   /* nonzero offset */

    if (is_exp) {
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
    Expr* kernel = factors[kidx];   /* borrowed: original Cos/Sin/Exp/Power[E,.] factor */
    Expr* u = kernel_arg(kernel);   /* args[0] for Cos/Sin/Exp; the exponent for Power[E,.] */
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

    /* Owned pole copies, PowerExpand-cleaned under symbolic-parameter mode so a
     * radical location (Sqrt[-4 a^2]/2) becomes the explicit I a: the residue at
     * an explicit purely-imaginary pole reduces to a clean real form, and the
     * Conjugate closure below can conjugate it via I -> -I. */
    size_t nk = keep.n + realp.n;
    Expr** poles = malloc(nk * sizeof(*poles));
    int* weight = malloc(nk * sizeof(*weight));
    size_t kk = 0;
    for (size_t i = 0; i < keep.n; i++)  { poles[kk] = res_powerclean(expr_copy(keep.v[i]));  weight[kk] = 2; kk++; }
    for (size_t i = 0; i < realp.n; i++) { poles[kk] = res_powerclean(expr_copy(realp.v[i])); weight[kk] = 1; kk++; }

    Expr* S = sum_residues(g, x, poles, weight, nk);
    for (size_t i = 0; i < nk; i++) expr_free(poles[i]);
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
                           mk_fn2("Plus", expr_copy(J), res_conjugate(J))));
    } else { /* KERN_SIN */
        Expr* diff = mk_fn2("Plus", expr_copy(J),
                            mk_fn2("Times", mk_int(-1), res_conjugate(J)));
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

static Expr* residue_family_trig(Expr* f, Expr* x, Expr* a, Expr* b,
                                 bool* diverges) {
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

    /* A pole ON the unit circle |z|=1 is a pole on the real x-axis where the
     * trig denominator vanishes: the periodic integral does not converge.  Flag
     * it conclusively (dispatcher emits Integrate::idiv); a symbolic root is a
     * plain undecidable bail. */
    bool circle_pole = false, undecidable = false;
    for (size_t i = 0; i < roots.n; i++) {
        double re, im;
        if (!res_reim(roots.v[i], &re, &im)) { undecidable = true; break; }
        double mag = sqrt(re * re + im * im);
        if (mag < 1.0 - RES_TOL)      ev_push(&keep, expr_copy(roots.v[i]));  /* inside disk */
        else if (mag > 1.0 + RES_TOL) { /* outside: skip */ }
        else { circle_pole = true; break; }   /* on the unit circle -> divergent */
    }
    ev_free(&roots);
    expr_free(Q);
    if (circle_pole) {
        if (diverges) *diverges = true;
        expr_free(z); expr_free(Fz); ev_free(&keep); return NULL;
    }
    if (undecidable || keep.n == 0) { expr_free(z); expr_free(Fz); ev_free(&keep); return NULL; }

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
    bool ok = pz && pz->type == EXPR_SYMBOL && pz->data.symbol.name == SYM_True;
    if (pz) expr_free(pz);
    return ok;
}

static Expr* residue_half_line(Expr* f, Expr* x, Expr* a, Expr* b,
                               bool* diverges) {
    if (!is_zero_pos_infinity(a, b)) return NULL;
    if (!is_even_in(f, x)) return NULL;

    /* -Infinity as DirectedInfinity[-1] via -1 * Infinity. */
    Expr* neg = eval_take(mk_fn2("Times", mk_int(-1), mk_sym(SYM_Infinity)));
    Expr* pos = mk_sym(SYM_Infinity);
    /* Recurse over the full line: a divergent full-line even integrand means the
     * half-line [0, Inf) diverges too, so propagate the flag. */
    Expr* full = integrate_residue_try(f, x, neg, pos, NULL, diverges);
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
 * Symbolic-parameter instantiation (built from Assumptions).
 * ---------------------------------------------------------------------- */

/* Real machine value of a concrete numeric leaf/expression (via N), real part
 * only.  Returns false for anything that does not reduce to a finite real -- in
 * particular a genuinely symbolic parameter (which is what we want to detect). */
static bool numeric_double(Expr* e, double* out) {
    Expr* n = ev1("N", expr_copy(e));
    if (!n) return false;
    bool ok = res_real_double(n, out);
    expr_free(n);
    return ok;
}

/* Collect the distinct free-parameter symbols of `e` (symbols in argument
 * position that are neither the integration variable `x`, the imaginary unit,
 * nor a numeric constant like Pi/E) into `pb` (up to `cap`).  Function heads are
 * skipped.  Returns the count. */
static size_t collect_params(const Expr* e, const Expr* x, ParamBound* pb,
                             size_t cap, size_t n) {
    if (!e) return n;
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol.name == x->data.symbol.name) return n;
        if (e->data.symbol.name == SYM_I) return n;
        for (size_t i = 0; i < n; i++)
            if (pb[i].sym == e->data.symbol.name) return n;   /* already seen */
        double tmp;
        if (numeric_double((Expr*)e, &tmp)) return n;    /* a numeric constant */
        if (n < cap) { pb[n].sym = e->data.symbol.name; pb[n].lo = -HUGE_VAL; pb[n].hi = HUGE_VAL; n++; }
        return n;
    }
    if (e->type != EXPR_FUNCTION) return n;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        n = collect_params(e->data.function.args[i], x, pb, cap, n);
    return n;
}

/* Tighten the bound of the parameter named by `var` given `var op const`.
 * `dir` = +1 for a lower bound (var > c / var >= c), -1 for an upper bound. */
static void bound_apply(ParamBound* pb, size_t np, const char* var, double c, int dir) {
    for (size_t i = 0; i < np; i++) {
        if (pb[i].sym != var) continue;
        if (dir > 0) { if (c > pb[i].lo) pb[i].lo = c; }
        else         { if (c < pb[i].hi) pb[i].hi = c; }
        return;
    }
}

/* Absorb a single ordered relation `L op R` (op one of Less/LessEqual/Greater/
 * GreaterEqual) into the parameter bounds: whichever side is a bare parameter
 * and the other a numeric constant tightens that parameter's interval. */
static void relation_apply(ParamBound* pb, size_t np, const char* op,
                           Expr* L, Expr* R) {
    bool lt = (strcmp(op, "Less") == 0 || strcmp(op, "LessEqual") == 0);
    bool gt = (strcmp(op, "Greater") == 0 || strcmp(op, "GreaterEqual") == 0);
    if (!lt && !gt) return;
    double c;
    if (L->type == EXPR_SYMBOL && numeric_double(R, &c))
        bound_apply(pb, np, L->data.symbol.name, c, lt ? -1 : +1);   /* L<c: upper; L>c: lower */
    else if (R->type == EXPR_SYMBOL && numeric_double(L, &c))
        bound_apply(pb, np, R->data.symbol.name, c, lt ? +1 : -1);   /* c<R: lower; c>R: upper */
}

/* Walk an assumption fact expression, tightening parameter bounds.  Handles
 * And / List conjunctions, chained Inequality[e,op,e,op,e,...], and the four
 * ordered binary relations. */
static void absorb_fact(ParamBound* pb, size_t np, Expr* fact) {
    if (!fact || fact->type != EXPR_FUNCTION) return;
    if (fact->data.function.head->type != EXPR_SYMBOL) return;
    const char* h = fact->data.function.head->data.symbol.name;
    size_t ac = fact->data.function.arg_count;
    if (strcmp(h, "And") == 0 || strcmp(h, "List") == 0) {
        for (size_t i = 0; i < ac; i++) absorb_fact(pb, np, fact->data.function.args[i]);
        return;
    }
    if (strcmp(h, "Inequality") == 0 && ac >= 3) {
        /* e0 op1 e1 op2 e2 ...: each (e_{2i}, op_{2i+1}, e_{2i+2}) is a relation. */
        for (size_t i = 0; i + 2 < ac; i += 2) {
            Expr* opsym = fact->data.function.args[i + 1];
            if (opsym->type != EXPR_SYMBOL) continue;
            relation_apply(pb, np, opsym->data.symbol.name,
                           fact->data.function.args[i], fact->data.function.args[i + 2]);
        }
        return;
    }
    if (ac == 2) relation_apply(pb, np, h,
                                fact->data.function.args[0], fact->data.function.args[1]);
}

/* A generic representative strictly inside [lo, hi], distinct per index to avoid
 * accidental coincidences (a value where two poles collide or a residue happens
 * to vanish).  Deterministic (no RNG): reproducible instantiations. */
static double pick_representative(double lo, double hi, size_t idx) {
    static const double FRACS[] = { 0.4142136, 0.6180340, 0.3183099, 0.5857864,
                                    0.2679492, 0.7320508, 0.5411961 };
    static const double SEEDS[] = { 1.3172, 2.7391, 0.6180, 3.1490,
                                    1.4213, 2.2360, 0.8177 };
    const size_t NF = sizeof(FRACS) / sizeof(FRACS[0]);
    double frac = FRACS[idx % NF], seed = SEEDS[idx % NF];
    bool lo_fin = lo > -HUGE_VAL, hi_fin = hi < HUGE_VAL;
    if (lo_fin && hi_fin) return lo + (hi - lo) * frac;
    if (lo_fin)           return lo + seed;
    if (hi_fin)           return hi - seed;
    return seed;                 /* unbounded: a generic positive value */
}

/* Build the parameter instantiation `List[Rule[p, val], ...]` from the integrand
 * and the assumptions.  The value chosen for each parameter is a generic point
 * of the interval the assumptions pin it to; it is used ONLY to read off the
 * signs (pole half-plane, kernel frequency) that decide which residues the
 * contour encloses -- the residue arithmetic itself stays symbolic, so the
 * closed form is correct by construction, with no numeric crosscheck.
 *
 * Returns NULL (declining symbolic-parameter mode, so the integral stays
 * unevaluated) when the integrand has no free parameters, OR when any parameter
 * is left two-sided unbounded by the assumptions: an unconstrained parameter
 * does not determine the sign of a pole that depends on it, so the
 * classification would be a guess.  A one-sided or interval constraint pins the
 * relevant sign on a connected region, which a single interior point reads off
 * correctly. */
static Expr* build_instantiation(Expr* f, Expr* x, Expr* assumptions) {
    ParamBound pb[16];
    size_t np = collect_params(f, x, pb, 16, 0);
    if (np == 0) return NULL;
    if (assumptions) absorb_fact(pb, np, assumptions);

    for (size_t i = 0; i < np; i++)
        if (pb[i].lo <= -HUGE_VAL && pb[i].hi >= HUGE_VAL) return NULL;  /* unconstrained */

    /* All-positive parameters licence a PowerExpand-based simplification of
     * radical pole locations (Sqrt[-4 a^2] -> 2 I a) and real-parameter
     * conjugation, which need the positivity branch. */
    g_all_pos = true;
    for (size_t i = 0; i < np; i++)
        if (pb[i].lo < 0.0) { g_all_pos = false; break; }

    /* Record the guaranteed intervals for the convergence/applicability gates. */
    g_nbounds = np;
    for (size_t i = 0; i < np; i++) g_bounds[i] = pb[i];

    Expr** rules = malloc(np * sizeof(*rules));
    for (size_t i = 0; i < np; i++) {
        double v = pick_representative(pb[i].lo, pb[i].hi, i);
        rules[i] = mk_fn2("Rule", mk_sym(pb[i].sym), expr_new_real(v));
    }
    Expr* lst = expr_new_function(mk_sym("List"), rules, np);
    free(rules);
    return lst;
}

/* Guaranteed [lo, hi] of `e` over the assumed region: the recorded interval for
 * a parameter symbol, the exact value for a numeric constant, else unbounded.
 * A family uses this (not the single instantiation point) to confirm a
 * convergence gate holds everywhere the assumptions permit. */
static void param_interval(Expr* e, double* lo, double* hi) {
    *lo = -HUGE_VAL; *hi = HUGE_VAL;
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < g_nbounds; i++)
            if (g_bounds[i].sym == e->data.symbol.name) {
                *lo = g_bounds[i].lo; *hi = g_bounds[i].hi; return;
            }
    }
    double v;
    if (numeric_double(e, &v)) { *lo = v; *hi = v; }
}

/* -------------------------------------------------------------------------
 * New contour families (implemented in the phases below).
 * ---------------------------------------------------------------------- */
static Expr* residue_family_rectangular(Expr* f, Expr* x, Expr* a, Expr* b);
static Expr* residue_family_mellin(Expr* f, Expr* x, Expr* a, Expr* b);
static Expr* residue_family_sector(Expr* f, Expr* x, Expr* a, Expr* b);

/* -------------------------------------------------------------------------
 * Keyhole / Mellin core -- Integrate[v^p R(v), {v, 0, Infinity}].
 *
 * For a branch power v^p (p non-integer, so s = p+1 is not an integer) times a
 * rational R(v) with poles off the positive real axis, the keyhole contour that
 * wraps the branch cut along [0, Inf) gives the Mellin transform
 *
 *      Int_0^Inf v^(s-1) R(v) dv = -Pi Sum_k [ z_k^(s-1) Res(R, z_k) ] e^(-i Pi s)
 *                                              / Sin(Pi s),          s = p + 1,
 *
 * where z_k ranges over the poles of R and z_k^(s-1) uses the branch arg in
 * (0, 2Pi) (so a pole at z_k is |z_k|^(s-1) e^(i theta_k (s-1)), theta_k in
 * (0,2Pi)).  The e^(-i Pi s)/Sin(Pi s) prefactor is the exact reduction of the
 * keyhole's 1/(1 - e^(2 Pi i s)) jump, which lands the value on a Sin (so a
 * numeric s reduces to an algebraic multiple of Pi).  Poles of any order are
 * summed via mellin_pole_contribution (the residue of the FULL integrand
 * z^(s-1) R(z), not z_k^(s-1) Res(R,z_k) which holds for simple poles only).
 * ---------------------------------------------------------------------- */

/* Split F = v^p R(v) with R rational in v and p a v-free NON-integer exponent.
 * On success *p_out and *R_out are owned; returns false (nothing owned) if there
 * is no such branch power or the remaining factor is not rational in v. */
static bool mellin_split(Expr* F, Expr* v, Expr** p_out, Expr** R_out) {
    *p_out = *R_out = NULL;
    Expr** fac; size_t nf; Expr* one[1];
    if (head_name_is(F, "Times")) { fac = F->data.function.args; nf = F->data.function.arg_count; }
    else { one[0] = F; fac = one; nf = 1; }

    Expr* p = NULL;
    ExprVec rest; ev_init(&rest);
    for (size_t i = 0; i < nf; i++) {
        Expr* fe = fac[i];
        bool is_branch = false;
        if (head_name_is(fe, "Power") && fe->data.function.arg_count == 2) {
            Expr* base = fe->data.function.args[0];
            Expr* e    = fe->data.function.args[1];
            if (base->type == EXPR_SYMBOL && base->data.symbol.name == v->data.symbol.name &&
                !contains_symbol(e, v) && e->type != EXPR_INTEGER) {
                if (p) { expr_free(p); ev_free(&rest); return false; }  /* two branch powers */
                p = expr_copy(e);
                is_branch = true;
            }
        }
        if (!is_branch) ev_push(&rest, expr_copy(fe));
    }
    if (!p) { ev_free(&rest); return false; }

    Expr* R;
    if (rest.n == 0) R = mk_int(1);
    else if (rest.n == 1) { R = rest.v[0]; rest.v[0] = NULL; }
    else {
        Expr** ra = malloc(rest.n * sizeof(*ra));
        for (size_t i = 0; i < rest.n; i++) { ra[i] = rest.v[i]; rest.v[i] = NULL; }
        R = eval_take(expr_new_function(mk_sym("Times"), ra, rest.n));
        free(ra);
    }
    ev_free(&rest);
    if (!R) { expr_free(p); return false; }

    /* R must be rational in v. */
    Expr* P; Expr* Q;
    if (!res_num_den(R, &P, &Q) || !res_polyq(P, v) || !res_polyq(Q, v)) {
        if (P) { expr_free(P); } if (Q) { expr_free(Q); }
        expr_free(p); expr_free(R); return false;
    }
    expr_free(P); expr_free(Q);
    *p_out = p; *R_out = R;
    return true;
}

/* The branch value z^(s-1) with arg in (0, 2Pi): |z|^(s-1) Exp[I theta (s-1)],
 * theta = Arg[z] shifted by 2Pi when it is <= 0 (lands a pole below/on the
 * negative axis on the upper lip).  `pm1` is the exponent s-1 = p. */
static Expr* mellin_branch(Expr* z, Expr* pm1) {
    /* theta = Arg[z], lifted into (0, 2Pi]. */
    double re, im;
    Expr* theta = ev1("Arg", expr_copy(z));
    if (!theta) return NULL;
    if (res_reim(z, &re, &im) && atan2(im, re) <= RES_TOL)
        theta = eval_take(mk_fn2("Plus", theta,
                     mk_fn2("Times", mk_int(2), mk_sym(SYM_Pi))));
    Expr* absz = ev1("Abs", expr_copy(z));
    Expr* mag  = mk_fn2("Power", absz, expr_copy(pm1));       /* |z|^(s-1) */
    Expr* phase = mk_fn1("Exp", mk_fn2("Times", mk_sym(SYM_I),
                       mk_fn2("Times", theta, expr_copy(pm1))));
    return eval_take(mk_fn2("Times", mag, phase));
}

/* Res[ z^(s-1) R(z), z_k ] for a pole z_k of ANY order, off the positive real
 * axis.  The keyhole formula sums the residues of the FULL integrand z^(s-1)R(z);
 * for a simple pole this is z_k^(s-1) Res(R,z_k), but for order >= 2 the two
 * differ (z^(s-1) contributes derivative terms) -- summing z_k^(s-1) Res(R,z_k)
 * there is WRONG (a pure double pole has Res(R)=0, silently zeroing the answer).
 *
 * Shift w = z - z_k:  z^(s-1) = z_k^(s-1) (1 + w/z_k)^(s-1).  The branch value
 * z_k^(s-1) (keyhole arg in (0,2Pi)) is mellin_branch; (1 + w/z_k)^(s-1) is
 * analytic at w=0 (regular binomial series, value 1 there -- NOT a Puiseux
 * expansion in w), so the contribution is
 *      z_k^(s-1) * Res_{w=0}[ (1 + w/z_k)^(s-1) R(w + z_k) ],
 * computed by the ordinary residue engine, correct for any pole order (the
 * simple-pole case falls out as m=1).  pm1 = s-1 = p.  Borrows R, v, z, pm1. */
static Expr* mellin_pole_contribution(Expr* R, Expr* v, Expr* z, Expr* pm1) {
    Expr* B = mellin_branch(z, pm1);                          /* z_k^(s-1), keyhole */
    if (!B) return NULL;
    /* Rshift = R /. v -> v + z  (moves the pole to w = 0). */
    Expr* Rshift = eval_take(mk_fn2("ReplaceAll", expr_copy(R),
                       mk_fn2("Rule", expr_copy(v),
                              mk_fn2("Plus", expr_copy(v), expr_copy(z)))));
    /* analytic = (1 + v/z)^pm1, regular at v = 0. */
    Expr* analytic = mk_fn2("Power",
        mk_fn2("Plus", mk_int(1),
            mk_fn2("Times", expr_copy(v), mk_fn2("Power", expr_copy(z), mk_int(-1)))),
        expr_copy(pm1));
    Expr* prod = eval_take(mk_fn2("Times", analytic, Rshift));
    Expr* zero = mk_int(0);
    Expr* r0 = prod ? residue_compute(prod, v, zero) : NULL;  /* any-order residue */
    if (prod) expr_free(prod);
    expr_free(zero);
    if (!r0) { expr_free(B); return NULL; }
    return eval_take(mk_fn2("Times", B, r0));
}

static Expr* mellin_core(Expr* F, Expr* v) {
    Expr* p; Expr* R;
    if (!mellin_split(F, v, &p, &R)) return NULL;

    /* s = p + 1.  Reject an integer s (Sin[Pi s] = 0: the keyhole degenerates,
     * and an integer power is a job for the even/rational half-line families). */
    Expr* s = eval_take(mk_fn2("Plus", expr_copy(p), mk_int(1)));
    if (!s) { expr_free(p); expr_free(R); return NULL; }
    { double sre, sim;
      if (res_reim(s, &sre, &sim) && fabs(sim) < RES_TOL &&
          fabs(sre - floor(sre + 0.5)) < RES_TOL) {
          expr_free(p); expr_free(R); expr_free(s); return NULL;
      } }

    /* Poles of R: roots of its denominator, excluding v = 0.  Convergence needs
     * 0 < Re(s) < deg(den R) - deg(num R): integrable at 0 (Re s > 0) and decay
     * at Infinity (v^(s-1) v^(dNum-dDen) integrable). */
    Expr* P; Expr* Q;
    if (!res_num_den(R, &P, &Q)) { expr_free(p); expr_free(R); expr_free(s); return NULL; }
    int dNum = res_degree(P, v), dDen = res_degree(Q, v);
    expr_free(P);
    /* Convergence over the whole assumed region: 0 < Re(s) < deg(den)-deg(num),
     * checked against the guaranteed interval (a one-sided assumption like a > 0
     * does not bound s above and must be refused).  s must also be real. */
    { double slo, shi, sre, sim;
      param_interval(s, &slo, &shi);
      if (!res_reim(s, &sre, &sim) || fabs(sim) > RES_TOL ||
          slo < -RES_TOL || shi > (double)(dDen - dNum) + RES_TOL) {
          expr_free(Q); expr_free(p); expr_free(R); expr_free(s); return NULL;
      } }
    ExprVec roots;
    if (!solve_roots(Q, v, &roots)) { expr_free(Q); expr_free(p); expr_free(R); expr_free(s); return NULL; }
    /* dQ = Q'(v), for the simple-pole test (Q'(z_k) != 0  <=>  z_k is a simple
     * pole).  Consumes Q. */
    Expr* dQ = eval_take(mk_fn2("D", Q, expr_copy(v)));

    /* Sigma = Sum_k Res[ z^(s-1) R(z), z_k ] over poles z_k off the positive real
     * axis (a pole on (0,Inf) would sit on the contour -> not handled).  A simple
     * pole takes the fast closed form z_k^(s-1) Res(R, z_k); an order >= 2 pole
     * (where that form is WRONG) takes the general full-integrand residue, which
     * is heavier (a Series of the shifted fractional-power integrand) so it is
     * reserved for exactly the poles that need it. */
    Expr* Sigma = mk_int(0);
    bool bad = false;
    for (size_t i = 0; i < roots.n && !bad; i++) {
        Expr* z = roots.v[i];
        double zre, zim;
        if (!res_reim(z, &zre, &zim)) { bad = true; break; }        /* symbolic pole */
        if (zim > -RES_TOL && zim < RES_TOL && zre > RES_TOL) { bad = true; break; } /* on (0,Inf) */
        /* Simple iff Q'(z) != 0 (numeric z, dQ a polynomial -> numeric value). */
        Expr* dQz = dQ ? eval_take(mk_fn2("ReplaceAll", expr_copy(dQ),
                            mk_fn2("Rule", expr_copy(v), expr_copy(z)))) : NULL;
        double dre, dim;
        bool simple = dQz && res_reim(dQz, &dre, &dim) &&
                      (fabs(dre) > RES_TOL || fabs(dim) > RES_TOL);
        if (dQz) expr_free(dQz);
        Expr* contrib;
        if (simple) {
            Expr* br = mellin_branch(z, p);
            Expr* rr = residue_compute(R, v, z);
            if (!br || !rr) { if (br) expr_free(br); if (rr) expr_free(rr); bad = true; break; }
            contrib = eval_take(mk_fn2("Times", br, rr));
        } else {
            contrib = mellin_pole_contribution(R, v, z, p);
        }
        if (!contrib) { bad = true; break; }
        Sigma = eval_take(mk_fn2("Plus", Sigma, contrib));
        if (!Sigma) { bad = true; break; }
    }
    ev_free(&roots);
    if (dQ) expr_free(dQ);
    if (bad || !Sigma) { if (Sigma) expr_free(Sigma); expr_free(p); expr_free(R); expr_free(s); return NULL; }
    expr_free(R);

    /* value = -Pi Sigma Exp[-I Pi s] / Sin[Pi s]. */
    Expr* raw = mk_fn2("Times",
        mk_fn2("Times", mk_int(-1), mk_sym(SYM_Pi)),
        mk_fn2("Times",
            mk_fn2("Times", Sigma, mk_fn1("Exp",
                mk_fn2("Times", mk_int(-1), mk_fn2("Times", mk_sym(SYM_I),
                    mk_fn2("Times", mk_sym(SYM_Pi), expr_copy(s)))))),
            mk_fn2("Power", mk_fn1("Sin", mk_fn2("Times", mk_sym(SYM_Pi), expr_copy(s))), mk_int(-1))));
    expr_free(s); expr_free(p);

    /* Close: symbolic exponent -> Simplify (lands on Csc[Pi s]); numeric s makes
     * the whole thing an algebraic multiple of Pi -> factor Pi and RootReduce. */
    Expr* value;
    if (g_inst) {
        value = ev1("Simplify", raw);
    } else {
        Expr* alg = ev1("RootReduce",
                        mk_fn2("Times", raw, mk_fn2("Power", mk_sym(SYM_Pi), mk_int(-1))));
        value = alg ? ev1("Simplify", mk_fn2("Times", mk_sym(SYM_Pi), alg)) : NULL;
    }
    if (!value) return NULL;

    double vv, re, im;
    bool ok = res_is_real_scalar(value, v, &vv) || res_is_finite_scalar(value, v, &re, &im);
    if (!ok) { expr_free(value); return NULL; }
    return value;
}

/* Half-line [0, Inf) branch-power integrand v^p R(v): the keyhole/Mellin core. */
static Expr* residue_family_mellin(Expr* f, Expr* x, Expr* a, Expr* b) {
    if (!is_zero_pos_infinity(a, b)) return NULL;
    return mellin_core(f, x);
}

/* Whole-line quasi-periodic integrand f(x) = Exp[c x] R(Exp[x]): the rectangular
 * contour of height 2 Pi i.  Reduced to the keyhole/Mellin core by w = Exp[x]:
 *   Int_{-Inf}^{Inf} f(x) dx = Int_0^Inf f(Log w) / w dw,
 * and f(Log w)/w = w^(c-1) R(w) is exactly the Mellin integrand.  This keeps a
 * single residue engine for both the strip and the branch cut. */
static Expr* residue_family_rectangular(Expr* f, Expr* x, Expr* a, Expr* b) {
    if (!is_neg_pos_infinity(a, b)) return NULL;
    /* Only worth the substitution when Exp[x] actually occurs. */
    if (!contains_head(f, "Exp") &&
        !(contains_symbol(f, x)))  /* cheap guard; real check is the split below */
        return NULL;

    Expr* w = mk_sym("IntegrateResidue`$w");
    /* g(w) = (f /. x -> Log[w]) / w. */
    Expr* fl = eval_take(mk_fn2("ReplaceAll", expr_copy(f),
                                mk_fn2("Rule", expr_copy(x), mk_fn1("Log", expr_copy(w)))));
    if (!fl) { expr_free(w); return NULL; }
    Expr* g = eval_take(mk_fn2("Times", fl, mk_fn2("Power", expr_copy(w), mk_int(-1))));
    if (!g) { expr_free(w); return NULL; }

    /* The substituted integrand must be a genuine branch-power * rational form
     * (a stray Log[w] left over means f was not of the Exp[c x] R(Exp[x]) type). */
    if (contains_head(g, "Log")) { expr_free(g); expr_free(w); return NULL; }
    Expr* value = mellin_core(g, w);
    expr_free(g); expr_free(w);
    return value;
}

/* -------------------------------------------------------------------------
 * Sector contour -- Integrate[x^m / (c + x^n), {x, 0, Infinity}], n possibly a
 * symbolic parameter (n > m + 1 for convergence).  The wedge of angle 2Pi/n maps
 * the ray back to itself scaled by Exp[2 Pi i/n], and the single enclosed pole
 * gives the classical closed form
 *
 *      Int_0^Inf x^(s-1)/(c + x^n) dx = (Pi/n) c^(s/n - 1) / Sin(Pi s/n),
 *      s = m + 1,   c > 0,   0 < Re(s) < n.
 *
 * This is the one family that admits a SYMBOLIC exponent n (the keyhole/Mellin
 * residue sum cannot enumerate n poles), so it is what powers Integrate[1/(1 +
 * x^n), ...].  A monomial numerator C x^m carries a constant factor C through. */

/* Split a monomial C * v^m (C free of v, m a nonnegative integer): on success
 * *C_out is owned, *m_out is the integer exponent.  Returns false otherwise. */
static bool monomial_split(Expr* num, Expr* v, Expr** C_out, int* m_out) {
    if (!contains_symbol(num, v)) { *C_out = expr_copy(num); *m_out = 0; return true; }
    if (num->type == EXPR_SYMBOL && num->data.symbol.name == v->data.symbol.name) {
        *C_out = mk_int(1); *m_out = 1; return true;
    }
    if (head_name_is(num, "Power") && num->data.function.arg_count == 2 &&
        num->data.function.args[0]->type == EXPR_SYMBOL &&
        num->data.function.args[0]->data.symbol.name == v->data.symbol.name &&
        num->data.function.args[1]->type == EXPR_INTEGER) {
        *C_out = mk_int(1);
        *m_out = (int)num->data.function.args[1]->data.integer;
        return true;
    }
    if (head_name_is(num, "Times")) {
        Expr* C = mk_int(1); int m = -1;
        for (size_t i = 0; i < num->data.function.arg_count; i++) {
            Expr* fe = num->data.function.args[i];
            if (!contains_symbol(fe, v)) { C = eval_take(mk_fn2("Times", C, expr_copy(fe))); continue; }
            if (m >= 0) { expr_free(C); return false; }            /* two v-factors */
            if (fe->type == EXPR_SYMBOL && fe->data.symbol.name == v->data.symbol.name) m = 1;
            else if (head_name_is(fe, "Power") &&
                     fe->data.function.args[0]->type == EXPR_SYMBOL &&
                     fe->data.function.args[0]->data.symbol.name == v->data.symbol.name &&
                     fe->data.function.args[1]->type == EXPR_INTEGER)
                m = (int)fe->data.function.args[1]->data.integer;
            else { expr_free(C); return false; }
        }
        if (m < 0) m = 0;
        *C_out = C; *m_out = m; return true;
    }
    return false;
}

static Expr* residue_family_sector(Expr* f, Expr* x, Expr* a, Expr* b) {
    if (!is_zero_pos_infinity(a, b)) return NULL;

    Expr* num; Expr* den;
    if (!res_num_den(f, &num, &den)) return NULL;

    /* Denominator must be exactly c + x^n: a Plus of an x-free constant c and a
     * power x^n (n free of x). */
    if (!head_name_is(den, "Plus") || den->data.function.arg_count != 2) {
        expr_free(num); expr_free(den); return NULL;
    }
    Expr* c = NULL; Expr* n = NULL;
    for (int i = 0; i < 2; i++) {
        Expr* t = den->data.function.args[i];
        if (!contains_symbol(t, x)) c = t;
        else if (head_name_is(t, "Power") && t->data.function.arg_count == 2 &&
                 t->data.function.args[0]->type == EXPR_SYMBOL &&
                 t->data.function.args[0]->data.symbol.name == x->data.symbol.name &&
                 !contains_symbol(t->data.function.args[1], x))
            n = t->data.function.args[1];
    }
    if (!c || !n) { expr_free(num); expr_free(den); return NULL; }

    Expr* C; int m;
    if (!monomial_split(num, x, &C, &m)) { expr_free(num); expr_free(den); return NULL; }

    /* Convergence / positivity gates that must hold over the WHOLE assumed
     * region (guaranteed interval bounds, not the instantiation point): c > 0
     * and n > m + 1 (so 0 < s = m+1 < n).  n must also be real. */
    double clo, chi, nlo, nhi, nim, nre;
    param_interval(c, &clo, &chi);
    param_interval(n, &nlo, &nhi);
    bool okc = clo > RES_TOL;                                   /* c > 0 guaranteed */
    bool okn = nlo >= (double)(m + 1) &&                        /* n > m+1 guaranteed */
               res_reim(n, &nre, &nim) && fabs(nim) < RES_TOL;  /* n real */
    if (!okc || !okn) { expr_free(C); expr_free(num); expr_free(den); return NULL; }

    /* value = C (Pi/n) c^(s/n - 1) Csc[Pi s/n], s = m + 1. */
    Expr* s   = mk_int(m + 1);
    Expr* son = mk_fn2("Times", expr_copy(s), mk_fn2("Power", expr_copy(n), mk_int(-1)));  /* s/n */
    Expr* cpow = mk_fn2("Power", expr_copy(c),
                        mk_fn2("Plus", expr_copy(son), mk_int(-1)));                        /* c^(s/n-1) */
    Expr* csc = mk_fn2("Power", mk_fn1("Sin", mk_fn2("Times", mk_sym(SYM_Pi), expr_copy(son))), mk_int(-1));
    Expr* raw = mk_fn2("Times", C,
                    mk_fn2("Times", mk_fn2("Times", mk_sym(SYM_Pi), mk_fn2("Power", expr_copy(n), mk_int(-1))),
                           mk_fn2("Times", cpow, csc)));
    expr_free(s); expr_free(son);
    expr_free(num); expr_free(den);

    Expr* value = ev1("Simplify", raw);
    double vv, re, im;
    bool ok = value && (res_is_real_scalar(value, x, &vv) || res_is_finite_scalar(value, x, &re, &im));
    if (!ok) { if (value) expr_free(value); return NULL; }
    return value;
}

/* -------------------------------------------------------------------------
 * Master entry + builtin.
 * ---------------------------------------------------------------------- */

Expr* integrate_residue_try(Expr* f, Expr* x, Expr* a, Expr* b,
                            Expr* assumptions, bool* diverges) {
    if (diverges) *diverges = false;
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;

    /* Enter symbolic-parameter mode when assumptions are supplied and we are the
     * outermost call (g_inst not already set by an enclosing call -- the even
     * half-line family recurses through here).  Only the frame that built the
     * instantiation frees it. */
    bool built_here = false;
    if (assumptions && !g_inst) {
        g_inst = build_instantiation(f, x, assumptions);
        built_here = (g_inst != NULL);
    }

    Expr* value = NULL;
    if (is_full_period(a, b)) {
        /* Full period (0,2Pi)/(-Pi,Pi): unit-circle trig-substitution family. */
        value = residue_family_trig(f, x, a, b, diverges);
    } else if (is_neg_pos_infinity(a, b)) {
        /* Whole line: Fourier/Jordan if a trig/exp kernel is present, then the
         * quasi-periodic (rectangular-contour) family, else rational. */
        value = residue_family_fourier(f, x, a, b);
        if (!value) value = residue_family_rectangular(f, x, a, b);
        if (!value) value = residue_family_rational(f, x, a, b, diverges);
    } else if (is_zero_pos_infinity(a, b)) {
        /* Half-line [0,Inf): even symmetry first (exact), then the keyhole/Mellin
         * (branch-power / log) and sector (1/(1+x^n)) families. */
        value = residue_half_line(f, x, a, b, diverges);
        if (!value) value = residue_family_mellin(f, x, a, b);
        if (!value) value = residue_family_sector(f, x, a, b);
    }

    if (built_here) { expr_free(g_inst); g_inst = NULL; g_nbounds = 0; }
    return value;
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
    return integrate_residue_try(f, x, a, b, NULL, NULL);
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
        "half-line [0,Infinity) variants.  Sums the enclosed residues; the value "
        "is correct by construction (no numeric crosscheck) and returned "
        "unevaluated when no family applies or the residue sum does not close to "
        "a scalar.  Symbolic-parameter integrals reach the keyhole/Mellin, sector "
        "and rectangular families through Integrate[f, spec, Assumptions -> ...].");
}
