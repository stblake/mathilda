/* integrate_goursat.c
 *
 * Goursat's pseudo-elliptic algorithm and its cube-/fourth-root
 * generalisations.  See integrate_goursat.h for the mathematical overview and
 * GoursatAppendix.wl for the reference implementation this file transcribes.
 *
 * ---------------------------------------------------------------------
 * Strategy
 * ---------------------------------------------------------------------
 * The reductions live one tier above the arithmetic primitives: this file is a
 * C *orchestration* of Mathilda's existing high-level builtins -- Solve (with
 * Cubics/Quartics radical solving), Together/Cancel over algebraic-number
 * extensions (Extension -> Automatic), PossibleZeroQ, Numerator/Denominator,
 * Expand -- via the eval_take(mk_fn...) idiom already used by
 * integrate_linratiorad.c.  A verbatim port of the WL appendix is impossible
 * (it relies on Association, which Mathilda lacks), so the control flow is
 * re-expressed in C while every algebraic step delegates to a builtin.
 *
 * 1. Recognition: split f into a rational cofactor F(t) and a single radical
 *    R(t)^q with q in {-1/2,-1/3,-2/3,-1/4,-3/4}; R must be a polynomial.
 * 2. Solve R == 0 for the roots (radicals via Cubics/Quartics -> True).
 * 3. Build the order-(den p) Mobius automorphism cycling those roots, project
 *    the integrand into eigencomponents, and test the elementarity criterion
 *    of the relevant theorem with PossibleZeroQ.
 * 4. On success, descend each non-obstructive eigenpiece to a genus-0 curve,
 *    integrate the resulting rational function recursively, and back-substitute.
 *    On failure (obstructed / roots not in radicals / a piece does not close)
 *    return NULL so the Integrate cascade continues.
 *
 * Correct by construction: like the sibling radical integrators, no
 * differentiate-back verification is performed -- the eigenspace criterion plus
 * recursive closure are the guards.
 *
 * Memory: builtins take ownership of `res`; helpers below own every Expr they
 * construct and free intermediates explicitly.  eval_take / *_ext / canonic /
 * subst_eval consume their primary argument.  We never expr_free(res).  An
 * Infinity root / fixed point is represented by a NULL Expr* sentinel.
 */

#include "integrate_goursat.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "internal.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "poly.h"
#include "zero_test.h"
#include "print.h"

#include <stdbool.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------- */
/* Small builders / evaluation helpers (mirror integrate_chebychev.c)     */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* name) { return expr_new_symbol(name); }

static Expr* mk_fn1(const char* name, Expr* a) {
    return expr_new_function(mk_sym(name), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(name), (Expr*[]){ a, b }, 2);
}
static Expr* mk_fn3(const char* name, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(mk_sym(name), (Expr*[]){ a, b, c }, 3);
}
static Expr* mk_fnv(const char* name, Expr** args, size_t n) {
    return expr_new_function(mk_sym(name), args, n);
}

static Expr* mk_pow(Expr* base, Expr* exp) { return mk_fn2("Power", base, exp); }
static Expr* mk_inv(Expr* e) { return mk_pow(e, mk_int(-1)); }
static Expr* mk_neg(Expr* e) { return mk_fn2("Times", mk_int(-1), e); }
static Expr* mk_rat(int64_t p, int64_t q) { return make_rational(p, q); }
static Expr* mk_pow_rat(Expr* base, int64_t n, int64_t d) {
    return mk_fn2("Power", base, make_rational(n, d));
}
static Expr* mk_pow_int(Expr* base, int64_t n) {
    return mk_fn2("Power", base, mk_int(n));
}
static Expr* mk_sqrt_expr(Expr* e) { return mk_pow_rat(e, 1, 2); }
static Expr* mk_rule(Expr* lhs, Expr* rhs) { return mk_fn2("Rule", lhs, rhs); }

/* Readable arithmetic builders (each consumes its operands). */
static Expr* gadd(Expr* a, Expr* b) { return mk_fn2("Plus", a, b); }
static Expr* gsub(Expr* a, Expr* b) { return mk_fn2("Plus", a, mk_neg(b)); }
static Expr* gmul(Expr* a, Expr* b) { return mk_fn2("Times", a, b); }
static Expr* gdiv(Expr* a, Expr* b) { return mk_fn2("Times", a, mk_inv(b)); }

/* Evaluate `call` to a fixed point, freeing `call`. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* True if `f` contains no subexpression structurally equal to `x`. */
static bool expr_free_of(const Expr* f, const Expr* x) {
    if (expr_eq((Expr*)f, (Expr*)x)) return false;
    if (f->type == EXPR_FUNCTION) {
        if (!expr_free_of(f->data.function.head, x)) return false;
        for (size_t i = 0; i < f->data.function.arg_count; i++)
            if (!expr_free_of(f->data.function.args[i], x)) return false;
    }
    return true;
}

/* True if `e` contains any unevaluated Integrate[...] call. */
static bool contains_unintegrated(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Integrate)) return true;
    if (contains_unintegrated(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_unintegrated(e->data.function.args[i])) return true;
    return false;
}

/* Together[e, Extension -> Automatic]; consumes e. */
static Expr* together_ext(Expr* e) {
    if (!e) return NULL;
    Expr* opt = mk_rule(mk_sym(SYM_Extension), mk_sym(SYM_Automatic));
    return eval_take(internal_together((Expr*[]){ e, opt }, 2));
}
/* Cancel[e, Extension -> Automatic]; consumes e. */
static Expr* cancel_ext(Expr* e) {
    if (!e) return NULL;
    Expr* opt = mk_rule(mk_sym(SYM_Extension), mk_sym(SYM_Automatic));
    return eval_take(internal_cancel((Expr*[]){ e, opt }, 2));
}
/* canonic[e] = Cancel[Together[e, Extension -> Automatic], Extension ->
 * Automatic]; consumes e. */
static Expr* canonic(Expr* e) {
    return cancel_ext(together_ext(e));
}

/* Expand[e]; consumes e. */
static Expr* expand_e(Expr* e) {
    if (!e) return NULL;
    return eval_take(internal_expand((Expr*[]){ e }, 1));
}
/* Simplify[e]; consumes e.  Needed to collapse roots-of-unity identities
 * (e.g. 1 + (-1)^(2/3) - (-1)^(1/3) -> 0) that Cancel/Together over the
 * Extension field leaves intact. */
static Expr* simplify_e(Expr* e) {
    if (!e) return NULL;
    return eval_take(mk_fn1("Simplify", e));
}
/* PowerExpand[e]; consumes e.  Collapses Sqrt of a perfect square (e.g.
 * Sqrt[k^2] -> k) under the positive-branch convention the integrator already
 * uses for radicands.  Applied to the V4 fixed points so a parameter-dependent
 * involution such as t -> 1/(k^2 t), whose fixed points Solve returns as the
 * spurious nested radicals +-Sqrt[k^2]/k^2, reduces to the intended +-1/k --
 * otherwise Sqrt[k^2] rides through lc/Q/gx as a SECOND algebraic generator and
 * the genus-0 integrand never closes.  A no-op on genuinely nested radicals
 * (Sqrt[1-k^2]) and on numeric constants; the final answer is still checked by
 * diff_back_ok. */
static Expr* powerexpand_e(Expr* e) {
    if (!e) return NULL;
    return eval_take(mk_fn1("PowerExpand", e));
}
/* Numerator[e] / Denominator[e]; borrow e. */
static Expr* numer_of(Expr* e) {
    return eval_take(internal_numerator((Expr*[]){ expr_copy(e) }, 1));
}
static Expr* denom_of(Expr* e) {
    return eval_take(internal_denominator((Expr*[]){ expr_copy(e) }, 1));
}

/* ReplaceAll[e, from -> to] then evaluate; consumes e and to, borrows from.
 * (Pass expr_copy(to) when the replacement value is reused.) */
static Expr* subst_eval(Expr* e, const Expr* from, Expr* to) {
    if (!e) { if (to) expr_free(to); return NULL; }
    Expr* rule = mk_rule(expr_copy((Expr*)from), to);
    return eval_take(internal_replace_all((Expr*[]){ e, rule }, 2));
}

/* Re-express every power of the monic radicand Q in G as a power of the true
 * radicand H = lc Q, keeping lc as a constant factor:
 *      G /. Power[Q, p_] :> Power[lc Q, p] Power[lc, -p].
 * Algebraically an identity (lc^p Q^p lc^-p = Q^p), so it leaves G a valid
 * antiderivative -- but it changes the FORM.  We integrate gx/Sqrt[Q] with a
 * MONIC Q (keeps the recursive Euler reduction rational and fast); naively
 * restoring the leading coefficient via a Sqrt[lc] prefactor is WRONG because
 * Sqrt[lc] Sqrt[Q] != Sqrt[lc Q] when Q < 0 (a branch-cut sign flip yielding a
 * bogus antiderivative).  Rewriting instead makes the radical the true H, which
 * after back-substitution resolves to the clean real Sqrt[R]*rational, and lc
 * survives only as a constant where it cancels the prefactor's 1/Sqrt[lc].
 * Consumes G; borrows Q, lc. */
static Expr* rebase_radical(Expr* G, Expr* Q, Expr* lc) {
    if (!G) return NULL;
    Expr* blank = mk_fnv("Blank", NULL, 0);
    Expr* ppat  = mk_fn2("Pattern", mk_sym("p"), blank);
    Expr* lhs   = mk_pow(expr_copy(Q), ppat);
    Expr* rhs   = mk_fn2("Times",
        mk_pow(mk_fn2("Times", expr_copy(lc), expr_copy(Q)), mk_sym("p")),
        mk_pow(expr_copy(lc), mk_neg(mk_sym("p"))));
    Expr* rule  = mk_fn2("RuleDelayed", lhs, rhs);
    return eval_take(internal_replace_all((Expr*[]){ G, rule }, 2));
}

/* PossibleZeroQ via the internal three-valued decision: strict TRUE only. */
static bool is_zero(const Expr* e) {
    if (!e) return false;
    return zero_test_decide(e) == ZERO_TEST_TRUE;
}

/* Recurse into the full integrator: Integrate[g, u].  Returns the closed
 * antiderivative, or NULL if it does not close.  Borrows g and u. */
static Expr* integrate_in(const Expr* g, const Expr* u) {
    Expr* r = eval_take(mk_fn2("Integrate", expr_copy((Expr*)g),
                               expr_copy((Expr*)u)));
    if (!r) return NULL;
    if (contains_unintegrated(r)) { expr_free(r); return NULL; }
    return r;
}

/* True iff e is the symbol True. */
static bool expr_is_true(Expr* e) {
    if (!e) return false;
    Expr* t = mk_sym("True");
    bool r = expr_eq(e, t);
    expr_free(t);
    return r;
}

/* ---------------------------------------------------------------------- */
/* Debug tracing -- controlled by the global Integrate`GoursatDebug        */
/* ---------------------------------------------------------------------- */
/*
 * When the user sets Integrate`GoursatDebug = True the descent narrates its
 * progress to stderr: whether the integrand matches the expected pseudo-elliptic
 * form, which involution / eigenspace criterion is tested and whether it holds,
 * and the recursive genus-0 reductions.  The flag is read ONCE at the outermost
 * entry (gs_set_debug, called from gs_guarded at gs_depth == 0) into the static
 * gs_debug so the recursive descents -- which run inside a separate
 * TimeConstrained sub-evaluation -- share the setting without re-evaluating the
 * symbol per log site.  Off by default: every log site is a single int test.
 */
static int gs_debug = 0;
static int gs_depth;   /* recursion depth, defined below; borrowed for indenting */

/* Read the current value of Integrate`GoursatDebug (True -> 1, else 0). */
static void gs_set_debug(void) {
    Expr* v = evaluate(mk_sym("Integrate`GoursatDebug"));
    gs_debug = expr_is_true(v) ? 1 : 0;
    if (v) expr_free(v);
}

/* Emit a formatted trace line (indented by the current recursion depth). */
static void gs_log(const char* fmt, ...) {
    if (!gs_debug) return;
    fprintf(stderr, "[Goursat]%*s ", 2 * gs_depth, "");
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* Emit a trace line "label: <expr>" with the expression pretty-printed. */
static void gs_log_expr(const char* label, Expr* e) {
    if (!gs_debug) return;
    char* s = e ? expr_to_string(e) : NULL;
    gs_log("%s: %s", label, s ? s : (e ? "<unprintable>" : "<null>"));
    if (s) free(s);
}

/* Differentiate-back verification: is D[result, x] - f identically zero?
 *
 * The ideal check is PossibleZeroQ[D[result,x] - f], and verification must NOT
 * use Simplify.  But PossibleZeroQ currently has two failure modes on this
 * subsystem's outputs (see POSSIBLE_ZEROQ_FAILURES.md): it MISFIRES (returns a
 * false negative) on the residual cyclotomic-tower nested radicals of the
 * period-3 higher-symmetry cases across a branch cut, and it is too SLOW on the
 * very large parametric antiderivatives to fit the descent's CPU budget.  Until
 * PossibleZeroQ handles both, verify numerically: pin the free parameters
 * (variables other than x) to fixed rationals, then sample x and require the
 * magnitude of D[result,x] - f to vanish (Abs collapses the R<0 complex branch).
 * Accept iff enough sampled points are ~0 and none is plainly non-zero.  Borrows
 * result, x, f. */
static bool diff_back_ok(Expr* result, Expr* x, Expr* f) {
    Expr* d = eval_take(mk_fn2("D", expr_copy(result), expr_copy(x)));
    if (!d) return false;
    Expr* diff = gsub(d, expr_copy(f));        /* D[result,x] - f */
    if (!diff) return false;

    /* Pin free parameters (e.g. the k in R = x(1-x)(1-k^2 x)) so the check is a
     * decisive single-variable numeric test rather than a symbolic one. */
    Expr* dd = expr_copy(diff);
    {
        Expr* vars = eval_take(mk_fn1("Variables", expr_copy(diff)));
        if (vars && head_is(vars, SYM_List)) {
            static const long an[8] = { 12, 17, 23, 29, 31, 37, 41, 43 };
            static const long ad[8] = {  7,  5,  9, 11, 13, 10,  6,  8 };
            size_t j = 0;
            for (size_t i = 0; dd && i < vars->data.function.arg_count; i++) {
                Expr* v = vars->data.function.args[i];
                if (expr_eq(v, x)) continue;
                dd = subst_eval(dd, v, make_rational(an[j % 8], ad[j % 8])); j++;
            }
        }
        if (vars) expr_free(vars);
    }
    expr_free(diff);
    if (!dd) return false;

    static const long pn[8] = { 17, 23, 31, 29, 37, 41, 19, 43 };
    static const long pd[8] = {  5,  7,  9, 11, 10, 12,  6,  8 };
    Expr* tol_lo = make_rational(1, 1000000000000L);   /* 1e-12: accept below */
    Expr* tol_hi = make_rational(1, 1000000L);         /* 1e-6 : reject above */
    int good = 0; bool bad = false;
    for (int i = 0; i < 8 && good < 3 && !bad; i++) {
        Expr* at  = subst_eval(expr_copy(dd), x, make_rational(pn[i], pd[i]));
        Expr* mag = at ? eval_take(mk_fn2("N", mk_fn1("Abs", at), mk_int(30))) : NULL;
        if (mag) {
            Expr* lo = eval_take(mk_fn2("Less",    expr_copy(mag), expr_copy(tol_lo)));
            Expr* hi = eval_take(mk_fn2("Greater", mag,            expr_copy(tol_hi)));
            if (expr_is_true(lo)) good++;
            else if (expr_is_true(hi)) bad = true;
            if (lo) expr_free(lo);
            if (hi) expr_free(hi);
        }
    }
    expr_free(tol_lo); expr_free(tol_hi); expr_free(dd);
    return good >= 2 && !bad;
}

/* Cheap numeric "is e definitely NOT identically zero in x?" test.
 *
 * The eigenspace elementarity criteria require certain projections of F to
 * vanish.  Combining those projections symbolically (canonic = Together/Cancel
 * over the splitting field of R) blows up super-polynomially when R has
 * cyclotomic roots, so the obstructed (non-elementary) integrands -- which are
 * exactly the ones that SHOULD be rejected -- are the slowest to reject.
 *
 * Decide it numerically on the *uncombined* expression instead: sample e at a
 * spread of fixed rational points and N-evaluate |e| (Abs collapses the complex
 * values that the cyclotomic radicals in e take to a real magnitude).  Return
 * true as soon as one finite sample is decisively above the noise floor -- proof
 * that e is not the zero function, so the criterion fails and the caller
 * declines.  A false return ("every sampled point is ~0") only authorises the
 * descent, whose antiderivative is independently verified by diff_back_ok, so a
 * numeric false-positive-for-zero cannot produce a wrong result.  Borrows e, x. */
static bool sample_clearly_nonzero(Expr* e, Expr* x) {
    if (!e) return false;
    static const long pn[8] = {  7, 13, 31, 23, 37, 19, 41, 29 };
    static const long pd[8] = { 10,  5, 10,  7, 10,  6, 12, 11 };
    Expr* tol_lo = make_rational(1, 1000000L);          /* 1e-6 : above noise floor */
    Expr* tol_hi = mk_int(1000000000L);                 /* 1e9  : reject pole blow-ups */
    bool nonzero = false;
    for (int i = 0; i < 8 && !nonzero; i++) {
        Expr* pt  = make_rational(pn[i], pd[i]);
        Expr* at  = subst_eval(expr_copy(e), x, pt);   /* consumes pt */
        Expr* mag = at ? eval_take(mk_fn2("N", mk_fn1("Abs", at), mk_int(30))) : NULL;
        if (mag) {
            /* Count a point as evidence of non-vanishing only when |e| is both
             * above the noise floor AND finite/bounded -- a sample that lands on
             * (or near) a pole of some f_j numericalises to Infinity even though
             * the symbolic projection is identically zero, so such points must be
             * skipped rather than mistaken for a genuine non-zero. */
            Expr* lo = eval_take(mk_fn2("Greater", expr_copy(mag), expr_copy(tol_lo)));
            Expr* hi = eval_take(mk_fn2("Less",    expr_copy(mag), expr_copy(tol_hi)));
            if (expr_is_true(lo) && expr_is_true(hi)) nonzero = true;
            if (lo) expr_free(lo);
            if (hi) expr_free(hi);
            expr_free(mag);
        }
    }
    expr_free(tol_lo); expr_free(tol_hi);
    return nonzero;
}

/* ---------------------------------------------------------------------- */
/* Recursion guard + fresh-symbol counter                                 */
/* ---------------------------------------------------------------------- */

#define GS_MAX_DEPTH 6
static int gs_depth = 0;
static unsigned long gs_sym_counter = 0;

static Expr* fresh_var(void) {
    char b[96];
    snprintf(b, sizeof(b), "Integrate`GoursatAlgebraic`u$%lu", gs_sym_counter++);
    return expr_new_symbol(b);
}

/* ---------------------------------------------------------------------- */
/* Root solving                                                           */
/* ---------------------------------------------------------------------- */

/* Solve[eq, t, Cubics -> True, Quartics -> True] and collect the distinct root
 * values (RHS of each Rule[t, val]) into a freshly-malloc'd, owned array.
 * Borrows eq, t.  Returns NULL with *n = 0 on no solutions / failure. */
static Expr** solve_eq_roots(Expr* eq, Expr* t, size_t* n) {
    *n = 0;
    Expr* call = mk_fnv("Solve", (Expr*[]){
        expr_copy(eq), expr_copy(t),
        mk_rule(mk_sym(SYM_Cubics),   mk_sym(SYM_True)),
        mk_rule(mk_sym(SYM_Quartics), mk_sym(SYM_True)) }, 4);
    Expr* sol = eval_take(call);
    if (!sol || !head_is(sol, SYM_List)) { if (sol) expr_free(sol); return NULL; }

    size_t cap = sol->data.function.arg_count;
    Expr** out = cap ? (Expr**)malloc(sizeof(Expr*) * cap) : NULL;
    size_t cnt = 0;
    for (size_t i = 0; i < sol->data.function.arg_count; i++) {
        Expr* sub = sol->data.function.args[i];   /* List[Rule[t, val]] */
        if (!head_is(sub, SYM_List) || sub->data.function.arg_count != 1) continue;
        Expr* rule = sub->data.function.args[0];
        if (!head_is(rule, SYM_Rule) || rule->data.function.arg_count != 2) continue;
        Expr* val = rule->data.function.args[1];
        bool dup = false;
        for (size_t k = 0; k < cnt; k++)
            if (expr_eq(out[k], val)) { dup = true; break; }
        if (!dup) out[cnt++] = expr_copy(val);
    }
    expr_free(sol);
    *n = cnt;
    return out;
}

/* Distinct roots of polynomial R in t. */
static Expr** poly_roots(Expr* R, Expr* t, size_t* n) {
    Expr* eq = mk_fn2("Equal", expr_copy(R), mk_int(0));
    Expr** rs = solve_eq_roots(eq, t, n);
    expr_free(eq);
    return rs;
}

/* Fixed points of a Mobius map S(t): solve S == t.  Writes alpha (always
 * owned) and beta (owned, or NULL meaning the second fixed point is at
 * Infinity -- i.e. S == t has a single finite solution).  Returns false on
 * an unexpected (non Mobius) solution count.  Borrows S, t. */
static bool fixed_points(Expr* S, Expr* t, Expr** alpha, Expr** beta) {
    Expr* eq = mk_fn2("Equal", expr_copy(S), expr_copy(t));
    size_t nf = 0;
    Expr** fps = solve_eq_roots(eq, t, &nf);
    expr_free(eq);
    *alpha = NULL; *beta = NULL;
    if (nf == 1) { *alpha = fps[0]; *beta = NULL; free(fps); return true; }
    if (nf >= 2) {
        *alpha = fps[0]; *beta = fps[1];
        for (size_t i = 2; i < nf; i++) expr_free(fps[i]);
        free(fps);
        return true;
    }
    if (fps) free(fps);
    return false;
}

/* ---------------------------------------------------------------------- */
/* eigenpiece -> rational function of x = z^m                             */
/* ---------------------------------------------------------------------- */

/* ToFunctionOf{Square,Cube,Fourth}: H is a rational function of z whose
 * numerator and denominator involve only multiples-of-m powers of z; return
 * the rational function of x obtained by z^m -> x.  Returns Integer 0 when H is
 * zero.  Borrows H, z; returns a freshly-owned, evaluated expression. */
static Expr* to_function_of_power(Expr* H, Expr* z, Expr* x, int m) {
    Expr* Hc = canonic(expr_copy(H));
    if (!Hc) return NULL;
    if (is_zero(Hc)) { expr_free(Hc); return mk_int(0); }

    Expr* num = expand_e(numer_of(Hc));
    Expr* den = expand_e(denom_of(Hc));
    expr_free(Hc);
    if (!num || !den) { if (num) expr_free(num); if (den) expr_free(den); return NULL; }

    int dN = get_degree_poly(num, z);
    int dD = get_degree_poly(den, z);

    /* numerator polynomial in x: sum_k Coefficient[num, z, m k] x^k */
    Expr* np = mk_int(0);
    for (int k = 0; k * m <= dN; k++) {
        Expr* ck = get_coeff(num, z, m * k);
        if (!ck) ck = mk_int(0);
        np = mk_fn2("Plus", np, mk_fn2("Times", ck, mk_pow_int(expr_copy(x), k)));
    }
    Expr* dp = mk_int(0);
    for (int k = 0; k * m <= dD; k++) {
        Expr* ck = get_coeff(den, z, m * k);
        if (!ck) ck = mk_int(0);
        dp = mk_fn2("Plus", dp, mk_fn2("Times", ck, mk_pow_int(expr_copy(x), k)));
    }
    expr_free(num); expr_free(den);
    return canonic(mk_fn2("Times", np, mk_inv(dp)));
}

/* ---------------------------------------------------------------------- */
/* Square-root case: Goursat's V4 theorem (p = 1/2)                       */
/* ---------------------------------------------------------------------- */

/* The Mobius involution swapping {a,b} with {c,d}.  Any slot may be NULL
 * (Infinity).  Borrows a,b,c,d,t; returns a freshly-owned canonic rational
 * function of t.  (Appendix MobiusInvolution.) */
static Expr* mobius_involution(Expr* a, Expr* b, Expr* c, Expr* d, Expr* t) {
    if (a == NULL) return mobius_involution(c, d, b, a, t);
    if (b == NULL) return mobius_involution(c, d, a, b, t);
    if (c == NULL) return mobius_involution(a, b, d, c, t);
    if (d == NULL) {
        /* (c t + a b - c (a+b)) / (t - c) */
        Expr* num = mk_fn3("Plus",
            mk_fn2("Times", expr_copy(c), expr_copy(t)),
            mk_fn2("Times", expr_copy(a), expr_copy(b)),
            mk_fn3("Times", mk_int(-1), expr_copy(c),
                   mk_fn2("Plus", expr_copy(a), expr_copy(b))));
        Expr* den = mk_fn2("Plus", expr_copy(t), mk_neg(expr_copy(c)));
        return canonic(mk_fn2("Times", num, mk_inv(den)));
    }
    /* ((ab-cd) t + (a+b) c d - (c+d) a b) / (((a+b)-(c+d)) t - (ab-cd)) */
    Expr* ab = mk_fn2("Times", expr_copy(a), expr_copy(b));
    Expr* cd = mk_fn2("Times", expr_copy(c), expr_copy(d));
    Expr* apb = mk_fn2("Plus", expr_copy(a), expr_copy(b));
    Expr* cpd = mk_fn2("Plus", expr_copy(c), expr_copy(d));
    Expr* ab_m_cd = mk_fn2("Plus", expr_copy(ab), mk_neg(expr_copy(cd)));

    Expr* num = mk_fn3("Plus",
        mk_fn2("Times", expr_copy(ab_m_cd), expr_copy(t)),
        mk_fn2("Times", expr_copy(apb), expr_copy(cd)),
        mk_neg(mk_fn2("Times", expr_copy(cpd), expr_copy(ab))));
    Expr* den = mk_fn2("Plus",
        mk_fn2("Times",
            mk_fn2("Plus", expr_copy(apb), mk_neg(expr_copy(cpd))),
            expr_copy(t)),
        mk_neg(expr_copy(ab_m_cd)));
    expr_free(ab); expr_free(cd); expr_free(apb); expr_free(cpd); expr_free(ab_m_cd);
    return canonic(mk_fn2("Times", num, mk_inv(den)));
}

/* ---------------------------------------------------------------------- */
/* Cyclotomic-tower fixed-point denesting                                  */
/* ---------------------------------------------------------------------- */
/*
 * `Solve[S == t]` returns the V4 involution's fixed points via the quadratic
 * formula, so they carry a nested radical `Sqrt[disc]` whose discriminant `disc`
 * is a CYCLOTOMIC CONSTANT (S's coefficients live in the splitting field of R).
 * Left opaque, that `Sqrt[disc]` is an out-of-field tower element: substituting
 * such a fixed point into `t(u)` leaves `R(t(u))·(1-u)^4` un-reducible (the
 * radical reads as an independent variable), so the descended `lc`/`Q`/`gx` keep
 * spurious `(1-u)` denominators and `Sqrt[cyclotomic]` coefficients, and the
 * genus-0 piece never closes -- `integrate_in` re-recognises a cubic-radical
 * Goursat problem and the descent recurses until the depth/time guard fires.
 *
 * But `disc` is, by construction, a rational multiple of a root of unity (e.g.
 * `12 (-1)^(2/3)`), whose square root IS a tower element expressible over the
 * generators (`Sqrt[12 (-1)^(2/3)] = 2 Sqrt[3] (-1)^(1/3)`).  Denest it: detect
 * `disc = c (-1)^(p/q)` numerically (verified exactly), then `PowerExpand` the
 * monomial.  The fixed points become explicit cyclotomic-tower elements and the
 * whole descent stays rational over that field.
 */

/* True iff e contains an inexact (Real / MPFR) leaf. */
static bool contains_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef EXPR_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_FUNCTION) {
        if (contains_inexact(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (contains_inexact(e->data.function.args[i])) return true;
    }
    return false;
}

/* If constant b (free of var) equals a rational multiple of a root of unity,
 * c*(-1)^(p/q), return that exact monomial form; else NULL.  Borrows b, var. */
static Expr* monomialize_root_of_unity(Expr* b, Expr* var) {
    if (!expr_free_of(b, var)) return NULL;
    Expr* nb = eval_take(mk_fn1("N", expr_copy(b)));
    if (!nb) return NULL;
    /* This denesting is valid ONLY for a genuine numeric cyclotomic constant.
     * With a symbolic PARAMETER present (e.g. the k in R = x(1-x)(1-k^2 x)) the
     * base is free of the integration variable yet still depends on k, so N[b]
     * retains Abs[..]/Arg[..]/k and the Rationalize[Arg[b]/Pi] heuristic below
     * fabricates a spurious rational approximation of 1/Pi (113/355) times
     * Arg[k^2] -- injecting Arg/Abs junk into the descent that then never
     * closes.  Require N[b] to reduce to an actual number before proceeding. */
    Expr* isnum = eval_take(mk_fn1("NumberQ", expr_copy(nb)));
    bool numeric = expr_is_true(isnum);
    if (isnum) expr_free(isnum);
    if (!numeric) { expr_free(nb); return NULL; }
    Expr* tol = make_rational(1, 1000000);
    Expr* c  = eval_take(mk_fn2("Rationalize", mk_fn1("Abs", expr_copy(nb)),
                                expr_copy(tol)));
    Expr* qq = eval_take(mk_fn2("Rationalize",
                   gdiv(mk_fn1("Arg", nb), mk_sym("Pi")), tol));  /* consumes nb, tol */
    if (!c || !qq) { if (c) expr_free(c); if (qq) expr_free(qq); return NULL; }
    Expr* cand = eval_take(mk_fn2("Times", c, mk_pow(mk_int(-1), qq)));  /* consumes c, qq */
    if (!cand) return NULL;
    /* Require an EXACT candidate (Rationalize must have produced rationals) and
     * exact equality with b -- never inject an inexact fixed point. */
    if (contains_inexact(cand)) { expr_free(cand); return NULL; }
    Expr* chk = gsub(expr_copy(b), expr_copy(cand));
    bool ok = is_zero(chk);
    expr_free(chk);
    if (!ok) { expr_free(cand); return NULL; }
    return cand;
}

/* Collect distinct Power[base, m/2] (m odd, base free of var) subexpressions. */
static void collect_const_halfpowers(Expr* e, Expr* var,
                                     Expr*** out, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        bool half = false;
        for (int m = 1; m <= 9 && !half; m += 2) {
            Expr* r = make_rational(m, 2);
            half = expr_eq(exp, r);
            expr_free(r);
        }
        if (half && expr_free_of(e->data.function.args[0], var)) {
            bool dup = false;
            for (size_t i = 0; i < *n; i++) if (expr_eq((*out)[i], e)) { dup = true; break; }
            if (!dup) {
                if (*n == *cap) {
                    *cap = *cap ? *cap * 2 : 4;
                    *out = (Expr**)realloc(*out, *cap * sizeof(Expr*));
                }
                (*out)[(*n)++] = e;
            }
        }
    }
    collect_const_halfpowers(e->data.function.head, var, out, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect_const_halfpowers(e->data.function.args[i], var, out, n, cap);
}

/* Denest every constant cyclotomic half-power in e (see block comment).
 * Returns a new owned tree, or NULL when nothing was denested.  Borrows e, var. */
static Expr* denest_const_radicals(Expr* e, Expr* var) {
    Expr** nodes = NULL; size_t nn = 0, cap = 0;
    collect_const_halfpowers(e, var, &nodes, &nn, &cap);
    Expr* out = NULL;
    for (size_t i = 0; i < nn; i++) {
        Expr* p = nodes[i];                       /* Power[base, exp] */
        Expr* mono = monomialize_root_of_unity(p->data.function.args[0], var);
        if (!mono) continue;
        Expr* repl = eval_take(mk_fn1("PowerExpand",
                         mk_pow(mono, expr_copy(p->data.function.args[1]))));
        if (!repl) continue;
        Expr* rule = mk_rule(expr_copy(p), repl);   /* consumes repl */
        Expr* base_e = out ? out : e;
        Expr* nx = eval_take(internal_replace_all(
                       (Expr*[]){ expr_copy(base_e), rule }, 2));
        if (out) expr_free(out);
        out = nx;
    }
    free(nodes);
    return out;
}

/* Reduce one V4 eigenpiece Fj (anti-invariant under involution S) to its
 * elementary antiderivative.  Returns the owned antiderivative (in t), or NULL
 * if the rational reduction does not close.  Borrows Fj, R, t, S. */
static Expr* reduce_v4_piece(Expr* Fj, Expr* R, Expr* t, Expr* S) {
    Expr* alpha = NULL; Expr* beta = NULL;
    if (!fixed_points(S, t, &alpha, &beta)) return NULL;
    /* Denest cyclotomic-tower nested radicals in the fixed points so the descent
     * stays rational over the splitting field (else it recurses; see block
     * comment above). */
    { Expr* ad = denest_const_radicals(alpha, t);
      if (ad) { expr_free(alpha); alpha = ad; } }
    if (beta) { Expr* bd = denest_const_radicals(beta, t);
                if (bd) { expr_free(beta); beta = bd; } }
    /* Collapse Sqrt[perfect-square] fixed points (e.g. Sqrt[k^2]/k^2 -> 1/k) so a
     * parameter-dependent involution's descent stays over Q(param) with the single
     * radical Sqrt[Q] instead of carrying Sqrt[k^2] as a spurious 2nd generator. */
    { Expr* ap = powerexpand_e(expr_copy(alpha));
      if (ap) { expr_free(alpha); alpha = ap; } }
    if (beta) { Expr* bp = powerexpand_e(expr_copy(beta));
                if (bp) { expr_free(beta); beta = bp; } }
    bool binf = (beta == NULL);
    gs_log_expr("    reduce_v4_piece: alpha", alpha);
    if (beta) gs_log_expr("    reduce_v4_piece: beta", beta);
    else gs_log("    reduce_v4_piece: beta = Infinity");

    Expr* u = fresh_var();     /* intermediate Mobius variable */
    Expr* x = fresh_var();     /* output square variable        */
    Expr* result = NULL;

    /* tu : t as a function of u. */
    Expr* tu = binf
        ? mk_fn2("Plus", expr_copy(alpha), expr_copy(u))
        : canonic(mk_fn2("Times",
            mk_fn2("Plus", expr_copy(alpha),
                   mk_neg(mk_fn2("Times", expr_copy(beta), expr_copy(u)))),
            mk_inv(mk_fn2("Plus", mk_int(1), mk_neg(expr_copy(u))))));

    /* Rfact = R(t(u)) [ (1-u)^4 ] expanded. */
    Expr* Rsub = subst_eval(expr_copy(R), t, expr_copy(tu));
    Expr* Rfact;
    if (binf) {
        Rfact = expand_e(Rsub);
    } else {
        Expr* om4 = mk_pow_int(mk_fn2("Plus", mk_int(1), mk_neg(expr_copy(u))), 4);
        Rfact = expand_e(canonic(mk_fn2("Times", Rsub, om4)));
    }
    if (!Rfact) goto done;

    Expr* lc = get_coeff(Rfact, u, 4);
    if (!lc || is_zero(lc)) { if (lc) expr_free(lc); expr_free(Rfact); goto done; }

    /* Q(x) = sum_{k=0}^{2} (Coefficient[Rfact,u,2k]/lc) x^k -- the MONIC even
     * part of Rfact as a polynomial in x = u^2.  Integrating gx/Sqrt[Q] with a
     * monic radicand keeps the recursive Euler reduction rational and fast.  The
     * true radicand is H = lc Q; Sqrt[lc] is restored after integration by
     * rebase_radical (NOT factored into the prefactor -- that is branch-unsafe,
     * since Sqrt[lc] Sqrt[Q] != Sqrt[lc Q] when Q < 0). */
    Expr* Q = mk_int(0);
    for (int k = 0; k <= 2; k++) {
        Expr* c2k = get_coeff(Rfact, u, 2 * k);
        if (!c2k) c2k = mk_int(0);
        Q = mk_fn2("Plus", Q,
            mk_fn3("Times", c2k, mk_inv(expr_copy(lc)), mk_pow_int(expr_copy(x), k)));
    }
    Q = canonic(Q);
    expr_free(Rfact);

    /* gu = Fj(t(u)) / u ; gx = ToFunctionOfSquare[gu]. */
    Expr* Fu = canonic(subst_eval(expr_copy(Fj), t, expr_copy(tu)));
    Expr* gu = canonic(mk_fn2("Times", Fu, mk_inv(expr_copy(u))));
    Expr* gx = gu ? to_function_of_power(gu, u, x, 2) : NULL;
    if (gu) expr_free(gu);
    if (!gx || !Q) {
        if (gx) expr_free(gx); if (Q) expr_free(Q); expr_free(lc); goto done;
    }

    /* prefactor = (binf ? 1 : alpha-beta) / (2 Sqrt[lc]). */
    Expr* pre_num = binf ? mk_int(1)
                         : mk_fn2("Plus", expr_copy(alpha), mk_neg(expr_copy(beta)));
    Expr* pre = mk_fn2("Times", pre_num,
        mk_inv(mk_fn2("Times", mk_int(2), mk_sqrt_expr(expr_copy(lc)))));

    /* Integrate gx / Sqrt[Q] in x (Q monic -> exact, fast rational reduction). */
    Expr* Qkeep = expr_copy(Q);
    Expr* integrand = mk_fn2("Times", gx, mk_pow_rat(Q, -1, 2));  /* consumes gx,Q */
    gs_depth++;
    Expr* G = integrate_in(integrand, x);
    gs_depth--;
    gs_log("    reduce_v4_piece: genus-0 integrate_in %s", G ? "CLOSED" : "did NOT close");
    expr_free(integrand);
    if (!G) { expr_free(pre); expr_free(Qkeep); expr_free(lc); goto done; }

    /* Rewrite Sqrt[Q] -> Sqrt[lc Q]/Sqrt[lc] so the radical is the true H = lc Q
     * (Sqrt[R]*rational after back-substitution) and Sqrt[lc] cancels the
     * prefactor's 1/Sqrt[lc] as the real constant lc -- never straddling the
     * branch cut on the variable radicand. */
    G = rebase_radical(G, Qkeep, lc);
    expr_free(Qkeep); expr_free(lc);
    if (!G) { expr_free(pre); goto done; }

    /* Back-substitute x -> (binf ? (t-alpha)^2 : ((t-alpha)/(t-beta))^2). */
    Expr* ta = mk_fn2("Plus", expr_copy(t), mk_neg(expr_copy(alpha)));
    Expr* back = binf
        ? mk_pow_int(ta, 2)
        : mk_pow_int(mk_fn2("Times", ta,
              mk_inv(mk_fn2("Plus", expr_copy(t), mk_neg(expr_copy(beta))))), 2);
    Expr* Gs = subst_eval(G, x, back);   /* consumes back */
    result = eval_take(mk_fn2("Times", pre, Gs));

done:
    if (alpha) expr_free(alpha);
    if (beta)  expr_free(beta);
    if (tu)    expr_free(tu);
    expr_free(u); expr_free(x);
    return result;
}

/* Count expression nodes (a cheap structural-complexity proxy). */
static long node_count(const Expr* e) {
    if (!e) return 0;
    long n = 1;
    if (e->type == EXPR_FUNCTION) {
        n += node_count(e->data.function.head);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            n += node_count(e->data.function.args[i]);
    }
    return n;
}

/* Score an involution by the cleanliness of its (denested) fixed points: each
 * V4 projection is anti-invariant under TWO involutions, and the descent over
 * whichever we pick must stay in that involution's fixed-point field.  A
 * complex (cyclotomic-tower) fixed point drags the reduction into a larger
 * tower whose nested radicals do not all denest -- so heavily prefer REAL fixed
 * points, then break ties by structural size.  Lower is better.  Borrows S, t. */
static long involution_score(Expr* S, Expr* t) {
    Expr* al = NULL; Expr* be = NULL;
    if (!fixed_points(S, t, &al, &be)) return 1000000000L;
    long score = 0;
    Expr* fps[2] = { al, be };
    for (int i = 0; i < 2; i++) {
        if (!fps[i]) continue;
        Expr* d = denest_const_radicals(fps[i], t);
        Expr* f = d ? d : fps[i];
        Expr* im  = eval_take(mk_fn2("N", mk_fn1("Im", expr_copy(f)), mk_int(20)));
        Expr* big = im ? eval_take(mk_fn2("Greater", mk_fn1("Abs", im),
                                   make_rational(1, 1000000000L))) : NULL;
        if (big && expr_is_true(big)) score += 1000000;   /* non-real: penalise */
        if (big) expr_free(big);
        score += node_count(f);
        if (d) expr_free(d);
    }
    if (al) expr_free(al);
    if (be) expr_free(be);
    return score;
}

/* Goursat V4 (p = 1/2): F rational, R cubic/quartic with simple roots. */
static Expr* goursat_v4(Expr* F, Expr* R, Expr* t) {
    size_t nr = 0;
    Expr** roots = poly_roots(R, t, &nr);
    int dR = get_degree_poly(R, t);

    /* Need four ramification points; a cubic contributes Infinity as the
     * fourth.  Require simple roots (distinct-root count matches degree). */
    Expr* r[4] = { NULL, NULL, NULL, NULL };
    bool ok = false;
    if (dR == 4 && nr == 4) {
        for (int i = 0; i < 4; i++) r[i] = roots[i];
        ok = true;
    } else if (dR == 3 && nr == 3) {
        for (int i = 0; i < 3; i++) r[i] = roots[i];
        r[3] = NULL;  /* Infinity */
        ok = true;
    }
    if (!ok) {
        for (size_t i = 0; i < nr; i++) expr_free(roots[i]);
        if (roots) free(roots);
        return NULL;
    }

    /* Three root-pairing involutions. */
    Expr* S1 = mobius_involution(r[0], r[1], r[2], r[3], t);
    Expr* S2 = mobius_involution(r[0], r[2], r[1], r[3], t);
    Expr* S3 = mobius_involution(r[0], r[3], r[1], r[2], t);
    Expr* S[3] = { S1, S2, S3 };

    /* V4 projections of F.  Build the four substituted copies f_j = F(S_j(t))
     * but DEFER the costly canonic: combining four rational functions over the
     * cyclotomic splitting field of R via Together explodes super-polynomially
     * (e.g. R = x^3 - 1, whose roots are the cube roots of unity), and the
     * trivial-character projection P0 is only needed for a zero test. */
    Expr* f0 = expr_copy(F);
    Expr* f1 = subst_eval(expr_copy(F), t, expr_copy(S1));
    Expr* f2 = subst_eval(expr_copy(F), t, expr_copy(S2));
    Expr* f3 = subst_eval(expr_copy(F), t, expr_copy(S3));

    Expr* P[4] = { NULL, NULL, NULL, NULL };
    Expr* result = NULL;

    /* Trivial-character projection P0 = (f0+f1+f2+f3)/4 must vanish, else the
     * integral is obstructed (non-elementary).  Decide this on the *uncombined*
     * sum by numeric sampling: a clearly non-zero sample declines instantly,
     * sidestepping the cyclotomic Together blowup that a symbolic is_zero(P0)
     * would trigger.  A passing (zero-ish) sample only authorises the descent,
     * whose result is independently checked by the diff_back_ok guard in
     * gs_core -- so a numeric false-positive cannot yield a wrong answer. */
    gs_log("  V4 criterion: trivial-character projection P0 = (f0+f1+f2+f3)/4 must vanish");
    if (f1 && f2 && f3) {
        Expr* P0raw = mk_fn3("Plus", mk_fn2("Plus", expr_copy(f0), expr_copy(f1)),
                                     expr_copy(f2), expr_copy(f3));
        bool obstructed = sample_clearly_nonzero(P0raw, t);
        expr_free(P0raw);
        if (obstructed) {
            gs_log("  P0 is decisively non-zero -- OBSTRUCTED for V4, declining");
            expr_free(f0); expr_free(f1); expr_free(f2); expr_free(f3);
            goto cleanup;
        }
    }
    gs_log("  P0 looks vanishing -- descending the non-trivial V4 projections");

    /* Goursat Theorem 1 fast path.  When F is ALREADY anti-invariant under one
     * of the three involutions (F + F(S_k) == 0), a single direct reduction
     * reduce_v4_piece(F, S_k) yields the antiderivative -- the character split
     * (Theorem 2) is unnecessary and actively harmful here: it decomposes such
     * an F into two more-complex pieces F_j = (F -/+ F(S_i))/2 whose genus-0
     * integrands are non-constant, so the two reductions bloat and (in this
     * integrator, lacking a full Simplify) never recombine.  Reducing F
     * directly keeps the genus-0 integrand as simple as the integrand permits
     * (a constant, for the parametric (1+k x)/((k x-1) Sqrt[R]) family), giving
     * the compact single-term answer.  f0=F, f1=F(S1), f2=F(S2), f3=F(S3) are
     * still live here; the anti-invariance test is a pure rational zero test
     * (no radicals), so is_zero is reliable.  Fall back to the Theorem-2 split
     * if no single involution works or the direct reduction fails to close. */
    {
        Expr* fS[3] = { f1, f2, f3 };
        int t1 = -1; long t1sc = 0;
        for (int k = 0; k < 3; k++) {
            if (!fS[k]) continue;
            Expr* s = canonic(mk_fn2("Plus", expr_copy(f0), expr_copy(fS[k])));
            bool anti = s && is_zero(s);
            if (s) expr_free(s);
            if (anti) {
                long sc = involution_score(S[k], t);
                gs_log("  Theorem 1: F is anti-invariant under S[%d] (score=%ld)", k, sc);
                if (t1 < 0 || sc < t1sc) { t1 = k; t1sc = sc; }
            }
        }
        if (t1 >= 0) {
            gs_log("  Theorem 1: reducing F directly under S[%d] -- no character split", t1);
            Expr* piece = reduce_v4_piece(F, R, t, S[t1]);
            if (piece) {
                gs_log_expr("  Theorem-1 piece antiderivative", piece);
                result = eval_take(piece);
                expr_free(f0); expr_free(f1); expr_free(f2); expr_free(f3);
                goto cleanup;
            }
            gs_log("  Theorem-1 direct reduction did not close -- falling back to Theorem-2 split");
        }
    }

    /* P0 looks vanishing -- now canonicalise the three non-trivial projections
     * that the descent actually integrates. */
    P[1] = canonic(mk_fn3("Plus", mk_fn2("Plus", expr_copy(f0), expr_copy(f1)),
                                  mk_neg(expr_copy(f2)), mk_neg(expr_copy(f3))));
    P[2] = canonic(mk_fn3("Plus", mk_fn2("Plus", expr_copy(f0), mk_neg(expr_copy(f1))),
                                  expr_copy(f2), mk_neg(expr_copy(f3))));
    P[3] = canonic(mk_fn3("Plus", mk_fn2("Plus", expr_copy(f0), mk_neg(expr_copy(f1))),
                                  mk_neg(expr_copy(f2)), expr_copy(f3)));
    for (int i = 1; i <= 3; i++)
        if (P[i]) P[i] = canonic(mk_fn2("Times", mk_rat(1, 4), P[i]));
    expr_free(f0); expr_free(f1); expr_free(f2); expr_free(f3);

    /* Sum the reductions of the non-zero projections P[1],P[2],P[3].  P[j] is
     * anti-invariant under the two involutions S[k], k != j-1; prefer one whose
     * S == t has a single finite fixed point (Infinity in the other slot). */
    Expr* total = mk_int(0);
    for (int j = 1; j <= 3; j++) {
        if (!P[j] || is_zero(P[j])) { gs_log("  projection P[%d] is zero -- skip", j); continue; }
        gs_log_expr("  nonzero projection", P[j]);
        /* candidate anti-involutions: indices in {0,1,2} \ {j-1}.  Pick the one
         * with the cleanest (preferably real) fixed points, so the descent stays
         * in the smallest cyclotomic tower and its nested radicals all denest. */
        int chosen = -1; long best = 0;
        for (int k = 0; k < 3; k++) {
            if (k == j - 1) continue;
            long sc = involution_score(S[k], t);
            gs_log("    candidate involution S[%d] score=%ld", k, sc);
            if (chosen < 0 || sc < best) { chosen = k; best = sc; }
        }
        gs_log("  chose involution S[%d] (score=%ld) for P[%d]", chosen, best, j);
        gs_log_expr("    chosen S", S[chosen]);
        Expr* piece = reduce_v4_piece(P[j], R, t, S[chosen]);
        if (piece) gs_log_expr("  piece antiderivative", piece);
        else gs_log("  piece FAILED to close (reduce_v4_piece -> NULL)");
        if (!piece) { expr_free(total); total = NULL; break; }
        total = mk_fn2("Plus", total, piece);
    }
    if (total) result = eval_take(total);

cleanup:
    for (int i = 0; i < 4; i++) if (P[i]) expr_free(P[i]);
    expr_free(S1); expr_free(S2); expr_free(S3);
    for (size_t i = 0; i < nr; i++) expr_free(roots[i]);
    if (roots) free(roots);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Cube-root case: order-3 Mobius eigendescent (p = 1/3 and 2/3)          */
/* ---------------------------------------------------------------------- */

/* omega = Exp[2 Pi I/3], a primitive cube root of unity. */
static Expr* mk_omega(void) {
    return mk_fn1("Exp", mk_fn3("Times", mk_rat(2, 3), mk_sym("Pi"), mk_sym("I")));
}

/* The order-3 Mobius transformation cycling the three roots r[0]->r[1]->r[2].
 * A NULL entry denotes Infinity.  Borrows r[*], t; returns owned canonic S(t).
 * (Appendix CyclicMobius.) */
static Expr* cyclic_mobius(Expr** r, Expr* t) {
    int pinf = -1;
    for (int i = 0; i < 3; i++) if (r[i] == NULL) pinf = i;

    if (pinf >= 0) {
        /* Two finite roots, in cyclic order starting after Infinity. */
        Expr* r1 = r[(pinf + 1) % 3];
        Expr* r2 = r[(pinf + 2) % 3];
        /* S = (r1 t - (r1^2 - r1 r2 + r2^2)) / (t - r2) */
        Expr* quad = gsub(gadd(mk_pow_int(expr_copy(r1), 2),
                               mk_pow_int(expr_copy(r2), 2)),
                          gmul(expr_copy(r1), expr_copy(r2)));
        Expr* num = gsub(gmul(expr_copy(r1), expr_copy(t)), quad);
        Expr* den = gsub(expr_copy(t), expr_copy(r2));
        return canonic(gdiv(num, den));
    }

    /* Three finite roots: solve {(A r_i + B)/(C r_i + 1) == r_{i+1}}. */
    Expr* A = fresh_var(); Expr* B = fresh_var(); Expr* C = fresh_var();
    Expr* eqs[3];
    for (int i = 0; i < 3; i++) {
        Expr* ri  = r[i];
        Expr* rip = r[(i + 1) % 3];
        /* A ri + B == r_{i+1} (C ri + 1) */
        Expr* lhs = gadd(gmul(expr_copy(A), expr_copy(ri)), expr_copy(B));
        Expr* rhs = gmul(expr_copy(rip),
                         gadd(gmul(expr_copy(C), expr_copy(ri)), mk_int(1)));
        eqs[i] = mk_fn2("Equal", lhs, rhs);
    }
    Expr* eqlist = mk_fnv("List", eqs, 3);
    Expr* vars = mk_fn3("List", expr_copy(A), expr_copy(B), expr_copy(C));
    Expr* solset = eval_take(mk_fn2("Solve", eqlist, vars));

    Expr* S = NULL;
    if (solset && head_is(solset, SYM_List) && solset->data.function.arg_count >= 1) {
        Expr* sol = solset->data.function.args[0];   /* List[Rule[A,_],...] */
        Expr* Sexpr = gdiv(gadd(gmul(expr_copy(A), expr_copy(t)), expr_copy(B)),
                           gadd(gmul(expr_copy(C), expr_copy(t)), mk_int(1)));
        Expr* Ssub = eval_take(internal_replace_all(
            (Expr*[]){ Sexpr, expr_copy(sol) }, 2));
        S = canonic(Ssub);
    }
    if (solset) expr_free(solset);
    expr_free(A); expr_free(B); expr_free(C);
    return S;
}

/* P_k projection under z -> omega z: (H + w^-k H(wz) + w^-2k H(w^2 z))/3.
 * Borrows H, z; returns owned canonic.  (Appendix EigenProjection.) */
static Expr* eigenproj3(Expr* H, Expr* z, int k) {
    Expr* w  = mk_omega();
    Expr* wz  = gmul(expr_copy(w), expr_copy(z));
    Expr* w2z = gmul(mk_pow_int(expr_copy(w), 2), expr_copy(z));
    Expr* Hw  = subst_eval(expr_copy(H), z, wz);    /* consumes wz  */
    Expr* Hw2 = subst_eval(expr_copy(H), z, w2z);   /* consumes w2z */
    Expr* t0 = expr_copy(H);
    Expr* t1 = gmul(mk_pow_int(expr_copy(w), -k), Hw);
    Expr* t2 = gmul(mk_pow_int(expr_copy(w), -2 * k), Hw2);
    expr_free(w);
    /* Simplify before canonic: collapses 1 + omega + omega^2 -> 0 so that a
     * vanishing eigencomponent is structurally zero, not a fake nonzero
     * constant that breaks the z^3-coefficient extraction downstream. */
    return canonic(simplify_e(gdiv(mk_fn3("Plus", t0, t1, t2), mk_int(3))));
}

/* Raw (evaluate-only, NO canonic/Simplify) order-3 eigenprojection, used solely
 * by the numeric criterion gate: the costly canonic/Simplify over the cyclotomic
 * field is exactly the blowup we want to avoid before deciding the integral is
 * obstructed, and numeric sampling does not need a canonical form.  Borrows H,z. */
static Expr* eigenproj3_raw(Expr* H, Expr* z, int k) {
    Expr* w  = mk_omega();
    Expr* wz  = gmul(expr_copy(w), expr_copy(z));
    Expr* w2z = gmul(mk_pow_int(expr_copy(w), 2), expr_copy(z));
    Expr* Hw  = subst_eval(expr_copy(H), z, wz);
    Expr* Hw2 = subst_eval(expr_copy(H), z, w2z);
    Expr* t0 = expr_copy(H);
    Expr* t1 = gmul(mk_pow_int(expr_copy(w), -k), Hw);
    Expr* t2 = gmul(mk_pow_int(expr_copy(w), -2 * k), Hw2);
    expr_free(w);
    return eval_take(gdiv(mk_fn3("Plus", t0, t1, t2), mk_int(3)));
}

/* Integrate one J-piece in uOut and back-substitute uOut -> backval.  Returns
 * Integer 0 for a zero piece, NULL on non-closure, else the owned contribution
 * in t.  Borrows J, uOut, backval. */
static Expr* integ_backsub(Expr* J, Expr* uOut, Expr* backval) {
    if (!J || is_zero(J)) return mk_int(0);
    gs_depth++;
    Expr* G = integrate_in(J, uOut);
    gs_depth--;
    if (!G) return NULL;
    return subst_eval(G, uOut, expr_copy(backval));   /* backval borrowed */
}

/* Cube-root Goursat (p = pnum/3, pnum in {1,2}): R cubic (or quadratic). */
static Expr* goursat_cubic(Expr* F, Expr* R, Expr* t, int pnum) {
    size_t nr = 0;
    Expr** roots = poly_roots(R, t, &nr);
    int dR = get_degree_poly(R, t);

    Expr* r[3] = { NULL, NULL, NULL };
    bool ok = false;
    if (dR == 3 && nr == 3) {
        for (int i = 0; i < 3; i++) r[i] = roots[i];
        ok = true;
    } else if (dR == 2 && nr == 2) {
        r[0] = roots[0]; r[1] = roots[1]; r[2] = NULL;  /* Infinity */
        ok = true;
    }
    if (!ok) {
        for (size_t i = 0; i < nr; i++) expr_free(roots[i]);
        if (roots) free(roots);
        return NULL;
    }

    Expr* result = NULL;
    Expr* z = NULL; Expr* x = NULL; Expr* uOut = NULL;
    Expr* alpha = NULL; Expr* beta = NULL;
    Expr* tz = NULL; Expr* H = NULL;
    Expr* H0 = NULL; Expr* H1 = NULL; Expr* H2 = NULL;
    Expr* cval = NULL; Expr* K = NULL;

    Expr* S = cyclic_mobius(r, t);
    if (!S) goto cleanup;
    if (!fixed_points(S, t, &alpha, &beta)) { expr_free(S); goto cleanup; }
    bool binf = (beta == NULL);

    z = fresh_var(); x = fresh_var(); uOut = fresh_var();

    /* t(z). */
    tz = binf
        ? gadd(expr_copy(alpha), expr_copy(z))
        : canonic(gdiv(gsub(expr_copy(alpha), gmul(expr_copy(beta), expr_copy(z))),
                       gsub(mk_int(1), expr_copy(z))));

    /* R(t(z)) [ (1-z)^3 ]. */
    Expr* Rsub = subst_eval(expr_copy(R), t, expr_copy(tz));
    Expr* Rz = binf ? canonic(Rsub)
                    : canonic(gmul(Rsub,
                          mk_pow_int(gsub(mk_int(1), expr_copy(z)), 3)));
    Expr* Rze = expand_e(Rz);
    if (!Rze) { expr_free(S); goto cleanup; }
    cval = get_coeff(Rze, z, 3);
    Expr* c0 = get_coeff(Rze, z, 0);
    expr_free(Rze);
    if (!cval || !c0 || is_zero(cval)) {
        if (cval) { expr_free(cval); cval = NULL; }
        if (c0) expr_free(c0);
        expr_free(S); goto cleanup;
    }
    K = canonic(gdiv(mk_neg(c0), expr_copy(cval)));   /* -c0/cval */

    /* H (p=1/3) or Htilde (p=2/3).  Build the raw (un-canonicalised) integrand
     * source `base` first; canonic is deferred until after the numeric gate. */
    Expr* base = NULL;
    {
        Expr* Ftz = subst_eval(expr_copy(F), t, expr_copy(tz));
        Expr* cpow = (pnum == 1) ? mk_pow_rat(expr_copy(cval), -1, 3)
                                 : mk_pow_rat(expr_copy(cval), -2, 3);
        if (pnum == 1) {
            base = binf
                ? gmul(Ftz, cpow)
                : gmul(gmul(gsub(expr_copy(alpha), expr_copy(beta)), Ftz),
                       gmul(cpow, mk_inv(gsub(mk_int(1), expr_copy(z)))));
        } else {
            base = binf
                ? gmul(Ftz, cpow)
                : gmul(gsub(expr_copy(alpha), expr_copy(beta)), gmul(Ftz, cpow));
        }
    }
    expr_free(S);
    if (!base) goto cleanup;

    /* Numeric elementarity gate (p=1/3 needs H1==0, p=2/3 needs H0==0): sample
     * the RAW criterion eigenprojection (no canonic/Simplify) and decline at
     * once when it is decisively non-zero -- this is the cyclotomic-field
     * Together/Simplify blowup that the obstructed cases would otherwise hit.
     * A passing sample falls through to the exact same symbolic path as before;
     * gs_core's diff_back_ok independently verifies any antiderivative. */
    gs_log("  cube criterion: omega-eigencomponent H%d of H must vanish",
           (pnum == 1) ? 1 : 0);
    {
        Expr* crit = eigenproj3_raw(base, z, (pnum == 1) ? 1 : 0);
        bool obstructed = sample_clearly_nonzero(crit, z);
        if (crit) expr_free(crit);
        if (obstructed) {
            gs_log("  numeric gate: H%d is decisively non-zero -- OBSTRUCTED (non-elementary), declining",
                   (pnum == 1) ? 1 : 0);
            expr_free(base); goto cleanup;
        }
    }
    gs_log("  numeric gate: criterion eigencomponent looks vanishing -- proceeding to exact descent");

    H = canonic(base);
    if (!H) goto cleanup;

    H0 = eigenproj3(H, z, 0);
    H1 = eigenproj3(H, z, 1);
    H2 = eigenproj3(H, z, 2);

    /* Elementarity criterion: p=1/3 needs H1==0, p=2/3 needs H0==0. */
    if (pnum == 1) {
        if (!H1 || !is_zero(H1)) {
            gs_log("  exact criterion: H1 != 0 -- criterion FAILS, declining");
            goto cleanup;
        }
    } else {
        if (!H0 || !is_zero(H0)) {
            gs_log("  exact criterion: H0 != 0 -- criterion FAILS, declining");
            goto cleanup;
        }
    }
    gs_log("  exact criterion HOLDS -- descending eigenpieces to genus-0");

    /* Common back-substitution atoms. */
    Expr* cube_c  = mk_pow_rat(expr_copy(cval), 1, 3);          /* cval^(1/3) */
    Expr* cube_R  = mk_pow_rat(expr_copy(R), 1, 3);             /* R^(1/3)    */
    Expr* amb     = binf ? NULL : gsub(expr_copy(alpha), expr_copy(beta));

    Expr* total = NULL;
    if (pnum == 1) {
        /* phi0 = ToFunctionOfCube[H0]; psi2 = ToFunctionOfCube[H2/z^2]. */
        Expr* phi0 = to_function_of_power(H0, z, x, 3);
        Expr* H2z2 = canonic(gdiv(expr_copy(H2), mk_pow_int(expr_copy(z), 2)));
        Expr* psi2 = H2z2 ? to_function_of_power(H2z2, z, x, 3) : NULL;
        if (H2z2) expr_free(H2z2);

        /* J0 = (phi0 /. x->K/(1-u^3)) u/(1-u^3);  J2 = (psi2 /. x->u^3+K) u. */
        Expr* om3 = gsub(mk_int(1), mk_pow_int(expr_copy(uOut), 3));   /* 1-u^3 */
        Expr* J0 = phi0 ? canonic(gmul(
            subst_eval(phi0, x, gdiv(expr_copy(K), expr_copy(om3))),
            gdiv(expr_copy(uOut), expr_copy(om3)))) : NULL;
        expr_free(om3);
        Expr* J2 = NULL;
        if (psi2) {
            Expr* u3K = gadd(mk_pow_int(expr_copy(uOut), 3), expr_copy(K));
            J2 = canonic(gmul(subst_eval(psi2, x, u3K), expr_copy(uOut)));
        }

        /* Back-substitutions. */
        Expr* J0back = binf
            ? gdiv(expr_copy(cube_R), gmul(expr_copy(cube_c),
                   gsub(expr_copy(t), expr_copy(alpha))))
            : gdiv(gmul(expr_copy(cube_R), expr_copy(amb)),
                   gmul(expr_copy(cube_c), gsub(expr_copy(t), expr_copy(alpha))));
        Expr* J2back = binf
            ? gdiv(expr_copy(cube_R), expr_copy(cube_c))
            : gdiv(gmul(expr_copy(cube_R), expr_copy(amb)),
                   gmul(expr_copy(cube_c), gsub(expr_copy(t), expr_copy(beta))));

        Expr* p0 = integ_backsub(J0, uOut, J0back);
        Expr* p2 = (p0 != NULL) ? integ_backsub(J2, uOut, J2back) : NULL;
        if (J0) expr_free(J0);
        if (J2) expr_free(J2);
        expr_free(J0back); expr_free(J2back);
        if (p0 && p2) total = gadd(p0, p2);
        else { if (p0) expr_free(p0); if (p2) expr_free(p2); }
    } else {
        /* phi1 = ToFunctionOfCube[H1/z]; phi2 = ToFunctionOfCube[H2/z^2]. */
        Expr* H1z = canonic(gdiv(expr_copy(H1), expr_copy(z)));
        Expr* phi1 = H1z ? to_function_of_power(H1z, z, x, 3) : NULL;
        if (H1z) expr_free(H1z);
        Expr* H2z2 = canonic(gdiv(expr_copy(H2), mk_pow_int(expr_copy(z), 2)));
        Expr* phi2 = H2z2 ? to_function_of_power(H2z2, z, x, 3) : NULL;
        if (H2z2) expr_free(H2z2);

        /* J1 = -(phi1 /. x->K u^3/(u^3-1)) u/(u^3-1);  J2 = phi2 /. x->u^3+K. */
        Expr* u3m1 = gsub(mk_pow_int(expr_copy(uOut), 3), mk_int(1));   /* u^3-1 */
        Expr* J1 = phi1 ? canonic(mk_neg(gmul(
            subst_eval(phi1, x, gdiv(gmul(expr_copy(K), mk_pow_int(expr_copy(uOut), 3)),
                                     expr_copy(u3m1))),
            gdiv(expr_copy(uOut), expr_copy(u3m1))))) : NULL;
        expr_free(u3m1);
        Expr* J2 = NULL;
        if (phi2) {
            Expr* u3K2 = gadd(mk_pow_int(expr_copy(uOut), 3), expr_copy(K));
            J2 = canonic(subst_eval(phi2, x, u3K2));
        }

        Expr* J1back = binf
            ? gdiv(gmul(expr_copy(cube_c), gsub(expr_copy(t), expr_copy(alpha))),
                   expr_copy(cube_R))
            : gdiv(gmul(expr_copy(cube_c), gsub(expr_copy(t), expr_copy(alpha))),
                   gmul(expr_copy(cube_R), expr_copy(amb)));
        Expr* J2back = binf
            ? gdiv(expr_copy(cube_R), expr_copy(cube_c))
            : gdiv(gmul(expr_copy(cube_R), expr_copy(amb)),
                   gmul(expr_copy(cube_c), gsub(expr_copy(t), expr_copy(beta))));

        Expr* p1 = integ_backsub(J1, uOut, J1back);
        Expr* p2 = (p1 != NULL) ? integ_backsub(J2, uOut, J2back) : NULL;
        if (J1) expr_free(J1);
        if (J2) expr_free(J2);
        expr_free(J1back); expr_free(J2back);
        if (p1 && p2) total = gadd(p1, p2);
        else { if (p1) expr_free(p1); if (p2) expr_free(p2); }
    }

    expr_free(cube_c); expr_free(cube_R); if (amb) expr_free(amb);
    if (total) result = eval_take(total);

cleanup:
    if (z) expr_free(z); if (x) expr_free(x); if (uOut) expr_free(uOut);
    if (alpha) expr_free(alpha); if (beta) expr_free(beta);
    if (tz) expr_free(tz); if (H) expr_free(H);
    if (H0) expr_free(H0); if (H1) expr_free(H1); if (H2) expr_free(H2);
    if (cval) expr_free(cval); if (K) expr_free(K);
    for (size_t i = 0; i < nr; i++) expr_free(roots[i]);
    if (roots) free(roots);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Third-kind cube-root logarithmic part (Goursat cube-root, p = 1/3)      */
/* ---------------------------------------------------------------------- */
/*
 * When the order-3 eigendescent criterion FAILS (the omega^1 eigencomponent
 * H1 != 0) the integrand F(t)/R(t)^(1/3) is still elementary if the obstruction
 * differential is of the THIRD kind -- i.e. F has a pole at a NON-branch point
 * of the cube-root cover (a pole not at a root of R).  Its antiderivative is the
 * logarithmic part, a sum of logs of the "linear-on-the-curve" functions
 * y - omega^j kappa t (y = R^(1/3), kappa = (lead R)^(1/3), omega = e^(2 pi i/3)):
 *
 *     Integral F/R^(1/3) dt  =  C * Sum_{j=0..2} omega^(2 j) Log[R^(1/3) - omega^j kappa t].
 *
 * The three log arguments are canonical: their product (norm over the sheets)
 * is  y^3 - kappa^3 t^3 = R(t) - (lead R) t^3 = N(t),  a polynomial of degree
 * <= 2 whose finite roots are one branch point (a root of R) and one NON-branch
 * point -- the latter forced to be F's pole.  The coefficient direction
 * (1, omega^2, omega) = (omega^(2 j)) is fixed: it is the unique vector killing
 * both the trace (V0) and the omega^1 (Y^1) components of the differential, so
 * only the overall scale C is free.  Reducing d/dt of the ansatz modulo
 * y^3 = R gives a Y-polynomial p0 + p1 Y + p2 Y^2 (rational in t); the two
 * conditions p0 == 0, p1 == 0 are the elementarity test, and matching the
 * remaining Y^2 slot to the integrand fixes  C = F * N / (R * p2)  -- which must
 * be free of t.  This is a Trager-Rothstein logarithmic-part computation
 * specialised to the cube-root-of-cubic curve; see the report in
 * goursat/goursat_cube_root_preprint_*.tex (third-kind remark) and
 * [[project_integrate_goursat]].
 *
 * Correct-by-verification: the eigenspace criterion does not certify this case,
 * so the answer is checked with the differentiate-back guard (diff_back_ok) in
 * gs_core.  Handles the m = 0 canonical linear form kappa t (F's pole a root of
 * N); a general affine form kappa t + m is deferred.  Borrows F, R, t.
 */
static Expr* goursat_cubic_thirdkind(Expr* F, Expr* R, Expr* t, int pnum) {
    if (pnum != 1) return NULL;                 /* p = 1/3 only for now */
    if (get_degree_poly(R, t) != 3) return NULL;

    Expr* lc = get_coeff(R, t, 3);              /* leading coefficient */
    if (!lc || is_zero(lc)) { if (lc) expr_free(lc); return NULL; }

    Expr* result = NULL;
    Expr* Y = fresh_var();
    Expr* w = mk_omega();
    Expr* kappa = mk_pow_rat(expr_copy(lc), 1, 3);          /* (lead R)^(1/3) */
    Expr* cubeR = mk_pow_rat(expr_copy(R), 1, 3);           /* R^(1/3)        */
    Expr* a = gmul(expr_copy(kappa), expr_copy(t));         /* kappa t        */
    Expr* Rp = eval_take(mk_fn2("D", expr_copy(R), expr_copy(t)));   /* R'   */
    Expr* Ypc = canonic(gdiv(expr_copy(Rp), gmul(mk_int(3), expr_copy(R))));  /* R'/(3R) */
    expr_free(Rp);

    /* omega^j a for j = 0, 1, 2. */
    Expr* wa[3] = {
        expr_copy(a),
        gmul(expr_copy(w), expr_copy(a)),
        gmul(mk_pow_int(expr_copy(w), 2), expr_copy(a))
    };
    /* omega^j kappa for j = 0, 1, 2 (the derivative of omega^j a is omega^j kappa). */
    Expr* wk[3] = {
        expr_copy(kappa),
        gmul(expr_copy(w), expr_copy(kappa)),
        gmul(mk_pow_int(expr_copy(w), 2), expr_copy(kappa))
    };
    /* Coefficient direction omega^(2 j). */
    Expr* w2[3] = {
        mk_int(1),
        mk_pow_int(expr_copy(w), 2),
        mk_pow_int(expr_copy(w), 4)
    };

    /* P(Y) = Sum_j omega^(2 j) (Ypc*Y - omega^j kappa) * Prod_{i != j}(Y - omega^i a). */
    Expr* P = mk_int(0);
    for (int j = 0; j < 3; j++) {
        int i1 = (j + 1) % 3, i2 = (j + 2) % 3;
        Expr* prod2 = gmul(gsub(expr_copy(Y), expr_copy(wa[i1])),
                           gsub(expr_copy(Y), expr_copy(wa[i2])));
        Expr* lin = gsub(gmul(expr_copy(Ypc), expr_copy(Y)), expr_copy(wk[j]));
        Expr* term = gmul(expr_copy(w2[j]), gmul(lin, prod2));
        P = gadd(P, term);
    }
    Expr* Pexp = expand_e(P);                   /* expand in Y (and t) */
    if (!Pexp) goto tk_cleanup;

    /* Reduce mod Y^3 = R:  p0' = p0 + p3 R,  keep p1, p2. */
    Expr* p0 = get_coeff(Pexp, Y, 0);
    Expr* p1 = get_coeff(Pexp, Y, 1);
    Expr* p2 = get_coeff(Pexp, Y, 2);
    Expr* p3 = get_coeff(Pexp, Y, 3);
    expr_free(Pexp);
    if (!p0) p0 = mk_int(0);
    if (!p1) p1 = mk_int(0);
    if (!p2) p2 = mk_int(0);
    if (!p3) p3 = mk_int(0);

    /* Reduce mod Y^3 = R: the Y^0 slot becomes p0 + p3 R.  canonic (Cancel over
     * Q(kappa, omega)) is needed before the zero test -- the raw coefficients
     * carry uncollapsed cyclotomic sums (1 + omega^2 - omega, ...) that
     * zero_test_decide cannot see through numerically. */
    Expr* p0r = canonic(gadd(p0, gmul(p3, expr_copy(R))));   /* consumes p0, p3 */
    Expr* p1r = canonic(p1);
    Expr* p2r = canonic(p2);

    /* Elementarity: the Y^0 and Y^1 slots of the reduced differential vanish. */
    bool ok = (p0r && is_zero(p0r) && p1r && is_zero(p1r)
               && p2r && !is_zero(p2r));
    if (ok) {
        /* N(t) = R - (lead R) t^3  (degree <= 2). */
        Expr* N = expand_e(gsub(expr_copy(R),
                        gmul(expr_copy(lc), mk_pow_int(expr_copy(t), 3))));
        Expr* Craw = gdiv(gmul(expr_copy(F), N),
                          gmul(expr_copy(R), expr_copy(p2r)));
        /* The overall scale C = F N / (R p2r) is a CONSTANT; obtain it by
         * substituting t -> t0
         * (which removes t) then canonic-ing the resulting parameter-only
         * algebraic number.  Sampling first is essential for a PARAMETRIC
         * radicand: canonic on the full t-bearing Craw over Q(param)(kappa, omega)
         * blows up, whereas the t-free sample is a cheap number-field Cancel.
         * diff_back_ok is the correctness gate, so a sample landing on a pole
         * (rejected here via ComplexInfinity / Indeterminate) or a spurious
         * value simply fails verification and declines. */
        Expr* C = NULL;
        {
            static const long sn[5] = { 1, 3, 2, 4, 5 };
            static const long sd[5] = { 7, 7, 9, 9, 11 };
            Expr* cinf = mk_sym("ComplexInfinity");
            Expr* indet = mk_sym("Indeterminate");
            for (int s = 0; s < 5 && !C; s++) {
                Expr* cs = canonic(subst_eval(expr_copy(Craw), t,
                                make_rational(sn[s], sd[s])));
                if (cs && expr_free_of(cs, cinf) && expr_free_of(cs, indet))
                    C = cs;
                else if (cs) expr_free(cs);
            }
            expr_free(cinf); expr_free(indet);
        }
        expr_free(Craw);
        if (C) {
            /* G = C * Sum_j omega^(2 j) Log[R^(1/3) - omega^j kappa t]. */
            Expr* G0 = mk_int(0);
            for (int j = 0; j < 3; j++) {
                Expr* arg = gsub(expr_copy(cubeR), expr_copy(wa[j]));
                Expr* lg = mk_fn1("Log", arg);
                G0 = gadd(G0, gmul(expr_copy(w2[j]), lg));
            }
            Expr* G = eval_take(gmul(C, G0));
            /* integrand f = F R^(-1/3) for the differentiate-back guard. */
            Expr* f = eval_take(gmul(expr_copy(F),
                            mk_pow_rat(expr_copy(R), -1, 3)));
            if (G && f && diff_back_ok(G, t, f)) result = expr_copy(G);
            if (G) expr_free(G);
            if (f) expr_free(f);
        }
    }
    if (p0r) expr_free(p0r);
    if (p1r) expr_free(p1r);
    if (p2r) expr_free(p2r);

tk_cleanup:
    for (int j = 0; j < 3; j++) { expr_free(wa[j]); expr_free(wk[j]); expr_free(w2[j]); }
    expr_free(Y); expr_free(w); expr_free(kappa); expr_free(cubeR);
    expr_free(a); expr_free(Ypc); expr_free(lc);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Period-3 case: sqrt(cubic) higher symmetry (Goursat 1887, Section 4)    */
/*                                                                         */
/* When R is a cubic linearly equivalent to t^3 - 1 there is an order-3    */
/* Mobius S fixing one of the four ramification points {roots, Infinity}   */
/* and cycling the other three.  If F is a non-trivial period-3 character  */
/* of <S> (F(S) = alpha F with alpha = Exp[2 Pi I/3]) then                 */
/* Int[F(t)/Sqrt[R(t)], t] is pseudo-elliptic.  This is the case the V4    */
/* Theorems 1-2 (goursat_v4) do NOT cover: the V4 trivial projection is    */
/* non-zero, but the period-3 trivial projection F + F(S) + F(S^2) = 0.    */
/*                                                                         */
/* Reduction (verified against the Section-4 worked example                */
/*   (1/3) dx/Sqrt[x(1-x)] = (t-1)/(t+2) dt/Sqrt[t^3-1],  x=((t-1)/(t+2))^3 */
/* ): with z = (t - a)/(t - a1) for the two fixed points a, a1 of S, the   */
/* map becomes z -> alpha z, the differential dt/Sqrt[R] becomes           */
/* dz/Sqrt[z(A z^3 + B)], and the alpha-character F1(z) = z phi(z^3).       */
/* Setting x = z^3 yields the elementary (1/3) phi(x)/Sqrt[x(A x + B)] dx.  */
/* ---------------------------------------------------------------------- */

/* True iff g(z) is (numerically) invariant under z -> omega z, i.e. g is a
 * function of z^3 alone -- the structural pre-condition for to_function_of_power
 * to extract phi(x) faithfully (it silently drops non-multiple-of-3 powers).
 * Used to detect the wrong fixed-point orientation (which turns an alpha-
 * character into an alpha^2-character, making Fz/z carry a stray z). Borrows. */
static bool is_cube_function(Expr* g, Expr* z) {
    Expr* w  = mk_omega();
    Expr* gw = subst_eval(expr_copy(g), z, gmul(w, expr_copy(z)));  /* g(omega z) */
    Expr* diff = gsub(expr_copy(g), gw);
    bool nonconst = sample_clearly_nonzero(diff, z);
    expr_free(diff);
    return !nonconst;
}

/* Reduce a pure alpha-character piece Fk (Fk(S)=alpha Fk; S the order-3 Mobius
 * with fixed points {a, a1}, z=(t-a)/(t-a1) sending S to z->alpha z) to its
 * elementary antiderivative in t, or NULL if Fk is not of this form for this
 * fixed-point orientation / the reduced rational integral does not close.
 * a or a1 may be NULL (Infinity).  Borrows Fk, R, t, a, a1. */
static Expr* period3_reduce(Expr* Fk, Expr* R, Expr* t, Expr* a, Expr* a1) {
    bool ainf = (a == NULL), a1inf = (a1 == NULL);
    if (ainf && a1inf) return NULL;

    Expr* z = fresh_var();
    Expr* x = fresh_var();
    Expr* result = NULL;
    Expr* tz = NULL, *Fz = NULL, *phi = NULL, *Qz = NULL, *AB = NULL, *G = NULL;
    Expr* wscale = NULL;   /* leading-coeff rescale: integration variable w = wscale * z^3 */

    /* t(z): z=(t-a)/(t-a1) sends a->0, a1->Infinity. */
    if (a1inf)     tz = gadd(expr_copy(a), expr_copy(z));              /* z = t - a    */
    else if (ainf) tz = gadd(expr_copy(a1), mk_inv(expr_copy(z)));     /* z = 1/(t-a1) */
    else           tz = canonic(gdiv(gsub(expr_copy(a),
                              gmul(expr_copy(a1), expr_copy(z))),
                              gsub(mk_int(1), expr_copy(z))));
    if (!tz) goto done;

    /* Fz = Fk(t(z)) must be z * phi(z^3).  Verify Fz/z is a function of z^3
     * (else this orientation gives the alpha^2-character -- decline it).  Do the
     * check on the RAW (un-canonicalised) substitution: it is a fast numeric
     * sample, and skipping the costly cyclotomic canonic on the common
     * non-matching orientation/integrand is what keeps decline fast. */
    Fz = subst_eval(expr_copy(Fk), t, expr_copy(tz));
    if (!Fz) goto done;
    Expr* FzOverZraw = gdiv(expr_copy(Fz), expr_copy(z));
    bool cube_ok = is_cube_function(FzOverZraw, z);
    if (!cube_ok) { expr_free(FzOverZraw); goto done; }
    Expr* FzOverZ = canonic(FzOverZraw);   /* matched orientation only */
    phi = FzOverZ ? to_function_of_power(FzOverZ, z, x, 3) : NULL;
    if (FzOverZ) expr_free(FzOverZ);
    if (!phi) goto done;

    /* Qz = R(t(z)) (dz/dt)^2 must be z (A z^3 + B) for the finite-fixed-point
     * reduction, i.e. Qz/z is a function of z^3.  (For a fixed point that is not
     * a ramification point -- e.g. the t^3=x / fix-Infinity reduction -- Qz has
     * a different shape that this branch does not handle; reject it fast on the
     * raw form rather than tripping a div-by-zero in to_function_of_power.) */
    {
        Expr* dtz = eval_take(mk_fn2("D", expr_copy(tz), expr_copy(z)));
        Expr* Rtz = subst_eval(expr_copy(R), t, expr_copy(tz));
        Qz = gdiv(Rtz, mk_pow_int(dtz, 2));         /* raw, no canonic yet */
    }
    if (!Qz) goto done;
    Expr* QzOverZraw = gdiv(expr_copy(Qz), expr_copy(z));
    if (!is_cube_function(QzOverZraw, z)) { expr_free(QzOverZraw); goto done; }
    Expr* QzOverZ = canonic(QzOverZraw);
    AB = QzOverZ ? to_function_of_power(QzOverZ, z, x, 3) : NULL;
    if (QzOverZ) expr_free(QzOverZ);
    if (!AB) goto done;

    /* integrand_x = (1/3) phi(x) / Sqrt[x (A x + B)].
     *
     * rad = x(A x + B) = A x^2 + B x has a (generally cyclotomic) LEADING
     * coefficient lc = A.  Mathilda's Sqrt-of-quadratic integrator handles an
     * algebraic-number constant on the linear/constant term but NOT on the
     * leading term (the Euler substitution needs Sqrt[lc] and chokes on an
     * algebraic lc).  Rescale x -> w/lc so the radicand becomes the MONIC
     * w^2 + B w (B = original linear coeff, kept intact on the linear term --
     * NOT divided by lc, which would require the cyclotomic constant-ratio
     * B/lc that Cancel/Simplify cannot reduce, e.g. (-1+zeta)/(1-zeta) -> -1):
     *      (1/3) phi(x)/Sqrt[rad]  dx
     *    = (1/(3 Sqrt[lc])) phi(w/lc)/Sqrt[w^2 + B w]  dw.
     * The 1/Sqrt[lc] is a CONSTANT (branch-safe -- never straddles the
     * w-dependent cut).  Integrate in w, then back-substitute w -> lc z(t)^3
     * (since x = w/lc and x = z(t)^3 give w = lc z^3). */
    {
        Expr* rad = expand_e(canonic(gmul(expr_copy(x), expr_copy(AB))));
        if (!rad) goto done;
        Expr* lc = get_coeff(rad, x, 2);                 /* leading coeff A */
        Expr* B  = get_coeff(rad, x, 1);                 /* linear coeff (kept) */
        expr_free(rad);
        if (!lc || is_zero(lc)) { if (lc) expr_free(lc); if (B) expr_free(B); goto done; }
        if (!B) B = mk_int(0);
        wscale = expr_copy(lc);                          /* w = lc * z^3 */
        Expr* monicrad = canonic(gadd(mk_pow_int(expr_copy(x), 2),
                                      gmul(B, expr_copy(x))));   /* w^2 + B w, consumes B */
        Expr* phi_sub = subst_eval(expr_copy(phi), x,
                                   gdiv(expr_copy(x), expr_copy(lc)));   /* phi(w/lc) */
        if (!monicrad || !phi_sub) {
            if (monicrad) expr_free(monicrad); if (phi_sub) expr_free(phi_sub);
            expr_free(lc); goto done;
        }
        /* Pull the constant prefactor 1/(3 Sqrt[lc]) OUT of the integral: a
         * nested-radical constant Sqrt[lc] left inside the integrand is mistaken
         * for a second variable radical by the radical integrator and hangs it.
         * Integrate the bare phi(w/lc)/Sqrt[w^2+Bw], then multiply the constant. */
        Expr* integrand = gmul(phi_sub, mk_pow_rat(monicrad, -1, 2));  /* consumes phi_sub, monicrad */
        gs_depth++;
        G = integrate_in(integrand, x);
        gs_depth--;
        expr_free(integrand);
        if (G) G = canonic(gmul(
            gmul(mk_rat(1, 3), mk_inv(mk_sqrt_expr(expr_copy(lc)))), G));
        expr_free(lc);
    }
    if (!G) goto done;

    /* back-substitute w -> wscale * z(t)^3  (w = lc x, x = z(t)^3). */
    {
        Expr* back;
        if (a1inf)     back = mk_pow_int(gsub(expr_copy(t), expr_copy(a)), 3);
        else if (ainf) back = mk_pow_int(mk_inv(gsub(expr_copy(t), expr_copy(a1))), 3);
        else           back = mk_pow_int(gdiv(gsub(expr_copy(t), expr_copy(a)),
                                              gsub(expr_copy(t), expr_copy(a1))), 3);
        if (wscale) back = gmul(expr_copy(wscale), back);   /* w = wscale * z^3 */
        result = subst_eval(G, x, back);   /* consumes G, back */
        G = NULL;
    }

done:
    if (wscale) expr_free(wscale);
    if (tz) expr_free(tz);
    if (Fz) expr_free(Fz);
    if (phi) expr_free(phi);
    if (Qz) expr_free(Qz);
    if (AB) expr_free(AB);
    if (G) expr_free(G);
    expr_free(z); expr_free(x);
    return result;
}

/* Reduce a present character piece by trying both fixed-point orientations
 * (one is the alpha-character orientation, the other alpha^2).  Returns the
 * antiderivative or NULL.  Borrows all. */
static Expr* period3_reduce_either(Expr* Fk, Expr* R, Expr* t, Expr* a, Expr* a1) {
    Expr* r = period3_reduce(Fk, R, t, a, a1);
    if (!r) r = period3_reduce(Fk, R, t, a1, a);
    return r;
}

/* Period-3 Goursat (p = 1/2, R a cubic with the t^3-1 higher symmetry). */
static Expr* goursat_period3(Expr* F, Expr* R, Expr* t) {
    size_t nr = 0;
    Expr** roots = poly_roots(R, t, &nr);
    int dR = get_degree_poly(R, t);
    if (dR != 3 || nr != 3) {
        for (size_t i = 0; i < nr; i++) expr_free(roots[i]);
        if (roots) free(roots);
        return NULL;
    }

    Expr* result = NULL;

    /* Four ramification points; try each as the fixed point of an order-3
     * Mobius cycling the other three. */
    Expr* bp[4] = { roots[0], roots[1], roots[2], NULL /* Infinity */ };
    for (int fix = 0; fix < 4 && !result; fix++) {
        Expr* cyc[3]; int ci = 0;
        for (int i = 0; i < 4; i++) if (i != fix) cyc[ci++] = bp[i];
        Expr* S = cyclic_mobius(cyc, t);
        if (!S) continue;

        /* Period-3 trivial projection F + F(S) + F(S^2) must vanish (numeric
         * gate on the *uncombined* sum -- a symbolic character projection in t
         * would trigger the cyclotomic Together blowup).  Compose F(S^2) as
         * (F(S))(S) rather than forming the cyclotomic S^2 = canonic(S.S). */
        gs_log("  period-3 (fixed point #%d): trivial projection F + F(S) + F(S^2) must vanish", fix);
        Expr* FS  = subst_eval(expr_copy(F), t, expr_copy(S));
        Expr* FS2 = FS ? subst_eval(expr_copy(FS), t, expr_copy(S)) : NULL;
        Expr* F0raw = mk_fn3("Plus", expr_copy(F), FS, FS2);
        bool trivnz = sample_clearly_nonzero(F0raw, t);
        expr_free(F0raw);
        if (trivnz) {
            gs_log("  period-3: trivial projection non-zero for this fixed point -- trying next");
            expr_free(S); continue;
        }

        Expr* a = NULL; Expr* a1 = NULL;
        /* Clean fixed points.  For the centered pure cubic R = lc((t-s)^3 + c)
         * (the t^3-1 family) the order-3 Mobius fixing a finite root r has fixed
         * points {r, 3s-2r} exactly (Goursat's ((t-r)/(t+2r))^3 substitutions,
         * here s = sum-of-roots/3).  Supplying these directly -- instead of
         * Solve[S==t] + Simplify, which for the cyclotomic roots returns nested
         * radicals Sqrt[-9(-1)^(1/3)] that leave roots-of-unity uncollapsed and
         * make to_function_of_power divide by zero (ComplexInfinity) -- keeps the
         * reduction inside the pure cyclotomic field Q(zeta_3), where canonic is
         * fast and complete.  s computed from the t^2 coefficient. */
        if (fix < 3) {
            Expr* c3 = get_coeff(R, t, 3);
            Expr* c2 = get_coeff(R, t, 2);
            Expr* s  = canonic(gdiv(mk_neg(c2 ? c2 : mk_int(0)),
                                    gmul(mk_int(3), c3 ? c3 : mk_int(1))));
            a  = expr_copy(bp[fix]);
            a1 = canonic(gsub(gmul(mk_int(3), expr_copy(s)),
                              gmul(mk_int(2), expr_copy(bp[fix]))));   /* 3s - 2r */
            expr_free(s);
            expr_free(S);
        } else {
            bool okfp = fixed_points(S, t, &a, &a1);
            expr_free(S);
            if (!okfp) { if (a) expr_free(a); if (a1) expr_free(a1); continue; }
            a  = simplify_e(a);
            a1 = simplify_e(a1);
        }

        /* Reduce F directly: substituting the (simple) ORIGINAL F into t(z)
         * collapses it BEFORE any cyclotomic projection (e.g. F(t(z)) = z for
         * (t-1)/(t+2)), and period3_reduce's is_cube_function check verifies it
         * is a pure alpha-character for this orientation -- trying both
         * fixed-point orders covers alpha vs alpha^2. */
        result = period3_reduce_either(F, R, t, a, a1);

        if (a) expr_free(a); if (a1) expr_free(a1);
    }

    for (size_t i = 0; i < nr; i++) expr_free(roots[i]);
    if (roots) free(roots);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Fourth-root case: order-4 Mobius eigendescent (p = 1/4 and 3/4)        */
/* ---------------------------------------------------------------------- */

/* Cross-ratio (r1, r3; r2, r4) = ((r1-r2)(r3-r4)) / ((r1-r4)(r3-r2)).
 * Borrows; returns owned (un-canonicalised). */
static Expr* cross_ratio(Expr* r1, Expr* r2, Expr* r3, Expr* r4) {
    Expr* num = gmul(gsub(expr_copy(r1), expr_copy(r2)),
                     gsub(expr_copy(r3), expr_copy(r4)));
    Expr* den = gmul(gsub(expr_copy(r1), expr_copy(r4)),
                     gsub(expr_copy(r3), expr_copy(r2)));
    return gdiv(num, den);
}

/* Find a harmonic cyclic ordering of four roots: an arrangement
 * (r1,r2,r3,r4) with (r1,r3;r2,r4) = -1, which is exactly the condition for an
 * order-4 Mobius cycling r1->r2->r3->r4 to exist (fourth-root paper Lemma 2.1).
 * Fixes r1 = roots[0] and tries the 6 orderings of the rest.  On success fills
 * `out` with borrowed pointers (into roots) and returns true. */
static bool harmonic_order(Expr** roots, Expr** out) {
    static const int perm[6][3] = {
        {1,2,3},{1,3,2},{2,1,3},{2,3,1},{3,1,2},{3,2,1}
    };
    for (int p = 0; p < 6; p++) {
        Expr* r1 = roots[0];
        Expr* r2 = roots[perm[p][0]];
        Expr* r3 = roots[perm[p][1]];
        Expr* r4 = roots[perm[p][2]];
        Expr* cr = cross_ratio(r1, r2, r3, r4);
        Expr* test = gadd(cr, mk_int(1));       /* cr + 1 == 0 ? */
        bool harm = is_zero(test);
        expr_free(test);
        if (harm) { out[0]=r1; out[1]=r2; out[2]=r3; out[3]=r4; return true; }
    }
    return false;
}

/* Order-4 Mobius cycling r[0]->r[1]->r[2]->r[3].  Solves
 * {(A r_i + B)/(C r_i + 1) == r_{i+1} : i = 0,1,2} for {A,B,C}.  Borrows
 * r[*], t; returns owned canonic S(t) or NULL.  (Appendix QuarticMobius.) */
static Expr* quartic_mobius(Expr** r, Expr* t) {
    Expr* A = fresh_var(); Expr* B = fresh_var(); Expr* C = fresh_var();
    Expr* eqs[3];
    for (int i = 0; i < 3; i++) {
        Expr* lhs = gadd(gmul(expr_copy(A), expr_copy(r[i])), expr_copy(B));
        Expr* rhs = gmul(expr_copy(r[i + 1]),
                         gadd(gmul(expr_copy(C), expr_copy(r[i])), mk_int(1)));
        eqs[i] = mk_fn2("Equal", lhs, rhs);
    }
    Expr* eqlist = mk_fnv("List", eqs, 3);
    Expr* vars = mk_fn3("List", expr_copy(A), expr_copy(B), expr_copy(C));
    Expr* solset = eval_take(mk_fn2("Solve", eqlist, vars));

    Expr* S = NULL;
    if (solset && head_is(solset, SYM_List) && solset->data.function.arg_count >= 1) {
        Expr* sol = solset->data.function.args[0];
        Expr* Sexpr = gdiv(gadd(gmul(expr_copy(A), expr_copy(t)), expr_copy(B)),
                           gadd(gmul(expr_copy(C), expr_copy(t)), mk_int(1)));
        S = canonic(eval_take(internal_replace_all(
            (Expr*[]){ Sexpr, expr_copy(sol) }, 2)));
    }
    if (solset) expr_free(solset);
    expr_free(A); expr_free(B); expr_free(C);
    return S;
}

/* P_k projection under z -> i z: (H + i^-k H(iz) + i^-2k H(-z) + i^-3k H(-iz))/4.
 * Borrows H, z; returns owned Simplify+canonic.  (Appendix EigenProjection4.) */
static Expr* eigenproj4(Expr* H, Expr* z, int k) {
    Expr* I = mk_sym("I");
    Expr* Hiz  = subst_eval(expr_copy(H), z, gmul(expr_copy(I), expr_copy(z)));
    Expr* Hmz  = subst_eval(expr_copy(H), z, mk_neg(expr_copy(z)));
    Expr* Hmiz = subst_eval(expr_copy(H), z,
                            mk_neg(gmul(expr_copy(I), expr_copy(z))));
    Expr* t0 = expr_copy(H);
    Expr* t1 = gmul(mk_pow_int(expr_copy(I), -k), Hiz);
    Expr* t2 = gmul(mk_pow_int(expr_copy(I), -2 * k), Hmz);
    Expr* t3 = gmul(mk_pow_int(expr_copy(I), -3 * k), Hmiz);
    expr_free(I);
    Expr* sum = expr_new_function(mk_sym("Plus"),
                    (Expr*[]){ t0, t1, t2, t3 }, 4);
    return canonic(simplify_e(gdiv(sum, mk_int(4))));
}

/* Raw (evaluate-only, NO canonic/Simplify) order-4 eigenprojection for the
 * numeric criterion gate.  Borrows H,z.  See eigenproj3_raw. */
static Expr* eigenproj4_raw(Expr* H, Expr* z, int k) {
    Expr* I = mk_sym("I");
    Expr* Hiz  = subst_eval(expr_copy(H), z, gmul(expr_copy(I), expr_copy(z)));
    Expr* Hmz  = subst_eval(expr_copy(H), z, mk_neg(expr_copy(z)));
    Expr* Hmiz = subst_eval(expr_copy(H), z,
                            mk_neg(gmul(expr_copy(I), expr_copy(z))));
    Expr* t0 = expr_copy(H);
    Expr* t1 = gmul(mk_pow_int(expr_copy(I), -k), Hiz);
    Expr* t2 = gmul(mk_pow_int(expr_copy(I), -2 * k), Hmz);
    Expr* t3 = gmul(mk_pow_int(expr_copy(I), -3 * k), Hmiz);
    expr_free(I);
    Expr* sum = expr_new_function(mk_sym("Plus"),
                    (Expr*[]){ t0, t1, t2, t3 }, 4);
    return eval_take(gdiv(sum, mk_int(4)));
}

/* Fourth-root Goursat (p = pnum/4, pnum in {1,3}): R quartic, harmonic roots. */
static Expr* goursat_quartic(Expr* F, Expr* R, Expr* t, int pnum) {
    size_t nr = 0;
    Expr** roots = poly_roots(R, t, &nr);
    int dR = get_degree_poly(R, t);
    if (dR != 4 || nr != 4) {
        for (size_t i = 0; i < nr; i++) expr_free(roots[i]);
        if (roots) free(roots);
        return NULL;
    }

    Expr* result = NULL;
    Expr* z = NULL; Expr* x = NULL; Expr* uOut = NULL;
    Expr* alpha = NULL; Expr* beta = NULL;
    Expr* tz = NULL; Expr* H = NULL;
    Expr* H0 = NULL; Expr* H1 = NULL; Expr* H2 = NULL; Expr* H3 = NULL;
    Expr* cval = NULL; Expr* K = NULL; Expr* S = NULL;

    /* Roots must be harmonic; otherwise this theory does not apply. */
    Expr* ord[4];
    if (!harmonic_order(roots, ord)) {
        gs_log("  quartic: the four roots of R are NOT harmonic (cross-ratio != -1) -- declining");
        goto cleanup;
    }
    gs_log("  quartic: roots are harmonic -- building order-4 Mobius automorphism");

    S = quartic_mobius(ord, t);
    if (!S) goto cleanup;
    if (!fixed_points(S, t, &alpha, &beta)) goto cleanup;
    bool binf = (beta == NULL);

    z = fresh_var(); x = fresh_var(); uOut = fresh_var();

    tz = binf
        ? gadd(expr_copy(alpha), expr_copy(z))
        : canonic(gdiv(gsub(expr_copy(alpha), gmul(expr_copy(beta), expr_copy(z))),
                       gsub(mk_int(1), expr_copy(z))));

    /* R(t(z)) [ (1-z)^4 ]. */
    Expr* Rsub = subst_eval(expr_copy(R), t, expr_copy(tz));
    Expr* Rz = binf ? canonic(Rsub)
                    : canonic(gmul(Rsub,
                          mk_pow_int(gsub(mk_int(1), expr_copy(z)), 4)));
    Expr* Rze = expand_e(Rz);
    if (!Rze) goto cleanup;
    cval = get_coeff(Rze, z, 4);
    Expr* c0 = get_coeff(Rze, z, 0);
    expr_free(Rze);
    if (!cval || !c0 || is_zero(cval)) {
        if (cval) { expr_free(cval); cval = NULL; }
        if (c0) expr_free(c0);
        goto cleanup;
    }
    K = canonic(gdiv(mk_neg(c0), expr_copy(cval)));

    /* H_p = (alpha-beta)(1-z)^(4p-2)(F/.t->tz)/cval^p  [binf: drop the
     * (alpha-beta) and (1-z) factors].  4p-2 = -1 (p=1/4) or +1 (p=3/4).
     * Build the raw source first; canonic is deferred until after the gate. */
    Expr* base = NULL;
    {
        Expr* Ftz = subst_eval(expr_copy(F), t, expr_copy(tz));
        Expr* cpow = (pnum == 1) ? mk_pow_rat(expr_copy(cval), -1, 4)
                                 : mk_pow_rat(expr_copy(cval), -3, 4);
        if (binf) {
            base = gmul(Ftz, cpow);
        } else {
            int e = (pnum == 1) ? -1 : 1;   /* (1-z)^(4p-2) */
            Expr* omz = mk_pow_int(gsub(mk_int(1), expr_copy(z)), e);
            base = gmul(gmul(gsub(expr_copy(alpha), expr_copy(beta)), Ftz),
                        gmul(cpow, omz));
        }
    }
    if (!base) goto cleanup;

    /* Numeric elementarity gate (p=1/4 needs H1==H2==0; p=3/4 needs H0==H1==0):
     * sample the two RAW criterion eigenprojections and decline at once if
     * either is decisively non-zero, sidestepping the cyclotomic-field
     * canonic/Simplify blowup on the obstructed cases.  A passing sample falls
     * through to the identical symbolic path; gs_core's diff_back_ok verifies
     * any antiderivative. */
    {
        int ka = (pnum == 1) ? 1 : 0;
        int kb = (pnum == 1) ? 2 : 1;
        gs_log("  quartic criterion: i-eigencomponents H%d and H%d of H must both vanish",
               ka, kb);
        Expr* ca = eigenproj4_raw(base, z, ka);
        bool obstructed = sample_clearly_nonzero(ca, z);
        if (ca) expr_free(ca);
        if (!obstructed) {
            Expr* cb = eigenproj4_raw(base, z, kb);
            obstructed = sample_clearly_nonzero(cb, z);
            if (cb) expr_free(cb);
        }
        if (obstructed) {
            gs_log("  numeric gate: an eigencomponent is decisively non-zero -- OBSTRUCTED, declining");
            expr_free(base); goto cleanup;
        }
    }
    gs_log("  numeric gate: criterion eigencomponents look vanishing -- proceeding to exact descent");

    H = canonic(base);
    if (!H) goto cleanup;

    H0 = eigenproj4(H, z, 0);
    H1 = eigenproj4(H, z, 1);
    H2 = eigenproj4(H, z, 2);
    H3 = eigenproj4(H, z, 3);

    /* Criterion: p=1/4 needs H1==0 and H2==0; p=3/4 needs H0==0 and H1==0. */
    if (pnum == 1) {
        if (!H1 || !H2 || !is_zero(H1) || !is_zero(H2)) {
            gs_log("  exact criterion: H1==0 && H2==0 FAILS -- declining");
            goto cleanup;
        }
    } else {
        if (!H0 || !H1 || !is_zero(H0) || !is_zero(H1)) {
            gs_log("  exact criterion: H0==0 && H1==0 FAILS -- declining");
            goto cleanup;
        }
    }
    gs_log("  exact criterion HOLDS -- descending eigenpieces to genus-0");

    Expr* q4_c = mk_pow_rat(expr_copy(cval), 1, 4);    /* cval^(1/4) */
    Expr* q4_R = mk_pow_rat(expr_copy(R), 1, 4);        /* R^(1/4)    */
    Expr* amb  = binf ? NULL : gsub(expr_copy(alpha), expr_copy(beta));

    Expr* total = NULL;
    if (pnum == 1) {
        /* phi0 = ToFunctionOfFourth[H0]; phi3 = ToFunctionOfFourth[H3/z^3]. */
        Expr* phi0 = to_function_of_power(H0, z, x, 4);
        Expr* H3z3 = canonic(gdiv(expr_copy(H3), mk_pow_int(expr_copy(z), 3)));
        Expr* phi3 = H3z3 ? to_function_of_power(H3z3, z, x, 4) : NULL;
        if (H3z3) expr_free(H3z3);

        /* J0 = (phi0 /. x->K/(1-w^4)) w^2/(1-w^4); J3 = (phi3 /. x->u^4+K) u^2. */
        Expr* om4 = gsub(mk_int(1), mk_pow_int(expr_copy(uOut), 4));
        Expr* J0 = phi0 ? canonic(gmul(
            subst_eval(phi0, x, gdiv(expr_copy(K), expr_copy(om4))),
            gdiv(mk_pow_int(expr_copy(uOut), 2), expr_copy(om4)))) : NULL;
        expr_free(om4);
        Expr* J3 = NULL;
        if (phi3) {
            Expr* u4K = gadd(mk_pow_int(expr_copy(uOut), 4), expr_copy(K));
            J3 = canonic(gmul(subst_eval(phi3, x, u4K),
                              mk_pow_int(expr_copy(uOut), 2)));
        }

        Expr* J0back = binf
            ? gdiv(expr_copy(q4_R), gmul(expr_copy(q4_c),
                   gsub(expr_copy(t), expr_copy(alpha))))
            : gdiv(gmul(expr_copy(q4_R), expr_copy(amb)),
                   gmul(expr_copy(q4_c), gsub(expr_copy(t), expr_copy(alpha))));
        Expr* J3back = binf
            ? gdiv(expr_copy(q4_R), expr_copy(q4_c))
            : gdiv(gmul(expr_copy(q4_R), expr_copy(amb)),
                   gmul(expr_copy(q4_c), gsub(expr_copy(t), expr_copy(beta))));

        Expr* p0 = integ_backsub(J0, uOut, J0back);
        Expr* p3 = (p0 != NULL) ? integ_backsub(J3, uOut, J3back) : NULL;
        if (J0) expr_free(J0);
        if (J3) expr_free(J3);
        expr_free(J0back); expr_free(J3back);
        if (p0 && p3) total = gadd(p0, p3);
        else { if (p0) expr_free(p0); if (p3) expr_free(p3); }
    } else {
        /* phi2 = ToFunctionOfFourth[H2/z^2]; phi3 = ToFunctionOfFourth[H3/z^3]. */
        Expr* H2z2 = canonic(gdiv(expr_copy(H2), mk_pow_int(expr_copy(z), 2)));
        Expr* phi2 = H2z2 ? to_function_of_power(H2z2, z, x, 4) : NULL;
        if (H2z2) expr_free(H2z2);
        Expr* H3z3 = canonic(gdiv(expr_copy(H3), mk_pow_int(expr_copy(z), 3)));
        Expr* phi3 = H3z3 ? to_function_of_power(H3z3, z, x, 4) : NULL;
        if (H3z3) expr_free(H3z3);

        /* J2 = -(phi2 /. x->K s^4/(s^4-1)) s^2/(s^4-1); J3 = phi3 /. x->v^4+K. */
        Expr* s4m1 = gsub(mk_pow_int(expr_copy(uOut), 4), mk_int(1));
        Expr* J2 = phi2 ? canonic(mk_neg(gmul(
            subst_eval(phi2, x, gdiv(gmul(expr_copy(K), mk_pow_int(expr_copy(uOut), 4)),
                                     expr_copy(s4m1))),
            gdiv(mk_pow_int(expr_copy(uOut), 2), expr_copy(s4m1))))) : NULL;
        expr_free(s4m1);
        Expr* J3 = NULL;
        if (phi3) {
            Expr* v4K = gadd(mk_pow_int(expr_copy(uOut), 4), expr_copy(K));
            J3 = canonic(subst_eval(phi3, x, v4K));
        }

        Expr* J2back = binf
            ? gdiv(gmul(expr_copy(q4_c), gsub(expr_copy(t), expr_copy(alpha))),
                   expr_copy(q4_R))
            : gdiv(gmul(expr_copy(q4_c), gsub(expr_copy(t), expr_copy(alpha))),
                   gmul(expr_copy(q4_R), expr_copy(amb)));
        Expr* J3back = binf
            ? gdiv(expr_copy(q4_R), expr_copy(q4_c))
            : gdiv(gmul(expr_copy(q4_R), expr_copy(amb)),
                   gmul(expr_copy(q4_c), gsub(expr_copy(t), expr_copy(beta))));

        Expr* p2 = integ_backsub(J2, uOut, J2back);
        Expr* p3 = (p2 != NULL) ? integ_backsub(J3, uOut, J3back) : NULL;
        if (J2) expr_free(J2);
        if (J3) expr_free(J3);
        expr_free(J2back); expr_free(J3back);
        if (p2 && p3) total = gadd(p2, p3);
        else { if (p2) expr_free(p2); if (p3) expr_free(p3); }
    }

    expr_free(q4_c); expr_free(q4_R); if (amb) expr_free(amb);
    if (total) result = eval_take(total);

cleanup:
    if (S) expr_free(S);
    if (z) expr_free(z); if (x) expr_free(x); if (uOut) expr_free(uOut);
    if (alpha) expr_free(alpha); if (beta) expr_free(beta);
    if (tz) expr_free(tz); if (H) expr_free(H);
    if (H0) expr_free(H0); if (H1) expr_free(H1);
    if (H2) expr_free(H2); if (H3) expr_free(H3);
    if (cval) expr_free(cval); if (K) expr_free(K);
    for (size_t i = 0; i < nr; i++) expr_free(roots[i]);
    if (roots) free(roots);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Recognition: f = F(t) * R(t)^q,  q in {-1/2,-1/3,-2/3,-1/4,-3/4}       */
/* ---------------------------------------------------------------------- */

/* On success returns true with owned *F (rational cofactor) and *R (polynomial
 * radicand), and the radical order p = pnum/pden (one of 1/2,1/3,2/3,1/4,3/4).
 * On false nothing is allocated.  Borrows f, x. */
static bool recognise(Expr* f, Expr* x, Expr** F, Expr** R,
                      int* pnum, int* pden) {
    Expr** factors; size_t nfac;
    if (head_is(f, SYM_Times)) {
        factors = f->data.function.args;
        nfac    = f->data.function.arg_count;
    } else {
        factors = &f;
        nfac    = 1;
    }

    Expr* cof = mk_int(1);
    Expr* rad = NULL;
    int pn = 0, pd = 0;

    for (size_t i = 0; i < nfac; i++) {
        Expr* g = factors[i];
        if (head_is(g, SYM_Power) && g->data.function.arg_count == 2 &&
            !expr_free_of(g->data.function.args[0], x)) {
            Expr* base = g->data.function.args[0];
            int64_t qn, qd;
            if (is_rational(g->data.function.args[1], &qn, &qd) && qd != 0) {
                /* Reduce defensively (Rational nodes are already normalised, but
                 * be safe) and force a positive denominator. */
                if (qd < 0) { qn = -qn; qd = -qd; }
                int64_t gg = gcd(qn < 0 ? -qn : qn, qd);
                if (gg > 1) { qn /= gg; qd /= gg; }
                /* Any rational exponent n/d with d in {2,3,4} is a radical.  Write
                 * R^(n/d) = R^k * R^(-p) where p = pnum/d in {1/2,1/3,2/3,1/4,3/4}
                 * (pnum = (-n) mod d in (0,d)) and the integer k = (n+pnum)/d folds
                 * into the rational cofactor F as R^k.  This admits POSITIVE-power
                 * radicals such as (1-x^3)^(1/3)/x = (1-x^3)/x * (1-x^3)^(-2/3),
                 * not just radicals already in the denominator. */
                if (qd == 2 || qd == 3 || qd == 4) {
                    if (rad) { goto fail; }      /* only one radical allowed */
                    int64_t pnum_l = ((-qn) % qd + qd) % qd;   /* in (0, qd) */
                    int64_t k = (qn + pnum_l) / qd;            /* integer */
                    rad = expr_copy(base);
                    pn = (int)pnum_l; pd = (int)qd;
                    if (k != 0)
                        cof = mk_fn2("Times", cof, mk_pow_int(expr_copy(base), k));
                    continue;
                }
            }
        }
        /* Everything else folds into the cofactor (validated rational below). */
        cof = mk_fn2("Times", cof, expr_copy(g));
    }

    if (!rad) goto fail;

    /* Radicand must be a polynomial in x. */
    Expr* Rexp = expand_e(expr_copy(rad));
    Expr* vars[1] = { x };
    if (!Rexp || !is_polynomial(Rexp, vars, 1)) { if (Rexp) expr_free(Rexp); goto fail; }

    /* Cofactor must be a rational function of x. */
    Expr* cc = canonic(expr_copy(cof));
    Expr* cnum = cc ? expand_e(numer_of(cc)) : NULL;
    Expr* cden = cc ? expand_e(denom_of(cc)) : NULL;
    bool crat = cnum && cden && is_polynomial(cnum, vars, 1)
                            && is_polynomial(cden, vars, 1);
    if (cnum) expr_free(cnum);
    if (cden) expr_free(cden);
    if (!crat) { if (cc) expr_free(cc); expr_free(Rexp); goto fail; }

    *F = cc ? cc : cof;
    if (cc) expr_free(cof);
    *R = Rexp;
    expr_free(rad);
    *pnum = pn; *pden = pd;
    return true;

fail:
    expr_free(cof);
    if (rad) expr_free(rad);
    return false;
}

/* ---------------------------------------------------------------------- */
/* Output normalisation: re-express nested radicals over Sqrt[R]           */
/* ---------------------------------------------------------------------- */
/*
 * The square-root descent integrates in a Mobius-substituted coordinate, so its
 * antiderivative carries a NESTED radical Sqrt[g] -- g a rational function of x
 * that, by construction, equals R times a perfect rational square (R is the
 * integrand's own radicand).  Sqrt[g] is therefore (rational)*Sqrt[R], but the
 * generic Cancel/Together that builds the result never extracts the square
 * factor (doing so is branch-cut-sensitive in general, so Simplify rightly
 * refuses it without assumptions).  The upshot: the returned antiderivative
 * carries a DIFFERENT radical from the integrand, and D[result,x] // Simplify
 * cannot reduce to the integrand (it would have to relate two distinct nested
 * radicals).  Re-express Sqrt[g] -> (rational)*Sqrt[R] here, where the integrator
 * KNOWS g/R is a perfect square and fixes the sole branch (a global +/- sign)
 * numerically -- a correct-by-construction normalisation of the *result*, not a
 * change to Simplify.  The rewrite is adopted only if the differentiate-back
 * guard still passes on the rewritten form.
 */

/* True iff e contains an unresolved Power[_, 1/2] (a square root). */
static bool has_sqrt(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is((Expr*)e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* half = make_rational(1, 2);
        bool m = expr_eq(e->data.function.args[1], half);
        expr_free(half);
        if (m) return true;
    }
    if (has_sqrt(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_sqrt(e->data.function.args[i])) return true;
    return false;
}

/* If exponent e is a rational with denominator 2 (numerator necessarily odd),
 * report 2*e (the odd numerator, carrying the sign) in *two_e and return true;
 * else false.  Recognises the half-integer powers Power[base, m/2] that carry a
 * single square-root generator.  Borrows e. */
static bool is_half_power(const Expr* e, int64_t* two_e) {
    if (!e || !head_is((Expr*)e, SYM_Rational)
        || e->data.function.arg_count != 2) return false;
    Expr* p = e->data.function.args[0];
    Expr* q = e->data.function.args[1];
    if (p->type != EXPR_INTEGER || q->type != EXPR_INTEGER) return false;
    if (q->data.integer != 2) return false;   /* reduced -> numerator is odd */
    *two_e = p->data.integer;
    return true;
}

/* Collect distinct radicands g of Power[g, 1/2] subexpressions of e.  Pointers
 * are borrowed into e; *out is realloc'd (caller free()s the array only). */
static void collect_sqrt_radicands(Expr* e, Expr*** out, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* half = make_rational(1, 2);
        bool m = expr_eq(e->data.function.args[1], half);
        expr_free(half);
        if (m) {
            Expr* g = e->data.function.args[0];
            bool dup = false;
            for (size_t i = 0; i < *n; i++)
                if (expr_eq((*out)[i], g)) { dup = true; break; }
            if (!dup) {
                if (*n == *cap) {
                    *cap = *cap ? *cap * 2 : 4;
                    *out = (Expr**)realloc(*out, *cap * sizeof(Expr*));
                }
                (*out)[(*n)++] = g;
            }
        }
    }
    collect_sqrt_radicands(e->data.function.head, out, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect_sqrt_radicands(e->data.function.args[i], out, n, cap);
}

/* Return true if the real part of cand/Sqrt[g] sampled at a fixed rational
 * point is decisively -1, i.e. cand must be negated to equal +Sqrt[g] (the two
 * differ at most by a global sign, the only branch freedom).  Borrows all. */
static bool radical_sign_negative(Expr* cand, Expr* g, Expr* x) {
    static const long pn[6] = { 17, 23, 31, 37, 19, 43 };
    static const long pd[6] = {  5,  7,  9, 10,  6,  8 };
    Expr* ratio = mk_fn2("Times", expr_copy(cand),
                         mk_inv(mk_sqrt_expr(expr_copy(g))));
    Expr* half = make_rational(1, 2);
    bool neg = false;
    for (int i = 0; i < 6; i++) {
        Expr* pt = make_rational(pn[i], pd[i]);
        Expr* at = subst_eval(expr_copy(ratio), x, pt);   /* consumes pt */
        Expr* re = at ? eval_take(mk_fn2("N", mk_fn1("Re", at), mk_int(20))) : NULL;
        if (re) {
            Expr* big = eval_take(mk_fn2("Greater",
                                  mk_fn1("Abs", expr_copy(re)), expr_copy(half)));
            if (expr_is_true(big)) {                     /* exact ratio is +/-1 */
                Expr* lt0 = eval_take(mk_fn2("Less", expr_copy(re), mk_int(0)));
                neg = expr_is_true(lt0);
                if (lt0) expr_free(lt0);
                expr_free(big); expr_free(re);
                break;
            }
            expr_free(big); expr_free(re);
        }
    }
    expr_free(half); expr_free(ratio);
    return neg;
}

/* If the radicand g (a rational function of x) is R times a perfect rational
 * square, return the single-radical form s*Sqrt[R] equal to Sqrt[g] (with the
 * global sign fixed numerically); else NULL.  Borrows g, R, x. */
static Expr* sqrt_over_R(Expr* g, Expr* R, Expr* x) {
    /* h = Factor[Cancel[Together[g/R]]] -- factoring exposes the square part. */
    Expr* h = eval_take(mk_fn1("Factor",
                  canonic(gdiv(expr_copy(g), expr_copy(R)))));
    if (!h) return NULL;
    /* s = PowerExpand[Sqrt[h]] is radical-free iff h is a perfect rational
     * square (times a constant whose own root PowerExpand resolves). */
    Expr* s = eval_take(mk_fn1("PowerExpand", mk_sqrt_expr(expr_copy(h))));
    bool ok = s && !has_sqrt(s);
    if (ok) {
        Expr* chk = gsub(mk_pow_int(expr_copy(s), 2), expr_copy(h));
        ok = is_zero(chk);
        expr_free(chk);
    }
    expr_free(h);
    if (!ok) { if (s) expr_free(s); return NULL; }
    Expr* cand = canonic(mk_fn2("Times", s, mk_sqrt_expr(expr_copy(R))));
    if (!cand) return NULL;
    if (radical_sign_negative(cand, g, x)) cand = eval_take(mk_neg(cand));
    return cand;
}

/* Re-express the antiderivative `res` so every nested Sqrt[g] that is a rational
 * multiple of Sqrt[R] is written (rational)*Sqrt[R].  Returns a new owned tree,
 * or NULL when nothing is rewritten.  Borrows res, R, x. */
static Expr* reexpress_over_radical(Expr* res, Expr* R, Expr* x) {
    Expr** rad = NULL; size_t nr = 0, cap = 0;
    collect_sqrt_radicands(res, &rad, &nr, &cap);
    Expr* out = NULL;
    for (size_t i = 0; i < nr; i++) {
        Expr* g = rad[i];
        if (expr_eq(g, R)) continue;            /* already over Sqrt[R] */
        Expr* cand = sqrt_over_R(g, R, x);
        if (!cand) continue;
        Expr* rule = mk_rule(mk_sqrt_expr(expr_copy(g)), cand);   /* consumes cand */
        Expr* base = out ? out : res;
        Expr* nx = eval_take(internal_replace_all(
                       (Expr*[]){ expr_copy(base), rule }, 2));
        if (out) expr_free(out);
        out = nx;
    }
    free(rad);
    return out;
}

/* ---------------------------------------------------------------------- */
/* Split-radical recombination                                            */
/*                                                                        */
/* reexpress_over_radical above heals an ISOLATED Sqrt[g] node.  The       */
/* eigendescent can instead leave the single radical Sqrt[R] split across  */
/* a Times as a PRODUCT/QUOTIENT of factor-wise roots -- e.g.              */
/*   Sqrt[R] / (Sqrt[x] Sqrt[(1-x)(1-k^2 x)])                              */
/* for R = x(1-x)(1-k^2 x).  Each factor alone is not a rational multiple  */
/* of Sqrt[R] (x/R and (1-x)(1-k^2 x)/R are not perfect squares), so the   */
/* per-node pass cannot touch it -- yet the COMBINED radicand is R^2 and    */
/* the whole term collapses to a pure rational.  Left split it is both a    */
/* spurious branch point (an apparent Sqrt[x] singularity the true         */
/* antiderivative does not have) and several independent radical generators */
/* that make a downstream Simplify blow up on multivariate GCDs.  The pass  */
/* below merges the x-dependent half-power factors of each Times into one   */
/* radicand and reduces it over R, restoring the canonical single-generator */
/* form.  Correct-by-construction: an algebraic identity on the integrand's */
/* domain, numeric global-sign fixed here and re-verified by diff_back_ok.  */
/* ---------------------------------------------------------------------- */

/* Numeric global-sign of a relative to b: pin the free parameters (variables
 * other than x) to fixed rationals, sample x, and return +1 if a/b is decisively
 * +1, -1 if decisively -1, else 0 (indeterminate -- caller must not rewrite).
 * The only branch freedom after merging roots is a global sign, so this pins it.
 * Borrows a, b, x. */
static int numeric_sign_ratio(Expr* a, Expr* b, Expr* x) {
    Expr* ratio = gdiv(expr_copy(a), expr_copy(b));
    Expr* vars = eval_take(mk_fn1("Variables", expr_copy(ratio)));
    if (vars && head_is(vars, SYM_List)) {
        static const long an[8] = { 12, 17, 23, 29, 31, 37, 41, 43 };
        static const long ad[8] = {  7,  5,  9, 11, 13, 10,  6,  8 };
        size_t j = 0;
        for (size_t i = 0; ratio && i < vars->data.function.arg_count; i++) {
            Expr* v = vars->data.function.args[i];
            if (expr_eq(v, x)) continue;
            ratio = subst_eval(ratio, v, make_rational(an[j % 8], ad[j % 8])); j++;
        }
    }
    if (vars) expr_free(vars);
    if (!ratio) return 0;

    static const long pn[6] = { 17, 23, 31, 37, 19, 43 };
    static const long pd[6] = {  5,  7,  9, 10,  6,  8 };
    Expr* eps = make_rational(1, 1000000L);
    int sign = 0;
    for (int i = 0; i < 6 && sign == 0; i++) {
        Expr* at = subst_eval(expr_copy(ratio), x, make_rational(pn[i], pd[i]));
        if (!at) continue;
        Expr* dp = eval_take(mk_fn2("N", mk_fn1("Abs",
                       gsub(expr_copy(at), mk_int(1))), mk_int(20)));
        Expr* dm = eval_take(mk_fn2("N", mk_fn1("Abs",
                       gadd(expr_copy(at), mk_int(1))), mk_int(20)));
        expr_free(at);
        if (dp) {
            Expr* lt = eval_take(mk_fn2("Less", expr_copy(dp), expr_copy(eps)));
            if (expr_is_true(lt)) sign = 1;
            expr_free(lt);
        }
        if (sign == 0 && dm) {
            Expr* lt = eval_take(mk_fn2("Less", expr_copy(dm), expr_copy(eps)));
            if (expr_is_true(lt)) sign = -1;
            expr_free(lt);
        }
        if (dp) expr_free(dp);
        if (dm) expr_free(dm);
    }
    expr_free(eps); expr_free(ratio);
    return sign;
}

/* Reduce Sqrt[rho] (rho a rational function of x) to canonical form over R:
 *   - a rational,          if rho is a perfect rational square;
 *   - a rational * Sqrt[R], if rho/R is a perfect rational square;
 *   - else NULL.
 * Result equals +Sqrt[rho] up to the global sign the caller fixes numerically.
 * Borrows rho, R, x. */
static Expr* reduce_sqrt_over_R(Expr* rho, Expr* R, Expr* x) {
    /* Perfect square?  Factor exposes the square part; PowerExpand takes the
     * root, which is radical-free exactly when rho is a rational square. */
    Expr* h = eval_take(mk_fn1("Factor", expr_copy(rho)));
    if (h) {
        Expr* s = powerexpand_e(mk_sqrt_expr(expr_copy(h)));
        if (s && !has_sqrt(s)) {
            Expr* chk = gsub(mk_pow_int(expr_copy(s), 2), expr_copy(h));
            bool ok = is_zero(chk);
            expr_free(chk);
            if (ok) { expr_free(h); return s; }
        }
        if (s) expr_free(s);
        expr_free(h);
    }
    /* Otherwise fall to the square*R case (sqrt_over_R returns rational*Sqrt[R],
     * already global-sign fixed against Sqrt[rho]). */
    return sqrt_over_R(rho, R, x);
}

/* Merge the x-dependent square-root factors of a single Times node t into the
 * one generator Sqrt[R].  Gather every Power[base, m/2] factor whose base
 * involves x; with two or more, form the combined rational radicand rho (bases
 * of +m/2 into the numerator, -m/2 into the denominator), reduce Sqrt[rho] over
 * R, reattach the integer-power leftovers and all non-radical (and constant-
 * radical) factors, then numeric-sign-fix the whole summand against t.  Returns a
 * NEW node equal to t in canonical single-generator form, or NULL when nothing
 * merges or the merged radicand does not reduce over R.  Borrows t, R, x. */
static Expr* combine_times_radicals(Expr* t, Expr* R, Expr* x) {
    if (!head_is(t, SYM_Times)) return NULL;
    size_t n = t->data.function.arg_count;

    Expr* num  = mk_int(1);   /* product of bases carrying a +1/2 part */
    Expr* den  = mk_int(1);   /* product of bases carrying a -1/2 part */
    Expr* rest = mk_int(1);   /* non-radical factors + integer-power leftovers */
    int nrad = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* f = t->data.function.args[i];
        int64_t two_e = 0;
        if (head_is(f, SYM_Power) && f->data.function.arg_count == 2
            && is_half_power(f->data.function.args[1], &two_e)
            && !expr_free_of(f->data.function.args[0], x)) {
            Expr* base = f->data.function.args[0];
            int64_t sgn = (two_e > 0) ? 1 : -1;
            if (sgn > 0) num = gmul(num, expr_copy(base));
            else         den = gmul(den, expr_copy(base));
            int64_t k = (two_e - sgn) / 2;            /* integer leftover power */
            if (k != 0) rest = gmul(rest, mk_pow_int(expr_copy(base), k));
            nrad++;
        } else {
            rest = gmul(rest, expr_copy(f));
        }
    }
    if (nrad < 2) { expr_free(num); expr_free(den); expr_free(rest); return NULL; }

    Expr* rho  = canonic(gdiv(num, den));             /* consumes num, den */
    Expr* repl = rho ? reduce_sqrt_over_R(rho, R, x) : NULL;
    if (rho) expr_free(rho);
    if (!repl) { expr_free(rest); return NULL; }

    Expr* cand = eval_take(gmul(rest, repl));          /* consumes rest, repl */
    if (!cand) return NULL;
    int sgn = numeric_sign_ratio(cand, t, x);
    if (sgn == 0) { expr_free(cand); return NULL; }    /* unverifiable -> leave */
    if (sgn < 0)  cand = eval_take(mk_neg(cand));
    return cand;
}

/* Functional rebuild of e that merges split square-root products at every Times
 * node (combine_times_radicals) into the single generator Sqrt[R].  Returns a
 * NEW owned tree (a copy when nothing merges).  Borrows e, R, x. */
static Expr* combine_radicals(Expr* e, Expr* R, Expr* x) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    Expr* head = combine_radicals(e->data.function.head, R, x);
    size_t n = e->data.function.arg_count;
    Expr** args = (n > 0) ? (Expr**)malloc(n * sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < n; i++)
        args[i] = combine_radicals(e->data.function.args[i], R, x);
    Expr* node = expr_new_function(head, args, n);     /* adopts head + args[i] */
    if (args) free(args);
    if (head_is(node, SYM_Times)) {
        Expr* merged = combine_times_radicals(node, R, x);
        if (merged) { expr_free(node); node = merged; }
    }
    return node;
}

/* ---------------------------------------------------------------------- */
/* Core driver                                                            */
/* ---------------------------------------------------------------------- */

static Expr* gs_core(Expr* f, Expr* x) {
    if (x->type != EXPR_SYMBOL)   return NULL;
    if (expr_free_of(f, x))       { gs_log("decline: integrand is free of the variable"); return NULL; }
    if (gs_depth >= GS_MAX_DEPTH) { gs_log("decline: recursion depth limit (%d) reached", GS_MAX_DEPTH); return NULL; }

    gs_log("------------------------------------------------------------");
    gs_log_expr("integrand f", f);
    gs_log_expr("variable", x);

    Expr* F = NULL; Expr* R = NULL;
    int pnum = 0, pden = 0;
    if (!recognise(f, x, &F, &R, &pnum, &pden)) {
        gs_log("FORM: no match -- no single radical R(x)^(n/d) with reduced d in {2,3,4}; declining");
        return NULL;
    }
    gs_log("FORM: matched -- reduced to canonical F(x)*R(x)^(-p), p = %d/%d", pnum, pden);
    gs_log("       (any rational radical exponent n/d, radical in numerator OR denominator, folds here:");
    gs_log("        R^(n/d) = R^k * R^(-p), the integer R^k absorbed into F)");
    gs_log_expr("  cofactor F", F);
    gs_log_expr("  radicand R", R);

    Expr* result = NULL;
    if (pden == 2) {
        gs_log("BRANCH: square-root case -- Goursat V4 (order-2 Mobius, Theorems 1-2)");
        result = goursat_v4(F, R, x);
        /* When V4 (Theorems 1-2) declines a square-root-of-cubic integrand, try
         * the period-3 higher-symmetry reduction (Goursat 1887 Section 4), which
         * covers the t^3-1-symmetric cases V4 misses (e.g.
         * (t-1)/((t+2) Sqrt[t^3-1])). */
        if (!result) {
            gs_log("V4 declined -- trying period-3 higher-symmetry reduction (Goursat 1887 Sec 4)");
            result = goursat_period3(F, R, x);
        }
    } else if (pden == 3) {
        gs_log("BRANCH: cube-root case (p = %d/3)", pnum);
        /* Try the constructive THIRD-KIND logarithmic-part reduction first: it is
         * cheap (no Mobius/Solve over the splitting field) and self-verifying
         * (diff_back_ok), and it covers the F-with-a-non-branch-pole cases the
         * order-3 eigendescent obstructs (H1 != 0).  Doing it first also avoids
         * the eigendescent's expensive -- for a parametric radicand, budget-
         * exhausting -- cyclic-Mobius canonicalisation on integrands it cannot
         * close anyway.  Fall back to the eigendescent (Theorem: H1 == 0) when
         * the log-sum ansatz does not verify. */
        gs_log("trying third-kind logarithmic-part reduction");
        result = goursat_cubic_thirdkind(F, R, x, pnum);
        if (!result) {
            gs_log("third-kind declined -- trying order-3 Mobius eigendescent");
            result = goursat_cubic(F, R, x, pnum);
        }
    } else if (pden == 4) {
        gs_log("BRANCH: fourth-root case -- order-4 Mobius eigendescent (p = %d/4)", pnum);
        result = goursat_quartic(F, R, x, pnum);
    }
    gs_log(result ? "descent produced a candidate antiderivative"
                  : "descent did not close -- no candidate");

    /* Normalise the result's nested radical over the integrand's own Sqrt[R]
     * (square-root case only), so D[result,x] // Simplify reduces to f.  Adopt
     * the rewrite only if the differentiate-back guard still passes on it. */
    if (result && pden == 2) {
        gs_log("reexpress_over_radical: START");
        Expr* improved = reexpress_over_radical(result, R, x);
        gs_log("reexpress_over_radical: DONE (%s)", improved ? "rewrote" : "no change");
        if (improved) {
            if (diff_back_ok(improved, x, f)) { expr_free(result); result = improved; }
            else expr_free(improved);
        }

        /* reexpress_over_radical re-expresses each nested radical as
         * rational*Sqrt[R] one node at a time; when the evaluator has already
         * distributed a factored radicand across a term, that leaves the single
         * radical Sqrt[R] SPLIT as a product of factor roots
         * (Sqrt[x] Sqrt[(1-x)(1-k^2 x)] Sqrt[R] / ... instead of Sqrt[R]).  That
         * split is a spurious branch point AND several independent radical
         * generators that make a downstream Simplify blow up on multivariate
         * GCDs.  It is only visible on the fully-evaluated form, so evaluate and
         * then recombine the x-dependent factor roots of each Times back into the
         * single generator Sqrt[R].  Correct by construction (an identity on the
         * integrand's domain, numeric global-sign fixed); adopt only if the
         * differentiate-back guard still passes. */
        gs_log("combine_radicals: START");
        Expr* flat = evaluate(expr_copy(result));
        Expr* merged = flat ? combine_radicals(flat, R, x) : NULL;
        gs_log("combine_radicals: DONE (%s)",
               (merged && flat && !expr_eq(merged, flat)) ? "rewrote" : "no change");
        if (merged) {
            if (flat && !expr_eq(merged, flat) && diff_back_ok(merged, x, f)) {
                expr_free(result); result = merged;
            } else expr_free(merged);
        }
        if (flat) expr_free(flat);
    }

    expr_free(F);
    expr_free(R);

    /* Safety guard (differentiate-back).  The eigenspace elementarity criterion
     * relies on PossibleZeroQ, which can misfire on deeply nested radical roots
     * (e.g. an irreducible quartic solved via Ferrari) -- catastrophic
     * cancellation can numericalise a genuinely nonzero projection to ~0, so the
     * criterion wrongly passes and an incorrect antiderivative (often the
     * degenerate 0) is emitted.  Verify D[result, x] == f before returning (see
     * diff_back_ok for why this is a numeric rather than a PossibleZeroQ check);
     * on failure, decline so the cascade continues. */
    if (result && !diff_back_ok(result, x, f)) {
        gs_log("VERIFY: differentiate-back guard FAILED -- discarding candidate, declining");
        expr_free(result); result = NULL;
    } else if (result) {
        gs_log("VERIFY: differentiate-back guard passed");
        gs_log_expr("RESULT", result);
    }
    return result;
}

/* ---------------------------------------------------------------------- */
/* Termination guard                                                      */
/* ---------------------------------------------------------------------- */
/*
 * The cube-/fourth-root descents canonicalise rational functions over the
 * splitting field of R.  When R has cyclotomic roots (e.g. 1 - x^3, whose
 * roots are the cube roots of unity) the core Together/Cancel/Simplify over
 * that algebraic-number field can blow up super-polynomially for unlucky
 * cofactors: Sqrt[1 - x^3]/x and (1 - x^3)^(1/3)/(1 + x) both spin for minutes
 * inside Cancel/Simplify even though the same machinery dispatches
 * (1 - x^3)^(1/3)/x in a tenth of a second.  Because try_goursat sits in the
 * *default* Integrate cascade, an unbounded attempt would hang every
 * Integrate[] of such an integrand -- not only explicit
 * Method -> "GoursatAlgebraic" calls.
 *
 * Bound the whole outermost attempt with a CPU-time budget -- exactly the
 * remedy the TimeConstrained docs prescribe for an open-ended computation --
 * and decline (leave the integral for the rest of the cascade) on exhaustion.
 * TimeConstrained's SIGPROF/ITIMER_PROF abort can interrupt the in-flight C
 * canonicaliser (a purely cooperative check could not), which is what makes
 * termination a hard guarantee on Linux/macOS.  Only the OUTERMOST attempt
 * (gs_depth == 0) is wrapped; the recursive genus-0 descents run inside that
 * single timer.  The documented price is a bounded leak of the abandoned
 * computation on the rare timeout path.
 */
#define GS_TIME_BUDGET_SEC 6
#define GS_RUN_HEAD "Integrate`GoursatAlgebraic`Run"

/* Internal, cascade-invisible builtin that runs the unguarded core so the
 * public entries can evaluate it under TimeConstrained.  Returns the symbol
 * $Failed (never a valid antiderivative) when the core declines, so the guard
 * sees one sentinel for both "declined" and "timed out". */
static Expr* builtin_gs_run(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return mk_sym("$Failed");
    Expr* r = gs_core(res->data.function.args[0], res->data.function.args[1]);
    return r ? r : mk_sym("$Failed");
}

/* Run gs_core under a CPU-time budget at the outermost level; nested recursive
 * calls (gs_depth > 0) run directly under the already-armed timer.  Borrows
 * f, x; returns the owned antiderivative or NULL (decline / timeout). */
static Expr* gs_guarded(Expr* f, Expr* x) {
    if (gs_depth > 0) return gs_core(f, x);
    gs_set_debug();   /* latch Integrate`GoursatDebug for this whole attempt */
    Expr* call = mk_fn3("TimeConstrained",
        mk_fn2(GS_RUN_HEAD, expr_copy(f), expr_copy(x)),
        mk_int(GS_TIME_BUDGET_SEC),
        mk_sym("$Failed"));
    Expr* r = eval_take(call);
    if (!r) return NULL;
    Expr* failed = mk_sym("$Failed");
    bool decline = expr_eq(r, failed);
    expr_free(failed);
    if (decline) { expr_free(r); return NULL; }
    return r;
}

/* ---------------------------------------------------------------------- */
/* Public entry points                                                    */
/* ---------------------------------------------------------------------- */

Expr* integrate_goursat_try(Expr* f, Expr* x) {
    return gs_guarded(f, x);
}

Expr* builtin_integrate_goursat(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;
    return gs_guarded(res->data.function.args[0], res->data.function.args[1]);
}

void integrate_goursat_init(void) {
    symtab_add_builtin(GS_RUN_HEAD, builtin_gs_run);
    symtab_get_def(GS_RUN_HEAD)->attributes |=
        ATTR_HOLDALL | ATTR_PROTECTED | ATTR_READPROTECTED;

    /* User-settable trace switch.  Default False; when set True the descent
     * narrates its progress (form recognition, criterion tests, reductions) to
     * stderr.  Left unprotected so the user can assign it. */
    {
        Expr* pat = expr_new_symbol("Integrate`GoursatDebug");
        Expr* val = expr_new_symbol(SYM_False);
        symtab_add_own_value("Integrate`GoursatDebug", pat, val);
        expr_free(pat); expr_free(val);
    }
    symtab_set_docstring("Integrate`GoursatDebug",
        "Integrate`GoursatDebug is a Boolean switch (default False); when set True,\n"
        "Integrate`GoursatAlgebraic prints a trace of form recognition, the\n"
        "involution/eigenspace criterion tests, and the recursive reductions to stderr.");

    symtab_add_builtin("Integrate`GoursatAlgebraic", builtin_integrate_goursat);
    symtab_get_def("Integrate`GoursatAlgebraic")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`GoursatAlgebraic",
        "Integrate`GoursatAlgebraic[f, x] integrates pseudo-elliptic integrands\n"
        "F(x) R(x)^(n/d) (F rational, R a polynomial, n/d any rational of reduced\n"
        "denominator d in {2,3,4}) by Goursat's algorithm and its cube-/fourth-root\n"
        "generalisations.  The radical may be in the numerator or the denominator:\n"
        "R^(n/d) is reduced to R^k R^(-p) with p in {1/2,1/3,2/3,1/4,3/4} and the\n"
        "integer power R^k folded into F.  A Mobius automorphism cycling the roots of\n"
        "R splits the integrand into eigencomponents that descend to genus-0 curves\n"
        "when the elementarity criterion holds.  Strict: returns unevaluated when f\n"
        "is not of this form, the integral is non-elementary, or the roots are not\n"
        "solvable in radicals.");
}
