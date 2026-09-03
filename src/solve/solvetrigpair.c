/*
 * solvetrigpair.c
 *
 * Argument-pair reducer for generalized two-argument trig / hyperbolic
 * equations.  See solvetrigpair.h for the algorithm narrative.
 *
 * Reuses the solveinv specialist to solve each single-argument atom, so
 * the periodic families (ArcXxx + 2 Pi C[k], the C[k] parameter, the
 * ConditionalExpression / Element[C[k], Integers] wrapping) are produced
 * by the same code path the single-head inverter already uses.
 */

#include "solvetrigpair.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "parse.h"
#include "solveinv.h"
#include "solvepoly.h"
#include "solvetrig.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "symtab.h"
#include "trig_canon.h"
#include "zero_test.h"

/* ------------------------------------------------------------------ *
 *  Tiny construction helpers (mirrors solvetrig.c style).            *
 * ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }
static Expr* mk_fn1(const char* head, Expr* a) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* mk_pow(Expr* base, Expr* exp) {
    return mk_fn2("Power", base, exp);
}
static Expr* mk_neg(Expr* e) {
    return mk_fn2("Times", mk_int(-1), e);
}

static bool head_is_sym(const Expr* e, const char* interned_head) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == interned_head;
}

static bool var_in(const Expr* e, const Expr* var) {
    if (!e || !var) return false;
    if (e->type == EXPR_SYMBOL && var->type == EXPR_SYMBOL
        && e->data.symbol.name == var->data.symbol.name) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (var_in(e->data.function.head, var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (var_in(e->data.function.args[i], var)) return true;
    return false;
}

/* Primary trig head (Sin/Cos/Sinh/Cosh) -- the only heads that appear in
 * the Together denominator once the reciprocal heads are normalized, and
 * thus the only pole factors we need to test. */
static bool is_primary_trig_head(const char* h) {
    return h == SYM_Sin || h == SYM_Cos || h == SYM_Sinh || h == SYM_Cosh;
}

/* ------------------------------------------------------------------ *
 *  Cached rule tables (parsed once in _init).                        *
 * ------------------------------------------------------------------ */

static Expr* pair_to_sincos   = NULL;  /* Tan/Sec/... -> Sin/Cos ratios  */
static Expr* pair_factor      = NULL;  /* sum-to-product / reverse-angle */
static Expr* pair_expand_exp  = NULL;  /* Power[E, e] -> Power[E, Expand[e]] */

/* Apply a cached rule list to `e` (consumed) via ReplaceAll; evaluate. */
static Expr* apply_rules(Expr* e, Expr* rules) {
    if (!rules) return e;
    return eval_and_free(mk_fn2("ReplaceAll", e, expr_copy(rules)));
}

/* ------------------------------------------------------------------ *
 *  Concatenate two solution Lists (a and b consumed).                *
 * ------------------------------------------------------------------ */

static Expr* concat_solutions(Expr* a, Expr* b) {
    if (!a) return b;
    if (!b) return a;
    if (!head_is_sym(a, SYM_List)) { expr_free(b); return a; }
    if (!head_is_sym(b, SYM_List)) { expr_free(b); return a; }
    size_t na = a->data.function.arg_count;
    size_t nb = b->data.function.arg_count;
    size_t n  = na + nb;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < na; i++) args[i]      = expr_copy(a->data.function.args[i]);
    for (size_t i = 0; i < nb; i++) args[na + i] = expr_copy(b->data.function.args[i]);
    Expr* out = expr_new_function(mk_sym("List"), args, n);
    free(args);
    expr_free(a); expr_free(b);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Solve the factored numerator: each var-bearing atom H[arg] == 0   *
 *  goes to the single-head inverse specialist; union the branches.   *
 * ------------------------------------------------------------------ */

static Expr* solve_atom(Expr* atom, Expr* var, Expr* dom,
                        const SolveInvOpts* opts) {
    Expr* eq = eval_and_free(mk_fn2("Equal", expr_copy(atom), mk_int(0)));
    /* Single invertible head first (H[g(var)] == 0), then a polynomial in
     * a single kernel (e.g. 2 Cosh[x] Cosh[y] - 2, or 1 + Cosh[A]^2 from a
     * Sech + Cosh collapse) which the isolator cannot peel. */
    Expr* s = solveinv_solve_inverse_equality(eq, var, dom, opts);
    if (!s) s = solvetrig_solve_poly_in_kernel(eq, var, dom, opts);
    expr_free(eq);
    return s;
}

static Expr* solve_numer(Expr* Nf, Expr* var, Expr* dom,
                         const SolveInvOpts* opts) {
    if (head_is_sym(Nf, SYM_Times)) {
        Expr* agg = NULL;
        bool any = false;
        for (size_t i = 0; i < Nf->data.function.arg_count; i++) {
            Expr* f = Nf->data.function.args[i];
            if (!var_in(f, var)) continue;   /* drop var-free factors */
            any = true;
            Expr* s = solve_atom(f, var, dom, opts);
            if (!s) { if (agg) expr_free(agg); return NULL; }
            agg = concat_solutions(agg, s);
        }
        if (!any) return NULL;
        return agg;
    }
    /* Single atom, or a Plus like (1 - Sin[A]) from a reciprocal trap --
     * hand the whole numerator to the isolator, which peels the lone
     * var-bearing head out of the surrounding var-free constants. */
    return solve_atom(Nf, var, dom, opts);
}

/* ------------------------------------------------------------------ *
 *  Pole / spurious gate.                                             *
 * ------------------------------------------------------------------ */

/* Collect primary-trig subexpressions over `var` from denominator `D`. */
static void collect_pole_factors(const Expr* D, const Expr* var,
                                 Expr*** arr, size_t* n, size_t* cap) {
    if (!D || D->type != EXPR_FUNCTION) return;
    if (D->data.function.head->type == EXPR_SYMBOL
        && is_primary_trig_head(D->data.function.head->data.symbol.name)
        && D->data.function.arg_count == 1
        && var_in(D->data.function.args[0], var)) {
        if (*n == *cap) {
            *cap = *cap ? *cap * 2 : 4;
            *arr = (Expr**)realloc(*arr, sizeof(Expr*) * (*cap));
        }
        (*arr)[(*n)++] = (Expr*)D;
        return;   /* no trig nests inside a primary-trig arg we care about */
    }
    for (size_t i = 0; i < D->data.function.arg_count; i++)
        collect_pole_factors(D->data.function.args[i], var, arr, n, cap);
}

/* Value carried by a solution RHS (unwrap ConditionalExpression). */
static Expr* cond_value(Expr* rhs) {
    if (head_is_sym(rhs, SYM_ConditionalExpression)
        && rhs->data.function.arg_count >= 1)
        return rhs->data.function.args[0];
    return rhs;
}

/* True iff some pole factor vanishes identically on this solution family.
 * Representative: set every generated parameter C[k] -> 0 (the denominator
 * is periodic with the family period, so one representative decides). */
static bool family_hits_pole(Expr* value, Expr** poles, size_t np,
                             const Expr* var, const char* param_head) {
    /* rep = value /. param_head[_] -> 0 */
    Expr* blank   = expr_new_function(mk_sym("Blank"), NULL, 0);
    Expr* pat     = expr_new_function(mk_sym(param_head),
                                      (Expr*[]){ blank }, 1);
    Expr* zrule   = mk_fn2("Rule", pat, mk_int(0));
    Expr* rep = eval_and_free(mk_fn2("ReplaceAll", expr_copy(value), zrule));

    bool polar = false;
    for (size_t i = 0; i < np && !polar; i++) {
        Expr* vrule = mk_fn2("Rule", expr_copy((Expr*)var), expr_copy(rep));
        Expr* sub   = eval_and_free(
            mk_fn2("ReplaceAll", expr_copy(poles[i]), vrule));
        if (zero_test_decide(sub) == ZERO_TEST_TRUE) polar = true;
        expr_free(sub);
    }
    expr_free(rep);
    return polar;
}

static Expr* pole_gate(Expr* sols, Expr* D, Expr* var,
                       const char* param_head) {
    if (!sols || !head_is_sym(sols, SYM_List)) return sols;
    Expr** poles = NULL; size_t np = 0, cap = 0;
    collect_pole_factors(D, var, &poles, &np, &cap);
    if (np == 0) { free(poles); return sols; }

    size_t nrows = sols->data.function.arg_count;
    Expr** kept = (Expr**)malloc(sizeof(Expr*) * (nrows ? nrows : 1));
    size_t nk = 0;
    for (size_t i = 0; i < nrows; i++) {
        Expr* row = sols->data.function.args[i];
        Expr* value = NULL;
        if (head_is_sym(row, SYM_List) && row->data.function.arg_count == 1) {
            Expr* rule = row->data.function.args[0];
            if (head_is_sym(rule, SYM_Rule)
                && rule->data.function.arg_count == 2)
                value = cond_value(rule->data.function.args[1]);
        }
        bool drop = value
            && family_hits_pole(value, poles, np, var, param_head);
        if (!drop) kept[nk++] = expr_copy(row);
    }
    Expr* out = expr_new_function(mk_sym("List"), kept, nk);
    free(kept);
    free(poles);
    expr_free(sols);
    return out;
}

/* ================================================================== *
 *  Engine 2: single rational exponential generator.                  *
 *                                                                    *
 *  For residuals that do NOT factor into single-argument atoms       *
 *  (nonzero-constant RHS like Tan[A]-Tan[B]==1, or mixed Sinh/Cosh),  *
 *  rewrite via TrigToExp and substitute a single generator           *
 *  u = E^(sigma*var), sigma = (I or 1) * g, g = the rational gcd of   *
 *  the var-coefficients of the exponentials.  The var-free part of    *
 *  each exponent (e.g. E^(I x)) folds into the polynomial             *
 *  coefficients, so the result is a Laurent polynomial in u; solve    *
 *  it and unwind each root through E^(sigma*var) == u0 (peel_exp).    *
 *  Only numeric var-coefficients admit a generator -- symbolic ones   *
 *  (Sinh[a x+b y]+...) decline here (Engine 1's exclusive domain).    *
 * ------------------------------------------------------------------ */

#define GEN_USYM "Solve`trigpairU"

static int64_t igcd64(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}
static int64_t ilcm64(int64_t a, int64_t b) {
    if (a == 0 || b == 0) return 0;
    int64_t g = igcd64(a, b);
    a /= g;
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    return a * b;
}

/* Rational value of a numeric Expr (Integer or Rational[p,q]). */
static bool expr_rational(const Expr* e, int64_t* num, int64_t* den) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *num = e->data.integer; *den = 1; return true; }
    if (head_is_sym(e, SYM_Rational) && e->data.function.arg_count == 2) {
        const Expr* p = e->data.function.args[0];
        const Expr* q = e->data.function.args[1];
        if (p->type == EXPR_INTEGER && q->type == EXPR_INTEGER && q->data.integer != 0) {
            *num = p->data.integer; *den = q->data.integer; return true;
        }
    }
    return false;
}

/* Decode a var-coefficient kappa into (is_imag, rational num/den).
 * kappa is Integer, Rational, or Complex[0, Integer|Rational]. */
static bool kappa_decode(const Expr* k, bool* is_imag, int64_t* num, int64_t* den) {
    if (head_is_sym(k, SYM_Complex) && k->data.function.arg_count == 2) {
        const Expr* re = k->data.function.args[0];
        if (!(re->type == EXPR_INTEGER && re->data.integer == 0)) return false;
        *is_imag = true;
        return expr_rational(k->data.function.args[1], num, den);
    }
    *is_imag = false;
    return expr_rational(k, num, den);
}

typedef struct {
    int64_t gnum, gden;   /* running rational gcd num/den (0/1 = none yet)   */
    bool    is_imag, have; /* uniform imaginary flag                          */
    bool    bad;          /* var occurred bare, or a non-numeric coefficient  */
} GenScan;

/* Reduce a rational in place. */
static void rat_reduce(int64_t* n, int64_t* d) {
    if (*d < 0) { *n = -*n; *d = -*d; }
    int64_t g = igcd64(*n, *d);
    if (g > 1) { *n /= g; *d /= g; }
}

/* Fold rational gcd: g = gcd(gnum,num)/lcm(gden,den). */
static void gcd_accumulate(GenScan* s, int64_t num, int64_t den) {
    rat_reduce(&num, &den);
    if (!s->have) { s->gnum = num < 0 ? -num : num; s->gden = den; s->have = true; return; }
    int64_t nn = igcd64(s->gnum, num);
    int64_t dd = ilcm64(s->gden, den);
    s->gnum = nn; s->gden = dd;
}

/* First pass: validate structure and accumulate the generator step. */
static void scan_generator(const Expr* e, const Expr* var, GenScan* s) {
    if (s->bad) return;
    if (e->type == EXPR_SYMBOL && var->type == EXPR_SYMBOL
        && e->data.symbol.name == var->data.symbol.name) { s->bad = true; return; }
    if (e->type != EXPR_FUNCTION) return;
    if (head_is_sym(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp_ = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == SYM_E
            && var_in(exp_, var)) {
            Expr* k = eval_and_free(mk_fn2("Coefficient",
                expr_copy((Expr*)exp_), expr_copy((Expr*)var)));
            bool imag; int64_t nu, de;
            bool ok = kappa_decode(k, &imag, &nu, &de);
            expr_free(k);
            if (!ok || nu == 0) { s->bad = true; return; }
            /* beta = exp_ /. var -> 0 must be var-free (exponent linear). */
            Expr* z = mk_fn2("Rule", expr_copy((Expr*)var), mk_int(0));
            Expr* beta = eval_and_free(mk_fn2("ReplaceAll", expr_copy((Expr*)exp_), z));
            if (var_in(beta, var)) { expr_free(beta); s->bad = true; return; }
            expr_free(beta);
            if (s->have && s->is_imag != imag) { s->bad = true; return; }
            s->is_imag = imag;
            gcd_accumulate(s, nu, de);
            return;   /* admissible: do not recurse into the exponent */
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        scan_generator(e->data.function.args[i], var, s);
}

/* Second pass: substitute E^(kappa var + beta) -> E^beta * U^m,
 * m = (kappa / g).  Reuses `usym` for U.  Tracks min exponent. */
static Expr* subst_generator(const Expr* e, const Expr* var, const Expr* usym,
                             int64_t gnum, int64_t gden, int64_t* min_m) {
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (head_is_sym(e, SYM_Power) && e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp_ = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == SYM_E
            && var_in(exp_, var)) {
            Expr* k = eval_and_free(mk_fn2("Coefficient",
                expr_copy((Expr*)exp_), expr_copy((Expr*)var)));
            bool imag; int64_t nu, de;
            kappa_decode(k, &imag, &nu, &de);
            expr_free(k);
            rat_reduce(&nu, &de);
            /* m = (nu/de) / (gnum/gden) = (nu*gden)/(de*gnum). */
            int64_t m = (nu * gden) / (de * gnum);
            if (m < *min_m) *min_m = m;
            Expr* z = mk_fn2("Rule", expr_copy((Expr*)var), mk_int(0));
            Expr* beta = eval_and_free(mk_fn2("ReplaceAll", expr_copy((Expr*)exp_), z));
            Expr* ebeta = eval_and_free(mk_fn2("Power", mk_sym("E"), beta));
            Expr* upow  = mk_pow(expr_copy((Expr*)usym), mk_int(m));
            return mk_fn2("Times", ebeta, upow);
        }
    }
    size_t n = e->data.function.arg_count;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        args[i] = subst_generator(e->data.function.args[i], var, usym, gnum, gden, min_m);
    Expr* out = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return out;
}

static Expr* try_exp_generator(Expr* residual, Expr* var, Expr* dom,
                               const SolveInvOpts* opts) {
    /* TrigToExp -> Together -> Numerator (all under suppression is not
     * needed: TrigToExp already emits Sin/Cos-free exponentials). */
    Expr* expf = eval_and_free(mk_fn1("TrigToExp", expr_copy(residual)));
    Expr* frac = eval_and_free(internal_together((Expr*[]){ expf }, 1));
    Expr* numer = eval_and_free(internal_numerator((Expr*[]){ expr_copy(frac) }, 1));
    expr_free(frac);
    /* Together leaves exponents unexpanded (e.g. x - 2y + 2(x+y), whose net
     * y-coefficient is 0); expand each so the coefficient scan is exact. */
    numer = eval_and_free(mk_fn2("ReplaceAll", numer, expr_copy(pair_expand_exp)));
    if (!var_in(numer, var)) { expr_free(numer); return NULL; }

    /* Determine the generator step g and orientation. */
    GenScan s = { 0, 1, false, false, false };
    scan_generator(numer, var, &s);
    if (s.bad || !s.have || s.gnum == 0) { expr_free(numer); return NULL; }

    Expr* usym = mk_sym(GEN_USYM);
    int64_t min_m = 0;
    Expr* inu = subst_generator(numer, var, usym, s.gnum, s.gden, &min_m);
    expr_free(numer);

    if (min_m < 0) {
        inu = eval_and_free(mk_fn2("Times", inu,
            mk_pow(expr_copy(usym), mk_int(-min_m))));
    }
    inu = eval_and_free(internal_expand((Expr*[]){ inu }, 1));
    if (!var_in(inu, usym)) { expr_free(inu); expr_free(usym); return NULL; }

    /* Solve the polynomial in u. */
    Expr* u_eq = eval_and_free(mk_fn2("Equal", inu, mk_int(0)));
    SolvePolyOpts polyopts = { false, false };
    Expr* usol = solvepoly_solve_polynomial_equality(u_eq, usym, NULL, &polyopts);
    expr_free(u_eq);
    if (!usol || !head_is_sym(usol, SYM_List)) {
        if (usol) expr_free(usol);
        expr_free(usym); return NULL;
    }

    /* sigma = (I if imaginary else 1) * (gnum/gden). */
    Expr* g_rat = (s.gden == 1) ? mk_int(s.gnum)
                                : mk_fn2("Rational", mk_int(s.gnum), mk_int(s.gden));
    Expr* sigma = s.is_imag ? mk_fn2("Times", mk_sym("I"), g_rat) : g_rat;

    Expr* agg = NULL;
    for (size_t i = 0; i < usol->data.function.arg_count; i++) {
        Expr* sol = usol->data.function.args[i];
        if (!head_is_sym(sol, SYM_List) || sol->data.function.arg_count != 1) continue;
        Expr* rule = sol->data.function.args[0];
        if (!head_is_sym(rule, SYM_Rule) || rule->data.function.arg_count != 2) continue;
        Expr* u0 = rule->data.function.args[1];
        /* E^(sigma var) == u0  ->  peel_exp via solveinv. */
        Expr* exp_eq = eval_and_free(mk_fn2("Equal",
            mk_pow(mk_sym("E"), mk_fn2("Times", expr_copy(sigma), expr_copy(var))),
            expr_copy(u0)));
        Expr* branch = solveinv_solve_inverse_equality(exp_eq, var, dom, opts);
        expr_free(exp_eq);
        if (branch) agg = concat_solutions(agg, branch);
    }
    expr_free(sigma);
    expr_free(usol);
    expr_free(usym);
    /* The peel_exp unwind already yields the canonical Log form; a Simplify
     * pass here would only strip complex logs for the rare clean-factoring
     * case yet triggers a pre-existing Simplify per-call leak, so it is
     * intentionally omitted (the exp-generator branch is the messy-output
     * fallback anyway). */
    return agg;
}

/* ------------------------------------------------------------------ *
 *  Public entry.                                                     *
 * ------------------------------------------------------------------ */

bool solvetrigpair_looks_applicable(const Expr* expr, const Expr* var) {
    return solvetrig_has_trig(expr, var);
}

Expr* solvetrigpair_solve(Expr* equation, Expr* var, Expr* dom,
                          const SolveInvOpts* opts) {
    if (!equation || !var) return NULL;
    if (!head_is_sym(equation, SYM_Equal)
        || equation->data.function.arg_count != 2) return NULL;
    if (!solvetrig_has_trig(equation, var)) return NULL;

    Expr* lhs = equation->data.function.args[0];
    Expr* rhs = equation->data.function.args[1];
    Expr* residual = eval_and_free(mk_fn2("Plus",
        expr_copy(lhs), mk_neg(expr_copy(rhs))));
    /* Kept for the Engine-2 (exp-generator) fallback, which needs the
     * original trig-form residual. */
    Expr* resid2 = expr_copy(residual);

    /* 1. Reciprocal-normalize, Together, numerator/denominator split, and
     *    factor -- all under trig_canon suppression, otherwise the very
     *    first step (Tan -> Sin/Cos) is undone immediately by the Times-
     *    level canonicaliser (Sin/Cos collapses straight back to Tan,
     *    1/Cos back to Sec). */
    trig_canon_suppress_inc();
    Expr* sincos = apply_rules(residual, pair_to_sincos);
    Expr* frac = eval_and_free(internal_together((Expr*[]){ sincos }, 1));
    Expr* N = eval_and_free(internal_numerator((Expr*[]){ expr_copy(frac) }, 1));
    Expr* D = eval_and_free(internal_denominator((Expr*[]){ frac }, 1));

    /* 2. Factor the numerator into single-argument atoms. */
    Expr* Nf = apply_rules(N, pair_factor);
    trig_canon_suppress_dec();

    const char* ph = (opts && opts->param_head) ? opts->param_head : intern_symbol("C");

    /* 3. Degenerate cases.  Factoring collapses a parity identity to a
     * literal 0 (via reverse-angle Sin[0]/Sinh[0]) and a contradiction to a
     * constant (via the Pythagorean rules), so both are decided here without
     * calling Simplify (which has a pre-existing per-call leak on trig-of-
     * fractional-argument inputs); the zero test is leak-free. */
    if (zero_test_decide(Nf) == ZERO_TEST_TRUE) {           /* -> True ({{}}) */
        expr_free(Nf); expr_free(D); expr_free(resid2);
        Expr* empty = expr_new_function(mk_sym("List"), NULL, 0);
        return expr_new_function(mk_sym("List"), (Expr*[]){ empty }, 1);
    }
    if (!var_in(Nf, var)) {                        /* nonzero const -> False */
        expr_free(Nf); expr_free(D); expr_free(resid2);
        return expr_new_function(mk_sym("List"), NULL, 0);
    }

    /* 4. Solve each atom via the inverse-function specialist. */
    Expr* sols = solve_numer(Nf, var, dom, opts);
    expr_free(Nf);

    /* 4b. Engine 2 fallback: exponential-generator polynomial for the
     * non-factoring residuals (nonzero-constant RHS, mixed Sinh/Cosh). */
    if (!sols) sols = try_exp_generator(resid2, var, dom, opts);
    expr_free(resid2);
    if (!sols) { expr_free(D); return NULL; }

    /* 5. Pole gate (Engine 1's trig-form denominator carries the poles). */
    sols = pole_gate(sols, D, var, ph);
    expr_free(D);
    return sols;
}

/* ------------------------------------------------------------------ *
 *  Qualified-builtin entry: Solve`SolveTrigPair                      *
 * ------------------------------------------------------------------ */

Expr* builtin_solve_trig_pair(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    Expr* equation = res->data.function.args[0];
    Expr* var      = res->data.function.args[1];
    Expr* dom      = (argc >= 3) ? res->data.function.args[2] : NULL;
    SolveInvOpts opts = { true, intern_symbol("C") };
    return solvetrigpair_solve(equation, var, dom, &opts);
}

void solvetrigpair_init(void) {
    /* Reciprocal normalization (same table trigsimp uses). */
    pair_expand_exp = parse_expression("{ Power[E, e_] :> Power[E, Expand[e]] }");

    pair_to_sincos = parse_expression(
        "{ Tan[x_]  :> Sin[x] / Cos[x], Cot[x_]  :> Cos[x] / Sin[x],"
        "  Sec[x_]  :> 1 / Cos[x],      Csc[x_]  :> 1 / Sin[x],"
        "  Tanh[x_] :> Sinh[x] / Cosh[x], Coth[x_] :> Cosh[x] / Sinh[x],"
        "  Sech[x_] :> 1 / Cosh[x],     Csch[x_] :> 1 / Sinh[x] }");

    /* Factorization: Pythagorean collapse (for contradiction -> const),
     * reverse angle-addition (Tan/Cot/mixed numerators), same-head
     * sum-to-product, and mixed Sin/Cos (Sec/Csc) phase-shift forms. */
    pair_factor = parse_expression(
        "{ "
        /* Pythagorean (same argument -> constant). */
        "Cos[x_]^2 + Sin[x_]^2 + r___ :> 1 + r, "
        "Cosh[x_]^2 - Sinh[x_]^2 + r___ :> 1 + r, "
        "-Cosh[x_]^2 + Sinh[x_]^2 + r___ :> -1 + r, "
        /* Reverse angle-addition (circular). */
        "Sin[a_] Cos[b_] + Cos[a_] Sin[b_] + r___ :> Sin[a + b] + r, "
        "-Sin[a_] Cos[b_] - Cos[a_] Sin[b_] + r___ :> -Sin[a + b] + r, "
        "Sin[a_] Cos[b_] - Cos[a_] Sin[b_] + r___ :> Sin[a - b] + r, "
        "-Sin[a_] Cos[b_] + Cos[a_] Sin[b_] + r___ :> -Sin[a - b] + r, "
        "Cos[a_] Cos[b_] - Sin[a_] Sin[b_] + r___ :> Cos[a + b] + r, "
        "-Cos[a_] Cos[b_] + Sin[a_] Sin[b_] + r___ :> -Cos[a + b] + r, "
        "Cos[a_] Cos[b_] + Sin[a_] Sin[b_] + r___ :> Cos[a - b] + r, "
        /* Reverse angle-addition (hyperbolic). */
        "Sinh[a_] Cosh[b_] + Cosh[a_] Sinh[b_] + r___ :> Sinh[a + b] + r, "
        "-Sinh[a_] Cosh[b_] - Cosh[a_] Sinh[b_] + r___ :> -Sinh[a + b] + r, "
        "Sinh[a_] Cosh[b_] - Cosh[a_] Sinh[b_] + r___ :> Sinh[a - b] + r, "
        "-Sinh[a_] Cosh[b_] + Cosh[a_] Sinh[b_] + r___ :> -Sinh[a - b] + r, "
        "Cosh[a_] Cosh[b_] + Sinh[a_] Sinh[b_] + r___ :> Cosh[a + b] + r, "
        "-Cosh[a_] Cosh[b_] - Sinh[a_] Sinh[b_] + r___ :> -Cosh[a + b] + r, "
        "Cosh[a_] Cosh[b_] - Sinh[a_] Sinh[b_] + r___ :> Cosh[a - b] + r, "
        /* Same-head sum-to-product (circular). */
        "Sin[a_] + Sin[b_] + r___ :> 2 Sin[(a + b)/2] Cos[(a - b)/2] + r, "
        "Sin[a_] - Sin[b_] + r___ :> 2 Cos[(a + b)/2] Sin[(a - b)/2] + r, "
        "Cos[a_] + Cos[b_] + r___ :> 2 Cos[(a + b)/2] Cos[(a - b)/2] + r, "
        "Cos[a_] - Cos[b_] + r___ :> -2 Sin[(a + b)/2] Sin[(a - b)/2] + r, "
        /* Same-head sum-to-product (hyperbolic). */
        "Sinh[a_] + Sinh[b_] + r___ :> 2 Sinh[(a + b)/2] Cosh[(a - b)/2] + r, "
        "Sinh[a_] - Sinh[b_] + r___ :> 2 Cosh[(a + b)/2] Sinh[(a - b)/2] + r, "
        "Cosh[a_] + Cosh[b_] + r___ :> 2 Cosh[(a + b)/2] Cosh[(a - b)/2] + r, "
        "Cosh[a_] - Cosh[b_] + r___ :> 2 Sinh[(a + b)/2] Sinh[(a - b)/2] + r, "
        /* Mixed Sin/Cos (Sec-vs-Csc numerators): Cos[b] = Sin[Pi/2 - b). */
        "Sin[a_] + Cos[b_] + r___ :> 2 Sin[(a + Pi/2 - b)/2] Cos[(a - Pi/2 + b)/2] + r, "
        "Sin[a_] - Cos[b_] + r___ :> 2 Cos[(a + Pi/2 - b)/2] Sin[(a - Pi/2 + b)/2] + r, "
        "-Sin[a_] + Cos[b_] + r___ :> -2 Cos[(a + Pi/2 - b)/2] Sin[(a - Pi/2 + b)/2] + r, "
        "-Sin[a_] - Cos[b_] + r___ :> -2 Sin[(a + Pi/2 - b)/2] Cos[(a - Pi/2 + b)/2] + r "
        "}");

    symtab_add_builtin("Solve`SolveTrigPair", builtin_solve_trig_pair);
    SymbolDef* def = symtab_get_def("Solve`SolveTrigPair");
    if (def) def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Solve`SolveTrigPair",
        "Solve`SolveTrigPair[lhs == rhs, var]\n"
        "Solve`SolveTrigPair[lhs == rhs, var, dom]\n"
        "\tThe argument-pair reducer used by Solve for generalized two-\n"
        "\targument trig / hyperbolic equations f[A(var)] +/- g[B(var)] ==\n"
        "\tc.  Reciprocal-normalizes to Sin/Cos (Sinh/Cosh), combines over\n"
        "\ta common denominator, factors the numerator into single-argument\n"
        "\tatoms via sum-to-product / reverse-angle-addition identities,\n"
        "\tsolves each atom through the inverse-function specialist, and\n"
        "\tdrops solution families that hit a denominator pole.");
}
