#include "simp.h"
#include "simp_internal.h"
#include "arithmetic.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expand.h"
#include "facpoly.h"
#include "numeric.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "expr.h"
#include "rationalize.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "trigrat.h"
#include "qa.h"
#include "qafactor.h"
#include "simp_log.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif


/* ----------------------------------------------------------------------- */
/* Algebraic-extension reduction: simp_algebraic                           */
/* ----------------------------------------------------------------------- */

/*
 * simp_algebraic rewrites an expression that contains one or more
 * distinct square-root sub-expressions Sqrt[u_i] by treating each
 * Sqrt[u_i] as a generator g_i of the algebraic extension
 *   K(vars)[g_1, ..., g_n] / (g_1^2 - u_1, ..., g_n^2 - u_n).
 * Standard rational arithmetic (Together) followed by reduction modulo
 * the relation ideal and successive rationalisation of the denominator
 * collapses identities that ordinary Together / Cancel cannot see, e.g.
 *
 *   (x/Sqrt[x^2+1] + 1) / ((Sqrt[x^2+1] + x)^2 + 1)  ->  1/(2 + 2 x^2)
 *   (x/(Sqrt[x^2+6] - Sqrt[6]))(1/Sqrt[x^2+6]
 *      - (Sqrt[x^2+6] - Sqrt[6])/x^2)                ->  Sqrt[6]/(x Sqrt[x^2+6])
 *   2/(Sqrt[(x+1)/(1-x)] - 1/Sqrt[(x+1)/(1-x)])      ->  Sqrt[(x+1)/(1-x)] (1-x)/x
 *
 * Algorithm:
 *   1. Walk the expression collecting every distinct surd argument
 *      u_i where Power[u_i, p/q] appears with q != 1. Bail if any q != 2
 *      (cube roots etc.), if any u_i contains an explicit complex
 *      literal, or if more than ALG_MAX_SURDS distinct bases appear.
 *   2. Substitute Sqrt[u_i] -> g_i for fresh distinct generator symbols.
 *      After substitution, the expression is a rational function in
 *      (vars, g_1, ..., g_n).
 *   3. Together  ->  N / D, both polynomials in (vars, g_1, ..., g_n).
 *   4. Reduce both N and D modulo the relation ideal {g_i^2 - u_i}_i
 *      via successive CoefficientList[..., g_i] decomposition. After
 *      one sweep across all generators the polynomial is multilinear
 *      in {g_i}: every g_i appears at degree 0 or 1.
 *   5. For i = 1..n, rationalise the i-th generator out of the
 *      denominator: multiply numerator and denominator by sigma_i(D)
 *      (D with g_i sign-flipped), then reduce again. After each step
 *      g_i has been eliminated from the denominator. The product
 *      D * sigma_i(D) lies in K[g_1, ..., g_{i-1}, g_{i+1}, ..., g_n]
 *      because the g_i terms in (a + b g_i)(a - b g_i) = a^2 - b^2 u_i
 *      cancel.
 *   6. Substitute g_i -> Sqrt[u_i] back, run Together / Cancel for
 *      cleanup, and accept the result iff its complexity score is
 *      strictly lower than the input.
 *
 * Principal-branch concern: the substitution Sqrt[u_i]^2 = u_i is only
 * sound where u_i lies in the principal branch's domain. We accept the
 * Mathematica-style convention (Simplify treats this as an identity on
 * the natural domain where the input is real) but skip when any u_i
 * contains an explicit complex literal (Complex[..,..] or the symbol I)
 * so we never produce a result that swallows a sign-of-imaginary-part
 * change silently.
 */

#define ALG_MAX_SURDS 4

/* Walk e collecting distinct surd bases. The walker enforces:
 *   - every Power[base, p/q] with q != 1 has q == 2,
 *   - distinct bases (by structural equality, expr_eq) accumulate into
 *     bases[0..*n_bases-1] up to max_n,
 *   - returns false on q != 2 or when bases would overflow max_n.
 *
 * Borrowed pointers into `e`. */
static bool alg_collect_sqrt_bases(const Expr* e, const Expr** bases,
                                   size_t* n_bases, size_t max_n) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        const Expr* base = e->data.function.args[0];
        Expr* exp        = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && q != 1) {
            if (q != 2) return false;
            bool seen = false;
            for (size_t i = 0; i < *n_bases; i++) {
                if (expr_eq((Expr*)base, (Expr*)bases[i])) { seen = true; break; }
            }
            if (!seen) {
                if (*n_bases >= max_n) return false;
                bases[(*n_bases)++] = base;
            }
        }
    }
    if (!alg_collect_sqrt_bases(e->data.function.head, bases, n_bases, max_n))
        return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (!alg_collect_sqrt_bases(e->data.function.args[i], bases, n_bases, max_n))
            return false;
    }
    return true;
}

/* Returns true if any sub-expression has head Complex or contains the
 * symbol I. Used to gate simp_algebraic off explicit-complex inputs
 * whose Sqrt[]^2 = arg identity could mask a branch flip. */
bool contains_explicit_complex(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL && e->data.symbol == SYM_I) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Complex) return true;
    if (contains_explicit_complex(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_explicit_complex(e->data.function.args[i])) return true;
    }
    return false;
}

/* For every Power[bases[i], p/2] in e, replace with bases[i]^floor(p/2)
 * * gens[i]^(p mod 2) (computed via floor-division so negative p is
 * handled correctly). Bases that don't appear are passed through. */
static Expr* alg_subst_sqrt_to_gens(const Expr* e, const Expr** bases,
                                    const char** gens, size_t n) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && q == 2) {
            const Expr* base = e->data.function.args[0];
            for (size_t i = 0; i < n; i++) {
                if (expr_eq((Expr*)base, (Expr*)bases[i])) {
                    int64_t m = p / 2;
                    int64_t r = p - 2 * m;
                    if (r < 0) { m -= 1; r += 2; }
                    Expr* base_pow = (m == 0)
                        ? expr_new_integer(1)
                        : eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){expr_copy((Expr*)base), expr_new_integer(m)}, 2));
                    Expr* g_pow = (r == 0)
                        ? expr_new_integer(1)
                        : expr_new_symbol(gens[i]);
                    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                              (Expr*[]){base_pow, g_pow}, 2));
                }
            }
        }
    }
    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = alg_subst_sqrt_to_gens(e->data.function.args[i], bases, gens, n);
    }
    Expr* new_head = alg_subst_sqrt_to_gens(e->data.function.head, bases, gens, n);
    Expr* result = eval_and_free(expr_new_function(new_head, new_args, count));
    free(new_args);
    return result;
}

/* For each generator symbol gens[i], replace with Sqrt[bases[i]]. */
static Expr* alg_subst_gens_to_sqrt(const Expr* e, const char** gens,
                                    const Expr** bases, size_t n) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < n; i++) {
            if (e->data.symbol == gens[i]) {
                return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){expr_copy((Expr*)bases[i]), make_rational(1, 2)}, 2));
            }
        }
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = alg_subst_gens_to_sqrt(e->data.function.args[i], gens, bases, n);
    }
    Expr* new_head = alg_subst_gens_to_sqrt(e->data.function.head, gens, bases, n);
    Expr* result = eval_and_free(expr_new_function(new_head, new_args, count));
    free(new_args);
    return result;
}

/* Replace every occurrence of generator gi_sym with -gi_sym in e. Used
 * to compute sigma_i(den) for rationalisation. */
static Expr* alg_sigma_negate(const Expr* e, const char* gi_sym) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL && e->data.symbol == gi_sym) {
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){expr_new_integer(-1), expr_copy((Expr*)e)}, 2));
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t count = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (count ? count : 1));
    for (size_t i = 0; i < count; i++) {
        new_args[i] = alg_sigma_negate(e->data.function.args[i], gi_sym);
    }
    Expr* new_head = alg_sigma_negate(e->data.function.head, gi_sym);
    Expr* result = eval_and_free(expr_new_function(new_head, new_args, count));
    free(new_args);
    return result;
}

/* Reduce poly modulo gi_sym^2 - u_i: returns A + B*gi where
 *   A = sum_{k even} a_k u^(k/2)
 *   B = sum_{k odd}  a_k u^((k-1)/2)
 * with a_k extracted via CoefficientList[poly, gi_sym]. The caller is
 * expected to have Expand-ed `poly` first. Returns NULL if
 * CoefficientList didn't yield a List. */
static Expr* alg_reduce_one_gen(const Expr* poly, const char* gi_sym,
                                const Expr* ui) {
    Expr* cl_args[2] = { expr_copy((Expr*)poly), expr_new_symbol(gi_sym) };
    Expr* cl_call = expr_new_function(expr_new_symbol("CoefficientList"),
                                      cl_args, 2);
    Expr* coefs = evaluate(cl_call);
    expr_free(cl_call);
    if (!coefs || coefs->type != EXPR_FUNCTION ||
        coefs->data.function.head->type != EXPR_SYMBOL ||
        coefs->data.function.head->data.symbol != SYM_List) {
        if (coefs) expr_free(coefs);
        return NULL;
    }

    size_t n = coefs->data.function.arg_count;
    Expr** evens = (Expr**)malloc(sizeof(Expr*) * (n + 1));
    Expr** odds  = (Expr**)malloc(sizeof(Expr*) * (n + 1));
    size_t ne = 0, no = 0;
    for (size_t k = 0; k < n; k++) {
        Expr* ck = coefs->data.function.args[k];
        int64_t exp_u = (int64_t)(k / 2);
        Expr* upow = (exp_u == 0)
            ? expr_new_integer(1)
            : eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                  (Expr*[]){expr_copy((Expr*)ui), expr_new_integer(exp_u)}, 2));
        Expr* term = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){expr_copy(ck), upow}, 2));
        if ((k & 1) == 0) evens[ne++] = term;
        else              odds[no++]  = term;
    }
    expr_free(coefs);

    Expr* a_sum;
    if (ne == 0)      a_sum = expr_new_integer(0);
    else if (ne == 1) a_sum = evens[0];
    else              a_sum = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                                                              evens, ne));
    Expr* b_sum;
    if (no == 0)      b_sum = expr_new_integer(0);
    else if (no == 1) b_sum = odds[0];
    else              b_sum = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                                                              odds, no));
    free(evens);
    free(odds);

    /* Combine A + B*gi into a single expression. */
    Expr* b_gi = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){b_sum, expr_new_symbol(gi_sym)}, 2));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                  (Expr*[]){a_sum, b_gi}, 2));
}

/* Polynomial-divide `poly` by `u` repeatedly (treating both as
 * polynomials in `var`) until the division has a non-zero remainder.
 * Returns the residual quotient and writes the multiplicity to *k_out
 * so that `poly = u^(*k_out) * residual` modulo non-divisibility.
 *
 * Used so that an implicit u_i^k factor inside den_r (e.g. x^4 hiding
 * (x^2)^2 when u = x^2) can be lifted into Power[g_i, 2k] -- once the
 * u-power is expressed in terms of the generator, polynomial GCD over
 * Q[vars, g_1, ..., g_n] cancels it against any g_i factors carried by
 * the multilinear numerator. */
static Expr* alg_extract_u_power(const Expr* poly, const Expr* u,
                                 const Expr* var, int* k_out) {
    int k = 0;
    Expr* cur = expr_copy((Expr*)poly);
    for (;;) {
        Expr* qa[3] = { expr_copy(cur), expr_copy((Expr*)u), expr_copy((Expr*)var) };
        Expr* qcall = expr_new_function(expr_new_symbol("PolynomialQuotient"), qa, 3);
        Expr* q = evaluate(qcall);
        expr_free(qcall);

        Expr* ra[3] = { expr_copy(cur), expr_copy((Expr*)u), expr_copy((Expr*)var) };
        Expr* rcall = expr_new_function(expr_new_symbol("PolynomialRemainder"), ra, 3);
        Expr* r = evaluate(rcall);
        expr_free(rcall);

        bool zero = (r && r->type == EXPR_INTEGER && r->data.integer == 0);
        if (r) expr_free(r);
        if (!zero) { expr_free(q); break; }
        expr_free(cur);
        cur = q;
        k++;
        /* Defensive cap: prevent runaway when PolynomialQuotient
         * misbehaves (e.g. floating-point coefficients sneaking in). */
        if (k > 100) break;
    }
    *k_out = k;
    return cur;
}

/* Returns true iff u is a polynomial in its own variables -- i.e.,
 * every Power[base, exp] in u has a non-negative integer exp. Rational
 * u (e.g. (x+1)/(1-x)) is rejected so the polynomial-division u-power
 * extraction never tries to divide by a non-polynomial divisor. */
static bool alg_u_is_polynomial(const Expr* u) {
    if (!u) return false;
    if (u->type != EXPR_FUNCTION) return true;   /* leaf is always polynomial */
    if (u->data.function.head &&
        u->data.function.head->type == EXPR_SYMBOL &&
        u->data.function.head->data.symbol == SYM_Power &&
        u->data.function.arg_count == 2) {
        Expr* exp = u->data.function.args[1];
        if (exp->type != EXPR_INTEGER && exp->type != EXPR_BIGINT) return false;
        if (exp->type == EXPR_INTEGER && exp->data.integer < 0) return false;
        if (exp->type == EXPR_BIGINT && mpz_sgn(exp->data.bigint) < 0) return false;
    }
    if (!alg_u_is_polynomial(u->data.function.head)) return false;
    for (size_t i = 0; i < u->data.function.arg_count; i++) {
        if (!alg_u_is_polynomial(u->data.function.args[i])) return false;
    }
    return true;
}

/* Pick the first variable in Variables[u]. Returns NULL when u is
 * variable-free (numeric / constant), in which case alg_extract_u_power
 * is undefined and the caller should skip the u-power-extraction step. */
static Expr* alg_pick_var(const Expr* u) {
    Expr* vars = call_unary_copy("Variables", u);
    if (!vars || vars->type != EXPR_FUNCTION ||
        vars->data.function.head->type != EXPR_SYMBOL ||
        vars->data.function.head->data.symbol != SYM_List ||
        vars->data.function.arg_count == 0) {
        if (vars) expr_free(vars);
        return NULL;
    }
    Expr* v = expr_copy(vars->data.function.args[0]);
    expr_free(vars);
    return v;
}

/* Reduce poly modulo all relations {gi_sym^2 - u_i}_i by sweeping each
 * generator once. The result is multilinear in (g_1, ..., g_n). The
 * input is Expand-ed before each generator pass so CoefficientList sees
 * the canonical polynomial form. Returns NULL on any inner failure. */
static Expr* alg_reduce_all(const Expr* poly, const char** gens,
                            const Expr** us, size_t n) {
    Expr* cur = call_unary_copy("Expand", poly);
    for (size_t i = 0; i < n; i++) {
        Expr* nxt = alg_reduce_one_gen(cur, gens[i], us[i]);
        expr_free(cur);
        if (!nxt) return NULL;
        cur = call_unary_owned("Expand", nxt);
    }
    return cur;
}

static Expr* simp_algebraic_impl(const Expr* e) {
    bool dbg = simp_debug_enabled();
    clock_t t0 = dbg ? clock() : 0;

    /* Cheap precondition: nothing to do without a half-integer Power. */
    if (!has_non_integer_power(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("Algebraic", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    /* Early denester short-circuit. For nested-radical inputs (e.g.
     * Sqrt[A + b Sqrt[C]] or 1/Sqrt[...]^k), the half-sum denester is
     * O(small) and the algebraic-tower Together-with-Extension below is
     * O(big). When the denester strictly wins on complexity, taking that
     * win here saves the 200+ ms it costs to compute and then reject
     * the Together[expr, Extension -> α] form for the same input.
     *
     * The complexity gate is strict (<) so we don't churn on inputs the
     * denester touches without reducing -- those still fall through to
     * the extension path which may legitimately help. */
    {
        Expr* denested = simp_denest_sqrt(e, NULL);
        if (denested && simp_default_complexity(denested)
                          < simp_default_complexity(e)) {
            if (dbg) simp_debug_log("Algebraic", e, denested,
                                    simp_debug_elapsed_ms(t0));
            return denested;
        }
        if (denested) expr_free(denested);
    }

    /* Phase G9 (cube-root and higher): when the input has exactly one
     * rational-base radical generator, route through
     * `Together[expr, Extension -> α]`.  Together's extension path uses
     * the qaupoly substrate (qaupoly_gcd / qaupoly_divrem), which handles
     * any q ≥ 2 natively — no Sqrt-only special case.  This is the
     * shortcut that lets Simplify simplify expressions involving
     * `Power[c, p/q]` for q > 2, which the older multi-Sqrt path
     * (alg_sigma_negate sign-flip rationalisation) below cannot.
     *
     * Multi-generator cases (n ≥ 2) fall through to the Sqrt-only
     * path because that path's general rationalisation by sign-flip
     * conjugates handles n ≥ 2 over Q(Sqrt[u_i]) directly.
     *
     * Same Layer-0 prefilter as builtin_simplify's alg_top branch:
     * inputs with nested radicals and no free polynomial variable end
     * up in the multi-generator (n ≥ 2) branch anyway because the
     * nested generator counts as one — autodetect's primitive-element
     * compositum work is ~tens-of-ms wasted before falling through. */
    if (!contains_explicit_complex(e)
        && !expr_has_nested_radical_radicand(e)) {
        QATower* qa_t = extension_autodetect(e);
        if (qa_t && qa_t->n == 1) {
            Expr* alpha = expr_copy(qa_t->alpha_renders[0]);
            Expr* tog = expr_new_function(
                expr_new_symbol("Together"),
                (Expr*[]){
                    expr_copy((Expr*)e),
                    expr_new_function(expr_new_symbol(SYM_Rule),
                        (Expr*[]){expr_new_symbol(SYM_Extension), alpha}, 2)
                }, 2);
            Expr* qa_out = evaluate(tog);
            expr_free(tog);
            qa_tower_free(qa_t);

            /* Accept only when the qaupoly path produced something at
             * least as compact as the input.  The Together-with-Extension
             * path occasionally factors a denominator unnecessarily —
             * letting simp_search's later round-loop pick the better
             * form is the safe move. */
            if (qa_out && simp_default_complexity(qa_out)
                            <= simp_default_complexity(e)) {
                if (dbg) simp_debug_log("Algebraic", e, qa_out,
                                        simp_debug_elapsed_ms(t0));
                return qa_out;
            }
            if (qa_out) expr_free(qa_out);
        } else if (qa_t) {
            qa_tower_free(qa_t);
        }
    }

    const Expr* bases[ALG_MAX_SURDS];
    size_t n = 0;
    if (!alg_collect_sqrt_bases(e, bases, &n, ALG_MAX_SURDS) || n == 0) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("Algebraic", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }
    /* Each surd's argument must itself be surd-free and contain no
     * explicit complex literals. */
    for (size_t i = 0; i < n; i++) {
        if (has_non_integer_power(bases[i]) ||
            contains_explicit_complex(bases[i])) {
            Expr* out = expr_copy((Expr*)e);
            if (dbg) simp_debug_log("Algebraic", e, out,
                                    simp_debug_elapsed_ms(t0));
            return out;
        }
    }
    if (contains_explicit_complex(e)) {
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("Algebraic", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    /* Allocate a fresh interned generator symbol per surd. The names
     * are $-prefixed so they won't collide with user symbols. */
    const char* gens[ALG_MAX_SURDS];
    static const char* gen_names[ALG_MAX_SURDS] = {
        "$pc_alggen0$", "$pc_alggen1$", "$pc_alggen2$", "$pc_alggen3$"
    };
    for (size_t i = 0; i < n; i++) gens[i] = intern_symbol(gen_names[i]);

    /* Step 2-3: substitute and Together. */
    Expr* sub = alg_subst_sqrt_to_gens(e, bases, gens, n);
    Expr* tg  = call_unary_owned("Together", sub);
    Expr* num = call_unary_copy("Numerator",   tg);
    Expr* den = call_unary_copy("Denominator", tg);
    expr_free(tg);

    /* Step 4: reduce both modulo the relation ideal. */
    Expr* num_r = alg_reduce_all(num, gens, bases, n);
    Expr* den_r = alg_reduce_all(den, gens, bases, n);
    expr_free(num); expr_free(den);
    if (!num_r || !den_r) {
        if (num_r) expr_free(num_r);
        if (den_r) expr_free(den_r);
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("Algebraic", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    /* Step 5: rationalise each generator out of the denominator in turn. */
    for (size_t i = 0; i < n; i++) {
        Expr* sig = alg_sigma_negate(den_r, gens[i]);
        Expr* num_mul = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){num_r, expr_copy(sig)}, 2));
        Expr* den_mul = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){den_r, sig}, 2));
        Expr* num_next = alg_reduce_all(num_mul, gens, bases, n);
        Expr* den_next = alg_reduce_all(den_mul, gens, bases, n);
        expr_free(num_mul); expr_free(den_mul);
        if (!num_next || !den_next) {
            if (num_next) expr_free(num_next);
            if (den_next) expr_free(den_next);
            Expr* out = expr_copy((Expr*)e);
            if (dbg) simp_debug_log("Algebraic", e, out,
                                    simp_debug_elapsed_ms(t0));
            return out;
        }
        num_r = num_next;
        den_r = den_next;
    }

    /* Step 6 (pre): pull each implicit u_i^k factor out of num_r and
     * den_r, replacing it with g_i^(2k) so that Cancel over
     * Q[vars, g_1, ..., g_n] sees the cancellation between the
     * multilinear g_i factor in the numerator and an implicit u_i^k
     * factor in the denominator. Without this step, x^4 in the
     * denominator and Sqrt[x^2] in the numerator look like coprime
     * polynomial atoms to Cancel even though x^4 = u^2 = g^4 modulo
     * the algebraic relation g^2 = u. */
    for (size_t i = 0; i < n; i++) {
        if (!alg_u_is_polynomial(bases[i])) continue;  /* rational u: skip */
        Expr* var = alg_pick_var(bases[i]);
        if (!var) continue;     /* numeric u_i: no polynomial division */

        int kn = 0, kd = 0;
        Expr* num_resid = alg_extract_u_power(num_r, bases[i], var, &kn);
        Expr* den_resid = alg_extract_u_power(den_r, bases[i], var, &kd);
        expr_free(var);

        if (kn == 0 && kd == 0) {
            expr_free(num_resid); expr_free(den_resid);
            continue;
        }
        /* num_r = num_resid * Power[g_i, 2*kn]
         * den_r = den_resid * Power[g_i, 2*kd] */
        if (kn > 0) {
            Expr* g_pow = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){expr_new_symbol(gens[i]), expr_new_integer(2 * kn)}, 2));
            expr_free(num_r);
            num_r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){num_resid, g_pow}, 2));
        } else {
            expr_free(num_r);
            num_r = num_resid;
        }
        if (kd > 0) {
            Expr* g_pow = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){expr_new_symbol(gens[i]), expr_new_integer(2 * kd)}, 2));
            expr_free(den_r);
            den_r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                      (Expr*[]){den_resid, g_pow}, 2));
        } else {
            expr_free(den_r);
            den_r = den_resid;
        }
    }

    /* Step 6: assemble num_r / den_r, substitute generators back, clean.
     *
     * Apply Factor to the polynomial-in-(vars, g_1..g_n) numerator and
     * denominator before substituting g_i -> Sqrt[u_i]. Without this
     * step, Cancel sees expanded polynomials whose common (u_i)^k
     * factors share denominators with Sqrt[u_i]^(2k); Factor exposes
     * the (u_i)^k structure so Cancel can combine
     * Power[u_i, k] * Power[u_i, 1/2] / Power[u_i, j]
     * into a single Power[u_i, ...] term. */
    Expr* num_factored = call_unary_owned("Factor", num_r);
    Expr* den_factored = call_unary_owned("Factor", den_r);

    Expr* den_inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                  (Expr*[]){den_factored, expr_new_integer(-1)}, 2));
    Expr* quot = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){num_factored, den_inv}, 2));

    Expr* with_sqrt = alg_subst_gens_to_sqrt(quot, gens, bases, n);
    expr_free(quot);
    Expr* result = call_unary_owned("Cancel", with_sqrt);

    /* Complexity gate: accept any form whose complexity score is no
     * greater than the input. The strict ">=" rejection used by the
     * simp_search round loop is too tight here because rationalisation
     * often hits a tied score (e.g. 1/(Sqrt[a]+Sqrt[b]) trades the
     * Power[..,-1] head for a single Times[-1, ...] term while keeping
     * two Sqrt leaves -- equal complexity but the rationalised form is
     * the conventionally preferred shape). simp_search's later round
     * loop will still pick the strictly-better form when one exists. */
    if (simp_default_complexity(result) > simp_default_complexity(e)) {
        expr_free(result);
        Expr* out = expr_copy((Expr*)e);
        if (dbg) simp_debug_log("Algebraic", e, out,
                                simp_debug_elapsed_ms(t0));
        return out;
    }

    if (dbg) simp_debug_log("Algebraic", e, result,
                            simp_debug_elapsed_ms(t0));
    return result;
}

Expr* simp_algebraic(const Expr* e) {
    return simp_memo_wrap(e, "$Algebraic", simp_algebraic_impl);
}

