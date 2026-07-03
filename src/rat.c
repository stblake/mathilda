#include "rat.h"
#include "common.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "poly.h"
#include "flint_bridge.h"
#include "expand.h"
#include "rationalize.h"
#include "sym_names.h"
#include "options.h"
#include "qafactor.h"
#include "facpoly.h"
#include "core.h"
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/* Together's recursive walker.  Cancel and Together both share the
 * algebraic-generator (radical) pass exported from poly.c, but
 * together_recursive is internal because it embeds the Plus-merging
 * logic specific to combining over a common denominator. */
static Expr* together_recursive(Expr* e);

static Expr* negate_expr(Expr* e) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), expr_copy(e)}, 2));
}

static bool is_superficially_negative(Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer < 0;
    if (e->type == EXPR_REAL) return e->data.real < 0.0;
    int64_t n, d;
    if (is_rational(e, &n, &d)) return n < 0;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL && e->data.function.head->data.symbol == SYM_Times) {
        if (e->data.function.arg_count > 0) {
            Expr* first = e->data.function.args[0];
            if (first->type == EXPR_INTEGER) return first->data.integer < 0;
            if (first->type == EXPR_REAL) return first->data.real < 0.0;
            if (is_rational(first, &n, &d)) return n < 0;
        }
    }
    return false;
}

/* den_has_negative_lead: true when the post-cancellation denominator
 * carries a globally negative sign that should be pushed up into the
 * numerator. PolynomialGCD's sign convention can leave both num and den
 * with negative leading coefficients (e.g. (x-3y)/(x+y) coming back as
 * (-x+3y)/(-x-y)); the leaf-count tiebreak in Simplify then picks the
 * uglier form. We treat two shapes as "negative lead":
 *   1. is_superficially_negative — a literal negative number, or
 *      a Times whose leading factor is negative.
 *   2. A Plus whose every term is superficially negative — i.e. the
 *      whole sum can be factored as -1 * (positive-lead sum).
 * Mixed-sign Plus (e.g. -x + y) is left alone since both orderings are
 * equally complex and either choice is arbitrary. */
static bool den_has_negative_lead(Expr* e) {
    if (!e) return false;
    if (is_superficially_negative(e)) return true;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Plus &&
        e->data.function.arg_count > 0) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (!is_superficially_negative(e->data.function.args[i])) return false;
        }
        return true;
    }
    return false;
}

static void extract_num_den(Expr* expr, Expr** num_out, Expr** den_out) {
    int64_t n, d;
    if (is_rational(expr, &n, &d)) {
        *num_out = expr_new_integer(n);
        *den_out = expr_new_integer(d);
        return;
    }

    Expr* re; Expr* im;
    if (is_complex(expr, &re, &im)) {
        int64_t rn = 0, rd = 1, in = 0, id = 1;
        is_rational(re, &rn, &rd);
        if (re->type == EXPR_INTEGER) { rn = re->data.integer; rd = 1; }
        is_rational(im, &in, &id);
        if (im->type == EXPR_INTEGER) { in = im->data.integer; id = 1; }

        int64_t common_den = lcm(rd, id);
        
        Expr* d_expr = expr_new_integer(common_den);
        *den_out = d_expr;
        *num_out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_copy(expr), expr_copy(d_expr)}, 2));
        return;
    }

    if (expr->type == EXPR_FUNCTION && expr->data.function.head->type == EXPR_SYMBOL && (expr->data.function.head->data.symbol == SYM_Power || expr->data.function.head->data.symbol == SYM_Exp)) {
        bool is_exp = expr->data.function.head->data.symbol == SYM_Exp;
        Expr* base = is_exp ? expr_new_symbol(SYM_E) : expr->data.function.args[0];
        Expr* exp = is_exp ? expr->data.function.args[0] : expr->data.function.args[1];

        if (exp->type == EXPR_FUNCTION && exp->data.function.head->type == EXPR_SYMBOL && exp->data.function.head->data.symbol == SYM_Plus) {
            size_t count = exp->data.function.arg_count;
            Expr** num_args = malloc(sizeof(Expr*) * count);
            Expr** den_args = malloc(sizeof(Expr*) * count);
            size_t n_c = 0, d_c = 0;
            for (size_t i = 0; i < count; i++) {
                Expr* arg = exp->data.function.args[i];
                if (is_superficially_negative(arg)) {
                    den_args[d_c++] = negate_expr(arg);
                } else {
                    num_args[n_c++] = expr_copy(arg);
                }
            }
            if (n_c == 0) *num_out = expr_new_integer(1);
            else if (n_c == 1) *num_out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(base), num_args[0]}, 2));
            else {
                Expr* p = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), num_args, n_c));
                *num_out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(base), p}, 2));
            }
            if (d_c == 0) *den_out = expr_new_integer(1);
            else if (d_c == 1) *den_out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(base), den_args[0]}, 2));
            else {
                Expr* p = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), den_args, d_c));
                *den_out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(base), p}, 2));
            }
            free(num_args);
            free(den_args);
            if (is_exp) expr_free(base);
            return;
        } else {
            if (is_superficially_negative(exp)) {
                *num_out = expr_new_integer(1);
                *den_out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(base), negate_expr(exp)}, 2));
            } else {
                if (is_exp) {
                    *num_out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(base), expr_copy(exp)}, 2));
                } else {
                    *num_out = expr_copy(expr);
                }
                *den_out = expr_new_integer(1);
            }
            if (is_exp) expr_free(base);
            return;
        }
    }

    if (expr->type == EXPR_FUNCTION && expr->data.function.head->type == EXPR_SYMBOL && expr->data.function.head->data.symbol == SYM_Times) {
        size_t count = expr->data.function.arg_count;
        Expr** n_args = malloc(sizeof(Expr*) * count);
        Expr** d_args = malloc(sizeof(Expr*) * count);
        size_t n_c = 0, d_c = 0;
        for (size_t i = 0; i < count; i++) {
            Expr* n_out; Expr* d_out;
            extract_num_den(expr->data.function.args[i], &n_out, &d_out);
            if (!(n_out->type == EXPR_INTEGER && n_out->data.integer == 1)) n_args[n_c++] = n_out;
            else expr_free(n_out);
            if (!(d_out->type == EXPR_INTEGER && d_out->data.integer == 1)) d_args[d_c++] = d_out;
            else expr_free(d_out);
        }
        if (n_c == 0) *num_out = expr_new_integer(1);
        else if (n_c == 1) *num_out = n_args[0];
        else *num_out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), n_args, n_c));

        if (d_c == 0) *den_out = expr_new_integer(1);
        else if (d_c == 1) *den_out = d_args[0];
        else *den_out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), d_args, d_c));

        free(n_args); free(d_args);
        return;
    }

    *num_out = expr_copy(expr);
    *den_out = expr_new_integer(1);
}

Expr* builtin_numerator(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* n; Expr* d;
    extract_num_den(res->data.function.args[0], &n, &d);
    expr_free(d);
    return n;
}

Expr* builtin_denominator(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* n; Expr* d;
    extract_num_den(res->data.function.args[0], &n, &d);
    expr_free(n);
    return d;
}

/* Strict variant of cancel_exact_div: returns NULL when exact_poly_div
 * cannot perform the division in Q[vars] (i.e. the divisor doesn't
 * actually divide the dividend in the working multivariate polynomial
 * ring). Used by together_recursive's Plus combine to detect when the
 * iterative PolynomialLCM produced a result that is only valid in some
 * algebraic-extension ring -- in that case the would-be quotient is a
 * symbolic Times[..., Power[dens[i], -1]] which, if propagated, leaves
 * a Power[Plus[...], -1] inside the combined numerator. The downstream
 * cancel_recursive -> PolynomialGCD path then re-enters multivariate
 * Euclid on a rational expression and explodes (case-13 hang).
 *
 * Soundness: returning NULL is always safe -- the caller falls back to
 * leaving the Plus uncombined, which is correctness-preserving. The
 * later simp_algebraic seed in simp_search handles the algebraic
 * extension properly when the multi-radical structure is amenable. */
static Expr* cancel_exact_div_strict(Expr* num, Expr* den) {
    if (is_zero_poly(num)) return expr_new_integer(0);
    if (den->type == EXPR_INTEGER && den->data.integer == 1) return expr_expand(num);

    Expr* exp_num = expr_expand(num);
    Expr* exp_den = expr_expand(den);

    size_t v_count = 0, v_cap = 16;
    Expr** vars = malloc(sizeof(Expr*) * v_cap);
    collect_variables(exp_num, &vars, &v_count, &v_cap);
    collect_variables(exp_den, &vars, &v_count, &v_cap);
    if (v_count > 0) qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);

    Expr* res = exact_poly_div(exp_num, exp_den, vars, v_count);

    for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
    free(vars);

    expr_free(exp_num);
    expr_free(exp_den);
    return res;  /* NULL if not exactly divisible */
}

static Expr* cancel_exact_div_wrapper(Expr* num, Expr* den) {
    if (is_zero_poly(num)) return expr_new_integer(0);
    if (den->type == EXPR_INTEGER && den->data.integer == 1) return expr_expand(num);

    Expr* exp_num = expr_expand(num);
    Expr* exp_den = expr_expand(den);

    size_t v_count = 0, v_cap = 16;
    Expr** vars = malloc(sizeof(Expr*) * v_cap);
    collect_variables(exp_num, &vars, &v_count, &v_cap);
    collect_variables(exp_den, &vars, &v_count, &v_cap);
    if (v_count > 0) qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);

    Expr* res = exact_poly_div(exp_num, exp_den, vars, v_count);

    for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
    free(vars);

#ifdef USE_FLINT
    /* The classical exact_poly_div treats an algebraic generator (Sqrt[k],
     * zeta_n, …) as a free variable and so cannot divide e.g. x^2-2 by x-Sqrt[2]
     * (which needs Sqrt[2]^2 = 2). When it fails, try exact division in the
     * actual extension field; only an exact quotient is returned. */
    if (!res) {
        Expr* fq = flint_extension_divexact(exp_num, exp_den);
        if (fq) {
            expr_free(exp_num);
            expr_free(exp_den);
            return fq;
        }
    }
#endif

    if (res) {
        expr_free(exp_num);
        expr_free(exp_den);
        return res;
    } else {
        Expr* t = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){exp_den, expr_new_integer(-1)}, 2));
        Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){exp_num, t}, 2));
        return r;
    }
}

/* Walk e looking for any Power[X, k] subterm where k is a negative
 * integer (or rational) AND X is itself a compound expression (Plus,
 * Times with non-numeric leaves, etc.) -- i.e. a non-trivial rational
 * sub-expression embedded inside e. Returns true if any is found.
 *
 * Used as a soundness gate before PolynomialGCD: when num or den has
 * such a sub-expression, the inputs aren't pure polynomials in the
 * working ring, and the multivariate GCD can run away (case-13 hang).
 * In that case the caller should leave the input uncancelled. */
static bool has_embedded_rational_subterm(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        bool neg_exp = false;
        if (exp->type == EXPR_INTEGER && exp->data.integer < 0) neg_exp = true;
        else {
            int64_t pn, qn;
            if (is_rational(exp, &pn, &qn) && pn < 0) neg_exp = true;
        }
        if (neg_exp) {
            /* Numeric base is fine -- represents a rational literal like
             * 1/2 = Power[2,-1]. Compound base means embedded rational. */
            if (base->type == EXPR_FUNCTION) return true;
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (has_embedded_rational_subterm(e->data.function.args[i])) return true;
    }
    return false;
}

/* Build a flat list of the multiplicative factors of `e`. For Times[a, b, ...]
 * the result is [a, b, ...]; for any other expression it is [e]. Each entry
 * is an expr_copy() of the source factor; the caller owns the entries and
 * the array.
 *
 * Used by the symbolic content extractor below to compute the structural
 * factor multiset of each Plus summand.                                      */
static void rat_factor_list_copy(const Expr* e, Expr*** out, size_t* n_out) {
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Times) {
        size_t n = e->data.function.arg_count;
        Expr** arr = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
        for (size_t i = 0; i < n; i++) arr[i] = expr_copy(e->data.function.args[i]);
        *out = arr;
        *n_out = n;
    } else {
        Expr** arr = malloc(sizeof(Expr*));
        arr[0] = expr_copy((Expr*)e);
        *out = arr;
        *n_out = 1;
    }
}

/* Filter a factor list down to *algebraic-constant atoms* — concretely,
 * Power[integer_literal, rational_non_integer_exponent] such as Sqrt[2],
 * Sqrt[3], CubeRoot[5], 2^(2/3).  These are exactly the factors that the
 * existing my_number_gcd / poly_content path cannot detect (it treats
 * Sqrt[2] as having integer content 1, blinding it to the structural
 * match in the coefficient ring).
 *
 * Everything else is dropped:
 *   - Integer / bigint / real literals → handled by the integer-content
 *     recursion.
 *   - Symbols and Power[symbol, integer] → polynomial-like, the standard
 *     multivariate PolynomialGCD picks these up.
 *   - Compound Power bases (Sqrt[1+x], etc.) → leaving these in would
 *     cause Cancel to rearrange intermediates that the integration
 *     dispatcher pattern-matches against (CRC-corpus DIFF NONZERO
 *     regressions in the Sqrt[(1+x)(2+3x)] family).
 *
 * The narrowness is the point: this pass exists *only* to close the
 * Cancel[Sqrt[2] / (Sqrt[2] + Sqrt[2] x^4)] gap, not to second-guess the
 * downstream polynomial path.                                            */
static void rat_strip_numeric_factors(Expr** arr, size_t* n_io) {
    size_t w = 0;
    for (size_t i = 0; i < *n_io; i++) {
        bool keep = false;
        Expr* f = arr[i];
        if (f->type == EXPR_FUNCTION && f->data.function.head &&
            f->data.function.head->type == EXPR_SYMBOL &&
            f->data.function.head->data.symbol == SYM_Power &&
            f->data.function.arg_count == 2) {
            Expr* base = f->data.function.args[0];
            Expr* exp  = f->data.function.args[1];
            bool base_is_literal_int =
                (base->type == EXPR_INTEGER || base->type == EXPR_BIGINT);
            int64_t p, q;
            bool exp_is_non_integer_rational =
                is_rational(exp, &p, &q) && q != 1;
            if (base_is_literal_int && exp_is_non_integer_rational) {
                keep = true;
            }
        }
        if (keep) {
            arr[w++] = f;
        } else {
            expr_free(f);
        }
    }
    *n_io = w;
}

/* Greedy multiset intersection of two factor lists, using expr_eq for
 * structural equality.  Modifies *a_io* in place to hold the intersection
 * (kept entries are the originals from a, unmatched a-entries are freed).
 * Everything in b is freed and *b_io* is set to NULL/0.                   */
static void rat_factor_intersect(Expr*** a_io, size_t* na_io,
                                 Expr*** b_io, size_t* nb_io) {
    Expr** a = *a_io;
    size_t na = *na_io;
    Expr** b = *b_io;
    size_t nb = *nb_io;
    bool* matched_b = nb ? calloc(nb, sizeof(bool)) : NULL;
    size_t w = 0;
    for (size_t i = 0; i < na; i++) {
        bool found = false;
        for (size_t j = 0; j < nb; j++) {
            if (!matched_b[j] && expr_eq(a[i], b[j])) {
                matched_b[j] = true;
                found = true;
                break;
            }
        }
        if (found) {
            a[w++] = a[i];
        } else {
            expr_free(a[i]);
        }
    }
    *na_io = w;
    for (size_t j = 0; j < nb; j++) expr_free(b[j]);
    free(b);
    free(matched_b);
    *b_io = NULL;
    *nb_io = 0;
}

/* Compute the symbolic multiplicative content of `e` — the product of
 * non-numeric factors common to every summand when `e` is a Plus, or the
 * non-numeric factor list of `e` itself when it is not.
 *
 * Examples
 *   Sqrt[2] + Sqrt[2]*x^4          -> Sqrt[2]
 *   Sqrt[2]                        -> Sqrt[2]
 *   2 + 3 x                        -> 1                 (numeric content only)
 *   Sqrt[2] x + Sqrt[2] y          -> Sqrt[2]
 *   Sqrt[2] x + Sqrt[3] y          -> 1                 (no common atom)
 *
 * Numeric content (integer GCD across coefficients) is intentionally NOT
 * extracted here; the existing my_number_gcd / poly_content path handles
 * that branch.  The caller owns the returned expression.                  */
static Expr* rat_symbolic_content(const Expr* e) {
    bool is_plus = (e->type == EXPR_FUNCTION &&
                    e->data.function.head &&
                    e->data.function.head->type == EXPR_SYMBOL &&
                    e->data.function.head->data.symbol == SYM_Plus);

    Expr** running = NULL;
    size_t n_running = 0;

    if (!is_plus) {
        rat_factor_list_copy(e, &running, &n_running);
        rat_strip_numeric_factors(running, &n_running);
    } else {
        size_t ns = e->data.function.arg_count;
        if (ns == 0) return expr_new_integer(1);

        rat_factor_list_copy(e->data.function.args[0], &running, &n_running);
        rat_strip_numeric_factors(running, &n_running);

        for (size_t i = 1; i < ns && n_running > 0; i++) {
            Expr** other = NULL;
            size_t n_other = 0;
            rat_factor_list_copy(e->data.function.args[i], &other, &n_other);
            rat_strip_numeric_factors(other, &n_other);
            rat_factor_intersect(&running, &n_running, &other, &n_other);
        }
    }

    if (n_running == 0) {
        free(running);
        return expr_new_integer(1);
    }
    if (n_running == 1) {
        Expr* r = running[0];
        free(running);
        return r;
    }
    /* expr_new_function copies the `running` array (memcpy), adopting the
     * element pointers but not the buffer — free our buffer, as the
     * n_running == 0 / == 1 branches above do. */
    Expr* r = expr_new_function(expr_new_symbol(SYM_Times), running, n_running);
    free(running);
    return r;
}

/* Distribute the divisor through a Plus so each summand is divided
 * independently — this lets the Times-evaluator collapse common factors
 * like Sqrt[2] * Power[Sqrt[2], -1] -> 1 *inside* every summand rather
 * than wrapping the whole Plus in an opaque Power[..., -1] that the
 * downstream PolynomialGCD step cannot see through.  Non-Plus inputs are
 * just multiplied by the inverse.                                         */
static Expr* rat_div_distribute(Expr* e, const Expr* divisor) {
    if (divisor->type == EXPR_INTEGER && divisor->data.integer == 1) {
        return expr_copy(e);
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Plus) {
        size_t n = e->data.function.arg_count;
        Expr** new_args = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
        for (size_t i = 0; i < n; i++) {
            Expr* inv = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Power),
                (Expr*[]){expr_copy((Expr*)divisor), expr_new_integer(-1)}, 2));
            new_args[i] = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Times),
                (Expr*[]){expr_copy(e->data.function.args[i]), inv}, 2));
        }
        Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), new_args, n));
        free(new_args);  /* expr_new_function copies the buffer */
        return r;
    }
    Expr* inv = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Power),
        (Expr*[]){expr_copy((Expr*)divisor), expr_new_integer(-1)}, 2));
    return eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Times),
        (Expr*[]){expr_copy(e), inv}, 2));
}

/* Pre-cancellation pass: factor out symbolic atoms (Sqrt[2], etc.) that
 * appear in every summand of num *and* every summand of den, then divide
 * them out.  Closes the gap left by my_number_gcd, which only sees
 * integer coefficients — without this pass, Cancel[Sqrt[2] /
 * (Sqrt[2] + Sqrt[2] x^4)] is rejected because PolynomialGCD([Sqrt[2]],
 * [Sqrt[2], Sqrt[2] x^4]) returns 1 (the structural Sqrt[2] match in the
 * coefficient ring is invisible to the integer-content recursion).
 *
 * On entry *num and *den are owned by the caller.  On exit, if a non-
 * trivial common factor was extracted, the originals are freed and
 * replaced with their divided forms (and downstream PolynomialGCD runs
 * on the smaller, primitive expressions).                                */
static void rat_strip_symbolic_common(Expr** num_io, Expr** den_io) {
    Expr* num_sym = rat_symbolic_content(*num_io);
    Expr* den_sym = rat_symbolic_content(*den_io);

    Expr** num_factors = NULL;
    size_t num_nf = 0;
    rat_factor_list_copy(num_sym, &num_factors, &num_nf);
    rat_strip_numeric_factors(num_factors, &num_nf);

    Expr** den_factors = NULL;
    size_t den_nf = 0;
    rat_factor_list_copy(den_sym, &den_factors, &den_nf);
    rat_strip_numeric_factors(den_factors, &den_nf);

    rat_factor_intersect(&num_factors, &num_nf, &den_factors, &den_nf);

    expr_free(num_sym);
    expr_free(den_sym);

    if (num_nf == 0) {
        free(num_factors);
        return;
    }

    Expr* common;
    if (num_nf == 1) {
        common = num_factors[0];
        free(num_factors);
    } else {
        /* expr_new_function copies the buffer (adopting the elements); free
         * our `num_factors` buffer as the num_nf == 1 branch does. */
        common = expr_new_function(expr_new_symbol(SYM_Times), num_factors, num_nf);
        free(num_factors);
    }

    Expr* new_num = rat_div_distribute(*num_io, common);
    Expr* new_den = rat_div_distribute(*den_io, common);
    expr_free(common);
    expr_free(*num_io);
    expr_free(*den_io);
    *num_io = new_num;
    *den_io = new_den;
}

/* Counter incremented around cancel_recursive's PolynomialGCD call.
 * Consumed by builtin_times to disable the eps=+1 Sqrt-coefficient
 * canonicalisation while we're inside a PolynomialGCD pass: that pass
 * multiplies polynomial coefficients many times, and canonicalising
 * c * Sqrt[a/b] -> c' * Sqrt[a'/b'] at each multiplication creates
 * moving-target coefficients that PolynomialGCD's Euclidean iteration
 * cannot drive to zero -- producing the QR symbolic hang on the
 * rank-deficient matrices once the Mathematica-style canonicalisation
 * is enabled.  Leaving the form alone inside PolynomialGCD keeps the
 * Euclidean run stable while still applying the canonical form at
 * user-facing evaluation. */
int cancel_recursive_inside_gcd = 0;

static Expr* cancel_recursive(Expr* e) {
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* head = e->data.function.head->data.symbol;
        if (head == SYM_List || head == SYM_Plus ||
            head == SYM_Equal || head == SYM_Less ||
            head == SYM_LessEqual || head == SYM_Greater ||
            head == SYM_GreaterEqual || head == SYM_And ||
            head == SYM_Or || head == SYM_Not) {

            size_t count = e->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = cancel_recursive(e->data.function.args[i]);
            }
            Expr* ret = eval_and_free(expr_new_function(expr_copy(e->data.function.head), args, count));
            free(args);
            return ret;
        }
    }

    Expr* num; Expr* den;
    extract_num_den(e, &num, &den);

    if (den->type == EXPR_INTEGER && den->data.integer == 1) {
        expr_free(den);
        return num;
    }

    rat_strip_symbolic_common(&num, &den);

    if (den->type == EXPR_INTEGER && den->data.integer == 1) {
        expr_free(den);
        return num;
    }

    /* Soundness gate: if num or den still carries an embedded rational
     * (a Power[Plus/Times, negative]), this isn't a clean polynomial
     * fraction. PolynomialGCD on it would re-enter multivariate Euclid
     * on a rational input and may explode (case-13/14/20 hang).
     * Leave the input uncancelled -- subsequent Simplify passes
     * (simp_algebraic, etc.) handle the algebraic structure soundly. */
    if (has_embedded_rational_subterm(num) || has_embedded_rational_subterm(den)) {
        expr_free(num);
        expr_free(den);
        return expr_copy(e);
    }

    cancel_recursive_inside_gcd++;
    Expr* g = NULL;
#ifdef USE_FLINT
    /* Fast path: the numerator/denominator here are ordinary polynomials
     * over Q (rational subterms were rejected by the soundness gate above).
     * FLINT's packed fmpq_mpoly GCD is dramatically faster than the generic
     * Expr subresultant Euclid the PolynomialGCD builtin falls back to on
     * multivariate input — the dominant cost of Cancel/Together/Simplify on
     * multi-generator algebraic fractions (e.g. simp_algebraic's Sqrt->g_i
     * substituted polynomials). The result is normalised to the same
     * primitive-integer, positive-leading associate the builtin returns, so
     * the surrounding cancellation and sign logic is unaffected. NULL
     * (non-polynomial / numeric input) falls through to the builtin. */
    g = flint_multivariate_gcd_normalized(num, den);
    if (g) g = eval_and_free(g);
#endif
    if (!g)
        g = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialGCD), (Expr*[]){expr_copy(num), expr_copy(den)}, 2));
    cancel_recursive_inside_gcd--;
    
    Expr* new_num = cancel_exact_div_wrapper(num, g);
    Expr* new_den = cancel_exact_div_wrapper(den, g);
    
    if (new_num && new_den && den_has_negative_lead(new_den)) {
        Expr* t1 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), new_num}, 2));
        Expr* t2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), new_den}, 2));
        new_num = expr_expand(t1);
        new_den = expr_expand(t2);
        expr_free(t1);
        expr_free(t2);
    }
    
    Expr* res;
    if (new_den && new_den->type == EXPR_INTEGER && new_den->data.integer == 1) {
        res = new_num;
        expr_free(new_den);
    } else if (new_num && new_den) {
        Expr* inv_den = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){new_den, expr_new_integer(-1)}, 2));
        res = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){new_num, inv_den}, 2));
    } else {
        if (new_num) expr_free(new_num);
        if (new_den) expr_free(new_den);
        res = expr_copy(e);
    }
    
    expr_free(g);
    expr_free(num);
    expr_free(den);
    
    return res;
}

/* Forward declarations for the Extension-aware combine path. */
static Expr* cancel_with_extension(const Expr* arg, const Expr* alpha);

/* together_recursive_ext: extension-aware Plus combining over Q(α).
 *
 * Mirrors together_recursive but threads `Extension -> α` through the
 * underlying PolynomialLCM and PolynomialQuotient calls so the inner
 * polynomial arithmetic stays in the Q(α)[x] substrate (qaupoly_*)
 * rather than dropping into the multivariate Q[α, x] subresultant-PRS
 * path, which is exponentially slow on Sqrt[α]-laden coefficients.
 * Re-enabled now that PolynomialQuotient/Remainder accept Extension. */
static Expr* together_recursive_ext(Expr* e, const Expr* alpha) {
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* head = e->data.function.head->data.symbol;

        if (head == SYM_Plus) {
            size_t count = e->data.function.arg_count;
            Expr** terms = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                terms[i] = together_recursive_ext(e->data.function.args[i], alpha);
            }

            Expr** nums = malloc(sizeof(Expr*) * count);
            Expr** dens = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                extract_num_den(terms[i], &nums[i], &dens[i]);
            }

            /* Compute LCM of denominators with Extension. */
            Expr* lcm_den = count > 0 ? expr_copy(dens[0]) : expr_new_integer(1);
            for (size_t i = 1; i < count; i++) {
                Expr* lcm_call = expr_new_function(
                    expr_new_symbol(SYM_PolynomialLCM),
                    (Expr*[]){
                        expr_copy(lcm_den),
                        expr_copy(dens[i]),
                        expr_new_function(
                            expr_new_symbol(SYM_Rule),
                            (Expr*[]){
                                expr_new_symbol(SYM_Extension),
                                expr_copy((Expr*)alpha)
                            }, 2)
                    }, 3);
                Expr* new_lcm = evaluate(lcm_call);
                expr_free(lcm_call);
                expr_free(lcm_den);
                lcm_den = new_lcm;
            }

            /* Pick a polynomial variable shared by all denominators (for
             * the PolynomialQuotient calls).  First non-α free symbol,
             * fall back to a sentinel `$Dummy` if everything is constant
             * (in which case the quotient is just lcm_den / dens[i] as
             * a number). */
            Expr* var = NULL;
            {
                size_t vc = 0, vcap = 8;
                Expr** vars = malloc(sizeof(Expr*) * vcap);
                collect_variables(lcm_den, &vars, &vc, &vcap);
                for (size_t i = 0; i < count; i++) {
                    collect_variables(dens[i], &vars, &vc, &vcap);
                }
                for (size_t k = 0; k < vc; k++) {
                    if (alpha && expr_eq(vars[k], (Expr*)alpha)) continue;
                    if (!var) var = expr_copy(vars[k]);
                }
                for (size_t k = 0; k < vc; k++) expr_free(vars[k]);
                free(vars);
                if (!var) var = expr_new_symbol("$Dummy");
            }

            /* Multiply each numerator by lcm_den/dens[i] and sum.
             * Pass Extension -> α to PolynomialQuotient so the division
             * runs in Q(α)[x] (qaupoly_divrem) rather than the slow
             * multivariate Q[α, x] path. */
            Expr** new_nums = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                Expr* q_call = expr_new_function(
                    expr_new_symbol(SYM_PolynomialQuotient),
                    (Expr*[]){
                        expr_copy(lcm_den),
                        expr_copy(dens[i]),
                        expr_copy(var),
                        expr_new_function(
                            expr_new_symbol(SYM_Rule),
                            (Expr*[]){
                                expr_new_symbol(SYM_Extension),
                                expr_copy((Expr*)alpha)
                            }, 2)
                    }, 4);
                Expr* q = evaluate(q_call);
                expr_free(q_call);
                new_nums[i] = eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Times),
                    (Expr*[]){nums[i], q}, 2));
            }
            expr_free(var);

            Expr* combined_num = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Plus), new_nums, count));
            Expr* inv_den = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Power),
                (Expr*[]){expr_copy(lcm_den), expr_new_integer(-1)}, 2));
            Expr* combined = eval_and_free(expr_new_function(
                expr_new_symbol(SYM_Times),
                (Expr*[]){combined_num, inv_den}, 2));

            Expr* result = cancel_with_extension(combined, alpha);
            if (!result) result = combined;
            else expr_free(combined);

            for (size_t i = 0; i < count; i++) {
                expr_free(terms[i]);
                expr_free(dens[i]);
            }
            free(terms);
            free(nums);
            free(dens);
            free(new_nums);
            expr_free(lcm_den);
            return result;
        }

        if (head == SYM_List || head == SYM_Equal || head == SYM_Less ||
            head == SYM_LessEqual || head == SYM_Greater ||
            head == SYM_GreaterEqual || head == SYM_And ||
            head == SYM_Or || head == SYM_Not) {
            size_t count = e->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = together_recursive_ext(e->data.function.args[i], alpha);
            }
            Expr* ret = eval_and_free(expr_new_function(
                expr_copy(e->data.function.head), args, count));
            free(args);
            return ret;
        }

        /* Other heads: recurse into args, then cancel-with-extension at
         * this level.  Mirrors together_recursive's catch-all branch. */
        size_t count = e->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            args[i] = together_recursive_ext(e->data.function.args[i], alpha);
        }
        Expr* ret = eval_and_free(expr_new_function(
            expr_copy(e->data.function.head), args, count));
        free(args);
        Expr* cancelled = cancel_with_extension(ret, alpha);
        if (!cancelled) return ret;
        expr_free(ret);
        return cancelled;
    }

    return cancel_with_extension(e, alpha);
}

/* Extension-aware Cancel: split into Numerator / Denominator (which
 * for a single-fraction input gives the polynomial parts directly,
 * and for a Plus expression term-threads structurally), compute
 * g = PolynomialGCD[num, den, Extension -> α] via the standard
 * builtin (which dispatches to qaupoly_gcd through poly.c), divide
 * both sides by g, and reassemble.  Returns NULL on any structural
 * failure (caller falls back to the no-extension path).
 *
 * Combining a sum-of-fractions argument first via Together is *not*
 * done here: that combine step would re-enter cancel_recursive over Q,
 * which on inputs with algebraic-number coefficients is the slow path
 * we are trying to avoid.  Callers (like Together-with-Extension) are
 * expected to pre-combine when they want sum-of-fractions handling. */
static Expr* cancel_with_extension(const Expr* arg, const Expr* alpha) {
    /* num = Numerator[arg], den = Denominator[arg]. */
    Expr* num_call = expr_new_function(
        expr_new_symbol(SYM_Numerator),
        (Expr*[]){expr_copy((Expr*)arg)}, 1);
    Expr* num = evaluate(num_call);
    expr_free(num_call);

    Expr* den_call = expr_new_function(
        expr_new_symbol(SYM_Denominator),
        (Expr*[]){expr_copy((Expr*)arg)}, 1);
    Expr* den = evaluate(den_call);
    expr_free(den_call);

    /* If denominator is 1 there's nothing to cancel. */
    if (den->type == EXPR_INTEGER && den->data.integer == 1) {
        expr_free(den);
        return num;
    }

    /* Strip algebraic-constant atoms (Sqrt[2], Sqrt[3], ...) that appear
     * in every summand of *both* num and den.  Mirrors the prepass in
     * cancel_recursive; without it the Extension path also misses these
     * cases, because PolynomialGCD[Sqrt[2], Sqrt[2] + Sqrt[2] x^4,
     * Extension -> α] still routes through the integer-content recursion
     * for its leading coefficient handling and returns 1.  The strip is
     * an alpha-independent algebraic identity (a/a = 1 in any ring), so
     * it is safe to run regardless of which extension the caller picked. */
    rat_strip_symbolic_common(&num, &den);
    if (den->type == EXPR_INTEGER && den->data.integer == 1) {
        expr_free(den);
        return num;
    }

    /* Expand num and den into genuine polynomials before the GCD.
     *
     * Numerator[arg]/Denominator[arg] frequently return UNEXPANDED products
     * -- e.g. when `arg` is a product of algebraic-coefficient fractions
     * (the cube-root Goursat descent's R(t(z))(1-z)^3 with the order-3
     * Mobius's Q(Sqrt[-3]) fixed points is exactly this shape).  Fed an
     * unexpanded product, PolynomialGCD[.,.,Extension -> alpha] only
     * partially reduces (its qaupoly path needs a coefficient list, so a
     * Times of factors -- especially one whose factors carry inner rational
     * or algebraic denominators -- yields a spuriously low-degree gcd or a
     * bare 1).  Expanding first turns both sides into proper Q(alpha)[x]
     * polynomials, after which the gcd and the two quotient divisions run in
     * the intended substrate.  (Cancel[Together[e, Extension -> Automatic],
     * Extension -> Automatic] previously left such an input unreduced.) */
    {
        Expr* num_call = expr_new_function(
            expr_new_symbol(SYM_Expand), (Expr*[]){ num }, 1);
        num = evaluate(num_call);   /* num_call adopts old num */
        expr_free(num_call);        /* frees the wrapper and the old num */

        Expr* den_call = expr_new_function(
            expr_new_symbol(SYM_Expand), (Expr*[]){ den }, 1);
        den = evaluate(den_call);
        expr_free(den_call);

        if (den->type == EXPR_INTEGER && den->data.integer == 1) {
            expr_free(den);
            return num;
        }
    }

    /* Compute g = PolynomialGCD[num, den, Extension -> alpha]. */
    Expr* gcd_call = expr_new_function(
        expr_new_symbol(SYM_PolynomialGCD),
        (Expr*[]){
            expr_copy(num),
            expr_copy(den),
            expr_new_function(
                expr_new_symbol(SYM_Rule),
                (Expr*[]){
                    expr_new_symbol(SYM_Extension),
                    expr_copy((Expr*)alpha)
                }, 2)
        }, 3);
    Expr* g = evaluate(gcd_call);
    expr_free(gcd_call);

    /* If GCD is trivial (1 or unevaluated), return the input as-is to
     * preserve current behavior. */
    bool trivial = false;
    if (g->type == EXPR_INTEGER && g->data.integer == 1) trivial = true;
    if (head_is(g, SYM_PolynomialGCD)) {
        trivial = true;
    }
    if (trivial) {
        expr_free(g);
        Expr* result = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times),
            (Expr*[]){
                num,
                eval_and_free(expr_new_function(
                    expr_new_symbol(SYM_Power),
                    (Expr*[]){den, expr_new_integer(-1)}, 2))
            }, 2));
        return result;
    }

    /* num/g and den/g via PolynomialQuotient.  We auto-detect the variable
     * inside the builtin; here we need to pass it explicitly.  Use the
     * first symbol that appears in the input (excluding alpha-render
     * symbols).  PolynomialQuotient with multivariate inputs threads the
     * division correctly via the existing implementation. */
    /* Pick a polynomial variable from `arg`. */
    Expr* var = NULL;
    {
        size_t vc = 0, vcap = 8;
        Expr** vars = malloc(sizeof(Expr*) * vcap);
        collect_variables((Expr*)arg, &vars, &vc, &vcap);
        for (size_t i = 0; i < vc; i++) {
            if (alpha && expr_eq(vars[i], (Expr*)alpha)) continue;
            if (!var) var = expr_copy(vars[i]);
        }
        for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
        free(vars);
    }
    if (!var) {
        expr_free(num); expr_free(den); expr_free(g);
        return NULL;
    }

    /* Pass Extension -> α so the division runs in Q(α)[x] rather than
     * dropping to the multivariate Q[α, x] path. */
    Expr* ext_rule_a = expr_new_function(
        expr_new_symbol(SYM_Rule),
        (Expr*[]){expr_new_symbol(SYM_Extension), expr_copy((Expr*)alpha)}, 2);
    Expr* new_num_call = expr_new_function(
        expr_new_symbol(SYM_PolynomialQuotient),
        (Expr*[]){num, expr_copy(g), expr_copy(var), ext_rule_a}, 4);
    Expr* new_num = evaluate(new_num_call);
    expr_free(new_num_call);

    Expr* ext_rule_b = expr_new_function(
        expr_new_symbol(SYM_Rule),
        (Expr*[]){expr_new_symbol(SYM_Extension), expr_copy((Expr*)alpha)}, 2);
    Expr* new_den_call = expr_new_function(
        expr_new_symbol(SYM_PolynomialQuotient),
        (Expr*[]){den, g, var, ext_rule_b}, 4);
    Expr* new_den = evaluate(new_den_call);
    expr_free(new_den_call);

    /* Bail when the Q(α)[x] quotient collapses to literal 0.  This is
     * the wrong-answer guard: it happens when the GCD over Q(α)[x] is
     * not a true divisor of `den` in the user's full algebraic-number
     * ring -- e.g. when `den` carries opaque non-α algebraic factors
     * (Sqrt[3], a nested Sqrt[radicand-with-α-inside], ...) that have
     * hidden relations to α the substrate doesn't see.  Building
     * Power[0, -1] here would silently inject ComplexInfinity into the
     * result (Simplify[D[Integrate[a x/(x^3+2), x], x]] reproduces this
     * if the bail is removed).  Returning NULL signals the caller to
     * fall back to the un-cancelled form -- correct-but-unsimplified,
     * which is the contract for cancel_with_extension. */
    if ((new_den->type == EXPR_INTEGER && new_den->data.integer == 0)
        || (new_den->type == EXPR_REAL && new_den->data.real == 0.0)) {
        expr_free(new_num);
        expr_free(new_den);
        return NULL;
    }

    /* Result = new_num / new_den. */
    Expr* inv_den = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Power),
        (Expr*[]){new_den, expr_new_integer(-1)}, 2));
    Expr* result = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Times),
        (Expr*[]){new_num, inv_den}, 2));
    return result;
}

/* Extension -> Automatic Cancel via PolynomialGCD + PolynomialQuotient.
 *
 * A robust fallback for the multi-generator tower case.  When the input's
 * splitting field needs several algebraic generators (e.g. Q(I, Sqrt[3]) =
 * Q(zeta_12), the field an order-3 Mobius map's fixed points live in),
 * qa_cancel_with_tower's gamma-substitution can decline (returns NULL) even
 * though PolynomialGCD[num, den, Extension -> Automatic] and
 * PolynomialQuotient[.., Extension -> Automatic] both reduce the input
 * correctly -- those builtins re-run extension_autodetect internally and pick
 * a working primitive element, whereas a single-generator cancel_with_extension
 * over just one of the generators cannot see the cross relations.
 *
 * This helper mirrors cancel_with_extension but threads Extension -> Automatic
 * (not a fixed alpha) through the inner polynomial calls, so it works for any
 * tower.  It expects `arg` to already be a single fraction P(x)/Q(x) (the
 * routing sends Plus inputs through together_recursive_ext first).  Returns the
 * reduced fraction, or NULL when the gcd is trivial / a quotient collapses /
 * no polynomial variable is present (caller falls back to the un-cancelled
 * form).  Borrows arg. */
static Expr* cancel_auto_gcd_quotient(const Expr* arg) {
    /* num = Expand[Numerator[arg]], den = Expand[Denominator[arg]]. */
    Expr* num_call = expr_new_function(expr_new_symbol(SYM_Numerator),
                                       (Expr*[]){expr_copy((Expr*)arg)}, 1);
    Expr* num0 = evaluate(num_call); expr_free(num_call);
    Expr* den_call = expr_new_function(expr_new_symbol(SYM_Denominator),
                                       (Expr*[]){expr_copy((Expr*)arg)}, 1);
    Expr* den0 = evaluate(den_call); expr_free(den_call);

    Expr* nx = expr_new_function(expr_new_symbol(SYM_Expand), (Expr*[]){num0}, 1);
    Expr* num = evaluate(nx); expr_free(nx);
    Expr* dx = expr_new_function(expr_new_symbol(SYM_Expand), (Expr*[]){den0}, 1);
    Expr* den = evaluate(dx); expr_free(dx);

    if (den->type == EXPR_INTEGER && den->data.integer == 1) {
        expr_free(den);
        return num;
    }

    /* Pick a polynomial variable (first non-alpha-render symbol in arg). */
    Expr* var = NULL;
    {
        size_t vc = 0, vcap = 8;
        Expr** vars = malloc(sizeof(Expr*) * vcap);
        collect_variables((Expr*)arg, &vars, &vc, &vcap);
        /* Skip the algebraic-atom render symbols (I, Sqrt[..], Power[c,p/q]);
         * collect_variables only returns bare EXPR_SYMBOLs, and the algebraic
         * atoms surface as Sqrt[]/Power[] functions or the symbol I, so the
         * first plain non-I symbol is the polynomial variable. */
        for (size_t i = 0; i < vc; i++) {
            if (vars[i]->type == EXPR_SYMBOL
                && vars[i]->data.symbol == SYM_I) continue;
            if (!var) var = expr_copy(vars[i]);
        }
        for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
        free(vars);
    }
    if (!var) { expr_free(num); expr_free(den); return NULL; }

    /* g = PolynomialGCD[num, den, Extension -> Automatic]. */
    Expr* gcd_call = expr_new_function(expr_new_symbol(SYM_PolynomialGCD),
        (Expr*[]){ expr_copy(num), expr_copy(den),
            expr_new_function(expr_new_symbol(SYM_Rule),
                (Expr*[]){ expr_new_symbol(SYM_Extension),
                           expr_new_symbol(SYM_Automatic) }, 2) }, 3);
    Expr* g = evaluate(gcd_call); expr_free(gcd_call);

    bool trivial = (g->type == EXPR_INTEGER && g->data.integer == 1)
                   || head_is(g, SYM_PolynomialGCD);
    if (trivial) { expr_free(g); expr_free(num); expr_free(den); expr_free(var); return NULL; }

    /* num/g and den/g via PolynomialQuotient[.., Extension -> Automatic]. */
    Expr* nq_call = expr_new_function(expr_new_symbol(SYM_PolynomialQuotient),
        (Expr*[]){ num, expr_copy(g), expr_copy(var),
            expr_new_function(expr_new_symbol(SYM_Rule),
                (Expr*[]){ expr_new_symbol(SYM_Extension),
                           expr_new_symbol(SYM_Automatic) }, 2) }, 4);
    Expr* nq = evaluate(nq_call); expr_free(nq_call);

    Expr* dq_call = expr_new_function(expr_new_symbol(SYM_PolynomialQuotient),
        (Expr*[]){ den, g, expr_copy(var),
            expr_new_function(expr_new_symbol(SYM_Rule),
                (Expr*[]){ expr_new_symbol(SYM_Extension),
                           expr_new_symbol(SYM_Automatic) }, 2) }, 4);
    Expr* dq = evaluate(dq_call); expr_free(dq_call);
    expr_free(var);

    if ((dq->type == EXPR_INTEGER && dq->data.integer == 0)
        || (dq->type == EXPR_REAL && dq->data.real == 0.0)) {
        expr_free(nq); expr_free(dq); return NULL;
    }

    Expr* inv_den = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ dq, expr_new_integer(-1) }, 2));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ nq, inv_den }, 2));
}

/* Count subnodes of `e` that mention any of `target`'s leaves. Used by
 * pick_best_tower_generator to rank tower generators by how much of
 * the input they "touch". The metric is "subnodes where target appears
 * as a structural sub-expression" — a conservative proxy for how much
 * a single-α Cancel/Together over Q(target) can simplify. */
static int subnode_mentions(const Expr* e, const Expr* target) {
    if (!e) return 0;
    if (expr_eq((Expr*)e, (Expr*)target)) return 1;
    if (e->type != EXPR_FUNCTION) return 0;
    int sum = 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        sum += subnode_mentions(e->data.function.args[i], target);
    }
    return sum;
}

/* Layer-3 guided fallback: instead of trying all N tower generators
 * sequentially when qa_cancel_with_tower declines, pick the one that
 * occurs most often in the input. Returns the index of the chosen
 * generator, or -1 if none have any occurrences (in which case the
 * fallback would have produced nothing useful anyway).
 *
 * Rationale: in a tower with n >= 2 generators, the single-α path can
 * only collapse fractions whose denominators are polynomials in THAT
 * α; generators that don't appear in any denominator give no
 * cancellation. Ranking by occurrence count picks the generator most
 * likely to be load-bearing. */
static int pick_best_tower_generator(const Expr* arg, const QATower* t) {
    int best_idx = -1;
    int best_count = 0;
    for (int i = 0; i < t->n; i++) {
        const Expr* gen = t->alpha_renders[i];
        if (!gen) continue;
        int c = subnode_mentions(arg, gen);
        if (c > best_count) {
            best_count = c;
            best_idx = i;
        }
    }
    return best_idx;
}

#ifdef USE_FLINT
/* True if `e` contains an algebraic-number atom: I, Sqrt[..], Complex[..], or
 * Power[base, p/q] with a non-integer rational exponent.  Distinguishes an
 * over-an-extension fraction from a plain rational one. */
static bool rat_has_algebraic_atom(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol == SYM_I;
    if (e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h && h->type == EXPR_SYMBOL) {
        if (h->data.symbol == SYM_Sqrt || h->data.symbol == SYM_Complex)
            return true;
        if (h->data.symbol == SYM_Power && e->data.function.arg_count == 2) {
            int64_t p, q;
            if (is_rational(e->data.function.args[1], &p, &q) && q != 1)
                return true;
        }
    }
    if (rat_has_algebraic_atom(h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rat_has_algebraic_atom(e->data.function.args[i])) return true;
    return false;
}

/* Cancel a single fraction rigorously over a detected algebraic-extension field
 * (Q(sqrt d), Q(zeta_n), radical tower) using the FLINT engine: g = gcd(num,den)
 * then num/g, den/g by exact division. Returns NULL — deferring to the classical
 * extension path — when no algebraic generator is present, the fraction is
 * already reduced (g = 1), or an exact division does not apply. */
static Expr* flint_cancel_fraction(Expr* arg) {
    /* Parametric radical field Q(t..)(sqrt k): normalise the whole rational
     * function via fmpz_mpoly_q. This covers a Plus of fractions (a sum whose
     * Numerator/Denominator does not combine, so the extract-num/den path below
     * bails and drops to the slow QA extension path — the LogToReal A^2+B^2
     * bottleneck). */
    {
        Expr* fn = flint_parametric_field_normalize(arg);
        if (fn) return fn;
    }
    Expr* num; Expr* den;
    extract_num_den(arg, &num, &den);
    if (den->type == EXPR_INTEGER && den->data.integer == 1) {
        expr_free(num); expr_free(den); return NULL;
    }
    /* Parametric radical Q(t..)(sqrt k): the classical fallback does NOT
     * terminate for a symbolic radicand, so once FLINT engages this field we
     * finish the reduction ourselves even when num, den are already coprime
     * (g == 1). Try it explicitly first so we can distinguish this regime. */
    Expr* pg = flint_parametric_sqrt_gcd(num, den);
    int parametric = (pg != NULL);
    Expr* g;
    if (parametric) {
        g = eval_and_free(pg);
    } else {
        /* Number field / cyclotomic / tower: the classical extension path is
         * fast AND does more (e.g. the roots-of-unity simplification the
         * cube-root Goursat descent relies on). Only take over when there is an
         * actual common factor; hand a coprime pair (g == 1) back to it. */
        g = flint_extension_gcd(num, den);
        if (!g) { expr_free(num); expr_free(den); return NULL; }
        g = eval_and_free(g);   /* the bridge returns a raw tree; normalise it */
    }
    if (g->type == EXPR_INTEGER && g->data.integer == 1) {
        expr_free(g);
        if (!parametric) { expr_free(num); expr_free(den); return NULL; }
        /* Parametric coprime: reduced form is num/den itself; returning NULL
         * would drop into the non-terminating classical path. */
        Expr* inv1 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){den, expr_new_integer(-1)}, 2));
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){num, inv1}, 2));
    }
    Expr* nn = flint_extension_divexact(num, g);
    Expr* nd = flint_extension_divexact(den, g);
    expr_free(g);
    expr_free(num);
    expr_free(den);
    if (!nn || !nd) {
        if (nn) expr_free(nn);
        if (nd) expr_free(nd);
        return NULL;
    }
    if (nd->type == EXPR_INTEGER && nd->data.integer == 1) {
        expr_free(nd);
        return eval_and_free(nn);
    }
    /* Capture whether the divided-out denominator is a numeric constant before
     * it is consumed below. */
    bool nd_numeric = (nd->type == EXPR_INTEGER) || head_is(nd, SYM_Rational);
    Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){nd, expr_new_integer(-1)}, 2));
    Expr* prod = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){nn, inv}, 2));
    /* A numeric denominator (e.g. the unit -1 from a GCD that FLINT normalised
     * monic in its own variable order rather than in x) leaves the result as
     * Times[scalar, Plus[...]], which does not auto-distribute — expand so the
     * sign/scale folds in and the surface form matches the classical path
     * (Sqrt[k] + x, not -(-Sqrt[k] - x)). */
    if (nd_numeric) {
        Expr* out = expr_expand(prod);
        expr_free(prod);
        return out;
    }
    return prod;
}
#endif

static Expr* builtin_cancel_compute(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;

#ifdef USE_FLINT
    /* Rigorous algebraic-extension cancellation via FLINT, ahead of the QA
     * extension path. Returns NULL for plain-rational / non-extension inputs,
     * which then take the existing code below. */
    if (res->data.function.arg_count == 1) {
        Expr* fc = flint_cancel_fraction(res->data.function.args[0]);
        if (fc) return fc;
    }
#endif

    /* Strip a trailing Extension -> α option, if any.  When the user
     * passed `Extension -> Automatic`, run extension_autodetect on the
     * polynomial argument and route through the existing α-path when a
     * single-generator extension is found (tier-1 limitation; towers
     * with n ≥ 2 fall back to the no-extension path). */
    size_t argc = res->data.function.arg_count;
    bool auto_flag = false;
    const Expr* alpha = extract_extension_option_full(res, &argc, &auto_flag);
    if (argc != 1) {
        if (argc == 0) return NULL;
        /* Wrong arg count after stripping options. */
        if (argc != 1) return NULL;
    }
    Expr* arg = res->data.function.args[0];

#ifdef USE_FLINT
    /* Proven reduction first for a single algebraic-coefficient fraction:
     * PolynomialGCD + PolynomialQuotient with Extension -> Automatic reduce it
     * correctly because they each re-detect the whole splitting field, whereas
     * both FLINT's number-field GCD (flint_cancel_fraction below) and the QA
     * tower path can leave such a fraction un-cancelled -- notably when the
     * numerator carries algebraic coefficient denominators like 1/Sqrt[3], the
     * shape the cube-root Goursat descent's Cancel[Together[R(t(z))(1-z)^3,
     * Extension], Extension] produces over Q(I, Sqrt[3]).  Gated on an actual
     * algebraic atom (plain-rational Cancel is untouched) and taken only on a
     * strict leaf-count shrink, so already-reduced inputs fall through to the
     * FLINT normalisation below. */
    if (auto_flag && rat_has_algebraic_atom(arg)) {
        bool is_sum = (arg->type == EXPR_FUNCTION && arg->data.function.head
                       && arg->data.function.head->type == EXPR_SYMBOL
                       && arg->data.function.head->data.symbol == SYM_Plus);
        if (!is_sum) {
            Expr* red = cancel_auto_gcd_quotient(arg);
            if (red) {
                int64_t in_sz  = leaf_count_internal((Expr*)arg, true);
                int64_t out_sz = leaf_count_internal(red, true);
                if (out_sz < in_sz) return red;
                expr_free(red);
            }
        }
    }

    /* Also route the Cancel[e, Extension -> ...] form (2-arg, after option
     * stripping) through FLINT — the QA autodetect path below does not detect a
     * symbolic radicand Sqrt[k] (parametric Q(a,b,k)(Sqrt k)) and the multivariate
     * Together/Cancel it falls back to hangs. The 1-arg fast path handled the
     * bare Cancel[e] above; this covers the Extension-carrying calls (e.g. the
     * Goursat descent's canonic = Cancel[Together[e, Ext], Ext]). */
    {
        Expr* fc = flint_cancel_fraction(arg);
        if (fc) return fc;
    }
#endif

    Expr* alpha_auto = NULL;
    QATower* auto_tower = NULL;
    /* Layer-0 input-level prefilter: extension_autodetect itself is
     * expensive on inputs containing nested radicals (the primitive-
     * element compositum construction runs over the whole tower). When
     * a nested radical sits inside a Power[Plus[...], p/q] (q ≥ 2) we
     * can predict from the input alone that the tower path will be
     * rejected — the γ-substitution leaves the outer Power opaque to
     * the no-extension Together (see qa_tower_has_nested_radical doc).
     * Skip the autodetect call entirely and fall straight through to
     * the no-extension path. */
    if (!alpha && auto_flag && !expr_has_nested_radical_radicand(arg)) {
        auto_tower = extension_autodetect(arg);
        if (auto_tower && auto_tower->n == 1) {
            alpha_auto = expr_copy(auto_tower->alpha_renders[0]);
            alpha = alpha_auto;
        } else if (auto_tower && auto_tower->n >= 2) {
            /* Multi-generator tower: route through qa_cancel_with_tower
             * which substitutes each α_i with its Q(γ) form, lifts to
             * QAUPoly over Q(γ), and runs qaupoly_gcd-based cancellation
             * end-to-end.  Returns NULL on lift failure; the caller
             * falls back to the no-extension path. */
            if (qa_tower_has_nested_radical(auto_tower)) {
                /* Layer-0 prefilter: when the tower has a nested-radical
                 * α_i (surfaced from a Power[Plus[...], p/q] node in the
                 * input), Step 1's γ-substitution leaves the surrounding
                 * Power[..., p/q] opaque to Step 4's no-extension
                 * Together — the result inflates past the leaf-count
                 * gate and gets rejected. Same gate also skips the
                 * N×fallback loop because together_recursive_ext on a
                 * nested-radical α_i hits the same opaque-Power wall.
                 * Drop straight to the no-extension path; the
                 * surrounding Simplify search picks up the denesting
                 * via simp_denest_sqrt on individual subnodes. */
                qa_tower_free(auto_tower);
                auto_tower = NULL;
            } else {
                Expr* tower_result = qa_cancel_with_tower(arg, auto_tower);
                if (tower_result) {
                    qa_tower_free(auto_tower);
                    return tower_result;
                }
                /* Layer-3 guided fallback. The previous linear N-try
                 * loop tried every tower generator (up to 4) and kept
                 * the smallest result, paying for each
                 * together_recursive_ext call sequentially. The
                 * single-α path's effectiveness is proportional to how
                 * often the chosen α appears in the input — generators
                 * that don't appear in any denominator can't collapse
                 * anything. Picking the most-used α gives the same
                 * result as the linear loop in most cases at
                 * ~1/N the cost. If the top pick doesn't beat the
                 * input, fall through (the linear loop's other tries
                 * would have failed the leaf-count gate anyway). */
                {
                    bool is_sum = (arg->type == EXPR_FUNCTION
                                   && arg->data.function.head
                                   && arg->data.function.head->type == EXPR_SYMBOL
                                   && arg->data.function.head->data.symbol == SYM_Plus);
                    int idx = pick_best_tower_generator(arg, auto_tower);
                    Expr* best = NULL;
                    if (idx >= 0) {
                        const Expr* gen = auto_tower->alpha_renders[idx];
                        Expr* cand = is_sum
                            ? together_recursive_ext(arg, gen)
                            : cancel_with_extension(arg, gen);
                        if (cand) {
                            Expr* tg = expr_new_function(expr_new_symbol(SYM_Together),
                                (Expr*[]){cand}, 1);
                            Expr* folded = evaluate(tg);
                            expr_free(tg);
                            if (folded) {
                                int64_t in_size = leaf_count_internal((Expr*)arg, true);
                                int64_t lc = leaf_count_internal(folded, true);
                                if (lc < in_size) {
                                    best = folded;
                                } else {
                                    expr_free(folded);
                                }
                            }
                        }
                    }
                    qa_tower_free(auto_tower);
                    auto_tower = NULL;
                    if (best) return best;
                }
            }
        }
    }

    if (alpha) {
        /* For Plus / sum-of-fractions inputs the user typically expects
         * "combine into a single fraction over Q(α), then cancel".
         * together_recursive_ext threads Extension -> α through the
         * underlying PolynomialLCM and PolynomialQuotient calls, so the
         * inner polynomial arithmetic stays in Q(α)[x] (no multivariate
         * Q[α, x] GCD blowup). */
        bool is_sum = (arg->type == EXPR_FUNCTION
                       && arg->data.function.head
                       && arg->data.function.head->type == EXPR_SYMBOL
                       && arg->data.function.head->data.symbol == SYM_Plus);
        if (is_sum) {
            Expr* result = together_recursive_ext(arg, alpha);
            if (result) {
                if (alpha_auto) expr_free(alpha_auto);
                if (auto_tower) qa_tower_free(auto_tower);
                return result;
            }
        } else {
            Expr* result = cancel_with_extension(arg, alpha);
            if (result) {
                if (alpha_auto) expr_free(alpha_auto);
                if (auto_tower) qa_tower_free(auto_tower);
                return result;
            }
        }
        /* Fall back to non-extension cancel on failure. */
        if (alpha_auto) { expr_free(alpha_auto); alpha_auto = NULL; }
        if (auto_tower) { qa_tower_free(auto_tower); auto_tower = NULL; }
    }

    /* General Extension -> Automatic fallback for a single fraction P(x)/Q(x).
     *
     * The per-alpha paths above can miss the full splitting field: the
     * single-generator autodetect drops a generator (notably I), and the
     * multi-generator qa_cancel_with_tower can decline outright.  But
     * PolynomialGCD and PolynomialQuotient with Extension -> Automatic each
     * re-run extension_autodetect internally and reduce over a working
     * primitive element -- so cancelling a single fraction directly through
     * them succeeds where the routing above returned the input unreduced
     * (e.g. Cancel[P/Q, Extension -> Automatic] with P carrying both I and
     * Sqrt[3] coefficients, the shape the cube-root Goursat descent produces).
     * Accept only when it strictly shrinks the input, matching the guided
     * fallback's contract; on no-change, fall through to the legacy paths. */
    if (auto_flag) {
        bool is_sum = (arg->type == EXPR_FUNCTION && arg->data.function.head
                       && arg->data.function.head->type == EXPR_SYMBOL
                       && arg->data.function.head->data.symbol == SYM_Plus);
        if (!is_sum) {
            Expr* red = cancel_auto_gcd_quotient(arg);
            if (red) {
                int64_t in_size = leaf_count_internal((Expr*)arg, true);
                int64_t lc = leaf_count_internal(red, true);
                if (lc < in_size) return red;
                expr_free(red);
            }
        }
    }

    /* Phase E: single-generator polynomial-radicand path.  Triggered
     * when Extension -> Automatic is requested and the standard
     * autodetect path didn't produce a tower (typically: the input
     * contains a Sqrt[poly] / Power[poly, 1/q] whose radicand has free
     * symbols, which qa_resolve_nested_radical's expr_collect_
     * atomic_algebraics rejects).  See qa_cancel_with_poly_radical
     * for the full algorithm + the documented Cardano-style multi-
     * radical limitation. */
    if (auto_flag && !alpha) {
        Expr* result = qa_cancel_with_poly_radical(arg);
        if (result) {
            if (alpha_auto) expr_free(alpha_auto);
            if (auto_tower) qa_tower_free(auto_tower);
            return result;
        }
    }

    /* Inexact coefficients block the rational-arithmetic cancellation
     * machinery (e.g. (x^2/9 - y^2/25.) — the 25. defeats Polynomial GCD).
     * Force-rationalise on the way in, numericalise on the way out. */
    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_cancel);
    }

    /* Algebraic-generator pass: if the input contains a sub-expression
     * u with fractional rational exponents (u^(p/q), q != 1), substitute
     * u -> g^m, run Together's recursive normalisation in g, then
     * substitute back.  together_recursive is preferred over
     * cancel_recursive directly because the substituted form may have
     * Plus terms with different g-exponents that need fraction
     * combination before GCD cancellation can fire. */
    {
        Expr* base = NULL;
        Expr* atom = NULL;
        int64_t m = 1;
        if (poly_find_radical_gen(arg, &base, &atom, &m)) {
            char* gen = poly_make_fresh_gen(arg);
            Expr* substituted = poly_subst_radical_to_gen(arg, base, atom, m, gen);
            Expr* combined = together_recursive(substituted);
            expr_free(substituted);
            Expr* final = poly_subst_radical_from_gen(combined, base, atom, m, gen);
            expr_free(combined);
            expr_free(base);
            expr_free(atom);
            free(gen);
            return final;
        }
    }

    return cancel_recursive(arg);
}

static Expr* together_recursive(Expr* e) {
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* head = e->data.function.head->data.symbol;
        
        if (head == SYM_Plus) {
            size_t count = e->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = together_recursive(e->data.function.args[i]);
            }
            
            Expr** nums = malloc(sizeof(Expr*) * count);
            Expr** dens = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                extract_num_den(args[i], &nums[i], &dens[i]);
            }
            
            Expr* lcm_den = count > 0 ? expr_copy(dens[0]) : expr_new_integer(1);
            for (size_t i = 1; i < count; i++) {
                Expr* call_lcm = expr_new_function(expr_new_symbol(SYM_PolynomialLCM), (Expr*[]){expr_copy(lcm_den), expr_copy(dens[i])}, 2);
                Expr* new_lcm = evaluate(call_lcm);
                expr_free(call_lcm);
                expr_free(lcm_den);
                lcm_den = new_lcm;
            }
            
            /* Strict-quotient combine: each lcm_den/dens[i] must divide
             * exactly in Q[vars]. If any quotient is not exact, the
             * iterative PolynomialLCM was unsound for the working
             * multivariate ring (typically: dens[i] divides lcm_den only
             * in some algebraic-extension ring, e.g. multi-radical sums
             * over Q like Sqrt[2]a+Sqrt[3]+Sqrt[5]+Sqrt[7]). Propagating
             * a symbolic Times[lcm_den, Power[dens[i], -1]] leaves a
             * Power[Plus[...], -1] inside the combined numerator, which
             * sends the downstream cancel_recursive -> PolynomialGCD
             * path into multivariate Euclid on a rational input -- this
             * is the case-13 hang. Bail cleanly: the Plus stays
             * uncombined and Simplify's simp_algebraic seed handles the
             * extension structure when applicable. */
            Expr** new_nums = malloc(sizeof(Expr*) * count);
            bool exact = true;
            for (size_t i = 0; i < count; i++) {
                Expr* q = cancel_exact_div_strict(lcm_den, dens[i]);
                if (!q) {
                    exact = false;
                    /* Free the partial quotients we already wrapped. */
                    for (size_t j = 0; j < i; j++) expr_free(new_nums[j]);
                    new_nums[i] = NULL;
                    break;
                }
                new_nums[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_copy(nums[i]), q}, 2));
            }

            if (!exact) {
                /* Cleanup and return the input Plus unchanged (the
                 * already-together'd summands are owned in args[]). */
                Expr* ret;
                if (count == 1) {
                    ret = expr_copy(args[0]);
                } else {
                    Expr** plus_args = malloc(sizeof(Expr*) * count);
                    for (size_t i = 0; i < count; i++) plus_args[i] = expr_copy(args[i]);
                    ret = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), plus_args, count));
                    free(plus_args);
                }
                for (size_t i = 0; i < count; i++) {
                    expr_free(args[i]);
                    expr_free(nums[i]);
                    expr_free(dens[i]);
                }
                free(args); free(nums); free(dens); free(new_nums);
                expr_free(lcm_den);
                return ret;
            }

            /* All quotients exact -- proceed with the standard combine. */
            Expr* combined_num = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), new_nums, count));
            Expr* inv_den = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_copy(lcm_den), expr_new_integer(-1)}, 2));
            Expr* combined_expr = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){combined_num, inv_den}, 2));

            Expr* res = cancel_recursive(combined_expr);
            expr_free(combined_expr);

            for (size_t i = 0; i < count; i++) {
                expr_free(args[i]);
                expr_free(nums[i]);
                expr_free(dens[i]);
            }
            free(args);
            free(nums);
            free(dens);
            free(new_nums);

            expr_free(lcm_den);

            return res;
        }
        
        if (head == SYM_List ||
            head == SYM_Equal || head == SYM_Less ||
            head == SYM_LessEqual || head == SYM_Greater ||
            head == SYM_GreaterEqual || head == SYM_And ||
            head == SYM_Or || head == SYM_Not) {
            
            size_t count = e->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = together_recursive(e->data.function.args[i]);
            }
            Expr* ret = eval_and_free(expr_new_function(expr_copy(e->data.function.head), args, count));
            free(args);
            return ret;
        }
        
        size_t count = e->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            args[i] = together_recursive(e->data.function.args[i]);
        }
        Expr* ret = eval_and_free(expr_new_function(expr_copy(e->data.function.head), args, count));
        free(args);
        
        Expr* cancelled = cancel_recursive(ret);
        expr_free(ret);
        return cancelled;
    }
    
    return cancel_recursive(e);
}

static Expr* builtin_together_compute(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;

    /* Strip a trailing Extension -> α option, if any.  `Extension ->
     * Automatic` triggers extension_autodetect on the polynomial
     * argument; a single-generator detection routes through the
     * existing α path. */
    size_t argc = res->data.function.arg_count;
    bool auto_flag = false;
    const Expr* alpha = extract_extension_option_full(res, &argc, &auto_flag);
    if (argc != 1) return NULL;
    Expr* arg = res->data.function.args[0];

#ifdef USE_FLINT
    /* Together = combine into a single fraction + cancel common factors, which
     * is exactly flint_cancel_fraction (extract num/den, GCD, exact-divide).
     * Route the algebraic-extension cases (incl. the parametric Q(a,b,k)(Sqrt k)
     * the QA autodetect below cannot see) through FLINT; NULL for plain-rational
     * or nothing-to-do inputs falls through to the classical path. Fixes the
     * Goursat descent's together_ext over the radical tower/parametric ring. */
    {
        Expr* fc = flint_cancel_fraction(arg);
        if (fc) return fc;
    }
#endif

    Expr* alpha_auto = NULL;
    QATower* auto_tower = NULL;
    /* Layer-0 input-level prefilter: see builtin_cancel for the
     * rationale. extension_autodetect's primitive-element construction
     * is the bulk of the cost for nested-radical inputs, and the tower
     * path is going to be rejected by qa_tower_has_nested_radical
     * downstream anyway. Skip both. */
    if (!alpha && auto_flag && !expr_has_nested_radical_radicand(arg)) {
        auto_tower = extension_autodetect(arg);
        if (auto_tower && auto_tower->n == 1) {
            alpha_auto = expr_copy(auto_tower->alpha_renders[0]);
            alpha = alpha_auto;
        } else if (auto_tower && auto_tower->n >= 2
                   && !internal_args_contain_inexact(res)) {
            /* Multi-generator tower: route through qa_cancel_with_tower.
             * For Together this gives "combine + cancel over Q(γ)" in
             * one shot — the function combines fractions internally
             * before lifting. */
            if (qa_tower_has_nested_radical(auto_tower)) {
                /* Layer-0 prefilter: nested-radical α_i (originating from
                 * Power[Plus[...], p/q] in the input) cannot be reduced
                 * by the γ-substitution path — Step 1 leaves the
                 * Power[non_int_base, p/q] structurally opaque to
                 * Step 4's no-extension Together. Tower call would
                 * inflate, hit the leaf-count gate, and return NULL.
                 * The N×fallback hits the same wall. Drop to the
                 * no-extension path; surrounding Simplify handles
                 * denesting via simp_denest_sqrt. */
                qa_tower_free(auto_tower);
                auto_tower = NULL;
            } else {
                Expr* tower_result = qa_cancel_with_tower(arg, auto_tower);
                if (tower_result) {
                    qa_tower_free(auto_tower);
                    return tower_result;
                }
                /* Layer-3 guided fallback (see builtin_cancel for full
                 * rationale). Pick the generator that occurs most often
                 * in the input — that's the α whose single-α Together
                 * is most likely to collapse fractions. The remaining
                 * generators pass through as opaque coefficients in
                 * together_recursive_ext, and the final no-extension
                 * Together fold-up combines like-coefficient terms. */
                {
                    int idx = pick_best_tower_generator(arg, auto_tower);
                    Expr* best = NULL;
                    if (idx >= 0) {
                        const Expr* gen = auto_tower->alpha_renders[idx];
                        Expr* cand = together_recursive_ext(arg, gen);
                        if (cand) {
                            Expr* tg = expr_new_function(expr_new_symbol(SYM_Together),
                                (Expr*[]){cand}, 1);
                            Expr* folded = evaluate(tg);
                            expr_free(tg);
                            if (folded) {
                                int64_t in_size = leaf_count_internal((Expr*)arg, true);
                                int64_t lc = leaf_count_internal(folded, true);
                                if (lc < in_size) {
                                    best = folded;
                                } else {
                                    expr_free(folded);
                                }
                            }
                        }
                    }
                    qa_tower_free(auto_tower);
                    auto_tower = NULL;
                    if (best) return best;
                }
            }
        }
    }

    /* Extension-aware path: combine + cancel in Q(α)[x] directly, with
     * Extension -> α threaded through PolynomialLCM and PolynomialQuotient
     * so the inner arithmetic stays on the qaupoly substrate. */
    if (alpha && !internal_args_contain_inexact(res)) {
        Expr* result = together_recursive_ext(arg, alpha);
        if (result) {
            if (alpha_auto) expr_free(alpha_auto);
            if (auto_tower) qa_tower_free(auto_tower);
            return result;
        }
        if (alpha_auto) { expr_free(alpha_auto); alpha_auto = NULL; }
        if (auto_tower) { qa_tower_free(auto_tower); auto_tower = NULL; }
    }

    /* Phase E: single-generator polynomial-radicand path.  Same wiring
     * as builtin_cancel_compute; see qa_cancel_with_poly_radical for
     * the algorithm. */
    if (auto_flag && !alpha) {
        Expr* result = qa_cancel_with_poly_radical(arg);
        if (result) {
            if (alpha_auto) expr_free(alpha_auto);
            if (auto_tower) qa_tower_free(auto_tower);
            return result;
        }
    }

    /* Inexact coefficients can't be combined by exact polynomial-LCM
     * machinery; rationalise inputs, run, then re-numericalise. */
    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_together);
    }

    /* Algebraic-generator pass (see builtin_cancel for the rationale). */
    {
        Expr* base = NULL;
        Expr* atom = NULL;
        int64_t m = 1;
        if (poly_find_radical_gen(arg, &base, &atom, &m)) {
            char* gen = poly_make_fresh_gen(arg);
            Expr* substituted = poly_subst_radical_to_gen(arg, base, atom, m, gen);
            Expr* combined = together_recursive(substituted);
            expr_free(substituted);
            Expr* final = poly_subst_radical_from_gen(combined, base, atom, m, gen);
            expr_free(combined);
            expr_free(base);
            expr_free(atom);
            free(gen);
            return final;
        }
    }

    return together_recursive(arg);
}

/* Layer-5 memoisation wrapper for Cancel.  builtin_simplify pushes a
 * FactorMemo for the duration of one Simplify call (see
 * simp.c:builtin_simplify), and the per-subnode descent in simp_search
 * re-asks the same `Cancel[arg, options...]` repeatedly across rounds.
 * Keying on the full `res` expression (which already encodes the args
 * and any Extension -> α option) hits cleanly when the same call is
 * issued from different code paths in the search.  Skip the cache
 * when no Simplify is active (factor_memo_active() returns NULL) so
 * standalone Cancel calls see no overhead. */
Expr* builtin_cancel(Expr* res) {
    FactorMemo* memo = factor_memo_active();
    if (memo) {
        const Expr* hit = factor_memo_lookup(memo, res);
        if (hit) return expr_copy((Expr*)hit);
    }
    Expr* out = builtin_cancel_compute(res);
    if (out && memo) factor_memo_store(memo, res, out);
    return out;
}

/* Layer-5 memoisation wrapper for Together (see builtin_cancel above). */
Expr* builtin_together(Expr* res) {
    FactorMemo* memo = factor_memo_active();
    if (memo) {
        const Expr* hit = factor_memo_lookup(memo, res);
        if (hit) return expr_copy((Expr*)hit);
    }
    Expr* out = builtin_together_compute(res);
    if (out && memo) factor_memo_store(memo, res, out);
    return out;
}

void rat_init(void) {
    symtab_add_builtin("Numerator", builtin_numerator);
    symtab_get_def("Numerator")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_add_builtin("Denominator", builtin_denominator);
    symtab_get_def("Denominator")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_add_builtin("Cancel", builtin_cancel);
    symtab_get_def("Cancel")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_add_builtin("Together", builtin_together);
    symtab_get_def("Together")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
}
