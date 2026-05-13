#include "rat.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "poly.h"
#include "expand.h"
#include "rationalize.h"
#include "sym_names.h"
#include "options.h"
#include "qafactor.h"
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
    return eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), expr_copy(e)}, 2));
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
        *num_out = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(expr), expr_copy(d_expr)}, 2));
        return;
    }

    if (expr->type == EXPR_FUNCTION && expr->data.function.head->type == EXPR_SYMBOL && (expr->data.function.head->data.symbol == SYM_Power || expr->data.function.head->data.symbol == SYM_Exp)) {
        bool is_exp = expr->data.function.head->data.symbol == SYM_Exp;
        Expr* base = is_exp ? expr_new_symbol("E") : expr->data.function.args[0];
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
            else if (n_c == 1) *num_out = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), num_args[0]}, 2));
            else {
                Expr* p = eval_and_free(expr_new_function(expr_new_symbol("Plus"), num_args, n_c));
                *num_out = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), p}, 2));
            }
            if (d_c == 0) *den_out = expr_new_integer(1);
            else if (d_c == 1) *den_out = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), den_args[0]}, 2));
            else {
                Expr* p = eval_and_free(expr_new_function(expr_new_symbol("Plus"), den_args, d_c));
                *den_out = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), p}, 2));
            }
            free(num_args);
            free(den_args);
            if (is_exp) expr_free(base);
            return;
        } else {
            if (is_superficially_negative(exp)) {
                *num_out = expr_new_integer(1);
                *den_out = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), negate_expr(exp)}, 2));
            } else {
                if (is_exp) {
                    *num_out = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), expr_copy(exp)}, 2));
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
        else *num_out = eval_and_free(expr_new_function(expr_new_symbol("Times"), n_args, n_c));

        if (d_c == 0) *den_out = expr_new_integer(1);
        else if (d_c == 1) *den_out = d_args[0];
        else *den_out = eval_and_free(expr_new_function(expr_new_symbol("Times"), d_args, d_c));

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

    if (res) {
        expr_free(exp_num);
        expr_free(exp_den);
        return res;
    } else {
        Expr* t = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){exp_den, expr_new_integer(-1)}, 2));
        Expr* r = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){exp_num, t}, 2));
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

    Expr* g = eval_and_free(expr_new_function(expr_new_symbol("PolynomialGCD"), (Expr*[]){expr_copy(num), expr_copy(den)}, 2));
    
    Expr* new_num = cancel_exact_div_wrapper(num, g);
    Expr* new_den = cancel_exact_div_wrapper(den, g);
    
    if (new_num && new_den && den_has_negative_lead(new_den)) {
        Expr* t1 = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), new_num}, 2));
        Expr* t2 = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), new_den}, 2));
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
        Expr* inv_den = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){new_den, expr_new_integer(-1)}, 2));
        res = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){new_num, inv_den}, 2));
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
                    expr_new_symbol("PolynomialLCM"),
                    (Expr*[]){
                        expr_copy(lcm_den),
                        expr_copy(dens[i]),
                        expr_new_function(
                            expr_new_symbol("Rule"),
                            (Expr*[]){
                                expr_new_symbol("Extension"),
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
                    expr_new_symbol("PolynomialQuotient"),
                    (Expr*[]){
                        expr_copy(lcm_den),
                        expr_copy(dens[i]),
                        expr_copy(var),
                        expr_new_function(
                            expr_new_symbol("Rule"),
                            (Expr*[]){
                                expr_new_symbol("Extension"),
                                expr_copy((Expr*)alpha)
                            }, 2)
                    }, 4);
                Expr* q = evaluate(q_call);
                expr_free(q_call);
                new_nums[i] = eval_and_free(expr_new_function(
                    expr_new_symbol("Times"),
                    (Expr*[]){nums[i], q}, 2));
            }
            expr_free(var);

            Expr* combined_num = eval_and_free(expr_new_function(
                expr_new_symbol("Plus"), new_nums, count));
            Expr* inv_den = eval_and_free(expr_new_function(
                expr_new_symbol("Power"),
                (Expr*[]){expr_copy(lcm_den), expr_new_integer(-1)}, 2));
            Expr* combined = eval_and_free(expr_new_function(
                expr_new_symbol("Times"),
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
        expr_new_symbol("Numerator"),
        (Expr*[]){expr_copy((Expr*)arg)}, 1);
    Expr* num = evaluate(num_call);
    expr_free(num_call);

    Expr* den_call = expr_new_function(
        expr_new_symbol("Denominator"),
        (Expr*[]){expr_copy((Expr*)arg)}, 1);
    Expr* den = evaluate(den_call);
    expr_free(den_call);

    /* If denominator is 1 there's nothing to cancel. */
    if (den->type == EXPR_INTEGER && den->data.integer == 1) {
        expr_free(den);
        return num;
    }

    /* Compute g = PolynomialGCD[num, den, Extension -> alpha]. */
    Expr* gcd_call = expr_new_function(
        expr_new_symbol("PolynomialGCD"),
        (Expr*[]){
            expr_copy(num),
            expr_copy(den),
            expr_new_function(
                expr_new_symbol("Rule"),
                (Expr*[]){
                    expr_new_symbol("Extension"),
                    expr_copy((Expr*)alpha)
                }, 2)
        }, 3);
    Expr* g = evaluate(gcd_call);
    expr_free(gcd_call);

    /* If GCD is trivial (1 or unevaluated), return the input as-is to
     * preserve current behavior. */
    bool trivial = false;
    if (g->type == EXPR_INTEGER && g->data.integer == 1) trivial = true;
    if (g->type == EXPR_FUNCTION
        && g->data.function.head
        && g->data.function.head->type == EXPR_SYMBOL
        && g->data.function.head->data.symbol
        && strcmp(g->data.function.head->data.symbol, "PolynomialGCD") == 0) {
        trivial = true;
    }
    if (trivial) {
        expr_free(g);
        Expr* result = eval_and_free(expr_new_function(
            expr_new_symbol("Times"),
            (Expr*[]){
                num,
                eval_and_free(expr_new_function(
                    expr_new_symbol("Power"),
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
        expr_new_symbol("Rule"),
        (Expr*[]){expr_new_symbol("Extension"), expr_copy((Expr*)alpha)}, 2);
    Expr* new_num_call = expr_new_function(
        expr_new_symbol("PolynomialQuotient"),
        (Expr*[]){num, expr_copy(g), expr_copy(var), ext_rule_a}, 4);
    Expr* new_num = evaluate(new_num_call);
    expr_free(new_num_call);

    Expr* ext_rule_b = expr_new_function(
        expr_new_symbol("Rule"),
        (Expr*[]){expr_new_symbol("Extension"), expr_copy((Expr*)alpha)}, 2);
    Expr* new_den_call = expr_new_function(
        expr_new_symbol("PolynomialQuotient"),
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
        expr_new_symbol("Power"),
        (Expr*[]){new_den, expr_new_integer(-1)}, 2));
    Expr* result = eval_and_free(expr_new_function(
        expr_new_symbol("Times"),
        (Expr*[]){new_num, inv_den}, 2));
    return result;
}

Expr* builtin_cancel(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;

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

    Expr* alpha_auto = NULL;
    QATower* auto_tower = NULL;
    if (!alpha && auto_flag) {
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
            Expr* tower_result = qa_cancel_with_tower(arg, auto_tower);
            if (tower_result) {
                qa_tower_free(auto_tower);
                return tower_result;
            }
            /* Tower path declined.  Try each tower generator in turn as
             * a single-α extension and pick the smallest result that
             * beats the input.  See builtin_together for rationale. */
            {
                bool is_sum = (arg->type == EXPR_FUNCTION
                               && arg->data.function.head
                               && arg->data.function.head->type == EXPR_SYMBOL
                               && arg->data.function.head->data.symbol == SYM_Plus);
                int64_t in_size = leaf_count_internal((Expr*)arg, true);
                Expr* best = NULL;
                int64_t best_size = -1;
                for (int i = 0; i < auto_tower->n; i++) {
                    const Expr* gen = auto_tower->alpha_renders[i];
                    if (!gen) continue;
                    Expr* cand = is_sum
                        ? together_recursive_ext(arg, gen)
                        : cancel_with_extension(arg, gen);
                    if (!cand) continue;
                    Expr* tg = expr_new_function(expr_new_symbol("Together"),
                        (Expr*[]){cand}, 1);
                    Expr* folded = evaluate(tg);
                    expr_free(tg);
                    if (!folded) continue;
                    int64_t lc = leaf_count_internal(folded, true);
                    if (lc < in_size && (best_size < 0 || lc < best_size)) {
                        if (best) expr_free(best);
                        best = folded;
                        best_size = lc;
                    } else {
                        expr_free(folded);
                    }
                }
                qa_tower_free(auto_tower);
                auto_tower = NULL;
                if (best) return best;
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
                Expr* call_lcm = expr_new_function(expr_new_symbol("PolynomialLCM"), (Expr*[]){expr_copy(lcm_den), expr_copy(dens[i])}, 2);
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
                new_nums[i] = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(nums[i]), q}, 2));
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
                    ret = eval_and_free(expr_new_function(expr_new_symbol("Plus"), plus_args, count));
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
            Expr* combined_num = eval_and_free(expr_new_function(expr_new_symbol("Plus"), new_nums, count));
            Expr* inv_den = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(lcm_den), expr_new_integer(-1)}, 2));
            Expr* combined_expr = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){combined_num, inv_den}, 2));

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

Expr* builtin_together(Expr* res) {
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

    Expr* alpha_auto = NULL;
    QATower* auto_tower = NULL;
    if (!alpha && auto_flag) {
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
            Expr* tower_result = qa_cancel_with_tower(arg, auto_tower);
            if (tower_result) {
                qa_tower_free(auto_tower);
                return tower_result;
            }
            /* Tower path declined (typically: Step 1 substitution into
             * γ-polynomial form blew up the leaf count past the safety
             * gate).  Before falling back to the generic no-extension
             * path — which can't see the algebraic relations between
             * generators and produces a combined-but-not-cancelled
             * fraction — try each tower generator in turn as a single-α
             * extension.  together_recursive_ext threads that one α
             * through PolynomialLCM/Quotient; the remaining generators
             * pass through as opaque coefficients, but the resulting
             * single-α cancellation is often enough to collapse the
             * fraction.  A final no-extension Together fold-up combines
             * like-coefficient terms that the single-α path leaves
             * separated. */
            {
                int64_t in_size = leaf_count_internal((Expr*)arg, true);
                Expr* best = NULL;
                int64_t best_size = -1;
                for (int i = 0; i < auto_tower->n; i++) {
                    const Expr* gen = auto_tower->alpha_renders[i];
                    if (!gen) continue;
                    Expr* cand = together_recursive_ext(arg, gen);
                    if (!cand) continue;
                    /* Final no-extension fold-up: combine like coefficient
                     * terms that the single-α path treats as polynomial-
                     * variable coefficients. */
                    Expr* tg = expr_new_function(expr_new_symbol("Together"),
                        (Expr*[]){cand}, 1);
                    Expr* folded = evaluate(tg);
                    expr_free(tg);
                    if (!folded) continue;
                    int64_t lc = leaf_count_internal(folded, true);
                    if (lc < in_size && (best_size < 0 || lc < best_size)) {
                        if (best) expr_free(best);
                        best = folded;
                        best_size = lc;
                    } else {
                        expr_free(folded);
                    }
                }
                qa_tower_free(auto_tower);
                auto_tower = NULL;
                if (best) return best;
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
