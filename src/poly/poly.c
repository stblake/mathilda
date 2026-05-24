/* ====================================================================
 * poly.c -- Polynomial algebra primitives.
 *
 * This module implements Mathilda's symbolic polynomial layer:
 *
 *   PolynomialQ        -- structural polynomial test
 *   Variables          -- collect free variables of an expression
 *   Coefficient        -- pattern-matching coefficient extraction
 *   CoefficientList    -- multi-variable coefficient table
 *   PolynomialGCD/LCM  -- multivariate gcd/lcm via subresultant PRS
 *   PolynomialQuotient -- (Euclidean) division quotient
 *   PolynomialRemainder-- (Euclidean) division remainder
 *   PolynomialMod      -- reduction mod a polynomial / integer / list
 *   PolynomialExtendedGCD -- extended Euclidean algorithm
 *   Collect            -- group terms by powers of a variable
 *   Decompose          -- functional decomposition of univariate poly
 *   HornerForm         -- nested Horner representation
 *   Resultant          -- Sylvester-matrix resultant
 *   Discriminant       -- (-1)^(n(n-1)/2)/a_n * Res(p, p')
 *
 * The file is organised top-down:
 *   1. Base/exponent decomposition helpers (BPList) used by Coefficient
 *      and Collect to break terms apart.
 *   2. Fast coefficient extraction helpers used internally to avoid
 *      the full evaluator pipeline (the hottest path in this file).
 *   3. The general Coefficient[] builtin (slow path with full pattern
 *      matching) and Variables/PolynomialQ.
 *   4. Polynomial division and content/primitive-part computations.
 *   5. Multivariate GCD/LCM driven by the recursive subresultant PRS.
 *   6. Collect, CoefficientList, Decompose, HornerForm, Resultant,
 *      Discriminant, PolynomialMod, PolynomialExtendedGCD.
 *
 * Memory convention (matches the rest of Mathilda):
 *   Every helper that returns Expr* hands fresh ownership to the
 *   caller. Built-in entry points (`builtin_*`) leave the input `res`
 *   alive -- the evaluator (or `internal_call_impl`) frees it after a
 *   non-NULL return.
 * ==================================================================== */

#include "internal.h"
#include "poly.h"
#include "expand.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "sort.h"
#include "replace.h"
#include "print.h"
#include "match.h"
#include "rationalize.h"
#include "sym_names.h"
#include "options.h"
#include "qa.h"
#include "qaupoly.h"
#include "qafactor.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

bool contains_any_symbol_from(Expr* expr, Expr* var);

/* Forward declaration for size-budget guards in poly_gcd_internal /
 * exact_poly_div etc. Defined further down with the subresultant
 * polynomial-remainder sequence helpers. */
static int64_t subres_leaf_count(Expr* e);

/* ------------------------------------------------------------------ */
/* BasePower decomposition                                            */
/*                                                                    */
/* A `BPList` represents a product as a flat list of (base, exponent) */
/* pairs. e.g. 3 * x^2 * y^a becomes                                   */
/*    {(3,1), (x,2), (y,a)}                                            */
/* This is the canonical form used by Coefficient[] and Collect[] to  */
/* identify which factors are "the variable's powers" versus the      */
/* numeric/symbolic coefficient.                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    Expr* base;
    Expr* exp;
} BasePower;

typedef struct {
    BasePower* data;
    size_t count;
    size_t capacity;
} BPList;

static void bp_init(BPList* l) {
    l->capacity = 8;
    l->count = 0;
    l->data = malloc(sizeof(BasePower) * l->capacity);
}

static void bp_add(BPList* l, Expr* base, Expr* exp) {
    for (size_t i = 0; i < l->count; i++) {
        if (expr_eq(l->data[i].base, base)) {
            Expr* new_exp = internal_plus((Expr*[]){expr_copy(l->data[i].exp), expr_copy(exp)}, 2);
            expr_free(l->data[i].exp);
            l->data[i].exp = new_exp;
            return;
        }
    }
    if (l->count == l->capacity) {
        l->capacity *= 2;
        l->data = realloc(l->data, sizeof(BasePower) * l->capacity);
    }
    l->data[l->count].base = expr_copy(base);
    l->data[l->count].exp = expr_copy(exp);
    l->count++;
}

static void bp_free(BPList* l) {
    for (size_t i = 0; i < l->count; i++) {
        expr_free(l->data[i].base);
        expr_free(l->data[i].exp);
    }
    free(l->data);
}

static void decompose_to_bp(Expr* e, BPList* l) {
    if (!e) return;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* head = e->data.function.head->data.symbol;
        if (head == SYM_Times) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) decompose_to_bp(e->data.function.args[i], l);
            return;
        }
        if (head == SYM_Power && e->data.function.arg_count == 2) {
            Expr* base = e->data.function.args[0];
            Expr* exp = e->data.function.args[1];
            if (exp->type == EXPR_INTEGER || is_rational(exp, NULL, NULL)) {
                bp_add(l, base, exp);
                return;
            }
            if (exp->type == EXPR_FUNCTION && exp->data.function.head->type == EXPR_SYMBOL &&
                exp->data.function.head->data.symbol == SYM_Plus) {
                for (size_t i = 0; i < exp->data.function.arg_count; i++) {
                    Expr* sub_exp = exp->data.function.args[i];
                    Expr* sub_p = internal_power((Expr*[]){expr_copy(base), expr_copy(sub_exp)}, 2);
                    decompose_to_bp(sub_p, l);
                    expr_free(sub_p);
                }
                return;
            }
        }
    }
    Expr* one = expr_new_integer(1);
    bp_add(l, e, one);
    expr_free(one);
}

static Expr* rebuild_from_bp(BPList* l) {
    if (l->count == 0) return expr_new_integer(1);
    Expr** args = malloc(sizeof(Expr*) * l->count);
    size_t count = 0;
    for (size_t i = 0; i < l->count; i++) {
        if (l->data[i].exp->type == EXPR_INTEGER) {
            int64_t exp = l->data[i].exp->data.integer;
            if (exp == 0) continue;
            if (exp == 1) args[count++] = expr_copy(l->data[i].base);
            else args[count++] = internal_power((Expr*[]){expr_copy(l->data[i].base), expr_copy(l->data[i].exp)}, 2);
        } else {
            args[count++] = internal_power((Expr*[]){expr_copy(l->data[i].base), expr_copy(l->data[i].exp)}, 2);
        }
    }
    if (count == 0) { free(args); return expr_new_integer(1); }
    if (count == 1) { Expr* res = args[0]; free(args); return res; }
    Expr* res = internal_times(args, count);
    free(args);
    return res;
}

/* ------------------------------------------------------------------ */
/* Fast coefficient extraction                                        */
/*                                                                    */
/* The general Coefficient[expr, var, n] builtin (`builtin_coefficient`*/
/* below) decomposes every term into a base/power list and uses       */
/* `get_k` to support compound monomial targets such as `x*y` or      */
/* `Sin[x]^2`. That generality is overkill for the polynomial code,   */
/* which always asks for the coefficient of an atomic variable raised */
/* to a non-negative integer power -- and which usually already has   */
/* the input fully expanded.                                          */
/*                                                                    */
/* `get_coeff_expanded` walks an already-expanded polynomial once,    */
/* computing the coefficient of var^n directly. This avoids three     */
/* layers of overhead per call:                                       */
/*   1. The evaluator pipeline (head/attr/listable/...) invoked by    */
/*      `evaluate(Coefficient[...])`.                                 */
/*   2. A redundant `expr_expand` call inside `builtin_coefficient`.  */
/*   3. The base/power bookkeeping in `decompose_to_bp` / `get_k`.    */
/*                                                                    */
/* Because polynomial division / GCD / coefficient list / resultant   */
/* repeatedly query coefficients of the same expanded polynomial,     */
/* this is one of the largest hot spots in the file.                  */
/* ------------------------------------------------------------------ */

/* True when `var` is a building block the fast path knows how to     */
/* recognise as a factor: any leaf, or a function call whose head is  */
/* not Plus/Times/Power. Compound monomials like `x*y` or `x^2` need  */
/* the slow path's pattern-matching logic.                             */
static bool var_is_atomic(Expr* var) {
    if (!var) return false;
    if (var->type == EXPR_SYMBOL || var->type == EXPR_INTEGER ||
        var->type == EXPR_REAL   || var->type == EXPR_BIGINT ||
        var->type == EXPR_STRING) return true;
    if (var->type == EXPR_FUNCTION && var->data.function.head->type == EXPR_SYMBOL) {
        const char* h = var->data.function.head->data.symbol;
        return h != SYM_Plus && h != SYM_Times && h != SYM_Power;
    }
    return false;
}

/* For a single multiplicand `f` and an atomic `var`, return the      */
/* exponent of var contributed by `f`:                                 */
/*    f == var               -> 1                                      */
/*    f == Power[var, k]     -> k       (k a non-negative integer)     */
/*    otherwise              -> 0       (f is treated as constant      */
/*                                       w.r.t. var)                   */
/* Returns -1 to signal "factor involves var with a non-integer or     */
/* negative exponent" -- the term is not polynomial in var, and the    */
/* caller should drop it from the coefficient sum.                     */
static int factor_var_exp(Expr* f, Expr* var) {
    if (expr_eq(f, var)) return 1;
    if (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL &&
        f->data.function.head->data.symbol == SYM_Power &&
        f->data.function.arg_count == 2 &&
        expr_eq(f->data.function.args[0], var)) {
        Expr* exp = f->data.function.args[1];
        if (exp->type == EXPR_INTEGER) {
            int64_t e = exp->data.integer;
            if (e < 0 || e > INT_MAX) return -1;
            return (int)e;
        }
        return -1;
    }
    return 0;
}

/* Walk `term` to compute its var-degree and accumulate every factor   */
/* that does not involve var into `*kept` (a growable array of fresh   */
/* copies). Recurses into nested Times nodes -- expanded expressions   */
/* sometimes carry Times[const, Times[...]] structure when expansion   */
/* finishes before the FLAT attribute has been re-applied. Returns     */
/* the accumulated var-degree, or -1 if `term` contains a non-integer  */
/* or negative power of var (i.e. is not polynomial in var).           */
static int collect_term_factors(Expr* term, Expr* var,
                                Expr*** kept, size_t* kept_count, size_t* kept_cap) {
    if (term->type == EXPR_FUNCTION && term->data.function.head->type == EXPR_SYMBOL &&
        term->data.function.head->data.symbol == SYM_Times) {
        int total_deg = 0;
        for (size_t j = 0; j < term->data.function.arg_count; j++) {
            int sub = collect_term_factors(term->data.function.args[j], var,
                                           kept, kept_count, kept_cap);
            if (sub < 0) return -1;
            total_deg += sub;
        }
        return total_deg;
    }
    int e = factor_var_exp(term, var);
    if (e < 0) return -1;
    if (e == 0) {
        if (*kept_count == *kept_cap) {
            *kept_cap = *kept_cap ? *kept_cap * 2 : 4;
            *kept = realloc(*kept, sizeof(Expr*) * (*kept_cap));
        }
        (*kept)[(*kept_count)++] = expr_copy(term);
    }
    return e;
}

/* Compute the contribution of a single expanded-polynomial term to    */
/* the coefficient of var^n. Returns NULL if the term contributes      */
/* nothing (its var-degree is not n, or it contains a non-polynomial   */
/* power of var). Otherwise returns a fresh Expr* the caller owns.     */
static Expr* term_coeff_or_null(Expr* term, Expr* var, int n) {
    Expr** kept = NULL;
    size_t kept_count = 0, kept_cap = 0;
    int deg = collect_term_factors(term, var, &kept, &kept_count, &kept_cap);
    if (deg < 0 || deg != n) {
        for (size_t k = 0; k < kept_count; k++) expr_free(kept[k]);
        free(kept);
        return NULL;
    }
    Expr* coeff;
    if (kept_count == 0)      coeff = expr_new_integer(1);
    else if (kept_count == 1) coeff = kept[0];
    else                       coeff = internal_times(kept, kept_count);
    free(kept);
    return coeff;
}

/* Return the coefficient of var^n in `expanded`, where `expanded` is  */
/* the result of a prior `expr_expand` call. `var` must be atomic      */
/* (caller's responsibility -- typically a symbol).                    */
static Expr* get_coeff_expanded(Expr* expanded, Expr* var, int n) {
    if (!expanded) return expr_new_integer(0);

    bool is_plus = (expanded->type == EXPR_FUNCTION &&
                    expanded->data.function.head->type == EXPR_SYMBOL &&
                    expanded->data.function.head->data.symbol == SYM_Plus);
    Expr** terms = is_plus ? expanded->data.function.args : &expanded;
    size_t term_count = is_plus ? expanded->data.function.arg_count : 1;

    Expr** matches = malloc(sizeof(Expr*) * term_count);
    size_t match_count = 0;
    for (size_t i = 0; i < term_count; i++) {
        Expr* c = term_coeff_or_null(terms[i], var, n);
        if (c) matches[match_count++] = c;
    }

    Expr* result;
    if (match_count == 0)      result = expr_new_integer(0);
    else if (match_count == 1) result = matches[0];
    else                        result = internal_plus(matches, match_count);
    free(matches);
    return result;
}

/* Bulk-extract every coefficient of `expanded` in var, for indices    */
/* 0..max_deg. Single pass over the terms (O(terms) tree walks) versus */
/* the (max_deg+1) walks performed by repeatedly calling               */
/* `get_coeff_expanded`. Each output slot is a freshly-allocated Expr* */
/* the caller must free; missing degrees are filled with Integer 0.    */
/* Returns false if `var` is not atomic.                                */
static bool get_all_coeffs_expanded(Expr* expanded, Expr* var, int max_deg,
                                    Expr*** out_coeffs) {
    if (!var_is_atomic(var) || max_deg < 0) return false;

    bool is_plus = (expanded && expanded->type == EXPR_FUNCTION &&
                    expanded->data.function.head->type == EXPR_SYMBOL &&
                    expanded->data.function.head->data.symbol == SYM_Plus);
    Expr** terms = is_plus ? expanded->data.function.args : &expanded;
    size_t term_count = expanded ? (is_plus ? expanded->data.function.arg_count : 1) : 0;

    /* Per-degree growable bucket of contributions. */
    Expr*** buckets = calloc((size_t)(max_deg + 1), sizeof(Expr**));
    size_t* bucket_count = calloc((size_t)(max_deg + 1), sizeof(size_t));
    size_t* bucket_cap   = calloc((size_t)(max_deg + 1), sizeof(size_t));

    for (size_t i = 0; i < term_count; i++) {
        Expr** kept = NULL;
        size_t kept_count = 0, kept_cap = 0;
        int deg = collect_term_factors(terms[i], var, &kept, &kept_count, &kept_cap);
        if (deg < 0 || deg > max_deg) {
            for (size_t k = 0; k < kept_count; k++) expr_free(kept[k]);
            free(kept);
            continue;
        }
        Expr* coeff;
        if (kept_count == 0)      coeff = expr_new_integer(1);
        else if (kept_count == 1) coeff = kept[0];
        else                       coeff = internal_times(kept, kept_count);
        free(kept);

        if (bucket_count[deg] == bucket_cap[deg]) {
            bucket_cap[deg] = bucket_cap[deg] ? bucket_cap[deg] * 2 : 4;
            buckets[deg] = realloc(buckets[deg], sizeof(Expr*) * bucket_cap[deg]);
        }
        buckets[deg][bucket_count[deg]++] = coeff;
    }

    Expr** out = malloc(sizeof(Expr*) * (size_t)(max_deg + 1));
    for (int d = 0; d <= max_deg; d++) {
        size_t bc = bucket_count[d];
        if (bc == 0)      out[d] = expr_new_integer(0);
        else if (bc == 1) out[d] = buckets[d][0];
        else              out[d] = internal_plus(buckets[d], bc);
        free(buckets[d]);
    }
    free(buckets);
    free(bucket_count);
    free(bucket_cap);

    *out_coeffs = out;
    return true;
}

/* ------------------------------------------------------------------ */
/* Slow-path Coefficient[] (full pattern-matching support)            */
/* ------------------------------------------------------------------ */

/* For a target monomial decomposed into its base/power list, count    */
/* how many copies (k) the term contains. Returns 0 when the term has  */
/* no coefficient contribution (target factors missing or matched a    */
/* fractional/negative number of times).                               */
static int64_t get_k(BPList* term_bp, BPList* target_bp) {
    if (target_bp->count == 0) return 0;
    int64_t k_val = -1;
    for (size_t i = 0; i < target_bp->count; i++) {
        Expr* t_base = target_bp->data[i].base;
        Expr* t_exp = target_bp->data[i].exp;
        bool found = false;
        for (size_t j = 0; j < term_bp->count; j++) {
            if (expr_eq(t_base, term_bp->data[j].base)) {
                found = true;
                Expr* trm_exp = term_bp->data[j].exp;
                if (t_exp->type == EXPR_INTEGER && trm_exp->type == EXPR_INTEGER) {
                    if (trm_exp->data.integer % t_exp->data.integer != 0) return 0;
                    int64_t ratio = trm_exp->data.integer / t_exp->data.integer;
                    if (ratio <= 0) return 0;
                    if (k_val == -1) k_val = ratio;
                    else if (k_val != ratio) return 0;
                } else if (expr_eq(t_exp, trm_exp)) {
                    if (k_val == -1) k_val = 1;
                    else if (k_val != 1) return 0;
                } else {
                    /* Symbolic-exponent ratio test (general algorithm):
                     * compute Cancel[trm_exp / t_exp] and accept the term
                     * iff the ratio is a positive integer >= 2. Handles
                     * Collect[a x^(2c) + b x^(2c), x^c] (ratio 2) and
                     * Collect[a x^(c+1) + b x^(2c+2), x^(c+1)] (ratio 2)
                     * uniformly without touching the integer fast path. */
                    Expr* div_args[2];
                    div_args[0] = expr_copy(trm_exp);
                    Expr* inv_args[2] = { expr_copy(t_exp), expr_new_integer(-1) };
                    div_args[1] = expr_new_function(expr_new_symbol("Power"), inv_args, 2);
                    Expr* div = expr_new_function(expr_new_symbol("Times"), div_args, 2);
                    Expr* canc_args[1] = { div };
                    Expr* ratio = internal_cancel(canc_args, 1);
                    bool ok = (ratio && ratio->type == EXPR_INTEGER &&
                               ratio->data.integer >= 2);
                    int64_t r_val = ok ? ratio->data.integer : 0;
                    if (ratio) expr_free(ratio);
                    if (!ok) return 0;
                    if (k_val == -1) k_val = r_val;
                    else if (k_val != r_val) return 0;
                }
                break;
            }
        }
        if (!found) return 0;
    }
    return k_val == -1 ? 0 : k_val;
}

/* Coefficient[expr, form, n=1]                                        */
/*                                                                     */
/* Mathematica-compatible coefficient extraction supporting compound   */
/* monomial forms (e.g. `Coefficient[expr, x*y]` or                    */
/* `Coefficient[expr, Sin[x]^2]`). Algorithm: expand once, decompose   */
/* the target into a base/power list, then for each summand check that */
/* every target factor appears with the same multiplier `k`. The       */
/* matching multiplier `k` is computed by `get_k`.                     */
/*                                                                     */
/* For the much more common atomic-variable case, internal callers     */
/* should use `get_coeff_expanded` directly to skip the BP machinery.  */
Expr* builtin_coefficient(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count < 2 || res->data.function.arg_count > 3)) return NULL;
    Expr* original_expr = res->data.function.args[0];
    Expr* form = res->data.function.args[1];
    int64_t n = (res->data.function.arg_count == 3 && res->data.function.args[2]->type == EXPR_INTEGER) ? res->data.function.args[2]->data.integer : 1;

    Expr* expanded = expr_expand(original_expr);
    if (!expanded) return NULL;

    BPList target_bp;
    bp_init(&target_bp);
    decompose_to_bp(form, &target_bp);

    size_t term_count = 1;
    Expr** terms = &expanded;
    if (expanded->type == EXPR_FUNCTION && expanded->data.function.head->type == EXPR_SYMBOL && expanded->data.function.head->data.symbol == SYM_Plus) {
        term_count = expanded->data.function.arg_count;
        terms = expanded->data.function.args;
    }

    Expr** coeffs = malloc(sizeof(Expr*) * term_count);
    size_t c_count = 0;

    for (size_t i = 0; i < term_count; i++) {
        Expr* term = terms[i];
        BPList term_bp;
        bp_init(&term_bp);
        decompose_to_bp(term, &term_bp);

        int64_t k = get_k(&term_bp, &target_bp);

        if (k == n) {
            if (n == 0) {
                coeffs[c_count++] = expr_copy(term);
            } else {
                size_t rem_count = 0;
                Expr** rem_args = malloc(sizeof(Expr*) * term_bp.count);
                for (size_t j = 0; j < term_bp.count; j++) {
                    Expr* base = term_bp.data[j].base;
                    Expr* exp = term_bp.data[j].exp;
                    bool matched = false;
                    for (size_t tj = 0; tj < target_bp.count; tj++) {
                        if (expr_eq(base, target_bp.data[tj].base)) {
                            matched = true; break;
                        }
                    }
                    if (!matched) {
                        if (exp->type == EXPR_INTEGER && exp->data.integer == 1) {
                            rem_args[rem_count++] = expr_copy(base);
                        } else {
                            rem_args[rem_count++] = internal_power((Expr*[]){expr_copy(base), expr_copy(exp)}, 2);
                        }
                    }
                }
                if (rem_count == 0) coeffs[c_count++] = expr_new_integer(1);
                else if (rem_count == 1) coeffs[c_count++] = rem_args[0];
                else coeffs[c_count++] = internal_times(rem_args, rem_count);
                free(rem_args);
            }
        }
        bp_free(&term_bp);
    }

    bp_free(&target_bp);
    expr_free(expanded);

    Expr* final_res;
    if (c_count == 0) final_res = expr_new_integer(0);
    else if (c_count == 1) final_res = coeffs[0];
    else final_res = internal_plus(coeffs, c_count);

    free(coeffs);
    return final_res;
}

/* ------------------------------------------------------------------ */
/* Constants & Predicates                                             */
/*                                                                    */
/* These helpers classify an expression as "constant w.r.t. some set  */
/* of variables", which drives both Variables[] (collect free atoms)  */
/* and PolynomialQ[] (structural polynomial test).                    */
/* ------------------------------------------------------------------ */

static bool contains_symbol(Expr* expr, const char* sym) {
    if (!expr) return false;
    if (expr->type == EXPR_SYMBOL) return strcmp(expr->data.symbol, sym) == 0;
    if (expr->type == EXPR_FUNCTION) {
        if (contains_symbol(expr->data.function.head, sym)) return true;
        for (size_t i = 0; i < expr->data.function.arg_count; i++)
            if (contains_symbol(expr->data.function.args[i], sym)) return true;
    }
    return false;
}

bool contains_any_symbol_from(Expr* expr, Expr* var) {
    if (!expr || !var) return false;
    if (var->type == EXPR_SYMBOL) return contains_symbol(expr, var->data.symbol);
    if (var->type == EXPR_FUNCTION) {
        if (contains_any_symbol_from(expr, var->data.function.head)) return true;
        for (size_t i = 0; i < var->data.function.arg_count; i++)
            if (contains_any_symbol_from(expr, var->data.function.args[i])) return true;
    }
    return false;
}

static bool is_constant_symbol(Expr* e) {
    if (e->type != EXPR_SYMBOL) return false;
    const char* s = e->data.symbol;
    return (s == SYM_Pi || s == SYM_E || s == SYM_I || 
            s == SYM_Infinity || s == SYM_ComplexInfinity);
}

static bool is_number(Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL) return true;
    if (is_rational(e, NULL, NULL)) return true;
    if (is_complex(e, NULL, NULL)) return true;
    return is_constant_symbol(e);
}

void collect_variables(Expr* e, Expr*** vars_ptr, size_t* count, size_t* capacity) {
    if (!e || is_number(e)) return;
    if (e->type == EXPR_FUNCTION) {
        const char* head = (e->data.function.head->type == EXPR_SYMBOL) ? e->data.function.head->data.symbol : "";
        if (strcmp(head, "Plus") == 0 || strcmp(head, "Times") == 0 || strcmp(head, "List") == 0) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) collect_variables(e->data.function.args[i], vars_ptr, count, capacity);
            return;
        }
        if (strcmp(head, "Power") == 0 && e->data.function.arg_count == 2) {
            Expr* base = e->data.function.args[0];
            if (is_number(base)) return;
            if (e->data.function.args[1]->type == EXPR_INTEGER || is_rational(e->data.function.args[1], NULL, NULL)) {
                collect_variables(base, vars_ptr, count, capacity); return;
            }
        }
    }
    for (size_t i = 0; i < *count; i++) if (expr_eq((*vars_ptr)[i], e)) return;
    if (*count == *capacity) { *capacity *= 2; *vars_ptr = realloc(*vars_ptr, sizeof(Expr*) * (*capacity)); }
    (*vars_ptr)[(*count)++] = expr_copy(e);
}

int compare_expr_ptrs(const void* a, const void* b) { return expr_compare(*(Expr**)a, *(Expr**)b); }

/* ====================================================================
 * Algebraic-generator substitution
 *
 * Detect a (base B, exponent atom A) pair such that the input is a
 * polynomial / rational function in Power[B, A] up to a common rational
 * scaling of exponents.  Substitute Power[B, A/m] -> g (g a fresh
 * symbol, m = lcm of the rational scalings' denominators) so all
 * matching Power sites become integer powers of g, let the caller run
 * a polynomial operation in g, then back-substitute g -> Power[B, A/m].
 *
 * Two flavours fall under this scheme:
 *
 *   * Radical case (A = 1, m > 1).  Power[B, p/q] -> g^(p*m/q).
 *     A bare occurrence of B (not wrapped in Power) counts as
 *     Power[B, 1] and substitutes to g^m.
 *
 *   * Exponential case (A != 1, multiple matching Power sites).
 *     For Power[B, c*A] (c rational), substitutes to g^(c*m).
 *
 * ==================================================================== */

/* Decompose an exponent expression `exp` into a leading rational
 * coefficient (c_n / c_d) and an "atom" sub-expression A, such that
 * exp = (c_n / c_d) * A.  Concretely:
 *   - If exp is itself rational, c = exp and A = 1.
 *   - If exp = Times[r, rest...] with r rational, c = r and
 *     A = Times[rest...]  (or the single rest arg if there is one).
 *   - Otherwise c = 1 and A = exp.
 * The returned A is always a fresh deep copy that the caller frees.
 */
static void decompose_exponent(Expr* exp, int64_t* c_n, int64_t* c_d, Expr** A_out) {
    int64_t p, q;
    if (is_rational(exp, &p, &q)) {
        *c_n = p; *c_d = q;
        *A_out = expr_new_integer(1);
        return;
    }
    if (exp->type == EXPR_FUNCTION &&
        exp->data.function.head->type == EXPR_SYMBOL &&
        exp->data.function.head->data.symbol == SYM_Times &&
        exp->data.function.arg_count >= 2 &&
        is_rational(exp->data.function.args[0], &p, &q)) {
        *c_n = p; *c_d = q;
        size_t n = exp->data.function.arg_count - 1;
        if (n == 1) {
            *A_out = expr_copy(exp->data.function.args[1]);
        } else {
            Expr** rest = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) {
                rest[i] = expr_copy(exp->data.function.args[i + 1]);
            }
            *A_out = eval_and_free(expr_new_function(expr_new_symbol("Times"), rest, n));
            free(rest);
        }
        return;
    }
    *c_n = 1; *c_d = 1;
    *A_out = expr_copy(exp);
}

/* True iff `e` is the integer 1. */
static bool atom_is_one(Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 1;
}

/* Find the first Power[B, exp] subterm with a non-trivial decomposition:
 *   - exp is rational with non-trivial denominator (radical case), or
 *   - exp has a non-trivial atom (exponential case).
 * Returns a fresh deep-copy of B in *base_out and of the atom in
 * *atom_out (always non-null on success). */
static bool walk_find_radical_base(Expr* e, Expr** base_out, Expr** atom_out) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        int64_t cn, cd;
        Expr* A = NULL;
        decompose_exponent(exp, &cn, &cd, &A);
        bool A_one = atom_is_one(A);
        if ((A_one && cd != 1) || !A_one) {
            *base_out = expr_copy(e->data.function.args[0]);
            *atom_out = A;
            return true;
        }
        expr_free(A);
    }
    if (walk_find_radical_base(e->data.function.head, base_out, atom_out)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (walk_find_radical_base(e->data.function.args[i], base_out, atom_out)) return true;
    }
    return false;
}

/* Test whether `e` is a "target Power" for the (B, A) pair: either
 *   - e == Power[B, c*A] with c rational, or
 *   - A is integer 1 and e == B (bare base, treated as c = 1).
 * On success, populate *c_n / *c_d with the rational coefficient. */
static bool is_target_power(Expr* e, Expr* B, Expr* A, int64_t* c_n, int64_t* c_d) {
    if (!e) return false;
    bool A_one = atom_is_one(A);
    if (A_one && expr_eq(e, B)) {
        *c_n = 1; *c_d = 1;
        return true;
    }
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2 &&
        expr_eq(e->data.function.args[0], B)) {
        Expr* exp = e->data.function.args[1];
        int64_t cn, cd;
        Expr* A_actual = NULL;
        decompose_exponent(exp, &cn, &cd, &A_actual);
        bool match = expr_eq(A_actual, A);
        expr_free(A_actual);
        if (match) {
            *c_n = cn;
            *c_d = cd;
            return true;
        }
    }
    return false;
}

/* Walk e: at every (B, A)-matching Power site update *m_out with lcm of
 * c-denominators, increment *count_out, and set *nontrivial_out true if
 * the c coefficient there is anything other than 1.  Also track
 * *varied_out: true if at least two matching sites have *different*
 * (cn, cd) pairs.  When varied_out stays false the substitution would
 * just rename one repeated power -> a single g^k and gain nothing,
 * while regressing the GCD path from Q[B^(1/m)][x] (B^(1/m) opaque) to
 * Q[g, x] bivariate (which can blow up via subresultant PRS). */
static void walk_gather(Expr* e, Expr* B, Expr* A,
                        int64_t* m_out, size_t* count_out,
                        bool* nontrivial_out,
                        bool* seen_out, int64_t* first_cn, int64_t* first_cd,
                        bool* varied_out) {
    if (!e) return;
    int64_t cn, cd;
    if (is_target_power(e, B, A, &cn, &cd)) {
        int64_t qa = cd < 0 ? -cd : cd;
        *m_out = lcm(*m_out, qa);
        (*count_out)++;
        if (!(cn == 1 && cd == 1)) *nontrivial_out = true;
        if (!*seen_out) {
            *seen_out = true;
            *first_cn = cn;
            *first_cd = cd;
        } else if (cn != *first_cn || cd != *first_cd) {
            *varied_out = true;
        }
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    walk_gather(e->data.function.head, B, A, m_out, count_out, nontrivial_out,
                seen_out, first_cn, first_cd, varied_out);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        walk_gather(e->data.function.args[i], B, A, m_out, count_out,
                    nontrivial_out, seen_out, first_cn, first_cd, varied_out);
    }
}

bool poly_find_radical_gen(Expr* e, Expr** base_out, Expr** atom_out, int64_t* m_out) {
    Expr* base = NULL;
    Expr* atom = NULL;
    if (!walk_find_radical_base(e, &base, &atom)) return false;
    int64_t m = 1;
    size_t count = 0;
    bool nontrivial_c = false;
    bool seen = false;
    int64_t first_cn = 0, first_cd = 0;
    bool varied = false;
    walk_gather(e, base, atom, &m, &count, &nontrivial_c,
                &seen, &first_cn, &first_cd, &varied);
    bool A_one = atom_is_one(atom);

    /* Detect whether the substitution would create a *bivariate* (or
     * higher) polynomial problem.  If `e` already has polynomial
     * variables besides the radical's base, replacing B^(1/m) with a
     * fresh symbol g introduces a second polynomial variable that
     * appears alongside the existing one(s) in every coefficient.
     * The downstream multivariate GCD can then enter the classical
     * subresultant-PRS coefficient explosion (#hang).  When the
     * substitution merely renames a single repeated power -> g^k
     * (varied=false) we gain nothing from making g a real variable —
     * better to leave B^(1/m) as an opaque coefficient so the GCD
     * stays in K[x] with K = Q[B^(1/m)]. */
    bool has_other_var = false;
    {
        size_t v_count = 0, v_cap = 8;
        Expr** vars = malloc(sizeof(Expr*) * v_cap);
        collect_variables(e, &vars, &v_count, &v_cap);
        for (size_t i = 0; i < v_count; i++) {
            if (!expr_eq(vars[i], base)) { has_other_var = true; break; }
        }
        for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
        free(vars);
    }

    /* Trigger conditions:
     *   - Radical case (A == 1): m > 1.  When the substitution would
     *     create a bivariate problem (has_other_var) AND no two sites
     *     actually combine (`varied` false), skip — substitution is a
     *     pure rename and only hurts.  Otherwise (univariate, or genuine
     *     combination) the substitution is beneficial.
     *   - Exponential case (A != 1): at least one Power[B, c*A] with
     *     c != 1, so g^c gives a non-trivial polynomial term.  When all
     *     matching sites have c = 1 the substitution g = Power[B, A]
     *     just renames Power[B, A] -> g and the operation gains nothing.
     */
    bool trigger = A_one ? (m > 1 && (!has_other_var || varied))
                         : (count >= 1 && nontrivial_c);
    if (!trigger) {
        expr_free(base); expr_free(atom);
        return false;
    }
    *base_out = base;
    *atom_out = atom;
    *m_out = m;
    return true;
}

Expr* poly_subst_radical_to_gen(Expr* e, Expr* base, Expr* atom, int64_t m, const char* gen) {
    if (!e) return NULL;
    int64_t cn, cd;
    if (is_target_power(e, base, atom, &cn, &cd)) {
        int64_t k = cn * (m / cd);
        return eval_and_free(expr_new_function(expr_new_symbol("Power"),
                  (Expr*[]){expr_new_symbol(gen), expr_new_integer(k)}, 2));
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    size_t count = e->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * count);

    /* Power[base, exp]: substitute the base only.  Exponents are
     * polynomial-degree literals (or symbolic) — never the place to
     * inject an algebraic-extension generator.  Without this guard
     * a polynomial like `-2 x^2 + 3 Sqrt[2] x` (radical-gen finds
     * base=2, atom=1) would have its bare integer `2` inside
     * Power[x, 2]'s exponent rewritten to gen^2, corrupting x^2 into
     * x^(gen^2). */
    bool is_power_node = (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Power
        && count == 2);

    for (size_t i = 0; i < count; i++) {
        if (is_power_node && i == 1) {
            new_args[i] = expr_copy(e->data.function.args[i]);
        } else {
            new_args[i] = poly_subst_radical_to_gen(e->data.function.args[i], base, atom, m, gen);
        }
    }
    Expr* new_head = poly_subst_radical_to_gen(e->data.function.head, base, atom, m, gen);
    Expr* result = eval_and_free(expr_new_function(new_head, new_args, count));
    free(new_args);
    return result;
}

/* Build Power[B, A * (k/m)] for given integer k. */
static Expr* build_back_power(Expr* base, Expr* atom, int64_t k, int64_t m) {
    if (atom_is_one(atom)) {
        Expr* exp = make_rational(k, m);
        return eval_and_free(expr_new_function(expr_new_symbol("Power"),
                  (Expr*[]){expr_copy(base), exp}, 2));
    }
    /* Power[B, (k/m)*A]. */
    Expr* coeff = make_rational(k, m);
    Expr* new_exp = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                  (Expr*[]){coeff, expr_copy(atom)}, 2));
    return eval_and_free(expr_new_function(expr_new_symbol("Power"),
              (Expr*[]){expr_copy(base), new_exp}, 2));
}

Expr* poly_subst_radical_from_gen(Expr* e, Expr* base, Expr* atom, int64_t m, const char* gen) {
    if (!e) return NULL;
    /* Bare gen: g -> Power[B, A/m]. */
    if (e->type == EXPR_SYMBOL && strcmp(e->data.symbol, gen) == 0) {
        return build_back_power(base, atom, 1, m);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    /* Power[gen, k] -> Power[B, k*A/m]. */
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_SYMBOL &&
        strcmp(e->data.function.args[0]->data.symbol, gen) == 0) {
        Expr* exp = e->data.function.args[1];
        int64_t pp, qq;
        if (is_rational(exp, &pp, &qq)) {
            /* g^(pp/qq) -> Power[B, A * (pp / (qq * m))]. */
            if (atom_is_one(atom)) {
                Expr* new_exp = make_rational(pp, qq * m);
                return eval_and_free(expr_new_function(expr_new_symbol("Power"),
                          (Expr*[]){expr_copy(base), new_exp}, 2));
            }
            Expr* coeff = make_rational(pp, qq * m);
            Expr* new_exp = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                      (Expr*[]){coeff, expr_copy(atom)}, 2));
            return eval_and_free(expr_new_function(expr_new_symbol("Power"),
                      (Expr*[]){expr_copy(base), new_exp}, 2));
        }
        /* g^exp where exp is symbolic — Power[B, exp*A/m]. */
        Expr* one_over_m = make_rational(1, m);
        Expr** factors;
        size_t nf;
        if (atom_is_one(atom)) {
            factors = malloc(sizeof(Expr*) * 2);
            factors[0] = expr_copy(exp);
            factors[1] = one_over_m;
            nf = 2;
        } else {
            factors = malloc(sizeof(Expr*) * 3);
            factors[0] = expr_copy(exp);
            factors[1] = one_over_m;
            factors[2] = expr_copy(atom);
            nf = 3;
        }
        Expr* new_exp = eval_and_free(expr_new_function(expr_new_symbol("Times"), factors, nf));
        free(factors);
        return eval_and_free(expr_new_function(expr_new_symbol("Power"),
                  (Expr*[]){expr_copy(base), new_exp}, 2));
    }
    size_t count = e->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * count);
    for (size_t i = 0; i < count; i++) {
        new_args[i] = poly_subst_radical_from_gen(e->data.function.args[i], base, atom, m, gen);
    }
    Expr* new_head = poly_subst_radical_from_gen(e->data.function.head, base, atom, m, gen);
    Expr* result = eval_and_free(expr_new_function(new_head, new_args, count));
    free(new_args);
    return result;
}

static bool walk_has_symbol_name(Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol, name) == 0;
    if (e->type == EXPR_FUNCTION) {
        if (walk_has_symbol_name(e->data.function.head, name)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (walk_has_symbol_name(e->data.function.args[i], name)) return true;
        }
    }
    return false;
}

char* poly_make_fresh_gen(Expr* e) {
    char buf[64];
    for (int n = 0; n < 1000; n++) {
        snprintf(buf, sizeof(buf), "$pc_radgen%d$", n);
        if (!walk_has_symbol_name(e, buf)) {
            size_t len = strlen(buf);
            char* out = malloc(len + 1);
            memcpy(out, buf, len + 1);
            return out;
        }
    }
    /* Fallback (collision-prone but extremely unlikely). */
    const char* fb = "$pc_radgen$";
    size_t len = strlen(fb);
    char* out = malloc(len + 1);
    memcpy(out, fb, len + 1);
    return out;
}

Expr* builtin_variables(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    size_t count = 0, capacity = 16;
    Expr** vars = malloc(sizeof(Expr*) * capacity);
    collect_variables(res->data.function.args[0], &vars, &count, &capacity);
    if (count > 0) qsort(vars, count, sizeof(Expr*), compare_expr_ptrs);
    Expr* list = expr_new_function(expr_new_symbol("List"), vars, count);
    free(vars); return list;
}

bool is_polynomial(Expr* e, Expr** vars, size_t var_count) {
    if (!e) return false;
    
    // Check if e is one of the variables
    for (size_t i = 0; i < var_count; i++) {
        if (expr_eq(e, vars[i])) return true;
    }
    
    // Check if e contains any of the variables.
    // If it DOES NOT contain any variable, it is a constant (polynomial of degree 0).
    bool contains_var = false;
    for (size_t i = 0; i < var_count; i++) {
        if (contains_any_symbol_from(e, vars[i])) {
            contains_var = true;
            break;
        }
    }
    if (!contains_var) return true;

    // If it contains a variable but didn't match expr_eq, we check its structure.
    if (e->type == EXPR_FUNCTION) {
        const char* head = (e->data.function.head->type == EXPR_SYMBOL) ? e->data.function.head->data.symbol : "";
        
        if (strcmp(head, "Plus") == 0 || strcmp(head, "Times") == 0) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (!is_polynomial(e->data.function.args[i], vars, var_count)) return false;
            }
            return true;
        }
        
        if (strcmp(head, "Power") == 0 && e->data.function.arg_count == 2) {
            Expr* base = e->data.function.args[0];
            Expr* exp = e->data.function.args[1];
            
            // Power[base, exp] is a polynomial if:
            // 1. exp is a non-negative integer AND base is a polynomial.
            if (exp->type == EXPR_INTEGER && exp->data.integer >= 0) {
                return is_polynomial(base, vars, var_count);
            }
            // Note: if base and exp are both constants (free of vars), 
            // it would have been caught by the !contains_var check above.
            return false;
        }
    }
    
    // If we reach here, it contains a variable but isn't a simple Plus/Times/Power 
    // and didn't match expr_eq. Thus it's not a polynomial in those variables.
    // e.g. Sin[x] is not a polynomial in x.
    return false;
}

Expr* builtin_polynomialq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* vars_expr = res->data.function.args[1];
    Expr** vars = NULL; size_t var_count = 0; bool free_vars = false;
    if (vars_expr->type == EXPR_FUNCTION && vars_expr->data.function.head->type == EXPR_SYMBOL && vars_expr->data.function.head->data.symbol == SYM_List) {
        var_count = vars_expr->data.function.arg_count;
        if (var_count > 0) { vars = malloc(sizeof(Expr*) * var_count); for (size_t i = 0; i < var_count; i++) vars[i] = vars_expr->data.function.args[i]; free_vars = true; }
    } else { vars = &res->data.function.args[1]; var_count = 1; }
    bool poly = is_polynomial(res->data.function.args[0], vars, var_count);
    if (free_vars) free(vars);
    return expr_new_symbol(poly ? "True" : "False");
}

/* ------------------------------------------------------------------ */
/* Polynomial GCD machinery                                           */
/*                                                                    */
/* The multivariate GCD is built bottom-up via the recursive          */
/* "subresultant PRS" algorithm:                                      */
/*                                                                    */
/*   poly_gcd_internal(A, B, vars, k) computes gcd(A, B) over the     */
/*   ring K[vars[0..k-1]], treating vars[k-1] as the main variable.   */
/*   Coefficient gcds in K[vars[0..k-2]] are computed by recursive    */
/*   calls. Eventually k==0 reduces to integer gcd.                   */
/*                                                                    */
/* Each level                                                         */
/*    1. Splits A and B into content (gcd of coefficients) and        */
/*       primitive part (poly / content).                             */
/*    2. Reduces the primitive parts with `pseudo_rem` until the      */
/*       remainder is zero.                                           */
/*    3. Multiplies the resulting gcd by the recursive content gcd.   */
/*                                                                    */
/* Pseudo-remainder is used (rather than the field-division remainder */
/* used by `poly_div_rem`) because it stays inside the coefficient    */
/* ring -- no rationals appear, even when the leading coefficient is  */
/* a complicated polynomial in the remaining variables.               */
/* ------------------------------------------------------------------ */

/* True if `e` is structurally / arithmetically the zero polynomial.   */
/* Cheap path covers literal 0; expanded path handles `(x-x)`-style    */
/* cancellations; deep path falls back to coefficient-list inspection. */
/* Depth-bounded helper: stops the recursion before overflowing the    */
/* C stack on inputs where CoefficientList does not strip a variable   */
/* and the recursive call sees the same expression again (e.g. when    */
/* `vars[0]` is an algebraic constant like Sqrt[5] and the polynomial  */
/* involves several Sqrt[..] mixed with x). Each recursive descent     */
/* genuinely strips one variable, so a small bound suffices for any    */
/* polynomial; an exhausted bound conservatively reports non-zero. */
#define IS_ZERO_POLY_MAX_DEPTH 32
static bool is_zero_poly_depth(Expr* e, int depth) {
    if (!e) return true;
    if (e->type == EXPR_INTEGER && e->data.integer == 0) return true;
    if (e->type == EXPR_REAL && e->data.real == 0.0) return true;
    if (depth >= IS_ZERO_POLY_MAX_DEPTH) return false;

    Expr* expanded = expr_expand(e);
    bool res = false;
    if (expanded->type == EXPR_INTEGER && expanded->data.integer == 0) res = true;
    else if (expanded->type == EXPR_REAL && expanded->data.real == 0.0) res = true;

    if (!res) {
        size_t v_count = 0, v_cap = 16;
        Expr** vars = malloc(sizeof(Expr*) * v_cap);
        collect_variables(expanded, &vars, &v_count, &v_cap);
        if (v_count > 0) {
            if (is_polynomial(expanded, vars, v_count)) {
                Expr* var = vars[0];
                Expr* clist = internal_coefficientlist((Expr*[]){expr_copy(expanded), expr_copy(var)}, 2);
                if (clist && clist->type == EXPR_FUNCTION &&
                    clist->data.function.head->type == EXPR_SYMBOL &&
                    clist->data.function.head->data.symbol == SYM_List) {
                    bool all_zero = true;
                    for (size_t i = 0; i < clist->data.function.arg_count; i++) {
                        if (!is_zero_poly_depth(clist->data.function.args[i], depth + 1)) {
                            all_zero = false;
                            break;
                        }
                    }
                    if (all_zero) res = true;
                }
                if (clist) expr_free(clist);
            }
        }
        for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
        free(vars);
    }
    expr_free(expanded);
    return res;
}

bool is_zero_poly(Expr* e) {
    return is_zero_poly_depth(e, 0);
}

static bool is_negative(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER && e->data.integer < 0) return true;
    if (e->type == EXPR_REAL && e->data.real < 0.0) return true;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        if (e->data.function.head->data.symbol == SYM_Times && e->data.function.arg_count > 0) {
            return is_negative(e->data.function.args[0]);
        }
        if (e->data.function.head->data.symbol == SYM_Rational && e->data.function.arg_count == 2) {
            return is_negative(e->data.function.args[0]);
        }
    }
    return false;
}

static int64_t get_int_content(Expr* e) {
    if (!e) return 1;
    if (e->type == EXPR_INTEGER) return llabs(e->data.integer);
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Complex && e->data.function.arg_count == 2) {
        return gcd(get_int_content(e->data.function.args[0]), get_int_content(e->data.function.args[1]));
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Times) {
        int64_t c = 1;
        for (size_t i=0; i<e->data.function.arg_count; i++) {
            c *= get_int_content(e->data.function.args[i]);
        }
        return c;
    }
    /* Plus content: GCD of summand contents.  Without this, the numG
     * accumulation in builtin_polynomialgcd treats the entire Plus as
     * a single non-numeric base with content 1, so e.g.
     *   PolynomialGCD[4, 4 + 4 x]
     * returned 1 instead of 4 -- and Cancel[4/(4+4x)] failed to
     * reduce. */
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Plus) {
        int64_t g = 0;  /* gcd-identity */
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            int64_t ci = get_int_content(e->data.function.args[i]);
            g = gcd(g, ci);
            if (g == 1) return 1;  /* short-circuit: cannot grow */
        }
        return g == 0 ? 1 : g;
    }
    return 1;
}

static Expr* my_number_gcd(Expr* a, Expr* b) {
    int64_t ca = get_int_content(a);
    int64_t cb = get_int_content(b);
    int64_t g = gcd(ca, cb);
    if (g == 0) g = 1;
    return expr_new_integer(g);
}

/* Degree of `e` as a polynomial in `var`. Walks the expression tree   */
/* once; returns 0 for constants, max-of-summands for Plus, sum-of-    */
/* factors for Times, and the integer exponent for `var^k`.            */
int get_degree_poly(Expr* e, Expr* var) {
    if (!e) return 0;
    if (expr_eq(e, var)) return 1;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* head = e->data.function.head->data.symbol;
        if (head == SYM_Power && e->data.function.arg_count == 2) {
            if (expr_eq(e->data.function.args[0], var)) {
                if (e->data.function.args[1]->type == EXPR_INTEGER) {
                    return (int)e->data.function.args[1]->data.integer;
                }
            }
        } else if (head == SYM_Plus) {
            int max_deg = 0;
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                int d = get_degree_poly(e->data.function.args[i], var);
                if (d > max_deg) max_deg = d;
            }
            return max_deg;
        } else if (head == SYM_Times) {
            int sum_deg = 0;
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                sum_deg += get_degree_poly(e->data.function.args[i], var);
            }
            return sum_deg;
        }
    }
    return 0;
}

/* Public coefficient extractor used across the polynomial code and    */
/* by external modules (facpoly.c, parfrac.c, simp.c, ...). The fast   */
/* path expands once and walks the result; the slow path defers to     */
/* the full Coefficient[] builtin via evaluate() so that compound      */
/* monomial targets (`x*y`, `Sin[x]^2`, ...) keep working.             */
Expr* get_coeff(Expr* e, Expr* var, int d) {
    if (var_is_atomic(var)) {
        Expr* expanded = expr_expand(e);
        if (expanded) {
            Expr* result = get_coeff_expanded(expanded, var, d);
            expr_free(expanded);
            return result;
        }
    }
    Expr* call = expr_new_function(expr_new_symbol("Coefficient"),
                                   (Expr*[]){expr_copy(e), expr_copy(var), expr_new_integer(d)}, 3);
    Expr* res = evaluate(call);
    expr_free(call);
    return res;
}

/* True iff e is rational-like (Integer, BigInt, Rational[n,d]) or
 * Complex[r, i] with both components rational-like. Field elements of
 * Q or Q[i] qualify; algebraic non-field atoms like Sqrt[2] do not. */
static bool is_rational_or_gaussian(const Expr* e) {
    if (!e) return false;
    if (is_rational_like((Expr*)e)) return true;
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Complex &&
        e->data.function.arg_count == 2) {
        return is_rational_like(e->data.function.args[0]) &&
               is_rational_like(e->data.function.args[1]);
    }
    return false;
}

/* Returns Q such that A = B * Q exactly (in K[vars]); otherwise NULL. */
/* Used by `poly_gcd_internal` to extract primitive parts and to stage */
/* GCD candidate division. The base case (var_count == 0) reduces to   */
/* exact integer / rational division -- everything else recurses on    */
/* the leading-coefficient ring.                                       */
Expr* exact_poly_div(Expr* A, Expr* B, Expr** vars, size_t var_count) {
    if (expr_eq(A, B)) return expr_new_integer(1);
    if (var_count == 0) {
        if (is_zero_poly(B)) return NULL;

        /* B is the multiplicative identity: A/1 = A. Always exact for
         * any A, including non-field atoms like Sqrt[2]. */
        if (B->type == EXPR_INTEGER && B->data.integer == 1) {
            return expr_copy(A);
        }
        /* A is zero: 0/B = 0 always (B != 0 already checked). */
        if (A->type == EXPR_INTEGER && A->data.integer == 0) {
            return expr_new_integer(0);
        }

        if (A->type == EXPR_BIGINT && B->type == EXPR_BIGINT) {
            mpz_t a, b, r, rem;
            expr_to_mpz(A, a);
            expr_to_mpz(B, b);
            mpz_init(r);
            mpz_init(rem);
            mpz_tdiv_qr(r, rem, a, b);
            bool exact = (mpz_cmp_ui(rem, 0) == 0);
            mpz_clear(a); mpz_clear(b); mpz_clear(rem);
            if (exact) {
                Expr* res = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
                mpz_clear(r);
                return res;
            }
            mpz_clear(r);
            // Fall through to rational quotient
        }

        /* Rational/Gaussian fallback: only when BOTH operands lie in a
         * field (Q or Q[i]). Then A/B is an exact field element and the
         * symbolic Times[A, B^{-1}] auto-evaluates to a clean Rational
         * or Complex[Rational, Rational]. For non-rational atoms
         * (Sqrt[2], symbolic radicals, etc.) we are NOT in a field —
         * Q[Sqrt[2], Sqrt[3], ...] is a polynomial ring, not a field —
         * and a "fallback" symbolic Times[A, 1/B] is unsound: it claims
         * an exact division that does not exist in the working ring,
         * and propagating it through the GCD/LCM machinery triggers
         * multivariate Euclidean coefficient explosion (case-13 Together
         * hang). Return NULL to signal "no exact division". */
        if (is_rational_or_gaussian(A) && is_rational_or_gaussian(B)) {
            return internal_times((Expr*[]){expr_copy(A), internal_power((Expr*[]){expr_copy(B), expr_new_integer(-1)}, 2)}, 2);
        }
        return NULL;
    }
    
    Expr* x = vars[var_count - 1];
    Expr* expandedA = expr_expand(A);
    Expr* expandedB = expr_expand(B);
    int degA = get_degree_poly(expandedA, x);
    int degB = get_degree_poly(expandedB, x);
    
    if (degA < degB || is_zero_poly(expandedB)) {
        expr_free(expandedA); expr_free(expandedB); return NULL;
    }
    
    Expr* Q = expr_new_integer(0);
    Expr* R = expandedA;
    /* expandedA / expandedB are already expanded; skip the redundant      */
    /* expr_expand inside get_coeff by calling the direct helper.          */
    Expr* lcB = get_coeff_expanded(expandedB, x, degB);

    while (true) {
        int degR = get_degree_poly(R, x);

        if (degR < degB || is_zero_poly(R)) break;

        Expr* lcR = get_coeff_expanded(R, x, degR);
        int d = degR - degB;

        Expr* q_coeff = exact_poly_div(lcR, lcB, vars, var_count - 1);
        if (!q_coeff) {
            expr_free(lcR); break;
        }
        /* If the quotient coefficient is a pure integer/bigint, the      */
        /* subtraction below cannot introduce new fractions and we can    */
        /* skip the costly together()+cancel() unification.                */
        bool q_coeff_pure = (q_coeff->type == EXPR_INTEGER || q_coeff->type == EXPR_BIGINT);

        Expr* x_pow = internal_power((Expr*[]){expr_copy(x), expr_new_integer(d)}, 2);
        Expr* term = internal_times((Expr*[]){q_coeff, x_pow}, 2);

        Expr* new_Q = internal_plus((Expr*[]){expr_copy(Q), expr_copy(term)}, 2);
        expr_free(Q);
        Q = new_Q;

        Expr* term_B = internal_times((Expr*[]){term, expr_copy(expandedB)}, 2);
        Expr* neg_term_B = internal_times((Expr*[]){expr_new_integer(-1), term_B}, 2);
        Expr* diff = internal_plus((Expr*[]){expr_copy(R), neg_term_B}, 2);
        Expr* new_R;
        if (q_coeff_pure) {
            new_R = expr_expand(diff);
            expr_free(diff);
        } else {
            Expr* together = internal_together((Expr*[]){diff}, 1);
            Expr* cancelled = internal_cancel((Expr*[]){together}, 1);
            new_R = expr_expand(cancelled);
            expr_free(cancelled);
        }

        expr_free(R);
        R = new_R;
        expr_free(lcR);
    }
    expr_free(lcB);
    expr_free(expandedB);
    
    if (!is_zero_poly(R)) {
        expr_free(R); expr_free(Q); return NULL;
    }
    expr_free(R);
    return Q;
}

/* Subresultant-style pseudo-remainder used by the multivariate GCD     */
/* loop. Computes lc(B)^k * A mod B in the polynomial ring without     */
/* needing leading-coefficient division, so it stays inside the        */
/* coefficient ring (no rationals introduced) and works over symbolic  */
/* coefficients.                                                       */
static Expr* pseudo_rem(Expr* A, Expr* B, Expr* x) {
    Expr* R = expr_expand(A);
    Expr* expandedB = expr_expand(B);
    int degB = get_degree_poly(expandedB, x);
    /* expandedB is already expanded -- direct fast path is safe. */
    Expr* lcB = get_coeff_expanded(expandedB, x, degB);

    while (true) {
        int degR = get_degree_poly(R, x);

        if (degR < degB || is_zero_poly(R)) break;

        Expr* lcR = get_coeff_expanded(R, x, degR);
        int d = degR - degB;
        
        Expr* t1 = internal_times((Expr*[]){expr_copy(lcB), R}, 2);
        Expr* x_pow = internal_power((Expr*[]){expr_copy(x), expr_new_integer(d)}, 2);
        Expr* t2 = internal_times((Expr*[]){lcR, x_pow, expr_copy(expandedB)}, 3);
        Expr* neg_t2 = internal_times((Expr*[]){expr_new_integer(-1), t2}, 2);
        Expr* diff = internal_plus((Expr*[]){t1, neg_t2}, 2);
        
        R = expr_expand(diff);
        expr_free(diff);
    }
    expr_free(lcB);
    expr_free(expandedB);
    return R;
}

/* Univariate polynomial long division of `p` by `q` in K[x] (K is the */
/* field generated by the symbolic coefficients). Writes the quotient  */
/* and remainder via the out parameters; both are NULL on failure      */
/* (e.g. divisor is zero). Used as the building block for both         */
/* PolynomialQuotient[] and PolynomialRemainder[].                     */
static void poly_div_rem(Expr* p, Expr* q, Expr* x, Expr** out_Q, Expr** out_R) {
    Expr* expandedA = expr_expand(p);
    Expr* expandedB = expr_expand(q);
    int degB = get_degree_poly(expandedB, x);

    if (is_zero_poly(expandedB)) {
        *out_Q = NULL;
        *out_R = NULL;
        expr_free(expandedA);
        expr_free(expandedB);
        return;
    }

    if (degB == 0) {
        Expr* invB = internal_power((Expr*[]){expr_copy(expandedB), expr_new_integer(-1)}, 2);
        *out_Q = internal_expand((Expr*[]){internal_times((Expr*[]){expr_copy(expandedA), invB}, 2)}, 1);
        *out_R = expr_new_integer(0);
        expr_free(expandedA);
        expr_free(expandedB);
        return;
    }

    Expr* Q = expr_new_integer(0);
    Expr* R = expandedA;
    /* Fast path: expandedB is already expanded -- skip the duplicate    */
    /* expr_expand inside get_coeff.                                     */
    Expr* lcB = get_coeff_expanded(expandedB, x, degB);

    while (true) {
        int degR = get_degree_poly(R, x);
        if (degR < degB || is_zero_poly(R)) break;

        Expr* lcR = get_coeff_expanded(R, x, degR);
        int d = degR - degB;

        /* Compute the next quotient coefficient = lcR / lcB.             */
        /* `q_coeff_pure` tracks whether the coefficient is a pure        */
        /* integer / bigint -- in that case the subtraction step cannot   */
        /* introduce new fractions and we can skip the expensive          */
        /* together()+cancel() unification of denominators.                */
        Expr* q_coeff;
        bool q_coeff_pure = false;
        if (expr_eq(lcR, lcB)) {
            q_coeff = expr_new_integer(1);
            q_coeff_pure = true;
        } else if (lcR->type == EXPR_BIGINT && lcB->type == EXPR_BIGINT) {
            mpz_t a, b, r, rem;
            expr_to_mpz(lcR, a);
            expr_to_mpz(lcB, b);
            mpz_init(r);
            mpz_init(rem);
            mpz_tdiv_qr(r, rem, a, b);
            bool exact = (mpz_cmp_ui(rem, 0) == 0);
            mpz_clear(a); mpz_clear(b); mpz_clear(rem);
            if (exact) {
                q_coeff = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
                mpz_clear(r);
                q_coeff_pure = true;
            } else {
                mpz_clear(r);
                Expr* lcB_inv = internal_power((Expr*[]){expr_copy(lcB), expr_new_integer(-1)}, 2);
                q_coeff = internal_times((Expr*[]){expr_copy(lcR), lcB_inv}, 2);
            }
        } else {
            Expr* lcB_inv = internal_power((Expr*[]){expr_copy(lcB), expr_new_integer(-1)}, 2);
            q_coeff = internal_times((Expr*[]){expr_copy(lcR), lcB_inv}, 2);
        }
        expr_free(lcR);

        Expr* x_pow = internal_power((Expr*[]){expr_copy(x), expr_new_integer(d)}, 2);
        Expr* term = internal_times((Expr*[]){q_coeff, x_pow}, 2);

        Expr* new_Q = internal_plus((Expr*[]){expr_copy(Q), expr_copy(term)}, 2);
        expr_free(Q);
        Q = new_Q;

        Expr* term_B = internal_times((Expr*[]){term, expr_copy(expandedB)}, 2);
        Expr* neg_term_B = internal_times((Expr*[]){expr_new_integer(-1), term_B}, 2);
        Expr* diff = internal_plus((Expr*[]){expr_copy(R), neg_term_B}, 2);
        Expr* new_R;
        if (q_coeff_pure) {
            /* No new rational structure introduced; expand alone is      */
            /* enough to combine like terms.                              */
            new_R = expr_expand(diff);
            expr_free(diff);
        } else {
            /* Symbolic / rational q_coeff -- run together()+cancel()      */
            /* to unify denominators introduced by Power[lcB,-1].          */
            Expr* together = internal_together((Expr*[]){diff}, 1);
            Expr* cancelled = internal_cancel((Expr*[]){together}, 1);
            new_R = expr_expand(cancelled);
            expr_free(cancelled);
        }

        expr_free(R);
        R = new_R;
    }
    expr_free(lcB);
    expr_free(expandedB);

    *out_Q = Q;
    *out_R = R;
}

/* PolynomialQuotient[a, b, x, Extension -> α] /                        */
/* PolynomialRemainder[a, b, x, Extension -> α].                        */
/*                                                                      */
/* Lifts `a` and `b` into Q(α)[x] via qa_expr_to_qaupoly and runs       */
/* qaupoly_divrem.  Returns the quotient (which==0) or remainder        */
/* (which==1) rendered back as a Mathilda Expr, or NULL if the lift      */
/* fails (caller falls back to the standard Q-treats-α-as-opaque path). */
/*                                                                      */
/* Without this path, PolynomialQuotient[(x - Sqrt[2])(x + Sqrt[2]), …] */
/* runs poly_div_rem with multivariate Q[Sqrt[2], x] coefficients,      */
/* whose Together+Cancel reconciliation step inside the inner loop is   */
/* exponentially slow on hard inputs (the slow path that previously     */
/* stalled together_recursive_ext).                                      */
static Expr* polynomialdivrem_with_extension(Expr* p, Expr* q, Expr* x,
                                             const Expr* alpha, int which) {
    if (!alpha) return NULL;
    if (x->type != EXPR_SYMBOL) return NULL;

    Expr* alpha_render = NULL;
    QAExt* ext = qa_resolve_extension(alpha, &alpha_render);
    if (!ext) return NULL;

    QAUPoly* pp = qa_expr_to_qaupoly(p, x, alpha_render, ext);
    QAUPoly* qq = pp ? qa_expr_to_qaupoly(q, x, alpha_render, ext) : NULL;

    Expr* result = NULL;
    if (pp && qq && !qaupoly_is_zero(qq)) {
        QAUPoly* Q = NULL;
        QAUPoly* R = NULL;
        if (qaupoly_divrem(pp, qq, &Q, &R)) {
            QAUPoly* out = (which == 0) ? Q : R;
            if (qaupoly_is_zero(out)) {
                result = expr_new_integer(0);
            } else {
                result = qaupoly_to_expr_alpha(out, x->data.symbol,
                                               alpha_render);
            }
        }
        if (Q) qaupoly_free(Q);
        if (R) qaupoly_free(R);
    }

    if (pp) qaupoly_free(pp);
    if (qq) qaupoly_free(qq);
    if (alpha_render) expr_free(alpha_render);
    qaext_free(ext);
    return result;
}

Expr* builtin_polynomialquotient(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 3) return NULL;

    /* Strip a trailing Extension -> α option, if any.  When the user
     * passed `Extension -> Automatic`, run extension_autodetect on the
     * two polynomial arguments (single-generator only — tier-1). */
    size_t poly_argc = res->data.function.arg_count;
    bool auto_flag = false;
    const Expr* alpha = extract_extension_option_full(res, &poly_argc, &auto_flag);
    if (poly_argc != 3) return NULL;

    Expr* p = res->data.function.args[0];
    Expr* q = res->data.function.args[1];
    Expr* x = res->data.function.args[2];

    Expr* alpha_auto = NULL;
    QATower* auto_tower = NULL;
    if (!alpha && auto_flag) {
        Expr* gen_args[2] = { p, q };
        auto_tower = extension_autodetect_args(gen_args, 2);
        if (auto_tower && auto_tower->n == 1) {
            alpha_auto = expr_copy(auto_tower->alpha_renders[0]);
            alpha = alpha_auto;
        }
    }

    if (alpha) {
        Expr* ext_result = polynomialdivrem_with_extension(p, q, x, alpha, 0);
        if (ext_result) {
            if (alpha_auto) expr_free(alpha_auto);
            if (auto_tower) qa_tower_free(auto_tower);
            return ext_result;
        }
        /* Fall through to non-extension path on lift failure (e.g.
         * coefficients outside Q(α), unrecognised α, multivariate). */
        if (alpha_auto) { expr_free(alpha_auto); alpha_auto = NULL; }
        if (auto_tower) { qa_tower_free(auto_tower); auto_tower = NULL; }
    }

    Expr *Q, *R;
    poly_div_rem(p, q, x, &Q, &R);
    if (!Q) return NULL;
    expr_free(R);
    Expr* expanded_Q = expr_expand(Q);
    expr_free(Q);
    return expanded_Q;
}

Expr* builtin_polynomialremainder(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 3) return NULL;

    /* Strip a trailing Extension -> α option, if any. */
    size_t poly_argc = res->data.function.arg_count;
    bool auto_flag = false;
    const Expr* alpha = extract_extension_option_full(res, &poly_argc, &auto_flag);
    if (poly_argc != 3) return NULL;

    Expr* p = res->data.function.args[0];
    Expr* q = res->data.function.args[1];
    Expr* x = res->data.function.args[2];

    Expr* alpha_auto = NULL;
    QATower* auto_tower = NULL;
    if (!alpha && auto_flag) {
        Expr* gen_args[2] = { p, q };
        auto_tower = extension_autodetect_args(gen_args, 2);
        if (auto_tower && auto_tower->n == 1) {
            alpha_auto = expr_copy(auto_tower->alpha_renders[0]);
            alpha = alpha_auto;
        }
    }

    if (alpha) {
        Expr* ext_result = polynomialdivrem_with_extension(p, q, x, alpha, 1);
        if (ext_result) {
            if (alpha_auto) expr_free(alpha_auto);
            if (auto_tower) qa_tower_free(auto_tower);
            return ext_result;
        }
        if (alpha_auto) { expr_free(alpha_auto); alpha_auto = NULL; }
        if (auto_tower) { qa_tower_free(auto_tower); auto_tower = NULL; }
    }

    Expr *Q, *R;
    poly_div_rem(p, q, x, &Q, &R);
    if (!R) return NULL;
    expr_free(Q);
    return R;
}

/* PolynomialQuotientRemainder[p, q, x] returns {Q, R} such that          */
/* p == Q*q + R with deg(R) < deg(q).  Two-output companion to            */
/* PolynomialQuotient / PolynomialRemainder; required by the              */
/* BronsteinRational pipeline (ExtendedEuclidean).                         */
Expr* builtin_polynomialquotientremainder(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 3) return NULL;

    size_t poly_argc = res->data.function.arg_count;
    bool auto_flag = false;
    const Expr* alpha = extract_extension_option_full(res, &poly_argc, &auto_flag);
    if (poly_argc != 3) return NULL;

    Expr* p = res->data.function.args[0];
    Expr* q = res->data.function.args[1];
    Expr* x = res->data.function.args[2];

    Expr* alpha_auto = NULL;
    QATower* auto_tower = NULL;
    if (!alpha && auto_flag) {
        Expr* gen_args[2] = { p, q };
        auto_tower = extension_autodetect_args(gen_args, 2);
        if (auto_tower && auto_tower->n == 1) {
            alpha_auto = expr_copy(auto_tower->alpha_renders[0]);
            alpha = alpha_auto;
        }
    }

    if (alpha) {
        /* For the extension path we run the divrem twice -- the slow        */
        /* path keeps the implementation simple and still avoids the         */
        /* multivariate stall that the no-Extension path would trigger.     */
        Expr* qext = polynomialdivrem_with_extension(p, q, x, alpha, 0);
        Expr* rext = polynomialdivrem_with_extension(p, q, x, alpha, 1);
        if (qext && rext) {
            Expr* list = expr_new_function(expr_new_symbol("List"),
                                           (Expr*[]){qext, rext}, 2);
            if (alpha_auto) expr_free(alpha_auto);
            if (auto_tower) qa_tower_free(auto_tower);
            return list;
        }
        if (qext) expr_free(qext);
        if (rext) expr_free(rext);
        if (alpha_auto) { expr_free(alpha_auto); alpha_auto = NULL; }
        if (auto_tower) { qa_tower_free(auto_tower); auto_tower = NULL; }
        /* Fall through to plain path on lift failure. */
    }

    Expr *Q, *R;
    poly_div_rem(p, q, x, &Q, &R);
    if (!Q || !R) {
        if (Q) expr_free(Q);
        if (R) expr_free(R);
        return NULL;
    }
    Expr* expanded_Q = expr_expand(Q);
    expr_free(Q);
    Expr* list = expr_new_function(expr_new_symbol("List"),
                                   (Expr*[]){expanded_Q, R}, 2);
    return list;
}

/* SubresultantPolynomialRemainders[a, b, x] returns the polynomial-       */
/* remainder chain {a, b, R_2, R_3, ...} in K(coeffs)[x], iterating         */
/* pseudo-remainder until a constant or zero remainder is reached.          */
/*                                                                          */
/* Note: this is *not* the Lazard-scaled subresultant chain.  For the      */
/* BronsteinRational / Lazard-Rioboo-Trager use case the algorithm only     */
/* consumes (i) the degree of each chain element in x and (ii) the         */
/* primitive part (in the auxiliary variable t) of each element.  Both      */
/* properties are invariant under content scaling, so the simpler          */
/* pseudo-rem chain is a correct substrate without the Lazard arithmetic.  */
Expr* builtin_subresultantpolynomialremainders(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;

    Expr* A = res->data.function.args[0];
    Expr* B = res->data.function.args[1];
    Expr* x = res->data.function.args[2];

    if (x->type != EXPR_SYMBOL) return NULL;

    Expr* a0 = expr_expand(A);
    Expr* b0 = expr_expand(B);

    /* Ensure deg(a0) >= deg(b0) per the standard chain orientation.
     * Mathematica's SubresultantPolynomialRemainders documents that the
     * first argument has higher degree; if not, swap so the chain
     * starts properly. */
    int dA = get_degree_poly(a0, x);
    int dB = get_degree_poly(b0, x);
    if (dA < dB) {
        Expr* t = a0; a0 = b0; b0 = t;
    }

    size_t cap = 8, n = 0;
    Expr** chain = (Expr**)malloc(sizeof(Expr*) * cap);
    chain[n++] = a0;
    chain[n++] = b0;

    while (true) {
        Expr* prev = chain[n - 2];
        Expr* cur  = chain[n - 1];
        if (is_zero_poly(cur)) break;
        int dcur = get_degree_poly(cur, x);
        if (dcur <= 0) break;
        Expr* r = pseudo_rem(prev, cur, x);
        if (is_zero_poly(r)) { expr_free(r); break; }
        if (n >= cap) { cap *= 2; chain = (Expr**)realloc(chain, sizeof(Expr*) * cap); }
        chain[n++] = r;
    }

    Expr* list = expr_new_function(expr_new_symbol("List"), chain, n);
    free(chain);
    return list;
}

Expr* poly_gcd_internal(Expr* A, Expr* B, Expr** vars, size_t var_count);

/* The polynomial content -- the GCD of A's coefficients with respect  */
/* to the highest-indexed variable in `vars`. Bulk-extracts every       */
/* coefficient in a single pass over the expanded terms (much faster    */
/* than asking get_coeff for each degree, which would re-walk the       */
/* whole polynomial each time).                                         */
Expr* poly_content(Expr* A, Expr** vars, size_t var_count) {
    if (var_count == 0) return expr_copy(A);
    Expr* x = vars[var_count - 1];
    Expr* expandedA = expr_expand(A);
    int degA = get_degree_poly(expandedA, x);

    Expr** coeffs = NULL;
    bool bulk = (degA >= 0) && get_all_coeffs_expanded(expandedA, x, degA, &coeffs);

    Expr* g = expr_new_integer(0);
    for (int i = 0; i <= degA; i++) {
        Expr* c = bulk ? coeffs[i] : get_coeff_expanded(expandedA, x, i);
        if (!is_zero_poly(c)) {
            if (is_zero_poly(g)) {
                expr_free(g);
                g = c;
            } else {
                Expr* new_g = poly_gcd_internal(g, c, vars, var_count - 1);
                expr_free(g);
                expr_free(c);
                g = new_g;
            }
        } else {
            expr_free(c);
        }
    }
    free(coeffs);
    expr_free(expandedA);
    if (is_zero_poly(g)) {
        expr_free(g); return expr_new_integer(0);
    }
    return g;
}

/* Recursive multivariate gcd(A, B). `vars[k-1]` is the main variable; */
/* coefficient gcds in K[vars[0..k-2]] are obtained by recursive       */
/* descent. Outline:                                                   */
/*    A = cont(A) * pp(A)                                              */
/*    B = cont(B) * pp(B)                                              */
/*    g_cont = gcd(cont(A), cont(B))           (one variable smaller)  */
/*    g_pp   = result of pseudo-remainder reduction on pp(A), pp(B)    */
/*    return g_cont * g_pp                                             */
/* The leading coefficient is then sign-normalised so the result has   */
/* a positive lead.                                                    */
Expr* poly_gcd_internal(Expr* A, Expr* B, Expr** vars, size_t var_count) {
    if (is_zero_poly(A)) return expr_expand(B);
    if (is_zero_poly(B)) return expr_expand(A);

    if (var_count == 0) {
        return my_number_gcd(A, B);
    }

    Expr* x = vars[var_count - 1];

    Expr* contA = poly_content(A, vars, var_count);
    Expr* ppA = exact_poly_div(A, contA, vars, var_count);
    if (!ppA) ppA = expr_expand(A);

    Expr* contB = poly_content(B, vars, var_count);
    Expr* ppB = exact_poly_div(B, contB, vars, var_count);
    if (!ppB) ppB = expr_expand(B);

    Expr* contGCD = poly_gcd_internal(contA, contB, vars, var_count - 1);
    expr_free(contA); expr_free(contB);

    Expr* U = ppA;
    Expr* V = ppB;

    if (get_degree_poly(U, x) < get_degree_poly(V, x)) {
        Expr* tmp = U; U = V; V = tmp;
    }

    /* Size budget: pseudo-remainder coefficients in the multivariate
     * Euclidean step can grow exponentially when the coefficient ring
     * has many algebraic generators (e.g. Q[a, Sqrt[3], Sqrt[5], ...]
     * coming from a multi-radical denominator). Bail out when an
     * intermediate U or V exceeds the budget — returning the trivial
     * GCD `contGCD` (the content GCD of the inputs) is sound: the GCD
     * is always a divisor of the true GCD, and missing a cancellation
     * is a correctness-preserving pessimism. The subresultant PRS path
     * uses the same shape — see `size_budget` at poly.c:3340. */
    int64_t input_size = subres_leaf_count(A) + subres_leaf_count(B);
    /* Size budget: pseudo-remainder coefficients in the multivariate
     * Euclidean step can grow exponentially when the coefficient ring
     * carries algebraic generators (Sqrt[2], Sqrt[3], ... pseudo-vars
     * coming from multi-radical denominators). Bail out when an
     * intermediate U or V exceeds the budget — returning the trivial
     * GCD `contGCD` is sound: a smaller-than-true GCD is always a
     * valid divisor, and missing a cancellation is correctness-
     * preserving pessimism. The Plus combine in cancel_recursive
     * detects an inexact quotient via cancel_exact_div_strict (added
     * for the same case-13 hang) and falls back cleanly when the GCD
     * was only the content rather than the full polynomial GCD.
     *
     * Budget shape: max(input_size, 2000). The 2000 floor keeps the
     * common univariate cases (degree-100 polynomials with hundreds
     * of leaves) well within budget. The 1× input_size cap on larger
     * inputs catches multi-radical coefficient explosions early — the
     * subres-PRS path has the same shape for the same reason
     * (poly.c:3340). */
    int64_t size_budget = input_size;
    if (size_budget < 2000) size_budget = 2000;
    bool budget_blown = false;
    int iter_count = 0;

    while (!is_zero_poly(V)) {
        int64_t lU = subres_leaf_count(U);
        int64_t lV = subres_leaf_count(V);
        if (lU > size_budget || lV > size_budget) {
            budget_blown = true;
            break;
        }
        /* Hard iteration cap: even when sizes look bounded, a degenerate
         * pseudo-remainder sequence can spin for many iterations on
         * multivariate inputs whose `is_zero_poly` check is itself O(n²)
         * in coefficient-list extraction. 50 iterations is well above
         * any well-formed univariate / bivariate Euclidean run. */
        if (iter_count > 50) {
            budget_blown = true;
            break;
        }
        iter_count++;
        Expr* R = pseudo_rem(U, V, x);
        expr_free(U);
        U = V;
        if (is_zero_poly(R)) {
            expr_free(R);
            V = NULL;
            break;
        }
        Expr* contR = poly_content(R, vars, var_count);
        V = exact_poly_div(R, contR, vars, var_count);
        if (!V) V = expr_copy(R);
        expr_free(R); expr_free(contR);
    }
    if (V) expr_free(V);

    /* On budget exhaustion, return the content GCD alone — a valid
     * (over-estimate-safe) divisor that lets the cancel pipeline
     * proceed without coefficient explosion. */
    if (budget_blown) {
        expr_free(U);
        return expr_expand(contGCD);
    }
    
    // Normalize U to have positive leading coefficient
    Expr* lc = get_coeff(U, x, get_degree_poly(U, x));
    if (is_negative(lc)) {
        Expr* negU = internal_times((Expr*[]){expr_new_integer(-1), expr_copy(U)}, 2);
        expr_free(U);
        U = expr_expand(negU);
        expr_free(negU);
    }
    expr_free(lc);
    
    Expr* res = internal_times((Expr*[]){expr_copy(contGCD), expr_copy(U)}, 2);
    expr_free(contGCD); expr_free(U);
    Expr* expanded_res = expr_expand(res);
    expr_free(res);
    return expanded_res;
}

/* Extension-aware polynomial GCD: lift each input to QAUPoly[x] over Q(α)
 * and fold via qaupoly_gcd.  Inputs must be univariate in the same
 * variable (after the alpha-render symbols are excluded).  Returns NULL
 * on any structural failure — the caller is expected to fall back to the
 * standard (over-Q) GCD path or to leave the call unevaluated.
 *
 * Multi-arg form folds left-to-right: gcd(a, b, c) = gcd(gcd(a, b), c).
 */
static Expr* polynomialgcd_with_extension(Expr** argv, size_t argc,
                                          const Expr* alpha) {
    if (argc < 1 || !alpha) return NULL;

    Expr* alpha_render = NULL;
    QAExt* ext = qa_resolve_extension(alpha, &alpha_render);
    if (!ext) return NULL;

    /* Determine the polynomial variable.  Collect free symbols across
     * every input, drop the alpha-render symbol(s), require exactly one
     * remaining variable. */
    size_t vc = 0, vcap = 8;
    Expr** vars = malloc(sizeof(Expr*) * vcap);
    for (size_t i = 0; i < argc; i++) {
        collect_variables(argv[i], &vars, &vc, &vcap);
    }

    Expr* poly_var = NULL;
    size_t live = 0;
    for (size_t i = 0; i < vc; i++) {
        if (alpha_render && expr_eq(vars[i], alpha_render)) continue;
        poly_var = vars[i];
        live++;
    }

    if (live > 1) {
        /* Multivariate: not supported on the extension path. */
        for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
        free(vars);
        if (alpha_render) expr_free(alpha_render);
        qaext_free(ext);
        return NULL;
    }

    /* live == 0 means every input is a constant (in Q(α)).  Treat as a
     * "GCD of constants"; we delegate to the standard path. */
    if (live == 0) {
        for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
        free(vars);
        if (alpha_render) expr_free(alpha_render);
        qaext_free(ext);
        return NULL;
    }

    /* Lift each input. */
    QAUPoly** ps = malloc(sizeof(QAUPoly*) * argc);
    for (size_t i = 0; i < argc; i++) ps[i] = NULL;
    bool ok = true;
    for (size_t i = 0; i < argc; i++) {
        ps[i] = qa_expr_to_qaupoly(argv[i], poly_var, alpha_render, ext);
        if (!ps[i]) { ok = false; break; }
    }

    Expr* result = NULL;
    if (ok) {
        QAUPoly* g = qaupoly_copy(ps[0]);
        for (size_t i = 1; i < argc && g; i++) {
            QAUPoly* next = qaupoly_gcd(g, ps[i]);
            qaupoly_free(g);
            g = next;
        }
        if (g && !qaupoly_is_zero(g)) {
            result = qaupoly_to_expr_alpha(g, poly_var->data.symbol,
                                           alpha_render);
        }
        if (g) qaupoly_free(g);
    }

    for (size_t i = 0; i < argc; i++) {
        if (ps[i]) qaupoly_free(ps[i]);
    }
    free(ps);
    for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
    free(vars);
    if (alpha_render) expr_free(alpha_render);
    qaext_free(ext);

    return result;
}

/* Build a fresh function-call Expr with `head_sym` as the head and the
 * given args (deep-copied).  Used by the extension-stripping wrappers
 * that need to dispatch to the standard (no-option) builtin path with a
 * shorter argument list. */
static Expr* expr_rebuild_call(const Expr* original, Expr** args,
                               size_t argc) {
    Expr** new_args = malloc(sizeof(Expr*) * argc);
    for (size_t i = 0; i < argc; i++) new_args[i] = expr_copy(args[i]);
    Expr* call = expr_new_function(
        expr_copy(original->data.function.head), new_args, argc);
    free(new_args);
    return call;
}

/* PolynomialGCD[p1, p2, ...]                                          */
/*                                                                     */
/* Pre-processes the inputs by extracting:                              */
/*   - the integer content (numeric gcd of all literal coefficients);   */
/*   - any factors that appear with positive integer exponent in every  */
/*     argument (handles `Sqrt[x]`-style irrational generators);        */
/* Then defers to `poly_gcd_internal` for the symbolic remainder.       */
/* The final result is `numG * common_factors * symbolic_gcd`.          */
/*                                                                     */
/* Option `Extension -> α`: factor over Q(α) using the QAUPoly         */
/* machinery (see polynomialgcd_with_extension above).  `Extension ->   */
/* None` is accepted and treated as the default (no extension).         */
/* `Extension -> Automatic` triggers `extension_autodetect_args` and    */
/* routes through the single-generator α-path when one is found.        */
Expr* builtin_polynomialgcd(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;

    /* Strip a trailing Extension -> α option, if any.  Extension ->
     * Automatic triggers extension_autodetect_args and routes through
     * the single-generator α-path when a single generator is found. */
    size_t poly_argc = res->data.function.arg_count;
    bool auto_flag = false;
    const Expr* alpha = extract_extension_option_full(res, &poly_argc, &auto_flag);
    Expr* alpha_auto = NULL;
    QATower* auto_tower = NULL;
    if (!alpha && auto_flag && poly_argc >= 1) {
        auto_tower = extension_autodetect_args(res->data.function.args,
                                               poly_argc);
        if (auto_tower && auto_tower->n == 1) {
            alpha_auto = expr_copy(auto_tower->alpha_renders[0]);
            alpha = alpha_auto;
        }
    }
    if (poly_argc != res->data.function.arg_count) {
        if (poly_argc == 0) {
            if (alpha_auto) expr_free(alpha_auto);
            if (auto_tower) qa_tower_free(auto_tower);
            return NULL;
        }
        if (alpha) {
            Expr* ext_result = polynomialgcd_with_extension(
                res->data.function.args, poly_argc, alpha);
            if (alpha_auto) { expr_free(alpha_auto); alpha_auto = NULL; }
            if (auto_tower) { qa_tower_free(auto_tower); auto_tower = NULL; }
            if (ext_result) return ext_result;
            /* Fall through to non-extension path on failure (multivariate,
             * unrecognised α, lift failure, etc.). */
        } else if (auto_tower && auto_tower->n >= 2
                   && !internal_args_contain_inexact(res)) {
            /* Phase C: multi-α PolynomialGCD over Q(γ).  The single-α
             * branch above already consumed n==1 towers; here we use
             * the compositum to handle multi-generator inputs (e.g.
             * `PolynomialGCD[..., ..., Extension -> Automatic]` where
             * the inputs share both Sqrt[2] and Sqrt[3]). */
            Expr* tower_result = qa_polynomialgcd_with_tower(
                res->data.function.args, poly_argc, auto_tower);
            qa_tower_free(auto_tower); auto_tower = NULL;
            if (tower_result) return tower_result;
        }
        /* Cleanup any residual auto-detect state before recursion. */
        if (alpha_auto) { expr_free(alpha_auto); alpha_auto = NULL; }
        if (auto_tower) { qa_tower_free(auto_tower); auto_tower = NULL; }
        /* Rebuild the call with options removed and recurse so the rest
         * of this function works against the trimmed argument list. */
        Expr* trimmed = expr_rebuild_call(res, res->data.function.args,
                                          poly_argc);
        Expr* result = builtin_polynomialgcd(trimmed);
        expr_free(trimmed);
        return result;
    }
    /* No options were stripped, but auto-detect may still have produced
     * an α to try.  When the standard path is about to run, give the
     * extension path a chance first. */
    if (alpha) {
        Expr* ext_result = polynomialgcd_with_extension(
            res->data.function.args, poly_argc, alpha);
        if (alpha_auto) { expr_free(alpha_auto); alpha_auto = NULL; }
        if (auto_tower) { qa_tower_free(auto_tower); auto_tower = NULL; }
        if (ext_result) return ext_result;
    } else if (auto_tower && auto_tower->n >= 2
               && !internal_args_contain_inexact(res)) {
        /* Phase C: multi-α tower fast path before the no-extension
         * decomposition runs (see options-stripped branch above). */
        Expr* tower_result = qa_polynomialgcd_with_tower(
            res->data.function.args, poly_argc, auto_tower);
        qa_tower_free(auto_tower); auto_tower = NULL;
        if (tower_result) return tower_result;
    } else {
        if (alpha_auto) expr_free(alpha_auto);
        if (auto_tower) qa_tower_free(auto_tower);
    }

    /* Inexact coefficients (e.g. PolynomialGCD[x^2 - 1.0, x - 1.0]) cannot
     * be reasoned about by the rational-arithmetic GCD machinery below.
     * Force-rationalise the inputs, run the exact algorithm, and
     * numericalise the result so the caller observes inexact-in /
     * inexact-out semantics. */
    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_polynomialgcd);
    }

    size_t count = res->data.function.arg_count;
    if (count == 0) return NULL;
    BPList* bps = malloc(sizeof(BPList) * count);
    for (size_t i = 0; i < count; i++) {
        bp_init(&bps[i]);
        decompose_to_bp(res->data.function.args[i], &bps[i]);
    }

    Expr* numG = expr_new_integer(0);
    for (size_t i = 0; i < count; i++) {
        Expr* num_i = expr_new_integer(1);
        for (size_t j = 0; j < bps[i].count; j++) {
            if (bps[i].data[j].exp->type == EXPR_INTEGER && bps[i].data[j].exp->data.integer == 1) {
                if (is_number(bps[i].data[j].base)) {
                    Expr* next = internal_times((Expr*[]){num_i, expr_copy(bps[i].data[j].base)}, 2);
                    num_i = next;
                    expr_free(bps[i].data[j].exp);
                    bps[i].data[j].exp = expr_new_integer(0);
                } else if (bps[i].data[j].base->type == EXPR_FUNCTION
                    && bps[i].data[j].base->data.function.head->type == EXPR_SYMBOL
                    && bps[i].data[j].base->data.function.head->data.symbol == SYM_Plus) {
                    /* Plus base has no Times-factor numeric coefficient
                     * to peel off, but it can carry an integer content
                     * (the GCD of its summand coefficients).  Extract
                     * that content into num_i and replace the base with
                     * its primitive part so the downstream
                     * poly_gcd_internal does not also extract the same
                     * content via its g_cont path -- otherwise the
                     * combined result would multiply numG by the same
                     * factor twice.  Without this branch,
                     *   PolynomialGCD[4, 4 + 4 x]   returned 1
                     * (and Cancel[4/(4 + 4 x)] failed to reduce). */
                    int64_t c = get_int_content(bps[i].data[j].base);
                    if (c > 1) {
                        Expr* base = bps[i].data[j].base;
                        size_t n = base->data.function.arg_count;
                        Expr** new_args = malloc(sizeof(Expr*) * n);
                        for (size_t k = 0; k < n; k++) {
                            new_args[k] = eval_and_free(internal_divide(
                                (Expr*[]){expr_copy(base->data.function.args[k]),
                                          expr_new_integer(c)}, 2));
                        }
                        Expr* prim = eval_and_free(expr_new_function(
                            expr_copy(base->data.function.head),
                            new_args, n));
                        free(new_args);
                        Expr* next = internal_times(
                            (Expr*[]){num_i, expr_new_integer(c)}, 2);
                        num_i = next;
                        expr_free(bps[i].data[j].base);
                        bps[i].data[j].base = prim;
                    }
                }
            }
        }
        if (i == 0) {
            expr_free(numG);
            numG = num_i;
        } else {
            Expr* next = my_number_gcd(numG, num_i);
            expr_free(numG);
            expr_free(num_i);
            numG = next;
        }
    }
    
    Expr** common_args = malloc(sizeof(Expr*) * bps[0].count);
    size_t common_count = 0;
    
    for (size_t i = 0; i < bps[0].count; i++) {
        Expr* base = bps[0].data[i].base;
        if (bps[0].data[i].exp->type != EXPR_INTEGER || bps[0].data[i].exp->data.integer == 0) continue;
        
        int64_t min_exp = bps[0].data[i].exp->data.integer;
        bool in_all = true;
        for (size_t k = 1; k < count; k++) {
            bool found = false;
            for (size_t j = 0; j < bps[k].count; j++) {
                if (expr_eq(bps[k].data[j].base, base) && bps[k].data[j].exp->type == EXPR_INTEGER) {
                    found = true;
                    if (bps[k].data[j].exp->data.integer < min_exp) {
                        min_exp = bps[k].data[j].exp->data.integer;
                    }
                    break;
                }
            }
            if (!found) { in_all = false; break; }
        }
        
        if (in_all && min_exp > 0) {
            if (min_exp == 1) {
                common_args[common_count++] = expr_copy(base);
            } else {
                common_args[common_count++] = internal_power((Expr*[]){expr_copy(base), expr_new_integer(min_exp)}, 2);
            }
            for (size_t k = 0; k < count; k++) {
                for (size_t j = 0; j < bps[k].count; j++) {
                    if (expr_eq(bps[k].data[j].base, base)) {
                        int64_t v = bps[k].data[j].exp->data.integer - min_exp;
                        expr_free(bps[k].data[j].exp);
                        bps[k].data[j].exp = expr_new_integer(v);
                        break;
                    }
                }
            }
        }
    }

    Expr** rems = malloc(sizeof(Expr*) * count);
    for (size_t i = 0; i < count; i++) {
        rems[i] = rebuild_from_bp(&bps[i]);
        bp_free(&bps[i]);
    }
    free(bps);

    size_t v_count = 0, v_cap = 16;
    Expr** vars = malloc(sizeof(Expr*) * v_cap);
    for (size_t i = 0; i < count; i++) {
        collect_variables(rems[i], &vars, &v_count, &v_cap);
    }
    if (v_count > 0) qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);


    // Check if arguments are polynomials
    bool all_poly = true;
    for (size_t i = 0; i < count; i++) {
        if (!is_polynomial(rems[i], vars, v_count)) {
            all_poly = false;
            break;
        }
    }
    if (!all_poly) {
        for (size_t i = 0; i < count; i++) expr_free(rems[i]);
        free(rems);
        for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
        free(vars);
        expr_free(numG);
        for (size_t i = 0; i < common_count; i++) expr_free(common_args[i]);
        free(common_args);
        return NULL;
    }
    
    Expr* cur_gcd = expr_copy(rems[0]);
    for (size_t i = 1; i < count; i++) {
        Expr* next_gcd = poly_gcd_internal(cur_gcd, rems[i], vars, v_count);
        expr_free(cur_gcd);
        cur_gcd = next_gcd;
    }
    
    size_t final_count = 0;
    if (!(numG->type == EXPR_INTEGER && numG->data.integer == 1)) final_count++;
    final_count += common_count;
    if (!(cur_gcd->type == EXPR_INTEGER && cur_gcd->data.integer == 1)) final_count++;
    
    if (final_count == 0) {
        free(common_args);
        for(size_t i=0; i<v_count; i++) expr_free(vars[i]);
        free(vars);
        expr_free(numG);
        { for(size_t i=0; i<count; i++) expr_free(rems[i]); free(rems); }
        /* cur_gcd is the Integer 1 (final_count == 0 implies it). Reuse
         * it instead of allocating a fresh one — saves a malloc/free pair
         * and, more importantly, drops the reference we would otherwise
         * leak. */
        return cur_gcd;
    }
    
    Expr** final_args = malloc(sizeof(Expr*) * final_count);
    size_t idx = 0;
    if (!(numG->type == EXPR_INTEGER && numG->data.integer == 1)) final_args[idx++] = numG;
    else expr_free(numG);
    for (size_t i = 0; i < common_count; i++) final_args[idx++] = common_args[i];
    if (!(cur_gcd->type == EXPR_INTEGER && cur_gcd->data.integer == 1)) final_args[idx++] = cur_gcd;
    else expr_free(cur_gcd);
    
    Expr* result;
    if (idx == 1) result = final_args[0];
    else result = internal_times(final_args, idx);
    
    free(common_args); 
    for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
    free(vars);
    { for(size_t i=0; i<count; i++) expr_free(rems[i]); free(rems); }
    free(final_args);
    return result;
}

static Expr* my_number_lcm(Expr* a, Expr* b) {
    int64_t ca = get_int_content(a);
    int64_t cb = get_int_content(b);
    int64_t l = lcm(ca, cb);
    if (l == 0) l = 1;
    return expr_new_integer(l);
}

/* PolynomialLCM[p1, p2, ...]                                          */
/*                                                                     */
/* Same flavour as PolynomialGCD: split off integer LCM and the        */
/* maximum-exponent symbolic factors that appear across all inputs,    */
/* then combine the polynomial parts via lcm(a,b) = a*b / gcd(a,b).    */
/* The final result also folds in the largest negative exponents of    */
/* any denominator generators, matching Mathematica's handling of      */
/* rational expressions.                                                */
/* Extension-aware polynomial LCM via QAUPoly: lcm(a, b) = a*b / gcd(a, b).
 * Multi-arg form folds left-to-right. */
static Expr* polynomiallcm_with_extension(Expr** argv, size_t argc,
                                          const Expr* alpha) {
    if (argc < 1 || !alpha) return NULL;

    Expr* alpha_render = NULL;
    QAExt* ext = qa_resolve_extension(alpha, &alpha_render);
    if (!ext) return NULL;

    size_t vc = 0, vcap = 8;
    Expr** vars = malloc(sizeof(Expr*) * vcap);
    for (size_t i = 0; i < argc; i++) {
        collect_variables(argv[i], &vars, &vc, &vcap);
    }
    Expr* poly_var = NULL;
    size_t live = 0;
    for (size_t i = 0; i < vc; i++) {
        if (alpha_render && expr_eq(vars[i], alpha_render)) continue;
        poly_var = vars[i];
        live++;
    }
    if (live != 1) {
        for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
        free(vars);
        if (alpha_render) expr_free(alpha_render);
        qaext_free(ext);
        return NULL;
    }

    QAUPoly** ps = malloc(sizeof(QAUPoly*) * argc);
    for (size_t i = 0; i < argc; i++) ps[i] = NULL;
    bool ok = true;
    for (size_t i = 0; i < argc; i++) {
        ps[i] = qa_expr_to_qaupoly(argv[i], poly_var, alpha_render, ext);
        if (!ps[i]) { ok = false; break; }
    }

    Expr* result = NULL;
    if (ok) {
        QAUPoly* L = qaupoly_copy(ps[0]);
        for (size_t i = 1; i < argc && L; i++) {
            QAUPoly* g = qaupoly_gcd(L, ps[i]);
            QAUPoly* prod = qaupoly_mul(L, ps[i]);
            QAUPoly *q = NULL, *r = NULL;
            if (g && qaupoly_divrem(prod, g, &q, &r)) {
                qaupoly_free(L);
                L = q;
                qaupoly_free(r);
            } else {
                /* Should not happen — gcd divides product exactly. */
                if (q) qaupoly_free(q);
                if (r) qaupoly_free(r);
                qaupoly_free(L);
                L = NULL;
            }
            if (g) qaupoly_free(g);
            qaupoly_free(prod);
        }
        if (L && !qaupoly_is_zero(L)) {
            QAUPoly* monic = qaupoly_make_monic(L);
            if (monic) {
                result = qaupoly_to_expr_alpha(monic, poly_var->data.symbol,
                                               alpha_render);
                qaupoly_free(monic);
            }
        }
        if (L) qaupoly_free(L);
    }

    for (size_t i = 0; i < argc; i++) if (ps[i]) qaupoly_free(ps[i]);
    free(ps);
    for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
    free(vars);
    if (alpha_render) expr_free(alpha_render);
    qaext_free(ext);
    return result;
}

Expr* builtin_polynomiallcm(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;

    /* Strip a trailing Extension -> α option.  Extension -> Automatic
     * triggers extension_autodetect_args on the polynomial arguments. */
    size_t poly_argc = res->data.function.arg_count;
    bool auto_flag = false;
    const Expr* alpha = extract_extension_option_full(res, &poly_argc, &auto_flag);
    Expr* alpha_auto = NULL;
    QATower* auto_tower = NULL;
    if (!alpha && auto_flag && poly_argc >= 1) {
        auto_tower = extension_autodetect_args(res->data.function.args,
                                               poly_argc);
        if (auto_tower && auto_tower->n == 1) {
            alpha_auto = expr_copy(auto_tower->alpha_renders[0]);
            alpha = alpha_auto;
        }
    }
    if (poly_argc != res->data.function.arg_count) {
        if (poly_argc == 0) {
            if (alpha_auto) expr_free(alpha_auto);
            if (auto_tower) qa_tower_free(auto_tower);
            return NULL;
        }
        if (alpha) {
            Expr* ext_result = polynomiallcm_with_extension(
                res->data.function.args, poly_argc, alpha);
            if (alpha_auto) { expr_free(alpha_auto); alpha_auto = NULL; }
            if (auto_tower) { qa_tower_free(auto_tower); auto_tower = NULL; }
            if (ext_result) return ext_result;
        } else if (auto_tower && auto_tower->n >= 2
                   && !internal_args_contain_inexact(res)) {
            /* Phase C: multi-α PolynomialLCM over Q(γ).  Mirrors the
             * GCD wiring above. */
            Expr* tower_result = qa_polynomiallcm_with_tower(
                res->data.function.args, poly_argc, auto_tower);
            qa_tower_free(auto_tower); auto_tower = NULL;
            if (tower_result) return tower_result;
        }
        if (alpha_auto) { expr_free(alpha_auto); alpha_auto = NULL; }
        if (auto_tower) { qa_tower_free(auto_tower); auto_tower = NULL; }
        Expr* trimmed = expr_rebuild_call(res, res->data.function.args,
                                          poly_argc);
        Expr* result = builtin_polynomiallcm(trimmed);
        expr_free(trimmed);
        return result;
    }
    /* No options stripped but auto-detect may still have produced an α. */
    if (alpha) {
        Expr* ext_result = polynomiallcm_with_extension(
            res->data.function.args, poly_argc, alpha);
        if (alpha_auto) { expr_free(alpha_auto); alpha_auto = NULL; }
        if (auto_tower) { qa_tower_free(auto_tower); auto_tower = NULL; }
        if (ext_result) return ext_result;
    } else if (auto_tower && auto_tower->n >= 2
               && !internal_args_contain_inexact(res)) {
        /* Phase C: multi-α tower fast path. */
        Expr* tower_result = qa_polynomiallcm_with_tower(
            res->data.function.args, poly_argc, auto_tower);
        qa_tower_free(auto_tower); auto_tower = NULL;
        if (tower_result) return tower_result;
    } else {
        if (alpha_auto) expr_free(alpha_auto);
        if (auto_tower) qa_tower_free(auto_tower);
    }

    /* Force-rationalise inexact coefficients so the exact LCM algorithm
     * applies; numericalise the final result. See builtin_polynomialgcd
     * above for the rationale. */
    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_polynomiallcm);
    }

    size_t count = res->data.function.arg_count;
    if (count == 0) return NULL;
    BPList* bps = malloc(sizeof(BPList) * count);
    for (size_t i = 0; i < count; i++) {
        bp_init(&bps[i]);
        decompose_to_bp(res->data.function.args[i], &bps[i]);
    }

    Expr* numL = expr_new_integer(1);
    for (size_t i = 0; i < count; i++) {
        Expr* num_i = expr_new_integer(1);
        for (size_t j = 0; j < bps[i].count; j++) {
            if (bps[i].data[j].exp->type == EXPR_INTEGER && bps[i].data[j].exp->data.integer == 1) {
                if (is_number(bps[i].data[j].base)) {
                    Expr* next = internal_times((Expr*[]){num_i, expr_copy(bps[i].data[j].base)}, 2);
                    num_i = next;
                    expr_free(bps[i].data[j].exp);
                    bps[i].data[j].exp = expr_new_integer(0);
                }
            }
        }
        if (i == 0) {
            expr_free(numL);
            numL = num_i;
        } else {
            Expr* next = my_number_lcm(numL, num_i);
            expr_free(numL);
            expr_free(num_i);
            numL = next;
        }
    }
    
    Expr** common_args = malloc(sizeof(Expr*) * 1024);
    size_t common_count = 0;
    size_t common_cap = 1024;
    
    for (size_t i = 0; i < bps[0].count; i++) {
        Expr* base = bps[0].data[i].base;
        if (bps[0].data[i].exp->type != EXPR_INTEGER || bps[0].data[i].exp->data.integer <= 0) continue;
        
        int64_t min_exp = bps[0].data[i].exp->data.integer;
        bool in_all = true;
        for (size_t k = 1; k < count; k++) {
            bool found = false;
            for (size_t j = 0; j < bps[k].count; j++) {
                if (expr_eq(bps[k].data[j].base, base) && bps[k].data[j].exp->type == EXPR_INTEGER) {
                    found = true;
                    if (bps[k].data[j].exp->data.integer < min_exp) {
                        min_exp = bps[k].data[j].exp->data.integer;
                    }
                    break;
                }
            }
            if (!found) { in_all = false; break; }
        }
        
        if (in_all && min_exp > 0) {
            if (common_count == common_cap) { common_cap *= 2; common_args = realloc(common_args, sizeof(Expr*) * common_cap); }
            if (min_exp == 1) {
                common_args[common_count++] = expr_copy(base);
            } else {
                common_args[common_count++] = internal_power((Expr*[]){expr_copy(base), expr_new_integer(min_exp)}, 2);
            }
            for (size_t k = 0; k < count; k++) {
                for (size_t j = 0; j < bps[k].count; j++) {
                    if (expr_eq(bps[k].data[j].base, base)) {
                        int64_t v = bps[k].data[j].exp->data.integer - min_exp;
                        expr_free(bps[k].data[j].exp);
                        bps[k].data[j].exp = expr_new_integer(v);
                        break;
                    }
                }
            }
        }
    }

    // Handle denominators (negative exponents)
    Expr** den_bases = malloc(sizeof(Expr*) * 1024);
    size_t den_bases_count = 0;
    size_t den_bases_cap = 1024;
    
    for (size_t i = 0; i < count; i++) {
        for (size_t j = 0; j < bps[i].count; j++) {
            if (bps[i].data[j].exp->type != EXPR_INTEGER || bps[i].data[j].exp->data.integer >= 0) continue;
            Expr* base = bps[i].data[j].base;
            bool found = false;
            for (size_t k = 0; k < den_bases_count; k++) {
                if (expr_eq(den_bases[k], base)) { found = true; break; }
            }
            if (!found) {
                if (den_bases_count == den_bases_cap) { den_bases_cap *= 2; den_bases = realloc(den_bases, sizeof(Expr*) * den_bases_cap); }
                den_bases[den_bases_count++] = expr_copy(base);
            }
        }
    }
    
    for (size_t i = 0; i < den_bases_count; i++) {
        Expr* base = den_bases[i];
        int64_t max_exp = -INT64_MAX; 
        for (size_t k = 0; k < count; k++) {
            int64_t e = 0;
            for (size_t j = 0; j < bps[k].count; j++) {
                if (expr_eq(bps[k].data[j].base, base) && bps[k].data[j].exp->type == EXPR_INTEGER) {
                    e = bps[k].data[j].exp->data.integer;
                    break;
                }
            }
            if (e > max_exp) max_exp = e;
        }
        
        if (max_exp < 0) {
            if (common_count == common_cap) { common_cap *= 2; common_args = realloc(common_args, sizeof(Expr*) * common_cap); }
            if (max_exp == -1) {
                common_args[common_count++] = internal_power((Expr*[]){expr_copy(base), expr_new_integer(-1)}, 2);
            } else {
                common_args[common_count++] = internal_power((Expr*[]){expr_copy(base), expr_new_integer(max_exp)}, 2);
            }
        }
        
        for (size_t k = 0; k < count; k++) {
            for (size_t j = 0; j < bps[k].count; j++) {
                if (expr_eq(bps[k].data[j].base, base)) {
                    expr_free(bps[k].data[j].exp);
                    bps[k].data[j].exp = expr_new_integer(0);
                    break;
                }
            }
        }
    }

    for (size_t i = 0; i < den_bases_count; i++) expr_free(den_bases[i]);
    free(den_bases);

    Expr** rems = malloc(sizeof(Expr*) * count);
    for (size_t i = 0; i < count; i++) {
        rems[i] = rebuild_from_bp(&bps[i]);
        bp_free(&bps[i]);
    }
    free(bps);
    
    size_t v_count = 0, v_cap = 16;
    Expr** vars = malloc(sizeof(Expr*) * v_cap);
    for (size_t i = 0; i < count; i++) {
        collect_variables(rems[i], &vars, &v_count, &v_cap);
    }
    if (v_count > 0) qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);

    /* If any input is not a polynomial in the collected variables (e.g.   */
    /* contains rational/symbolic exponents like y^(1/3)), bail out and    */
    /* leave the expression unevaluated -- otherwise poly_gcd_internal /   */
    /* exact_poly_div will loop forever on non-polynomial inputs.          */
    bool all_poly = true;
    for (size_t i = 0; i < count; i++) {
        if (!is_polynomial(rems[i], vars, v_count)) { all_poly = false; break; }
    }
    if (!all_poly) {
        for (size_t i = 0; i < count; i++) expr_free(rems[i]);
        free(rems);
        for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
        free(vars);
        expr_free(numL);
        for (size_t i = 0; i < common_count; i++) expr_free(common_args[i]);
        free(common_args);
        return NULL;
    }

    Expr* cur_lcm = expr_copy(rems[0]);
    for (size_t i = 1; i < count; i++) {
        Expr* cur_gcd = poly_gcd_internal(cur_lcm, rems[i], vars, v_count);
        Expr* next_lcm;
        if (cur_gcd->type == EXPR_INTEGER && cur_gcd->data.integer == 1) {
            next_lcm = internal_times((Expr*[]){expr_copy(cur_lcm), expr_copy(rems[i])}, 2);
        } else {
            Expr* Q1 = exact_poly_div(cur_lcm, cur_gcd, vars, v_count);
            Expr* Q2 = exact_poly_div(rems[i], cur_gcd, vars, v_count);
            if (!Q1) Q1 = expr_copy(cur_lcm); 
            if (!Q2) Q2 = expr_copy(rems[i]); 
            
            int c1 = (Q1->type == EXPR_FUNCTION && Q1->data.function.head->type == EXPR_SYMBOL && Q1->data.function.head->data.symbol == SYM_Plus) ? Q1->data.function.arg_count : 1;
            int c2 = (Q2->type == EXPR_FUNCTION && Q2->data.function.head->type == EXPR_SYMBOL && Q2->data.function.head->data.symbol == SYM_Plus) ? Q2->data.function.arg_count : 1;
            
            if (c1 <= c2) {
                next_lcm = internal_times((Expr*[]){Q1, expr_copy(rems[i])}, 2);
                expr_free(Q2);
            } else {
                next_lcm = internal_times((Expr*[]){expr_copy(cur_lcm), Q2}, 2);
                expr_free(Q1);
            }
        }
        expr_free(cur_gcd);
        expr_free(cur_lcm); 
        cur_lcm = next_lcm;
    }
    
    size_t final_count = 0;
    if (!(numL->type == EXPR_INTEGER && numL->data.integer == 1)) final_count++;
    final_count += common_count;
    if (!(cur_lcm->type == EXPR_INTEGER && cur_lcm->data.integer == 1)) final_count++;
    
    if (final_count == 0) {
        free(common_args);
        for(size_t i=0; i<v_count; i++) expr_free(vars[i]);
        free(vars);
        expr_free(numL);
        { for(size_t i=0; i<count; i++) expr_free(rems[i]); free(rems); }
        /* cur_lcm is the Integer 1 (final_count == 0 implies it). Reuse
         * it instead of allocating a fresh one and leaking the existing. */
        return cur_lcm;
    }
    
    Expr** final_args = malloc(sizeof(Expr*) * final_count);
    size_t idx = 0;
    if (!(numL->type == EXPR_INTEGER && numL->data.integer == 1)) final_args[idx++] = numL;
    else expr_free(numL);
    for (size_t i = 0; i < common_count; i++) final_args[idx++] = common_args[i];
    if (!(cur_lcm->type == EXPR_INTEGER && cur_lcm->data.integer == 1)) final_args[idx++] = cur_lcm;
    else expr_free(cur_lcm);
    
    Expr* result;
    if (idx == 1) result = final_args[0];
    else result = internal_times(final_args, idx);
    
    free(common_args); 
    for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
    free(vars);
    { for(size_t i=0; i<count; i++) expr_free(rems[i]); free(rems); }
    free(final_args);
    return result;
}


static bool is_threadable_head(Expr* head) {
    if (head->type != EXPR_SYMBOL) return false;
    const char* s = head->data.symbol;
    return s == SYM_List || s == SYM_Equal || s == SYM_Unequal ||
           s == SYM_Less || s == SYM_LessEqual || s == SYM_Greater ||
           s == SYM_GreaterEqual || s == SYM_And || s == SYM_Or ||
           s == SYM_Not || s == SYM_SameQ || s == SYM_UnsameQ;
}

typedef struct {
    Expr* base;
    Expr* exp;
    Expr** coeffs;
    size_t coeff_count;
    size_t coeff_cap;
} CollectGroup;

/* Compute the residual BP list = `term_bp` minus k copies of every
 * factor in `var_bp`. Caller pre-validated that get_k(term_bp,
 * var_bp) == k, so each var_bp factor is present in term_bp with
 * sufficient exponent. Returns a fresh BPList with the leftover
 * factors. Uses integer-exp arithmetic; non-integer exponents are
 * preserved unchanged when their structural copy matches (those would
 * have produced k==1 from get_k anyway). */
static void bp_compute_residual(BPList* term_bp, BPList* var_bp,
                                int64_t k, BPList* out) {
    bp_init(out);
    /* Build a "consumed exponent per term_bp index" map. */
    int64_t* consumed = calloc(term_bp->count, sizeof(int64_t));
    for (size_t vi = 0; vi < var_bp->count; vi++) {
        Expr* vb = var_bp->data[vi].base;
        Expr* ve = var_bp->data[vi].exp;
        for (size_t tj = 0; tj < term_bp->count; tj++) {
            if (expr_eq(vb, term_bp->data[tj].base)) {
                if (ve->type == EXPR_INTEGER) {
                    consumed[tj] += k * ve->data.integer;
                } else {
                    /* Non-integer var exp: get_k returned 1, exponent
                     * matches structurally — consume the whole. */
                    consumed[tj] = -1; /* sentinel: drop entirely */
                }
                break;
            }
        }
    }
    for (size_t tj = 0; tj < term_bp->count; tj++) {
        Expr* tb = term_bp->data[tj].base;
        Expr* te = term_bp->data[tj].exp;
        if (consumed[tj] == -1) continue;  /* dropped */
        if (te->type == EXPR_INTEGER && consumed[tj] > 0) {
            int64_t rem = te->data.integer - consumed[tj];
            if (rem == 0) continue;
            Expr* re = expr_new_integer(rem);
            bp_add(out, tb, re);
            expr_free(re);
        } else {
            bp_add(out, tb, te);
        }
    }
    free(consumed);
}

static Expr* collect_internal(Expr* expr, Expr** vars, size_t num_vars, size_t var_idx, Expr* h) {
    if (expr->type == EXPR_FUNCTION && is_threadable_head(expr->data.function.head)) {
        size_t count = expr->data.function.arg_count;
        Expr** new_args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            new_args[i] = collect_internal(expr->data.function.args[i], vars, num_vars, var_idx, h);
        }
        return eval_and_free(expr_new_function(expr_copy(expr->data.function.head), new_args, count));
    }

    if (var_idx == num_vars) {
        if (h) {
            return eval_and_free(expr_new_function(expr_copy(h), (Expr*[]){expr_copy(expr)}, 1));
        }
        return expr_copy(expr);
    }

    Expr* var = vars[var_idx];

    /* Plus keys must be treated as opaque atoms: expanding the input
     * with respect to a Plus pattern would distribute it across Times
     * factors and destroy the Plus subterm we are trying to collect
     * on. Skip expansion in that case. Times / atomic / Power keys
     * still benefit from expansion to surface the var^k structure.
     *
     * MMA's `Collect[a (c+s) + b (c+s), c+s]` returns `(a+b)(c+s)`,
     * not the distributed `a c + a s + b c + b s`. */
    bool key_is_plus = (var->type == EXPR_FUNCTION
                        && var->data.function.head
                        && var->data.function.head->type == EXPR_SYMBOL
                        && var->data.function.head->data.symbol == SYM_Plus);

    Expr* expanded = key_is_plus ? expr_copy(expr) : expr_expand_patt(expr, var);
    if (!expanded) return expr_copy(expr);

    size_t term_count = 1;
    Expr** terms = &expanded;
    if (expanded->type == EXPR_FUNCTION && expanded->data.function.head->type == EXPR_SYMBOL && expanded->data.function.head->data.symbol == SYM_Plus) {
        term_count = expanded->data.function.arg_count;
        terms = expanded->data.function.args;
    }

    /* Decompose the key. Two paths:
     *
     * 1. Single-base key (var_bp.count == 1) -- atoms `x`, powers
     *    `Power[x, e]`, Plus keys `Plus[c, s]`. Group each term by
     *    `term_exp / var_exp`, the rational/symbolic exponent at which
     *    the key's base appears in the term. This handles non-integer
     *    exponents like `Sqrt[x] = x^(1/2)` correctly.
     *
     * 2. Multi-factor key (var_bp.count > 1) -- Times monomials
     *    `x*y*z`, etc. Use `get_k` to find the integer multiplicity at
     *    which `var` divides the term. For these keys, fractional
     *    "k" doesn't have a coherent meaning (you can't have half a
     *    `x*y`), so the integer-only ratio is what matches MMA's
     *    behaviour. */
    BPList var_bp;
    bp_init(&var_bp);
    decompose_to_bp(var, &var_bp);

    typedef struct { Expr* k_expr; int64_t k_int; bool int_form; Expr** coeffs; size_t cc, ccap; } CGrp;
    size_t gcap = 8, gcount = 0;
    CGrp* groups = malloc(sizeof(CGrp) * gcap);

    bool single_base = (var_bp.count == 1);
    Expr* sb_base = single_base ? var_bp.data[0].base : NULL;
    Expr* sb_exp  = single_base ? var_bp.data[0].exp  : NULL;

    for (size_t i = 0; i < term_count; i++) {
        Expr* term = terms[i];
        BPList term_bp;
        bp_init(&term_bp);
        decompose_to_bp(term, &term_bp);

        Expr* k_expr = NULL;       /* group key for single_base path */
        int64_t k_int = 0;          /* group key for multi-factor path */
        bool int_form = !single_base;
        Expr* coeff = NULL;

        if (single_base) {
            int matched_idx = -1;
            for (size_t j = 0; j < term_bp.count; j++) {
                if (expr_eq(sb_base, term_bp.data[j].base)) {
                    matched_idx = (int)j;
                    break;
                }
            }
            /* Power-of-power fallback: Collect target is `Power[B, e_t]`
             * (an atomic Power decomposition because e_t is symbolic /
             * Times) and a term factor is some other power of the same
             * primitive base B. Compute k such that the term factor
             * equals (target)^k by exponent ratio: e_term_eff / e_t.
             *
             * Examples handled here that the expr_eq match cannot:
             *   target x^c, term factor x^(2c)              -> k = 2
             *   target Power[x, c], term factor Power[x, 3c] -> k = 3
             *   target x^(c+1), term factor x^(2c+2)         -> k = 2
             *
             * Only matches positive integer ratios (matches the integer-
             * integer sb_exp branch's semantics). */
            int pop_match_k_int = 0;
            if (matched_idx == -1 &&
                sb_base->type == EXPR_FUNCTION &&
                sb_base->data.function.head->type == EXPR_SYMBOL &&
                sb_base->data.function.head->data.symbol == SYM_Power &&
                sb_base->data.function.arg_count == 2) {
                Expr* B   = sb_base->data.function.args[0];
                Expr* e_t = sb_base->data.function.args[1];
                for (size_t j = 0; j < term_bp.count; j++) {
                    Expr* tb = term_bp.data[j].base;
                    Expr* tx = term_bp.data[j].exp;
                    /* Effective B-exponent contributed by this entry. */
                    Expr* e_term_eff = NULL;
                    if (expr_eq(B, tb)) {
                        e_term_eff = expr_copy(tx);
                    } else if (tb->type == EXPR_FUNCTION &&
                               tb->data.function.head->type == EXPR_SYMBOL &&
                               tb->data.function.head->data.symbol == SYM_Power &&
                               tb->data.function.arg_count == 2 &&
                               expr_eq(B, tb->data.function.args[0])) {
                        e_term_eff = internal_times((Expr*[]){
                            expr_copy(tb->data.function.args[1]),
                            expr_copy(tx) }, 2);
                    }
                    if (!e_term_eff) continue;

                    Expr* inv_e_t = internal_power((Expr*[]){
                        expr_copy(e_t), expr_new_integer(-1) }, 2);
                    Expr* k_cand = internal_times((Expr*[]){
                        e_term_eff, inv_e_t }, 2);
                    if (k_cand->type == EXPR_INTEGER && k_cand->data.integer >= 1) {
                        matched_idx = (int)j;
                        pop_match_k_int = (int)k_cand->data.integer;
                    }
                    expr_free(k_cand);
                    if (matched_idx != -1) break;
                }
            }
            if (matched_idx == -1) {
                k_expr = expr_new_integer(0);
                coeff = rebuild_from_bp(&term_bp);
            } else if (pop_match_k_int > 0) {
                /* Power-of-power match: k_expr is the positive integer k.
                 * Coefficient drops the matched entry (the only B-bearing
                 * entry; we matched against a single term factor). */
                k_expr = expr_new_integer(pop_match_k_int);
                BPList rest;
                bp_init(&rest);
                for (size_t j = 0; j < term_bp.count; j++) {
                    if ((int)j == matched_idx) continue;
                    bp_add(&rest, term_bp.data[j].base, term_bp.data[j].exp);
                }
                coeff = rebuild_from_bp(&rest);
                bp_free(&rest);
            } else {
                Expr* term_exp = term_bp.data[matched_idx].exp;
                /* k = term_exp / var_exp. */
                if (sb_exp->type == EXPR_INTEGER && sb_exp->data.integer == 1) {
                    k_expr = expr_copy(term_exp);
                } else if (sb_exp->type == EXPR_INTEGER && term_exp->type == EXPR_INTEGER) {
                    int64_t a = term_exp->data.integer;
                    int64_t b = sb_exp->data.integer;
                    if (b != 0 && a % b == 0) {
                        k_expr = expr_new_integer(a / b);
                    } else {
                        /* Non-integer ratio with Power key -> no clean
                         * fit; treat as no-match. */
                        k_expr = expr_new_integer(0);
                    }
                } else {
                    /* Symbolic / rational exponents: build the literal
                     * ratio so equal terms collect (group key compared
                     * by expr_eq). */
                    Expr* inv_var = internal_power((Expr*[]){expr_copy(sb_exp), expr_new_integer(-1)}, 2);
                    k_expr = internal_times((Expr*[]){expr_copy(term_exp), inv_var}, 2);
                }
                /* Coefficient: term without the matched BP entry. */
                bool zero_k = (k_expr->type == EXPR_INTEGER && k_expr->data.integer == 0);
                if (zero_k) {
                    coeff = rebuild_from_bp(&term_bp);
                } else {
                    /* Drop the matched entry. */
                    BPList rest;
                    bp_init(&rest);
                    for (size_t j = 0; j < term_bp.count; j++) {
                        if ((int)j == matched_idx) continue;
                        bp_add(&rest, term_bp.data[j].base, term_bp.data[j].exp);
                    }
                    coeff = rebuild_from_bp(&rest);
                    bp_free(&rest);
                }
            }
        } else {
            k_int = get_k(&term_bp, &var_bp);
            if (k_int == 0) {
                coeff = rebuild_from_bp(&term_bp);
            } else {
                BPList res_bp;
                bp_compute_residual(&term_bp, &var_bp, k_int, &res_bp);
                coeff = rebuild_from_bp(&res_bp);
                bp_free(&res_bp);
            }
        }
        bp_free(&term_bp);

        /* Find or create group with matching k. */
        int found = -1;
        for (size_t j = 0; j < gcount; j++) {
            if (int_form != groups[j].int_form) continue;
            if (int_form) {
                if (groups[j].k_int == k_int) { found = (int)j; break; }
            } else {
                if (expr_eq(k_expr, groups[j].k_expr)) { found = (int)j; break; }
            }
        }
        if (found == -1) {
            if (gcount == gcap) { gcap *= 2; groups = realloc(groups, sizeof(CGrp) * gcap); }
            groups[gcount].int_form = int_form;
            groups[gcount].k_int = k_int;
            groups[gcount].k_expr = single_base ? k_expr : NULL;
            groups[gcount].cc = 0;
            groups[gcount].ccap = 4;
            groups[gcount].coeffs = malloc(sizeof(Expr*) * 4);
            found = (int)gcount++;
        } else {
            if (single_base && k_expr) expr_free(k_expr);
        }
        if (groups[found].cc == groups[found].ccap) {
            groups[found].ccap *= 2;
            groups[found].coeffs = realloc(groups[found].coeffs, sizeof(Expr*) * groups[found].ccap);
        }
        groups[found].coeffs[groups[found].cc++] = coeff;
    }

    Expr** final_terms = malloc(sizeof(Expr*) * gcount);
    size_t final_count = 0;
    for (size_t i = 0; i < gcount; i++) {
        Expr* coeff_sum;
        if (groups[i].cc == 1) {
            coeff_sum = expr_copy(groups[i].coeffs[0]);
        } else {
            Expr** ca = malloc(sizeof(Expr*) * groups[i].cc);
            for (size_t j = 0; j < groups[i].cc; j++) ca[j] = expr_copy(groups[i].coeffs[j]);
            coeff_sum = internal_plus(ca, groups[i].cc);
            free(ca);
        }

        Expr* collected_coeff = collect_internal(coeff_sum, vars, num_vars, var_idx + 1, h);
        expr_free(coeff_sum);

        Expr* term;
        bool zero_group;
        if (groups[i].int_form) {
            zero_group = (groups[i].k_int == 0);
        } else {
            zero_group = (groups[i].k_expr->type == EXPR_INTEGER && groups[i].k_expr->data.integer == 0);
        }
        if (zero_group) {
            term = collected_coeff;
        } else {
            Expr* var_pow;
            if (groups[i].int_form) {
                var_pow = (groups[i].k_int == 1)
                    ? expr_copy(var)
                    : internal_power((Expr*[]){expr_copy(var), expr_new_integer(groups[i].k_int)}, 2);
            } else {
                /* Single-base key: rebuild Power[base, k_expr] (or
                 * just `base` when k_expr = 1). */
                bool k_is_one = (groups[i].k_expr->type == EXPR_INTEGER
                                 && groups[i].k_expr->data.integer == 1);
                if (k_is_one) {
                    /* var = base^var_exp; with k=1 this is just var. */
                    var_pow = expr_copy(var);
                } else {
                    /* k_expr is the multiplicity of `var` in the term:
                     * term contains var^k_expr. Build Power[var, k_expr]
                     * (NOT Power[sb_base, k_expr]) so the exponent on
                     * the actual collect target is faithful. The
                     * evaluator then folds Power[Power[x, e], k] to
                     * Power[x, e*k] automatically. */
                    var_pow = internal_power((Expr*[]){expr_copy(var), expr_copy(groups[i].k_expr)}, 2);
                }
            }
            term = internal_times((Expr*[]){collected_coeff, var_pow}, 2);
        }
        final_terms[final_count++] = term;

        if (!groups[i].int_form && groups[i].k_expr) expr_free(groups[i].k_expr);
        for (size_t j = 0; j < groups[i].cc; j++) expr_free(groups[i].coeffs[j]);
        free(groups[i].coeffs);
    }
    free(groups);
    bp_free(&var_bp);

    Expr* result;
    if (final_count == 0) result = expr_new_integer(0);
    else if (final_count == 1) result = final_terms[0];
    else result = internal_plus(final_terms, final_count);

    free(final_terms);
    expr_free(expanded);
    return result;
}

Expr* builtin_collect(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2 || res->data.function.arg_count > 3) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* vars_expr = res->data.function.args[1];
    Expr* h = (res->data.function.arg_count == 3) ? res->data.function.args[2] : NULL;

    size_t num_vars = 1;
    Expr** vars = &vars_expr;

    if (vars_expr->type == EXPR_FUNCTION && vars_expr->data.function.head->type == EXPR_SYMBOL && vars_expr->data.function.head->data.symbol == SYM_List) {
        num_vars = vars_expr->data.function.arg_count;
        vars = vars_expr->data.function.args;
    }

    return collect_internal(expr, vars, num_vars, 0, h);
}

/* Recursive worker for CoefficientList[expr, {v1, v2, ...}]. Builds a  */
/* nested list whose shape mirrors the variable order. Uses the bulk    */
/* coefficient-extraction helper so each level walks `expanded` once    */
/* instead of (degree+1) times.                                         */
static Expr* coeff_list_rec(Expr* expr, Expr** vars, int* max_degrees, size_t num_vars, size_t var_idx) {
    if (var_idx == num_vars) return expr_copy(expr);
    Expr* var = vars[var_idx];
    int d = max_degrees[var_idx];
    Expr* expanded = expr_expand(expr);
    if (!expanded) return expr_copy(expr);

    Expr** coeffs = NULL;
    bool bulk = get_all_coeffs_expanded(expanded, var, d, &coeffs);

    Expr** args = malloc(sizeof(Expr*) * (d + 1));
    for (int i = 0; i <= d; i++) {
        Expr* c = bulk ? coeffs[i] : get_coeff_expanded(expanded, var, i);
        args[i] = coeff_list_rec(c, vars, max_degrees, num_vars, var_idx + 1);
        expr_free(c);
    }
    free(coeffs);
    expr_free(expanded);

    Expr* list = expr_new_function(expr_new_symbol("List"), args, d + 1);
    free(args);
    return list;
}

Expr* builtin_coefficientlist(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* vars_expr = res->data.function.args[1];

    if (expr->type == EXPR_INTEGER && expr->data.integer == 0) {
        return expr_new_function(expr_new_symbol("List"), NULL, 0);
    }

    size_t num_vars = 1;
    Expr** vars = &vars_expr;

    if (vars_expr->type == EXPR_FUNCTION && vars_expr->data.function.head->type == EXPR_SYMBOL && vars_expr->data.function.head->data.symbol == SYM_List) {
        num_vars = vars_expr->data.function.arg_count;
        vars = vars_expr->data.function.args;
        if (num_vars == 0) return expr_copy(expr);
    }

    Expr* expanded = expr_expand(expr);
    if (!expanded) return expr_copy(expr);

    int* max_degrees = malloc(sizeof(int) * num_vars);
    for (size_t i = 0; i < num_vars; i++) {
        int d = get_degree_poly(expanded, vars[i]);
        max_degrees[i] = (d < 0) ? 0 : d;
    }

    Expr* result = coeff_list_rec(expanded, vars, max_degrees, num_vars, 0);
    free(max_degrees);
    expr_free(expanded);
    return result;
}

/* Functional decomposition of a univariate polynomial:                */
/* given f(x), find polynomials g, h with f = g(h(x)) and deg(h) >= 2. */
/* The recursive driver applies two reductions:                        */
/*   - Common-degree shortcut: if every nonzero monomial has degree    */
/*     divisible by `d`, substitute y = x^d.                            */
/*   - Trial composition: for each divisor s of n=deg(f), build the    */
/*     candidate inner polynomial of degree s using the expansion of   */
/*     a_n * (x^s)^r and check whether the outer polynomial divides    */
/*     evenly (matches the standard "Kozen-Landau" decomposition).     */
static Expr* decompose_recursive(Expr* f, Expr* x) {
    Expr* expanded = expr_expand(f);
    int n = get_degree_poly(expanded, x);
    if (n < 2) {
        Expr* res = expr_new_function(expr_new_symbol("List"), (Expr*[]){expr_copy(expanded)}, 1);
        expr_free(expanded);
        return res;
    }

    /* Compute the gcd `d` of all i for which the i-th coefficient is     */
    /* non-zero. If d > 1, every term has degree divisible by d and we    */
    /* can substitute y = x^d to obtain a smaller polynomial.             */
    Expr** all_coeffs = NULL;
    bool have_bulk = get_all_coeffs_expanded(expanded, x, n, &all_coeffs);

    int d = 0;
    for (int i = 1; i <= n; i++) {
        Expr* c = have_bulk ? all_coeffs[i] : get_coeff_expanded(expanded, x, i);
        if (!is_zero_poly(c)) {
            if (d == 0) d = i;
            else {
                int64_t tmp_d = gcd(d, i);
                d = (int)tmp_d;
            }
        }
        if (!have_bulk) expr_free(c);
    }

    if (d > 1) {
        Expr* H = internal_power((Expr*[]){expr_copy(x), expr_new_integer(d)}, 2);

        Expr** g_args = malloc(sizeof(Expr*) * (n/d + 1));
        int g_count = 0;
        for (int i = 0; i <= n; i += d) {
            Expr* c = have_bulk ? expr_copy(all_coeffs[i]) : get_coeff_expanded(expanded, x, i);
            if (!is_zero_poly(c)) {
                if (i == 0) {
                    g_args[g_count++] = c;
                } else if (i == d) {
                    Expr* t = internal_times((Expr*[]){c, expr_copy(x)}, 2);
                    g_args[g_count++] = t;
                } else {
                    Expr* xp = internal_power((Expr*[]){expr_copy(x), expr_new_integer(i/d)}, 2);
                    Expr* t = internal_times((Expr*[]){c, xp}, 2);
                    g_args[g_count++] = t;
                }
            } else {
                expr_free(c);
            }
        }
        Expr* g;
        if (g_count == 0) g = expr_new_integer(0);
        else if (g_count == 1) g = g_args[0];
        else g = internal_plus(g_args, g_count);
        free(g_args);

        if (have_bulk) {
            for (int i = 0; i <= n; i++) expr_free(all_coeffs[i]);
            free(all_coeffs);
        }
        expr_free(expanded);

        if (expr_eq(g, x)) {
            expr_free(g);
            return expr_new_function(expr_new_symbol("List"), (Expr*[]){H}, 1);
        }

        Expr* Lg = decompose_recursive(g, x);
        expr_free(g);
        
        size_t Lg_count = Lg->data.function.arg_count;
        Expr** L_args = malloc(sizeof(Expr*) * (Lg_count + 1));
        for (size_t i = 0; i < Lg_count; i++) L_args[i] = expr_copy(Lg->data.function.args[i]);
        L_args[Lg_count] = H;
        
        Expr* res = expr_new_function(expr_new_symbol("List"), L_args, Lg_count + 1);
        free(L_args);
        expr_free(Lg);
        return res;
    }

    Expr* a_n = have_bulk ? expr_copy(all_coeffs[n]) : get_coeff_expanded(expanded, x, n);
    if (have_bulk) {
        for (int i = 0; i <= n; i++) expr_free(all_coeffs[i]);
        free(all_coeffs);
        all_coeffs = NULL;
        have_bulk = false;
    }
    for (int s = 2; s < n; s++) {
        if (n % s != 0) continue;
        int r = n / s;
        
        Expr* H = internal_power((Expr*[]){expr_copy(x), expr_new_integer(s)}, 2);
        
        bool valid = true;
        for (int k = 1; k < s; k++) {
            Expr* temp_E = internal_power((Expr*[]){expr_copy(H), expr_new_integer(r)}, 2);
            Expr* E = expr_expand(temp_E);
            expr_free(temp_E);
            Expr* C = get_coeff(E, x, n - k);
            expr_free(E);
            
            Expr* a_nk = get_coeff(expanded, x, n - k);
            
            Expr* temp_num = internal_plus((Expr*[]){a_nk, internal_times((Expr*[]){expr_new_integer(-1), expr_copy(a_n), C}, 3)}, 2);
            Expr* num = expr_expand(temp_num);
            expr_free(temp_num);
            
            Expr* temp_den = internal_times((Expr*[]){expr_new_integer(r), expr_copy(a_n)}, 2);
            Expr* den = expr_expand(temp_den);
            expr_free(temp_den);
            
            size_t v_count = 0, v_cap = 16;
            Expr** vars = malloc(sizeof(Expr*) * v_cap);
            collect_variables(num, &vars, &v_count, &v_cap);
            collect_variables(den, &vars, &v_count, &v_cap);
            if (v_count > 0) qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);
            
            size_t vx_count = 0;
            Expr** vars_nox = malloc(sizeof(Expr*) * v_count);
            for (size_t i = 0; i < v_count; i++) {
                if (!expr_eq(vars[i], x)) vars_nox[vx_count++] = vars[i];
            }
            
            Expr* c_sk = exact_poly_div(num, den, vars_nox, vx_count);
            
            for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
            free(vars);
            free(vars_nox);
            expr_free(num);
            expr_free(den);
            
            if (!c_sk) {
                valid = false;
                break;
            }
            
            Expr* term;
            if (s - k == 1) {
                term = internal_times((Expr*[]){c_sk, expr_copy(x)}, 2);
            } else {
                Expr* xp = internal_power((Expr*[]){expr_copy(x), expr_new_integer(s - k)}, 2);
                term = internal_times((Expr*[]){c_sk, xp}, 2);
            }
            Expr* temp_H = internal_plus((Expr*[]){H, term}, 2);
            Expr* next_H = expr_expand(temp_H);
            expr_free(temp_H);
            H = next_H;
        }
        
        if (valid) {
            Expr* Q = expr_copy(expanded);
            Expr** g_terms = malloc(sizeof(Expr*) * (r + 1));
            int g_count = 0;
            for (int i = 0; i <= r; i++) {
                Expr *new_Q, *Rem;
                poly_div_rem(Q, H, x, &new_Q, &Rem);
                expr_free(Q);
                Q = new_Q;
                if (get_degree_poly(Rem, x) > 0) {
                    expr_free(Rem);
                    valid = false;
                    break;
                }
                g_terms[g_count++] = Rem;
            }
            
            if (valid && is_zero_poly(Q)) {
                expr_free(Q);
                
                Expr** g_args = malloc(sizeof(Expr*) * g_count);
                int actual_g_count = 0;
                for (int i = 0; i < g_count; i++) {
                    if (!is_zero_poly(g_terms[i])) {
                        if (i == 0) {
                            g_args[actual_g_count++] = expr_copy(g_terms[i]);
                        } else if (i == 1) {
                            g_args[actual_g_count++] = internal_times((Expr*[]){expr_copy(g_terms[i]), expr_copy(x)}, 2);
                        } else {
                            Expr* xp = internal_power((Expr*[]){expr_copy(x), expr_new_integer(i)}, 2);
                            g_args[actual_g_count++] = internal_times((Expr*[]){expr_copy(g_terms[i]), xp}, 2);
                        }
                    }
                    expr_free(g_terms[i]);
                }
                free(g_terms);
                
                Expr* g;
                if (actual_g_count == 0) g = expr_new_integer(0);
                else if (actual_g_count == 1) g = g_args[0];
                else g = internal_plus(g_args, actual_g_count);
                free(g_args);
                
                Expr* Lg = decompose_recursive(g, x);
                Expr* Lh = decompose_recursive(H, x);
                expr_free(g);
                expr_free(H);
                expr_free(expanded);
                expr_free(a_n);
                
                size_t c1 = Lg->data.function.arg_count;
                size_t c2 = Lh->data.function.arg_count;
                Expr** final_args = malloc(sizeof(Expr*) * (c1 + c2));
                for (size_t i = 0; i < c1; i++) final_args[i] = expr_copy(Lg->data.function.args[i]);
                for (size_t i = 0; i < c2; i++) final_args[c1 + i] = expr_copy(Lh->data.function.args[i]);
                
                Expr* res = expr_new_function(expr_new_symbol("List"), final_args, c1 + c2);
                free(final_args);
                expr_free(Lg);
                expr_free(Lh);
                return res;
            } else {
                for (int i = 0; i < g_count; i++) expr_free(g_terms[i]);
                free(g_terms);
                expr_free(Q);
            }
        }
        expr_free(H);
    }
    expr_free(a_n);
    
    Expr* res = expr_new_function(expr_new_symbol("List"), (Expr*[]){expr_copy(expanded)}, 1);
    expr_free(expanded);
    return res;
}

Expr* builtin_decompose(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* poly = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    
    Expr* pq = internal_polynomialq((Expr*[]){expr_copy(poly), expr_copy(x)}, 2);
    bool is_poly = (pq->type == EXPR_SYMBOL && pq->data.symbol == SYM_True);
    expr_free(pq);
    
    if (!is_poly) {
        return NULL;
    }
    
    return decompose_recursive(poly, x);
}

/* Recursive Horner-form construction:                                  */
/*   p(x) = c_0 + c_1 x + ... + c_n x^n                                 */
/*       => c_0 + x (c_1 + x (c_2 + ... + x c_n))                      */
/* For multivariate inputs, recurses into each coefficient using the   */
/* tail of the variable list.                                          */
static Expr* horner_form_rec(Expr* expr, Expr** vars, size_t num_vars) {
    if (num_vars == 0) return expr_copy(expr);
    Expr* v = vars[0];
    
    Expr* expanded = expr_expand(expr);
    
    Expr* pq = internal_polynomialq((Expr*[]){expr_copy(expanded), expr_copy(v)}, 2);
    bool is_poly = (pq->type == EXPR_SYMBOL && pq->data.symbol == SYM_True);
    expr_free(pq);
    
    if (!is_poly) {
        expr_free(expanded);
        return NULL;
    }
    
    Expr* cl = eval_and_free(expr_new_function(expr_new_symbol("CoefficientList"), (Expr*[]){expr_copy(expanded), expr_copy(v)}, 2));
    expr_free(expanded);
    
    if (!cl || cl->type != EXPR_FUNCTION || cl->data.function.head->data.symbol != SYM_List) {
        if (cl) expr_free(cl);
        return NULL;
    }
    
    size_t count = cl->data.function.arg_count;
    if (count == 0) {
        expr_free(cl);
        return expr_new_integer(0);
    }
    
    Expr* H = horner_form_rec(cl->data.function.args[count - 1], vars + 1, num_vars - 1);
    if (!H) {
        expr_free(cl);
        return NULL;
    }
    
    for (int i = (int)count - 2; i >= 0; i--) {
        Expr* c_i = horner_form_rec(cl->data.function.args[i], vars + 1, num_vars - 1);
        if (!c_i) {
            expr_free(cl);
            expr_free(H);
            return NULL;
        }
        
        bool h_zero = is_zero_poly(H);
        
        if (h_zero) {
            expr_free(H);
            H = c_i;
        } else {
            Expr* t = internal_times((Expr*[]){expr_copy(v), H}, 2);
            bool c_zero = is_zero_poly(c_i);
            if (c_zero) {
                expr_free(c_i);
                H = t;
            } else {
                H = internal_plus((Expr*[]){c_i, t}, 2);
            }
        }
    }
    
    expr_free(cl);
    return H;
}

Expr* builtin_hornerform(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1 || res->data.function.arg_count > 3) return NULL;
    
    Expr* expr = res->data.function.args[0];
    
    Expr* num = NULL;
    Expr* den = NULL;
    
    if (expr->type == EXPR_FUNCTION && expr->data.function.head->type == EXPR_SYMBOL && expr->data.function.head->data.symbol == SYM_Times) {
        size_t n_cap = 16, n_count = 0;
        size_t d_cap = 16, d_count = 0;
        Expr** n_args = malloc(sizeof(Expr*) * n_cap);
        Expr** d_args = malloc(sizeof(Expr*) * d_cap);
        
        for (size_t i = 0; i < expr->data.function.arg_count; i++) {
            Expr* arg = expr->data.function.args[i];
            if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && arg->data.function.head->data.symbol == SYM_Power && arg->data.function.arg_count == 2) {
                Expr* exp = arg->data.function.args[1];
                if ((exp->type == EXPR_INTEGER && exp->data.integer < 0) || 
                    (exp->type == EXPR_FUNCTION && exp->data.function.head->type == EXPR_SYMBOL && exp->data.function.head->data.symbol == SYM_Rational && exp->data.function.args[0]->data.integer < 0)) {
                    if (d_count == d_cap) { d_cap *= 2; d_args = realloc(d_args, sizeof(Expr*) * d_cap); }
                    if (exp->type == EXPR_INTEGER) {
                        if (exp->data.integer == -1) {
                            d_args[d_count++] = expr_copy(arg->data.function.args[0]);
                        } else {
                            d_args[d_count++] = internal_power((Expr*[]){expr_copy(arg->data.function.args[0]), expr_new_integer(-exp->data.integer)}, 2);
                        }
                    } else { 
                        Expr* new_rat = eval_and_free(expr_new_function(expr_new_symbol("Rational"), (Expr*[]){expr_new_integer(-exp->data.function.args[0]->data.integer), expr_copy(exp->data.function.args[1])}, 2));
                        d_args[d_count++] = internal_power((Expr*[]){expr_copy(arg->data.function.args[0]), new_rat}, 2);
                    }
                    continue;
                }
            }
            if (n_count == n_cap) { n_cap *= 2; n_args = realloc(n_args, sizeof(Expr*) * n_cap); }
            n_args[n_count++] = expr_copy(arg);
        }
        
        if (n_count == 0) num = expr_new_integer(1);
        else if (n_count == 1) num = n_args[0];
        else num = internal_times(n_args, n_count);
        
        if (d_count == 0) den = expr_new_integer(1);
        else if (d_count == 1) den = d_args[0];
        else den = internal_times(d_args, d_count);
        
        free(n_args); free(d_args);
    } else if (expr->type == EXPR_FUNCTION && expr->data.function.head->type == EXPR_SYMBOL && expr->data.function.head->data.symbol == SYM_Power && expr->data.function.arg_count == 2) {
        Expr* exp = expr->data.function.args[1];
        if (exp->type == EXPR_INTEGER && exp->data.integer < 0) {
            num = expr_new_integer(1);
            if (exp->data.integer == -1) {
                den = expr_copy(expr->data.function.args[0]);
            } else {
                den = internal_power((Expr*[]){expr_copy(expr->data.function.args[0]), expr_new_integer(-exp->data.integer)}, 2);
            }
        } else {
            num = expr_copy(expr);
            den = expr_new_integer(1);
        }
    } else {
        num = expr_copy(expr);
        den = expr_new_integer(1);
    }

    Expr* vars1_expr = NULL;
    Expr* vars2_expr = NULL;
    
    if (res->data.function.arg_count == 1) {
        vars1_expr = eval_and_free(expr_new_function(expr_new_symbol("Variables"), (Expr*[]){expr_copy(num)}, 1));
        vars2_expr = eval_and_free(expr_new_function(expr_new_symbol("Variables"), (Expr*[]){expr_copy(den)}, 1));
    } else if (res->data.function.arg_count == 2) {
        vars1_expr = expr_copy(res->data.function.args[1]);
        vars2_expr = expr_copy(res->data.function.args[1]);
    } else if (res->data.function.arg_count == 3) {
        vars1_expr = expr_copy(res->data.function.args[1]);
        vars2_expr = expr_copy(res->data.function.args[2]);
    }
    
    if (vars1_expr && (vars1_expr->type != EXPR_FUNCTION || vars1_expr->data.function.head->type != EXPR_SYMBOL || vars1_expr->data.function.head->data.symbol != SYM_List)) {
        vars1_expr = eval_and_free(expr_new_function(expr_new_symbol("List"), (Expr*[]){vars1_expr}, 1));
    }
    if (vars2_expr && (vars2_expr->type != EXPR_FUNCTION || vars2_expr->data.function.head->type != EXPR_SYMBOL || vars2_expr->data.function.head->data.symbol != SYM_List)) {
        vars2_expr = eval_and_free(expr_new_function(expr_new_symbol("List"), (Expr*[]){vars2_expr}, 1));
    }
    
    Expr** vars1 = vars1_expr ? vars1_expr->data.function.args : NULL;
    size_t num_vars1 = vars1_expr ? vars1_expr->data.function.arg_count : 0;
    
    Expr** vars2 = vars2_expr ? vars2_expr->data.function.args : NULL;
    size_t num_vars2 = vars2_expr ? vars2_expr->data.function.arg_count : 0;
    
    Expr* h_num = horner_form_rec(num, vars1, num_vars1);
    if (!h_num) {
        printf("HornerForm::poly: "); 
        char* s = expr_to_string(expr);
        printf("%s is not a polynomial.\n", s);
        free(s);
        if (vars1_expr) expr_free(vars1_expr);
        if (vars2_expr) expr_free(vars2_expr);
        expr_free(num);
        expr_free(den);
        return expr_copy(res);
    }
    
    Expr* h_den = NULL;
    if (den->type == EXPR_INTEGER && den->data.integer == 1) {
        h_den = expr_copy(den);
    } else {
        h_den = horner_form_rec(den, vars2, num_vars2);
        if (!h_den) {
            printf("HornerForm::poly: "); 
            char* s = expr_to_string(expr);
            printf("%s is not a polynomial.\n", s);
            free(s);
            if (vars1_expr) expr_free(vars1_expr);
            if (vars2_expr) expr_free(vars2_expr);
            expr_free(num);
            expr_free(den);
            expr_free(h_num);
            return expr_copy(res);
        }
    }
    
    Expr* result = NULL;
    if (h_den->type == EXPR_INTEGER && h_den->data.integer == 1) {
        result = h_num;
        expr_free(h_den);
    } else {
        Expr* inv_den = internal_power((Expr*[]){h_den, expr_new_integer(-1)}, 2);
        result = internal_times((Expr*[]){h_num, inv_den}, 2);
    }
    
    if (vars1_expr) expr_free(vars1_expr);
    if (vars2_expr) expr_free(vars2_expr);
    expr_free(num);
    expr_free(den);
    
    return result;
}

/* ------------------------------------------------------------------ */
/* Bronstein-style subresultant PRS for Resultant.                    */
/*   See "SubResultant" in Bronstein, Symbolic Integration I, p.24.   */
/*                                                                    */
/* This computes Resultant(A, B, x) without ever forming the          */
/* (n+m)x(n+m) Sylvester matrix.  The chain works in D[x] (D the      */
/* coefficient ring), making only O(min(n, m)) calls to a pseudo-     */
/* remainder, with a single *scalar* exact division per chain step    */
/* (R / β_i, β_i ∈ D, not D[x]).  In contrast, the Sylvester+Det path */
/* needs O(n^3) exact divisions of polynomial coefficients in Bareiss */
/* elimination, which collapses to O(n!) Laplace expansion when any   */
/* one of those divisions cannot be certified — most commonly over    */
/* algebraic-number coefficient rings like Q(α).                       */
/* ------------------------------------------------------------------ */

/* Local leaf-count helper for size-based bailout in the chain.       */
static int64_t subres_leaf_count(Expr* e) {
    if (!e) return 0;
    if (e->type != EXPR_FUNCTION) return 1;
    int64_t c = 1;  /* count the head */
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        c += subres_leaf_count(e->data.function.args[i]);
    }
    return c;
}

/* True if `e` contains any subterm Power[X, Rational[a, b]] with b > 1, */
/* i.e. an algebraic number (Sqrt[N], cube roots, etc.).  Pseudo-       */
/* remainder over D = Q(α)[t] generates Power[base, k/2] forms (e.g.    */
/* Sqrt[3]^3 → 3^(3/2)) that don't combine with Times[base, Sqrt[base]] */
/* via Plus alone, so the chain bloats geometrically.  We detect this   */
/* up front and let the Sylvester+Det path handle it instead.            */
static bool subres_has_algebraic(Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Power &&
        e->data.function.arg_count == 2) {
        Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_FUNCTION &&
            exp->data.function.head->type == EXPR_SYMBOL &&
            exp->data.function.head->data.symbol == SYM_Rational &&
            exp->data.function.arg_count == 2) {
            Expr* den = exp->data.function.args[1];
            if (den->type == EXPR_INTEGER && den->data.integer != 1) {
                return true;
            }
            if (den->type == EXPR_BIGINT) return true;
        }
    }
    if (e->data.function.head &&
        subres_has_algebraic(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (subres_has_algebraic(e->data.function.args[i])) return true;
    }
    return false;
}

/* Standard pseudo-remainder: returns lc(B)^(deg(A) - deg(B) + 1) * A */
/* mod B in D[x], assuming deg(A) >= deg(B).  Distinguished from the  */
/* looser pseudo_rem above (which multiplies by the minimum power of  */
/* lc(B) needed to drive the remainder below deg(B)) — Bronstein's    */
/* chain identities require the full d+1 power so that division by    */
/* β_i is exact in D.                                                  */
static Expr* pseudo_rem_standard(Expr* A, Expr* B, Expr* x) {
    Expr* expandedB = expr_expand(B);
    int dB = get_degree_poly(expandedB, x);
    Expr* expandedA = expr_expand(A);
    int dA = get_degree_poly(expandedA, x);

    if (is_zero_poly(expandedB) || dA < dB) {
        expr_free(expandedB);
        return expandedA;
    }

    Expr* lcB = get_coeff_expanded(expandedB, x, dB);
    int expected_iters = dA - dB + 1;
    int iters = 0;
    Expr* R = expandedA;

    while (true) {
        int degR = get_degree_poly(R, x);
        if (degR < dB || is_zero_poly(R)) break;

        Expr* lcR = get_coeff_expanded(R, x, degR);
        int d = degR - dB;

        Expr* t1 = internal_times((Expr*[]){expr_copy(lcB), R}, 2);
        Expr* x_pow = internal_power(
            (Expr*[]){expr_copy(x), expr_new_integer(d)}, 2);
        Expr* t2 = internal_times(
            (Expr*[]){lcR, x_pow, expr_copy(expandedB)}, 3);
        Expr* neg_t2 = internal_times(
            (Expr*[]){expr_new_integer(-1), t2}, 2);
        Expr* diff = internal_plus((Expr*[]){t1, neg_t2}, 2);

        R = expr_expand(diff);
        expr_free(diff);
        iters++;
    }

    /* Pad with extra lc(B) multiplications when an iteration step      */
    /* dropped the remainder degree by more than 1 (so we exited early).*/
    if (iters < expected_iters) {
        Expr* pad_pow = internal_power(
            (Expr*[]){expr_copy(lcB),
                      expr_new_integer(expected_iters - iters)}, 2);
        Expr* padded = internal_times((Expr*[]){R, pad_pow}, 2);
        R = expr_expand(padded);
        expr_free(padded);
    }

    expr_free(lcB);
    expr_free(expandedB);
    return R;
}

/* Coefficient-wise division of `poly` (a polynomial in `var`) by the */
/* scalar `denom`.  Builds the result by extracting coefficients,      */
/* dividing each by denom, and rebuilding.  Cancel handles rational   */
/* simplification.  Returns a fresh expression; the caller frees it.  */
static Expr* poly_divide_by_scalar(Expr* poly, Expr* denom, Expr* var) {
    if (is_zero_poly(poly)) return expr_new_integer(0);

    Expr* expanded = expr_expand(poly);
    int d = get_degree_poly(expanded, var);
    if (d < 0 || is_zero_poly(expanded)) {
        expr_free(expanded);
        return expr_new_integer(0);
    }

    Expr* inv_denom = internal_power(
        (Expr*[]){expr_copy(denom), expr_new_integer(-1)}, 2);

    Expr** coeffs = NULL;
    bool bulk = get_all_coeffs_expanded(expanded, var, d, &coeffs);

    Expr** terms = malloc(sizeof(Expr*) * (d + 1));
    int term_count = 0;
    for (int i = 0; i <= d; i++) {
        Expr* c = bulk ? coeffs[i] : get_coeff_expanded(expanded, var, i);
        if (is_zero_poly(c)) {
            expr_free(c);
            continue;
        }
        Expr* prod = internal_times((Expr*[]){c, expr_copy(inv_denom)}, 2);
        Expr* simp = internal_cancel((Expr*[]){prod}, 1);
        if (i == 0) {
            terms[term_count++] = simp;
        } else if (i == 1) {
            terms[term_count++] = internal_times(
                (Expr*[]){simp, expr_copy(var)}, 2);
        } else {
            Expr* xp = internal_power(
                (Expr*[]){expr_copy(var), expr_new_integer(i)}, 2);
            terms[term_count++] = internal_times(
                (Expr*[]){simp, xp}, 2);
        }
    }
    free(coeffs);
    expr_free(inv_denom);
    expr_free(expanded);

    Expr* sum;
    if (term_count == 0) {
        free(terms);
        return expr_new_integer(0);
    }
    if (term_count == 1) {
        sum = terms[0];
    } else {
        sum = internal_plus(terms, term_count);
    }
    free(terms);
    Expr* result = expr_expand(sum);
    expr_free(sum);
    return result;
}

/* Bronstein subresultant Resultant.  Returns NULL on hard failure so */
/* the caller can fall back to the Sylvester+Det path.  Both inputs   */
/* are required to have positive degree in `var` — degenerate cases   */
/* (constant inputs) are handled by resultant_internal before calling */
/* this.                                                               */
static Expr* resultant_subresultant(Expr* P, Expr* Q, Expr* var) {
    /* Skip Bronstein for inputs with algebraic numbers (Sqrt[N], cube     */
    /* roots, etc.).  Pseudo-remainder over Q(α)[t][x] generates           */
    /* Power[α, k/m] forms that our Plus doesn't combine with the          */
    /* equivalent Times[α^q, Sqrt[α]] forms, leading to geometric chain    */
    /* bloat.  Sylvester+Det is more reliable here even at O(n!) Laplace.  */
    if (subres_has_algebraic(P) || subres_has_algebraic(Q)) {
        return NULL;
    }

    int dP = get_degree_poly(P, var);
    int dQ = get_degree_poly(Q, var);

    /* Order chain so deg(R_0) >= deg(R_1).  Theorem 1.4.1: swapping  */
    /* arguments multiplies the resultant by (-1)^(deg(P)*deg(Q)).    */
    int swap_sign = 1;
    Expr *A_init, *B_init;
    int dA_init, dB_init;
    if (dP >= dQ) {
        A_init = expr_expand(P);
        B_init = expr_expand(Q);
        dA_init = dP; dB_init = dQ;
    } else {
        A_init = expr_expand(Q);
        B_init = expr_expand(P);
        dA_init = dQ; dB_init = dP;
        if ((dP & 1) && (dQ & 1)) swap_sign = -1;
    }

    /* Chain storage. */
    size_t cap = 8;
    Expr** R_chain = (Expr**)malloc(sizeof(Expr*) * cap);
    int*   degs    = (int*)  malloc(sizeof(int)   * cap);
    Expr** lcs     = (Expr**)malloc(sizeof(Expr*) * cap);
    Expr** betas   = (Expr**)malloc(sizeof(Expr*) * cap);
    for (size_t t = 0; t < cap; t++) {
        R_chain[t] = NULL; lcs[t] = NULL; betas[t] = NULL; degs[t] = -1;
    }

    R_chain[0] = A_init;
    R_chain[1] = B_init;
    degs[0] = dA_init;
    degs[1] = dB_init;
    lcs[0] = NULL;  /* never read */
    lcs[1] = get_coeff_expanded(B_init, var, dB_init);

    int delta = dA_init - dB_init;          /* δ_1 */
    Expr* gamma = expr_new_integer(-1);     /* γ_1 */
    /* β_1 = (-1)^(δ_1 + 1). */
    betas[1] = ((delta + 1) & 1)
                ? expr_new_integer(-1)
                : expr_new_integer(1);

    int prev_delta = delta;
    int i = 1;
    Expr* result = NULL;
    bool failed = false;

    /* Size budget: if a chain element exceeds this threshold relative to */
    /* the inputs, the cancellation in poly_divide_by_scalar isn't       */
    /* keeping pace (typically Q(α)[t] coefficient rings where           */
    /* internal_cancel can't fully reduce algebraic-number denominators). */
    /* Bail out so resultant_internal falls back to Sylvester+Det.        */
    int64_t input_size = subres_leaf_count(A_init) + subres_leaf_count(B_init);
    int64_t size_budget = input_size * 30;
    if (size_budget < 5000) size_budget = 5000;

    while (!is_zero_poly(R_chain[i])) {
        /* Pre-check: bail before doing the next (potentially expensive) */
        /* pseudo-remainder if the current chain head is already bloated. */
        if (subres_leaf_count(R_chain[i]) > size_budget) {
            failed = true;
            break;
        }

        Expr* pprem_val = pseudo_rem_standard(R_chain[i-1], R_chain[i], var);
        if (!pprem_val) { failed = true; break; }

        Expr* R_next = poly_divide_by_scalar(pprem_val, betas[i], var);
        expr_free(pprem_val);
        if (!R_next) { failed = true; break; }

        if (subres_leaf_count(R_next) > size_budget) {
            expr_free(R_next);
            failed = true;
            break;
        }

        /* Grow arrays if needed (need slot i+2 reachable below). */
        if ((size_t)(i + 2) >= cap) {
            size_t new_cap = cap * 2;
            R_chain = realloc(R_chain, sizeof(Expr*) * new_cap);
            degs    = realloc(degs,    sizeof(int)   * new_cap);
            lcs     = realloc(lcs,     sizeof(Expr*) * new_cap);
            betas   = realloc(betas,   sizeof(Expr*) * new_cap);
            for (size_t t = cap; t < new_cap; t++) {
                R_chain[t] = NULL; lcs[t] = NULL;
                betas[t] = NULL; degs[t] = -1;
            }
            cap = new_cap;
        }

        R_chain[i+1] = R_next;
        if (is_zero_poly(R_next)) {
            i = i + 1;
            break;
        }
        degs[i+1] = get_degree_poly(R_next, var);
        lcs[i+1] = get_coeff_expanded(R_next, var, degs[i+1]);

        /* Advance index, then update γ_i, δ_i, β_i. */
        i = i + 1;

        /* γ_i = (-r_{i-1})^(δ_{i-1}) · γ_{i-1}^(1 - δ_{i-1}) */
        Expr* neg_r = internal_times(
            (Expr*[]){expr_new_integer(-1), expr_copy(lcs[i-1])}, 2);
        Expr* term1 = internal_power(
            (Expr*[]){neg_r, expr_new_integer(prev_delta)}, 2);
        Expr* term2 = internal_power(
            (Expr*[]){expr_copy(gamma), expr_new_integer(1 - prev_delta)}, 2);
        Expr* gamma_new = internal_times((Expr*[]){term1, term2}, 2);
        Expr* gamma_simp = internal_cancel((Expr*[]){gamma_new}, 1);
        Expr* gamma_expanded = expr_expand(gamma_simp);
        expr_free(gamma_simp);
        expr_free(gamma);
        gamma = gamma_expanded;

        /* δ_i = deg(R_{i-1}) - deg(R_i). */
        int new_delta = degs[i-1] - degs[i];
        prev_delta = new_delta;

        /* β_i = -r_{i-1} · γ_i^(δ_i). */
        Expr* gamma_pow = internal_power(
            (Expr*[]){expr_copy(gamma), expr_new_integer(new_delta)}, 2);
        Expr* beta_new = internal_times(
            (Expr*[]){expr_new_integer(-1),
                      expr_copy(lcs[i-1]), gamma_pow}, 3);
        Expr* beta_simp = internal_cancel((Expr*[]){beta_new}, 1);
        betas[i] = expr_expand(beta_simp);
        expr_free(beta_simp);
    }

    int k = i - 1;

    if (!failed) {
        if (k == 0) {
            /* Should not occur for positive-degree inputs. */
            failed = true;
        } else if (degs[k] > 0) {
            /* Common factor of positive degree -> resultant is 0. */
            result = expr_new_integer(0);
        } else if (degs[k-1] == 1) {
            /* Fast path: resultant = R_k (a non-zero constant). */
            result = expr_copy(R_chain[k]);
        } else {
            /* General case: result = s · c · R_k^deg(R_{k-1}). */
            int s = 1;
            Expr* c = expr_new_integer(1);
            for (int j = 1; j <= k - 1; j++) {
                if ((degs[j-1] & 1) && (degs[j] & 1)) s = -s;
                int delta_j = degs[j-1] - degs[j];
                int exp_ratio = degs[j];
                int exp_rj = degs[j-1] - degs[j+1];

                Expr* rj_1pdj = internal_power(
                    (Expr*[]){expr_copy(lcs[j]),
                              expr_new_integer(1 + delta_j)}, 2);
                Expr* inv_rj_1pdj = internal_power(
                    (Expr*[]){rj_1pdj, expr_new_integer(-1)}, 2);
                Expr* ratio = internal_times(
                    (Expr*[]){expr_copy(betas[j]), inv_rj_1pdj}, 2);
                Expr* ratio_pow = internal_power(
                    (Expr*[]){ratio, expr_new_integer(exp_ratio)}, 2);
                Expr* rj_pow_2 = internal_power(
                    (Expr*[]){expr_copy(lcs[j]),
                              expr_new_integer(exp_rj)}, 2);
                Expr* incr = internal_times(
                    (Expr*[]){c, ratio_pow, rj_pow_2}, 3);
                Expr* simp = internal_cancel((Expr*[]){incr}, 1);
                c = expr_expand(simp);
                expr_free(simp);
            }

            Expr* Rk_pow = internal_power(
                (Expr*[]){expr_copy(R_chain[k]),
                          expr_new_integer(degs[k-1])}, 2);
            Expr* sc = (s == 1)
                ? c
                : internal_times((Expr*[]){expr_new_integer(-1), c}, 2);
            Expr* unsigned_res = internal_times((Expr*[]){sc, Rk_pow}, 2);
            result = expr_expand(unsigned_res);
            expr_free(unsigned_res);
        }
        if (result && swap_sign == -1) {
            Expr* neg = internal_times(
                (Expr*[]){expr_new_integer(-1), result}, 2);
            result = expr_expand(neg);
            expr_free(neg);
        }
    }

    /* Cleanup chain storage. */
    for (size_t t = 0; t < cap; t++) {
        if (R_chain[t]) expr_free(R_chain[t]);
        if (lcs[t])     expr_free(lcs[t]);
        if (betas[t])   expr_free(betas[t]);
    }
    free(R_chain);
    free(degs);
    free(lcs);
    free(betas);
    expr_free(gamma);

    return result;
}

/* Resultant of P, Q in `var`. We exploit two algebraic identities for */
/* a fast path before falling back to the Sylvester matrix:            */
/*   Res(P1*P2, Q) = Res(P1,Q) * Res(P2,Q)                              */
/*   Res(P^k, Q)   = Res(P,Q)^k                                         */
/* The general case prefers Bronstein's subresultant PRS                */
/* (resultant_subresultant) and falls back to constructing the          */
/* (n+m)x(n+m) Sylvester matrix and                                     */
/* takes its determinant (delegated to the linalg module).             */
static Expr* resultant_internal(Expr* P, Expr* Q, Expr* var) {
    if (P->type == EXPR_FUNCTION && P->data.function.head->type == EXPR_SYMBOL) {
        if (P->data.function.head->data.symbol == SYM_Times) {
            size_t count = P->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = resultant_internal(P->data.function.args[i], Q, var);
            }
            Expr* res = internal_times(args, count);
            free(args);
            return res;
        } else if (P->data.function.head->data.symbol == SYM_Power && P->data.function.arg_count == 2) {
            Expr* r = resultant_internal(P->data.function.args[0], Q, var);
            return internal_power((Expr*[]){r, expr_copy(P->data.function.args[1])}, 2);
        }
    }
    
    if (Q->type == EXPR_FUNCTION && Q->data.function.head->type == EXPR_SYMBOL) {
        if (Q->data.function.head->data.symbol == SYM_Times) {
            size_t count = Q->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = resultant_internal(P, Q->data.function.args[i], var);
            }
            Expr* res = internal_times(args, count);
            free(args);
            return res;
        } else if (Q->data.function.head->data.symbol == SYM_Power && Q->data.function.arg_count == 2) {
            Expr* r = resultant_internal(P, Q->data.function.args[0], var);
            return internal_power((Expr*[]){r, expr_copy(Q->data.function.args[1])}, 2);
        }
    }
    
    Expr* exp_P = expr_expand(P);
    Expr* exp_Q = expr_expand(Q);
    int n = get_degree_poly(exp_P, var);
    int m = get_degree_poly(exp_Q, var);
    
    if (n == 0 && m == 0) {
        expr_free(exp_P); expr_free(exp_Q);
        return expr_new_integer(1);
    }
    if (n == 0) {
        Expr* r = internal_power((Expr*[]){expr_copy(exp_P), expr_new_integer(m)}, 2);
        expr_free(exp_P); expr_free(exp_Q);
        return r;
    }
    if (m == 0) {
        Expr* r = internal_power((Expr*[]){expr_copy(exp_Q), expr_new_integer(n)}, 2);
        expr_free(exp_P); expr_free(exp_Q);
        return r;
    }

    /* Try Bronstein's subresultant PRS first.  Returns NULL for inputs  */
    /* with algebraic-number coefficients (Sqrt[N], etc.) so Sylvester+  */
    /* Det handles those, and on size-budget exhaustion in pathological  */
    /* cases.                                                             */
    {
        Expr* sub = resultant_subresultant(exp_P, exp_Q, var);
        if (sub) {
            expr_free(exp_P); expr_free(exp_Q);
            return sub;
        }
    }

    /* Both inputs are already expanded -- one bulk pass per polynomial   */
    /* gives us all coefficients in O(terms), instead of O(deg * terms)   */
    /* across (deg+1) separate get_coeff queries.                          */
    Expr** p_all = NULL;
    Expr** p_coeffs = malloc(sizeof(Expr*) * (n + 1));
    if (get_all_coeffs_expanded(exp_P, var, n, &p_all)) {
        for (int i = 0; i <= n; i++) p_coeffs[i] = p_all[n - i];
        free(p_all);
    } else {
        for (int i = 0; i <= n; i++) p_coeffs[i] = get_coeff_expanded(exp_P, var, n - i);
    }

    Expr** q_all = NULL;
    Expr** q_coeffs = malloc(sizeof(Expr*) * (m + 1));
    if (get_all_coeffs_expanded(exp_Q, var, m, &q_all)) {
        for (int i = 0; i <= m; i++) q_coeffs[i] = q_all[m - i];
        free(q_all);
    } else {
        for (int i = 0; i <= m; i++) q_coeffs[i] = get_coeff_expanded(exp_Q, var, m - i);
    }
    
    int dim = n + m;
    Expr** rows = malloc(sizeof(Expr*) * dim);
    for (int i = 0; i < m; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * dim);
        for (int j = 0; j < dim; j++) {
            if (j >= i && j - i <= n) row_elems[j] = expr_copy(p_coeffs[j - i]);
            else row_elems[j] = expr_new_integer(0);
        }
        rows[i] = expr_new_function(expr_new_symbol("List"), row_elems, dim);
        free(row_elems);
    }
    
    for (int i = 0; i < n; i++) {
        Expr** row_elems = malloc(sizeof(Expr*) * dim);
        for (int j = 0; j < dim; j++) {
            if (j >= i && j - i <= m) row_elems[j] = expr_copy(q_coeffs[j - i]);
            else row_elems[j] = expr_new_integer(0);
        }
        rows[m + i] = expr_new_function(expr_new_symbol("List"), row_elems, dim);
        free(row_elems);
    }
    
    Expr* matrix = expr_new_function(expr_new_symbol("List"), rows, dim);
    free(rows);
    
    Expr* det_call = expr_new_function(expr_new_symbol("Det"), (Expr*[]){matrix}, 1);
    Expr* evaluated_det = evaluate(det_call);
    expr_free(det_call);
    
    Expr* result = expr_expand(evaluated_det);
    expr_free(evaluated_det);
    
    for (int i = 0; i <= n; i++) expr_free(p_coeffs[i]);
    free(p_coeffs);
    for (int i = 0; i <= m; i++) expr_free(q_coeffs[i]);
    free(q_coeffs);
    expr_free(exp_P); expr_free(exp_Q);
    
    return result;
}

Expr* builtin_resultant(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    
    Expr* p1 = res->data.function.args[0];
    Expr* p2 = res->data.function.args[1];
    Expr* var = res->data.function.args[2];
    
    Expr* pq1 = internal_polynomialq((Expr*[]){expr_copy(p1), expr_copy(var)}, 2);
    bool is_poly1 = (pq1->type == EXPR_SYMBOL && pq1->data.symbol == SYM_True);
    expr_free(pq1);
    
    Expr* pq2 = internal_polynomialq((Expr*[]){expr_copy(p2), expr_copy(var)}, 2);
    bool is_poly2 = (pq2->type == EXPR_SYMBOL && pq2->data.symbol == SYM_True);
    expr_free(pq2);
    
    if (!is_poly1 || !is_poly2) {
        return NULL;
    }
    
    return resultant_internal(p1, p2, var);
}

/* d/dvar of a univariate polynomial in already-expanded form. The     */
/* coefficient ring is left untouched -- we just multiply each c_i by  */
/* its exponent. Used by the discriminant routine.                     */
static Expr* poly_derivative(Expr* exp_p, Expr* var) {
    int n = get_degree_poly(exp_p, var);
    if (n <= 0) return expr_new_integer(0);

    Expr** coeffs = NULL;
    bool bulk = get_all_coeffs_expanded(exp_p, var, n, &coeffs);

    Expr** args = malloc(sizeof(Expr*) * n);
    int count = 0;
    for (int i = 1; i <= n; i++) {
        Expr* c = bulk ? coeffs[i] : get_coeff_expanded(exp_p, var, i);
        if (!is_zero_poly(c)) {
            Expr* t1 = internal_times((Expr*[]){expr_new_integer(i), c}, 2);
            if (i == 1) {
                args[count++] = t1;
            } else if (i == 2) {
                Expr* t2 = internal_times((Expr*[]){t1, expr_copy(var)}, 2);
                args[count++] = t2;
            } else {
                Expr* xp = internal_power((Expr*[]){expr_copy(var), expr_new_integer(i-1)}, 2);
                Expr* t2 = internal_times((Expr*[]){t1, xp}, 2);
                args[count++] = t2;
            }
        } else {
            expr_free(c);
        }
    }
    if (bulk) {
        /* coeffs[0] is the constant term — we never consumed it; free it.
         * coeffs[1..n] were either consumed into `args` or freed above via
         * `c` when zero. The container itself must also be freed. */
        expr_free(coeffs[0]);
        free(coeffs);
    }
    if (count == 0) {
        free(args);
        return expr_new_integer(0);
    } else if (count == 1) {
        Expr* res = args[0];
        free(args);
        return res;
    } else {
        Expr* res = internal_plus(args, count);
        free(args);
        return res;
    }
}

/* Discriminant[p, x] = (-1)^(n*(n-1)/2) / a_n * Resultant(p, p', x).   */
/* The sign and 1/a_n factor are applied after the resultant call.     */
Expr* builtin_discriminant(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    
    Expr* poly = res->data.function.args[0];
    Expr* var = res->data.function.args[1];
    
    Expr* pq = internal_polynomialq((Expr*[]){expr_copy(poly), expr_copy(var)}, 2);
    bool is_poly = (pq->type == EXPR_SYMBOL && pq->data.symbol == SYM_True);
    expr_free(pq);
    
    if (!is_poly) return NULL;
    
    Expr* exp_poly = expr_expand(poly);
    if (!exp_poly) { return NULL; }
    int n = get_degree_poly(exp_poly, var);
    
    if (n < 0) {
        expr_free(exp_poly);
        return NULL;
    }
    if (n == 0 || n == 1) {
        expr_free(exp_poly);
        return expr_new_integer(0); // Discriminant of constant or linear is 0 in some conventions, or 1? Mathematica says 1 for linear.
    }
    
    Expr* a_n = get_coeff(exp_poly, var, n);
    Expr* deriv = poly_derivative(exp_poly, var);
    if (!deriv) { expr_free(exp_poly); expr_free(a_n); return NULL; }
    Expr* res_val = resultant_internal(exp_poly, deriv, var);
    if (!res_val) { expr_free(exp_poly); expr_free(a_n); expr_free(deriv); return NULL; }
    expr_free(deriv);
    expr_free(exp_poly);
    
    int64_t sign_pow = (int64_t)n * (n - 1) / 2;
    int64_t sign = (sign_pow % 2 != 0) ? -1 : 1;
    
    Expr* num = internal_times((Expr*[]){expr_new_integer(sign), res_val}, 2);
    Expr* den = internal_power((Expr*[]){a_n, expr_new_integer(-1)}, 2);
    
    Expr* final_res = internal_times((Expr*[]){num, den}, 2);
    Expr* ret = expr_expand(final_res);
    expr_free(final_res);
    
    return ret;
}

static Expr* apply_floor_to_coeffs(Expr* e) {
    if (!e) return NULL;
    if (e->type == EXPR_INTEGER) return expr_copy(e);
    if (e->type == EXPR_REAL) return expr_new_integer((int64_t)floor(e->data.real));
    if (e->type == EXPR_FUNCTION && e->data.function.head->data.symbol == SYM_Rational) {
        int64_t n = e->data.function.args[0]->data.integer;
        int64_t d = e->data.function.args[1]->data.integer;
        int64_t q = n / d;
        int64_t r = n % d;
        if (r != 0 && ((n < 0) ^ (d < 0))) q -= 1;
        return expr_new_integer(q);
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->data.symbol == SYM_Plus) {
        Expr** args = malloc(sizeof(Expr*) * e->data.function.arg_count);
        for(size_t i=0; i<e->data.function.arg_count; i++) args[i] = apply_floor_to_coeffs(e->data.function.args[i]);
        Expr* res = eval_and_free(expr_new_function(expr_copy(e->data.function.head), args, e->data.function.arg_count));
        free(args); return res;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->data.symbol == SYM_Times) {
        Expr** args = malloc(sizeof(Expr*) * e->data.function.arg_count);
        for(size_t i=0; i<e->data.function.arg_count; i++) {
            if (i == 0) args[i] = apply_floor_to_coeffs(e->data.function.args[i]); 
            else args[i] = expr_copy(e->data.function.args[i]);
        }
        Expr* res = eval_and_free(expr_new_function(expr_copy(e->data.function.head), args, e->data.function.arg_count));
        free(args); return res;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->data.symbol == SYM_Complex) {
        Expr* re = apply_floor_to_coeffs(e->data.function.args[0]);
        Expr* im = apply_floor_to_coeffs(e->data.function.args[1]);
        Expr* res = eval_and_free(expr_new_function(expr_new_symbol("Complex"), (Expr*[]){re, im}, 2));
        return res;
    }
    return expr_copy(e);
}

static Expr* integer_poly_div(Expr* A, Expr* B, Expr** vars, size_t var_count) {
    Expr* q = exact_poly_div(A, B, vars, var_count);
    if (q) {
        Expr* f = apply_floor_to_coeffs(q);
        expr_free(q);
        if (is_zero_poly(f)) {
            expr_free(f); return NULL;
        }
        return f;
    }
    return NULL;
}

/* Reduce `poly` modulo a single divisor `m`. If `m` is an integer,    */
/* every numeric coefficient is reduced into [0, |m|). Otherwise `m`   */
/* is treated as a polynomial in its highest-indexed variable; we      */
/* repeatedly subtract scaled copies of `m` to kill any term whose     */
/* main-variable degree is at least deg(m) -- the polynomial analogue  */
/* of integer modular reduction.                                       */
static Expr* polynomial_mod_single(Expr* poly, Expr* m, bool use_integer_div) {
    if (m->type == EXPR_INTEGER) {
        int64_t m_val = m->data.integer;
        if (m_val == 0) return expr_copy(poly);
        if (m_val < 0) m_val = -m_val;
        
        Expr* expanded = expr_expand(poly);
        if (expanded->type == EXPR_FUNCTION && expanded->data.function.head->data.symbol == SYM_Plus) {
            Expr** args = malloc(sizeof(Expr*) * expanded->data.function.arg_count);
            for(size_t i=0; i<expanded->data.function.arg_count; i++) {
                Expr* term = expanded->data.function.args[i];
                if (term->type == EXPR_INTEGER) {
                    int64_t c = term->data.integer % m_val;
                    if (c < 0) c += m_val;
                    args[i] = expr_new_integer(c);
                } else if (term->type == EXPR_FUNCTION && term->data.function.head->data.symbol == SYM_Times && term->data.function.args[0]->type == EXPR_INTEGER) {
                    int64_t c = term->data.function.args[0]->data.integer % m_val;
                    if (c < 0) c += m_val;
                    Expr** t_args = malloc(sizeof(Expr*) * term->data.function.arg_count);
                    t_args[0] = expr_new_integer(c);
                    for(size_t j=1; j<term->data.function.arg_count; j++) t_args[j] = expr_copy(term->data.function.args[j]);
                    args[i] = internal_times(t_args, term->data.function.arg_count);
                    free(t_args);
                } else if (term->type == EXPR_FUNCTION && term->data.function.head->data.symbol == SYM_Complex) {
                    int64_t r = term->data.function.args[0]->data.integer % m_val;
                    if (r < 0) r += m_val;
                    int64_t i_val = term->data.function.args[1]->data.integer % m_val;
                    if (i_val < 0) i_val += m_val;
                    args[i] = eval_and_free(expr_new_function(expr_new_symbol("Complex"), (Expr*[]){expr_new_integer(r), expr_new_integer(i_val)}, 2));
                } else {
                    int64_t c = 1 % m_val;
                    if (c < 0) c += m_val;
                    if (c == 1) args[i] = expr_copy(term);
                    else args[i] = internal_times((Expr*[]){expr_new_integer(c), expr_copy(term)}, 2);
                }
            }
            Expr* res = internal_plus(args, expanded->data.function.arg_count);
            free(args);
            expr_free(expanded);
            return res;
        } else {
            Expr* term = expanded;
            Expr* res = NULL;
            if (term->type == EXPR_INTEGER) {
                int64_t c = term->data.integer % m_val;
                if (c < 0) c += m_val;
                res = expr_new_integer(c);
            } else if (term->type == EXPR_FUNCTION && term->data.function.head->data.symbol == SYM_Times && term->data.function.args[0]->type == EXPR_INTEGER) {
                int64_t c = term->data.function.args[0]->data.integer % m_val;
                if (c < 0) c += m_val;
                Expr** t_args = malloc(sizeof(Expr*) * term->data.function.arg_count);
                t_args[0] = expr_new_integer(c);
                for(size_t j=1; j<term->data.function.arg_count; j++) t_args[j] = expr_copy(term->data.function.args[j]);
                res = internal_times(t_args, term->data.function.arg_count);
                free(t_args);
            } else if (term->type == EXPR_FUNCTION && term->data.function.head->data.symbol == SYM_Complex) {
                int64_t r = term->data.function.args[0]->data.integer % m_val;
                if (r < 0) r += m_val;
                int64_t i_val = term->data.function.args[1]->data.integer % m_val;
                if (i_val < 0) i_val += m_val;
                res = eval_and_free(expr_new_function(expr_new_symbol("Complex"), (Expr*[]){expr_new_integer(r), expr_new_integer(i_val)}, 2));
            } else {
                int64_t c = 1 % m_val;
                if (c < 0) c += m_val;
                if (c == 1) res = expr_copy(term);
                else res = internal_times((Expr*[]){expr_new_integer(c), expr_copy(term)}, 2);
            }
            expr_free(expanded);
            return res;
        }
    }
    
    size_t v_count = 0, v_cap = 16;
    Expr** vars = malloc(sizeof(Expr*) * v_cap);
    
    Expr* t_list = eval_and_free(expr_new_function(expr_new_symbol("List"), (Expr*[]){expr_copy(poly), expr_copy(m)}, 2));
    collect_variables(t_list, &vars, &v_count, &v_cap);
    expr_free(t_list);
    
    if (v_count == 0) {
        free(vars);
        return expr_copy(poly);
    }
    qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);
    
    size_t m_v_count = 0, m_v_cap = 16;
    Expr** m_vars = malloc(sizeof(Expr*) * m_v_cap);
    collect_variables(m, &m_vars, &m_v_count, &m_v_cap);
    qsort(m_vars, m_v_count, sizeof(Expr*), compare_expr_ptrs);
    
    if (m_v_count == 0) {
        free(vars); free(m_vars); return expr_copy(poly);
    }
    Expr* main_var = m_vars[0];
    
    Expr** temp_vars = malloc(sizeof(Expr*) * v_count);
    size_t idx = 0;
    for (size_t i = 0; i < v_count; i++) {
        if (expr_compare(vars[i], main_var) != 0) {
            temp_vars[idx++] = vars[i];
        }
    }
    temp_vars[idx++] = main_var;
    
    Expr* exp_m = expr_expand(m);
    int d = get_degree_poly(exp_m, main_var);
    if (d == 0) {
        free(vars); free(m_vars); free(temp_vars); expr_free(exp_m); return expr_copy(poly);
    }
    Expr* lc = get_coeff(exp_m, main_var, d);
    
    Expr* curr_poly = expr_expand(poly);
    
    bool changed = true;
    while (changed) {
        changed = false;
        
        size_t t_count = 0;
        Expr** terms = NULL;
        if (curr_poly->type == EXPR_FUNCTION && curr_poly->data.function.head->data.symbol == SYM_Plus) {
            t_count = curr_poly->data.function.arg_count;
            terms = curr_poly->data.function.args;
        } else {
            t_count = 1;
            terms = &curr_poly;
        }
        
        for (size_t i = 0; i < t_count; i++) {
            Expr* term = terms[i];
            int term_deg = get_degree_poly(term, main_var);
            if (term_deg >= d) {
                Expr* term_coeff = get_coeff(term, main_var, term_deg);
                Expr* q = use_integer_div ? integer_poly_div(term_coeff, lc, temp_vars, v_count - 1) : exact_poly_div(term_coeff, lc, temp_vars, v_count - 1);
                if (q) {
                    Expr* x_pow = (term_deg - d == 0) ? expr_new_integer(1) : ((term_deg - d == 1) ? expr_copy(main_var) : internal_power((Expr*[]){expr_copy(main_var), expr_new_integer(term_deg - d)}, 2));
                    Expr* mult = internal_times((Expr*[]){q, x_pow}, 2);
                    Expr* sub = internal_times((Expr*[]){mult, expr_copy(exp_m)}, 2);
                    Expr* neg_sub = internal_times((Expr*[]){expr_new_integer(-1), sub}, 2);
                    
                    Expr* t_eval = internal_plus((Expr*[]){expr_copy(curr_poly), neg_sub}, 2);
                    Expr* next_poly = expr_expand(t_eval);
                    expr_free(t_eval);
                    expr_free(curr_poly);
                    curr_poly = next_poly;
                    changed = true;
                    expr_free(term_coeff);
                    break;
                }
                expr_free(term_coeff);
            }
        }
    }
    
    expr_free(lc);
    expr_free(exp_m);
    for(size_t i=0; i<v_count; i++) expr_free(vars[i]);
    for(size_t i=0; i<m_v_count; i++) expr_free(m_vars[i]);
    free(vars); free(m_vars); free(temp_vars);
    return curr_poly;
}

Expr* builtin_polynomialmod(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    
    Expr* expr = res->data.function.args[0];
    Expr* m = res->data.function.args[1];
    
    if (expr->type == EXPR_FUNCTION) {
        const char* head = expr->data.function.head->type == EXPR_SYMBOL ? expr->data.function.head->data.symbol : "";
        if (strcmp(head, "List") == 0 || strcmp(head, "Equal") == 0 || strcmp(head, "Unequal") == 0 ||
            strcmp(head, "Less") == 0 || strcmp(head, "LessEqual") == 0 || strcmp(head, "Greater") == 0 || 
            strcmp(head, "GreaterEqual") == 0 || strcmp(head, "And") == 0 || strcmp(head, "Or") == 0 || 
            strcmp(head, "Not") == 0) {
            Expr** args = malloc(sizeof(Expr*) * expr->data.function.arg_count);
            for (size_t i = 0; i < expr->data.function.arg_count; i++) {
                Expr* ap_args[] = {expr_copy(expr->data.function.args[i]), expr_copy(m)};
                args[i] = internal_polynomialmod(ap_args, 2);
            }
            Expr* ret = eval_and_free(expr_new_function(expr_copy(expr->data.function.head), args, expr->data.function.arg_count));
            free(args);
            return ret;
        }
    }
    
    if (m->type == EXPR_FUNCTION && m->data.function.head->data.symbol == SYM_List) {
        Expr* curr = expr_copy(expr);
        bool has_integer = false;
        for (size_t i = 0; i < m->data.function.arg_count; i++) {
            if (m->data.function.args[i]->type == EXPR_INTEGER) {
                has_integer = true;
                break;
            }
        }
        for (size_t i = 0; i < m->data.function.arg_count; i++) {
            Expr* next = polynomial_mod_single(curr, m->data.function.args[i], has_integer);
            expr_free(curr);
            curr = next;
        }
        return curr;
    }
    
    return polynomial_mod_single(expr, m, false);
}
static bool is_constant_1(Expr* e) {
    if (!e) return false;
    Expr* ev = eval_and_free(expr_copy(e));
    bool res = (ev->type == EXPR_INTEGER && ev->data.integer == 1);
    expr_free(ev);
    return res;
}

static int64_t mod_inverse_int_poly(int64_t a, int64_t m) {
    int64_t m0 = m, t, q;
    int64_t x0 = 0, x1 = 1;
    if (m == 1) return 0;
    while (a > 1) {
        if (m == 0) return 0;
        q = a / m;
        t = m;
        m = a % m, a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    if (x1 < 0) x1 += m0;
    return x1;
}

/* PolynomialExtendedGCD[a, b, x] (optionally Modulus -> p)            */
/*                                                                     */
/* Standard extended Euclidean iteration:                               */
/*   r_{i+1} = r_{i-1} - q_i * r_i                                      */
/*   s_{i+1} = s_{i-1} - q_i * s_i                                      */
/*   t_{i+1} = t_{i-1} - q_i * t_i                                      */
/* On exit, r_0 is the gcd and (s_0, t_0) are the Bezout coefficients. */
/* The result is finally normalised so the gcd is monic in `x`.        */
Expr* builtin_polynomialextendedgcd(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 3 && res->data.function.arg_count != 4)) return NULL;
    Expr* A = res->data.function.args[0];
    Expr* B = res->data.function.args[1];
    Expr* x = res->data.function.args[2];
    Expr* mod_p = NULL;
    if (res->data.function.arg_count == 4) {
        Expr* rule = res->data.function.args[3];
        if (rule->type == EXPR_FUNCTION && rule->data.function.head->data.symbol == SYM_Rule &&
            rule->data.function.args[0]->type == EXPR_SYMBOL && rule->data.function.args[0]->data.symbol == SYM_Modulus) {
            mod_p = rule->data.function.args[1];
        } else {
            return NULL;
        }
    }

    Expr* r0 = expr_expand(A);
    Expr* r1 = expr_expand(B);
    if (mod_p) {
        Expr* next0 = internal_polynomialmod((Expr*[]){r0, expr_copy(mod_p)}, 2);
        r0 = next0;
        Expr* next1 = internal_polynomialmod((Expr*[]){r1, expr_copy(mod_p)}, 2);
        r1 = next1;
    }

    Expr* s0 = expr_new_integer(1);
    Expr* t0 = expr_new_integer(0);
    Expr* s1 = expr_new_integer(0);
    Expr* t1 = expr_new_integer(1);

    while (!is_zero_poly(r1)) {
        if (get_degree_poly(r1, x) == 0) {
            expr_free(r0); r0 = r1; r1 = expr_new_integer(0);
            expr_free(s0); s0 = s1; s1 = expr_new_integer(0);
            expr_free(t0); t0 = t1; t1 = expr_new_integer(0);
            break;
        }
        Expr *q = NULL, *r2 = NULL;
        poly_div_rem(r0, r1, x, &q, &r2);
        if (!q || !r2) {
            if (q) expr_free(q);
            if (r2) expr_free(r2);
            break;
        }

        // s2 = s0 - q * s1
        Expr* q_s1 = internal_times((Expr*[]){expr_copy(q), expr_copy(s1)}, 2);
        Expr* neg_q_s1 = internal_times((Expr*[]){expr_new_integer(-1), q_s1}, 2);
        Expr* s2 = internal_plus((Expr*[]){expr_copy(s0), neg_q_s1}, 2);

        // t2 = t0 - q * t1
        Expr* q_t1 = internal_times((Expr*[]){expr_copy(q), expr_copy(t1)}, 2);
        Expr* neg_q_t1 = internal_times((Expr*[]){expr_new_integer(-1), q_t1}, 2);
        Expr* t2 = internal_plus((Expr*[]){expr_copy(t0), neg_q_t1}, 2);

        s2 = internal_expand((Expr*[]){s2}, 1);
        t2 = internal_expand((Expr*[]){t2}, 1);
        r2 = internal_expand((Expr*[]){r2}, 1);

        if (!mod_p) {
            s2 = internal_cancel((Expr*[]){internal_together((Expr*[]){s2}, 1)}, 1);
            t2 = internal_cancel((Expr*[]){internal_together((Expr*[]){t2}, 1)}, 1);
            r2 = internal_cancel((Expr*[]){internal_together((Expr*[]){r2}, 1)}, 1);
        } else {
            s2 = internal_polynomialmod((Expr*[]){s2, expr_copy(mod_p)}, 2);
            t2 = internal_polynomialmod((Expr*[]){t2, expr_copy(mod_p)}, 2);
            r2 = internal_polynomialmod((Expr*[]){r2, expr_copy(mod_p)}, 2);
        }

        expr_free(r0); r0 = r1; r1 = r2;
        expr_free(s0); s0 = s1; s1 = s2;
        expr_free(t0); t0 = t1; t1 = t2;
        expr_free(q);
    }
    expr_free(r1);
    expr_free(s1);
    expr_free(t1);

    // Normalize so that GCD is monic
    int deg = get_degree_poly(r0, x);
    if (deg >= 0 && !is_zero_poly(r0)) {
        Expr* lc = get_coeff(r0, x, deg);
        if (!is_constant_1(lc)) {
            Expr* lc_inv;
            if (mod_p) {
                if (lc->type == EXPR_INTEGER && mod_p->type == EXPR_INTEGER) {
                    lc_inv = expr_new_integer(mod_inverse_int_poly(lc->data.integer, mod_p->data.integer));
                } else {
                    lc_inv = internal_power((Expr*[]){expr_copy(lc), expr_new_integer(-1)}, 2);
                    lc_inv = internal_polynomialmod((Expr*[]){lc_inv, expr_copy(mod_p)}, 2);
                }
            } else {
                lc_inv = internal_power((Expr*[]){expr_copy(lc), expr_new_integer(-1)}, 2);
            }

            Expr* nr0 = internal_expand((Expr*[]){internal_times((Expr*[]){expr_copy(r0), expr_copy(lc_inv)}, 2)}, 1);
            Expr* ns0 = internal_expand((Expr*[]){internal_times((Expr*[]){expr_copy(s0), expr_copy(lc_inv)}, 2)}, 1);
            Expr* nt0 = internal_expand((Expr*[]){internal_times((Expr*[]){expr_copy(t0), expr_copy(lc_inv)}, 2)}, 1);
            
            if (!mod_p) {
                nr0 = internal_cancel((Expr*[]){internal_together((Expr*[]){nr0}, 1)}, 1);
                ns0 = internal_cancel((Expr*[]){internal_together((Expr*[]){ns0}, 1)}, 1);
                nt0 = internal_cancel((Expr*[]){internal_together((Expr*[]){nt0}, 1)}, 1);
            } else {
                nr0 = internal_polynomialmod((Expr*[]){nr0, expr_copy(mod_p)}, 2);
                ns0 = internal_polynomialmod((Expr*[]){ns0, expr_copy(mod_p)}, 2);
                nt0 = internal_polynomialmod((Expr*[]){nt0, expr_copy(mod_p)}, 2);
            }

            expr_free(r0); r0 = nr0;
            expr_free(s0); s0 = ns0;
            expr_free(t0); t0 = nt0;
            expr_free(lc_inv);
        }
        expr_free(lc);
    }

    Expr* coef_list = expr_new_function(expr_new_symbol("List"), (Expr*[]){s0, t0}, 2);
    Expr* ret = expr_new_function(expr_new_symbol("List"), (Expr*[]){r0, coef_list}, 2);
    return ret;
}

void poly_init(void) {
    symtab_add_builtin("PolynomialQ", builtin_polynomialq);
    symtab_get_def("PolynomialQ")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Variables", builtin_variables);
    symtab_get_def("Variables")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Coefficient", builtin_coefficient);
    symtab_get_def("Coefficient")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("CoefficientList", builtin_coefficientlist);
    symtab_get_def("CoefficientList")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("PolynomialGCD", builtin_polynomialgcd);
    symtab_get_def("PolynomialGCD")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("PolynomialLCM", builtin_polynomiallcm);
    symtab_get_def("PolynomialLCM")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("PolynomialQuotient", builtin_polynomialquotient);
    symtab_get_def("PolynomialQuotient")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("PolynomialRemainder", builtin_polynomialremainder);
    symtab_get_def("PolynomialRemainder")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("PolynomialQuotientRemainder", builtin_polynomialquotientremainder);
    symtab_get_def("PolynomialQuotientRemainder")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PolynomialQuotientRemainder",
        "PolynomialQuotientRemainder[p, q, x] returns {Quotient, Remainder}\n"
        "such that p == Quotient*q + Remainder, with deg(Remainder) < deg(q)\n"
        "in x. Single-pass companion to PolynomialQuotient/PolynomialRemainder.\n"
        "Accepts an optional Extension -> alpha rule (default None) to perform\n"
        "the division over Q(alpha)[x] rather than the rational coefficient field.");
    symtab_add_builtin("SubresultantPolynomialRemainders", builtin_subresultantpolynomialremainders);
    symtab_get_def("SubresultantPolynomialRemainders")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("SubresultantPolynomialRemainders",
        "SubresultantPolynomialRemainders[a, b, x] gives the polynomial-remainder\n"
        "chain {a, b, R_2, R_3, ...} obtained by iterating pseudo-remainder over\n"
        "K(coeffs)[x] until a constant or zero remainder is reached. Used by the\n"
        "Lazard-Rioboo-Trager rational integration pipeline; the chain is correct\n"
        "modulo content scaling, which downstream consumers strip with primitive[].");
    symtab_add_builtin("PolynomialMod", builtin_polynomialmod);
    symtab_get_def("PolynomialMod")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Collect", builtin_collect);
    symtab_get_def("Collect")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Decompose", builtin_decompose);
    symtab_get_def("Decompose")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("HornerForm", builtin_hornerform);
    symtab_get_def("HornerForm")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Resultant", builtin_resultant);
    symtab_get_def("Resultant")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("PolynomialExtendedGCD", builtin_polynomialextendedgcd);
    symtab_get_def("PolynomialExtendedGCD")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("Discriminant", builtin_discriminant);
    symtab_get_def("Discriminant")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
}
