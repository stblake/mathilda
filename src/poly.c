/* ====================================================================
 * poly.c -- Polynomial algebra primitives.
 *
 * This module implements PicoCAS's symbolic polynomial layer:
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
 * Memory convention (matches the rest of PicoCAS):
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
#include <string.h>
#include <stdlib.h>
#include <limits.h>

bool contains_any_symbol_from(Expr* expr, Expr* var);

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
                } else {
                    if (expr_eq(t_exp, trm_exp)) {
                        if (k_val == -1) k_val = 1;
                        else if (k_val != 1) return 0;
                    } else return 0;
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
bool is_zero_poly(Expr* e) {
    if (!e) return true;
    if (e->type == EXPR_INTEGER && e->data.integer == 0) return true;
    if (e->type == EXPR_REAL && e->data.real == 0.0) return true;
    
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
                        if (!is_zero_poly(clist->data.function.args[i])) {
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

/* Returns Q such that A = B * Q exactly (in K[vars]); otherwise NULL. */
/* Used by `poly_gcd_internal` to extract primitive parts and to stage */
/* GCD candidate division. The base case (var_count == 0) reduces to   */
/* exact integer / rational division -- everything else recurses on    */
/* the leading-coefficient ring.                                       */
Expr* exact_poly_div(Expr* A, Expr* B, Expr** vars, size_t var_count) {
    if (expr_eq(A, B)) return expr_new_integer(1);
    if (var_count == 0) {
        if (is_zero_poly(B)) return NULL;
        
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
            // Fall through to Times[A, B^-1] if not exact integer division
        }
        
        return internal_times((Expr*[]){expr_copy(A), internal_power((Expr*[]){expr_copy(B), expr_new_integer(-1)}, 2)}, 2);
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

Expr* builtin_polynomialquotient(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* p = res->data.function.args[0];
    Expr* q = res->data.function.args[1];
    Expr* x = res->data.function.args[2];
    
    Expr *Q, *R;
    poly_div_rem(p, q, x, &Q, &R);
    if (!Q) return NULL;
    expr_free(R);
    Expr* expanded_Q = expr_expand(Q);
    expr_free(Q);
    return expanded_Q;
}

Expr* builtin_polynomialremainder(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* p = res->data.function.args[0];
    Expr* q = res->data.function.args[1];
    Expr* x = res->data.function.args[2];
    
    Expr *Q, *R;
    poly_div_rem(p, q, x, &Q, &R);
    if (!R) return NULL;
    expr_free(Q);
    return R;
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
    
    while (!is_zero_poly(V)) {
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

/* PolynomialGCD[p1, p2, ...]                                          */
/*                                                                     */
/* Pre-processes the inputs by extracting:                              */
/*   - the integer content (numeric gcd of all literal coefficients);   */
/*   - any factors that appear with positive integer exponent in every  */
/*     argument (handles `Sqrt[x]`-style irrational generators);        */
/* Then defers to `poly_gcd_internal` for the symbolic remainder.       */
/* The final result is `numG * common_factors * symbolic_gcd`.          */
Expr* builtin_polynomialgcd(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;

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
        return expr_new_integer(1);
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
Expr* builtin_polynomiallcm(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;

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
        return expr_new_integer(1);
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
    Expr* expanded = expr_expand_patt(expr, var);
    if (!expanded) return expr_copy(expr);

    size_t term_count = 1;
    Expr** terms = &expanded;
    if (expanded->type == EXPR_FUNCTION && expanded->data.function.head->type == EXPR_SYMBOL && expanded->data.function.head->data.symbol == SYM_Plus) {
        term_count = expanded->data.function.arg_count;
        terms = expanded->data.function.args;
    }

    size_t group_cap = 16;
    size_t group_count = 0;
    CollectGroup* groups = malloc(sizeof(CollectGroup) * group_cap);

    for (size_t i = 0; i < term_count; i++) {
        Expr* term = terms[i];
        BPList term_bp;
        bp_init(&term_bp);
        decompose_to_bp(term, &term_bp);

        Expr* match_base = NULL;
        Expr* match_exp = NULL;
        int matched_idx = -1;

        for (size_t j = 0; j < term_bp.count; j++) {
            MatchEnv* env = env_new();
            if (match(term_bp.data[j].base, var, env)) {
                match_base = expr_copy(term_bp.data[j].base);
                match_exp = expr_copy(term_bp.data[j].exp);
                matched_idx = j;
                env_free(env);
                break;
            }
            env_free(env);
        }

        Expr* coeff = NULL;
        if (matched_idx == -1) {
            coeff = rebuild_from_bp(&term_bp);
        } else {
            expr_free(term_bp.data[matched_idx].base);
            expr_free(term_bp.data[matched_idx].exp);
            for (size_t j = matched_idx; j < term_bp.count - 1; j++) {
                term_bp.data[j] = term_bp.data[j + 1];
            }
            term_bp.count--;
            coeff = rebuild_from_bp(&term_bp);
        }
        bp_free(&term_bp);

        if (!match_exp) match_exp = expr_new_integer(0);

        int found = -1;
        for (size_t j = 0; j < group_count; j++) {
            bool base_eq = (!match_base && !groups[j].base) || (match_base && groups[j].base && expr_eq(match_base, groups[j].base));
            bool exp_eq = expr_eq(match_exp, groups[j].exp);
            if (base_eq && exp_eq) {
                found = j;
                break;
            }
        }

        if (found == -1) {
            if (group_count == group_cap) {
                group_cap *= 2;
                groups = realloc(groups, sizeof(CollectGroup) * group_cap);
            }
            groups[group_count].base = match_base; 
            groups[group_count].exp = match_exp;   
            groups[group_count].coeff_cap = 4;
            groups[group_count].coeff_count = 0;
            groups[group_count].coeffs = malloc(sizeof(Expr*) * 4);
            found = group_count;
            group_count++;
        } else {
            if (match_base) expr_free(match_base);
            expr_free(match_exp);
        }

        if (groups[found].coeff_count == groups[found].coeff_cap) {
            groups[found].coeff_cap *= 2;
            groups[found].coeffs = realloc(groups[found].coeffs, sizeof(Expr*) * groups[found].coeff_cap);
        }
        groups[found].coeffs[groups[found].coeff_count++] = coeff;
    }

    Expr** final_terms = malloc(sizeof(Expr*) * group_count);
    size_t final_count = 0;

    for (size_t i = 0; i < group_count; i++) {
        Expr* coeff_sum = NULL;
        if (groups[i].coeff_count == 1) {
            coeff_sum = expr_copy(groups[i].coeffs[0]);
        } else {
            Expr** c_args = malloc(sizeof(Expr*) * groups[i].coeff_count);
            for (size_t j = 0; j < groups[i].coeff_count; j++) c_args[j] = expr_copy(groups[i].coeffs[j]);
            coeff_sum = internal_plus(c_args, groups[i].coeff_count);
            free(c_args);
        }

        Expr* collected_coeff = collect_internal(coeff_sum, vars, num_vars, var_idx + 1, h);
        expr_free(coeff_sum);

        Expr* term = NULL;
        if (!groups[i].base || (groups[i].exp->type == EXPR_INTEGER && groups[i].exp->data.integer == 0)) {
            term = collected_coeff;
        } else {
            Expr* power = internal_power((Expr*[]){expr_copy(groups[i].base), expr_copy(groups[i].exp)}, 2);
            term = internal_times((Expr*[]){collected_coeff, power}, 2);
        }

        final_terms[final_count++] = term;

        if (groups[i].base) expr_free(groups[i].base);
        expr_free(groups[i].exp);
        for (size_t j = 0; j < groups[i].coeff_count; j++) {
            expr_free(groups[i].coeffs[j]);
        }
        free(groups[i].coeffs);
    }
    free(groups);

    Expr* result = NULL;
    if (final_count == 0) {
        result = expr_new_integer(0);
    } else if (final_count == 1) {
        result = final_terms[0];
    } else {
        result = internal_plus(final_terms, final_count);
    }

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

/* Resultant of P, Q in `var`. We exploit two algebraic identities for */
/* a fast path before falling back to the Sylvester matrix:            */
/*   Res(P1*P2, Q) = Res(P1,Q) * Res(P2,Q)                              */
/*   Res(P^k, Q)   = Res(P,Q)^k                                         */
/* The general case builds the (n+m)x(n+m) Sylvester matrix and        */
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
