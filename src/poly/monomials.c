/* monomials.c — MonomialList, CoefficientRules, FromCoefficientRules
 *
 * A sparse {exponent-vector -> coefficient} view of a multivariate polynomial,
 * plus the list of its monomials and the inverse reconstruction. All three
 * heads share one core:
 *
 *   1. resolve the variable list (explicit, `All`, or default Variables[poly]);
 *   2. reduce coefficients modulo an optional `Modulus -> m`;
 *   3. Expand and split into additive terms;
 *   4. decompose each term into (integer exponent vector, coefficient) w.r.t.
 *      the variables;
 *   5. merge like monomials and drop zero coefficients;
 *   6. sort by a monomial order (six named orders + explicit weight matrix).
 *
 * MonomialList and CoefficientRules differ only in how each (expvec, coeff)
 * term is rendered. FromCoefficientRules is the inverse of CoefficientRules.
 *
 * Every named order is a special case of "descending lexicographic order of the
 * weighted exponent vectors w.v" (Wolfram's own model), so a single weight
 * matrix drives one comparator. For k variables (e_i = i-th unit row,
 * deg = all-ones row), greatest monomial first:
 *
 *   Lexicographic                    e_0, e_1, ..., e_{k-1}          (default)
 *   NegativeLexicographic           -e_0, ..., -e_{k-1}   (= Sort, ascending)
 *   DegreeLexicographic              deg, e_0, ..., e_{k-2}
 *   DegreeReverseLexicographic       deg, -e_{k-1}, ..., -e_1
 *   NegativeDegreeLexicographic     -deg, e_0, ..., e_{k-2}
 *   NegativeDegreeReverseLexicographic  -deg, -e_{k-1}, ..., -e_1
 *
 * These reproduce Wolfram's explicit-matrix spellings exactly (e.g.
 * DegreeLexicographic on {x,y} is {{1,1},{1,0}}; DegreeReverseLexicographic on
 * {x,y,z} is {{1,1,1},{0,0,-1},{0,-1,0}}).
 */

#include "poly.h"
#include "expr.h"
#include "sym_names.h"
#include "internal.h"
#include "expand.h"
#include "symtab.h"
#include "sparse/sparse.h"   /* CoefficientArrays builds SparseArrays */
#include "attr.h"
#include "groebner.h"   /* gb_build_order_matrix, gb_classify_named_order */
#include "flint_bridge.h" /* flint_field_monomials — field-coefficient fast path */
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "pack.h"      /* pack_unpack: materialise a visible NDArray argument */

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>

/* ------------------------------------------------------------------ */
/*  Small structural predicates                                        */
/* ------------------------------------------------------------------ */

static int head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == sym;
}
static int is_list_head(const Expr* e)  { return head_is(e, SYM_List); }
static int is_rule_head(const Expr* e)  {
    return (head_is(e, SYM_Rule) || head_is(e, SYM_RuleDelayed)) &&
           e->data.function.arg_count == 2;
}

/* ------------------------------------------------------------------ */
/*  One monomial: exponent vector + coefficient + weighted sort key    */
/* ------------------------------------------------------------------ */

typedef struct {
    int*     exps;   /* k non-negative integers (owned) */
    Expr*    coeff;  /* owned; NULL once moved into an output node */
    int64_t* key;    /* keylen weighted sort values (owned; filled before sort) */
} Mono;

/* qsort is context-free in C99, so the comparator reads the shared key length
 * set immediately before each sort. The evaluator is single-threaded here. */
static int g_keylen;

static int cmp_key_desc(const void* pa, const void* pb) {
    const Mono* a = (const Mono*)pa;
    const Mono* b = (const Mono*)pb;
    for (int i = 0; i < g_keylen; i++) {
        if (a->key[i] != b->key[i]) return a->key[i] > b->key[i] ? -1 : 1;
    }
    return 0;
}

static void monos_free(Mono* monos, size_t n) {
    if (!monos) return;
    for (size_t i = 0; i < n; i++) {
        free(monos[i].exps);
        free(monos[i].key);
        if (monos[i].coeff) expr_free(monos[i].coeff);
    }
    free(monos);
}

/* ------------------------------------------------------------------ */
/*  Weight matrix for a monomial order                                 */
/* ------------------------------------------------------------------ */

/* Fill an all-zero rows*k matrix, then set the requested pattern. `deg_sign`
 * is +1 for a degree-first order, -1 for its Negative variant, 0 for a pure
 * (non-degree) order. `rev` selects reverse-lex tail rows when true. */
static int64_t* named_matrix(int k, int deg_sign, int rev, int neg_lex,
                             int* rows_out) {
    /* Single source of truth lives in groebner.c so GroebnerBasis /
     * PolynomialReduce and this family build byte-identical matrices. */
    return gb_build_order_matrix(k, deg_sign, rev, neg_lex, rows_out);
}

/* Parse an explicit {{...},{...}} weight matrix into a fresh rows*k row-major
 * int64 array (caller frees). Every row must have length k and integer
 * entries. Returns NULL on any shape mismatch or non-integer entry. */
static int64_t* explicit_matrix(const Expr* v, int k, int* rows_out) {
    if (!is_list_head(v)) return NULL;
    size_t r = v->data.function.arg_count;
    if (r == 0) return NULL;
    int64_t* w = malloc(sizeof(int64_t) * r * (size_t)(k > 0 ? k : 1));
    for (size_t i = 0; i < r; i++) {
        const Expr* row = v->data.function.args[i];
        if (!is_list_head(row) || (int)row->data.function.arg_count != k) {
            free(w); return NULL;
        }
        for (int j = 0; j < k; j++) {
            const Expr* e = row->data.function.args[j];
            if (e->type == EXPR_INTEGER) {
                w[i * (size_t)k + j] = e->data.integer;
            } else {
                free(w); return NULL;
            }
        }
    }
    *rows_out = (int)r;
    return w;
}

/* Resolve the order argument (may be NULL -> default "Lexicographic") into a
 * weight matrix. Returns NULL on an unknown order string / malformed matrix. */
static int64_t* order_weight_matrix(const Expr* order_arg, int k, int* rows_out) {
    if (!order_arg) return named_matrix(k, 0, 0, 0, rows_out);   /* Lexicographic */
    /* A named order, as a string or a symbol:
     * "DegreeReverseLexicographic" or DegreeReverseLexicographic. */
    if (order_arg->type == EXPR_STRING || order_arg->type == EXPR_SYMBOL) {
        const char* name = order_arg->type == EXPR_STRING ? order_arg->data.string
                                                          : order_arg->data.symbol.name;
        int ds = 0, rv = 0, nl = 0; GBOrder hint;
        if (gb_classify_named_order(name, &ds, &rv, &nl, &hint))
            return named_matrix(k, ds, rv, nl, rows_out);
        return NULL;                                /* unknown order name */
    }
    if (is_list_head(order_arg))
        return explicit_matrix(order_arg, k, rows_out);
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Variable resolution                                                */
/* ------------------------------------------------------------------ */

/* Return a fresh array of k BORROWED variable pointers and set *k. If
 * vars_arg is NULL, computes Variables[poly] (stored in *owned, caller frees
 * with expr_free). A List borrows its elements; a bare expression is one
 * variable. Returns NULL on allocation trouble only. */
static Expr** resolve_vars(Expr* poly, Expr* vars_arg, int* k, Expr** owned) {
    *owned = NULL;
    Expr* list = NULL;         /* the List whose args are the variables */
    if (!vars_arg) {
        *owned = internal_variables((Expr*[]){ expr_copy(poly) }, 1);
        list = *owned;
    } else if (is_list_head(vars_arg)) {
        list = vars_arg;
    } else {
        /* single variable */
        Expr** vv = malloc(sizeof(Expr*));
        vv[0] = vars_arg;
        *k = 1;
        return vv;
    }
    int n = (list && is_list_head(list)) ? (int)list->data.function.arg_count : 0;
    Expr** vv = malloc(sizeof(Expr*) * (size_t)(n > 0 ? n : 1));
    for (int i = 0; i < n; i++) vv[i] = list->data.function.args[i];
    *k = n;
    return vv;
}

/* ------------------------------------------------------------------ */
/*  Term decomposition                                                 */
/* ------------------------------------------------------------------ */

/* Accumulate the coefficient factors of one term (everything not a power of a
 * variable) into a growable array of fresh copies, and add variable exponents
 * into exps[]. Recurses into nested Times. Returns 0 if the term contains a
 * variable raised to a negative or non-integer power (not polynomial). */
static int walk_term(Expr* f, Expr** vars, int k, int* exps,
                     Expr*** kept, size_t* nk, size_t* cap) {
    if (head_is(f, SYM_Times)) {
        for (size_t i = 0; i < f->data.function.arg_count; i++)
            if (!walk_term(f->data.function.args[i], vars, k, exps, kept, nk, cap))
                return 0;
        return 1;
    }
    /* factor == variable */
    for (int i = 0; i < k; i++) {
        if (expr_eq(f, vars[i])) { exps[i] += 1; return 1; }
    }
    /* factor == Power[variable, non-negative integer] */
    if (head_is(f, SYM_Power) && f->data.function.arg_count == 2) {
        Expr* base = f->data.function.args[0];
        Expr* e    = f->data.function.args[1];
        for (int i = 0; i < k; i++) {
            if (expr_eq(base, vars[i])) {
                if (e->type == EXPR_INTEGER && e->data.integer >= 0 &&
                    e->data.integer <= INT_MAX) {
                    exps[i] += (int)e->data.integer;
                    return 1;
                }
                return 0;   /* variable to a bad power: not polynomial */
            }
        }
    }
    /* everything else is part of the coefficient */
    if (*nk == *cap) { *cap *= 2; *kept = realloc(*kept, sizeof(Expr*) * *cap); }
    (*kept)[(*nk)++] = expr_copy(f);
    return 1;
}

/* Decompose one expanded term into exps[k] and a fresh coefficient. Returns 0
 * if not polynomial in the variables. */
static int term_to_expvec(Expr* term, Expr** vars, int k, int* exps,
                          Expr** out_coeff) {
    for (int i = 0; i < k; i++) exps[i] = 0;
    size_t cap = 4, nk = 0;
    Expr** kept = malloc(sizeof(Expr*) * cap);
    if (!walk_term(term, vars, k, exps, &kept, &nk, &cap)) {
        for (size_t i = 0; i < nk; i++) expr_free(kept[i]);
        free(kept);
        return 0;
    }
    Expr* coeff;
    if (nk == 0) coeff = expr_new_integer(1);
    else         coeff = internal_times(kept, nk);   /* consumes the copies */
    free(kept);
    *out_coeff = coeff;
    return 1;
}

/* ------------------------------------------------------------------ */
/*  Build the sorted, merged monomial list                             */
/* ------------------------------------------------------------------ */

/* On success returns a malloc'd array of `*n_out` monomials (possibly zero)
 * sorted by the weight matrix; on non-polynomial input returns NULL. `poly`,
 * `vars`, `W`, and `modulus` are borrowed. */
/* Sort by the weight key, merge byte-equal exponent vectors (adjacent after the
 * sort) summing their coefficients, and drop zero coefficients.  Shared by the
 * generic term loop and the field fast path so both finish identically.  Sets
 * *n_out and returns the (in-place compacted) `monos`. */
static Mono* finish_monomials(Mono* monos, size_t n, int keylen, int k,
                              size_t* n_out) {
    g_keylen = keylen;
    qsort(monos, n, sizeof(Mono), cmp_key_desc);

    size_t w = 0;
    for (size_t r = 0; r < n; ) {
        size_t s = r + 1;
        Expr* coeff = monos[r].coeff;
        monos[r].coeff = NULL;
        while (s < n &&
               memcmp(monos[r].exps, monos[s].exps, sizeof(int) * (size_t)k) == 0) {
            coeff = internal_plus((Expr*[]){ coeff, monos[s].coeff }, 2);
            monos[s].coeff = NULL;
            s++;
        }
        int is_zero = (coeff->type == EXPR_INTEGER && coeff->data.integer == 0);
        if (is_zero) {
            expr_free(coeff);
            free(monos[r].exps);
            free(monos[r].key);
        } else {
            monos[w].exps  = monos[r].exps;
            monos[w].key   = monos[r].key;
            monos[w].coeff = coeff;
            w++;
        }
        /* free the merged-away representatives (r+1 .. s-1) */
        for (size_t x = r + 1; x < s; x++) {
            free(monos[x].exps);
            free(monos[x].key);
        }
        r = s;
    }

    *n_out = w;
    return monos;
}

/* Native field-coefficient fast path for build_monomials.  When `poly` is a
 * polynomial in `vars` with AlgebraicNumber[theta,..] (one number field Q(theta))
 * coefficients, read the monomials straight out of the reduced FLINT polynomial
 * — skipping Expand + per-term internal_times/internal_plus over AlgebraicNumber
 * scalars — then finish through the shared sort/merge so the result is identical
 * to the generic path.  Returns the Mono array (via finish_monomials) on success,
 * or NULL to decline (caller uses the generic path). */
static Mono* build_monomials_field(Expr* poly, Expr** vars, int k,
                                   const int64_t* W, int rows, size_t* n_out) {
    int*   fexps  = NULL;
    Expr** fcoeff = NULL;
    size_t fn     = 0;
    int rc = flint_field_monomials(poly, vars, k, &fexps, &fcoeff, &fn);
    if (rc < 0) return NULL;

    int keylen = rows + k;
    Mono* monos = malloc(sizeof(Mono) * (fn > 0 ? fn : 1));
    for (size_t t = 0; t < fn; t++) {
        int* exps = malloc(sizeof(int) * (size_t)(k > 0 ? k : 1));
        for (int j = 0; j < k; j++) exps[j] = fexps[t * (size_t)k + j];
        int64_t* key = malloc(sizeof(int64_t) * (size_t)(keylen > 0 ? keylen : 1));
        for (int r = 0; r < rows; r++) {
            int64_t acc = 0;
            for (int j = 0; j < k; j++) acc += W[(size_t)r * k + j] * exps[j];
            key[r] = acc;
        }
        for (int j = 0; j < k; j++) key[rows + j] = exps[j];
        monos[t].exps  = exps;
        monos[t].coeff = fcoeff[t];      /* move ownership */
        monos[t].key   = key;
    }
    free(fexps);
    free(fcoeff);
    return finish_monomials(monos, fn, keylen, k, n_out);
}

/* Plain-Q counterpart of build_monomials_field: reads the (expvec, rational
 * coefficient) terms of a Q-polynomial straight out of fmpq_mpoly (via
 * flint_polynomial_monomials), skipping the generic Expr Expand + per-term walk
 * + Orderless re-sort.  Same finish_monomials tail, so the result is identical
 * to the generic path. */
static Mono* build_monomials_polyQ(Expr* poly, Expr** vars, int k,
                                   const int64_t* W, int rows, size_t* n_out) {
    int*   fexps  = NULL;
    Expr** fcoeff = NULL;
    size_t fn     = 0;
    int rc = flint_polynomial_monomials(poly, vars, k, &fexps, &fcoeff, &fn);
    if (rc < 0) return NULL;

    int keylen = rows + k;
    Mono* monos = malloc(sizeof(Mono) * (fn > 0 ? fn : 1));
    for (size_t t = 0; t < fn; t++) {
        int* exps = malloc(sizeof(int) * (size_t)(k > 0 ? k : 1));
        for (int j = 0; j < k; j++) exps[j] = fexps[t * (size_t)k + j];
        int64_t* key = malloc(sizeof(int64_t) * (size_t)(keylen > 0 ? keylen : 1));
        for (int r = 0; r < rows; r++) {
            int64_t acc = 0;
            for (int j = 0; j < k; j++) acc += W[(size_t)r * k + j] * exps[j];
            key[r] = acc;
        }
        for (int j = 0; j < k; j++) key[rows + j] = exps[j];
        monos[t].exps  = exps;
        monos[t].coeff = fcoeff[t];      /* move ownership */
        monos[t].key   = key;
    }
    free(fexps);
    free(fcoeff);
    return finish_monomials(monos, fn, keylen, k, n_out);
}

static Mono* build_monomials(Expr* poly, Expr** vars, int k,
                             const int64_t* W, int rows, Expr* modulus,
                             size_t* n_out) {
    *n_out = 0;

    /* Native read-off (no modulus): straight from FLINT, skipping the generic
     * Expand + per-term evaluator work + Orderless re-sort.  Try the field path
     * (algebraic coefficients) first, then the plain-Q path; each declines
     * (NULL) unless `poly` is a polynomial in exactly `vars`. */
    if (!modulus && flint_bridge_available()) {
        Mono* fm = build_monomials_field(poly, vars, k, W, rows, n_out);
        if (fm) return fm;
        Mono* qm = build_monomials_polyQ(poly, vars, k, W, rows, n_out);
        if (qm) return qm;
    }

    /* Modulus first (reduces integer coefficients into [0,m) and drops the
     * terms that vanish), then Expand. */
    Expr* reduced = NULL;
    if (modulus) {
        reduced = internal_polynomialmod(
            (Expr*[]){ expr_copy(poly), expr_copy(modulus) }, 2);
    }
    Expr* expanded = expr_expand(reduced ? reduced : poly);
    if (reduced) expr_free(reduced);
    if (!expanded) return NULL;

    /* Additive terms. */
    Expr** terms;
    size_t nterms;
    if (head_is(expanded, SYM_Plus)) {
        terms  = expanded->data.function.args;
        nterms = expanded->data.function.arg_count;
    } else {
        terms  = &expanded;
        nterms = 1;
    }

    Mono* monos = malloc(sizeof(Mono) * (nterms > 0 ? nterms : 1));
    size_t n = 0;
    int keylen = rows + k;

    for (size_t t = 0; t < nterms; t++) {
        int* exps = malloc(sizeof(int) * (size_t)(k > 0 ? k : 1));
        Expr* coeff;
        if (!term_to_expvec(terms[t], vars, k, exps, &coeff)) {
            free(exps);
            monos_free(monos, n);
            expr_free(expanded);
            return NULL;             /* not polynomial in the variables */
        }
        /* weighted sort key: rows weighted values, then the raw exponents as a
         * deterministic tie-break for a degenerate user matrix. */
        int64_t* key = malloc(sizeof(int64_t) * (size_t)(keylen > 0 ? keylen : 1));
        for (int r = 0; r < rows; r++) {
            int64_t acc = 0;
            for (int j = 0; j < k; j++) acc += W[(size_t)r * k + j] * exps[j];
            key[r] = acc;
        }
        for (int j = 0; j < k; j++) key[rows + j] = exps[j];
        monos[n].exps  = exps;
        monos[n].coeff = coeff;
        monos[n].key   = key;
        n++;
    }
    expr_free(expanded);

    return finish_monomials(monos, n, keylen, k, n_out);
}

/* ------------------------------------------------------------------ */
/*  Shared front end for MonomialList / CoefficientRules               */
/* ------------------------------------------------------------------ */

/* Parse (poly, [vars], [order]) + trailing Modulus->m, resolve variables and
 * the order matrix, and build the sorted monomials. Returns the Mono array (or
 * NULL to decline) and hands the caller the borrowed `vars`/`owned`/`k` it
 * needs to render output. */
static Mono* parse_and_build(Expr* res, Expr*** vars_out, Expr** owned_out,
                             int* k_out, size_t* n_out) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc  = res->data.function.arg_count;
    Expr** args  = res->data.function.args;

    /* strip trailing option rules; capture Modulus */
    Expr* modulus = NULL;
    size_t npos = argc;
    while (npos > 0 && is_rule_head(args[npos - 1])) {
        Expr* rule = args[npos - 1];
        Expr* key  = rule->data.function.args[0];
        if (key->type == EXPR_SYMBOL && key->data.symbol.name == SYM_Modulus) {
            Expr* m = rule->data.function.args[1];
            if (m->type == EXPR_INTEGER || m->type == EXPR_BIGINT) modulus = m;
        }
        npos--;
    }
    if (npos < 1 || npos > 3) return NULL;

    Expr* poly = args[0];
    /* Symbolic/structural head only: a bare List or NDArray is not a polynomial. */
    if (poly->type == EXPR_NDARRAY || is_list_head(poly)) return NULL;

    Expr* vars_arg  = (npos >= 2) ? args[1] : NULL;
    Expr* order_arg = (npos >= 3) ? args[2] : NULL;
    /* `All` (or omitted) means Variables[poly]. */
    if (vars_arg && vars_arg->type == EXPR_SYMBOL &&
        vars_arg->data.symbol.name == SYM_All) {
        vars_arg = NULL;
    }

    int k;
    Expr* owned = NULL;
    Expr** vars = resolve_vars(poly, vars_arg, &k, &owned);

    int rows;
    int64_t* W = order_weight_matrix(order_arg, k, &rows);
    if (!W) { free(vars); if (owned) expr_free(owned); return NULL; }

    size_t n;
    Mono* monos = build_monomials(poly, vars, k, W, rows, modulus, &n);
    free(W);
    if (!monos) { free(vars); if (owned) expr_free(owned); return NULL; }

    *vars_out  = vars;
    *owned_out = owned;
    *k_out     = k;
    *n_out     = n;
    return monos;
}

/* ------------------------------------------------------------------ */
/*  CoefficientRules                                                   */
/* ------------------------------------------------------------------ */

Expr* builtin_coefficientrules(Expr* res) {
    Expr** vars; Expr* owned; int k; size_t n;
    Mono* monos = parse_and_build(res, &vars, &owned, &k, &n);
    if (!monos) return NULL;

    Expr** out = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr** es = malloc(sizeof(Expr*) * (size_t)(k > 0 ? k : 1));
        for (int j = 0; j < k; j++) es[j] = expr_new_integer(monos[i].exps[j]);
        Expr* explist = expr_new_function(expr_new_symbol(SYM_List), es, (size_t)k);
        free(es);
        Expr* coeff = monos[i].coeff; monos[i].coeff = NULL;    /* move */
        out[i] = expr_new_function(expr_new_symbol(SYM_Rule),
                                   (Expr*[]){ explist, coeff }, 2);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, n);
    free(out);

    monos_free(monos, n);
    free(vars);
    if (owned) expr_free(owned);
    return result;
}

/* ------------------------------------------------------------------ */
/*  MonomialList                                                       */
/* ------------------------------------------------------------------ */

Expr* builtin_monomiallist(Expr* res) {
    Expr** vars; Expr* owned; int k; size_t n;
    Mono* monos = parse_and_build(res, &vars, &owned, &k, &n);
    if (!monos) return NULL;

    Expr** out = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr** fac = malloc(sizeof(Expr*) * (size_t)(k + 1));
        size_t nf = 0;
        fac[nf++] = monos[i].coeff; monos[i].coeff = NULL;      /* move */
        for (int j = 0; j < k; j++) {
            int e = monos[i].exps[j];
            if (e == 1) {
                fac[nf++] = expr_copy(vars[j]);
            } else if (e > 1) {
                fac[nf++] = internal_power(
                    (Expr*[]){ expr_copy(vars[j]), expr_new_integer(e) }, 2);
            }
        }
        out[i] = internal_times(fac, nf);   /* Times drops a leading 1 */
        free(fac);
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, n);
    free(out);

    monos_free(monos, n);
    free(vars);
    if (owned) expr_free(owned);
    return result;
}

/* ------------------------------------------------------------------ */
/*  CoefficientArrays                                                  */
/* ------------------------------------------------------------------ */
/* CoefficientArrays[polys, vars] decomposes a polynomial system by TOTAL DEGREE
 * in vars: element d + 1 of the result carries every degree-d coefficient, so
 * that for a system the first two elements are exactly the constant vector and
 * the coefficient matrix of a linear system --
 *
 *     CoefficientArrays[{x + 2 y - 3, 3 x - y}, {x, y}]
 *       -> {{-3, 0}, {{1, 2}, {3, -1}}}
 *
 * which is what a linear solver wants: `{b, M} = CoefficientArrays[eqs, vars]`.
 *
 * Index convention (Mathematica's). A degree-d coefficient lives at the index
 * tuple naming its variables in NON-DECREASING order, and every other
 * permutation of that tuple is 0 -- the array is upper-triangular, not
 * symmetric. So `CoefficientArrays[x y, {x, y}]` puts the 1 at {1, 2} and
 * leaves {2, 1} zero, and `x^2` goes to {1, 1}. A list of polynomials adds the
 * equation index as the LEADING axis, so element d + 1 has rank d + 1; a single
 * polynomial has rank d (element 1 is then a bare scalar, not a vector).
 *
 * Result type. Every array element is a SparseArray with
 * background 0 -- built straight from the monomials through a COO builder, so
 * storage is proportional to the number of terms, not to k^d. The one scalar is
 * the constant term of a SINGLE polynomial (rank 0). A system's constant vector
 * is a rank-1 SparseArray. Callers test `Head[ca[[2]]] === SparseArray` (a
 * symbolic-algebra package's collector aborts otherwise), and `Normal` recovers
 * the dense form.
 *
 * Rank beyond SP_MAX_RANK (a degree-33 element) has no SparseArray, and such an
 * element is cheap only when k is tiny, so it alone falls back to a dense
 * nested List. The density guard applies to those dense elements only: a
 * total over CA_MAX_ENTRIES is declined with ::dense rather than silently
 * exhausting memory. */

#define CA_MAX_ENTRIES 20000000L

/* Nested dense zero array of rank `rank` with every axis `k` long. */
static Expr* ca_zeros(int k, int rank) {
    if (rank == 0) return expr_new_integer(0);
    Expr** args = malloc(sizeof(Expr*) * (size_t)(k > 0 ? k : 1));
    for (int i = 0; i < k; i++) args[i] = ca_zeros(k, rank - 1);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), args, (size_t)k);
    free(args);
    return out;
}

/* Overwrite the entry of the nested dense array `arr` at idx[0..rank) (0-based)
 * with `val` (adopted). */
static void ca_set(Expr* arr, const int* idx, int rank, Expr* val) {
    for (int d = 0; d < rank; d++) {
        if (!is_list_head(arr) || (size_t)idx[d] >= arr->data.function.arg_count) {
            expr_free(val);
            return;
        }
        if (d == rank - 1) {
            expr_free(arr->data.function.args[idx[d]]);
            arr->data.function.args[idx[d]] = val;
            return;
        }
        arr = arr->data.function.args[idx[d]];
    }
    /* rank == 0: nothing addressable, caller handles the scalar case. */
    expr_free(val);
}

/* The constant term of `e` when no variable is named: `e` itself, reduced by
 * the modulus and expanded exactly as build_monomials treats a polynomial. */
static Expr* ca_constant(Expr* e, Expr* modulus) {
    Expr* reduced = modulus ? internal_polynomialmod(
        (Expr*[]){ expr_copy(e), expr_copy(modulus) }, 2) : NULL;
    Expr* out = expr_expand(reduced ? reduced : e);
    if (reduced) expr_free(reduced);
    return out ? out : expr_copy(e);
}

Expr* builtin_coefficientarrays(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    /* Strip trailing option rules; capture Modulus, as CoefficientRules does. */
    Expr* modulus = NULL;
    size_t npos = argc;
    while (npos > 0 && is_rule_head(args[npos - 1])) {
        Expr* rule = args[npos - 1];
        Expr* key  = rule->data.function.args[0];
        if (key->type == EXPR_SYMBOL && key->data.symbol.name == SYM_Modulus) {
            Expr* m = rule->data.function.args[1];
            if (m->type == EXPR_INTEGER || m->type == EXPR_BIGINT) modulus = m;
        }
        npos--;
    }
    if (npos < 1 || npos > 2) return NULL;

    Expr* polys = args[0];

    /* A VISIBLE NDArray is a system of constant polynomials, so materialise it
     * and go on. The transparency gate never touches a visible NDArray (only a
     * packed List), and this head is not on AWARE, so without this the call
     * would fall through to "one polynomial" -- and an NDArray is not a
     * polynomial, so the answer would be a silently unevaluated
     * CoefficientArrays[NDArray[...], vars] rather than {{1, 2, 3}}. */
    Expr* unpacked = NULL;
    if (polys->type == EXPR_NDARRAY) {
        unpacked = pack_unpack(polys);
        if (!unpacked) return NULL;
        polys = unpacked;
    }

    /* A List argument is a SYSTEM (equation index becomes the leading axis);
     * anything else is one polynomial and the arrays carry no leading axis. */
    int      is_system = is_list_head(polys);
    size_t   m         = is_system ? polys->data.function.arg_count : 1;
    Expr**   eqs       = malloc(sizeof(Expr*) * (m > 0 ? m : 1));
    for (size_t i = 0; i < m; i++)
        eqs[i] = is_system ? polys->data.function.args[i] : polys;

    /* Variables. `All` or omitted means Variables[] of the whole input -- taken
     * over the List as a unit so every equation shares one column order. */
    Expr* vars_arg = (npos >= 2) ? args[1] : NULL;
    if (vars_arg && vars_arg->type == EXPR_SYMBOL &&
        vars_arg->data.symbol.name == SYM_All) vars_arg = NULL;

    int k;
    Expr* owned = NULL;
    Expr** vars = resolve_vars(polys, vars_arg, &k, &owned);

    /* No variables (CoefficientArrays[p, {}], or a constant input): every term
     * is the constant term, so the answer is the single degree-0 element --
     * the expanded scalar {p} for one polynomial, the constant vector as a
     * rank-1 SparseArray for a system (the shape the k > 0 path gives its
     * first element). Handled here because the monomial builders and the
     * order matrices are written for k >= 1. */
    if (k == 0) {
        Expr* el = NULL;
        if (!is_system) {
            el = ca_constant(polys, modulus);
        } else {
            int64_t dim = (int64_t)m;
            SpCoo c0;
            sp_coo_init(&c0, 1, &dim, SPV_EXPR);
            bool pushed = true;
            for (size_t i = 0; i < m && pushed; i++) {
                int64_t at = (int64_t)i;
                pushed = sp_coo_push_expr(&c0, &at,
                             ca_constant(eqs[i], modulus));
            }
            SparseArray* sa = NULL;
            if (pushed) sa = sp_coo_finish(&c0, expr_new_integer(0), true);
            else        sp_coo_free(&c0);
            if (sa) el = sp_to_expr(sa);
        }
        free(eqs); free(vars);
        if (owned) expr_free(owned);
        if (unpacked) expr_free(unpacked);
        return el ? expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ el }, 1)
                  : NULL;
    }
    if (k < 0) {
        free(eqs); free(vars);
        if (owned) expr_free(owned);
        if (unpacked) expr_free(unpacked);
        return NULL;
    }

    /* The default (lexicographic) monomial order: the sort order is irrelevant
     * here, since every monomial is placed by its own exponent vector, but
     * build_monomials needs a weight matrix. */
    int rows = 0;
    int64_t* W = order_weight_matrix(NULL, k, &rows);
    if (!W) {
        free(eqs); free(vars);
        if (owned) expr_free(owned);
        if (unpacked) expr_free(unpacked);
        return NULL;
    }

    /* Decompose every equation, and find the maximum total degree. */
    Mono** all   = calloc(m ? m : 1, sizeof(Mono*));
    size_t* counts = calloc(m ? m : 1, sizeof(size_t));
    int maxdeg = 0;
    int ok = 1;
    for (size_t i = 0; i < m && ok; i++) {
        all[i] = build_monomials(eqs[i], vars, k, W, rows, modulus, &counts[i]);
        if (!all[i]) { ok = 0; break; }
        for (size_t t = 0; t < counts[i]; t++) {
            int d = 0;
            for (int j = 0; j < k; j++) d += all[i][t].exps[j];
            if (d > maxdeg) maxdeg = d;
        }
    }
    free(W);

    /* Element d + 1 has rank d (+ 1 for a system). Up to SP_MAX_RANK it is a
     * SparseArray; past that it is dense, holding m * k^d entries -- guard
     * the dense total only. */
    int sys_off = is_system ? 1 : 0;
    if (ok) {
        double total = 0.0, kd = 1.0;
        for (int d = 0; d <= maxdeg; d++) {
            if (d + sys_off > SP_MAX_RANK) total += (double)m * kd;
            if (total > (double)CA_MAX_ENTRIES) break;
            kd *= (double)k;
        }
        if (total > (double)CA_MAX_ENTRIES) {
            mth_message("CoefficientArrays", "dense",
                "the dense result for degree %d in %d variables exceeds %ld entries; "
                "arrays of rank above %d are dense and this size is declined.",
                maxdeg, k, CA_MAX_ENTRIES, SP_MAX_RANK);
            ok = 0;
        }
    }

    if (!ok) {
        for (size_t i = 0; i < m; i++) if (all[i]) monos_free(all[i], counts[i]);
        free(all); free(counts); free(eqs); free(vars);
        if (owned) expr_free(owned);
        if (unpacked) expr_free(unpacked);
        return NULL;
    }

    /* One builder per sparse element (degree d with rank d + sys_off in
     * 1..SP_MAX_RANK); a dense nested List for the single polynomial's scalar
     * slot (d == 0, initially 0) and for any element past SP_MAX_RANK. */
    size_t nel  = (size_t)(maxdeg >= 0 ? maxdeg : 0) + 1;
    Expr** out  = calloc(nel, sizeof(Expr*));
    SpCoo* coo  = calloc(nel, sizeof(SpCoo));
    bool*  spd  = calloc(nel, sizeof(bool));
    int64_t dims[SP_MAX_RANK];
    for (int d = 0; d <= maxdeg; d++) {
        int rank = d + sys_off;
        if (rank >= 1 && rank <= SP_MAX_RANK) {
            if (is_system) dims[0] = (int64_t)m;
            for (int r = sys_off; r < rank; r++) dims[r] = (int64_t)k;
            sp_coo_init(&coo[d], rank, dims, SPV_EXPR);
            spd[d] = true;
        } else if (is_system) {
            Expr** rowsv = malloc(sizeof(Expr*) * (m > 0 ? m : 1));
            for (size_t i = 0; i < m; i++) rowsv[i] = ca_zeros(k, d);
            out[d] = expr_new_function(expr_new_symbol(SYM_List), rowsv, m);
            free(rowsv);
        } else {
            out[d] = ca_zeros(k, d);
        }
    }

    /* Place each monomial at (equation,) its sorted variable-index tuple. */
    int*     idx = malloc(sizeof(int) * (size_t)(maxdeg > 0 ? maxdeg : 1));
    int64_t* co  = malloc(sizeof(int64_t) * nel);
    for (size_t i = 0; i < m && ok; i++) {
        for (size_t t = 0; t < counts[i]; t++) {
            int d = 0;
            for (int j = 0; j < k; j++)
                for (int r = 0; r < all[i][t].exps[j]; r++) idx[d++] = j;
            Expr* coeff = all[i][t].coeff; all[i][t].coeff = NULL;   /* move */
            if (spd[d]) {
                if (is_system) co[0] = (int64_t)i;
                for (int r = 0; r < d; r++) co[r + sys_off] = idx[r];
                if (!sp_coo_push_expr(&coo[d], co, coeff)) { ok = 0; break; }
                continue;
            }
            if (!is_system && d == 0) {        /* the scalar constant term */
                expr_free(out[0]);
                out[0] = coeff;
                continue;
            }
            ca_set(is_system ? out[d]->data.function.args[i] : out[d], idx, d, coeff);
        }
    }
    free(idx); free(co);

    for (int d = 0; d <= maxdeg; d++) {
        if (!spd[d]) continue;
        if (!ok) { sp_coo_free(&coo[d]); continue; }
        SparseArray* sa = sp_coo_finish(&coo[d], expr_new_integer(0), true);
        if (sa) out[d] = sp_to_expr(sa);
        else    ok = 0;
    }
    free(coo); free(spd);

    Expr* result = NULL;
    if (ok) {
        result = expr_new_function(expr_new_symbol(SYM_List), out, nel);
    } else {
        for (int d = 0; d <= maxdeg; d++) if (out[d]) expr_free(out[d]);
    }
    free(out);

    for (size_t i = 0; i < m; i++) monos_free(all[i], counts[i]);
    free(all); free(counts); free(eqs); free(vars);
    if (owned) expr_free(owned);
    if (unpacked) expr_free(unpacked);
    return result;
}

/* ------------------------------------------------------------------ */
/*  FromCoefficientRules                                               */
/* ------------------------------------------------------------------ */

Expr* builtin_fromcoefficientrules(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) {
        mth_message("FromCoefficientRules", "argrx", "FromCoefficientRules called with %zu arguments; 2 arguments are expected.", argc);
        return NULL;
    }
    Expr* rules    = res->data.function.args[0];
    Expr* vars_arg = res->data.function.args[1];
    if (!is_list_head(rules)) return NULL;

    /* variables: a List or a single variable */
    int k;
    Expr** vars;
    if (is_list_head(vars_arg)) {
        k = (int)vars_arg->data.function.arg_count;
        vars = malloc(sizeof(Expr*) * (size_t)(k > 0 ? k : 1));
        for (int i = 0; i < k; i++) vars[i] = vars_arg->data.function.args[i];
    } else {
        k = 1;
        vars = malloc(sizeof(Expr*));
        vars[0] = vars_arg;
    }

    size_t m = rules->data.function.arg_count;
    Expr** terms = malloc(sizeof(Expr*) * (m > 0 ? m : 1));
    size_t nt = 0;

    for (size_t i = 0; i < m; i++) {
        Expr* rule = rules->data.function.args[i];
        if (!is_rule_head(rule)) goto fail;
        Expr* explist = rule->data.function.args[0];
        Expr* coeff   = rule->data.function.args[1];
        if (!is_list_head(explist) ||
            (int)explist->data.function.arg_count != k) goto fail;

        Expr** fac = malloc(sizeof(Expr*) * (size_t)(k + 1));
        size_t nf = 0;
        fac[nf++] = expr_copy(coeff);
        for (int j = 0; j < k; j++) {
            Expr* ej = explist->data.function.args[j];
            if (ej->type != EXPR_INTEGER) {
                for (size_t x = 0; x < nf; x++) expr_free(fac[x]);
                free(fac);
                goto fail;
            }
            int64_t e = ej->data.integer;
            if (e == 0) continue;
            if (e == 1) fac[nf++] = expr_copy(vars[j]);
            else        fac[nf++] = internal_power(
                            (Expr*[]){ expr_copy(vars[j]), expr_copy(ej) }, 2);
        }
        terms[nt++] = internal_times(fac, nf);
        free(fac);
    }

    Expr* result = (nt == 0) ? expr_new_integer(0) : internal_plus(terms, nt);
    free(terms);
    free(vars);
    return result;

fail:
    for (size_t x = 0; x < nt; x++) expr_free(terms[x]);
    free(terms);
    free(vars);
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  GroebnerBasis`DistributedTermsList                                 */
/* ------------------------------------------------------------------ */
/* GroebnerBasis`DistributedTermsList[poly, vars, MonomialOrder -> o] gives the
 * distributed form {{{e1, ..., ek}, c}, ...}, vars} -- CoefficientRules with
 * {exps, coeff} pairs instead of rules, the term list paired with the
 * variables, in descending monomial order (Lexicographic by default). A List
 * of polynomials gives one term list per polynomial. An internal
 * helper, which symbolic packages call directly (a loop-integral package reads
 * Feynman-parameter polynomials this way); a generalised variable it would
 * adjoin (y^(1/2) in vars {x, y}) is declined. */

/* One polynomial's {{exps, coeff}, ...}, or NULL when it is not polynomial. */
static Expr* dtl_terms(Expr* poly, Expr** vars, int k, const int64_t* W, int rows) {
    size_t n;
    Mono* monos = build_monomials(poly, vars, k, W, rows, NULL, &n);
    if (!monos) return NULL;
    Expr** out = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
    for (size_t i = 0; i < n; i++) {
        Expr** es = malloc(sizeof(Expr*) * (size_t)(k > 0 ? k : 1));
        for (int j = 0; j < k; j++) es[j] = expr_new_integer(monos[i].exps[j]);
        Expr* explist = expr_new_function(expr_new_symbol(SYM_List), es, (size_t)k);
        free(es);
        Expr* coeff = monos[i].coeff; monos[i].coeff = NULL;    /* move */
        out[i] = expr_new_function(expr_new_symbol(SYM_List),
                                   (Expr*[]){ explist, coeff }, 2);
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), out, n);
    free(out);
    monos_free(monos, n);
    return r;
}

Expr* builtin_distributedtermslist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    Expr* order_arg = NULL;
    size_t npos = argc;
    while (npos > 0 && is_rule_head(args[npos - 1])) {
        Expr* key = args[npos - 1]->data.function.args[0];
        if (key->type == EXPR_SYMBOL && key->data.symbol.name == SYM_MonomialOrder)
            order_arg = args[npos - 1]->data.function.args[1];
        npos--;
    }
    if (npos != 2 || !is_list_head(args[1])) return NULL;

    int k;
    Expr* owned = NULL;
    Expr** vars = resolve_vars(args[0], args[1], &k, &owned);
    int rows;
    int64_t* W = order_weight_matrix(order_arg, k, &rows);
    if (!W) { free(vars); return NULL; }

    Expr* terms;
    if (is_list_head(args[0])) {
        size_t m = args[0]->data.function.arg_count;
        Expr** each = malloc(sizeof(Expr*) * (m > 0 ? m : 1));
        size_t done = 0;
        for (; done < m; done++) {
            each[done] = dtl_terms(args[0]->data.function.args[done], vars, k, W, rows);
            if (!each[done]) break;
        }
        if (done < m) {
            for (size_t i = 0; i < done; i++) expr_free(each[i]);
            terms = NULL;
        } else {
            terms = expr_new_function(expr_new_symbol(SYM_List), each, m);
        }
        free(each);
    } else {
        terms = dtl_terms(args[0], vars, k, W, rows);
    }
    free(W); free(vars);
    if (!terms) return NULL;
    return expr_new_function(expr_new_symbol(SYM_List),
                             (Expr*[]){ terms, expr_copy(args[1]) }, 2);
}

/* ------------------------------------------------------------------ */
/*  Registration                                                       */
/* ------------------------------------------------------------------ */

void monomials_init(void) {
    symtab_add_builtin("MonomialList", builtin_monomiallist);
    symtab_get_def("MonomialList")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("CoefficientRules", builtin_coefficientrules);
    symtab_get_def("CoefficientRules")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("CoefficientArrays", builtin_coefficientarrays);
    symtab_get_def("CoefficientArrays")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CoefficientArrays",
        "CoefficientArrays[polys, vars] gives the coefficient arrays of the "
        "polynomial system polys in the variables vars, grouped by total degree: "
        "element d + 1 holds every degree-d coefficient, so {b, m} = "
        "CoefficientArrays[eqs, vars] is the constant vector and coefficient "
        "matrix of a linear system. A degree-d coefficient sits at the index "
        "tuple naming its variables in non-decreasing order, every other "
        "permutation being 0; a list of polynomials carries the equation index "
        "as the leading axis. CoefficientArrays[poly, vars] treats a non-List "
        "first argument as a single polynomial. Option Modulus -> p reduces the "
        "coefficients. Each array is a SparseArray (background 0); the "
        "constant term of a single polynomial is a scalar.");
    symtab_add_builtin("GroebnerBasis`DistributedTermsList", builtin_distributedtermslist);
    symtab_get_def("GroebnerBasis`DistributedTermsList")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("GroebnerBasis`DistributedTermsList",
        "GroebnerBasis`DistributedTermsList[poly, vars] gives {{{e1, ..., ek}, c}, ...}, vars}: "
        "the terms of the polynomial poly in vars as exponent vectors with their "
        "coefficients, in descending monomial order. A list of polynomials gives "
        "one term list each. Option MonomialOrder -> o (default Lexicographic) "
        "selects the order.");
    symtab_add_builtin("FromCoefficientRules", builtin_fromcoefficientrules);
    symtab_get_def("FromCoefficientRules")->attributes |= ATTR_PROTECTED;
}
