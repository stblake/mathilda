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
#include "attr.h"

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
    /* Pure lexicographic / negative-lexicographic: k identity(-signed) rows. */
    if (deg_sign == 0) {
        int64_t* w = calloc((size_t)(k > 0 ? k : 1) * (size_t)(k > 0 ? k : 1),
                            sizeof(int64_t));
        for (int i = 0; i < k; i++) w[(size_t)i * k + i] = neg_lex ? -1 : 1;
        *rows_out = k;
        return w;
    }
    /* Degree-first orders: row 0 is the (signed) total-degree row, then k-1
     * lex or reverse-lex tail rows. */
    int rows = k;                         /* deg row + (k-1) tail rows */
    int64_t* w = calloc((size_t)(rows > 0 ? rows : 1) * (size_t)(k > 0 ? k : 1),
                        sizeof(int64_t));
    for (int j = 0; j < k; j++) w[j] = deg_sign;        /* total-degree row */
    for (int r = 1; r < k; r++) {
        if (rev) w[(size_t)r * k + (k - r)] = -1;        /* -e_{k-r} */
        else     w[(size_t)r * k + (r - 1)] =  1;        /*  e_{r-1} */
    }
    *rows_out = rows;
    return w;
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
    if (order_arg->type == EXPR_STRING) {
        const char* s = order_arg->data.string;
        if (strcmp(s, "Lexicographic") == 0)
            return named_matrix(k, 0, 0, 0, rows_out);
        if (strcmp(s, "NegativeLexicographic") == 0)
            return named_matrix(k, 0, 0, 1, rows_out);
        if (strcmp(s, "DegreeLexicographic") == 0)
            return named_matrix(k, +1, 0, 0, rows_out);
        if (strcmp(s, "DegreeReverseLexicographic") == 0)
            return named_matrix(k, +1, 1, 0, rows_out);
        if (strcmp(s, "NegativeDegreeLexicographic") == 0)
            return named_matrix(k, -1, 0, 0, rows_out);
        if (strcmp(s, "NegativeDegreeReverseLexicographic") == 0)
            return named_matrix(k, -1, 1, 0, rows_out);
        return NULL;                                /* unknown order string */
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
static Mono* build_monomials(Expr* poly, Expr** vars, int k,
                             const int64_t* W, int rows, Expr* modulus,
                             size_t* n_out) {
    *n_out = 0;

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

    /* Sort, then merge byte-equal exponent vectors (adjacent after the sort)
     * and drop zero coefficients. */
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
/*  FromCoefficientRules                                               */
/* ------------------------------------------------------------------ */

Expr* builtin_fromcoefficientrules(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) {
        fprintf(stderr,
                "FromCoefficientRules::argrx: FromCoefficientRules called with "
                "%zu arguments; 2 arguments are expected.\n", argc);
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
/*  Registration                                                       */
/* ------------------------------------------------------------------ */

void monomials_init(void) {
    symtab_add_builtin("MonomialList", builtin_monomiallist);
    symtab_get_def("MonomialList")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("CoefficientRules", builtin_coefficientrules);
    symtab_get_def("CoefficientRules")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("FromCoefficientRules", builtin_fromcoefficientrules);
    symtab_get_def("FromCoefficientRules")->attributes |= ATTR_PROTECTED;
}
