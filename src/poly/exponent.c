/*
 * exponent.c -- Exponent[expr, form] / Exponent[expr, form, h].
 *
 * Gives the maximum power (default h = Max) with which `form` appears in the
 * expanded form of `expr`.  Purely syntactic: `expr` is expanded first, but no
 * zero-coefficient recognition is attempted (Exponent[zero x^2 + x + 1, x] = 2
 * even when `zero` is numerically 0 but not in normal form).
 *
 * Algorithm
 * ---------
 *   1. expanded = Expand[expr].
 *   2. Split into additive terms (the args of a top-level Plus, else the whole
 *      expression as a single term).  The genuine zero polynomial has NO terms,
 *      so its exponent set is empty and h[] fires -- Max[] = -Infinity.
 *   3. For each term (a monomial), read off the exponent of `form`:
 *        - form decomposes into (base, fe) pairs (a symbol/kernel -> (form,1);
 *          a Power[b,e] -> (b,e); a product -> one pair per factor);
 *        - the exponent of a single base in a monomial is the power to which it
 *          is raised (0 if absent, symbolic/rational allowed);
 *        - for a product form the term's exponent is Min over the form's bases
 *          of (base-exponent / fe) -- the largest k with form^k dividing it.
 *   4. Collect the exponents into a sorted, de-duplicated set and return
 *      h @@ set (h defaults to Max).  Symbolic exponents flow through Max/Min
 *      unevaluated (Max[1/2, 1 + n]).
 *
 * Listable + Protected.  Listable makes Exponent[expr, {f1, f2, ...}] thread
 * into the per-form list of exponents for free.
 */

#include "exponent.h"

#include "expand.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "print.h"
#include "expr.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <gmp.h>

/* ===================================================================== */
/* Small helpers                                                         */
/* ===================================================================== */

static Expr* mk_int(long n) { return expr_new_integer(n); }

static bool head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

static bool is_rule_head(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rule ||
            e->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
           e->data.function.arg_count == 2;
}

/* True iff the expanded expression is the (numeric) zero polynomial. */
static bool expr_is_zero_number(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) == 0;
    if (e->type == EXPR_REAL)    return e->data.real == 0.0;
    return false;
}

/* ===================================================================== */
/* Diagnostics                                                           */
/* ===================================================================== */

static Expr* exponent_emit_argt(size_t argc) {
    fprintf(stderr,
            "Exponent::argt: Exponent called with %zu argument%s; "
            "2 or 3 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

static Expr* exponent_emit_nonopt(Expr* bad, size_t pos, Expr* res) {
    char* bad_str = expr_to_string(bad);
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "Exponent::nonopt: Options expected (instead of %s) beyond "
            "position %zu in %s. An option must be a rule or a list of rules.\n",
            bad_str ? bad_str : "?", pos, call_str ? call_str : "?");
    free(bad_str);
    free(call_str);
    return NULL;
}

/* ===================================================================== */
/* Exponent extraction                                                   */
/* ===================================================================== */

/* Exponent (owned Expr) of a single `base` in a monomial `term`.  For an
 * expanded monomial every base appears in at most one factor, but factors are
 * summed defensively.  Returns integer 0 when `base` does not appear. */
static Expr* base_exp_in_monomial(const Expr* term, const Expr* base) {
    if (expr_eq((Expr*)term, (Expr*)base)) return mk_int(1);

    if (head_is(term, "Power") && term->data.function.arg_count == 2 &&
        expr_eq(term->data.function.args[0], (Expr*)base)) {
        return expr_copy(term->data.function.args[1]);
    }

    if (head_is(term, "Times")) {
        size_t k = term->data.function.arg_count;
        Expr** parts = malloc(k * sizeof(Expr*));
        size_t np = 0;
        for (size_t i = 0; i < k; i++) {
            Expr* f = term->data.function.args[i];
            if (expr_eq(f, (Expr*)base)) {
                parts[np++] = mk_int(1);
            } else if (head_is(f, "Power") && f->data.function.arg_count == 2 &&
                       expr_eq(f->data.function.args[0], (Expr*)base)) {
                parts[np++] = expr_copy(f->data.function.args[1]);
            }
        }
        Expr* out;
        if (np == 0) { out = mk_int(0); }
        else if (np == 1) { out = parts[0]; }
        else {
            Expr* sum = expr_new_function(expr_new_symbol("Plus"), parts, np);
            out = eval_and_free(sum);
        }
        free(parts);
        return out;
    }

    return mk_int(0);
}

/* Decompose `form` into (base, fe) pairs.  Fills freshly-allocated parallel
 * arrays (caller frees the elements and the arrays).  A numeric form yields no
 * pairs (it never "appears"), so the exponent is 0. */
static void form_pairs(const Expr* form, Expr*** bases_out, Expr*** fes_out,
                       size_t* n_out) {
    size_t cap = 4, n = 0;
    Expr** bases = malloc(cap * sizeof(Expr*));
    Expr** fes = malloc(cap * sizeof(Expr*));

    #define PUSH(B, E) do { \
        if (n == cap) { cap *= 2; \
            bases = realloc(bases, cap * sizeof(Expr*)); \
            fes = realloc(fes, cap * sizeof(Expr*)); } \
        bases[n] = (B); fes[n] = (E); n++; } while (0)

    if (head_is(form, "Times")) {
        for (size_t i = 0; i < form->data.function.arg_count; i++) {
            Expr* f = form->data.function.args[i];
            if (f->type == EXPR_INTEGER || f->type == EXPR_BIGINT ||
                f->type == EXPR_REAL) {
                continue;  /* numeric factor: not a base */
            }
            if (head_is(f, "Power") && f->data.function.arg_count == 2) {
                PUSH(expr_copy(f->data.function.args[0]),
                     expr_copy(f->data.function.args[1]));
            } else {
                PUSH(expr_copy(f), mk_int(1));
            }
        }
    } else if (form->type == EXPR_INTEGER || form->type == EXPR_BIGINT ||
               form->type == EXPR_REAL) {
        /* numeric form: no bases */
    } else if (head_is(form, "Power") && form->data.function.arg_count == 2) {
        PUSH(expr_copy(form->data.function.args[0]),
             expr_copy(form->data.function.args[1]));
    } else {
        PUSH(expr_copy((Expr*)form), mk_int(1));
    }

    #undef PUSH
    *bases_out = bases;
    *fes_out = fes;
    *n_out = n;
}

/* Divide exponent `te` by the fixed form-exponent `fe` (owned inputs consumed);
 * returns owned result.  fe == 1 short-circuits. */
static Expr* ratio(Expr* te, Expr* fe) {
    if (fe->type == EXPR_INTEGER && fe->data.integer == 1) {
        expr_free(fe);
        return te;
    }
    Expr* inv = expr_new_function(expr_new_symbol("Power"),
                                  (Expr*[]){ fe, mk_int(-1) }, 2);
    Expr* prod = expr_new_function(expr_new_symbol("Times"),
                                   (Expr*[]){ te, inv }, 2);
    return eval_and_free(prod);
}

/* Exponent (owned Expr) of `form` in the monomial `term`. */
static Expr* exp_in_term(const Expr* term, const Expr* form) {
    Expr** bases; Expr** fes; size_t np;
    form_pairs(form, &bases, &fes, &np);

    Expr* out;
    if (np == 0) {
        out = mk_int(0);
    } else if (np == 1) {
        Expr* te = base_exp_in_monomial(term, bases[0]);
        out = ratio(te, fes[0]);  /* consumes te, fes[0] */
        fes[0] = NULL;
    } else {
        /* largest k with form^k | term  ==  Min_i (exp(base_i) / fe_i) */
        Expr** ratios = malloc(np * sizeof(Expr*));
        for (size_t i = 0; i < np; i++) {
            Expr* te = base_exp_in_monomial(term, bases[i]);
            ratios[i] = ratio(te, fes[i]);  /* consumes te, fes[i] */
            fes[i] = NULL;
        }
        Expr* mn = expr_new_function(expr_new_symbol("Min"), ratios, np);
        out = eval_and_free(mn);
        free(ratios);
    }

    for (size_t i = 0; i < np; i++) {
        expr_free(bases[i]);
        if (fes[i]) expr_free(fes[i]);
    }
    free(bases);
    free(fes);
    return out;
}

/* ===================================================================== */
/* Sorted, de-duplicated exponent set                                    */
/* ===================================================================== */

/* Insertion sort by expr_compare, then drop adjacent duplicates (expr_eq).
 * Small sets (one entry per distinct monomial degree); O(n^2) is fine. */
static void sort_dedup(Expr** a, size_t* n_io) {
    size_t n = *n_io;
    for (size_t i = 1; i < n; i++) {
        Expr* key = a[i];
        size_t j = i;
        while (j > 0 && expr_compare(a[j - 1], key) > 0) {
            a[j] = a[j - 1];
            j--;
        }
        a[j] = key;
    }
    size_t w = 0;
    for (size_t i = 0; i < n; i++) {
        if (w > 0 && expr_eq(a[w - 1], a[i])) {
            expr_free(a[i]);
        } else {
            a[w++] = a[i];
        }
    }
    *n_io = w;
}

/* ===================================================================== */
/* Entry point                                                           */
/* ===================================================================== */

Expr* builtin_exponent(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return exponent_emit_argt(argc);

    /* Positions 4+ must be options (rules); we implement none, so any
     * non-rule there is a nonopt error.  Report the last offender at pos 3. */
    if (argc > 3) {
        Expr* last_bad = NULL;
        for (size_t i = 3; i < argc; i++) {
            if (!is_rule_head(res->data.function.args[i]))
                last_bad = res->data.function.args[i];
        }
        if (last_bad) return exponent_emit_nonopt(last_bad, 3, res);
    }

    Expr* expr = res->data.function.args[0];
    Expr* form = res->data.function.args[1];
    Expr* h    = (argc >= 3) ? res->data.function.args[2] : NULL;

    Expr* expanded = expr_expand(expr);
    if (!expanded) return NULL;

    size_t cap = 8, ne = 0;
    Expr** exps = malloc(cap * sizeof(Expr*));

    if (!expr_is_zero_number(expanded)) {
        if (head_is(expanded, "Plus")) {
            size_t nt = expanded->data.function.arg_count;
            if (nt > cap) { cap = nt; exps = realloc(exps, cap * sizeof(Expr*)); }
            for (size_t i = 0; i < nt; i++)
                exps[ne++] = exp_in_term(expanded->data.function.args[i], form);
        } else {
            exps[ne++] = exp_in_term(expanded, form);
        }
    }
    expr_free(expanded);

    sort_dedup(exps, &ne);

    Expr* hhead = h ? expr_copy(h) : expr_new_symbol("Max");
    Expr* call = expr_new_function(hhead, exps, ne);  /* adopts exps elements */
    free(exps);
    return eval_and_free(call);
}

/* ===================================================================== */
/* Init                                                                  */
/* ===================================================================== */

void exponent_init(void) {
    symtab_add_builtin("Exponent", builtin_exponent);
    symtab_get_def("Exponent")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    /* Docstring lives in info.c. */
}
