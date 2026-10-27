/*
 * facpoly_list.c -- FactorList[poly] / FactorList[poly, opts].
 *
 * A thin wrapper over Factor: factor `poly` via `Factor[poly, opts...]` (all
 * options are forwarded verbatim -- GaussianIntegers, Extension, ...), then
 * split the resulting product into {factor, exponent} pairs.
 *
 *   FactorList[x^2 - 1]              -> {{1, 1}, {-1 + x, 1}, {1 + x, 1}}
 *   FactorList[2 x^3 + 2 x^2 - ...]  -> {{2, 1}, {-1 + x, 1}, {1 + x, 2}}
 *   FactorList[x^4 - 2, Extension -> Sqrt[2]]
 *                                    -> {{1, 1}, {Sqrt[2] + x^2, 1}, {-Sqrt[2] + x^2, 1}}
 *
 * The first element is always the overall numerical factor {c, 1} (which is
 * {1, 1} when there is no numerical factor).  Denominator factors of a
 * rational function carry negative exponents.
 *
 * Parsing rules on the Factor output R (a Times, a bare factor, or a number):
 *   - a number literal (Integer / Rational / Real / Complex / ...) multiplies
 *     into the overall numerical factor `c`;
 *   - Power[base, e] with an *integer* e is the pair {base, e} (this is a
 *     factor raised to a multiplicity, positive or negative);
 *   - anything else -- including Power[base, 1/2] = Sqrt[base], which is an
 *     irreducible factor in its own right -- is the pair {factor, 1}.
 */

#include "facpoly_list.h"
#include "facpoly.h"
#include "expsplit.h"

#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "print.h"
#include "expr.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ===================================================================== */
/* Small helpers                                                         */
/* ===================================================================== */

static Expr* mk_int(long n) { return expr_new_integer(n); }

static Expr* mk_list2(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ a, b }, 2);
}

static bool is_rule_head(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rule ||
            e->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
           e->data.function.arg_count == 2;
}

static bool head_name_is(const Expr* e, const char* interned) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == interned;
}

/* A pure numeric literal: it belongs to the overall numerical factor. */
static bool is_number_literal(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
        case EXPR_REAL:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        default:
            break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        (head_name_is(e, SYM_Rational) || head_name_is(e, SYM_Complex))) {
        return true;
    }
    return false;
}

static bool is_integer_literal(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT);
}

/* ===================================================================== */
/* Diagnostics                                                           */
/* ===================================================================== */

static Expr* factorlist_emit_argx(size_t argc) {
    mth_message("FactorList", "argx", "FactorList called with %zu arguments; 1 argument is expected.", argc);
    return NULL;
}

static Expr* factorlist_emit_nonopt(Expr* bad, size_t pos, Expr* res) {
    char* bad_str = expr_to_string(bad);
    char* call_str = expr_to_string(res);
    mth_message("FactorList", "nonopt", "Options expected (instead of %s) beyond position %zu in %s. An option must be a rule or a list of rules.", bad_str ? bad_str : "?", pos, call_str ? call_str : "?");
    free(bad_str);
    free(call_str);
    return NULL;
}

/* ===================================================================== */
/* Factor-output -> {factor, exponent} pairs                             */
/* ===================================================================== */

/* Append the pairs for one multiplicative factor `f` to (bases,exps); or, when
 * f is a number literal, fold it into *numeric (consumed & replaced). */
static void absorb_factor(const Expr* f, Expr** numeric,
                          Expr*** bases, Expr*** exps, size_t* n, size_t* cap) {
    if (is_number_literal(f)) {
        Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                                       (Expr*[]){ *numeric, expr_copy((Expr*)f) }, 2);
        *numeric = eval_and_free(prod);
        return;
    }

    Expr* base; Expr* exp;
    if (head_name_is(f, SYM_Power) && f->data.function.arg_count == 2 &&
        is_integer_literal(f->data.function.args[1])) {
        base = expr_copy(f->data.function.args[0]);
        exp  = expr_copy(f->data.function.args[1]);
    } else {
        base = expr_copy((Expr*)f);
        exp  = mk_int(1);
    }

    if (*n == *cap) {
        *cap *= 2;
        *bases = realloc(*bases, *cap * sizeof(Expr*));
        *exps  = realloc(*exps,  *cap * sizeof(Expr*));
    }
    (*bases)[*n] = base;
    (*exps)[*n]  = exp;
    (*n)++;
}

/* Split an evaluated factored form `fac` (a Times, a bare factor, or a number)
 * into the pair-list {{numeric, 1}, {base_i, exp_i}, ...}.  Takes ownership of
 * `fac` (frees it).  Shared by FactorList and FactorSquareFreeList, whose only
 * difference is which head produced the product. */
/* Does the sum p lead with a negative coefficient -- the first term of
 * CoefficientRules[p, Variables[p]], lexicographic as in the reference? Then
 * FactorList reports a -1 for it: FactorList[1 - x] is {{-1, 1}, {-1 + x, 1}},
 * FactorList[x - y] keeps {x - y, 1}. (Factor's own normal form, -(1 - x)
 * = -1 + x times -1, is gone once the product evaluates: -(1 + x) is
 * -1 - x.) */
static bool sum_leads_negative(const Expr* p) {
    Expr* vars = eval_and_free(expr_new_function(expr_new_symbol(SYM_Variables),
                                                 (Expr*[]){ expr_copy((Expr*)p) }, 1));
    Expr* cr = eval_and_free(expr_new_function(expr_new_symbol(SYM_CoefficientRules),
                                               (Expr*[]){ expr_copy((Expr*)p), vars }, 2));
    bool neg = false;
    if (head_name_is(cr, SYM_List) && cr->data.function.arg_count >= 1) {
        Expr* r = cr->data.function.args[0];
        if (is_rule_head(r)) {
            Expr* c = r->data.function.args[1];
            neg = (c->type == EXPR_INTEGER && c->data.integer < 0) ||
                  (c->type == EXPR_BIGINT && mpz_sgn(c->data.bigint) < 0) ||
                  (head_name_is(c, SYM_Rational) && c->data.function.args[0]->type == EXPR_INTEGER &&
                   c->data.function.args[0]->data.integer < 0);
        }
    }
    expr_free(cr);
    return neg;
}

static Expr* product_to_pair_list(Expr* fac) {
    Expr* numeric = mk_int(1);
    size_t cap = 8, n = 0;
    Expr** bases = malloc(cap * sizeof(Expr*));
    Expr** exps  = malloc(cap * sizeof(Expr*));

    if (head_name_is(fac, SYM_Times)) {
        for (size_t i = 0; i < fac->data.function.arg_count; i++)
            absorb_factor(fac->data.function.args[i], &numeric,
                          &bases, &exps, &n, &cap);
    } else {
        absorb_factor(fac, &numeric, &bases, &exps, &n, &cap);
    }
    expr_free(fac);

    /* A simple factor whose sign the evaluated product absorbed -- -(1 + x)
     * is -1 - x, and (-1/2) (1 + x) is 1/2 (-1 - x) -- gets its -1 back into
     * the numerical factor, so its leading coefficient is positive again. */
    for (size_t i = 0; i < n; i++) {
        if (!(exps[i]->type == EXPR_INTEGER && exps[i]->data.integer == 1)) continue;
        if (!head_name_is(bases[i], SYM_Plus) || !sum_leads_negative(bases[i])) continue;
        Expr* p = bases[i];
        size_t m = p->data.function.arg_count;
        Expr** terms = malloc(sizeof(Expr*) * (m ? m : 1));
        for (size_t k = 0; k < m; k++)
            terms[k] = expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ mk_int(-1), expr_copy(p->data.function.args[k]) }, 2);
        bases[i] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, m));
        free(terms);
        expr_free(p);
        numeric = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                                  (Expr*[]){ mk_int(-1), numeric }, 2));
    }

    /* Assemble {{numeric, 1}, {base_i, exp_i}, ...}.  A rational numeric
     * factor p/q is listed as {p, 1}, {q, -1}, the {p, 1} dropped when p is
     * 1: FactorList[-15/2 (x^2 - 1)] starts {{-15, 1}, {2, -1}, ...}. */
    size_t total = n + 2;
    Expr** items = malloc(total * sizeof(Expr*));
    size_t k = 0;
    if (head_name_is(numeric, SYM_Rational) && numeric->data.function.arg_count == 2) {
        Expr* pn = expr_copy(numeric->data.function.args[0]);
        Expr* qd = expr_copy(numeric->data.function.args[1]);
        expr_free(numeric);
        if (pn->type == EXPR_INTEGER && pn->data.integer == 1) expr_free(pn);
        else items[k++] = mk_list2(pn, mk_int(1));
        items[k++] = mk_list2(qd, mk_int(-1));
    } else {
        items[k++] = mk_list2(numeric, mk_int(1));
    }
    for (size_t i = 0; i < n; i++)
        items[k++] = mk_list2(bases[i], exps[i]);
    total = k;
    free(bases);
    free(exps);

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, total);
    free(items);
    return out;
}

/* ===================================================================== */
/* Entry point                                                           */
/* ===================================================================== */

Expr* builtin_factorlist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return factorlist_emit_argx(0);

    /* Positions 2+ must be options (rules).  Report the last offender. */
    Expr* last_bad = NULL;
    for (size_t i = 1; i < argc; i++) {
        if (!is_rule_head(res->data.function.args[i]))
            last_bad = res->data.function.args[i];
    }
    if (last_bad) return factorlist_emit_nonopt(last_bad, 1, res);
    if (!factor_options_valid(res, "FactorList", 0)) return NULL;

    /* FactorList does not thread: a list, an equation or a logical
     * combination is one (unfactored) factor of itself. */
    {
        Expr* a0 = res->data.function.args[0];
        if (a0->type == EXPR_FUNCTION && a0->data.function.head->type == EXPR_SYMBOL) {
            const char* h = a0->data.function.head->data.symbol.name;
            if (h == SYM_List || h == SYM_Equal || h == SYM_Unequal || h == SYM_Less
                || h == SYM_LessEqual || h == SYM_Greater || h == SYM_GreaterEqual
                || h == SYM_And || h == SYM_Or || h == SYM_Not)
                return expr_new_function(expr_new_symbol(SYM_List),
                    (Expr*[]){ mk_list2(mk_int(1), mk_int(1)), mk_list2(expr_copy(a0), mk_int(1)) }, 2);
        }
    }

    /* (x + y)^(-3 + 3 e) lists as {(x + y)^e, 3}, {x + y, -3}: powers with a
     * symbolic exponent factor over their generators (expsplit.c), which come
     * right after the constant, in order of appearance, as in the reference. */
    {
        ExpSplitTab st = {0};
        Expr* split = expsplit_abstract(res->data.function.args[0], &st);
        if (split) {
            Expr** fa = malloc(argc * sizeof(Expr*));
            fa[0] = split;
            for (size_t i = 1; i < argc; i++) fa[i] = expr_copy(res->data.function.args[i]);
            Expr* fl = eval_and_free(expr_new_function(expr_new_symbol(SYM_FactorList), fa, argc));
            free(fa);
            Expr* out = NULL;
            if (fl && fl->type == EXPR_FUNCTION && fl->data.function.head->type == EXPR_SYMBOL &&
                fl->data.function.head->data.symbol.name == SYM_List) {
                size_t np = fl->data.function.arg_count;
                Expr** pairs = malloc((np ? np : 1) * sizeof(Expr*));
                size_t k = 0;
                /* the constant, then the generators in table order, then the rest */
                for (int pass = 0; pass < 3; pass++)
                    for (size_t i = 0; i < np; i++) {
                        Expr* pr = fl->data.function.args[i];
                        Expr* f = (pr->type == EXPR_FUNCTION && pr->data.function.arg_count == 2)
                                  ? pr->data.function.args[0] : NULL;
                        int g = f ? expsplit_index(&st, f) : -1;
                        int cls = (i == 0 && f && expr_is_numeric_like(f)) ? 0 : (g >= 0 ? 1 : 2);
                        if (cls != pass) continue;
                        if (g >= 0)
                            pairs[k++] = expr_new_function(expr_new_symbol(SYM_List),
                                (Expr*[]){ expsplit_generator(&st, (size_t)g),
                                           expr_copy(pr->data.function.args[1]) }, 2);
                        else
                            pairs[k++] = expsplit_restore(expr_copy(pr), &st);
                    }
                /* The ordinary factors keep FactorList's own order: refactor
                 * their product (it has no splittable power left), so the
                 * generator symbols' names cannot reorder x[3] after x[2] + x[3]. */
                size_t r0 = 0;
                while (r0 < k && expsplit_index(&st, pairs[r0]) < 0) {
                    const Expr* f0 = pairs[r0]->data.function.args[0];
                    bool gen = false;
                    for (size_t g = 0; g < st.n && !gen; g++) {
                        Expr* gx = expsplit_generator(&st, g);
                        gen = expr_eq(gx, (Expr*)f0);
                        expr_free(gx);
                    }
                    bool cst = (r0 == 0 && expr_is_numeric_like(f0));
                    if (!gen && !cst) break;
                    r0++;
                }
                if (k - r0 >= 2) {
                    Expr** pw = malloc((k - r0) * sizeof(Expr*));
                    for (size_t i = r0; i < k; i++)
                        pw[i - r0] = expr_new_function(expr_new_symbol(SYM_Power),
                            (Expr*[]){ expr_copy(pairs[i]->data.function.args[0]),
                                       expr_copy(pairs[i]->data.function.args[1]) }, 2);
                    Expr* prod = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), pw, k - r0));
                    free(pw);
                    Expr* rl = eval_and_free(expr_new_function(expr_new_symbol(SYM_FactorList), &prod, 1));
                    if (rl && rl->type == EXPR_FUNCTION && rl->data.function.arg_count >= 1 &&
                        rl->data.function.head->type == EXPR_SYMBOL &&
                        rl->data.function.head->data.symbol.name == SYM_List) {
                        size_t m = rl->data.function.arg_count;
                        Expr* c0 = rl->data.function.args[0];
                        bool unit = c0->type == EXPR_FUNCTION && c0->data.function.arg_count == 2 &&
                                    c0->data.function.args[0]->type == EXPR_INTEGER &&
                                    c0->data.function.args[0]->data.integer == 1;
                        size_t skip = unit ? 1 : 0;
                        if (m - skip == k - r0) {           /* same factors, FactorList's order */
                            for (size_t i = r0; i < k; i++) expr_free(pairs[i]);
                            for (size_t i = skip; i < m; i++)
                                pairs[r0 + i - skip] = expr_copy(rl->data.function.args[i]);
                        }
                    }
                    expr_free(rl);
                }
                out = expr_new_function(expr_new_symbol(SYM_List), pairs, k);
                free(pairs);
            }
            expr_free(fl);
            expsplit_free(&st);
            if (out) return out;
        }
    }

    /* Factor[poly, opts...] -- forward every argument verbatim. */
    Expr** fargs = malloc(argc * sizeof(Expr*));
    for (size_t i = 0; i < argc; i++)
        fargs[i] = expr_copy(res->data.function.args[i]);
    Expr* fac = eval_and_free(
        expr_new_function(expr_new_symbol(SYM_Factor), fargs, argc));
    free(fargs);
    if (!fac) return NULL;

    return product_to_pair_list(fac);
}

/* FactorSquareFreeList[poly] / FactorSquareFreeList[poly, Extension -> ext].
 *
 * A thin wrapper over FactorSquareFree: it square-free-factors via
 * `FactorSquareFree[poly, opts...]` (every argument -- poly and the Extension
 * option -- forwarded verbatim, so all extension handling lives in
 * FactorSquareFree) and splits the resulting product into {factor, exponent}
 * pairs, with the leading pair being the overall numerical factor {c, 1}. */
Expr* builtin_factorsquarefreelist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return NULL;  /* leave FactorSquareFreeList[] unevaluated */
    if (!factor_options_valid(res, "FactorSquareFreeList", 1)) return NULL;
    {
        Expr* a0 = res->data.function.args[0];
        if (a0->type == EXPR_FUNCTION && a0->data.function.head->type == EXPR_SYMBOL) {
            const char* h = a0->data.function.head->data.symbol.name;
            if (h == SYM_List || h == SYM_Equal || h == SYM_Unequal || h == SYM_Less
                || h == SYM_LessEqual || h == SYM_Greater || h == SYM_GreaterEqual
                || h == SYM_And || h == SYM_Or || h == SYM_Not)
                return expr_new_function(expr_new_symbol(SYM_List),
                    (Expr*[]){ mk_list2(mk_int(1), mk_int(1)), mk_list2(expr_copy(a0), mk_int(1)) }, 2);
        }
    }

    /* FactorSquareFree[poly, opts...] -- forward every argument verbatim. */
    Expr** fargs = malloc(argc * sizeof(Expr*));
    for (size_t i = 0; i < argc; i++)
        fargs[i] = expr_copy(res->data.function.args[i]);
    Expr* fac = eval_and_free(
        expr_new_function(expr_new_symbol(SYM_FactorSquareFree), fargs, argc));
    free(fargs);
    if (!fac) return NULL;

    /* If FactorSquareFree declined (bad arity / unrecognised option), it comes
     * back as an unevaluated FactorSquareFree[...] expression -- propagate the
     * decline so FactorSquareFreeList[...] stays unevaluated too, rather than
     * wrapping the inert call as a spurious factor. */
    if (head_name_is(fac, SYM_FactorSquareFree)) {
        expr_free(fac);
        return NULL;
    }

    return product_to_pair_list(fac);
}

/* ===================================================================== */
/* Init                                                                  */
/* ===================================================================== */

void factorlist_init(void) {
    symtab_add_builtin("FactorList", builtin_factorlist);
    symtab_get_def("FactorList")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c. */

    symtab_add_builtin("FactorSquareFreeList", builtin_factorsquarefreelist);
    symtab_get_def("FactorSquareFreeList")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c. */
}
