/* groebnerbasis.c
 *
 * `GroebnerBasis[polys, vars]` and `GroebnerBasis[polys, mainVars,
 * elimVars]` -- the user-facing entry point that drives the Buchberger
 * core in groebner.c.  See SPEC.md and docs/spec/builtins/arithmetic-
 * and-algebra.md for the surface this implements.
 *
 * Supported options (others emit `GroebnerBasis::nimpl`):
 *   MonomialOrder      -> Lexicographic            (default)
 *                       | DegreeReverseLexicographic
 *                       | EliminationOrder         (forced by the 3-arg form)
 *   CoefficientDomain  -> Rationals | Automatic    (default)
 *   Method             -> "Buchberger" | Automatic (default)
 */

#include "groebnerbasis.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "eval.h"
#include "groebner.h"
#include "internal.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ */
/*  Diagnostic helpers                                                 */
/* ------------------------------------------------------------------ */

static void warn_once(const char* tag, const char* msg) {
    fprintf(stderr, "GroebnerBasis::%s: %s\n", tag, msg);
}

/* ------------------------------------------------------------------ */
/*  Argument shape detection                                           */
/* ------------------------------------------------------------------ */

static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_List;
}

static bool is_rule(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && (e->data.function.head->data.symbol == SYM_Rule
         || e->data.function.head->data.symbol == SYM_RuleDelayed);
}

static bool is_equal(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Equal;
}

/* ------------------------------------------------------------------ */
/*  Option extraction                                                  */
/* ------------------------------------------------------------------ */

/* In-out: `*n_pos` starts at the full argument count; the function
 * decrements it for each trailing Rule[] arg.  All three output options
 * default to Automatic-equivalent values.  Returns `false` on a setting
 * we cannot honour (caller chooses whether to bail or fall back). */
static bool extract_options(const Expr* res, size_t* n_pos,
                            GBOrder* order, bool* warned_nimpl) {
    size_t argc = res->data.function.arg_count;
    size_t cut = argc;
    while (cut > 0 && is_rule(res->data.function.args[cut - 1])) cut--;
    *n_pos = cut;

    /* Defaults. */
    *order = GB_ORDER_LEX;

    for (size_t i = cut; i < argc; i++) {
        Expr* rule = res->data.function.args[i];
        if (rule->data.function.arg_count != 2) continue;
        Expr* key = rule->data.function.args[0];
        Expr* val = rule->data.function.args[1];
        if (key->type != EXPR_SYMBOL) continue;

        if (key->data.symbol == SYM_MonomialOrder) {
            if (val->type == EXPR_SYMBOL) {
                if (val->data.symbol == SYM_Lexicographic
                 || val->data.symbol == SYM_Automatic) {
                    *order = GB_ORDER_LEX;
                } else if (val->data.symbol == SYM_DegreeReverseLexicographic) {
                    *order = GB_ORDER_GREVLEX;
                } else if (val->data.symbol == SYM_EliminationOrder) {
                    *order = GB_ORDER_ELIM;
                } else {
                    warn_once("nimpl", "unsupported MonomialOrder value; "
                                       "falling back to Lexicographic");
                    *warned_nimpl = true;
                    *order = GB_ORDER_LEX;
                }
            } else {
                warn_once("nimpl", "weight-matrix MonomialOrder not "
                                   "implemented; falling back to Lexicographic");
                *warned_nimpl = true;
            }
        } else if (key->data.symbol == SYM_CoefficientDomain) {
            bool ok = (val->type == EXPR_SYMBOL
                       && (val->data.symbol == SYM_Rationals
                        || val->data.symbol == SYM_Automatic));
            if (!ok) {
                warn_once("nimpl", "only CoefficientDomain -> Rationals is "
                                   "implemented");
                *warned_nimpl = true;
            }
        } else if (key->data.symbol == SYM_Method) {
            bool ok = false;
            if (val->type == EXPR_SYMBOL
                && val->data.symbol == SYM_Automatic) ok = true;
            if (val->type == EXPR_STRING && val->data.string
                && strcmp(val->data.string, "Buchberger") == 0) ok = true;
            if (!ok) {
                warn_once("nimpl", "only Method -> \"Buchberger\" is "
                                   "implemented; falling back");
                *warned_nimpl = true;
            }
        }
    }
    return true;
}

/* ------------------------------------------------------------------ */
/*  Polynomial-equation preprocessing                                  */
/* ------------------------------------------------------------------ */

/* Rewrite `Equal[a, b]` as `a - b`, then run `Expand` so products of
 * polynomials (e.g. `(x - 1) (x - 2)`) reach a Plus-of-monomials shape
 * that `gb_from_expr` can ingest.  Other shapes are expanded as well
 * since input polynomials may carry factored sub-expressions.  Caller
 * owns the return. */
static Expr* normalise_polynomial(const Expr* p) {
    Expr* base;
    if (is_equal(p) && p->data.function.arg_count == 2) {
        Expr* a = expr_copy(p->data.function.args[0]);
        Expr* b = expr_copy(p->data.function.args[1]);
        base = internal_subtract((Expr*[]){ a, b }, 2);
    } else {
        base = expr_copy((Expr*)p);
    }
    /* Expand factored products / collect like terms.  Expand is also a
     * cheap no-op on already-expanded inputs. */
    Expr* expanded = internal_expand((Expr*[]){ base }, 1);
    /* Finally evaluate to fix-point so Plus/Times normalise. */
    Expr* normalised = evaluate(expanded);
    expr_free(expanded);
    return normalised;
}

/* ------------------------------------------------------------------ */
/*  Builtin entry                                                      */
/* ------------------------------------------------------------------ */

Expr* builtin_groebner_basis(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) {
        if (argc == 0) {
            fprintf(stderr, "GroebnerBasis::argt: GroebnerBasis called with "
                            "0 arguments; 2 or 3 expected.\n");
        } else {
            fprintf(stderr, "GroebnerBasis::argt: GroebnerBasis called with "
                            "1 argument; 2 or 3 expected.\n");
        }
        return NULL;
    }

    size_t n_pos;
    GBOrder order;
    bool warned_nimpl = false;
    extract_options(res, &n_pos, &order, &warned_nimpl);

    if (n_pos < 2 || n_pos > 3) {
        fprintf(stderr, "GroebnerBasis::argt: GroebnerBasis takes 2 or 3 "
                        "positional arguments.\n");
        return NULL;
    }

    Expr* polys_list = res->data.function.args[0];
    Expr* vars_list  = res->data.function.args[1];
    Expr* elim_list  = (n_pos == 3) ? res->data.function.args[2] : NULL;

    if (!is_list(polys_list)) return NULL;

    /* Accept a single-symbol shorthand for the variable list:
     *   GroebnerBasis[polys, x]    -> GroebnerBasis[polys, {x}]
     *   GroebnerBasis[polys, m, x] -> GroebnerBasis[polys, m, {x}]
     * We wrap symbol args in a synthetic List for uniform downstream
     * processing.  The wrappers are freed at the end of the function. */
    Expr* wrap_vars = NULL;
    Expr* wrap_elim = NULL;
    if (!is_list(vars_list)) {
        if (vars_list->type != EXPR_SYMBOL) return NULL;
        Expr** wrapped = (Expr**)malloc(sizeof(Expr*));
        wrapped[0] = expr_copy(vars_list);
        wrap_vars = expr_new_function(expr_new_symbol("List"), wrapped, 1);
        free(wrapped);
        vars_list = wrap_vars;
    }
    if (n_pos == 3 && !is_list(elim_list)) {
        if (elim_list->type != EXPR_SYMBOL) {
            if (wrap_vars) expr_free(wrap_vars);
            return NULL;
        }
        Expr** wrapped = (Expr**)malloc(sizeof(Expr*));
        wrapped[0] = expr_copy(elim_list);
        wrap_elim = expr_new_function(expr_new_symbol("List"), wrapped, 1);
        free(wrapped);
        elim_list = wrap_elim;
    }

    /* Empty polynomial list -> empty basis. */
    size_t n_polys = polys_list->data.function.arg_count;
    if (n_polys == 0) {
        if (wrap_vars) expr_free(wrap_vars);
        if (wrap_elim) expr_free(wrap_elim);
        return expr_new_function(expr_new_symbol("List"), NULL, 0);
    }

    /* Build the joint variable array: elim variables first when in the
     * 3-arg form, then the main variables.  The elimination order then
     * lex-prefers monomials carrying elim-block exponents. */
    size_t n_main = vars_list->data.function.arg_count;
    size_t n_elim = elim_list ? elim_list->data.function.arg_count : 0;
    size_t n_vars = n_main + n_elim;
    if (n_vars == 0) {
        if (wrap_vars) expr_free(wrap_vars);
        if (wrap_elim) expr_free(wrap_elim);
        return NULL;
    }

    Expr** all_vars = (Expr**)malloc(sizeof(Expr*) * n_vars);
    for (size_t i = 0; i < n_elim; i++) {
        all_vars[i] = elim_list->data.function.args[i];
    }
    for (size_t i = 0; i < n_main; i++) {
        all_vars[n_elim + i] = vars_list->data.function.args[i];
    }
    int elim_pivot = (int)n_elim;

    /* The 3-arg form forces Lex (with elim variables placed first) so
     * the elimination theorem applies: any monomial involving any elim
     * variable is lex-larger than any monomial without one.  The user's
     * MonomialOrder option is informational in this case (Mathematica
     * behaves the same way). */
    GBOrder use_order = (n_elim > 0) ? GB_ORDER_LEX : order;

    /* Convert each input polynomial. */
    GBPoly** F = (GBPoly**)malloc(sizeof(GBPoly*) * n_polys);
    size_t nF = 0;
    bool failed = false;
    for (size_t i = 0; i < n_polys; i++) {
        Expr* norm = normalise_polynomial(polys_list->data.function.args[i]);
        GBPoly* p = gb_from_expr(norm, all_vars, (int)n_vars,
                                 use_order, elim_pivot);
        expr_free(norm);
        if (!p) { failed = true; break; }
        if (p->n_terms == 0) { gb_poly_free(p); continue; }
        F[nF++] = p;
    }
    if (failed) {
        for (size_t i = 0; i < nF; i++) gb_poly_free(F[i]);
        free(F);
        free(all_vars);
        if (wrap_vars) expr_free(wrap_vars);
        if (wrap_elim) expr_free(wrap_elim);
        return NULL;
    }

    if (nF == 0) {
        /* All inputs reduced to zero -> empty basis. */
        free(F);
        free(all_vars);
        if (wrap_vars) expr_free(wrap_vars);
        if (wrap_elim) expr_free(wrap_elim);
        return expr_new_function(expr_new_symbol("List"), NULL, 0);
    }

    /* Fast path: any input is a non-zero constant -> ideal = <1>. */
    for (size_t i = 0; i < nF; i++) {
        if (gb_poly_is_constant(F[i])) {
            for (size_t j = 0; j < nF; j++) gb_poly_free(F[j]);
            free(F);
            free(all_vars);
            if (wrap_vars) expr_free(wrap_vars);
            if (wrap_elim) expr_free(wrap_elim);
            Expr** one = (Expr**)malloc(sizeof(Expr*));
            one[0] = expr_new_integer(1);
            Expr* lst = expr_new_function(expr_new_symbol("List"), one, 1);
            free(one);
            return lst;
        }
    }

    size_t out_n = 0;
    GBPoly** G = gb_buchberger(F, nF, &out_n);

    /* Free input working set. */
    for (size_t i = 0; i < nF; i++) gb_poly_free(F[i]);
    free(F);

    /* Elimination filter: drop polynomials that still mention any elim
     * variable.  By the elimination theorem, the surviving polynomials
     * form a Gröbner basis of the elimination ideal. */
    if (n_elim > 0) {
        int* elim_idx = (int*)malloc(sizeof(int) * n_elim);
        for (size_t i = 0; i < n_elim; i++) elim_idx[i] = (int)i;
        size_t k = 0;
        for (size_t i = 0; i < out_n; i++) {
            if (gb_poly_free_of_vars(G[i], elim_idx, (int)n_elim)) {
                G[k++] = G[i];
            } else {
                gb_poly_free(G[i]);
            }
        }
        out_n = k;
        free(elim_idx);
    }

    /* If the basis collapses to {<non-zero constant>} -> {1}. */
    bool has_const = false;
    for (size_t i = 0; i < out_n; i++) {
        if (gb_poly_is_constant(G[i]) && !gb_poly_is_zero(G[i])) {
            has_const = true; break;
        }
    }
    if (has_const) {
        gb_basis_free(G, out_n);
        free(all_vars);
        if (wrap_vars) expr_free(wrap_vars);
        if (wrap_elim) expr_free(wrap_elim);
        Expr** one = (Expr**)malloc(sizeof(Expr*));
        one[0] = expr_new_integer(1);
        Expr* lst = expr_new_function(expr_new_symbol("List"), one, 1);
        free(one);
        return lst;
    }

    /* Convert each basis polynomial back to an Expr (using the elim-
     * filtered variable view: the elim variables are now unused, but
     * gb_to_expr will skip exponent-0 vars anyway). */
    Expr** items = (out_n > 0) ? (Expr**)malloc(sizeof(Expr*) * out_n) : NULL;
    for (size_t i = 0; i < out_n; i++) {
        items[i] = gb_to_expr(G[i], all_vars);
    }
    gb_basis_free(G, out_n);
    free(all_vars);
    if (wrap_vars) expr_free(wrap_vars);
    if (wrap_elim) expr_free(wrap_elim);

    Expr* lst = expr_new_function(expr_new_symbol("List"), items, out_n);
    free(items);
    return lst;
}

/* ------------------------------------------------------------------ */
/*  Initialisation                                                     */
/* ------------------------------------------------------------------ */

void groebner_init(void) {
    symtab_add_builtin("GroebnerBasis", builtin_groebner_basis);
    SymbolDef* def = symtab_get_def("GroebnerBasis");
    if (def) def->attributes |= ATTR_PROTECTED;
}
