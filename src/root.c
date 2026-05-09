/* root.c
 *
 * Mathematica-style symbolic Root and RootSum.  See root.h for the
 * representation contract.  Both heads are held forms — the C
 * builtins below intentionally return NULL for every input so the
 * evaluator leaves the call exactly as the caller wrote it.  The
 * useful work happens elsewhere:
 *
 *   - src/deriv.c   — D[RootSum[f1, f2], x] threads through the body.
 *   - src/intrat.c  — NaiveLogPart constructs RootSum nodes when the
 *                     LRT log part has no closed-form real expression.
 *
 * Splitting this module out of intrat.c keeps the held-symbolic
 * machinery available to any future caller (e.g. Solve, Reduce,
 * factorisation over algebraic extensions) that needs to name the
 * roots of a polynomial without committing to a radical form.
 */

#include "root.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"

#include <stdlib.h>
#include <string.h>

/* Held builtin: the evaluator already left the call unevaluated by
 * the time we get here, so there is nothing to compute.  Returning
 * NULL preserves the input verbatim. */
static Expr* builtin_root(Expr* res) {
    (void)res;
    return NULL;
}

/* ------------------------------------------------------------------
 * RootSum simplifier: Lagrange-interpolation closed form.
 *
 * For a squarefree polynomial d(t) with roots α_i and any polynomial
 * a(t) with deg(a) < deg(d), the Hermite/Lagrange interpolation
 * identity
 *
 *     Σ_i  a(α_i) / (d'(α_i) (x − α_i))  ==  a(x) / d(x)
 *
 * collapses any RootSum whose body matches the shape
 *     Function[ a(#) / (d'(#) (x − #)) ]
 * to the closed rational function a(x) / d(x).
 *
 * This is the post-D form produced by D[RootSum[Function[d],
 * Function[a Log[x − #]/d'(#)]], x] and is the identity Mathematica
 * uses to recognise that D[RootSum[1+#+#^6&, Log[x−#]/(1+6 #^5)&], x]
 * == 1/(1 + x + x^6).
 * ------------------------------------------------------------------ */

/* True iff e contains no subexpression structurally equal to target. */
static bool expr_free_of(const Expr* e, const Expr* target) {
    if (expr_eq(e, target)) return false;
    if (e->type == EXPR_FUNCTION) {
        if (!expr_free_of(e->data.function.head, target)) return false;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (!expr_free_of(e->data.function.args[i], target)) return false;
        }
    }
    return true;
}

/* Substitute every Slot[1] in `e` with a deep copy of `replacement`. */
static Expr* subst_slot1(const Expr* e, const Expr* replacement) {
    if (!e) return NULL;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "Slot") == 0
        && e->data.function.arg_count == 1
        && e->data.function.args[0]->type == EXPR_INTEGER
        && e->data.function.args[0]->data.integer == 1) {
        return expr_copy((Expr*)replacement);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        new_args[i] = subst_slot1(e->data.function.args[i], replacement);
    }
    Expr* new_head = subst_slot1(e->data.function.head, replacement);
    Expr* out = expr_new_function(new_head, new_args, n);
    free(new_args);
    return out;
}

/* Walk `e` looking for the literal subexpression Plus[x, Times[-1,
 * Slot[1]]] (i.e. `x - #1`) for some `x` free of Slot[1].  On match,
 * returns a fresh copy of `x`; otherwise NULL.  Picks the first
 * match found. */
static Expr* find_x_minus_slot1(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Plus
        && e->data.function.arg_count == 2) {
        /* Look for one arg = -Slot[1] (i.e. Times[-1, Slot[1]]) and
         * the other free of Slot[1]. */
        Expr* a = e->data.function.args[0];
        Expr* b = e->data.function.args[1];
        Expr* slot = expr_new_function(expr_new_symbol("Slot"),
            (Expr*[]){ expr_new_integer(1) }, 1);
        Expr* neg_slot = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(slot) }, 2);
        Expr* x_candidate = NULL;
        if (expr_eq(a, neg_slot) && expr_free_of(b, slot)) {
            x_candidate = expr_copy(b);
        } else if (expr_eq(b, neg_slot) && expr_free_of(a, slot)) {
            x_candidate = expr_copy(a);
        }
        expr_free(slot); expr_free(neg_slot);
        if (x_candidate) return x_candidate;
    }
    /* Recurse. */
    Expr* found = find_x_minus_slot1(e->data.function.head);
    if (found) return found;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        found = find_x_minus_slot1(e->data.function.args[i]);
        if (found) return found;
    }
    return NULL;
}

/* Compute Cancel[Together[e]]. */
static Expr* simplify_rational(Expr* e) {
    Expr* tg = internal_together((Expr*[]){ expr_copy(e) }, 1);
    Expr* tg_e = evaluate(tg); expr_free(tg);
    Expr* ca = internal_cancel((Expr*[]){ tg_e }, 1);
    return ca;
}

/* Try the Lagrange closed-form simplification.  Returns the
 * simplified Expr* on success, NULL otherwise.  Does not consume
 * its inputs. */
static Expr* rootsum_try_lagrange(Expr* poly_fn, Expr* body_fn) {
    /* Both args must be 1-arg Function nodes (Slot[1] form). */
    if (poly_fn->type != EXPR_FUNCTION
        || poly_fn->data.function.head->type != EXPR_SYMBOL
        || poly_fn->data.function.head->data.symbol != SYM_Function
        || poly_fn->data.function.arg_count != 1) return NULL;
    if (body_fn->type != EXPR_FUNCTION
        || body_fn->data.function.head->type != EXPR_SYMBOL
        || body_fn->data.function.head->data.symbol != SYM_Function
        || body_fn->data.function.arg_count != 1) return NULL;

    Expr* poly = poly_fn->data.function.args[0];
    Expr* body = body_fn->data.function.args[0];

    /* Find the variable x such that (x - #1) appears in the body. */
    Expr* x = find_x_minus_slot1(body);
    if (!x) return NULL;

    /* d_prime = D[poly, Slot[1]]. */
    Expr* slot = expr_new_function(expr_new_symbol("Slot"),
        (Expr*[]){ expr_new_integer(1) }, 1);
    Expr* d_call = expr_new_function(expr_new_symbol("D"),
        (Expr*[]){ expr_copy(poly), expr_copy(slot) }, 2);
    Expr* d_prime = evaluate(d_call); expr_free(d_call);

    /* a(#1) = body * (x - #1) * d'(#1).  Use the Lagrange identity
     * Σ_α body[α] = a(x)/d(x) iff body = a(t)/(d'(t)(x-t)). */
    Expr* x_minus_slot = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(x),
                   expr_new_function(expr_new_symbol("Times"),
                       (Expr*[]){ expr_new_integer(-1),
                                  expr_copy(slot) }, 2) }, 2);
    Expr* prod_raw = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_copy(body), x_minus_slot, d_prime }, 3);
    Expr* prod_eval = evaluate(prod_raw); expr_free(prod_raw);
    Expr* a_of_slot = simplify_rational(prod_eval); expr_free(prod_eval);

    /* If the simplified product still contains x, the body did not
     * fit the Lagrange shape — bail out. */
    if (!expr_free_of(a_of_slot, x)) {
        expr_free(a_of_slot); expr_free(slot); expr_free(x);
        return NULL;
    }

    /* a(x) = a_of_slot[Slot[1] -> x],  d(x) = poly[Slot[1] -> x]. */
    Expr* a_at_x = subst_slot1(a_of_slot, x);
    Expr* d_at_x = subst_slot1(poly, x);
    expr_free(a_of_slot); expr_free(slot);

    /* Build a(x)/d(x).  Run Together so the rational form prints
     * canonically. */
    Expr* quotient_raw = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ a_at_x,
                   expr_new_function(expr_new_symbol("Power"),
                       (Expr*[]){ d_at_x, expr_new_integer(-1) }, 2) }, 2);
    Expr* quotient_eval = evaluate(quotient_raw); expr_free(quotient_raw);
    Expr* tg = internal_together((Expr*[]){ quotient_eval }, 1);
    Expr* result = evaluate(tg); expr_free(tg);
    expr_free(x);
    return result;
}

static Expr* builtin_rootsum(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* fn1 = res->data.function.args[0];
    Expr* fn2 = res->data.function.args[1];
    Expr* simplified = rootsum_try_lagrange(fn1, fn2);
    if (simplified) return simplified;
    return NULL;
}

/* Replace every occurrence of `bvar` (a symbol) in `e` with Slot[1].
 * Walks the tree in place, returning a freshly allocated copy. */
static Expr* substitute_bvar_with_slot(Expr* e, const char* bvar_name) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL && strcmp(e->data.symbol, bvar_name) == 0) {
        return expr_new_function(expr_new_symbol("Slot"),
            (Expr*[]){ expr_new_integer(1) }, 1);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    Expr* head = substitute_bvar_with_slot(e->data.function.head, bvar_name);
    size_t n = e->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        args[i] = substitute_bvar_with_slot(e->data.function.args[i], bvar_name);
    }
    Expr* out = expr_new_function(head, args, n);
    free(args);
    return out;
}

Expr* root_make_rootsum(Expr* bvar, Expr* poly, Expr* body) {
    /* Convert to the Mathematica-canonical Slot form:
     *   RootSum[Function[poly_in_slot1], Function[body_in_slot1]]
     * where every occurrence of `bvar` in `poly` and `body` becomes
     * `Slot[1]`.  This matches `(... &)` syntax (1-arg Function) and
     * keeps the held output free of internal context-qualified
     * symbols like Integrate`Private`t.  The 2-arg `Function[t, expr]`
     * form is also accepted by D[RootSum, x] in src/deriv.c for back-
     * compat, but the Slot form is what we emit going forward. */
    if (bvar && bvar->type == EXPR_SYMBOL) {
        const char* name = bvar->data.symbol;
        Expr* poly_s = substitute_bvar_with_slot(poly, name);
        Expr* body_s = substitute_bvar_with_slot(body, name);
        expr_free(bvar); expr_free(poly); expr_free(body);

        Expr* fn1 = expr_new_function(expr_new_symbol("Function"),
            (Expr*[]){ poly_s }, 1);
        Expr* fn2 = expr_new_function(expr_new_symbol("Function"),
            (Expr*[]){ body_s }, 1);
        return expr_new_function(expr_new_symbol("RootSum"),
            (Expr*[]){ fn1, fn2 }, 2);
    }

    /* Defensive fallback for non-symbol bvar: original 2-arg form. */
    Expr* fn1_args[2] = { bvar, poly };
    Expr* fn1 = expr_new_function(expr_new_symbol("Function"), fn1_args, 2);

    Expr* fn2_args[2] = { expr_copy(bvar), body };
    Expr* fn2 = expr_new_function(expr_new_symbol("Function"), fn2_args, 2);

    Expr* rs_args[2] = { fn1, fn2 };
    return expr_new_function(expr_new_symbol("RootSum"), rs_args, 2);
}

void root_init(void) {
    symtab_add_builtin("Root", builtin_root);
    symtab_get_def("Root")->attributes |= ATTR_PROTECTED | ATTR_HOLDALL;
    symtab_set_docstring("Root",
        "Root[Function[t, p[t]], k]\n"
        "\tRepresents the k-th real root of the univariate polynomial p\n"
        "\tin the variable t.  Held symbolic form — the system makes no\n"
        "\tattempt to express the root as a closed-form radical.");

    symtab_add_builtin("RootSum", builtin_rootsum);
    symtab_get_def("RootSum")->attributes |= ATTR_PROTECTED | ATTR_HOLDALL;
    symtab_set_docstring("RootSum",
        "RootSum[Function[t, p[t]], Function[t, body[t]]]\n"
        "\tThe formal sum of body[\\[Alpha]] over the roots \\[Alpha] of\n"
        "\tp[\\[Alpha]] == 0.  Held symbolic form, used by the rational\n"
        "\tintegrator's NaiveLogPart fallback when the logarithmic part\n"
        "\tcannot be expressed in closed-form real elementary functions.\n"
        "\tDifferentiation threads through the body Function:\n"
        "\t  D[RootSum[f1, Function[t, body]], x]\n"
        "\t    == RootSum[f1, Function[t, D[body, x]]].");
}
