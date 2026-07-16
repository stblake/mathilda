/* integrate_symmetry.c
 *
 * Definite integration by origin-symmetry reduction.  See integrate_symmetry.h
 * for the contract.  The two reductions are:
 *
 *   odd  f(-x) == -f(x):  Int_{-c}^{c} f = 0        (half must converge)
 *   even f(-x) ==  f(x):  Int_{-c}^{c} f = 2 Int_0^c f   (finite c only)
 *
 * Correctness hinges on two things: (1) the parity test is a *proof*
 * (Simplify[f +/- f(-x)] === 0, not the heuristic PossibleZeroQ, which mistakes
 * a decaying trend for a zero -- see project_possiblezeroq_decay_false_positive),
 * and (2) the value is only claimed once the half integral is confirmed finite,
 * so a divergent principal value (e.g. Int x/(1+x^2) over the whole line) is
 * never reported as 0 but instead left for another mechanism.
 */

#include "integrate_symmetry.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"

#include <string.h>
#include <stdbool.h>

/* ---- small construction / evaluation helpers ---------------------------- */

static Expr* cp(const Expr* e) { return expr_copy((Expr*)e); }

static Expr* mk_fn(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
static Expr* mk_fn1(const char* head, Expr* a) {
    Expr* args[1] = { a };
    return mk_fn(head, args, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return mk_fn(head, args, 2);
}

/* Evaluate `call` (taking ownership) and return the result. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);   /* evaluate does not free its input */
    expr_free(call);
    return r;
}
/* Simplify[e], consuming e. */
static Expr* simplify_take(Expr* e) { return eval_take(mk_fn1("Simplify", e)); }

static bool head_name_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

static bool sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

/* Does `e` contain the variable `x` (a symbol) anywhere? */
static bool contains_symbol(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return x->type == EXPR_SYMBOL && e->data.symbol.name == x->data.symbol.name;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol(e->data.function.args[i], x)) return true;
    return false;
}

/* True iff `e` mentions any divergence / indeterminacy symbol. */
static bool mentions_nonfinite(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return sym_is(e, "Infinity") || sym_is(e, "ComplexInfinity") ||
               sym_is(e, "Indeterminate") || sym_is(e, "Undefined");
    if (head_name_is(e, "DirectedInfinity") || head_name_is(e, "Integrate"))
        return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (mentions_nonfinite(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (mentions_nonfinite(e->data.function.args[i])) return true;
    return false;
}

/* A value counts as a finite definite result when it is non-NULL, free of the
 * integration variable, and free of any divergence / unevaluated-Integrate
 * marker.  (A ConditionalExpression is accepted: its convergence strip is the
 * caller's concern, and doubling / zeroing it stays correct.) */
static bool is_finite_value(const Expr* v, const Expr* x) {
    return v && !contains_symbol(v, x) && !mentions_nonfinite(v);
}

/* f with x -> -x, evaluated: ReplaceAll[f, x -> -x]. */
static Expr* reflect(const Expr* f, const Expr* x) {
    Expr* negx = mk_fn2("Times", expr_new_integer(-1), cp(x));
    Expr* rule = mk_fn2("Rule", cp(x), negx);
    return eval_take(mk_fn2("ReplaceAll", cp(f), rule));
}

/* Proof that Simplify[e] === 0. */
static bool proves_zero(Expr* e) {
    Expr* s = simplify_take(e);
    bool z = s && ((s->type == EXPR_INTEGER && s->data.integer == 0) ||
                   (s->type == EXPR_REAL && s->data.real == 0.0));
    if (s) expr_free(s);
    return z;
}

/* Is [a, b] symmetric about the origin?  Either b == -a as a proven algebraic
 * identity (a + b simplifies to 0), or the two-sided infinite line
 * (a = -Infinity, b = +Infinity). */
static bool is_pos_infinity(const Expr* b) {
    if (sym_is(b, "Infinity")) return true;
    /* DirectedInfinity[1] */
    return head_name_is(b, "DirectedInfinity") &&
           b->data.function.arg_count == 1 &&
           b->data.function.args[0]->type == EXPR_INTEGER &&
           b->data.function.args[0]->data.integer == 1;
}
static bool is_neg_infinity(const Expr* a) {
    /* DirectedInfinity[-1] */
    if (head_name_is(a, "DirectedInfinity") && a->data.function.arg_count == 1 &&
        a->data.function.args[0]->type == EXPR_INTEGER &&
        a->data.function.args[0]->data.integer == -1)
        return true;
    /* -Infinity is stored as Times[-1, Infinity]. */
    if (head_name_is(a, "Times") && a->data.function.arg_count == 2) {
        Expr* p = a->data.function.args[0];
        Expr* q = a->data.function.args[1];
        if (p->type == EXPR_INTEGER && p->data.integer == -1 && is_pos_infinity(q))
            return true;
    }
    return false;
}
static bool interval_is_symmetric(const Expr* a, const Expr* b, bool* infinite) {
    *infinite = false;
    if (is_neg_infinity(a) && is_pos_infinity(b)) { *infinite = true; return true; }
    if (is_pos_infinity(a) || is_pos_infinity(b) ||
        is_neg_infinity(a) || is_neg_infinity(b)) return false;
    Expr* sum = mk_fn2("Plus", cp(a), cp(b));
    return proves_zero(sum);
}

/* Evaluate the half integral Integrate[f, {x, 0, c}] (carrying assumptions if
 * present), returning the evaluated value (possibly unevaluated). */
static Expr* half_integral(const Expr* f, const Expr* x, const Expr* c,
                           const Expr* assumptions) {
    Expr* spec = expr_new_function(
        expr_new_symbol("List"),
        (Expr*[]){ cp(x), expr_new_integer(0), cp(c) }, 3);
    Expr* call;
    if (assumptions) {
        Expr* as = mk_fn2("Rule", expr_new_symbol("Assumptions"), cp(assumptions));
        call = expr_new_function(expr_new_symbol("Integrate"),
                                 (Expr*[]){ cp(f), spec, as }, 3);
    } else {
        call = expr_new_function(expr_new_symbol("Integrate"),
                                 (Expr*[]){ cp(f), spec }, 2);
    }
    return eval_take(call);
}

/* ---- core --------------------------------------------------------------- */

Expr* integrate_symmetry_try(Expr* f, Expr* x, Expr* a, Expr* b,
                             Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!contains_symbol(f, x)) return NULL;

    bool infinite = false;
    if (!interval_is_symmetric(a, b, &infinite)) return NULL;

    Expr* fr = reflect(f, x);
    if (!fr) return NULL;

    /* Odd: f + f(-x) == 0.  Result is 0 provided the half converges. */
    Expr* odd_sum = mk_fn2("Plus", cp(f), cp(fr));
    bool is_odd = proves_zero(odd_sum);
    if (is_odd) {
        expr_free(fr);
        Expr* half = half_integral(f, x, b, assumptions);
        bool fin = is_finite_value(half, x);
        if (half) expr_free(half);
        if (fin) return expr_new_integer(0);
        return NULL;   /* half diverges / unknown: not our call to make */
    }

    /* Even: f - f(-x) == 0.  Reduce to 2 Integrate[f, {x, 0, c}].  This runs
     * only after residue has declined, so its clean closed forms are preserved;
     * here it mainly unlocks the half-line Ramanujan/Mellin path for a
     * non-rational even integrand (e.g. a Gaussian on the whole line). */
    (void)infinite;
    Expr* even_diff = mk_fn2("Plus", cp(f), mk_fn2("Times", expr_new_integer(-1), fr));
    bool is_even = proves_zero(even_diff);
    if (is_even) {
        Expr* half = half_integral(f, x, b, assumptions);
        if (is_finite_value(half, x))
            return simplify_take(mk_fn2("Times", expr_new_integer(2), half));
        if (half) expr_free(half);
        return NULL;
    }

    return NULL;
}

/* ---- builtin ------------------------------------------------------------ */

Expr* builtin_integrate_symmetry(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* f    = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!head_name_is(spec, "List") || spec->data.function.arg_count != 3)
        return NULL;
    Expr* x = spec->data.function.args[0];
    Expr* a = spec->data.function.args[1];
    Expr* b = spec->data.function.args[2];
    if (x->type != EXPR_SYMBOL) return NULL;

    Expr* assumptions = NULL;
    for (size_t t = 2; t < argc; t++) {
        Expr* opt = res->data.function.args[t];
        if (opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2 &&
            opt->data.function.head->type == EXPR_SYMBOL &&
            (opt->data.function.head->data.symbol.name == SYM_Rule ||
             opt->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
            Expr* lhs = opt->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL &&
                strcmp(lhs->data.symbol.name, "Assumptions") == 0) {
                assumptions = opt->data.function.args[1];
                continue;
            }
        }
        return NULL;
    }
    return integrate_symmetry_try(f, x, a, b, assumptions);
}

void integrate_symmetry_init(void) {
    symtab_add_builtin("Integrate`Symmetry", builtin_integrate_symmetry);
    symtab_get_def("Integrate`Symmetry")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`Symmetry",
        "Integrate`Symmetry[f, {x, a, b}] evaluates a definite integral over an "
        "interval symmetric about the origin by parity: an odd integrand "
        "(f(-x) == -f(x)) integrates to 0, and an even integrand (f(-x) == f(x)) "
        "over a finite [-c, c] reduces to 2 Integrate[f, {x, 0, c}].  The parity "
        "is proved by Simplify; the value is claimed only when the half integral "
        "converges, so a divergent principal value is never reported as 0.  "
        "Returns unevaluated when the interval is not symmetric, the integrand "
        "has no definite parity, or the half integral does not close.");
}
