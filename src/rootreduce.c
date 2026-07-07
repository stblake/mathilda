/* Mathilda — RootReduce implementation.
 *
 * RootReduce[expr] canonicalises an algebraic expression. It dispatches between
 * two rigorous FLINT engines depending on the shape of `expr`:
 *
 *   (1) Constant algebraic NUMBERS (no free symbol) — integers, rationals,
 *       radicals, roots of unity, the imaginary unit and Root[] objects
 *       combined by +,-,*,/,^ — are canonicalised via FLINT `qqbar`
 *       (src/poly/flint_qqbar.c) to a single representative: a rational, a
 *       quadratic radical expression, or a Root[Function[minpoly&], k] object.
 *       This is WL's central RootReduce behaviour.
 *
 *   (2) Algebraic FUNCTIONS over a tower Q(params)(radicals) — radicals whose
 *       radicand carries a free variable (e.g. the Goursat k^(1/3) towers) —
 *       are rationalised by flint_algebraic_field_canonical (src/poly/
 *       flint_bridge.c): the denominator is inverted in the field by an exact
 *       linear solve, no numeric oracle.
 *
 * RootReduce also threads over equations, inequalities and logic functions
 * (Equal, Less, And, ...), and for equations/inequalities of constant algebraic
 * numbers it decides the (in)equality exactly via `qqbar`. It is Listable, so
 * it threads over lists elementwise.
 *
 * Options: Method -> "Automatic" | "Recursive" | "NumberField" (see flint_qqbar).
 *
 * When `expr` carries no algebraic content (or the case is out of scope) it is
 * returned unchanged, matching WL. Ownership follows the builtin contract:
 * return a new tree or steal from `res`; never expr_free(res).
 */

#include "rootreduce.h"

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "symtab.h"
#include "sym_names.h"
#include "flint_bridge.h"
#include "flint_qqbar.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head && e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol, name) == 0;
}

static const char* head_name(const Expr* e) {
    return (e && e->type == EXPR_FUNCTION && e->data.function.head &&
            e->data.function.head->type == EXPR_SYMBOL)
           ? e->data.function.head->data.symbol : NULL;
}

/* An option rule Method -> value (Rule or RuleDelayed with lhs symbol Method). */
static int is_option_rule(const Expr* e) {
    if ((head_is(e, "Rule") || head_is(e, "RuleDelayed")) &&
        e->data.function.arg_count == 2) {
        const Expr* lhs = e->data.function.args[0];
        return lhs->type == EXPR_SYMBOL;   /* any symbolic-lhs rule = an option */
    }
    return 0;
}

/* Parse Method from the option rules. *bad set if an unknown value is given. */
static QQBarMethod parse_method(const Expr* res, int* bad) {
    QQBarMethod m = QQBAR_METHOD_AUTOMATIC;
    *bad = 0;
    size_t n = res->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        const Expr* a = res->data.function.args[i];
        if (!is_option_rule(a)) continue;
        const Expr* lhs = a->data.function.args[0];
        if (strcmp(lhs->data.symbol, "Method") != 0) continue;
        const Expr* v = a->data.function.args[1];
        const char* s = (v->type == EXPR_STRING) ? v->data.string
                      : (v->type == EXPR_SYMBOL) ? v->data.symbol : NULL;
        if (s && strcmp(s, "Automatic") == 0)        m = QQBAR_METHOD_AUTOMATIC;
        else if (s && strcmp(s, "Recursive") == 0)   m = QQBAR_METHOD_RECURSIVE;
        else if (s && strcmp(s, "NumberField") == 0) m = QQBAR_METHOD_NUMBERFIELD;
        else *bad = 1;
    }
    return m;
}

/* Heads RootReduce threads over. */
static int is_relational(const Expr* e) {
    const char* h = head_name(e);
    if (!h) return 0;
    static const char* heads[] = {
        "Equal", "Unequal", "Less", "LessEqual", "Greater", "GreaterEqual",
        "And", "Or", "Not", "Xor", "Implies", NULL };
    for (int i = 0; heads[i]; i++) if (strcmp(h, heads[i]) == 0) return 1;
    return 0;
}

static Expr* bool_expr(int v) { return expr_new_symbol(v ? "True" : "False"); }

/* Thread RootReduce over a relational/logical head. For a binary
 * (in)equality of constant algebraic numbers, decide it exactly via qqbar;
 * otherwise map RootReduce over the parts and re-evaluate. */
static Expr* thread_relational(const Expr* arg) {
    const char* h = head_name(arg);
    size_t n = arg->data.function.arg_count;

    if (n == 2) {
        const Expr* a = arg->data.function.args[0];
        const Expr* b = arg->data.function.args[1];
        if (strcmp(h, "Equal") == 0 || strcmp(h, "Unequal") == 0) {
            int r = flint_qqbar_equal(a, b);      /* 1 eq, 0 neq, -1 undecided */
            if (r >= 0) return bool_expr(strcmp(h, "Equal") == 0 ? (r == 1) : (r == 0));
        } else {
            int c = flint_qqbar_compare(a, b);    /* sign(a-b) or -2 undecided */
            if (c != -2) {
                int val = strcmp(h, "Less") == 0        ? (c < 0)
                        : strcmp(h, "LessEqual") == 0    ? (c <= 0)
                        : strcmp(h, "Greater") == 0      ? (c > 0)
                        : strcmp(h, "GreaterEqual") == 0 ? (c >= 0) : 0;
                return bool_expr(val);
            }
        }
    }

    /* Generic threading: head[ RootReduce[part_i], ... ], then evaluate. */
    Expr** parts = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        parts[i] = expr_new_function(expr_new_symbol("RootReduce"),
                       (Expr*[]){ expr_copy(arg->data.function.args[i]) }, 1);
    Expr* out = expr_new_function(expr_copy(arg->data.function.head), parts, n);
    free(parts);
    return eval_and_free(out);
}

/* RootReduce[expr, opts] — see file header. Not HoldAll: args are evaluated. */
Expr* builtin_rootreduce(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;

    /* Split positional args from trailing options. */
    size_t npos = 0, pos_idx = 0;
    for (size_t i = 0; i < n; i++)
        if (!is_option_rule(res->data.function.args[i])) { npos++; pos_idx = i; }

    int method_bad = 0;
    QQBarMethod method = parse_method(res, &method_bad);
    if (method_bad) {
        fprintf(stderr, "RootReduce::mtd: the Method option must be one of "
                        "\"Automatic\", \"NumberField\" or \"Recursive\".\n");
        return NULL;
    }
    if (npos != 1) {
        fprintf(stderr, "RootReduce::argx: RootReduce called with %zu argument%s; "
                        "1 argument is expected.\n", npos, npos == 1 ? "" : "s");
        return NULL;
    }

    Expr* arg = res->data.function.args[pos_idx];

    /* G4: thread over equations / inequalities / logic. */
    if (is_relational(arg)) return thread_relational(arg);

    /* G1/G2: constant algebraic number -> Root / quadratic radical / rational. */
    Expr* q = flint_qqbar_canonical(arg, method);
    if (q) return q;

    /* Parametric algebraic function: rationalise the denominator over the tower. */
#ifdef USE_FLINT
    Expr* r = flint_algebraic_field_canonical(arg);
    if (r) return r;
#endif

    /* Identity: steal the positional arg out of res (evaluator frees res). */
    res->data.function.args[pos_idx] = NULL;
    return arg;
}

void rootreduce_init(void) {
    symtab_add_builtin("RootReduce", builtin_rootreduce);
    symtab_get_def("RootReduce")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;

    /* Options[RootReduce] = {Method -> Automatic}. */
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol(SYM_Method), expr_new_symbol(SYM_Automatic) }, 2);
    Expr* opts = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ rule }, 1);
    symtab_set_options("RootReduce", opts);

    symtab_set_docstring("RootReduce",
        "RootReduce[expr] canonicalises an algebraic expression: a constant "
        "algebraic number becomes a rational, a quadratic radical, or a Root "
        "object; a rational function over a radical tower has its denominator "
        "rationalised. Threads over lists, equations, inequalities and logic. "
        "Option: Method -> \"Automatic\" | \"Recursive\" | \"NumberField\".");
}
