/*
 * solve.c
 *
 * The `Solve` router: classifies the input equation system, parses
 * options, and dispatches to a specialist solver.  The only specialist
 * wired up in this initial cut is Solve`SolvePolynomialEquality
 * (src/solvepoly.c) for a single polynomial equality in one variable.
 *
 * `Solve` has the HoldAll attribute -- the user's `vars` argument
 * (typically a bare symbol) must reach the router without OwnValue
 * substitution.  The `expr` argument is evaluated explicitly inside
 * the builtin to normalise nested calls (e.g. (x-1)(x-2) gets
 * canonicalised before we try to recognise the polynomial shape).
 */

#include "solve.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "solvepoly.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ *
 *  Option parsing.                                                    *
 * ------------------------------------------------------------------ */

typedef struct {
    SolvePolyOpts poly;
    Expr* dom;             /* borrowed; default = NULL ( = Complexes) */
} SolveOpts;

/* Recognised Solve option-name symbols. */
static bool is_known_option_name(const char* s) {
    return s == SYM_Cubics
        || s == SYM_Quartics
        || s == SYM_GeneratedParameters
        || s == SYM_VerifySolutions
        || s == SYM_Assumptions
        || s == SYM_InverseFunctions
        || s == SYM_Method
        || s == SYM_Modulus;
}

/* True iff `e` is Rule[opt, _] or RuleDelayed[opt, _] for a recognised
 * Solve option name.  Used to peel options off the end of the arg list. */
static bool is_option_arg(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    const Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL) return false;
    return is_known_option_name(lhs->data.symbol);
}

/* Returns true iff `e` is the symbol True. */
static bool is_true(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_True;
}

/* Apply a single option rule to `opts`.  Unknown values do not abort
 * (they are silently ignored for now) -- only unknown option *names*
 * are rejected, by is_known_option_name. */
static void apply_option(const Expr* rule, SolveOpts* opts) {
    const Expr* lhs = rule->data.function.args[0];
    const Expr* rhs = rule->data.function.args[1];
    const char* name = lhs->data.symbol;
    if (name == SYM_Cubics)   { opts->poly.cubics_radical = is_true(rhs); return; }
    if (name == SYM_Quartics) { opts->poly.quartics_radical = is_true(rhs); return; }
    /* GeneratedParameters / VerifySolutions / Assumptions /
     * InverseFunctions / Method / Modulus: parsed but not yet wired
     * into the polynomial specialist. */
}

/* Warn once per distinct unevaluated form about an unrecognised
 * option.  Mirrors the integrate.c:254-262 idiom. */
static void warn_bad_option(const Expr* res, const Expr* opt) {
    static uint64_t last_warned_hash = 0;
    uint64_t h = expr_hash(res);
    if (h == last_warned_hash) return;
    last_warned_hash = h;
    const Expr* lhs = (opt && opt->type == EXPR_FUNCTION
                       && opt->data.function.arg_count == 2)
        ? opt->data.function.args[0] : NULL;
    const char* name = (lhs && lhs->type == EXPR_SYMBOL)
        ? lhs->data.symbol : "?";
    fprintf(stderr,
        "Solve::optx: Unknown option %s in Solve.\n",
        name);
}

/* ------------------------------------------------------------------ *
 *  Argument classification.                                           *
 * ------------------------------------------------------------------ */

/* Returns the single variable contained in `vars` (a symbol or a
 * length-1 List of one symbol) and writes it to *var_out.  Returns
 * false if `vars` is not a supported shape (e.g. multivariate). */
static bool classify_single_var(Expr* vars, Expr** var_out) {
    if (vars->type == EXPR_SYMBOL) { *var_out = vars; return true; }
    if (vars->type == EXPR_FUNCTION
        && vars->data.function.head->type == EXPR_SYMBOL
        && vars->data.function.head->data.symbol == SYM_List
        && vars->data.function.arg_count == 1
        && vars->data.function.args[0]->type == EXPR_SYMBOL) {
        *var_out = vars->data.function.args[0];
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Builtin entry.                                                     *
 * ------------------------------------------------------------------ */

Expr* builtin_solve(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    /* Walk trailing args, peeling options.  Position of first option
     * = end of positional args. */
    size_t pos_end = argc;
    while (pos_end > 0) {
        Expr* a = res->data.function.args[pos_end - 1];
        if (a->type == EXPR_FUNCTION
            && a->data.function.head->type == EXPR_SYMBOL
            && (a->data.function.head->data.symbol == SYM_Rule
                || a->data.function.head->data.symbol == SYM_RuleDelayed)
            && a->data.function.arg_count == 2
            && a->data.function.args[0]->type == EXPR_SYMBOL) {
            const char* name = a->data.function.args[0]->data.symbol;
            if (is_known_option_name(name)) {
                pos_end--;
                continue;
            }
            /* Trailing Rule[] that is not a recognised option name is
             * a syntax error in Solve -- bail. */
            if (!is_option_arg(a)) {
                /* shape is right (Rule[sym, _]) but name is wrong */
                warn_bad_option(res, a);
                return NULL;
            }
        }
        break;
    }

    /* Positional args: expr [, vars [, dom]] */
    if (pos_end < 2 || pos_end > 3) {
        if (pos_end < 2) return NULL;
        /* pos_end > 3: too many positional args before options.  */
        return NULL;
    }

    /* Parse options. */
    SolveOpts opts = {{ false, false }, NULL};
    for (size_t i = pos_end; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (is_option_arg(a)) apply_option(a, &opts);
    }

    /* HoldAll deferred evaluation of expr; evaluate it now (vars is
     * intentionally left as a literal symbol so that an OwnValue on
     * `x` does not get substituted). */
    Expr* expr = eval_and_free(expr_copy(res->data.function.args[0]));
    Expr* vars = res->data.function.args[1];   /* borrowed, unevaluated */
    Expr* dom  = (pos_end >= 3) ? res->data.function.args[2] : NULL;

    Expr* var = NULL;
    if (!classify_single_var(vars, &var)) {
        expr_free(expr);
        return NULL;
    }

    /* Equation may have been simplified all the way to True / False
     * during the explicit evaluate() above (e.g. Equal[1, 0] → False).
     * Mathematica's convention is:
     *   True  → {{}}   (tautology: full-dimensional solution set)
     *   False → {}     (contradiction: no solutions)            */
    Expr* out = NULL;
    if (expr->type == EXPR_SYMBOL && expr->data.symbol == SYM_True) {
        Expr* empty = expr_new_function(expr_new_symbol("List"), NULL, 0);
        out = expr_new_function(expr_new_symbol("List"),
                                (Expr*[]){ empty }, 1);
    } else if (expr->type == EXPR_SYMBOL && expr->data.symbol == SYM_False) {
        out = expr_new_function(expr_new_symbol("List"), NULL, 0);
    } else if (expr->type == EXPR_FUNCTION
        && expr->data.function.head->type == EXPR_SYMBOL
        && expr->data.function.head->data.symbol == SYM_Equal
        && expr->data.function.arg_count == 2) {
        out = solvepoly_solve_polynomial_equality(expr, var, dom, &opts.poly);
    }

    expr_free(expr);
    return out;  /* evaluator frees res on non-NULL return */
}

/* ------------------------------------------------------------------ *
 *  Init.                                                              *
 * ------------------------------------------------------------------ */

void solve_init(void) {
    symtab_add_builtin("Solve", builtin_solve);
    SymbolDef* def = symtab_get_def("Solve");
    if (def) def->attributes |= ATTR_PROTECTED | ATTR_HOLDALL;
    symtab_set_docstring("Solve",
        "Solve[expr, vars]\n"
        "\tAttempts to solve the equation or system expr for the\n"
        "\tvariables vars.\n"
        "Solve[expr, vars, dom]\n"
        "\tSolves over the domain dom.  Default Complexes; Reals filters\n"
        "\tdown to real roots via per-degree discriminant tests.\n"
        "\n"
        "Options:\n"
        "    Cubics              -> False  (radical form for cubics)\n"
        "    Quartics            -> False  (radical form for quartics)\n"
        "    GeneratedParameters -> C       (reserved)\n"
        "    VerifySolutions     -> Automatic (reserved)\n"
        "\n"
        "Initial implementation handles single polynomial equalities in\n"
        "one variable; the polynomial specialist is also directly\n"
        "callable as Solve`SolvePolynomialEquality.");

    symtab_set_docstring("Cubics",
        "Cubics is an option for Solve that controls whether cubic\n"
        "\tequations are solved via explicit radical formulas\n"
        "\t(Cubics -> True) or returned as held Root[] objects\n"
        "\t(default Cubics -> False).");
    symtab_set_docstring("Quartics",
        "Quartics is an option for Solve that controls whether quartic\n"
        "\tequations are solved via explicit radical formulas\n"
        "\t(Quartics -> True) or returned as held Root[] objects\n"
        "\t(default Quartics -> False).");
    symtab_set_docstring("GeneratedParameters",
        "GeneratedParameters is an option for Solve specifying how\n"
        "\tnewly introduced parameters should be named.  Default: C.\n"
        "\tReserved -- not yet wired into the polynomial specialist.");
    symtab_set_docstring("VerifySolutions",
        "VerifySolutions is an option for Solve that decides whether to\n"
        "\tverify each returned solution by back-substitution.\n"
        "\tDefault: Automatic.  Reserved.");

    solvepoly_init();
}
