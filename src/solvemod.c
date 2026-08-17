/*
 * solvemod.c
 *
 * Modular solving pre-pass for Solve[poly == 0, x, Modulus -> p].
 * See solvemod.h for the contract and the refusal policy.
 */

#include "solvemod.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "poly.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ *
 *  Tiny construction helpers (mirror the sibling solve specialists). *
 * ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* mk_rule(Expr* lhs, Expr* rhs) { return mk_fn2("Rule", lhs, rhs); }
static Expr* mk_list(Expr** args, size_t n) {
    return expr_new_function(mk_sym("List"), args, n);
}

/* Extract a usable small modulus from `modulus`.  Returns the value in
 * [2, SOLVEMOD_P_MAX] on success, or 0 to refuse (symbolic, non-integer,
 * BigInt, or out of range). */
static int64_t extract_modulus(const Expr* modulus) {
    if (!modulus || modulus->type != EXPR_INTEGER) return 0;
    int64_t p = modulus->data.integer;
    if (p < 2 || p > SOLVEMOD_P_MAX) return 0;
    return p;
}

/* Resolve the single solve variable.  `vars` is the symbol-only spec:
 * a bare symbol, or a List of exactly one symbol.  Returns the symbol
 * (borrowed) or NULL to refuse (system / multivariable / malformed). */
static Expr* single_var(Expr* vars) {
    if (!vars) return NULL;
    if (vars->type == EXPR_SYMBOL) return vars;
    if (vars->type == EXPR_FUNCTION
        && vars->data.function.head->type == EXPR_SYMBOL
        && vars->data.function.head->data.symbol.name == SYM_List
        && vars->data.function.arg_count == 1
        && vars->data.function.args[0]->type == EXPR_SYMBOL) {
        return vars->data.function.args[0];
    }
    return NULL;
}

/* Is `residual mod p` the integer 0 when `var -> r` is substituted?
 * internal_mod handles rational coefficients (via modular inverse), so a
 * residual like x/3 - 1 is reduced correctly. */
static bool root_mod_p(Expr* residual, Expr* var, int64_t r, int64_t p) {
    Expr* rule = mk_rule(expr_copy(var), mk_int(r));
    Expr* sub = eval_and_free(internal_replace_all(
        (Expr*[]){ expr_copy(residual), rule }, 2));
    Expr* modded = eval_and_free(internal_mod(
        (Expr*[]){ sub, mk_int(p) }, 2));
    bool is_zero = (modded && modded->type == EXPR_INTEGER
                    && modded->data.integer == 0);
    if (modded) expr_free(modded);
    return is_zero;
}

Expr* solvemod_solve_modular(Expr* equation, Expr* vars, Expr* dom,
                             Expr* modulus) {
    (void)dom;  /* Modulus overrides any domain. */
    if (!equation || !vars) return NULL;

    int64_t p = extract_modulus(modulus);
    if (p == 0) return NULL;                 /* bad / out-of-range modulus */

    Expr* var = single_var(vars);
    if (!var) return NULL;                   /* system / multivariable */

    /* Require an explicit equation lhs == rhs. */
    if (equation->type != EXPR_FUNCTION
        || equation->data.function.head->type != EXPR_SYMBOL
        || equation->data.function.head->data.symbol.name != SYM_Equal
        || equation->data.function.arg_count != 2) {
        return NULL;
    }
    Expr* lhs = equation->data.function.args[0];
    Expr* rhs = equation->data.function.args[1];

    /* residual = lhs - rhs. */
    Expr* residual = eval_and_free(mk_fn2("Plus",
        expr_copy(lhs),
        mk_fn2("Times", mk_int(-1), expr_copy(rhs))));

    /* Only polynomial equations are solvable by residue enumeration. */
    if (!is_polynomial(residual, &var, 1)) {
        expr_free(residual);
        return NULL;
    }

    /* Constant residual (no var): 0 mod p -> tautology {{}}, else {}. */
    if (!contains_any_symbol_from(residual, var)) {
        Expr* modded = eval_and_free(internal_mod(
            (Expr*[]){ expr_copy(residual), mk_int(p) }, 2));
        bool zero = (modded && modded->type == EXPR_INTEGER
                     && modded->data.integer == 0);
        if (modded) expr_free(modded);
        expr_free(residual);
        if (zero) {
            Expr* empty = mk_list(NULL, 0);
            return mk_list((Expr*[]){ empty }, 1);
        }
        return mk_list(NULL, 0);
    }

    /* Enumerate residues 0 .. p-1 and collect the roots (ascending). */
    Expr** sols = (Expr**)malloc(sizeof(Expr*) * (size_t)p);
    size_t nsol = 0;
    for (int64_t r = 0; r < p; r++) {
        if (root_mod_p(residual, var, r, p)) {
            Expr* rule = mk_rule(expr_copy(var), mk_int(r));
            sols[nsol++] = mk_list((Expr*[]){ rule }, 1);
        }
    }
    expr_free(residual);

    Expr* result = mk_list(sols, nsol);
    free(sols);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Qualified-builtin entry: Solve`SolveModular[eqn, var, p]          *
 * ------------------------------------------------------------------ */

static Expr* builtin_solve_modular(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 3) return NULL;
    Expr* equation = res->data.function.args[0];
    Expr* var      = res->data.function.args[1];
    Expr* modulus  = res->data.function.args[2];
    return solvemod_solve_modular(equation, var, NULL, modulus);
}

void solvemod_init(void) {
    symtab_add_builtin("Solve`SolveModular", builtin_solve_modular);
    symtab_set_docstring("Solve`SolveModular",
        "Solve`SolveModular[eqn, var, p]\n"
        "\tInternal: solves a single-variable polynomial equation eqn\n"
        "\tfor var over Z/pZ by residue enumeration, returning\n"
        "\t{{var -> r}, ...} with r ascending in [0, p).  Refuses\n"
        "\t(returns unevaluated) for systems, non-polynomial residuals,\n"
        "\tor a modulus outside [2, 100000].");
}
