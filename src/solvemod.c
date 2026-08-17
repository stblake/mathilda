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
#include <stdlib.h>

#include "eval.h"
#include "expr.h"
#include "gbmod.h"
#include "groebner.h"
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

/* ------------------------------------------------------------------ *
 *  Modular systems: Solve[{eqns}, {vars}, Modulus -> p] over GF(p).   *
 *  Requires p PRIME (GF(p) is a field, so the Buchberger engine in     *
 *  gbmod.c applies).  Composite p is refused (single-variable composite *
 *  still works via the residue-enumeration path below).                *
 * ------------------------------------------------------------------ */

/* Runaway backstop: a positive-dimensional ideal over GF(p) can have up to
 * p^dim points; cap the collected set and refuse (unevaluated) rather than
 * truncate. */
#define MODSYS_MAX_SOLS 200000

typedef struct {
    GFpPoly** gb; int ngb;
    Expr**    var_arr; int nvar;
    uint64_t  p;
    Expr**    sols; int nsol, solcap;
    uint64_t* solvals;   /* nsol * nvar, parallel to sols, for canonical sort */
    bool      overflow;
} ModSysCtx;

/* Substitute the already-assigned higher variables (indices > vi) into g. */
static GFpPoly* modsys_subst_higher(const GFpPoly* g, const uint64_t* values,
                                    int vi, int nvar) {
    GFpPoly* cur = gfp_poly_copy(g);
    for (int idx = vi + 1; idx < nvar; idx++) {
        GFpPoly* nx = gfp_poly_subst(cur, idx, values[idx]);
        gfp_poly_free(cur);
        cur = nx;
    }
    return cur;
}

static void modsys_emit(ModSysCtx* c, const uint64_t* values) {
    if (c->nsol == c->solcap) {
        c->solcap = c->solcap ? c->solcap * 2 : 16;
        c->sols = (Expr**)realloc(c->sols, sizeof(Expr*) * (size_t)c->solcap);
        c->solvals = (uint64_t*)realloc(c->solvals,
            sizeof(uint64_t) * (size_t)c->solcap * (size_t)c->nvar);
    }
    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (size_t)c->nvar);
    for (int i = 0; i < c->nvar; i++) {
        rules[i] = mk_rule(expr_copy(c->var_arr[i]), mk_int((int64_t)values[i]));
        c->solvals[(size_t)c->nsol * (size_t)c->nvar + (size_t)i] = values[i];
    }
    c->sols[c->nsol++] = mk_list(rules, (size_t)c->nvar);
    free(rules);
}

/* Triangular back-substitution over GF(p): solve the last variable first
 * (its lex GB generator is univariate), branch over each residue root, and
 * recurse.  A variable with no univariate generator at its level is free and
 * ranges over all of GF(p); the leaf verifies every generator vanishes. */
static void modsys_rec(ModSysCtx* c, uint64_t* values, int depth) {
    if (c->overflow) return;
    if (c->nsol > MODSYS_MAX_SOLS) { c->overflow = true; return; }
    if (depth == c->nvar) {
        for (int g = 0; g < c->ngb; g++)
            if (gfp_eval_full(c->gb[g], values) != 0) return;
        modsys_emit(c, values);
        return;
    }
    int vi = c->nvar - 1 - depth;
    GFpPoly* best = NULL; int bestdeg = 0;
    for (int g = 0; g < c->ngb; g++) {
        GFpPoly* s = modsys_subst_higher(c->gb[g], values, vi, c->nvar);
        if (gfp_univariate_in(s, vi) && gfp_degree_in(s, vi) >= 1) {
            int dg = gfp_degree_in(s, vi);
            if (!best || dg < bestdeg) { if (best) gfp_poly_free(best); best = s; bestdeg = dg; continue; }
        }
        gfp_poly_free(s);
    }
    if (best) {
        for (uint64_t r = 0; r < c->p && !c->overflow; r++)
            if (gfp_eval_univariate(best, vi, r) == 0) {
                values[vi] = r;
                modsys_rec(c, values, depth + 1);
            }
        gfp_poly_free(best);
    } else {
        for (uint64_t r = 0; r < c->p && !c->overflow; r++) {   /* free variable */
            values[vi] = r;
            modsys_rec(c, values, depth + 1);
        }
    }
}

/* Canonical ascending order by the value vector (x0, x1, ...). */
static int modsys_val_cmp(ModSysCtx* c, int a, int b) {
    for (int i = 0; i < c->nvar; i++) {
        uint64_t va = c->solvals[(size_t)a * (size_t)c->nvar + (size_t)i];
        uint64_t vb = c->solvals[(size_t)b * (size_t)c->nvar + (size_t)i];
        if (va != vb) return va < vb ? -1 : 1;
    }
    return 0;
}

/* Solve a polynomial system over GF(p).  `equations` is List/And of Equal (or
 * a single Equal); `vars` is a List of >= 1 distinct symbols; p is prime.
 * Returns {{v -> r, ...}, ...} (r in [0, p)), {} if inconsistent, or NULL to
 * refuse (non-polynomial, denominator pole mod p, or overflow). */
static Expr* solvemod_solve_system(Expr* equations, Expr* vars, int64_t p) {
    if (!gfp_is_prime((uint64_t)p)) return NULL;   /* GF(p) needs a prime field */

    if (!vars || vars->type != EXPR_FUNCTION
        || vars->data.function.head->type != EXPR_SYMBOL
        || vars->data.function.head->data.symbol.name != SYM_List) return NULL;
    int n = (int)vars->data.function.arg_count;
    if (n < 1) return NULL;
    Expr** var_arr = vars->data.function.args;
    for (int j = 0; j < n; j++)
        if (var_arr[j]->type != EXPR_SYMBOL) return NULL;
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++)
            if (var_arr[a]->data.symbol.name == var_arr[b]->data.symbol.name) return NULL;

    /* Equation array. */
    Expr* single_holder[1];
    Expr** eq_arr; int m;
    if (equations->type == EXPR_FUNCTION
        && equations->data.function.head->type == EXPR_SYMBOL
        && (equations->data.function.head->data.symbol.name == SYM_List
            || equations->data.function.head->data.symbol.name == SYM_And)) {
        eq_arr = equations->data.function.args;
        m = (int)equations->data.function.arg_count;
    } else if (equations->type == EXPR_FUNCTION
        && equations->data.function.head->type == EXPR_SYMBOL
        && equations->data.function.head->data.symbol.name == SYM_Equal) {
        single_holder[0] = equations; eq_arr = single_holder; m = 1;
    } else {
        return NULL;
    }

    /* Build the GF(p) polynomial system: residual = Expand(lhs - rhs), to a
     * rational GBPoly (gb_from_expr), reduced into GF(p) (gfp_from_gbpoly). */
    GFpPoly** F = (GFpPoly**)malloc(sizeof(GFpPoly*) * (size_t)(m > 0 ? m : 1));
    int nF = 0; bool ok = true;
    for (int i = 0; i < m && ok; i++) {
        Expr* eq = eq_arr[i];
        if (eq->type != EXPR_FUNCTION
            || eq->data.function.head->type != EXPR_SYMBOL
            || eq->data.function.head->data.symbol.name != SYM_Equal
            || eq->data.function.arg_count != 2) { ok = false; break; }
        Expr* diff = mk_fn2("Plus", expr_copy(eq->data.function.args[0]),
                            mk_fn2("Times", mk_int(-1), expr_copy(eq->data.function.args[1])));
        Expr* resid = eval_and_free(internal_expand((Expr*[]){ diff }, 1));
        GBPoly* g = gb_from_expr(resid, var_arr, n, GB_ORDER_LEX, 0, NULL);
        expr_free(resid);
        if (!g) { ok = false; break; }                    /* non-polynomial */
        GFpPoly* gp = gfp_from_gbpoly(g, GB_ORDER_LEX, (uint64_t)p);
        gb_poly_free(g);
        if (!gp) { ok = false; break; }                   /* denominator pole mod p */
        if (!gfp_poly_is_zero(gp)) F[nF++] = gp; else gfp_poly_free(gp);
    }
    if (!ok) {
        for (int i = 0; i < nF; i++) gfp_poly_free(F[i]);
        free(F);
        return NULL;
    }

    /* All residuals vanished mod p -> every point is a solution: enumerate
     * GF(p)^n (guarded by the overflow cap). */
    size_t nG = 0;
    GFpPoly** G;
    if (nF == 0) { G = NULL; nG = 0; }
    else G = gfp_buchberger(F, (size_t)nF, &nG);
    for (int i = 0; i < nF; i++) gfp_poly_free(F[i]);
    free(F);

    /* Unit ideal (a nonzero constant generator) -> no solutions -> {}. */
    for (size_t i = 0; i < nG; i++) {
        if (gfp_poly_is_constant(G[i]) && !gfp_poly_is_zero(G[i])) {
            gfp_basis_free(G, nG);
            return mk_list(NULL, 0);
        }
    }

    ModSysCtx c;
    c.gb = G; c.ngb = (int)nG;
    c.var_arr = var_arr; c.nvar = n;
    c.p = (uint64_t)p;
    c.sols = NULL; c.nsol = 0; c.solcap = 0; c.solvals = NULL; c.overflow = false;

    uint64_t* values = (uint64_t*)calloc((size_t)n, sizeof(uint64_t));
    modsys_rec(&c, values, 0);
    free(values);
    gfp_basis_free(G, nG);

    if (c.overflow) {
        for (int i = 0; i < c.nsol; i++) expr_free(c.sols[i]);
        free(c.sols); free(c.solvals);
        return NULL;
    }

    /* Sort tuples into canonical ascending order (selection sort over a small
     * index set; solution counts here are modest). */
    int* order = (int*)malloc(sizeof(int) * (size_t)(c.nsol > 0 ? c.nsol : 1));
    for (int i = 0; i < c.nsol; i++) order[i] = i;
    for (int i = 0; i + 1 < c.nsol; i++) {
        int pick = i;
        for (int j = i + 1; j < c.nsol; j++)
            if (modsys_val_cmp(&c, order[j], order[pick]) < 0) pick = j;
        if (pick != i) { int t = order[i]; order[i] = order[pick]; order[pick] = t; }
    }
    Expr** sorted = (Expr**)malloc(sizeof(Expr*) * (size_t)(c.nsol > 0 ? c.nsol : 1));
    for (int i = 0; i < c.nsol; i++) sorted[i] = c.sols[order[i]];
    free(order);
    Expr* result = mk_list(sorted, (size_t)c.nsol);
    free(sorted);
    free(c.sols);
    free(c.solvals);
    return result;
}

Expr* solvemod_solve_modular(Expr* equation, Expr* vars, Expr* dom,
                             Expr* modulus) {
    (void)dom;  /* Modulus overrides any domain. */
    if (!equation || !vars) return NULL;

    int64_t p = extract_modulus(modulus);
    if (p == 0) return NULL;                 /* bad / out-of-range modulus */

    Expr* var = single_var(vars);
    if (!var) return solvemod_solve_system(equation, vars, p);  /* system over GF(p) */

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
