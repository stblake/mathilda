/*
 * dsolve_triangular.c — DSolve`TriangularSystem.
 *
 * Solves a system whose inter-function dependency graph is a DAG (a
 * lower-triangular / cascade system) by forward substitution.  Repeatedly find
 * an as-yet-unused equation that mentions exactly one still-unsolved function,
 * solve it for that function with the scalar engine, substitute the solution
 * (and its derivatives) forward into the remaining equations, and renumber the
 * generated constants so they do not collide.  This generalizes
 * DSolve`DecoupleSystem (which needs each equation to mention exactly one
 * function) to genuinely coupled-but-triangular systems, at ANY coefficient —
 * constant OR variable — that the constant-A matrix exponential cannot reach:
 *
 *     {y'==0, x'+y==0}       -> y=C[1], x=C[2]-C[1] x   (constant, triangular)
 *     {y'==x^2 y, x'==y}     -> variable-coefficient triangular
 *
 * When no equation isolates a single unsolved function the system is not
 * triangular and the method declines (returns NULL), leaving the later
 * LinearFirstOrderSystem (matrix exponential) to try irreducibly-coupled
 * constant systems.  Verification / condition-fitting / assembly are handled by
 * the shared substrate (dsolve_run_system).
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../eval.h"
#include "../internal.h"
#include <stdlib.h>

/* A solved function's arbitrary constants must NOT sit in the C[] namespace
 * while later functions are being solved: the scalar engine always emits its
 * own integration constants as C[1], C[2], ..., so a solved C[1] substituted
 * forward would merge with the next solve's fresh C[1].  We therefore park each
 * solved body's constants in the private head KHEAD, disjoint from C[], and
 * remap KHEAD[k] -> C[k] once at the very end. */
#define KHEAD "DSolve`sysK"

/* Rename src[j] -> dst[base+j] for j = 1..m in `body`; body consumed, owned. */
static Expr* rename_params(Expr* body, const char* src, const char* dst, int base, int m) {
    if (m <= 0) return body;
    Expr** rules = malloc((size_t)m * sizeof(Expr*));
    for (int j = 1; j <= m; j++) {
        Expr* from = expr_new_function(expr_new_symbol(src),
                         (Expr*[]){ expr_new_integer(j) }, 1);
        Expr* to   = expr_new_function(expr_new_symbol(dst),
                         (Expr*[]){ expr_new_integer(base + j) }, 1);
        rules[j - 1] = expr_new_function(expr_new_symbol(SYM_Rule),
                           (Expr*[]){ from, to }, 2);
    }
    Expr* rl = expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)m);
    free(rules);
    return eval_and_free(internal_replace_all((Expr*[]){ body, rl }, 2));
}

/* Substitute the solved function `fname -> body` (and Derivative[m][fname] ->
 * D[body,{x,m}] for m = 1..maxord) into `residual`.  `residual` is consumed and
 * `body` borrowed; the evaluated result is returned owned. */
static Expr* subst_solution(Expr* residual, const char* fname, const Expr* body,
                            int maxord, const char* xvar) {
    for (int m = maxord; m >= 1; m--) {
        Expr* dm = expr_copy((Expr*)body);
        for (int k = 0; k < m; k++) dm = ds_d(dm, expr_new_symbol(xvar));
        residual = ds_subst(residual, ds_make_funcapp(fname, m, xvar), dm);
    }
    return ds_subst(residual, ds_make_funcapp(fname, 0, xvar), expr_copy((Expr*)body));
}

Expr** dsolve_triangular_solve(DSolveProblem* P) {
    size_t n = P->nfun;
    const char* xvar = P->ind_names[0];
    if (P->neq != n) return NULL;                 /* square, one governing eqn each */

    Expr** work    = malloc(n * sizeof(Expr*));   /* residuals, solved funcs peeled */
    for (size_t e = 0; e < n; e++) work[e] = expr_copy(P->eq_residuals[e]);
    bool*  eq_used = calloc(n, sizeof(bool));
    bool*  solved  = calloc(n, sizeof(bool));
    Expr** bodies  = calloc(n, sizeof(Expr*));
    int    offset  = 0;
    size_t done    = 0;

    bool progress = true;
    while (done < n && progress) {
        progress = false;
        for (size_t e = 0; e < n && !progress; e++) {
            if (eq_used[e]) continue;
            /* count the unsolved functions this equation still mentions */
            long which = -1; int cnt = 0;
            for (size_t j = 0; j < n; j++) {
                if (solved[j]) continue;
                if (ds_contains(work[e], P->fun_names[j])) { which = (long)j; cnt++; }
            }
            if (cnt != 1) continue;
            size_t j = (size_t)which;

            /* solve  work[e] == 0  for function j via the scalar engine */
            Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                           (Expr*[]){ expr_copy(work[e]), expr_new_integer(0) }, 2);
            Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                             (Expr*[]){ eq, expr_new_symbol(P->fun_names[j]),
                                        expr_new_symbol(xvar) }, 3);
            Expr* r = eval_and_free(call);
            Expr* body = dsolve_extract_system_body(r, P->fun_names[j]);
            expr_free(r);
            if (!body) continue;                  /* try another candidate equation */

            /* park this solution's fresh C[1..m] in KHEAD[offset+1..], disjoint
             * from the C[] the next scalar solve will generate */
            bodies[j] = rename_params(body, "C", KHEAD, offset, P->max_order[j]);
            offset += P->max_order[j];
            solved[j] = true; eq_used[e] = true; done++; progress = true;

            /* substitute the solved function forward into the remaining equations */
            for (size_t f = 0; f < n; f++) {
                if (eq_used[f]) continue;
                work[f] = subst_solution(work[f], P->fun_names[j], bodies[j],
                                         P->max_order[j], xvar);
            }
        }
    }

    bool ok = (done == n);
    for (size_t e = 0; e < n; e++) expr_free(work[e]);
    free(work); free(eq_used); free(solved);
    if (!ok) {
        for (size_t j = 0; j < n; j++) if (bodies[j]) expr_free(bodies[j]);
        free(bodies);
        return NULL;
    }
    /* remap the private constants back to a contiguous C[1..offset] */
    for (size_t j = 0; j < n; j++)
        bodies[j] = rename_params(bodies[j], KHEAD, "C", 0, offset);
    return bodies;
}

void dsolve_triangular_init(void) {
    /* dispatched directly for nfun>1; no backtick builtin */
}
