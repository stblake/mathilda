/*
 * dsolve_linearizable.c — DSolve`Linearizable.
 *
 * First-order ODEs that become LINEAR or BERNOULLI in a new variable u = phi(y)
 * for an elementary phi.  Given y' == F(x,y), the substitution u = phi(y) gives
 *     u' = phi'(y) F(x,y),
 * and for the right phi the right-hand side, re-expressed through y = psi(u)
 * (psi = phi^{-1}), collapses to a function G(x,u) that is linear or Bernoulli
 * in u — which the scalar cascade already solves.  The recovered u = H(x,C) is
 * mapped back by y = psi(H).
 *
 * Candidate table (phi, psi):
 *     Log/Exp   — u = Log[y]  (y = Exp[u]);  u = Exp[y] (y = Log[u])
 *     Sin/Cos   — u = Sin[y]  (y = ArcSin[u]);  u = Cos[y] (y = ArcCos[u])
 *     Tan       — u = Tan[y]  (y = ArcTan[u])
 * The Sin/Cos/Tan cases introduce a Sqrt[1-u^2] / Sqrt[1+u^2] from
 * Cos[ArcSin[u]] etc.; the method keeps a candidate only when that radical (and
 * every inverse-function artifact) cancels, leaving G a clean rational/
 * elementary function of x and u.  Examples that reduce cleanly:
 *     y' Cos[y] = Cos[x] Sin[y]^2 + Sin[y]           -> u=Sin[y]: u' = u + Cos[x] u^2  (Bernoulli)
 *     y' - Tan[y]/(x+1) = (x+1)e^x Sec[y]            -> u=Sin[y]: u' - u/(x+1) = (x+1)e^x  (linear)
 *     y' = e^{x-y}(e^x - e^y)                        -> u=Exp[y]: u' + e^x u = e^{2x}  (linear)
 *     y' = -(x Log[y]+Log[y]-1)y/(x+1)               -> u=Log[y]: u' + u = 1/(x+1)  (linear)
 *
 * Gate: F must contain a transcendental function of y (Sin/Cos/.../Log/Exp of a
 * y-argument).  This both selects the intended equations and prevents recursion
 * — the reduced equation is rational in u, so it can never re-enter this method.
 * Runs after the y-linear / y-Bernoulli specialists and before the Lie backstop.
 *
 * Verification is the substrate's (dsolve_run back-substitutes y = psi(H) into
 * the original ODE), so a spurious candidate can only be dropped.
 */
#include "dsolve_common.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../eval.h"
#include <stdlib.h>
#include <string.h>

/* Does `e` apply a transcendental function to an argument containing `yname`? */
static bool lz_transc_of_y(const Expr* e, const char* yname) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* n = h->data.symbol.name;
        static const char* const T[] = { "Sin","Cos","Tan","Cot","Sec","Csc",
                                          "Sinh","Cosh","Tanh","Log", NULL };
        for (int i = 0; T[i]; i++) if (strcmp(n, T[i]) == 0) {
            for (size_t k = 0; k < e->data.function.arg_count; k++)
                if (ds_contains(e->data.function.args[k], yname)) return true;
        }
        /* Exp[...] is Power[E, arg]; catch a y in the exponent */
        if (n == SYM_Power && e->data.function.arg_count == 2) {
            const Expr* base = e->data.function.args[0];
            if (base->type == EXPR_SYMBOL && strcmp(base->data.symbol.name, "E") == 0
                && ds_contains(e->data.function.args[1], yname)) return true;
        }
    }
    if (lz_transc_of_y(h, yname)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (lz_transc_of_y(e->data.function.args[i], yname)) return true;
    return false;
}

/* A non-integer power of a usym-containing base — an uncancelled inverse-function
 * radical (Cos[ArcSin[u]] = Sqrt[1-u^2]).  A genuine Sqrt like this means the
 * substitution did not linearise, and Simplify cannot remove it, so it is a
 * cheap definitive reject. */
static bool lz_has_frac_power_of_u(const Expr* e, const char* usym) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power
        && e->data.function.arg_count == 2
        && e->data.function.args[1]->type != EXPR_INTEGER
        && ds_contains(e->data.function.args[0], usym))
        return true;
    if (lz_has_frac_power_of_u(h, usym)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (lz_has_frac_power_of_u(e->data.function.args[i], usym)) return true;
    return false;
}

/* An inverse-trig function (ArcSin/ArcCos/ArcTan) of a usym-containing argument
 * that plain evaluation left uncollapsed (e.g. Cos[2 ArcCos[u]]); a full Simplify
 * may still collapse it (to 2u^2-1), so this triggers the expensive pass rather
 * than a reject. */
static bool lz_has_arc_of_u(const Expr* e, const char* usym) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* n = h->data.symbol.name;
        if ((strcmp(n, "ArcSin") == 0 || strcmp(n, "ArcCos") == 0 ||
             strcmp(n, "ArcTan") == 0))
            for (size_t k = 0; k < e->data.function.arg_count; k++)
                if (ds_contains(e->data.function.args[k], usym)) return true;
    }
    if (lz_has_arc_of_u(h, usym)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (lz_has_arc_of_u(e->data.function.args[i], usym)) return true;
    return false;
}

/* Re-entry guard: the reduced equation is rational in u, but keep a hard stop so
 * a candidate that slips through the artifact check can never recurse into this
 * method (or spin the heuristic backstop) on its own sub-solve. */
static int g_lz_active = 0;

/* Cheap applicability pre-filter: does `e` contain the kernel a candidate needs,
 * applied to the dependent symbol `y`?  kind 0 = Log[y], 1 = Exp[y]
 * (Power[E, y]), 2 = a circular-trig function of y.  Skipping the others avoids a
 * costly Simplify on a substitution that cannot linearise (e.g. Log on a trig
 * ODE builds Sin[Exp[u]], on which Simplify spins). */
static bool lz_has_kernel(const Expr* e, int kind, const char* y) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* n = h->data.symbol.name;
        if (kind == 0 && strcmp(n, "Log") == 0 && ds_contains(e->data.function.args[0], y))
            return true;
        if (kind == 1 && n == SYM_Power && e->data.function.arg_count == 2) {
            const Expr* b = e->data.function.args[0];
            if (b->type == EXPR_SYMBOL && strcmp(b->data.symbol.name, "E") == 0
                && ds_contains(e->data.function.args[1], y)) return true;
        }
        if (kind == 2) {
            static const char* const TR[] = { "Sin","Cos","Tan","Cot","Sec","Csc", NULL };
            for (int i = 0; TR[i]; i++) if (strcmp(n, TR[i]) == 0
                && ds_contains(e->data.function.args[0], y)) return true;
        }
    }
    if (lz_has_kernel(h, kind, y)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (lz_has_kernel(e->data.function.args[i], kind, y)) return true;
    return false;
}

/* Candidate c -> kernel kind for the pre-filter. */
static int lz_cand_kind(int c) { return (c == 0) ? 0 : (c == 1) ? 1 : 2; }

Expr** dsolve_linearizable_try(DSolveProblem* P, size_t* nbranch) {
    if (P->nfun != 1 || P->neq != 1 || P->max_order[0] != 1) return NULL;
    if (g_lz_active) return NULL;
    const char* yname = P->fun_names[0];
    const char* xvar  = P->ind_names[0];

    Expr* F = dsolve_solve_top_derivative(P, 1);
    if (!F) return NULL;

    /* Gate: only genuinely transcendental-in-y right-hand sides. */
    if (!lz_transc_of_y(F, yname)) { expr_free(F); return NULL; }

    const char* Ysym = intern_symbol("DSolve`lzY");
    const char* usym = intern_symbol("DSolve`lzU");
    const char* ufun = intern_symbol("DSolve`lzUf");

    /* F with y[x] -> the plain symbol Ysym */
    Expr* Fy = ds_subst(F, ds_make_funcapp(yname, 0, xvar), expr_new_symbol(Ysym));

    static const char* const PHI[] = { "Log", "Exp", "Sin", "Cos", "Tan", NULL };
    static const char* const PSI[] = { "Exp", "Log", "ArcSin", "ArcCos", "ArcTan", NULL };

    Expr** out = NULL;
    for (int c = 0; PHI[c] && !out; c++) {
        if (!lz_has_kernel(Fy, lz_cand_kind(c), Ysym)) continue;   /* pre-filter */
        /* dphi = d/dY phi[Y] */
        Expr* phiY = ds_call1(PHI[c], expr_new_symbol(Ysym));
        Expr* dphi = ds_d(phiY, expr_new_symbol(Ysym));
        /* Gy = dphi * Fy, simplified in Y (canonicalises the trig product) */
        Expr* Gy = ds_simplify(eval_and_free(ds_call2(SYM_Times, dphi, expr_copy(Fy))));
        /* substitute Y -> psi[u], then just EVALUATE (cheap): this collapses the
         * inverse compositions (Sin[ArcSin u] -> u, Cos[ArcSin u] -> Sqrt[1-u^2]),
         * so a substitution that does not linearise shows a genuine Sqrt of u and
         * is rejected WITHOUT paying for a full Simplify. */
        Expr* psiU = ds_call1(PSI[c], expr_new_symbol(usym));
        Expr* G = eval_and_free(ds_subst(Gy, expr_new_symbol(Ysym), psiU));  /* consumes Gy */

        if (lz_has_frac_power_of_u(G, usym)) { expr_free(G); continue; }     /* genuine Sqrt: reject */
        if (lz_has_arc_of_u(G, usym)) {              /* inverse left uncollapsed: full Simplify may fix */
            G = ds_simplify(G);
            if (lz_has_frac_power_of_u(G, usym) || lz_has_arc_of_u(G, usym)) { expr_free(G); continue; }
        }
        if (!ds_free_of(G, Ysym)) { expr_free(G); continue; }

        /* Solve u' == G(x,u).  When G is LINEAR in u — u' == A(x) u + B(x) —
         * integrate directly (fast, no recursion).  Otherwise (Bernoulli, …)
         * recurse into the scalar cascade under the re-entry guard. */
        Expr* dGu = ds_d(expr_copy(G), expr_new_symbol(usym));       /* dG/du */
        Expr** ubodies = NULL;
        size_t nb = 0;
        if (ds_free_of(dGu, usym)) {                                 /* linear */
            Expr* B  = ds_subst(expr_copy(G), expr_new_symbol(usym), expr_new_integer(0));
            Expr* Pc = eval_and_free(ds_call2(SYM_Times, expr_new_integer(-1), expr_copy(dGu)));
            Expr* ub = dsolve_linear_factor_solve(Pc, B, xvar);      /* consumes Pc, B */
            if (ub) { ubodies = malloc(sizeof(Expr*)); ubodies[0] = ub; nb = 1; }
        } else {
            Expr* Gu = ds_subst(expr_copy(G), expr_new_symbol(usym),
                                ds_make_funcapp(ufun, 0, xvar));
            Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
                            (Expr*[]){ ds_make_funcapp(ufun, 1, xvar), Gu }, 2);
            Expr* call = expr_new_function(expr_new_symbol(SYM_DSolve),
                            (Expr*[]){ eqn, ds_make_funcapp(ufun, 0, xvar),
                                       expr_new_symbol(xvar) }, 3);
            g_lz_active++;
            Expr* r = eval_and_free(call);
            g_lz_active--;
            ubodies = dsolve_extract_applied_bodies(r, ufun, &nb);
            expr_free(r);
        }
        expr_free(dGu);
        expr_free(G);
        if (!ubodies || nb == 0) { if (ubodies) free(ubodies); continue; }

        /* y = psi(H) for each recovered branch H (u -> H) */
        Expr** yb = malloc(nb * sizeof(Expr*));
        for (size_t k = 0; k < nb; k++)
            yb[k] = eval_and_free(ds_call1(PSI[c], ubodies[k]));   /* consumes ubodies[k] */
        free(ubodies);
        out = yb;
        *nbranch = nb;
    }

    expr_free(Fy);
    return out;
}

static Expr* builtin_dsolve_linearizable(Expr* res) {
    return dsolve_method_builtin(res, dsolve_linearizable_try);
}

void dsolve_linearizable_init(void) {
    symtab_add_builtin("DSolve`Linearizable", builtin_dsolve_linearizable);
    symtab_get_def("DSolve`Linearizable")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DSolve`Linearizable",
        "DSolve`Linearizable[eqn, y, x] solves a first-order ODE that becomes "
        "linear or Bernoulli under a substitution u = phi(y) (phi in "
        "Log/Exp/Sin/Cos/Tan): it reduces to u' == G(x,u), solves that with the "
        "scalar cascade, and maps back y = phi^{-1}(u).");
}
