/* Mathilda — NDSolve Method-of-Lines front-end (PDE support, Phase 1).
 *
 * A partial differential equation in one temporal and one spatial variable is
 * solved by discretizing the spatial operator on a uniform grid: each grid
 * value U_j(t) = u(t, x_j) becomes an unknown of a large first-order ODE system
 *
 *     dU_j/dt = G(t, x_j, U, U_x|_j, U_xx|_j, ...)
 *
 * which the shared adaptive driver (nd_integrate) and every stepper already
 * solve.  Spatial derivatives are replaced by second-order central finite
 * differences; interior nodes are the unknowns and Dirichlet boundary values
 * fold into their neighbours' stencils.  The per-node right-hand side is built
 * symbolically (substituting the grid stencils into the solved-for temporal
 * derivative) and evaluated by the same Block-localized sampler used for ODEs,
 * so variable-coefficient and nonlinear PDEs work with no extra machinery.
 *
 * The accepted time nodes carry the whole spatial vector, so the solution is a
 * complete tensor grid (adaptive in t, uniform in x); it is handed to the
 * multidimensional Interpolation builtin to produce a 2-D InterpolatingFunction
 * over (t, x), applied as u[t, x].
 *
 * Phase 1 scope: one spatial dimension, one dependent function, temporal order
 * 1 or 2, spatial orders 1 and 2, Dirichlet boundary conditions, machine
 * precision.  Higher-order stencils (Fornberg), other boundary conditions, the
 * compiled linear operator, 2-D regions, MPFR and nonlinear stress cases are
 * later phases.  Nonlinear/variable-coefficient PDEs in 1-D already work here
 * through the symbolic sampler.
 */
#include "ndsolve_common.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../eval.h"
#include "../common.h"
#include "../numeric.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static void nd_mol_warn(const char* tag, const char* msg) {
    fprintf(stderr, "NDSolve::%s: %s\n", tag, msg);
}

/* ------------------------------------------------------------------ *
 *  PDE funcapp matching:  u[t,x]  or  Derivative[a,b][u][t,x]         *
 * ------------------------------------------------------------------ */
static bool nd_pde_match(const Expr* e, const char* fname,
                         int* a, int* b, const Expr** arg0, const Expr** arg1) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    const Expr* head = e->data.function.head;
    *arg0 = e->data.function.args[0];
    *arg1 = e->data.function.args[1];
    if (head->type == EXPR_SYMBOL) {
        if (head->data.symbol.name != fname) return false;
        *a = 0; *b = 0; return true;
    }
    /* Derivative[a,b][u] applied to (t,x). */
    if (head->type == EXPR_FUNCTION && head->data.function.arg_count == 1
        && head->data.function.args[0]->type == EXPR_SYMBOL
        && head->data.function.args[0]->data.symbol.name == fname) {
        const Expr* d = head->data.function.head;      /* Derivative[a,b] */
        if (d->type == EXPR_FUNCTION && d->data.function.arg_count == 2
            && d->data.function.head->type == EXPR_SYMBOL
            && d->data.function.head->data.symbol.name == SYM_Derivative
            && d->data.function.args[0]->type == EXPR_INTEGER
            && d->data.function.args[1]->type == EXPR_INTEGER) {
            *a = (int)d->data.function.args[0]->data.integer;
            *b = (int)d->data.function.args[1]->data.integer;
            return true;
        }
    }
    return false;
}

/* Build the literal  u[t,x]  (a==b==0) or  Derivative[a,b][u][t,x]. */
static Expr* nd_pde_lit(const char* fname, int a, int b,
                        const char* tvar, const char* xvar) {
    if (a == 0 && b == 0) {
        Expr* ar[2] = { expr_new_symbol(tvar), expr_new_symbol(xvar) };
        return expr_new_function(expr_new_symbol(fname), ar, 2);
    }
    Expr* di[2] = { expr_new_integer(a), expr_new_integer(b) };
    Expr* dhead = expr_new_function(expr_new_symbol(SYM_Derivative), di, 2);
    Expr* uu[1] = { expr_new_symbol(fname) };
    Expr* d2 = expr_new_function(dhead, uu, 1);
    Expr* ar[2] = { expr_new_symbol(tvar), expr_new_symbol(xvar) };
    return expr_new_function(d2, ar, 2);
}

/* NDSolve is HoldAll, so a derivative written with D — D[u[t,x],t],
 * D[u[t,x],{x,2}] — arrives unevaluated instead of as the Derivative[a,b][u]
 * form the matcher expects.  Return a structural copy of `e` in which every
 * D[...] node is evaluated (canonicalizing it to Derivative form) while the
 * surrounding structure (Equal, the boundary/initial values, ...) is left
 * intact — so an equation like Derivative[1,0][u][t,x] == ... is never
 * prematurely reduced. */
static Expr* nd_normalize_derivs(const Expr* e) {
    if (!e) return NULL;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_D)
        return eval_and_free(expr_copy((Expr*)e));
    if (e->type == EXPR_FUNCTION) {
        size_t n = e->data.function.arg_count;
        Expr* head = nd_normalize_derivs(e->data.function.head);
        Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++) args[i] = nd_normalize_derivs(e->data.function.args[i]);
        Expr* r = expr_new_function(head, args, n);
        free(args);
        return r;
    }
    return expr_copy((Expr*)e);
}

/* True if any subexpression of `e` is an application of `fname`. */
static bool nd_contains_func(const Expr* e, const char* fname) {
    if (!e) return false;
    int a, b; const Expr *g0, *g1;
    if (nd_pde_match(e, fname, &a, &b, &g0, &g1)) return true;
    if (e->type == EXPR_FUNCTION) {
        if (nd_contains_func(e->data.function.head, fname)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (nd_contains_func(e->data.function.args[i], fname)) return true;
    }
    return false;
}

/* Recursively track the temporal order, which spatial orders (1,2) appear, and
 * whether an unsupported term (mixed derivative or spatial order > 2) shows up. */
static void nd_scan_orders(const Expr* e, const char* fname, int* torder,
                           bool* has_s1, bool* has_s2, bool* unsupported) {
    if (!e) return;
    int a, b; const Expr *g0, *g1;
    if (nd_pde_match(e, fname, &a, &b, &g0, &g1)) {
        if (a > *torder) *torder = a;
        if (a == 0 && b == 1) *has_s1 = true;
        if (a == 0 && b == 2) *has_s2 = true;
        if (a > 0 && b > 0) *unsupported = true;   /* mixed space-time */
        if (b > 2) *unsupported = true;            /* spatial order > 2 */
    }
    if (e->type == EXPR_FUNCTION) {
        nd_scan_orders(e->data.function.head, fname, torder, has_s1, has_s2, unsupported);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            nd_scan_orders(e->data.function.args[i], fname, torder, has_s1, has_s2, unsupported);
    }
}

/* ------------------------------------------------------------------ *
 *  Arithmetic-tree builders (own their argument Exprs)                *
 * ------------------------------------------------------------------ */
static Expr* nd_scaled(double c, Expr* e) {          /* c * e */
    Expr* ar[2] = { expr_new_real(c), e };
    return expr_new_function(expr_new_symbol(SYM_Times), ar, 2);
}
static Expr* nd_plus2(Expr* a, Expr* b) {
    Expr* ar[2] = { a, b };
    return expr_new_function(expr_new_symbol(SYM_Plus), ar, 2);
}
static Expr* nd_plus3(Expr* a, Expr* b, Expr* c) {
    Expr* ar[3] = { a, b, c };
    return expr_new_function(expr_new_symbol(SYM_Plus), ar, 3);
}

/* Evaluate a one-variable expression (borrowed) at var==value to a double. */
static bool nd_eval_at(const Expr* e, const char* var, double value,
                       NumericSpec spec, double* out) {
    Expr* lit = expr_new_symbol(var);
    Expr* sub = expr_new_real(value);
    Expr* r = nd_replace_all(expr_copy((Expr*)e), &lit, &sub, 1);
    expr_free(lit); expr_free(sub);
    bool ok = r && nd_eval_to_double(r, spec, out);
    expr_free(r);
    return ok;
}

static bool nd_mol_is_option(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    if (h->data.symbol.name != SYM_Rule && h->data.symbol.name != SYM_RuleDelayed) return false;
    return e->data.function.arg_count == 2 && e->data.function.args[0]->type == EXPR_SYMBOL;
}

/* Recursively search an option value for an integer sub-option `name` (handles
 * Method -> {"MethodOfLines", "SpatialDiscretization" -> {..,"MinPoints"->n}}). */
static bool nd_scan_int_subopt(const Expr* e, const char* name, long* out) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if ((head_is((Expr*)e, SYM_Rule) || head_is((Expr*)e, SYM_RuleDelayed))
        && e->data.function.arg_count == 2) {
        const Expr* lhs = e->data.function.args[0];
        const char* key = NULL;
        if (lhs->type == EXPR_STRING) key = lhs->data.string;
        else if (lhs->type == EXPR_SYMBOL) key = lhs->data.symbol.name;
        if (key && strcmp(key, name) == 0) {
            Expr* v = eval_and_free(expr_copy((Expr*)e->data.function.args[1]));
            double d;
            bool ok = v && nd_to_double(v, &d) && d > 0;
            if (ok) *out = (long)(d + 0.5);
            expr_free(v);
            if (ok) return true;
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nd_scan_int_subopt(e->data.function.args[i], name, out)) return true;
    return false;
}

/* ------------------------------------------------------------------ *
 *  Second-order central spatial stencils at interior node `jj`         *
 *  (interior index 0..N-1 ↔ grid index jj+1; boundaries are Dirichlet)*
 * ------------------------------------------------------------------ */
static Expr* nd_neighbor(size_t jj, size_t N, int dir, int torder,
                         Expr** ysym, Expr* bc_left, Expr* bc_right) {
    /* dir = -1 (left) or +1 (right).  Grid boundary → Dirichlet value expr. */
    if (dir < 0 && jj == 0)     return expr_copy(bc_left);
    if (dir > 0 && jj == N - 1) return expr_copy(bc_right);
    size_t nb = (dir < 0) ? jj - 1 : jj + 1;
    return expr_copy(ysym[nb * (size_t)torder + 0]);
}

static Expr* nd_stencil1(size_t jj, size_t N, int torder, double h,
                         Expr** ysym, Expr* bc_left, Expr* bc_right) {
    Expr* left  = nd_neighbor(jj, N, -1, torder, ysym, bc_left, bc_right);
    Expr* right = nd_neighbor(jj, N, +1, torder, ysym, bc_left, bc_right);
    /* (right - left) / (2h) */
    return nd_plus2(nd_scaled(1.0 / (2.0 * h), right),
                    nd_scaled(-1.0 / (2.0 * h), left));
}

static Expr* nd_stencil2(size_t jj, size_t N, int torder, double h,
                         Expr** ysym, Expr* bc_left, Expr* bc_right) {
    Expr* left   = nd_neighbor(jj, N, -1, torder, ysym, bc_left, bc_right);
    Expr* right  = nd_neighbor(jj, N, +1, torder, ysym, bc_left, bc_right);
    Expr* center = expr_copy(ysym[jj * (size_t)torder + 0]);
    double c = 1.0 / (h * h);
    /* (left - 2 center + right) / h^2 */
    return nd_plus3(nd_scaled(c, left), nd_scaled(-2.0 * c, center),
                    nd_scaled(c, right));
}

/* ------------------------------------------------------------------ *
 *  Result assembly: 2-D InterpolatingFunction over (t, x)             *
 * ------------------------------------------------------------------ */
static Expr* nd_mol_build_result(const char* fname, bool applied,
                                 const char* tvar, const char* xvar,
                                 double xmin, double h, size_t nx, int torder,
                                 Expr* bc_left, Expr* bc_right,
                                 const NdSolution* sol, NumericSpec spec) {
    size_t n = sol->n, d = sol->d;
    if (n < 2) return NULL;
    size_t npts = n * nx;
    Expr** entries = malloc(sizeof(Expr*) * npts);
    size_t p = 0;
    for (size_t i = 0; i < n; i++) {
        double ti = sol->ts[i];
        for (size_t g = 0; g < nx; g++) {
            double xg = xmin + (double)g * h;
            double val;
            if (g == 0) {
                if (!nd_eval_at(bc_left, tvar, ti, spec, &val)) val = 0.0;
            } else if (g == nx - 1) {
                if (!nd_eval_at(bc_right, tvar, ti, spec, &val)) val = 0.0;
            } else {
                val = sol->Ys[i * d + (g - 1) * (size_t)torder + 0];
            }
            Expr* coord[2] = { expr_new_real(ti), expr_new_real(xg) };
            Expr* coordL = expr_new_function(expr_new_symbol(SYM_List), coord, 2);
            Expr* pair[2] = { coordL, expr_new_real(val) };
            entries[p++] = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
        }
    }
    Expr* data = expr_new_function(expr_new_symbol(SYM_List), entries, npts);
    free(entries);
    Expr* ifun = eval_and_free(expr_new_function(expr_new_symbol(SYM_Interpolation), &data, 1));
    if (!head_is(ifun, SYM_InterpolatingFunction)) { expr_free(ifun); return NULL; }

    Expr* lhs;
    if (applied) {
        Expr* ar[2] = { expr_new_symbol(tvar), expr_new_symbol(xvar) };
        lhs = expr_new_function(expr_new_symbol(fname), ar, 2);
        Expr* ar2[2] = { expr_new_symbol(tvar), expr_new_symbol(xvar) };
        ifun = expr_new_function(ifun, ar2, 2);
    } else {
        lhs = expr_new_symbol(fname);
    }
    Expr* rule = nd_call2(SYM_Rule, lhs, ifun);
    Expr* inner = expr_new_function(expr_new_symbol(SYM_List), &rule, 1);
    Expr* outer = expr_new_function(expr_new_symbol(SYM_List), &inner, 1);
    return outer;
}

/* ------------------------------------------------------------------ *
 *  Method of Lines solver (Phase 1)                                   *
 * ------------------------------------------------------------------ */
Expr* nd_mol_solve(Expr* res, const NdOpts* o0, const char* forced_method) {
    NdOpts o = *o0;
    if (forced_method && strcmp(forced_method, "MethodOfLines") != 0)
        o.method = intern_symbol(forced_method);
#ifdef USE_MPFR
    if (o0->spec.mode == NUMERIC_MODE_MPFR)
        nd_mol_warn("mppde", "arbitrary-precision PDEs are not yet supported; "
                             "using machine precision");
#endif
    o.spec = numeric_machine_spec();
    o.wp_bits = 53;
    NumericSpec spec = o.spec;

    Expr** A = res->data.function.args;
    size_t argc = res->data.function.arg_count;

    /* ---- positional range count ---- */
    size_t pos_end = 2;
    while (pos_end < argc && !nd_mol_is_option(A[pos_end])) pos_end++;
    size_t nranges = (pos_end >= 2) ? pos_end - 2 : 0;
    if (nranges < 2) return NULL;
    if (nranges > 2) {
        nd_mol_warn("pdedim", "only one spatial dimension is supported in this phase");
        return NULL;
    }

    /* ---- temporal + spatial ranges ---- */
    Expr* rt = A[2]; Expr* rx = A[3];
    if (!head_is(rt, SYM_List) || rt->data.function.arg_count != 3 ||
        rt->data.function.args[0]->type != EXPR_SYMBOL) return NULL;
    if (!head_is(rx, SYM_List) || rx->data.function.arg_count != 3 ||
        rx->data.function.args[0]->type != EXPR_SYMBOL) return NULL;
    const char* tvar = rt->data.function.args[0]->data.symbol.name;
    const char* xvar = rx->data.function.args[0]->data.symbol.name;
    double tmin, tmax, xmin, xmax;
    if (!nd_eval_to_double(rt->data.function.args[1], spec, &tmin)) return NULL;
    if (!nd_eval_to_double(rt->data.function.args[2], spec, &tmax)) return NULL;
    if (!nd_eval_to_double(rx->data.function.args[1], spec, &xmin)) return NULL;
    if (!nd_eval_to_double(rx->data.function.args[2], spec, &xmax)) return NULL;
    if (!(xmax > xmin)) return NULL;

    /* ---- dependent function (single) ---- */
    Expr* funcs = A[1];
    Expr* fitem = funcs;
    if (head_is(funcs, SYM_List)) {
        if (funcs->data.function.arg_count != 1) {
            nd_mol_warn("pdesys", "systems of PDEs are not yet supported");
            return NULL;
        }
        fitem = funcs->data.function.args[0];
    }
    const char* fname; bool applied;
    if (fitem->type == EXPR_SYMBOL) { fname = fitem->data.symbol.name; applied = false; }
    else if (fitem->type == EXPR_FUNCTION && fitem->data.function.head->type == EXPR_SYMBOL
             && fitem->data.function.arg_count == 2) {
        fname = fitem->data.function.head->data.symbol.name; applied = true;
    } else return NULL;

    /* ---- grid resolution ---- */
    long nx_l = 25;
    for (size_t i = pos_end; i < argc; i++) {
        if (!nd_mol_is_option(A[i])) continue;
        if (A[i]->data.function.args[0]->data.symbol.name == SYM_Method)
            nd_scan_int_subopt(A[i], "MinPoints", &nx_l);
    }
    if (nx_l < 5) nx_l = 5;
    if (nx_l > 20000) nx_l = 20000;
    size_t nx = (size_t)nx_l;
    double h = (xmax - xmin) / (double)(nx - 1);
    size_t N = nx - 2;                 /* interior unknowns (Dirichlet) */
    if (N < 1) return NULL;

    /* ---- normalize equations (flatten + canonicalize D -> Derivative) ---- */
    Expr* eqns = A[0];
    Expr** eqsrc; size_t neq;
    if (head_is(eqns, SYM_List)) { neq = eqns->data.function.arg_count; eqsrc = eqns->data.function.args; }
    else { neq = 1; eqsrc = &eqns; }
    Expr** eqitems = malloc(sizeof(Expr*) * neq);
    for (size_t e = 0; e < neq; e++) eqitems[e] = nd_normalize_derivs(eqsrc[e]);

    /* ---- scan orders ---- */
    int torder = 0; bool has_s1 = false, has_s2 = false, unsupported = false;
    for (size_t e = 0; e < neq; e++)
        nd_scan_orders(eqitems[e], fname, &torder, &has_s1, &has_s2, &unsupported);
    if (unsupported || torder < 1 || torder > 2) {
        nd_mol_warn(unsupported ? "pdeord" : "pdetime",
                    unsupported ? "only spatial orders 1 and 2 (unmixed) are supported "
                                  "in this phase"
                                : "temporal order must be 1 or 2 in this phase");
        for (size_t e = 0; e < neq; e++) expr_free(eqitems[e]);
        free(eqitems);
        return NULL;
    }
    (void)has_s1; (void)has_s2;

    /* ---- classify equations: evolution PDE, ICs, Dirichlet BCs ---- */
    Expr* pde_eq = NULL;
    Expr** ic = calloc((size_t)torder, sizeof(Expr*));   /* ic[a]: value in x   */
    Expr* bc_left = NULL; Expr* bc_right = NULL;          /* values in t         */
    double t0 = tmin; bool have_t0 = false;
    for (size_t e = 0; e < neq; e++) {
        Expr* eq = eqitems[e];
        if (!head_is(eq, SYM_Equal) || eq->data.function.arg_count != 2) continue;
        Expr* L = eq->data.function.args[0];
        Expr* R = eq->data.function.args[1];
        int a, b; const Expr *g0, *g1; Expr* val = NULL;
        if (nd_pde_match(L, fname, &a, &b, &g0, &g1) && !nd_contains_func(R, fname)) val = R;
        else if (nd_pde_match(R, fname, &a, &b, &g0, &g1) && !nd_contains_func(L, fname)) val = L;
        if (val) {
            double num;
            bool g0_num = nd_eval_to_double((Expr*)g0, spec, &num);
            bool g1_tvar = (g1->type == EXPR_SYMBOL && g1->data.symbol.name == tvar);
            bool g1_xvar = (g1->type == EXPR_SYMBOL && g1->data.symbol.name == xvar);
            bool g0_tvar = (g0->type == EXPR_SYMBOL && g0->data.symbol.name == tvar);
            (void)g1_tvar;
            double xb;
            bool g1_num = nd_eval_to_double((Expr*)g1, spec, &xb);
            if (b == 0 && g0_num && g1_xvar && a < torder) {
                /* initial condition: order-a time derivative at t0 */
                if (!ic[a]) ic[a] = expr_copy(val);
                t0 = num; have_t0 = true;
            } else if (a == 0 && b == 0 && g0_tvar && g1_num) {
                /* Dirichlet boundary condition */
                if (fabs(xb - xmin) <= 1e-9 * (fabs(xmin) + 1.0)) {
                    if (!bc_left) bc_left = expr_copy(val);
                } else if (fabs(xb - xmax) <= 1e-9 * (fabs(xmax) + 1.0)) {
                    if (!bc_right) bc_right = expr_copy(val);
                }
            }
            continue;
        }
        if (!pde_eq) pde_eq = eq;    /* evolution equation */
    }
    (void)have_t0;

    bool ok = (pde_eq != NULL) && (bc_left != NULL) && (bc_right != NULL);
    for (int a = 0; a < torder && ok; a++) if (!ic[a]) ok = false;
    if (!ok) {
        nd_mol_warn("pdeic", "insufficient initial/boundary conditions (Phase 1 "
                             "requires Dirichlet BCs at both ends and full initial data)");
        for (int a = 0; a < torder; a++) expr_free(ic[a]);
        free(ic); expr_free(bc_left); expr_free(bc_right);
        for (size_t e = 0; e < neq; e++) expr_free(eqitems[e]);
        free(eqitems);
        return NULL;
    }

    /* ---- solve the PDE for its top temporal derivative -> G ---- */
    Expr* R0 = nd_call2(SYM_Subtract, expr_copy(pde_eq->data.function.args[0]),
                                      expr_copy(pde_eq->data.function.args[1]));
    Expr* topLit = nd_pde_lit(fname, torder, 0, tvar, xvar);
    Expr* Psym = expr_new_symbol("NDSolve`Pt");
    Expr* Rp = nd_replace_all(R0, &topLit, &Psym, 1);            /* consumes R0 */
    Expr* dargs[2] = { expr_copy(Rp), expr_copy(Psym) };
    Expr* aE = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), dargs, 2));
    Expr* zero = expr_new_integer(0);
    Expr* bE = nd_replace_all(expr_copy(Rp), &Psym, &zero, 1);
    expr_free(zero);
    Expr* invE = nd_call2(SYM_Power, expr_copy(aE), expr_new_integer(-1));
    Expr* g3[3] = { expr_new_integer(-1), expr_copy(bE), invE };
    Expr* G = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), g3, 3));
    expr_free(topLit); expr_free(Psym); expr_free(Rp); expr_free(aE); expr_free(bE);

    if (!G) {
        nd_mol_warn("pdesolve", "could not solve the PDE for its highest temporal derivative");
        for (int a = 0; a < torder; a++) expr_free(ic[a]);
        free(ic); expr_free(bc_left); expr_free(bc_right);
        for (size_t e = 0; e < neq; e++) expr_free(eqitems[e]);
        free(eqitems);
        return NULL;
    }

    /* ---- reduced-state symbols  w[jj*torder + m] ---- */
    size_t d = N * (size_t)torder;
    Expr** ysym = malloc(sizeof(Expr*) * d);
    for (size_t idx = 0; idx < d; idx++) {
        char buf[48]; snprintf(buf, sizeof(buf), "NDSolve`w%zu", idx);
        ysym[idx] = expr_new_symbol(intern_symbol(buf));
    }

    /* ---- assemble the NdProblem ---- */
    NdProblem P; memset(&P, 0, sizeof(P));
    P.d = d; P.spec = spec; P.tvar = tvar;
    P.eval_monitor = o.eval_monitor;
    P.ysym = ysym;
    P.f = calloc(d, sizeof(Expr*));
    P.Y0 = calloc(d, sizeof(double));
    P.bind_y = malloc(sizeof(NdBind) * d);
    nd_bind_snapshot(&P.bind_t, tvar);
    for (size_t idx = 0; idx < d; idx++)
        nd_bind_snapshot(&P.bind_y[idx], ysym[idx]->data.symbol.name);

    /* Literals of the derivative forms substituted per node. */
    Expr* lit_u   = nd_pde_lit(fname, 0, 0, tvar, xvar);
    Expr* lit_ux  = nd_pde_lit(fname, 0, 1, tvar, xvar);
    Expr* lit_uxx = nd_pde_lit(fname, 0, 2, tvar, xvar);
    Expr* lit_xv  = expr_new_symbol(xvar);
    Expr** lit_ut = malloc(sizeof(Expr*) * (size_t)torder);    /* lit_ut[m] = D[m,0] */
    for (int m = 1; m < torder; m++) lit_ut[m] = nd_pde_lit(fname, m, 0, tvar, xvar);

    for (size_t jj = 0; jj < N && G; jj++) {
        size_t base = jj * (size_t)torder;
        double xj = xmin + (double)(jj + 1) * h;
        /* chain: dU^(m)/dt = U^(m+1) */
        for (int m = 0; m + 1 < torder; m++)
            P.f[base + (size_t)m] = expr_copy(ysym[base + (size_t)m + 1]);
        /* top: substitute stencils + node coordinate into G */
        size_t cnt = 0;
        Expr* lits[8]; Expr* subs[8];
        lits[cnt] = lit_u;   subs[cnt] = expr_copy(ysym[base + 0]);                cnt++;
        lits[cnt] = lit_ux;  subs[cnt] = nd_stencil1(jj, N, torder, h, ysym, bc_left, bc_right); cnt++;
        lits[cnt] = lit_uxx; subs[cnt] = nd_stencil2(jj, N, torder, h, ysym, bc_left, bc_right); cnt++;
        for (int m = 1; m < torder; m++) {
            lits[cnt] = lit_ut[m]; subs[cnt] = expr_copy(ysym[base + (size_t)m]); cnt++;
        }
        lits[cnt] = lit_xv;  subs[cnt] = expr_new_real(xj);                        cnt++;
        P.f[base + (size_t)torder - 1] = nd_replace_all(expr_copy(G), lits, subs, cnt);
        for (size_t s = 0; s < cnt; s++) expr_free(subs[s]);   /* lits are shared */
    }
    expr_free(lit_u); expr_free(lit_ux); expr_free(lit_uxx); expr_free(lit_xv);
    for (int m = 1; m < torder; m++) expr_free(lit_ut[m]);
    free(lit_ut);

    bool build_ok = true;
    for (size_t i = 0; i < d; i++) if (!P.f[i]) build_ok = false;

    /* ---- initial data ---- */
    for (size_t jj = 0; jj < N && build_ok; jj++) {
        double xj = xmin + (double)(jj + 1) * h;
        for (int m = 0; m < torder; m++) {
            double v;
            if (!nd_eval_at(ic[m], xvar, xj, spec, &v)) { build_ok = false; break; }
            P.Y0[jj * (size_t)torder + (size_t)m] = v;
        }
    }
    P.t0 = t0; P.tmin = tmin; P.tmax = tmax;

    /* ---- integrate ---- */
    Expr* result = NULL;
    if (build_ok) {
        const NdStepper* S = nd_lookup_stepper(o.method);
        if (!S) S = nd_default_stepper();
        NdSolution sol; nd_solution_init(&sol, d);
        NdStatus st = nd_integrate(&P, S, &o, &sol);
        nd_bind_restore(&P.bind_t);
        for (size_t i = 0; i < d; i++) nd_bind_restore(&P.bind_y[i]);
        if (st == ND_ERR_MAXSTEPS) nd_mol_warn("mxst", "maximum number of steps reached; returning partial solution");
        else if (st == ND_ERR_STEPSIZE) nd_mol_warn("ndsz", "step size effectively zero; stiffness suspected (try Method->\"BDF\")");
        else if (st == ND_ERR_NONCONV) nd_mol_warn("ndcf", "corrector failed to converge");
        else if (st == ND_ERR_SAMPLE) nd_mol_warn("nrnum", "spatial operator did not evaluate to a number");
        result = nd_mol_build_result(fname, applied, tvar, xvar, xmin, h, nx, torder,
                                     bc_left, bc_right, &sol, spec);
        nd_solution_free(&sol);
    } else {
        nd_bind_restore(&P.bind_t);
        for (size_t i = 0; i < d; i++) nd_bind_restore(&P.bind_y[i]);
    }

    /* ---- cleanup ---- */
    for (size_t i = 0; i < d; i++) { expr_free(P.f[i]); expr_free(ysym[i]); }
    free(P.f); free(ysym); free(P.Y0); free(P.bind_y);
    if (P.jac) {
        for (size_t i = 0; i < d; i++) {
            for (size_t j = 0; j < d; j++) expr_free(P.jac[i][j]);
            free(P.jac[i]);
        }
        free(P.jac);
    }
    expr_free(G);
    for (int a = 0; a < torder; a++) expr_free(ic[a]);
    free(ic); expr_free(bc_left); expr_free(bc_right);
    for (size_t e = 0; e < neq; e++) expr_free(eqitems[e]);
    free(eqitems);
    return result;
}
