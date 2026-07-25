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
#include "ndsolve_stencil.h"
#include "ndsolve_operator.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../eval.h"
#include "../common.h"
#include "../numeric.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define ND_MAX_SORDER 8    /* highest spatial derivative order handled */

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

/* Build  Derivative[0,b][u][tvar, pt]  (b==0 -> u[tvar, pt]) using an explicit
 * boundary-point Expr `pt` (copied), so the literal matches the exact form the
 * boundary condition uses (e.g. an integer 0, not a real 0.0). */
static Expr* nd_bc_lit(const char* fname, int b, const char* tvar, const Expr* pt) {
    Expr* t = expr_new_symbol(tvar);
    Expr* p = expr_copy((Expr*)pt);
    if (b == 0) { Expr* ar[2] = { t, p }; return expr_new_function(expr_new_symbol(fname), ar, 2); }
    Expr* di[2] = { expr_new_integer(0), expr_new_integer(b) };
    Expr* dhead = expr_new_function(expr_new_symbol(SYM_Derivative), di, 2);
    Expr* uu[1] = { expr_new_symbol(fname) };
    Expr* d2 = expr_new_function(dhead, uu, 1);
    Expr* ar[2] = { t, p };
    return expr_new_function(d2, ar, 2);
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

/* Recursively track the temporal order, which spatial orders appear
 * (has_s[1..ND_MAX_SORDER]), and whether an unsupported term (mixed
 * space-time derivative, or spatial order beyond ND_MAX_SORDER) shows up. */
static void nd_scan_orders(const Expr* e, const char* fname, int* torder,
                           bool* has_s, bool* unsupported) {
    if (!e) return;
    int a, b; const Expr *g0, *g1;
    if (nd_pde_match(e, fname, &a, &b, &g0, &g1)) {
        if (a > *torder) *torder = a;
        if (a == 0 && b >= 1) {
            if (b <= ND_MAX_SORDER) has_s[b] = true;
            else *unsupported = true;
        }
        if (a > 0 && b > 0) *unsupported = true;   /* mixed space-time */
    }
    if (e->type == EXPR_FUNCTION) {
        nd_scan_orders(e->data.function.head, fname, torder, has_s, unsupported);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            nd_scan_orders(e->data.function.args[i], fname, torder, has_s, unsupported);
    }
}

/* One occurrence of the dependent function in an equation. */
typedef struct { int a, b; const Expr* arg0; const Expr* arg1; } NdFApp;

static void nd_collect_fapps(const Expr* e, const char* fname,
                             NdFApp* out, int* n, int max) {
    if (!e) return;
    int a, b; const Expr *g0, *g1;
    if (nd_pde_match(e, fname, &a, &b, &g0, &g1) && *n < max) {
        out[*n].a = a; out[*n].b = b; out[*n].arg0 = g0; out[*n].arg1 = g1; (*n)++;
    }
    if (e->type == EXPR_FUNCTION) {
        nd_collect_fapps(e->data.function.head, fname, out, n, max);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            nd_collect_fapps(e->data.function.args[i], fname, out, n, max);
    }
}

/* ------------------------------------------------------------------ *
 *  Arithmetic-tree builders (own their argument Exprs)                *
 * ------------------------------------------------------------------ */
static Expr* nd_scaled(double c, Expr* e) {          /* c * e */
    Expr* ar[2] = { expr_new_real(c), e };
    return expr_new_function(expr_new_symbol(SYM_Times), ar, 2);
}
static Expr* nd_plus2any(Expr* a, Expr* b) {         /* a + b (consumes both) */
    Expr* ar[2] = { a, b };
    return expr_new_function(expr_new_symbol(SYM_Plus), ar, 2);
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

/* Value of grid node `g` as a symbolic Expr (owned): an interior reduced-state
 * symbol, a boundary value expression (Dirichlet/Neumann/Robin, at grid 0 or
 * nx-1), or — in the periodic case — the wrapped interior state symbol. */
static Expr* nd_node_value(int g, int nx, bool periodic, int torder,
                           Expr** ysym, Expr* bc_left, Expr* bc_right) {
    if (periodic) {
        int M = nx - 1;
        int gw = ((g % M) + M) % M;
        return expr_copy(ysym[(size_t)gw * (size_t)torder + 0]);
    }
    if (g <= 0)        return expr_copy(bc_left);
    if (g >= nx - 1)   return expr_copy(bc_right);
    return expr_copy(ysym[(size_t)(g - 1) * (size_t)torder + 0]);
}

/* ------------------------------------------------------------------ *
 *  Spatial-derivative stencil at grid node `jgrid` (Fornberg weights)  *
 *  as a symbolic combination of neighbouring node values.  Central and *
 *  boundary-shifted for Dirichlet/Neumann/Robin; symmetric and wrapped *
 *  cyclically when `periodic`.                                          *
 * ------------------------------------------------------------------ */
static Expr* nd_mol_stencil_expr(int jgrid, int nx, bool periodic, int deriv, int order,
                                 double h, int torder, Expr** ysym,
                                 Expr* bc_left, Expr* bc_right) {
    int idx[64]; double w[64]; int npts;
    if (periodic) {
        int M = nx - 1;
        int p = (deriv <= 2) ? (order + 1) : (deriv + order);
        if (p < deriv + 1) p = deriv + 1;
        if (p > M) p = M;
        if (p > 64) p = 64;
        int half = p / 2;
        double xs[64] = { 0 };
        for (int k = 0; k < p; k++) { idx[k] = jgrid - half + k; xs[k] = (double)(idx[k] - jgrid); }
        nd_fornberg_weights(deriv, xs, p, 0.0, w);
        double hd = pow(h, (double)deriv);
        for (int k = 0; k < p; k++) w[k] /= hd;
        npts = p;
    } else {
        nd_stencil_build(jgrid, nx, deriv, order, h, idx, w, &npts);
    }
    Expr** terms = malloc(sizeof(Expr*) * (size_t)npts);
    int nt = 0;
    for (int k = 0; k < npts; k++) {
        if (w[k] == 0.0) continue;
        Expr* val = nd_node_value(idx[k], nx, periodic, torder, ysym, bc_left, bc_right);
        terms[nt++] = nd_scaled(w[k], val);
    }
    Expr* r;
    if (nt == 0) { free(terms); r = expr_new_integer(0); }
    else if (nt == 1) { r = terms[0]; free(terms); }
    else { r = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)nt); free(terms); }
    return r;
}

/* Build the boundary-node value expression for a general linear boundary
 * condition  a*u + b*u_x + r = 0  at grid node `jbound` (0 or nx-1), by
 * eliminating u_x with a one-sided finite-difference stencil:
 *   u_x(jbound) = w0*u_bound + Ux_rest  (Ux_rest over interior nodes)
 *   =>  u_bound = -(r + b*Ux_rest) / (a + b*w0).
 * a, b, r are borrowed (copied in); the returned Expr is evaluated to fold
 * constants (Dirichlet -> the value r-term). */
static Expr* nd_bc_make_value(Expr* a, Expr* b, Expr* r, int jbound, int nx,
                              int dord, double h, int torder, Expr** ysym) {
    int idx[64]; double w[64]; int npts;
    nd_stencil_build(jbound, nx, 1, dord, h, idx, w, &npts);
    double w0 = 0.0;
    Expr** terms = malloc(sizeof(Expr*) * (size_t)npts);
    int nt = 0;
    for (int k = 0; k < npts; k++) {
        if (idx[k] == jbound) { w0 = w[k]; continue; }
        if (w[k] == 0.0) continue;
        int gi = idx[k];
        if (gi < 1 || gi > nx - 2) continue;      /* only interior states  */
        terms[nt++] = nd_scaled(w[k], expr_copy(ysym[(size_t)(gi - 1) * (size_t)torder + 0]));
    }
    Expr* ux_rest;
    if (nt == 0) { free(terms); ux_rest = expr_new_integer(0); }
    else if (nt == 1) { ux_rest = terms[0]; free(terms); }
    else { ux_rest = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)nt); free(terms); }

    Expr* denom = nd_plus2any(expr_copy(a), nd_scaled(w0, expr_copy(b)));
    Expr* numer = nd_scaled(-1.0, nd_plus2any(expr_copy(r),
                            expr_new_function(expr_new_symbol(SYM_Times),
                                (Expr*[]){ expr_copy(b), ux_rest }, 2)));
    Expr* inv = nd_call2(SYM_Power, denom, expr_new_integer(-1));
    Expr* val = expr_new_function(expr_new_symbol(SYM_Times),
                                  (Expr*[]){ numer, inv }, 2);
    return eval_and_free(val);
}

/* ------------------------------------------------------------------ *
 *  Result assembly: 2-D InterpolatingFunction over (t, x)             *
 * ------------------------------------------------------------------ */
/* Evaluate a (possibly state-dependent) boundary expression at the currently
 * bound point (t and reduced-state symbols set by the caller). */
static bool nd_bound_eval(Expr* bc, NumericSpec spec, double* out) {
    Expr* v = eval_and_free(expr_copy(bc));
    if (!v) return false;
    Expr* nv = numericalize(v, spec);
    expr_free(v);
    bool ok = nv && nd_to_double(nv, out) && isfinite(*out);
    expr_free(nv);
    return ok;
}

static Expr* nd_mol_build_result(NdProblem* P, const char* fname, bool applied,
                                 const char* tvar, const char* xvar,
                                 double xmin, double h, size_t nx, int torder,
                                 bool periodic, Expr* bc_left, Expr* bc_right,
                                 const NdSolution* sol, NumericSpec spec) {
    size_t n = sol->n, d = sol->d;
    if (n < 2) return NULL;
    size_t npts = n * nx;
    Expr** entries = malloc(sizeof(Expr*) * npts);
    size_t p = 0;
    for (size_t i = 0; i < n; i++) {
        double ti = sol->ts[i];
        /* Bind t and all reduced-state symbols so that boundary-value
         * expressions (Neumann/Robin depend on interior states) evaluate. */
        if (!periodic) {
            nd_bind_set(&P->bind_t, expr_new_real(ti));
            for (size_t k = 0; k < d; k++)
                nd_bind_set(&P->bind_y[k], expr_new_real(sol->Ys[i * d + k]));
        }
        for (size_t g = 0; g < nx; g++) {
            double xg = xmin + (double)g * h;
            double val;
            if (periodic) {
                size_t gg = (g == nx - 1) ? 0 : g;   /* last node = image of node 0 */
                val = sol->Ys[i * d + gg * (size_t)torder + 0];
            } else if (g == 0) {
                if (!nd_bound_eval(bc_left, spec, &val)) val = 0.0;
            } else if (g == nx - 1) {
                if (!nd_bound_eval(bc_right, spec, &val)) val = 0.0;
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

    /* ---- grid resolution + finite-difference order + Compiled ---- */
    long nx_l = 25;
    long dord_l = 4;                   /* DifferenceOrder default (WL-faithful) */
    bool compiled = true;              /* Compiled -> False forces symbolic RHS */
    const char* SYM_Compiled = intern_symbol("Compiled");
    for (size_t i = pos_end; i < argc; i++) {
        if (!nd_mol_is_option(A[i])) continue;
        const char* opt = A[i]->data.function.args[0]->data.symbol.name;
        if (opt == SYM_Method) {
            nd_scan_int_subopt(A[i], "MinPoints", &nx_l);
            nd_scan_int_subopt(A[i], "DifferenceOrder", &dord_l);
        } else if (opt == SYM_Compiled) {
            Expr* v = eval_and_free(expr_copy(A[i]->data.function.args[1]));
            if (v && v->type == EXPR_SYMBOL && v->data.symbol.name == intern_symbol("False"))
                compiled = false;
            expr_free(v);
        }
    }
    if (nx_l < 5) nx_l = 5;
    if (nx_l > 20000) nx_l = 20000;
    size_t nx = (size_t)nx_l;
    double h = (xmax - xmin) / (double)(nx - 1);
    if (dord_l < 1) dord_l = 1;
    if (dord_l > (long)nx - 1) dord_l = (long)nx - 1;   /* need <= nx nodes */
    int dord = (int)dord_l;

    /* ---- normalize equations (flatten + canonicalize D -> Derivative) ---- */
    Expr* eqns = A[0];
    Expr** eqsrc; size_t neq;
    if (head_is(eqns, SYM_List)) { neq = eqns->data.function.arg_count; eqsrc = eqns->data.function.args; }
    else { neq = 1; eqsrc = &eqns; }
    Expr** eqitems = malloc(sizeof(Expr*) * neq);
    for (size_t e = 0; e < neq; e++) eqitems[e] = nd_normalize_derivs(eqsrc[e]);

    /* ---- scan orders ---- */
    int torder = 0; bool unsupported = false;
    bool has_s[ND_MAX_SORDER + 1] = { 0 };
    for (size_t e = 0; e < neq; e++)
        nd_scan_orders(eqitems[e], fname, &torder, has_s, &unsupported);
    /* list of distinct spatial derivative orders present */
    int sord[ND_MAX_SORDER]; int nsord = 0;
    for (int b = 1; b <= ND_MAX_SORDER; b++) if (has_s[b]) sord[nsord++] = b;
    if (unsupported || torder < 1 || torder > 2) {
        nd_mol_warn(unsupported ? "pdeord" : "pdetime",
                    unsupported ? "mixed space-time derivatives (and spatial order "
                                  "above 8) are not supported"
                                : "temporal order must be 1 or 2 in this phase");
        for (size_t e = 0; e < neq; e++) expr_free(eqitems[e]);
        free(eqitems);
        return NULL;
    }

    /* ---- classify equations: evolution PDE, ICs, boundary conditions ----
     * A boundary condition is stored as the coefficients of the general linear
     * form  a*u + b*u_x + r == 0  at a boundary point (Dirichlet: b=0; Neumann:
     * a=0; Robin: both), which the value elimination below turns into a
     * boundary-node value expression. `u[t,xmin]==u[t,xmax]` marks periodicity. */
    Expr* pde_eq = NULL;
    Expr** ic = calloc((size_t)torder, sizeof(Expr*));   /* ic[a]: value in x   */
    Expr *aL = NULL, *bL = NULL, *rL = NULL;             /* left boundary spec  */
    Expr *aR = NULL, *bR = NULL, *rR = NULL;             /* right boundary spec */
    bool periodic = false;
    double t0 = tmin;
    for (size_t e = 0; e < neq; e++) {
        Expr* eq = eqitems[e];
        if (!head_is(eq, SYM_Equal) || eq->data.function.arg_count != 2) continue;
        Expr* L = eq->data.function.args[0];
        Expr* R = eq->data.function.args[1];

        NdFApp fa[16]; int nfa = 0;
        nd_collect_fapps(eq, fname, fa, &nfa, 16);
        bool interior = false, is_ic = false; int ic_ord = 0;
        int nbnd = 0; const Expr* bnd_pt = NULL; double bnd_x = 0.0;
        for (int i = 0; i < nfa; i++) {
            bool a0_t = fa[i].arg0->type == EXPR_SYMBOL && fa[i].arg0->data.symbol.name == tvar;
            bool a1_x = fa[i].arg1->type == EXPR_SYMBOL && fa[i].arg1->data.symbol.name == xvar;
            double n0, n1;
            bool a0_num = nd_eval_to_double((Expr*)fa[i].arg0, spec, &n0);
            bool a1_num = nd_eval_to_double((Expr*)fa[i].arg1, spec, &n1);
            if (a0_t && a1_x) interior = true;
            else if (a0_num && a1_x) { is_ic = true; ic_ord = fa[i].a; t0 = n0; }
            else if (a0_t && a1_num) { nbnd++; bnd_pt = fa[i].arg1; bnd_x = n1; }
        }

        if (interior) { if (!pde_eq) pde_eq = eq; continue; }

        if (is_ic) {
            int a, b; const Expr *g0, *g1; Expr* val = NULL;
            if (nd_pde_match(L, fname, &a, &b, &g0, &g1) && !nd_contains_func(R, fname)) val = R;
            else if (nd_pde_match(R, fname, &a, &b, &g0, &g1) && !nd_contains_func(L, fname)) val = L;
            if (val && ic_ord < torder && !ic[ic_ord]) ic[ic_ord] = expr_copy(val);
            continue;
        }

        if (nbnd > 0) {
            /* periodic: u[t,c1] == u[t,c2] (order 0, distinct boundaries) */
            int la, lb, ra, rb; const Expr *l0, *l1, *r0, *r1;
            bool Lm = nd_pde_match(L, fname, &la, &lb, &l0, &l1);
            bool Rm = nd_pde_match(R, fname, &ra, &rb, &r0, &r1);
            if (Lm && Rm && la == 0 && lb == 0 && ra == 0 && rb == 0) {
                double x1, x2;
                if (nd_eval_to_double((Expr*)l1, spec, &x1) &&
                    nd_eval_to_double((Expr*)r1, spec, &x2) && fabs(x1 - x2) > 1e-9) {
                    periodic = true; continue;
                }
            }
            /* general linear BC: extract a,b,r by substituting u,u_x -> symbols */
            Expr* res = nd_call2(SYM_Subtract, expr_copy(L), expr_copy(R));
            Expr* Ulit = nd_bc_lit(fname, 0, tvar, bnd_pt);
            Expr* Uxlit = nd_bc_lit(fname, 1, tvar, bnd_pt);
            Expr* Us = expr_new_symbol("NDSolve`Ub");
            Expr* Uxs = expr_new_symbol("NDSolve`Uxb");
            Expr* subL[2] = { Ulit, Uxlit }; Expr* subR[2] = { Us, Uxs };
            Expr* res2 = nd_replace_all(res, subL, subR, 2);   /* consumes res */
            Expr* da[2] = { expr_copy(res2), expr_copy(Us) };
            Expr* aco = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), da, 2));
            Expr* db[2] = { expr_copy(res2), expr_copy(Uxs) };
            Expr* bco = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), db, 2));
            Expr* z0 = expr_new_integer(0), *z1 = expr_new_integer(0);
            Expr* zl[2] = { Us, Uxs }; Expr* zs[2] = { z0, z1 };
            Expr* rco = nd_replace_all(expr_copy(res2), zl, zs, 2);
            expr_free(z0); expr_free(z1);
            expr_free(Ulit); expr_free(Uxlit); expr_free(Us); expr_free(Uxs); expr_free(res2);
            bool left = fabs(bnd_x - xmin) <= 1e-9 * (fabs(xmin) + 1.0);
            if (left && !aL) { aL = aco; bL = bco; rL = rco; }
            else if (!left && !aR) { aR = aco; bR = bco; rR = rco; }
            else { expr_free(aco); expr_free(bco); expr_free(rco); }
            continue;
        }
        if (!pde_eq) pde_eq = eq;    /* evolution equation */
    }

    bool have_bc = periodic || (aL && aR);
    bool ok = (pde_eq != NULL) && have_bc;
    for (int a = 0; a < torder && ok; a++) if (!ic[a]) ok = false;
    if (!ok) {
        nd_mol_warn("pdeic", "insufficient initial/boundary conditions (need a "
                             "condition at each end — Dirichlet/Neumann/Robin or "
                             "periodic — and full initial data)");
        for (int a = 0; a < torder; a++) expr_free(ic[a]);
        free(ic);
        expr_free(aL); expr_free(bL); expr_free(rL);
        expr_free(aR); expr_free(bR); expr_free(rR);
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
        free(ic);
        expr_free(aL); expr_free(bL); expr_free(rL);
        expr_free(aR); expr_free(bR); expr_free(rR);
        for (size_t e = 0; e < neq; e++) expr_free(eqitems[e]);
        free(eqitems);
        return NULL;
    }

    /* ---- unknown grid nodes ----
     * Dirichlet/Neumann/Robin: interior nodes 1..nx-2 (grid base 1).
     * Periodic: nodes 0..nx-2 (node nx-1 is the image of node 0). */
    size_t nunk = periodic ? (nx - 1) : (nx - 2);
    int gbase = periodic ? 0 : 1;
    size_t d = nunk * (size_t)torder;

    /* ---- reduced-state symbols ---- */
    Expr** ysym = malloc(sizeof(Expr*) * d);
    for (size_t idx = 0; idx < d; idx++) {
        char buf[48]; snprintf(buf, sizeof(buf), "NDSolve`w%zu", idx);
        ysym[idx] = expr_new_symbol(intern_symbol(buf));
    }

    /* ---- eliminate boundary nodes into value expressions (non-periodic) ---- */
    Expr* bc_left = NULL; Expr* bc_right = NULL;
    if (!periodic) {
        bc_left  = nd_bc_make_value(aL, bL, rL, 0, (int)nx, dord, h, torder, ysym);
        bc_right = nd_bc_make_value(aR, bR, rR, (int)nx - 1, (int)nx, dord, h, torder, ysym);
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

    /* Literals of the derivative forms substituted per node: u[t,x], each
     * spatial order present, each lower temporal order, and the coordinate x. */
    Expr* lit_u   = nd_pde_lit(fname, 0, 0, tvar, xvar);
    Expr* lit_xv  = expr_new_symbol(xvar);
    Expr** lit_sx = malloc(sizeof(Expr*) * (size_t)(nsord ? nsord : 1)); /* D[0,b] */
    for (int i = 0; i < nsord; i++) lit_sx[i] = nd_pde_lit(fname, 0, sord[i], tvar, xvar);
    Expr** lit_ut = malloc(sizeof(Expr*) * (size_t)torder);    /* lit_ut[m] = D[m,0] */
    for (int m = 1; m < torder; m++) lit_ut[m] = nd_pde_lit(fname, m, 0, tvar, xvar);

    for (size_t u = 0; u < nunk && G; u++) {
        size_t base = u * (size_t)torder;
        int jgrid = (int)u + gbase;
        double xj = xmin + (double)jgrid * h;
        /* chain: dU^(m)/dt = U^(m+1) */
        for (int m = 0; m + 1 < torder; m++)
            P.f[base + (size_t)m] = expr_copy(ysym[base + (size_t)m + 1]);
        /* top: substitute stencils + lower time derivatives + coordinate into G */
        size_t cnt = 0;
        Expr* lits[ND_MAX_SORDER + 4]; Expr* subs[ND_MAX_SORDER + 4];
        lits[cnt] = lit_u; subs[cnt] = expr_copy(ysym[base + 0]); cnt++;
        for (int i = 0; i < nsord; i++) {
            lits[cnt] = lit_sx[i];
            subs[cnt] = nd_mol_stencil_expr(jgrid, (int)nx, periodic, sord[i], dord, h,
                                            torder, ysym, bc_left, bc_right);
            cnt++;
        }
        for (int m = 1; m < torder; m++) {
            lits[cnt] = lit_ut[m]; subs[cnt] = expr_copy(ysym[base + (size_t)m]); cnt++;
        }
        lits[cnt] = lit_xv; subs[cnt] = expr_new_real(xj); cnt++;
        P.f[base + (size_t)torder - 1] = nd_replace_all(expr_copy(G), lits, subs, cnt);
        for (size_t s = 0; s < cnt; s++) expr_free(subs[s]);   /* lits are shared */
    }
    expr_free(lit_u); expr_free(lit_xv);
    for (int i = 0; i < nsord; i++) expr_free(lit_sx[i]);
    free(lit_sx);
    for (int m = 1; m < torder; m++) expr_free(lit_ut[m]);
    free(lit_ut);

    bool build_ok = true;
    for (size_t i = 0; i < d; i++) if (!P.f[i]) build_ok = false;

    /* ---- initial data ---- */
    for (size_t u = 0; u < nunk && build_ok; u++) {
        double xj = xmin + (double)((int)u + gbase) * h;
        for (int m = 0; m < torder; m++) {
            double v;
            if (!nd_eval_at(ic[m], xvar, xj, spec, &v)) { build_ok = false; break; }
            P.Y0[u * (size_t)torder + (size_t)m] = v;
        }
    }
    P.t0 = t0; P.tmin = tmin; P.tmax = tmax;

    /* ---- compile a linear operator (fast RHS / exact Jacobian) if possible ---- */
    if (build_ok && compiled) P.op = nd_operator_try_build(&P);

    /* Parabolic problems (a diffusion term with first-order time evolution) are
     * stiff; default them to BDF when the user didn't force a method. */
    bool parabolic = (torder == 1) &&
                     (has_s[2] || has_s[4] || has_s[6] || has_s[8]);

    /* ---- integrate ---- */
    Expr* result = NULL;
    if (build_ok) {
        /* "MethodOfLines" is the spatial controller, not a time stepper — treat
         * it (and Automatic) as unset so stiffness auto-selection can apply. */
        const char* tim = o.method;
        if (tim && strcmp(tim, "MethodOfLines") == 0) tim = NULL;
        const NdStepper* S = tim ? nd_lookup_stepper(tim)
                           : (parabolic ? nd_lookup_stepper("BDF") : nd_default_stepper());
        if (!S) S = nd_default_stepper();
        NdSolution sol; nd_solution_init(&sol, d);
        NdStatus st = nd_integrate(&P, S, &o, &sol);
        if (st == ND_ERR_MAXSTEPS) nd_mol_warn("mxst", "maximum number of steps reached; returning partial solution");
        else if (st == ND_ERR_STEPSIZE) nd_mol_warn("ndsz", "step size effectively zero; stiffness suspected (try Method->\"BDF\")");
        else if (st == ND_ERR_NONCONV) nd_mol_warn("ndcf", "corrector failed to converge");
        else if (st == ND_ERR_SAMPLE) nd_mol_warn("nrnum", "spatial operator did not evaluate to a number");
        /* build the result while the bindings are still live (boundary/periodic
         * node values may depend on the reduced-state symbols). */
        result = nd_mol_build_result(&P, fname, applied, tvar, xvar, xmin, h, nx, torder,
                                     periodic, bc_left, bc_right, &sol, spec);
        nd_solution_free(&sol);
    }
    nd_bind_restore(&P.bind_t);
    for (size_t i = 0; i < d; i++) nd_bind_restore(&P.bind_y[i]);

    /* ---- cleanup ---- */
    nd_operator_free(P.op);
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
    free(ic);
    expr_free(aL); expr_free(bL); expr_free(rL);
    expr_free(aR); expr_free(bR); expr_free(rR);
    expr_free(bc_left); expr_free(bc_right);
    for (size_t e = 0; e < neq; e++) expr_free(eqitems[e]);
    free(eqitems);
    return result;
}
