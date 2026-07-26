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
 * Scope: one spatial dimension, temporal order 1 or 2, spatial orders 1-8,
 * Dirichlet/Neumann/Robin/periodic boundary conditions, Fornberg stencils, the
 * compiled linear operator, complex (Schrödinger) PDEs, MPFR precision, and
 * nonlinear/variable-coefficient PDEs through the symbolic sampler; two spatial
 * dimensions in nd_mol_solve_2d.  Coupled *systems* of 1-D PDEs are handled by
 * nd_mol_solve_system (block-per-function reduced state), which also offers
 * donor-cell / Lax-Friedrichs upwind schemes for advective terms; the scalar
 * path here supports the same upwind options.
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

static Expr* nd_mol_solve_2d(Expr* res, const NdOpts* o0, const char* forced_method);
static Expr* nd_mol_solve_system(Expr* res, const NdOpts* o0, const char* forced_method);

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

/* ================================================================== *
 *  Complex-valued PDEs (Schrödinger): Re/Im realification            *
 * ================================================================== */

/* True if `e` involves the imaginary unit (a Complex[...] node or symbol I). */
bool nd_has_imaginary(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol.name, "I") == 0;
    if (e->type == EXPR_FUNCTION) {
        const Expr* h = e->data.function.head;
        if (h->type == EXPR_SYMBOL && strcmp(h->data.symbol.name, "Complex") == 0) return true;
        if (nd_has_imaginary(h)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (nd_has_imaginary(e->data.function.args[i])) return true;
    }
    return false;
}

/* Evaluate `e` (numericalized) and extract its real and imaginary parts. */
bool nd_split_reim(Expr* val, NumericSpec spec, double* re, double* im) {
    Expr* nv = numericalize(val, spec);
    if (!nv) return false;
    Expr* reE = eval_and_free(nd_call1(intern_symbol("Re"), expr_copy(nv)));
    Expr* imE = eval_and_free(nd_call1(intern_symbol("Im"), expr_copy(nv)));
    expr_free(nv);
    bool ok = reE && imE && nd_to_double(reE, re) && nd_to_double(imE, im)
              && isfinite(*re) && isfinite(*im);
    expr_free(reE); expr_free(imE);
    return ok;
}

/* Evaluate a one-variable expression at var==value to a complex (re, im). */
static bool nd_eval_complex_at(const Expr* e, const char* var, double value,
                               NumericSpec spec, double* re, double* im) {
    Expr* lit = expr_new_symbol(var);
    Expr* sub = expr_new_real(value);
    Expr* r = nd_replace_all(expr_copy((Expr*)e), &lit, &sub, 1);
    expr_free(lit); expr_free(sub);
    if (!r) return false;
    bool ok = nd_split_reim(r, spec, re, im);
    expr_free(r);
    return ok;
}

/* Evaluate a (bound) boundary expression to a complex (re, im). */
static bool nd_bound_eval_complex(Expr* bc, NumericSpec spec, double* re, double* im) {
    Expr* raw = eval_and_free(expr_copy(bc));
    if (!raw) return false;
    bool ok = nd_split_reim(raw, spec, re, im);
    expr_free(raw);
    return ok;
}

/* Transform the complex first-order system dZ/dt = f(Z) (Z = P->ysym, length
 * nunk) into the real system of the interleaved parts (2k = Re Z_k, 2k+1 = Im).
 * Each f_k is expanded with Z_j -> wR_j + I wI_j and split by ComplexExpand into
 * its Re/Im parts.  Replaces P->f/ysym/Y0/bind_y/d in place. */
bool nd_realify(NdProblem* P, size_t nunk, const double* Y0re, const double* Y0im) {
    size_t d = P->d;                 /* == nunk (temporal order 1) */
    size_t dR = 2 * nunk;
    Expr** ysymR = malloc(sizeof(Expr*) * dR);
    for (size_t i = 0; i < dR; i++) {
        /* keep the "NDSolve`w<index>" convention the operator builder parses */
        char b[48]; snprintf(b, sizeof b, "NDSolve`w%zu", i);
        ysymR[i] = expr_new_symbol(intern_symbol(b));
    }
    /* substitution Z_k -> wR_k + I wI_k */
    Expr** subs = malloc(sizeof(Expr*) * nunk);
    for (size_t k = 0; k < nunk; k++) {
        Expr* iw = nd_call2("Times", expr_new_symbol("I"), expr_copy(ysymR[2*k+1]));
        subs[k] = nd_call2("Plus", expr_copy(ysymR[2*k]), iw);
    }
    Expr** fR = calloc(dR, sizeof(Expr*));
    bool ok = true;
    for (size_t k = 0; k < nunk && ok; k++) {
        Expr* fsub = nd_replace_all(expr_copy(P->f[k]), P->ysym, subs, nunk);
        Expr* reP = nd_call1("Re", expr_copy(fsub));
        Expr* imP = nd_call1("Im", expr_copy(fsub));
        expr_free(fsub);
        fR[2*k]   = eval_and_free(nd_call1("ComplexExpand", reP));
        fR[2*k+1] = eval_and_free(nd_call1("ComplexExpand", imP));
        if (!fR[2*k] || !fR[2*k+1]) ok = false;
    }
    for (size_t k = 0; k < nunk; k++) expr_free(subs[k]);
    free(subs);
    if (!ok) {
        for (size_t i = 0; i < dR; i++) expr_free(fR[i]);
        free(fR);
        for (size_t k = 0; k < nunk; k++) expr_free(ysymR[k*2]), expr_free(ysymR[k*2+1]);
        free(ysymR);
        return false;
    }
    /* tear down the old complex state, install the real one */
    for (size_t i = 0; i < d; i++) nd_bind_restore(&P->bind_y[i]);
    for (size_t i = 0; i < d; i++) { expr_free(P->f[i]); expr_free(P->ysym[i]); }
    free(P->f); free(P->ysym); free(P->Y0); free(P->bind_y);
    P->d = dR; P->ysym = ysymR; P->f = fR;
    P->Y0 = malloc(sizeof(double) * dR);
    for (size_t k = 0; k < nunk; k++) { P->Y0[2*k] = Y0re[k]; P->Y0[2*k+1] = Y0im[k]; }
    P->bind_y = malloc(sizeof(NdBind) * dR);
    for (size_t i = 0; i < dR; i++) nd_bind_snapshot(&P->bind_y[i], ysymR[i]->data.symbol.name);
    return true;
}

/* Build one real InterpolatingFunction over (t,x) selecting either the Re
 * (part=0) or Im (part=1) component of a realified complex solution. */
static Expr* nd_build_cplx_if(NdProblem* P, const char* tvar, double xmin, double h,
                              size_t nx, size_t nunk, int gbase, bool periodic,
                              Expr* bc_left, Expr* bc_right, const NdSolution* sol,
                              NumericSpec spec, int part) {
    size_t n = sol->n, d = sol->d;
    size_t npts = n * nx;
    Expr** entries = malloc(sizeof(Expr*) * npts);
    size_t p = 0;
    for (size_t i = 0; i < n; i++) {
        double ti = sol->ts[i];
        if (!periodic) {
            nd_bind_set(&P->bind_t, expr_new_real(ti));
            for (size_t k = 0; k < d; k++) nd_bind_set(&P->bind_y[k], expr_new_real(sol->Ys[i*d + k]));
        }
        for (size_t g = 0; g < nx; g++) {
            double xg = xmin + (double)g * h;
            double val = 0.0;
            bool interior = periodic ? (g < nx - 1)
                                     : (g >= 1 && g + 1 < nx);
            if (interior) {
                size_t node = periodic ? g : (g - 1);   /* interior node index */
                val = sol->Ys[i*d + 2*node + (size_t)part];
            } else if (periodic) {
                val = sol->Ys[i*d + 0 + (size_t)part];   /* wrap to node 0 */
            } else {
                double re, im;
                Expr* bc = (g == 0) ? bc_left : bc_right;
                if (nd_bound_eval_complex(bc, spec, &re, &im)) val = part ? im : re;
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
    (void)tvar; (void)nunk; (void)gbase;
    return ifun;
}

/* Assemble the complex result: u -> Function[{t,x}, ifRe[t,x] + I ifIm[t,x]]. */
static Expr* nd_mol_build_result_complex(NdProblem* P, const char* fname, bool applied,
                                         const char* tvar, const char* xvar, double xmin,
                                         double h, size_t nx, size_t nunk, int gbase,
                                         bool periodic, Expr* bc_left, Expr* bc_right,
                                         const NdSolution* sol, NumericSpec spec) {
    if (sol->n < 2) return NULL;
    Expr* ifRe = nd_build_cplx_if(P, tvar, xmin, h, nx, nunk, gbase, periodic,
                                  bc_left, bc_right, sol, spec, 0);
    Expr* ifIm = nd_build_cplx_if(P, tvar, xmin, h, nx, nunk, gbase, periodic,
                                  bc_left, bc_right, sol, spec, 1);
    if (!ifRe || !ifIm) { expr_free(ifRe); expr_free(ifIm); return NULL; }
    /* body = ifRe[t,x] + I ifIm[t,x] */
    Expr* ra[2] = { expr_new_symbol(tvar), expr_new_symbol(xvar) };
    Expr* reApp = expr_new_function(ifRe, ra, 2);
    Expr* ia[2] = { expr_new_symbol(tvar), expr_new_symbol(xvar) };
    Expr* imApp = expr_new_function(ifIm, ia, 2);
    Expr* iTimes = nd_call2("Times", expr_new_symbol("I"), imApp);
    Expr* body = nd_call2("Plus", reApp, iTimes);
    Expr* lhs, *rhs;
    if (applied) {
        Expr* la[2] = { expr_new_symbol(tvar), expr_new_symbol(xvar) };
        lhs = expr_new_function(expr_new_symbol(fname), la, 2);
        rhs = body;
    } else {
        lhs = expr_new_symbol(fname);
        Expr* params[2] = { expr_new_symbol(tvar), expr_new_symbol(xvar) };
        Expr* plist = expr_new_function(expr_new_symbol(SYM_List), params, 2);
        Expr* fargs[2] = { plist, body };
        rhs = expr_new_function(expr_new_symbol("Function"), fargs, 2);
    }
    Expr* rule = nd_call2(SYM_Rule, lhs, rhs);
    Expr* inner = expr_new_function(expr_new_symbol(SYM_List), &rule, 1);
    return expr_new_function(expr_new_symbol(SYM_List), &inner, 1);
}

/* ================================================================== *
 *  Coupled systems of 1-D PDEs                                        *
 *                                                                     *
 *  {u1, u2, ...}(t, x) evolving under a coupled system of evolution   *
 *  equations.  Each dependent function gets its own block of interior *
 *  reduced-state unknowns in one global first-order ODE vector; the   *
 *  per-node right-hand side of function f substitutes the finite-     *
 *  difference stencils / node values of *every* function (so coupling *
 *  terms like (h u)_x or g h_x resolve), then the shared driver       *
 *  integrates the whole vector.  Each function yields its own 2-D      *
 *  InterpolatingFunction, returned as {{u1 -> if1, u2 -> if2, ...}}.   *
 *                                                                     *
 *  Scope: one spatial dimension, temporal order 1 or 2 per function,  *
 *  each evolution equation already solvable for exactly one function's *
 *  highest temporal derivative (no coupled mass matrix), real-valued, *
 *  machine precision.  Spatial scheme: centered Fornberg (default) or  *
 *  an upwind option (Lax-Friedrichs / donor-cell) for advective terms. *
 * ================================================================== */

/* Spatial discretization scheme for first-derivative (advective) terms. */
typedef enum { ND_SCHEME_CENTERED = 0, ND_SCHEME_UPWIND_LF, ND_SCHEME_UPWIND_DONOR } NdScheme;

/* Per-function descriptor for a coupled system. */
typedef struct {
    const char* fname;
    bool  applied;
    int   torder;                       /* highest temporal order            */
    bool  has_s[ND_MAX_SORDER + 1];     /* which spatial orders appear        */
    int   sord[ND_MAX_SORDER];          /* distinct spatial orders (sorted)   */
    int   nsord;
    Expr** ic;                          /* [torder] initial data in x (owned) */
    Expr *aL, *bL, *rL;                 /* left BC  a u + b u_x + r == 0       */
    Expr *aR, *bR, *rR;                 /* right BC                            */
    bool  periodic;
    Expr* pde_eq;                       /* borrowed evolution eq (in eqitems)  */
    Expr* G;                            /* solved top temporal derivative      */
    double t0;
    size_t nunk;                        /* interior unknown nodes             */
    int    gbase;                       /* first unknown grid index (0 or 1)   */
    size_t base;                        /* offset of this block in state vec   */
    Expr*  bc_left; Expr* bc_right;     /* eliminated boundary node values     */
} NdSysFunc;

/* Detect a Method-selected upwind spatial scheme.  DifferenceOrder -> 1 (or an
 * explicit "Upwind"/"DonorCell" spatial method) selects donor-cell upwinding for
 * scalar-separable advection; "LaxFriedrichs" selects artificial-viscosity
 * stabilization (robust for nonlinear systems).  Centered otherwise. */
/* Recursively test for a boolean flag sub-option `name`: a rule `name -> True`
 * (or a bare string element "name" inside the Method/SpatialDiscretization
 * lists).  Coexists with the standard SpatialDiscretization grid list. */
static bool nd_scan_flag_subopt(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_STRING) return strcmp(e->data.string, name) == 0;
    if (e->type != EXPR_FUNCTION) return false;
    if ((head_is((Expr*)e, SYM_Rule) || head_is((Expr*)e, SYM_RuleDelayed))
        && e->data.function.arg_count == 2) {
        const Expr* lhs = e->data.function.args[0];
        const char* key = NULL;
        if (lhs->type == EXPR_STRING) key = lhs->data.string;
        else if (lhs->type == EXPR_SYMBOL) key = lhs->data.symbol.name;
        if (key && strcmp(key, name) == 0) {
            const Expr* v = e->data.function.args[1];
            return !(v->type == EXPR_SYMBOL && v->data.symbol.name == intern_symbol("False"));
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (nd_scan_flag_subopt(e->data.function.args[i], name)) return true;
    return false;
}

/* First-derivative stencil at grid node `jgrid` for one function, honoring the
 * selected scheme.  Centered -> the shared Fornberg stencil.  Donor-cell upwind
 * -> a symbolic one-sided difference biased by the sign of `wind` (a constant
 * advection speed expr, borrowed) via Max[Sign[.],0]; the wind sign selects the
 * backward (wind>0) or forward (wind<0) neighbor.  Lax-Friedrichs is handled by
 * the caller as an added viscosity term, so it uses the centered stencil here. */
static Expr* nd_first_deriv_stencil(NdScheme scheme, const Expr* wind,
                                    int jgrid, int nx, bool periodic, int order,
                                    double h, int torder, Expr** ysym,
                                    Expr* bc_left, Expr* bc_right) {
    if (scheme != ND_SCHEME_UPWIND_DONOR)
        return nd_mol_stencil_expr(jgrid, nx, periodic, 1, order, h, torder, ysym,
                                   bc_left, bc_right);
    /* Donor-cell: backward = (u_j - u_{j-1})/h, forward = (u_{j+1} - u_j)/h. */
    Expr* uj  = nd_node_value(jgrid,     nx, periodic, torder, ysym, bc_left, bc_right);
    Expr* ujm = nd_node_value(jgrid - 1, nx, periodic, torder, ysym, bc_left, bc_right);
    Expr* ujp = nd_node_value(jgrid + 1, nx, periodic, torder, ysym, bc_left, bc_right);
    Expr* back = nd_scaled(1.0 / h, nd_plus2any(expr_copy(uj), nd_scaled(-1.0, ujm)));
    Expr* fwd  = nd_scaled(1.0 / h, nd_plus2any(ujp, nd_scaled(-1.0, expr_copy(uj))));
    expr_free(uj);
    /* weight_back = Max[Sign[wind], 0], weight_fwd = Max[-Sign[wind], 0]. */
    Expr* sgn  = nd_call1(intern_symbol("Sign"), expr_copy((Expr*)wind));
    Expr* wb   = eval_and_free(nd_call2(intern_symbol("Max"), expr_copy(sgn), expr_new_integer(0)));
    Expr* wf   = eval_and_free(nd_call2(intern_symbol("Max"),
                               nd_scaled(-1.0, sgn), expr_new_integer(0)));
    Expr* tb   = nd_call2(SYM_Times, wb, back);
    Expr* tf   = nd_call2(SYM_Times, wf, fwd);
    return nd_plus2any(tb, tf);
}

static Expr* nd_mol_solve_system(Expr* res, const NdOpts* o0, const char* forced_method) {
    NdOpts o = *o0;
    if (forced_method && strcmp(forced_method, "MethodOfLines") != 0)
        o.method = intern_symbol(forced_method);
    o.spec = numeric_machine_spec(); o.wp_bits = 53;   /* systems: machine precision */
    NumericSpec spec = o.spec;

    Expr** A = res->data.function.args;
    size_t argc = res->data.function.arg_count;

    size_t pos_end = 2;
    while (pos_end < argc && !nd_mol_is_option(A[pos_end])) pos_end++;
    size_t nranges = (pos_end >= 2) ? pos_end - 2 : 0;
    if (nranges != 2) {
        if (nranges >= 3) nd_mol_warn("pdesys", "systems are supported in one spatial "
                                                "dimension only");
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

    /* ---- dependent functions ---- */
    Expr* funcs = A[1];
    if (!head_is(funcs, SYM_List) || funcs->data.function.arg_count < 2) return NULL;
    size_t nfun = funcs->data.function.arg_count;
    NdSysFunc* F = calloc(nfun, sizeof(NdSysFunc));
    bool parse_ok = true;
    for (size_t fi = 0; fi < nfun; fi++) {
        Expr* it = funcs->data.function.args[fi];
        if (it->type == EXPR_SYMBOL) { F[fi].fname = it->data.symbol.name; F[fi].applied = false; }
        else if (it->type == EXPR_FUNCTION && it->data.function.head->type == EXPR_SYMBOL
                 && it->data.function.arg_count == 2) {
            F[fi].fname = it->data.function.head->data.symbol.name; F[fi].applied = true;
        } else { parse_ok = false; break; }
    }
    if (!parse_ok) { free(F); return NULL; }

    /* ---- grid + difference order + scheme ---- */
    long nx_l = 25, dord_l = 4;
    bool compiled = true;
    NdScheme scheme = ND_SCHEME_CENTERED;
    const char* SYM_Compiled = intern_symbol("Compiled");
    for (size_t i = pos_end; i < argc; i++) {
        if (!nd_mol_is_option(A[i])) continue;
        const char* opt = A[i]->data.function.args[0]->data.symbol.name;
        if (opt == SYM_Method) {
            nd_scan_int_subopt(A[i], "MinPoints", &nx_l);
            long dv = dord_l;
            if (nd_scan_int_subopt(A[i], "DifferenceOrder", &dv)) {
                dord_l = dv;
                if (dv == 1) scheme = ND_SCHEME_UPWIND_LF;   /* first-order = LF upwinding */
            }
            /* Systems use Lax-Friedrichs-family upwinding (donor-cell needs a
             * per-characteristic wind that is ill-defined for coupled fluxes). */
            if (nd_scan_flag_subopt(A[i], "Upwind") || nd_scan_flag_subopt(A[i], "DonorCell") ||
                nd_scan_flag_subopt(A[i], "LaxFriedrichs") || nd_scan_flag_subopt(A[i], "Rusanov"))
                scheme = ND_SCHEME_UPWIND_LF;
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
    if (dord_l > (long)nx - 1) dord_l = (long)nx - 1;
    int dord = (int)dord_l;
    int sten_ord = dord;
    if (scheme != ND_SCHEME_CENTERED && sten_ord < 2) sten_ord = 2;   /* keep centered base */
    if (sten_ord > (long)nx - 1) sten_ord = (int)nx - 1;

    /* ---- normalize equations ---- */
    Expr* eqns = A[0];
    Expr** eqsrc; size_t neq;
    if (head_is(eqns, SYM_List)) { neq = eqns->data.function.arg_count; eqsrc = eqns->data.function.args; }
    else { neq = 1; eqsrc = &eqns; }
    Expr** eqitems = malloc(sizeof(Expr*) * neq);
    for (size_t e = 0; e < neq; e++) eqitems[e] = nd_normalize_derivs(eqsrc[e]);

    /* ---- per-function order scan (over every equation) ---- */
    bool sys_unsupported = false;
    for (size_t fi = 0; fi < nfun; fi++) {
        int torder = 0; bool uns = false;
        for (size_t e = 0; e < neq; e++)
            nd_scan_orders(eqitems[e], F[fi].fname, &torder, F[fi].has_s, &uns);
        F[fi].torder = torder;
        F[fi].nsord = 0;
        for (int b = 1; b <= ND_MAX_SORDER; b++) if (F[fi].has_s[b]) F[fi].sord[F[fi].nsord++] = b;
        if (uns || torder < 1 || torder > 2) sys_unsupported = true;
        F[fi].ic = calloc((size_t)(torder > 0 ? torder : 1), sizeof(Expr*));
        F[fi].t0 = tmin;
    }
    if (sys_unsupported) {
        nd_mol_warn("pdeord", "each function must have temporal order 1 or 2 with no "
                              "mixed space-time derivatives");
        goto cleanup_early;
    }

    /* ---- classify every equation ---- */
    for (size_t e = 0; e < neq; e++) {
        Expr* eq = eqitems[e];
        if (!head_is(eq, SYM_Equal) || eq->data.function.arg_count != 2) continue;
        Expr* L = eq->data.function.args[0];
        Expr* R = eq->data.function.args[1];

        /* scan the fapps of each function to categorize this equation */
        bool any_tx = false, any_ic = false, any_bc = false;
        int ic_fi = -1, ic_ord = 0; double ic_t0 = tmin;
        int bc_fi = -1;
        for (size_t fi = 0; fi < nfun; fi++) {
            NdFApp fa[16]; int nfa = 0;
            nd_collect_fapps(eq, F[fi].fname, fa, &nfa, 16);
            for (int i = 0; i < nfa; i++) {
                bool a0_t = fa[i].arg0->type == EXPR_SYMBOL && fa[i].arg0->data.symbol.name == tvar;
                bool a1_x = fa[i].arg1->type == EXPR_SYMBOL && fa[i].arg1->data.symbol.name == xvar;
                double n0, n1;
                bool a0_num = nd_eval_to_double((Expr*)fa[i].arg0, spec, &n0);
                bool a1_num = nd_eval_to_double((Expr*)fa[i].arg1, spec, &n1);
                if (a0_t && a1_x) any_tx = true;
                else if (a0_num && a1_x) { any_ic = true; ic_fi = (int)fi; ic_ord = fa[i].a; ic_t0 = n0; }
                else if (a0_t && a1_num) { any_bc = true; bc_fi = (int)fi; }
            }
        }

        if (any_ic && ic_fi >= 0) {
            const char* fn = F[ic_fi].fname;
            int a, b; const Expr *g0, *g1; Expr* val = NULL;
            if (nd_pde_match(L, fn, &a, &b, &g0, &g1) && !nd_contains_func(R, fn)) val = R;
            else if (nd_pde_match(R, fn, &a, &b, &g0, &g1) && !nd_contains_func(L, fn)) val = L;
            if (val && ic_ord < F[ic_fi].torder && !F[ic_fi].ic[ic_ord]) {
                F[ic_fi].ic[ic_ord] = expr_copy(val);
                F[ic_fi].t0 = ic_t0;
            }
            continue;
        }

        if (any_bc && bc_fi >= 0) {
            const char* fn = F[bc_fi].fname;
            /* periodic: fn[t,c1] == fn[t,c2] with distinct boundaries */
            int la, lb, ra, rb; const Expr *l0, *l1, *r0, *r1;
            bool Lm = nd_pde_match(L, fn, &la, &lb, &l0, &l1);
            bool Rm = nd_pde_match(R, fn, &ra, &rb, &r0, &r1);
            if (Lm && Rm && la == 0 && lb == 0 && ra == 0 && rb == 0) {
                double x1, x2;
                if (nd_eval_to_double((Expr*)l1, spec, &x1) &&
                    nd_eval_to_double((Expr*)r1, spec, &x2) && fabs(x1 - x2) > 1e-9) {
                    F[bc_fi].periodic = true; continue;
                }
            }
            /* general linear BC: recover a,b,r at the boundary point */
            const Expr* bnd_pt = NULL; double bnd_x = 0.0;
            NdFApp fa[16]; int nfa = 0; nd_collect_fapps(eq, fn, fa, &nfa, 16);
            for (int i = 0; i < nfa; i++) {
                double n1;
                if (fa[i].arg0->type == EXPR_SYMBOL && fa[i].arg0->data.symbol.name == tvar
                    && nd_eval_to_double((Expr*)fa[i].arg1, spec, &n1)) { bnd_pt = fa[i].arg1; bnd_x = n1; }
            }
            if (!bnd_pt) continue;
            Expr* sub0 = nd_call2(SYM_Subtract, expr_copy(L), expr_copy(R));
            Expr* Ulit = nd_bc_lit(fn, 0, tvar, bnd_pt);
            Expr* Uxlit = nd_bc_lit(fn, 1, tvar, bnd_pt);
            Expr* Us = expr_new_symbol("NDSolve`Ub");
            Expr* Uxs = expr_new_symbol("NDSolve`Uxb");
            Expr* subL[2] = { Ulit, Uxlit }; Expr* subR[2] = { Us, Uxs };
            Expr* sub2 = nd_replace_all(sub0, subL, subR, 2);
            Expr* da[2] = { expr_copy(sub2), expr_copy(Us) };
            Expr* aco = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), da, 2));
            Expr* db[2] = { expr_copy(sub2), expr_copy(Uxs) };
            Expr* bco = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), db, 2));
            Expr* z0 = expr_new_integer(0), *z1 = expr_new_integer(0);
            Expr* zl[2] = { Us, Uxs }; Expr* zs[2] = { z0, z1 };
            Expr* rco = nd_replace_all(expr_copy(sub2), zl, zs, 2);
            expr_free(z0); expr_free(z1);
            expr_free(Ulit); expr_free(Uxlit); expr_free(Us); expr_free(Uxs); expr_free(sub2);
            bool left = fabs(bnd_x - xmin) <= 1e-9 * (fabs(xmin) + 1.0);
            if (left && !F[bc_fi].aL) { F[bc_fi].aL = aco; F[bc_fi].bL = bco; F[bc_fi].rL = rco; }
            else if (!left && !F[bc_fi].aR) { F[bc_fi].aR = aco; F[bc_fi].bR = bco; F[bc_fi].rR = rco; }
            else { expr_free(aco); expr_free(bco); expr_free(rco); }
            continue;
        }

        if (any_tx) {
            /* interior evolution equation: owner = the single function whose
             * highest temporal derivative appears here. */
            int owner = -1, nowner = 0;
            for (size_t fi = 0; fi < nfun; fi++) {
                NdFApp fa[16]; int nfa = 0; nd_collect_fapps(eq, F[fi].fname, fa, &nfa, 16);
                for (int i = 0; i < nfa; i++) {
                    bool a0_t = fa[i].arg0->type == EXPR_SYMBOL && fa[i].arg0->data.symbol.name == tvar;
                    bool a1_x = fa[i].arg1->type == EXPR_SYMBOL && fa[i].arg1->data.symbol.name == xvar;
                    if (a0_t && a1_x && fa[i].a == F[fi].torder && fa[i].b == 0) {
                        owner = (int)fi; nowner++; break;
                    }
                }
            }
            if (nowner == 1 && !F[owner].pde_eq) F[owner].pde_eq = eq;
            else if (nowner != 1) sys_unsupported = true;   /* coupled mass matrix */
        }
    }

    if (sys_unsupported) {
        nd_mol_warn("pdesys", "each evolution equation must be solvable for exactly one "
                              "function's highest temporal derivative");
        goto cleanup_early;
    }

    /* ---- validate completeness per function ---- */
    for (size_t fi = 0; fi < nfun; fi++) {
        bool have_bc = F[fi].periodic || (F[fi].aL && F[fi].aR);
        bool ok = F[fi].pde_eq && have_bc;
        for (int a = 0; a < F[fi].torder && ok; a++) if (!F[fi].ic[a]) ok = false;
        if (!ok) { sys_unsupported = true; break; }
    }
    if (sys_unsupported) {
        nd_mol_warn("pdeic", "each function needs an evolution equation, a boundary "
                             "condition at each end (or periodicity), and full initial data");
        goto cleanup_early;
    }

    /* ---- solve each evolution equation for its top temporal derivative ---- */
    bool solve_ok = true, solve_warned = false;
    for (size_t fi = 0; fi < nfun && solve_ok; fi++) {
        Expr* pde_eq = F[fi].pde_eq;
        Expr* R0 = nd_call2(SYM_Subtract, expr_copy(pde_eq->data.function.args[0]),
                                          expr_copy(pde_eq->data.function.args[1]));
        Expr* topLit = nd_pde_lit(F[fi].fname, F[fi].torder, 0, tvar, xvar);
        Expr* Psym = expr_new_symbol("NDSolve`Pt");
        Expr* Rp = nd_replace_all(R0, &topLit, &Psym, 1);
        Expr* dargs[2] = { expr_copy(Rp), expr_copy(Psym) };
        Expr* aE = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), dargs, 2));
        Expr* zero = expr_new_integer(0);
        Expr* bE = nd_replace_all(expr_copy(Rp), &Psym, &zero, 1);
        expr_free(zero);
        Expr* invE = nd_call2(SYM_Power, expr_copy(aE), expr_new_integer(-1));
        Expr* g3[3] = { expr_new_integer(-1), expr_copy(bE), invE };
        F[fi].G = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), g3, 3));
        expr_free(topLit); expr_free(Psym); expr_free(Rp); expr_free(aE); expr_free(bE);
        if (!F[fi].G) solve_ok = false;
        else if (nd_has_imaginary(F[fi].G)) {
            nd_mol_warn("pdecplx", "complex-valued systems are not supported");
            solve_ok = false; solve_warned = true;
        }
    }
    if (!solve_ok) {
        if (!solve_warned) nd_mol_warn("pdesolve", "could not solve an evolution "
                                                   "equation for its temporal derivative");
        goto cleanup_solved;
    }

    /* ---- global block-per-function reduced-state layout ---- */
    size_t d = 0;
    for (size_t fi = 0; fi < nfun; fi++) {
        F[fi].nunk  = F[fi].periodic ? (nx - 1) : (nx - 2);
        F[fi].gbase = F[fi].periodic ? 0 : 1;
        F[fi].base  = d;
        d += F[fi].nunk * (size_t)F[fi].torder;
    }
    Expr** ysym = malloc(sizeof(Expr*) * d);
    for (size_t idx = 0; idx < d; idx++) {
        char buf[48]; snprintf(buf, sizeof(buf), "NDSolve`w%zu", idx);
        ysym[idx] = expr_new_symbol(intern_symbol(buf));
    }

    /* ---- eliminate boundary nodes per function ---- */
    for (size_t fi = 0; fi < nfun; fi++) {
        if (F[fi].periodic) continue;
        F[fi].bc_left  = nd_bc_make_value(F[fi].aL, F[fi].bL, F[fi].rL, 0, (int)nx,
                                          sten_ord, h, F[fi].torder, &ysym[F[fi].base]);
        F[fi].bc_right = nd_bc_make_value(F[fi].aR, F[fi].bR, F[fi].rR, (int)nx - 1, (int)nx,
                                          sten_ord, h, F[fi].torder, &ysym[F[fi].base]);
    }

    /* ---- assemble the global NdProblem ---- */
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

    /* Precompute the substitution literal patterns (shared across all nodes):
     * for every function g -> its value, each spatial order, each lower temporal
     * order; plus the spatial coordinate x. */
    typedef struct { int kind; size_t g; int p; } SubSpec;  /* kind:0 val 1 spat 2 time 3 x */
    size_t maxlits = 1;
    for (size_t g = 0; g < nfun; g++)
        maxlits += 1 + (size_t)F[g].nsord + (size_t)(F[g].torder - 1);
    Expr** lits = malloc(sizeof(Expr*) * maxlits);
    SubSpec* spec_of = malloc(sizeof(SubSpec) * maxlits);
    size_t nlits = 0;
    for (size_t g = 0; g < nfun; g++) {
        lits[nlits] = nd_pde_lit(F[g].fname, 0, 0, tvar, xvar);
        spec_of[nlits] = (SubSpec){ 0, g, 0 }; nlits++;
        for (int i = 0; i < F[g].nsord; i++) {
            lits[nlits] = nd_pde_lit(F[g].fname, 0, F[g].sord[i], tvar, xvar);
            spec_of[nlits] = (SubSpec){ 1, g, F[g].sord[i] }; nlits++;
        }
        for (int m = 1; m < F[g].torder; m++) {
            lits[nlits] = nd_pde_lit(F[g].fname, m, 0, tvar, xvar);
            spec_of[nlits] = (SubSpec){ 2, g, m }; nlits++;
        }
    }
    lits[nlits] = expr_new_symbol(xvar);
    spec_of[nlits] = (SubSpec){ 3, 0, 0 }; nlits++;

    /* ---- build the per-node RHS for every function ---- */
    bool build_ok = true;
    for (size_t fi = 0; fi < nfun && build_ok; fi++) {
        int torder_f = F[fi].torder;
        /* Lax-Friedrichs viscosity coefficient (grid-scaled): applied to purely
         * hyperbolic functions (first spatial order present, no diffusion). */
        bool hyperbolic = F[fi].has_s[1] && !F[fi].has_s[2];
        double lf_nu = (scheme == ND_SCHEME_UPWIND_LF && hyperbolic) ? (0.5 * h) : 0.0;
        for (size_t u = 0; u < F[fi].nunk; u++) {
            size_t base = F[fi].base + u * (size_t)torder_f;
            int jgrid = (int)u + F[fi].gbase;
            double xj = xmin + (double)jgrid * h;
            for (int m = 0; m + 1 < torder_f; m++)
                P.f[base + (size_t)m] = expr_copy(ysym[base + (size_t)m + 1]);
            /* build substitutions for this node across all functions */
            Expr** subs = malloc(sizeof(Expr*) * nlits);
            for (size_t s = 0; s < nlits; s++) {
                SubSpec sp = spec_of[s];
                size_t g = sp.g;
                Expr** ys = &ysym[F[g].base];
                if (sp.kind == 0) {
                    subs[s] = nd_node_value(jgrid, (int)nx, F[g].periodic, F[g].torder, ys,
                                            F[g].bc_left, F[g].bc_right);
                } else if (sp.kind == 1) {
                    if (sp.p == 1)
                        subs[s] = nd_first_deriv_stencil(scheme, ys[0] /*unused for centered*/,
                                       jgrid, (int)nx, F[g].periodic, sten_ord, h,
                                       F[g].torder, ys, F[g].bc_left, F[g].bc_right);
                    else
                        subs[s] = nd_mol_stencil_expr(jgrid, (int)nx, F[g].periodic, sp.p,
                                       sten_ord, h, F[g].torder, ys, F[g].bc_left, F[g].bc_right);
                } else if (sp.kind == 2) {
                    int gnode = jgrid - F[g].gbase;
                    if (gnode >= 0 && gnode < (int)F[g].nunk)
                        subs[s] = expr_copy(ys[(size_t)gnode * (size_t)F[g].torder + (size_t)sp.p]);
                    else
                        subs[s] = expr_new_integer(0);
                } else {
                    subs[s] = expr_new_real(xj);
                }
            }
            Expr* rhs = nd_replace_all(expr_copy(F[fi].G), lits, subs, nlits);
            for (size_t s = 0; s < nlits; s++) expr_free(subs[s]);
            free(subs);
            /* Lax-Friedrichs artificial viscosity: + nu * u_xx (centered). */
            if (lf_nu > 0.0 && rhs) {
                Expr* visc = nd_mol_stencil_expr(jgrid, (int)nx, F[fi].periodic, 2, sten_ord,
                                                 h, torder_f, &ysym[F[fi].base],
                                                 F[fi].bc_left, F[fi].bc_right);
                rhs = nd_plus2any(rhs, nd_scaled(lf_nu, visc));
            }
            P.f[base + (size_t)torder_f - 1] = rhs;
            if (!rhs) build_ok = false;
        }
    }
    for (size_t s = 0; s < nlits; s++) expr_free(lits[s]);
    free(lits); free(spec_of);
    for (size_t i = 0; i < d; i++) if (!P.f[i]) build_ok = false;

    /* ---- initial data ---- */
    for (size_t fi = 0; fi < nfun && build_ok; fi++) {
        for (size_t u = 0; u < F[fi].nunk; u++) {
            int jgrid = (int)u + F[fi].gbase;
            double xj = xmin + (double)jgrid * h;
            for (int m = 0; m < F[fi].torder; m++) {
                double v;
                if (!nd_eval_at(F[fi].ic[m], xvar, xj, spec, &v)) { build_ok = false; break; }
                P.Y0[F[fi].base + u * (size_t)F[fi].torder + (size_t)m] = v;
            }
        }
    }
    P.t0 = tmin; P.tmin = tmin; P.tmax = tmax;

    /* ---- compile a linear operator if the whole system is linear ---- */
    if (build_ok && compiled) P.op = nd_operator_try_build(&P);

    /* parabolic (stiff) if any function is diffusion-driven at temporal order 1 */
    bool parabolic = false;
    for (size_t fi = 0; fi < nfun; fi++)
        if (F[fi].torder == 1 && (F[fi].has_s[2] || F[fi].has_s[4] || F[fi].has_s[6] || F[fi].has_s[8]))
            parabolic = true;

    /* ---- integrate ---- */
    Expr* result = NULL;
    if (build_ok) {
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

        /* one InterpolatingFunction per function -> {{u1->if1, u2->if2, ...}} */
        if (sol.n >= 2) {
            Expr** rules = malloc(sizeof(Expr*) * nfun);
            bool res_ok = true;
            for (size_t fi = 0; fi < nfun && res_ok; fi++) {
                size_t npts = sol.n * nx;
                Expr** entries = malloc(sizeof(Expr*) * npts);
                size_t pp = 0;
                for (size_t i = 0; i < sol.n; i++) {
                    double ti = sol.ts[i];
                    nd_bind_set(&P.bind_t, expr_new_real(ti));
                    for (size_t k = 0; k < d; k++) nd_bind_set(&P.bind_y[k], expr_new_real(sol.Ys[i * d + k]));
                    for (size_t g = 0; g < nx; g++) {
                        double xg = xmin + (double)g * h;
                        double val;
                        if (F[fi].periodic) {
                            size_t gg = (g == nx - 1) ? 0 : g;
                            val = sol.Ys[i * d + F[fi].base + gg * (size_t)F[fi].torder + 0];
                        } else if (g == 0) {
                            if (!nd_bound_eval(F[fi].bc_left, spec, &val)) val = 0.0;
                        } else if (g == nx - 1) {
                            if (!nd_bound_eval(F[fi].bc_right, spec, &val)) val = 0.0;
                        } else {
                            val = sol.Ys[i * d + F[fi].base + (g - 1) * (size_t)F[fi].torder + 0];
                        }
                        Expr* coord[2] = { expr_new_real(ti), expr_new_real(xg) };
                        Expr* coordL = expr_new_function(expr_new_symbol(SYM_List), coord, 2);
                        Expr* pair[2] = { coordL, expr_new_real(val) };
                        entries[pp++] = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
                    }
                }
                Expr* data = expr_new_function(expr_new_symbol(SYM_List), entries, npts);
                free(entries);
                Expr* ifun = eval_and_free(expr_new_function(expr_new_symbol(SYM_Interpolation), &data, 1));
                if (!head_is(ifun, SYM_InterpolatingFunction)) { expr_free(ifun); res_ok = false; break; }
                Expr* lhs;
                if (F[fi].applied) {
                    Expr* ar[2] = { expr_new_symbol(tvar), expr_new_symbol(xvar) };
                    lhs = expr_new_function(expr_new_symbol(F[fi].fname), ar, 2);
                    Expr* ar2[2] = { expr_new_symbol(tvar), expr_new_symbol(xvar) };
                    ifun = expr_new_function(ifun, ar2, 2);
                } else lhs = expr_new_symbol(F[fi].fname);
                rules[fi] = nd_call2(SYM_Rule, lhs, ifun);
            }
            if (res_ok) {
                Expr* inner = expr_new_function(expr_new_symbol(SYM_List), rules, nfun);
                result = expr_new_function(expr_new_symbol(SYM_List), &inner, 1);
            }
            free(rules);
        }
        nd_solution_free(&sol);
    }
    nd_bind_restore(&P.bind_t);
    for (size_t i = 0; i < d; i++) nd_bind_restore(&P.bind_y[i]);

    nd_operator_free(P.op);
    for (size_t i = 0; i < d; i++) { expr_free(P.f[i]); expr_free(P.ysym[i]); }
    free(P.f); free(P.ysym); free(P.Y0); free(P.bind_y);
    if (P.jac) {
        for (size_t i = 0; i < d; i++) { for (size_t j = 0; j < d; j++) expr_free(P.jac[i][j]); free(P.jac[i]); }
        free(P.jac);
    }

    for (size_t fi = 0; fi < nfun; fi++) { expr_free(F[fi].bc_left); expr_free(F[fi].bc_right); }
    /* fallthrough cleanup below */
    for (size_t fi = 0; fi < nfun; fi++) {
        expr_free(F[fi].G);
        for (int a = 0; a < F[fi].torder; a++) expr_free(F[fi].ic[a]);
        free(F[fi].ic);
        expr_free(F[fi].aL); expr_free(F[fi].bL); expr_free(F[fi].rL);
        expr_free(F[fi].aR); expr_free(F[fi].bR); expr_free(F[fi].rR);
    }
    for (size_t e = 0; e < neq; e++) expr_free(eqitems[e]);
    free(eqitems);
    free(F);
    return result;

cleanup_solved:
    for (size_t fi = 0; fi < nfun; fi++) expr_free(F[fi].G);
cleanup_early:
    for (size_t fi = 0; fi < nfun; fi++) {
        if (F[fi].ic) { for (int a = 0; a < F[fi].torder; a++) expr_free(F[fi].ic[a]); free(F[fi].ic); }
        expr_free(F[fi].aL); expr_free(F[fi].bL); expr_free(F[fi].rL);
        expr_free(F[fi].aR); expr_free(F[fi].bR); expr_free(F[fi].rR);
    }
    for (size_t e = 0; e < neq; e++) expr_free(eqitems[e]);
    free(eqitems);
    free(F);
    return NULL;
}

/* ------------------------------------------------------------------ *
 *  Method of Lines solver (Phase 1)                                   *
 * ------------------------------------------------------------------ */
Expr* nd_mol_solve(Expr* res, const NdOpts* o0, const char* forced_method) {
    NdOpts o = *o0;
    if (forced_method && strcmp(forced_method, "MethodOfLines") != 0)
        o.method = intern_symbol(forced_method);
#ifdef USE_MPFR
    bool use_mpfr = (o0->spec.mode == NUMERIC_MODE_MPFR);
#else
    bool use_mpfr = false;
#endif
    if (!use_mpfr) { o.spec = numeric_machine_spec(); o.wp_bits = 53; }
    NumericSpec spec = o.spec;

    Expr** A = res->data.function.args;
    size_t argc = res->data.function.arg_count;

    /* ---- positional range count ---- */
    size_t pos_end = 2;
    while (pos_end < argc && !nd_mol_is_option(A[pos_end])) pos_end++;
    size_t nranges = (pos_end >= 2) ? pos_end - 2 : 0;
    if (nranges < 2) return NULL;
    if (nranges == 3) return nd_mol_solve_2d(res, o0, forced_method);   /* two spatial dims */
    if (nranges > 3) {
        nd_mol_warn("pdedim", "at most two spatial dimensions are supported");
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
        if (funcs->data.function.arg_count == 0) return NULL;
        if (funcs->data.function.arg_count > 1)
            return nd_mol_solve_system(res, o0, forced_method);   /* coupled system */
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
    NdScheme scheme = ND_SCHEME_CENTERED;   /* opt-in upwinding for advective terms */
    const char* SYM_Compiled = intern_symbol("Compiled");
    for (size_t i = pos_end; i < argc; i++) {
        if (!nd_mol_is_option(A[i])) continue;
        const char* opt = A[i]->data.function.args[0]->data.symbol.name;
        if (opt == SYM_Method) {
            nd_scan_int_subopt(A[i], "MinPoints", &nx_l);
            long dv = dord_l;
            if (nd_scan_int_subopt(A[i], "DifferenceOrder", &dv)) {
                dord_l = dv;
                if (dv == 1) scheme = ND_SCHEME_UPWIND_DONOR;   /* first-order = upwind */
            }
            if (nd_scan_flag_subopt(A[i], "Upwind") || nd_scan_flag_subopt(A[i], "DonorCell"))
                scheme = ND_SCHEME_UPWIND_DONOR;
            else if (nd_scan_flag_subopt(A[i], "LaxFriedrichs") || nd_scan_flag_subopt(A[i], "Rusanov"))
                scheme = ND_SCHEME_UPWIND_LF;
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
    /* Upwinding is realized by the donor-cell first derivative or the LF
     * viscosity term; the underlying centered stencils still want an even
     * accuracy order (DifferenceOrder->1 selects the scheme, not one-sided base). */
    if (scheme != ND_SCHEME_CENTERED && dord < 2) dord = 2;
    if (dord > (long)nx - 1) dord = (int)nx - 1;

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

    /* Complex-valued PDE (e.g. Schrödinger): the solved RHS carries the
     * imaginary unit.  Handled by Re/Im realification below (temporal order 1);
     * the realified system is real, integrated at machine precision. */
    bool cplx = nd_has_imaginary(G);
    if (cplx) { use_mpfr = false; o.spec = numeric_machine_spec(); o.wp_bits = 53; spec = o.spec; }
    if (cplx && torder != 1) {
        nd_mol_warn("pdecplx", "complex PDEs are supported only at temporal order 1");
        expr_free(G);
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

    /* Donor-cell upwind wind field: wind = -dG/d(u_x), so the upwind direction can
     * switch with the (possibly state-dependent) local advection speed at runtime. */
    Expr* wind_lit = NULL;
    if (scheme == ND_SCHEME_UPWIND_DONOR && has_s[1] && G) {
        Expr* lit_ux = nd_pde_lit(fname, 0, 1, tvar, xvar);
        Expr* Q = expr_new_symbol("NDSolve`Wq");
        Expr* Gq = nd_replace_all(expr_copy(G), &lit_ux, &Q, 1);
        Expr* da[2] = { Gq, expr_copy(Q) };
        Expr* coef = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), da, 2));
        wind_lit = coef ? eval_and_free(nd_scaled(-1.0, coef)) : NULL;
        expr_free(lit_ux); expr_free(Q);
    }
    /* Lax-Friedrichs artificial-viscosity coefficient (grid-scaled), applied to a
     * purely hyperbolic (advective, non-diffusive) scalar equation. */
    double lf_nu = (scheme == ND_SCHEME_UPWIND_LF && has_s[1] && !has_s[2]) ? (0.5 * h) : 0.0;

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
            if (scheme == ND_SCHEME_UPWIND_DONOR && sord[i] == 1 && wind_lit) {
                Expr* wl[2] = { lit_u, lit_xv };
                Expr* ws[2] = { expr_copy(ysym[base + 0]), expr_new_real(xj) };
                Expr* wind_node = nd_replace_all(expr_copy(wind_lit), wl, ws, 2);
                expr_free(ws[0]); expr_free(ws[1]);
                subs[cnt] = nd_first_deriv_stencil(scheme, wind_node, jgrid, (int)nx, periodic,
                                                   dord, h, torder, ysym, bc_left, bc_right);
                expr_free(wind_node);
            } else {
                subs[cnt] = nd_mol_stencil_expr(jgrid, (int)nx, periodic, sord[i], dord, h,
                                                torder, ysym, bc_left, bc_right);
            }
            cnt++;
        }
        for (int m = 1; m < torder; m++) {
            lits[cnt] = lit_ut[m]; subs[cnt] = expr_copy(ysym[base + (size_t)m]); cnt++;
        }
        lits[cnt] = lit_xv; subs[cnt] = expr_new_real(xj); cnt++;
        Expr* rhs = nd_replace_all(expr_copy(G), lits, subs, cnt);
        for (size_t s = 0; s < cnt; s++) expr_free(subs[s]);   /* lits are shared */
        if (lf_nu > 0.0 && rhs) {
            Expr* visc = nd_mol_stencil_expr(jgrid, (int)nx, periodic, 2, dord, h,
                                             torder, ysym, bc_left, bc_right);
            rhs = nd_plus2any(rhs, nd_scaled(lf_nu, visc));
        }
        P.f[base + (size_t)torder - 1] = rhs;
    }
    expr_free(wind_lit);
    expr_free(lit_u); expr_free(lit_xv);
    for (int i = 0; i < nsord; i++) expr_free(lit_sx[i]);
    free(lit_sx);
    for (int m = 1; m < torder; m++) expr_free(lit_ut[m]);
    free(lit_ut);

    bool build_ok = true;
    for (size_t i = 0; i < d; i++) if (!P.f[i]) build_ok = false;

    /* ---- initial data ---- */
    double* Y0re = NULL; double* Y0im = NULL;   /* complex path */
    if (cplx) { Y0re = malloc(sizeof(double) * nunk); Y0im = malloc(sizeof(double) * nunk); }
    for (size_t u = 0; u < nunk && build_ok; u++) {
        double xj = xmin + (double)((int)u + gbase) * h;
        if (cplx) {
            if (!nd_eval_complex_at(ic[0], xvar, xj, spec, &Y0re[u], &Y0im[u])) build_ok = false;
        } else {
            for (int m = 0; m < torder; m++) {
                double v;
                if (!nd_eval_at(ic[m], xvar, xj, spec, &v)) { build_ok = false; break; }
                P.Y0[u * (size_t)torder + (size_t)m] = v;
            }
        }
    }
    P.t0 = t0; P.tmin = tmin; P.tmax = tmax;

    /* ---- realify a complex system into interleaved (Re, Im) real unknowns ---- */
    if (build_ok && cplx) {
        if (!nd_realify(&P, nunk, Y0re, Y0im)) build_ok = false;
        d = P.d;
    }
    free(Y0re); free(Y0im);

    /* ---- compile a linear operator (fast RHS / exact Jacobian) if possible ---- */
    if (build_ok && compiled && !use_mpfr) P.op = nd_operator_try_build(&P);

    /* Parabolic problems (a diffusion term with first-order time evolution) are
     * stiff; default them to BDF when the user didn't force a method.  Complex
     * (Schrödinger) systems have imaginary spectra — keep the explicit default. */
    bool parabolic = (torder == 1) && !cplx &&
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
#ifdef USE_MPFR
        if (use_mpfr && !cplx) {
            /* Arbitrary precision: integrate at MPFR precision (explicit) and
             * assemble the MPFR 2-D InterpolatingFunction. */
            NdMolGrid grid = { fname, applied, tvar, xvar, xmin, h, nx, torder,
                               periodic, bc_left, bc_right };
            result = nd_solve_mpfr_mol(&P, &o, S, &grid);
        } else
#endif
        {
        NdSolution sol; nd_solution_init(&sol, d);
        NdStatus st = nd_integrate(&P, S, &o, &sol);
        if (st == ND_ERR_MAXSTEPS) nd_mol_warn("mxst", "maximum number of steps reached; returning partial solution");
        else if (st == ND_ERR_STEPSIZE) nd_mol_warn("ndsz", "step size effectively zero; stiffness suspected (try Method->\"BDF\")");
        else if (st == ND_ERR_NONCONV) nd_mol_warn("ndcf", "corrector failed to converge");
        else if (st == ND_ERR_SAMPLE) nd_mol_warn("nrnum", "spatial operator did not evaluate to a number");
        /* build the result while the bindings are still live (boundary/periodic
         * node values may depend on the reduced-state symbols). */
        if (cplx)
            result = nd_mol_build_result_complex(&P, fname, applied, tvar, xvar, xmin, h,
                                                 nx, nunk, gbase, periodic, bc_left, bc_right,
                                                 &sol, spec);
        else
            result = nd_mol_build_result(&P, fname, applied, tvar, xvar, xmin, h, nx, torder,
                                         periodic, bc_left, bc_right, &sol, spec);
        nd_solution_free(&sol);
        }
    }
    nd_bind_restore(&P.bind_t);
    for (size_t i = 0; i < d; i++) nd_bind_restore(&P.bind_y[i]);

    /* ---- cleanup (P.ysym/P.f may have been replaced by realification) ---- */
    nd_operator_free(P.op);
    for (size_t i = 0; i < d; i++) { expr_free(P.f[i]); expr_free(P.ysym[i]); }
    free(P.f); free(P.ysym); free(P.Y0); free(P.bind_y);
    (void)ysym;
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

/* ================================================================== *
 *  Two spatial dimensions (Phase 4)                                   *
 *                                                                     *
 *  u(t, x, y) on a rectangle, discretized on an nx*ny tensor grid.    *
 *  Interior nodes are the unknowns; each edge carries a general linear *
 *  boundary condition (Dirichlet/Neumann/Robin) whose boundary node    *
 *  values are eliminated into the interior stencils.  Pure x- and      *
 *  y-derivatives use Fornberg stencils along each axis; the result is  *
 *  a 3-D InterpolatingFunction applied as u[t, x, y].  Scope: linear   *
 *  boundary conditions, unmixed spatial derivatives.                   *
 * ================================================================== */

/* Match  u[t,x,y]  or  Derivative[a,b,c][u][t,x,y]. */
static bool nd_pde3_match(const Expr* e, const char* fname, int* a, int* b, int* c,
                          const Expr** g0, const Expr** g1, const Expr** g2) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 3) return false;
    const Expr* head = e->data.function.head;
    *g0 = e->data.function.args[0];
    *g1 = e->data.function.args[1];
    *g2 = e->data.function.args[2];
    if (head->type == EXPR_SYMBOL) {
        if (head->data.symbol.name != fname) return false;
        *a = 0; *b = 0; *c = 0; return true;
    }
    if (head->type == EXPR_FUNCTION && head->data.function.arg_count == 1
        && head->data.function.args[0]->type == EXPR_SYMBOL
        && head->data.function.args[0]->data.symbol.name == fname) {
        const Expr* d = head->data.function.head;
        if (d->type == EXPR_FUNCTION && d->data.function.arg_count == 3
            && d->data.function.head->type == EXPR_SYMBOL
            && d->data.function.head->data.symbol.name == SYM_Derivative
            && d->data.function.args[0]->type == EXPR_INTEGER
            && d->data.function.args[1]->type == EXPR_INTEGER
            && d->data.function.args[2]->type == EXPR_INTEGER) {
            *a = (int)d->data.function.args[0]->data.integer;
            *b = (int)d->data.function.args[1]->data.integer;
            *c = (int)d->data.function.args[2]->data.integer;
            return true;
        }
    }
    return false;
}

/* Build  Derivative[a,b,c][u][t,x,y]  (a==b==c==0 -> u[t,x,y]). */
static Expr* nd_pde3_lit(const char* fname, int a, int b, int c,
                         const char* tv, const char* xv, const char* yv) {
    if (a == 0 && b == 0 && c == 0) {
        Expr* ar[3] = { expr_new_symbol(tv), expr_new_symbol(xv), expr_new_symbol(yv) };
        return expr_new_function(expr_new_symbol(fname), ar, 3);
    }
    Expr* di[3] = { expr_new_integer(a), expr_new_integer(b), expr_new_integer(c) };
    Expr* dhead = expr_new_function(expr_new_symbol(SYM_Derivative), di, 3);
    Expr* uu[1] = { expr_new_symbol(fname) };
    Expr* d2 = expr_new_function(dhead, uu, 1);
    Expr* ar[3] = { expr_new_symbol(tv), expr_new_symbol(xv), expr_new_symbol(yv) };
    return expr_new_function(d2, ar, 3);
}

/* Build Derivative[0,dx,dy][u][t, ax, ay] (dx==dy==0 -> u[t,ax,ay]); ax,ay copied. */
static Expr* nd_bc3_lit(const char* fname, int dx, int dy, const char* tv,
                        const Expr* ax, const Expr* ay) {
    Expr* t = expr_new_symbol(tv);
    Expr* px = expr_copy((Expr*)ax);
    Expr* py = expr_copy((Expr*)ay);
    if (dx == 0 && dy == 0) {
        Expr* ar[3] = { t, px, py };
        return expr_new_function(expr_new_symbol(fname), ar, 3);
    }
    Expr* di[3] = { expr_new_integer(0), expr_new_integer(dx), expr_new_integer(dy) };
    Expr* dhead = expr_new_function(expr_new_symbol(SYM_Derivative), di, 3);
    Expr* uu[1] = { expr_new_symbol(fname) };
    Expr* d2 = expr_new_function(dhead, uu, 1);
    Expr* ar[3] = { t, px, py };
    return expr_new_function(d2, ar, 3);
}

/* One occurrence of the dependent function in a 3-argument (t,x,y) equation. */
typedef struct { int a, b, c; const Expr *g0, *g1, *g2; } NdFApp3;
static void nd_collect_fapps3(const Expr* e, const char* fname,
                              NdFApp3* out, int* n, int max) {
    if (!e) return;
    int a, b, c; const Expr *g0, *g1, *g2;
    if (nd_pde3_match(e, fname, &a, &b, &c, &g0, &g1, &g2) && *n < max) {
        out[*n].a = a; out[*n].b = b; out[*n].c = c;
        out[*n].g0 = g0; out[*n].g1 = g1; out[*n].g2 = g2; (*n)++;
    }
    if (e->type == EXPR_FUNCTION) {
        nd_collect_fapps3(e->data.function.head, fname, out, n, max);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            nd_collect_fapps3(e->data.function.args[i], fname, out, n, max);
    }
}

static void nd_scan3(const Expr* e, const char* fname, int* torder,
                     int (*sp)[2], int* nsp, bool* unsupported) {
    if (!e) return;
    int a, b, c; const Expr *g0, *g1, *g2;
    if (nd_pde3_match(e, fname, &a, &b, &c, &g0, &g1, &g2)) {
        if (a > *torder) *torder = a;
        if (a == 0 && (b > 0 || c > 0)) {
            if (b > 0 && c > 0) *unsupported = true;          /* mixed x-y */
            else if (b > ND_MAX_SORDER || c > ND_MAX_SORDER) *unsupported = true;
            else {
                bool found = false;
                for (int i = 0; i < *nsp; i++) if (sp[i][0] == b && sp[i][1] == c) found = true;
                if (!found && *nsp < 2 * ND_MAX_SORDER) { sp[*nsp][0] = b; sp[*nsp][1] = c; (*nsp)++; }
            }
        }
        if (a > 0 && (b > 0 || c > 0)) *unsupported = true;   /* mixed space-time */
    }
    if (e->type == EXPR_FUNCTION) {
        nd_scan3(e->data.function.head, fname, torder, sp, nsp, unsupported);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            nd_scan3(e->data.function.args[i], fname, torder, sp, nsp, unsupported);
    }
}

/* 2-D grid context passed to the node-value / stencil builders. */
typedef struct {
    size_t nx, ny;
    double xmin, ymin, hx, hy;
    int torder;
    Expr** ysym;
    const char *tv, *xv, *yv;
    /* Each edge carries a general linear BC  a*u + b*u_n + r == 0  (n = the
     * outward-normal axis), coefficients as exprs in (t, tangential coord).
     * Dirichlet is the b==0 case; Neumann is a==0; Robin has both. */
    Expr *axl, *bxl, *rxl;   /* x = xmin edge (exprs in t, y) */
    Expr *axr, *bxr, *rxr;   /* x = xmax edge                */
    Expr *ayl, *byl, *ryl;   /* y = ymin edge (exprs in t, x) */
    Expr *ayr, *byr, *ryr;   /* y = ymax edge                */
    int dord;                /* difference order for the BC one-sided stencils */
    NumericSpec spec;
} Nd2D;

static Expr* nd_nodeval_2d(const Nd2D* G, int ix, int iy);

/* Eliminate a boundary node on `axis` (0=x, 1=y) at grid index `jbound` (0 or
 * n-1), transverse index `tj`, using the edge's linear BC  a*u + b*u_n + r == 0.
 * A one-sided first-derivative stencil gives
 *   u_n(jbound) = w0*u_bound + sum_k w_k*u_k,   so
 *   u_bound = -(r + b*sum_k w_k*u_k) / (a + b*w0).
 * a,b,r are exprs in (t, tangential coord); the tangential symbol is replaced by
 * its coordinate first.  Interior stencil neighbours resolve through
 * nd_nodeval_2d — an interior reduced-state symbol, or, when tj lands on a
 * boundary (a corner), the transverse edge (which references only interior
 * nodes, so the recursion is at most two deep). */
static Expr* nd_bc_eliminate_2d(const Nd2D* G, int axis, int jbound, int tj,
                                const Expr* a, const Expr* b, const Expr* r) {
    int n = axis == 0 ? (int)G->nx : (int)G->ny;
    double h = axis == 0 ? G->hx : G->hy;
    const char* tanv = axis == 0 ? G->yv : G->xv;
    double tanc = axis == 0 ? G->ymin + (double)tj * G->hy
                            : G->xmin + (double)tj * G->hx;
    Expr* tlit = expr_new_symbol(tanv);
    Expr* tsub = expr_new_real(tanc);
    Expr* as = eval_and_free(nd_replace_all(expr_copy((Expr*)a), &tlit, &tsub, 1));
    Expr* bs = eval_and_free(nd_replace_all(expr_copy((Expr*)b), &tlit, &tsub, 1));
    Expr* rs = eval_and_free(nd_replace_all(expr_copy((Expr*)r), &tlit, &tsub, 1));
    expr_free(tlit); expr_free(tsub);

    /* Dirichlet fast path (b == 0):  u_bound = -r/a. */
    if (bs->type == EXPR_INTEGER && bs->data.integer == 0) {
        expr_free(bs);
        Expr* inv = nd_call2(SYM_Power, as, expr_new_integer(-1));
        Expr* val = nd_scaled(-1.0, expr_new_function(expr_new_symbol(SYM_Times),
                              (Expr*[]){ rs, inv }, 2));
        return eval_and_free(val);
    }

    int idx[64]; double w[64]; int npts;
    nd_stencil_build(jbound, n, 1, G->dord, h, idx, w, &npts);
    double w0 = 0.0;
    Expr** terms = malloc(sizeof(Expr*) * (size_t)npts);
    int nt = 0;
    for (int k = 0; k < npts; k++) {
        if (idx[k] == jbound) { w0 = w[k]; continue; }
        if (w[k] == 0.0) continue;
        if (idx[k] <= 0 || idx[k] >= n - 1) continue;   /* skip same-axis edges */
        Expr* val = axis == 0 ? nd_nodeval_2d(G, idx[k], tj)
                              : nd_nodeval_2d(G, tj, idx[k]);
        terms[nt++] = nd_scaled(w[k], val);
    }
    Expr* un_rest;
    if (nt == 0) { free(terms); un_rest = expr_new_integer(0); }
    else if (nt == 1) { un_rest = terms[0]; free(terms); }
    else { un_rest = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)nt); free(terms); }

    Expr* denom = nd_plus2any(expr_copy(as), nd_scaled(w0, expr_copy(bs)));
    Expr* numer = nd_scaled(-1.0, nd_plus2any(expr_copy(rs),
                            expr_new_function(expr_new_symbol(SYM_Times),
                                (Expr*[]){ expr_copy(bs), un_rest }, 2)));
    Expr* inv = nd_call2(SYM_Power, denom, expr_new_integer(-1));
    Expr* val = expr_new_function(expr_new_symbol(SYM_Times),
                                  (Expr*[]){ numer, inv }, 2);
    expr_free(as); expr_free(bs); expr_free(rs);
    return eval_and_free(val);
}

/* Value at grid node (ix,iy): interior reduced-state symbol, or a boundary node
 * eliminated through its edge's general linear BC. */
static Expr* nd_nodeval_2d(const Nd2D* G, int ix, int iy) {
    int nx = (int)G->nx, ny = (int)G->ny;
    if (ix >= 1 && ix <= nx - 2 && iy >= 1 && iy <= ny - 2) {
        size_t iu = (size_t)(iy - 1) * (size_t)(nx - 2) + (size_t)(ix - 1);
        return expr_copy(G->ysym[iu * (size_t)G->torder + 0]);
    }
    /* x-edges take priority so corners resolve via the transverse (y) edge. */
    if (ix <= 0)       return nd_bc_eliminate_2d(G, 0, 0,      iy, G->axl, G->bxl, G->rxl);
    if (ix >= nx - 1)  return nd_bc_eliminate_2d(G, 0, nx - 1, iy, G->axr, G->bxr, G->rxr);
    if (iy <= 0)       return nd_bc_eliminate_2d(G, 1, 0,      ix, G->ayl, G->byl, G->ryl);
    return nd_bc_eliminate_2d(G, 1, ny - 1, ix, G->ayr, G->byr, G->ryr);
}

/* Stencil for the `deriv`-th derivative along axis (0=x,1=y) at node (ix,iy). */
static Expr* nd_stencil2d(const Nd2D* G, int ix, int iy, int axis, int deriv, int order) {
    int idx[64]; double w[64]; int npts;
    int j = (axis == 0) ? ix : iy;
    int n = (axis == 0) ? (int)G->nx : (int)G->ny;
    double h = (axis == 0) ? G->hx : G->hy;
    nd_stencil_build(j, n, deriv, order, h, idx, w, &npts);
    Expr** terms = malloc(sizeof(Expr*) * (size_t)npts);
    int nt = 0;
    for (int k = 0; k < npts; k++) {
        if (w[k] == 0.0) continue;
        Expr* val = (axis == 0) ? nd_nodeval_2d(G, idx[k], iy)
                                : nd_nodeval_2d(G, ix, idx[k]);
        terms[nt++] = nd_scaled(w[k], val);
    }
    Expr* r;
    if (nt == 0) { free(terms); r = expr_new_integer(0); }
    else if (nt == 1) { r = terms[0]; free(terms); }
    else { r = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)nt); free(terms); }
    return r;
}

static Expr* nd_mol_solve_2d(Expr* res, const NdOpts* o0, const char* forced_method) {
    NdOpts o = *o0;
    if (forced_method && strcmp(forced_method, "MethodOfLines") != 0)
        o.method = intern_symbol(forced_method);
    o.spec = numeric_machine_spec(); o.wp_bits = 53;
    NumericSpec spec = o.spec;

    Expr** A = res->data.function.args;
    size_t argc = res->data.function.arg_count;
    size_t pos_end = 2;
    while (pos_end < argc && !nd_mol_is_option(A[pos_end])) pos_end++;

    Expr* rt = A[2]; Expr* rx = A[3]; Expr* ry = A[4];
    Expr* rng[3] = { rt, rx, ry };
    for (int i = 0; i < 3; i++)
        if (!head_is(rng[i], SYM_List) || rng[i]->data.function.arg_count != 3
            || rng[i]->data.function.args[0]->type != EXPR_SYMBOL) return NULL;
    const char* tv = rt->data.function.args[0]->data.symbol.name;
    const char* xv = rx->data.function.args[0]->data.symbol.name;
    const char* yv = ry->data.function.args[0]->data.symbol.name;
    double tmin, tmax, xmin, xmax, ymin, ymax;
    if (!nd_eval_to_double(rt->data.function.args[1], spec, &tmin)) return NULL;
    if (!nd_eval_to_double(rt->data.function.args[2], spec, &tmax)) return NULL;
    if (!nd_eval_to_double(rx->data.function.args[1], spec, &xmin)) return NULL;
    if (!nd_eval_to_double(rx->data.function.args[2], spec, &xmax)) return NULL;
    if (!nd_eval_to_double(ry->data.function.args[1], spec, &ymin)) return NULL;
    if (!nd_eval_to_double(ry->data.function.args[2], spec, &ymax)) return NULL;
    if (!(xmax > xmin) || !(ymax > ymin)) return NULL;

    /* dependent function */
    Expr* funcs = A[1]; Expr* fitem = funcs;
    if (head_is(funcs, SYM_List)) {
        if (funcs->data.function.arg_count != 1) return NULL;
        fitem = funcs->data.function.args[0];
    }
    const char* fname; bool applied;
    if (fitem->type == EXPR_SYMBOL) { fname = fitem->data.symbol.name; applied = false; }
    else if (fitem->type == EXPR_FUNCTION && fitem->data.function.head->type == EXPR_SYMBOL
             && fitem->data.function.arg_count == 3) {
        fname = fitem->data.function.head->data.symbol.name; applied = true;
    } else return NULL;

    /* grid + difference order */
    long nx_l = 15, ny_l = 15, dord_l = 4;
    for (size_t i = pos_end; i < argc; i++) {
        if (!nd_mol_is_option(A[i])) continue;
        if (A[i]->data.function.args[0]->data.symbol.name == SYM_Method) {
            nd_scan_int_subopt(A[i], "MinPoints", &nx_l);
            nd_scan_int_subopt(A[i], "MinPoints", &ny_l);
            nd_scan_int_subopt(A[i], "DifferenceOrder", &dord_l);
        }
    }
    if (nx_l < 5) nx_l = 5;
    if (nx_l > 400) nx_l = 400;
    if (ny_l < 5) ny_l = 5;
    if (ny_l > 400) ny_l = 400;
    size_t nx = (size_t)nx_l, ny = (size_t)ny_l;
    double hx = (xmax - xmin) / (double)(nx - 1), hy = (ymax - ymin) / (double)(ny - 1);
    if (dord_l < 1) dord_l = 1;
    if (dord_l > (long)nx - 1) dord_l = (long)nx - 1;
    if (dord_l > (long)ny - 1) dord_l = (long)ny - 1;
    int dord = (int)dord_l;

    /* normalize equations */
    Expr* eqns = A[0]; Expr** eqsrc; size_t neq;
    if (head_is(eqns, SYM_List)) { neq = eqns->data.function.arg_count; eqsrc = eqns->data.function.args; }
    else { neq = 1; eqsrc = &eqns; }
    Expr** eqitems = malloc(sizeof(Expr*) * neq);
    for (size_t e = 0; e < neq; e++) eqitems[e] = nd_normalize_derivs(eqsrc[e]);

    /* scan orders */
    int torder = 0, nsp = 0; bool unsupported = false;
    int sp[2 * ND_MAX_SORDER][2];
    for (size_t e = 0; e < neq; e++)
        nd_scan3(eqitems[e], fname, &torder, sp, &nsp, &unsupported);
    if (unsupported || torder < 1 || torder > 2) {
        nd_mol_warn("pdeord", "2-D PDEs support unmixed spatial derivatives and "
                              "temporal order 1 or 2");
        for (size_t e = 0; e < neq; e++) expr_free(eqitems[e]);
        free(eqitems); return NULL;
    }

    /* classify: evolution PDE / IC / general linear edge BCs.  Each edge BC is
     * stored as the coefficients of  a*u + b*u_n + r == 0  (Dirichlet: b=0;
     * Neumann: a=0; Robin: both), which the boundary-node elimination turns into
     * a value expression.  Edge order: 0=x-min, 1=x-max, 2=y-min, 3=y-max. */
    Expr* pde_eq = NULL;
    Expr** ic = calloc((size_t)torder, sizeof(Expr*));    /* ic[a]: value in x,y  */
    Expr* ea[4] = { NULL, NULL, NULL, NULL };  /* per-edge a coefficient */
    Expr* eb[4] = { NULL, NULL, NULL, NULL };  /* per-edge b coefficient */
    Expr* er[4] = { NULL, NULL, NULL, NULL };  /* per-edge r term         */
    bool eset[4] = { false, false, false, false };
    bool bc_unsupported = false, periodic_req = false;
    double t0 = tmin;
    for (size_t e = 0; e < neq; e++) {
        Expr* eq = eqitems[e];
        if (!head_is(eq, SYM_Equal) || eq->data.function.arg_count != 2) continue;
        Expr* L = eq->data.function.args[0];
        Expr* R = eq->data.function.args[1];

        NdFApp3 fa[24]; int nfa = 0;
        nd_collect_fapps3(eq, fname, fa, &nfa, 24);
        bool interior = false, is_ic = false; int ic_ord = 0;
        int edge = -1; const Expr* bpt = NULL; bool multiedge = false;
        for (int i = 0; i < nfa; i++) {
            bool a0_t = fa[i].g0->type == EXPR_SYMBOL && fa[i].g0->data.symbol.name == tv;
            bool a1_x = fa[i].g1->type == EXPR_SYMBOL && fa[i].g1->data.symbol.name == xv;
            bool a2_y = fa[i].g2->type == EXPR_SYMBOL && fa[i].g2->data.symbol.name == yv;
            double n0, n1, n2;
            bool a0_n = nd_eval_to_double((Expr*)fa[i].g0, spec, &n0);
            bool a1_n = nd_eval_to_double((Expr*)fa[i].g1, spec, &n1);
            bool a2_n = nd_eval_to_double((Expr*)fa[i].g2, spec, &n2);
            int ee = -1; const Expr* ep = NULL;
            if (a0_t && a1_x && a2_y) interior = true;
            else if (a0_n && a1_x && a2_y) { is_ic = true; ic_ord = fa[i].a; t0 = n0; }
            else if (a0_t && a1_n && a2_y) {                    /* x-edge */
                ep = fa[i].g1;
                if (fabs(n1 - xmin) <= 1e-9 * (fabs(xmin) + 1.0)) ee = 0;
                else if (fabs(n1 - xmax) <= 1e-9 * (fabs(xmax) + 1.0)) ee = 1;
            } else if (a0_t && a1_x && a2_n) {                  /* y-edge */
                ep = fa[i].g2;
                if (fabs(n2 - ymin) <= 1e-9 * (fabs(ymin) + 1.0)) ee = 2;
                else if (fabs(n2 - ymax) <= 1e-9 * (fabs(ymax) + 1.0)) ee = 3;
            }
            if (ee >= 0) {
                if (edge >= 0 && edge != ee) multiedge = true;  /* couples two edges */
                edge = ee; bpt = ep;
            }
        }

        if (interior) { if (!pde_eq) pde_eq = eq; continue; }

        /* a single condition referencing two distinct edges (e.g.
         * u[t,xmin,y]==u[t,xmax,y]) is a periodic coupling — not yet supported. */
        if (multiedge) { periodic_req = true; continue; }

        if (is_ic) {
            int a, b, c; const Expr *g0, *g1, *g2; Expr* val = NULL;
            if (nd_pde3_match(L, fname, &a, &b, &c, &g0, &g1, &g2) && !nd_contains_func(R, fname)) val = R;
            else if (nd_pde3_match(R, fname, &a, &b, &c, &g0, &g1, &g2) && !nd_contains_func(L, fname)) val = L;
            if (val && ic_ord < torder && !ic[ic_ord]) ic[ic_ord] = expr_copy(val);
            continue;
        }

        if (edge >= 0 && bpt) {
            /* extract a, b, r for  a*u + b*u_n + r == 0  by substituting the
             * boundary value u and normal derivative u_n with fresh symbols. */
            int axis = (edge < 2) ? 0 : 1;
            Expr* ax = axis == 0 ? expr_copy((Expr*)bpt) : expr_new_symbol(xv);
            Expr* ay = axis == 0 ? expr_new_symbol(yv)   : expr_copy((Expr*)bpt);
            Expr* Ulit  = nd_bc3_lit(fname, 0, 0, tv, ax, ay);
            Expr* Unlit = axis == 0 ? nd_bc3_lit(fname, 1, 0, tv, ax, ay)
                                    : nd_bc3_lit(fname, 0, 1, tv, ax, ay);
            expr_free(ax); expr_free(ay);
            Expr* Us  = expr_new_symbol("NDSolve`Ub");
            Expr* Uns = expr_new_symbol("NDSolve`Unb");
            Expr* subL[2] = { Ulit, Unlit }; Expr* subR[2] = { Us, Uns };
            Expr* rr = nd_call2(SYM_Subtract, expr_copy(L), expr_copy(R));
            Expr* rr2 = nd_replace_all(rr, subL, subR, 2);      /* consumes rr */
            if (nd_contains_func(rr2, fname)) bc_unsupported = true;
            Expr* da[2] = { expr_copy(rr2), expr_copy(Us) };
            Expr* aco = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), da, 2));
            Expr* db[2] = { expr_copy(rr2), expr_copy(Uns) };
            Expr* bco = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), db, 2));
            Expr* z0 = expr_new_integer(0), *z1 = expr_new_integer(0);
            Expr* zl[2] = { Us, Uns }; Expr* zs[2] = { z0, z1 };
            Expr* rco = eval_and_free(nd_replace_all(expr_copy(rr2), zl, zs, 2));
            expr_free(z0); expr_free(z1);
            expr_free(Ulit); expr_free(Unlit); expr_free(Us); expr_free(Uns); expr_free(rr2);
            if (!eset[edge]) { ea[edge] = aco; eb[edge] = bco; er[edge] = rco; eset[edge] = true; }
            else { expr_free(aco); expr_free(bco); expr_free(rco); }
            continue;
        }
        if (!pde_eq) pde_eq = eq;
    }
    bool ok = pde_eq && eset[0] && eset[1] && eset[2] && eset[3] && !bc_unsupported && !periodic_req;
    for (int a = 0; a < torder && ok; a++) if (!ic[a]) ok = false;
    if (!ok) {
        if (periodic_req)
            nd_mol_warn("pdebc", "periodic boundary conditions are not yet supported for 2-D PDEs");
        else if (bc_unsupported)
            nd_mol_warn("pdebc", "2-D PDE boundary conditions must be linear "
                                 "(Dirichlet/Neumann/Robin) with unmixed derivatives");
        else
            nd_mol_warn("pdeic", "2-D PDEs need a boundary condition on each of the "
                                 "four edges and full initial data");
        for (int a = 0; a < torder; a++) expr_free(ic[a]);
        free(ic);
        for (int i = 0; i < 4; i++) { expr_free(ea[i]); expr_free(eb[i]); expr_free(er[i]); }
        for (size_t e = 0; e < neq; e++) expr_free(eqitems[e]);
        free(eqitems);
        return NULL;
    }

    /* solve for the top temporal derivative -> G */
    Expr* R0 = nd_call2(SYM_Subtract, expr_copy(pde_eq->data.function.args[0]),
                                      expr_copy(pde_eq->data.function.args[1]));
    Expr* topLit = nd_pde3_lit(fname, torder, 0, 0, tv, xv, yv);
    Expr* Psym = expr_new_symbol("NDSolve`Pt");
    Expr* Rp = nd_replace_all(R0, &topLit, &Psym, 1);
    Expr* dargs[2] = { expr_copy(Rp), expr_copy(Psym) };
    Expr* aE = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), dargs, 2));
    Expr* zero = expr_new_integer(0);
    Expr* bE = nd_replace_all(expr_copy(Rp), &Psym, &zero, 1);
    expr_free(zero);
    Expr* invE = nd_call2(SYM_Power, expr_copy(aE), expr_new_integer(-1));
    Expr* g3[3] = { expr_new_integer(-1), expr_copy(bE), invE };
    Expr* Gexpr = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), g3, 3));
    expr_free(topLit); expr_free(Psym); expr_free(Rp); expr_free(aE); expr_free(bE);

    size_t Nint = (nx - 2) * (ny - 2);
    size_t d = Nint * (size_t)torder;
    Expr** ysym = malloc(sizeof(Expr*) * d);
    for (size_t i = 0; i < d; i++) {
        char buf[48]; snprintf(buf, sizeof(buf), "NDSolve`w%zu", i);
        ysym[i] = expr_new_symbol(intern_symbol(buf));
    }

    Nd2D GC;
    GC.nx = nx; GC.ny = ny; GC.xmin = xmin; GC.ymin = ymin; GC.hx = hx; GC.hy = hy;
    GC.torder = torder; GC.ysym = ysym; GC.tv = tv; GC.xv = xv; GC.yv = yv;
    GC.axl = ea[0]; GC.bxl = eb[0]; GC.rxl = er[0];
    GC.axr = ea[1]; GC.bxr = eb[1]; GC.rxr = er[1];
    GC.ayl = ea[2]; GC.byl = eb[2]; GC.ryl = er[2];
    GC.ayr = ea[3]; GC.byr = eb[3]; GC.ryr = er[3];
    GC.dord = dord; GC.spec = spec;

    NdProblem P; memset(&P, 0, sizeof(P));
    P.d = d; P.spec = spec; P.tvar = tv;
    P.eval_monitor = o.eval_monitor;
    P.ysym = ysym;
    P.f = calloc(d, sizeof(Expr*));
    P.Y0 = calloc(d, sizeof(double));
    P.bind_y = malloc(sizeof(NdBind) * d);
    nd_bind_snapshot(&P.bind_t, tv);
    for (size_t i = 0; i < d; i++) nd_bind_snapshot(&P.bind_y[i], ysym[i]->data.symbol.name);

    Expr* lit_u  = nd_pde3_lit(fname, 0, 0, 0, tv, xv, yv);
    Expr* lit_xv = expr_new_symbol(xv);
    Expr* lit_yv = expr_new_symbol(yv);
    Expr** lit_sp = malloc(sizeof(Expr*) * (size_t)(nsp ? nsp : 1));
    for (int i = 0; i < nsp; i++) lit_sp[i] = nd_pde3_lit(fname, 0, sp[i][0], sp[i][1], tv, xv, yv);
    Expr** lit_ut = malloc(sizeof(Expr*) * (size_t)torder);
    for (int m = 1; m < torder; m++) lit_ut[m] = nd_pde3_lit(fname, m, 0, 0, tv, xv, yv);

    bool build_ok = (Gexpr != NULL);
    for (size_t iy = 1; iy + 1 < ny && build_ok; iy++)
        for (size_t ix = 1; ix + 1 < nx && build_ok; ix++) {
            size_t iu = (iy - 1) * (nx - 2) + (ix - 1);
            size_t base = iu * (size_t)torder;
            double xj = xmin + (double)ix * hx, yj = ymin + (double)iy * hy;
            for (int m = 0; m + 1 < torder; m++)
                P.f[base + (size_t)m] = expr_copy(ysym[base + (size_t)m + 1]);
            size_t cnt = 0;
            Expr* lits[2 * ND_MAX_SORDER + 5]; Expr* subs[2 * ND_MAX_SORDER + 5];
            lits[cnt] = lit_u; subs[cnt] = expr_copy(ysym[base + 0]); cnt++;
            for (int i = 0; i < nsp; i++) {
                int axis = (sp[i][0] > 0) ? 0 : 1;
                int deriv = (axis == 0) ? sp[i][0] : sp[i][1];
                lits[cnt] = lit_sp[i];
                subs[cnt] = nd_stencil2d(&GC, (int)ix, (int)iy, axis, deriv, dord);
                cnt++;
            }
            for (int m = 1; m < torder; m++) {
                lits[cnt] = lit_ut[m]; subs[cnt] = expr_copy(ysym[base + (size_t)m]); cnt++;
            }
            lits[cnt] = lit_xv; subs[cnt] = expr_new_real(xj); cnt++;
            lits[cnt] = lit_yv; subs[cnt] = expr_new_real(yj); cnt++;
            P.f[base + (size_t)torder - 1] = nd_replace_all(expr_copy(Gexpr), lits, subs, cnt);
            for (size_t s = 0; s < cnt; s++) expr_free(subs[s]);
        }
    expr_free(lit_u); expr_free(lit_xv); expr_free(lit_yv);
    for (int i = 0; i < nsp; i++) expr_free(lit_sp[i]);
    free(lit_sp);
    for (int m = 1; m < torder; m++) expr_free(lit_ut[m]);
    free(lit_ut);
    for (size_t i = 0; i < d; i++) if (!P.f[i]) build_ok = false;

    /* initial data */
    for (size_t iy = 1; iy + 1 < ny && build_ok; iy++)
        for (size_t ix = 1; ix + 1 < nx && build_ok; ix++) {
            size_t iu = (iy - 1) * (nx - 2) + (ix - 1);
            double xj = xmin + (double)ix * hx, yj = ymin + (double)iy * hy;
            for (int m = 0; m < torder; m++) {
                Expr* lit2[2] = { expr_new_symbol(xv), expr_new_symbol(yv) };
                Expr* sub2[2] = { expr_new_real(xj), expr_new_real(yj) };
                Expr* r = nd_replace_all(expr_copy(ic[m]), lit2, sub2, 2);
                expr_free(lit2[0]); expr_free(lit2[1]); expr_free(sub2[0]); expr_free(sub2[1]);
                double v; bool okv = r && nd_eval_to_double(r, spec, &v);
                expr_free(r);
                if (!okv) { build_ok = false; break; }
                P.Y0[iu * (size_t)torder + (size_t)m] = v;
            }
        }
    P.t0 = t0; P.tmin = tmin; P.tmax = tmax;

    if (build_ok) P.op = nd_operator_try_build(&P);
    bool parabolic = (torder == 1);   /* 2-D evolution PDEs are diffusion-like */

    Expr* result = NULL;
    if (build_ok) {
        const char* tim = o.method;
        if (tim && strcmp(tim, "MethodOfLines") == 0) tim = NULL;
        const NdStepper* S = tim ? nd_lookup_stepper(tim)
                           : (parabolic ? nd_lookup_stepper("BDF") : nd_default_stepper());
        if (!S) S = nd_default_stepper();
        NdSolution sol; nd_solution_init(&sol, d);
        NdStatus st = nd_integrate(&P, S, &o, &sol);
        if (st == ND_ERR_MAXSTEPS) nd_mol_warn("mxst", "maximum number of steps reached");
        else if (st == ND_ERR_STEPSIZE) nd_mol_warn("ndsz", "step size effectively zero");
        else if (st == ND_ERR_NONCONV) nd_mol_warn("ndcf", "corrector failed to converge");
        else if (st == ND_ERR_SAMPLE) nd_mol_warn("nrnum", "spatial operator did not evaluate");

        /* assemble the 3-D InterpolatingFunction over (t,x,y) */
        if (sol.n >= 2) {
            size_t npts = sol.n * nx * ny;
            Expr** entries = malloc(sizeof(Expr*) * npts);
            size_t p = 0;
            for (size_t i = 0; i < sol.n; i++) {
                double ti = sol.ts[i];
                nd_bind_set(&P.bind_t, expr_new_real(ti));
                /* bind interior states so Neumann/Robin boundary nodes, which
                 * are eliminated into interior values, evaluate correctly. */
                for (size_t k = 0; k < d; k++)
                    nd_bind_set(&P.bind_y[k], expr_new_real(sol.Ys[i * d + k]));
                for (size_t iy = 0; iy < ny; iy++)
                    for (size_t ix = 0; ix < nx; ix++) {
                        double xg = xmin + (double)ix * hx, yg = ymin + (double)iy * hy;
                        double val;
                        if (ix >= 1 && ix + 1 < nx && iy >= 1 && iy + 1 < ny) {
                            size_t iu = (iy - 1) * (nx - 2) + (ix - 1);
                            val = sol.Ys[i * d + iu * (size_t)torder + 0];
                        } else {
                            Expr* bv = nd_nodeval_2d(&GC, (int)ix, (int)iy);
                            if (!nd_bound_eval(bv, spec, &val)) val = 0.0;
                            expr_free(bv);
                        }
                        Expr* coord[3] = { expr_new_real(ti), expr_new_real(xg), expr_new_real(yg) };
                        Expr* coordL = expr_new_function(expr_new_symbol(SYM_List), coord, 3);
                        Expr* pair[2] = { coordL, expr_new_real(val) };
                        entries[p++] = expr_new_function(expr_new_symbol(SYM_List), pair, 2);
                    }
            }
            Expr* data = expr_new_function(expr_new_symbol(SYM_List), entries, npts);
            free(entries);
            Expr* ifun = eval_and_free(expr_new_function(expr_new_symbol(SYM_Interpolation), &data, 1));
            if (head_is(ifun, SYM_InterpolatingFunction)) {
                Expr* lhs;
                if (applied) {
                    Expr* ar[3] = { expr_new_symbol(tv), expr_new_symbol(xv), expr_new_symbol(yv) };
                    lhs = expr_new_function(expr_new_symbol(fname), ar, 3);
                    Expr* ar2[3] = { expr_new_symbol(tv), expr_new_symbol(xv), expr_new_symbol(yv) };
                    ifun = expr_new_function(ifun, ar2, 3);
                } else lhs = expr_new_symbol(fname);
                Expr* rule = nd_call2(SYM_Rule, lhs, ifun);
                Expr* inner = expr_new_function(expr_new_symbol(SYM_List), &rule, 1);
                result = expr_new_function(expr_new_symbol(SYM_List), &inner, 1);
            } else expr_free(ifun);
        }
        nd_solution_free(&sol);
    }
    nd_bind_restore(&P.bind_t);
    for (size_t i = 0; i < d; i++) nd_bind_restore(&P.bind_y[i]);

    nd_operator_free(P.op);
    for (size_t i = 0; i < d; i++) { expr_free(P.f[i]); expr_free(ysym[i]); }
    free(P.f); free(ysym); free(P.Y0); free(P.bind_y);
    if (P.jac) {
        for (size_t i = 0; i < d; i++) { for (size_t j = 0; j < d; j++) expr_free(P.jac[i][j]); free(P.jac[i]); }
        free(P.jac);
    }
    expr_free(Gexpr);
    for (int a = 0; a < torder; a++) expr_free(ic[a]);
    free(ic);
    for (int i = 0; i < 4; i++) { expr_free(ea[i]); expr_free(eb[i]); expr_free(er[i]); }
    for (size_t e = 0; e < neq; e++) expr_free(eqitems[e]);
    free(eqitems);
    return result;
}
