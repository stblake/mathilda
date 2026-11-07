/* Mathilda -- the finite element method for planar elliptic problems.
 *
 * See fem.h for the user-level surface. The pipeline:
 *
 *   region      Rectangle / Polygon / Triangle / Disk (sector) / Annulus
 *               -> boundary loops with straight and circular edges (FemGeom)
 *   mesh        Delaunay refinement to the requested MaxCellMeasure
 *               (fem_mesh.c), then quadratic nodes; midside nodes of curved
 *               boundary edges sit on the arc (isoparametric elements)
 *   operator    the residual R = op u - f of the equation is read as
 *                 -c11 u_xx - 2 c12 u_xy - c22 u_yy + bx u_x + by u_y + d u + r
 *               (coefficients from D[R, u_..]; c constant, b, d, r any
 *               functions of x and y), NeumannValue[g, pred] terms apart
 *   weak form   a(u, v) = Int c grad u . grad v + (b . grad u) v + d u v,
 *               plus s Int_edges g v for every NeumannValue term s NV[g, pred]
 *               on the boundary edges whose midpoint satisfies pred (g may be
 *               linear in u: Robin); every other boundary edge has the
 *               natural zero-flux condition
 *   quadrature  the degree-6 12-point Dunavant rule on each triangle, 5-point
 *               Gauss-Legendre on each boundary edge
 *   Dirichlet   DirichletCondition[u == h, pred] fixes every boundary node
 *               where pred holds (True: all of them); fixed values are moved
 *               to the right-hand side
 *   solve       K u = F by a sparse direct factorisation (CHOLMOD when K is
 *               symmetric positive definite, UMFPACK otherwise); the
 *               eigenproblem K x = lambda M x by shift-and-invert Krylov on
 *               (K - sigma M)^-1 M, every eigenvalue then refined as the
 *               Rayleigh quotient x.K.x / x.M.x
 *
 * Accuracy: P2 elements converge as h^3 in the solution and h^4 in the
 * eigenvalues on a smooth problem; a re-entrant corner of angle pi/alpha caps
 * both at h^(2 alpha). The default mesh (MaxCellMeasure -> area/1500, about
 * 3000 quadratic nodes) gives about 7 digits for a smooth eigenvalue.
 *
 * Coefficient and source fields are evaluated at the quadrature points in
 * one batch: compiled (autocompile) when they are in the compilable subset,
 * else by the evaluator. */
#include "fem.h"
#include "fem_mesh.h"
#include "fem_chol.h"
#include "ndsolve_common.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../symtab.h"
#include "../attr.h"
#include "../eval.h"
#include "../common.h"
#include "../numeric.h"
#include "../ndarray.h"
#include "../message.h"
#include "../compile/autocompile.h"
#include "../modular.h"           /* mth_block_enter / leave */
#include "../sparse/sparse_la.h"
#include "../linalg/krylov_schur.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const char* SYM_DirichletCondition;
static const char* SYM_NeumannValue;
static const char* SYM_FEMMesh;          /* NDSolve`FEM`Mesh[coords, elements, curved] */
static const char* SYM_FEMFunction;      /* NDSolve`FEM`MeshFunction[mesh, values]      */

/* ------------------------------------------------------------------ */
/*  Small helpers                                                      */
/* ------------------------------------------------------------------ */

static bool fem_num(Expr* e, double* out) {
    return nd_eval_to_double(e, numeric_machine_spec(), out) && isfinite(*out);
}

static bool fem_point(Expr* e, double* x, double* y) {
    if (!head_is(e, SYM_List)) {
        if (e && e->type == EXPR_NDARRAY && e->data.ndarray.rank == 1 && e->data.ndarray.dims[0] == 2) {
            Expr* a = ndarray_element_to_expr(e, 0), *b = ndarray_element_to_expr(e, 1);
            bool ok = fem_num(a, x) && fem_num(b, y);
            expr_free(a); expr_free(b);
            return ok;
        }
        return false;
    }
    if (e->data.function.arg_count != 2) return false;
    return fem_num(e->data.function.args[0], x) && fem_num(e->data.function.args[1], y);
}

/* The elements of a List or rank-1 packed array (fresh copies), or NULL. */
static Expr** fem_items(Expr* e, size_t* n) {
    if (head_is(e, SYM_List)) {
        *n = e->data.function.arg_count;
        Expr** v = malloc(sizeof(Expr*) * (*n ? *n : 1));
        for (size_t i = 0; i < *n; i++) v[i] = expr_copy(e->data.function.args[i]);
        return v;
    }
    if (e && e->type == EXPR_NDARRAY) {
        Expr* l = ndarray_to_nested_list(e);
        Expr** v = fem_items(l, n);
        expr_free(l);
        return v;
    }
    return NULL;
}
static void fem_items_free(Expr** v, size_t n) {
    if (!v) return;
    for (size_t i = 0; i < n; i++) expr_free(v[i]);
    free(v);
}

static Expr* fem_list(Expr** items, size_t n) {
    return expr_new_function(expr_new_symbol(SYM_List), items, n);
}

/* A packed real vector / matrix (the buffer is adopted). */
static Expr* fem_packed(double* buf, int rank, int64_t d0, int64_t d1) {
    int64_t dims[2] = { d0, d1 };
    Expr* e = expr_new_ndarray_raw(rank, dims, buf, NDT_FLOAT64);
    e->data.ndarray.present_as = NDA_HEAD_LIST;
    return e;
}

/* Read a real buffer of exactly n elements from a packed array or List. */
static bool fem_read_reals(const Expr* e, int64_t n, double* out) {
    if (e && e->type == EXPR_NDARRAY) {
        int64_t tot = 1;
        for (int k = 0; k < e->data.ndarray.rank; k++) tot *= e->data.ndarray.dims[k];
        if (tot != n || ndt_is_complex(e->data.ndarray.dtype)) return false;
        for (int64_t i = 0; i < n; i++) {
            double re, im;
            ndt_get(e->data.ndarray.data, (size_t)i, e->data.ndarray.dtype, &re, &im);
            out[i] = re;
        }
        return true;
    }
    if (!head_is((Expr*)e, SYM_List) || (int64_t)e->data.function.arg_count != n) return false;
    for (int64_t i = 0; i < n; i++)
        if (!nd_to_double(e->data.function.args[i], &out[i]) || !isfinite(out[i])) return false;
    return true;
}

static bool fem_contains_symbol(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == name;
    if (e->type != EXPR_FUNCTION) return false;
    if (fem_contains_symbol(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (fem_contains_symbol(e->data.function.args[i], name)) return true;
    return false;
}

/* Evaluate the field f(x, y) (borrowed) at n points. Constant fields are
 * broadcast; then the compiled program; then one evaluation of f on packed
 * coordinate vectors; then point by point. */
static bool fem_field(Expr* f, const char* xv, const char* yv,
                      const double* X, const double* Y, int64_t n, double* out) {
    if (!f) { for (int64_t i = 0; i < n; i++) out[i] = 0.0; return true; }
    bool hx = fem_contains_symbol(f, xv), hy = fem_contains_symbol(f, yv);
    if (!hx && !hy) {
        double c;
        if (!fem_num(f, &c)) return false;
        for (int64_t i = 0; i < n; i++) out[i] = c;
        return true;
    }
    Expr* vx = expr_new_symbol(xv), *vy = expr_new_symbol(yv);
    const Expr* vs[2] = { vx, vy };
    AutoCompiled* ac = autocompile_new(f, vs, 2);
    bool ok = true;
    if (ac) {
        for (int64_t i = 0; i < n && ok; i++) {
            double in[2] = { X[i], Y[i] };
            if (!autocompiled_eval_real(ac, in, &out[i]) || !isfinite(out[i])) {
                Expr* subs[2] = { expr_new_real(X[i]), expr_new_real(Y[i]) };
                Expr* lits[2] = { vx, vy };
                Expr* r = nd_replace_all(expr_copy(f), lits, subs, 2);
                ok = fem_num(r, &out[i]);
                expr_free(r); expr_free(subs[0]); expr_free(subs[1]);
            }
        }
        autocompiled_free(ac);
        expr_free(vx); expr_free(vy);
        return ok;
    }
    /* one batch evaluation over packed coordinates */
    double* bx = malloc(sizeof(double) * (size_t)n);
    double* by = malloc(sizeof(double) * (size_t)n);
    memcpy(bx, X, sizeof(double) * (size_t)n);
    memcpy(by, Y, sizeof(double) * (size_t)n);
    Expr* subs[2] = { fem_packed(bx, 1, n, 0), fem_packed(by, 1, n, 0) };
    Expr* lits[2] = { vx, vy };
    Expr* r = nd_replace_all(expr_copy(f), lits, subs, 2);
    Expr* rn = r ? eval_and_free(nd_call1(SYM_N, r)) : NULL;
    expr_free(subs[0]); expr_free(subs[1]);
    bool batch = rn && fem_read_reals(rn, n, out);
    expr_free(rn);
    if (!batch) {
        for (int64_t i = 0; i < n && ok; i++) {
            Expr* s2[2] = { expr_new_real(X[i]), expr_new_real(Y[i]) };
            Expr* ri = nd_replace_all(expr_copy(f), lits, s2, 2);
            ok = fem_num(ri, &out[i]);
            expr_free(ri); expr_free(s2[0]); expr_free(s2[1]);
        }
    }
    expr_free(vx); expr_free(vy);
    return ok;
}

/* A predicate pred(x, y) at one point: True, False, or -1 (not Boolean). */
static int fem_pred(Expr* pred, const char* xv, const char* yv, double x, double y) {
    if (pred->type == EXPR_SYMBOL) {
        if (pred->data.symbol.name == SYM_True) return 1;
        if (pred->data.symbol.name == SYM_False) return 0;
    }
    Expr* lits[2] = { expr_new_symbol(xv), expr_new_symbol(yv) };
    Expr* subs[2] = { expr_new_real(x), expr_new_real(y) };
    Expr* r = nd_replace_all(expr_copy(pred), lits, subs, 2);
    int v = -1;
    if (r && r->type == EXPR_SYMBOL) {
        if (r->data.symbol.name == SYM_True) v = 1;
        else if (r->data.symbol.name == SYM_False) v = 0;
    }
    expr_free(r);
    for (int k = 0; k < 2; k++) { expr_free(lits[k]); expr_free(subs[k]); }
    return v;
}

/* ------------------------------------------------------------------ */
/*  Regions                                                            */
/* ------------------------------------------------------------------ */

bool fem_is_region_spec(const Expr* e) {
    return head_is((Expr*)e, SYM_Element) && e->data.function.arg_count == 2
        && head_is(e->data.function.args[0], SYM_List);
}

/* Append the points of an arc (centre, r) from angle t0 to t1, excluding the
 * end point, nseg edges each on arc `ai`. */
static void arc_points(double* xy, int* arcs, int* n, double cx, double cy, double r,
                       double t0, double t1, int nseg, int ai) {
    for (int k = 0; k < nseg; k++) {
        double t = t0 + (t1 - t0) * (double)k / (double)nseg;
        xy[2 * *n] = cx + r * cos(t); xy[2 * *n + 1] = cy + r * sin(t);
        arcs[*n] = ai;
        (*n)++;
    }
}

static int arc_segments(double dtheta) {
    int n = (int)ceil(fabs(dtheta) / (M_PI / 8.0) - 1e-9);
    return n < 2 ? 2 : n;
}

/* The boundary of `reg` (evaluated) as a FemGeom; false with *why set. */
static bool fem_region(Expr* reg0, FemGeom* g, const char** why) {
    fem_geom_init(g);
    Expr* reg = eval_and_free(expr_copy(reg0));
    bool ok = false;
    *why = "the region is not a Rectangle, Polygon, Triangle, Disk or Annulus";
    if (!reg || reg->type != EXPR_FUNCTION || reg->data.function.head->type != EXPR_SYMBOL) goto done;
    const char* h = reg->data.function.head->data.symbol.name;
    size_t na = reg->data.function.arg_count;
    Expr** A = reg->data.function.args;
    if (h == SYM_Rectangle) {
        double x0 = 0, y0 = 0, x1 = 1, y1 = 1;
        if (na >= 1 && !fem_point(A[0], &x0, &y0)) goto done;
        if (na >= 2) { if (!fem_point(A[1], &x1, &y1)) goto done; }
        else if (na == 1) { x1 = x0 + 1; y1 = y0 + 1; }
        if (na > 2 || !(x1 > x0) || !(y1 > y0)) goto done;
        double xy[8] = { x0, y0, x1, y0, x1, y1, x0, y1 };
        ok = fem_geom_add_loop(g, xy, NULL, 4, false);
    } else if (h == SYM_Polygon || h == intern_symbol("Triangle")) {
        Expr** P = NULL; size_t np = 0;
        if (na == 0 && h != SYM_Polygon) {
            double xy[6] = { 0, 0, 1, 0, 0, 1 };
            ok = fem_geom_add_loop(g, xy, NULL, 3, false);
            goto done;
        }
        if (na != 1 || !(P = fem_items(A[0], &np)) || np < 3) { fem_items_free(P, np); goto done; }
        double* xy = malloc(sizeof(double) * 2 * np);
        bool pts = true;
        for (size_t i = 0; i < np && pts; i++) pts = fem_point(P[i], &xy[2 * i], &xy[2 * i + 1]);
        fem_items_free(P, np);
        /* a closing repeat of the first vertex */
        if (pts && np > 3 && xy[0] == xy[2 * (np - 1)] && xy[1] == xy[2 * (np - 1) + 1]) np--;
        if (pts) ok = fem_geom_add_loop(g, xy, NULL, (int)np, false);
        else *why = "the polygon has a non-numeric vertex";
        free(xy);
    } else if (h == SYM_Disk) {
        double cx = 0, cy = 0, r = 1, t0 = 0, t1 = 2 * M_PI;
        bool sector = false;
        if (na >= 1 && !fem_point(A[0], &cx, &cy)) goto done;
        if (na >= 2) {
            double rx, ry;
            if (fem_point(A[1], &rx, &ry)) {
                if (fabs(rx - ry) > 1e-14 * fabs(rx)) { *why = "an elliptical Disk is not supported"; goto done; }
                r = rx;
            } else if (!fem_num(A[1], &r)) goto done;
        }
        if (na >= 3) {
            if (!fem_point(A[2], &t0, &t1)) goto done;
            sector = fabs(t1 - t0) < 2 * M_PI - 1e-12;
        }
        if (na > 3 || !(r > 0)) goto done;
        int ai = fem_geom_add_arc(g, cx, cy, r);
        if (!sector) {
            double xy[2 * 64]; int arcs[64], n = 0;
            arc_points(xy, arcs, &n, cx, cy, r, 0.0, 2 * M_PI, 16, ai);
            ok = fem_geom_add_loop(g, xy, arcs, n, false);
        } else {
            if (t1 < t0) { double s = t0; t0 = t1; t1 = s; }
            int ns = arc_segments(t1 - t0);
            double* xy = malloc(sizeof(double) * 2 * (size_t)(ns + 2));
            int* arcs = malloc(sizeof(int) * (size_t)(ns + 2)), n = 0;
            arc_points(xy, arcs, &n, cx, cy, r, t0, t1, ns, ai);
            xy[2 * n] = cx + r * cos(t1); xy[2 * n + 1] = cy + r * sin(t1); arcs[n] = -1; n++;
            xy[2 * n] = cx; xy[2 * n + 1] = cy; arcs[n] = -1; n++;
            ok = fem_geom_add_loop(g, xy, arcs, n, false);
            free(xy); free(arcs);
        }
    } else if (h == SYM_Annulus) {
        double cx = 0, cy = 0, r0 = 0.5, r1 = 1, t0 = 0, t1 = 2 * M_PI;
        if (na >= 1 && !fem_point(A[0], &cx, &cy)) goto done;
        if (na >= 2 && !fem_point(A[1], &r0, &r1)) goto done;
        bool sector = false;
        if (na >= 3) {
            if (!fem_point(A[2], &t0, &t1)) goto done;
            sector = fabs(t1 - t0) < 2 * M_PI - 1e-12;
        }
        if (na > 3 || !(r0 > 0) || !(r1 > r0)) goto done;
        int a0 = fem_geom_add_arc(g, cx, cy, r0), a1 = fem_geom_add_arc(g, cx, cy, r1);
        if (!sector) {
            double xy[2 * 64]; int arcs[64], n = 0;
            arc_points(xy, arcs, &n, cx, cy, r1, 0.0, 2 * M_PI, 16, a1);
            ok = fem_geom_add_loop(g, xy, arcs, n, false);
            n = 0;
            arc_points(xy, arcs, &n, cx, cy, r0, 0.0, 2 * M_PI, 16, a0);
            ok = ok && fem_geom_add_loop(g, xy, arcs, n, true);
        } else {
            if (t1 < t0) { double s = t0; t0 = t1; t1 = s; }
            int ns = arc_segments(t1 - t0);
            double* xy = malloc(sizeof(double) * 4 * (size_t)(ns + 2));
            int* arcs = malloc(sizeof(int) * 2 * (size_t)(ns + 2)), n = 0;
            arc_points(xy, arcs, &n, cx, cy, r1, t0, t1, ns, a1);
            xy[2 * n] = cx + r1 * cos(t1); xy[2 * n + 1] = cy + r1 * sin(t1); arcs[n] = -1; n++;
            /* inner arc backwards from t1 to t0 */
            arc_points(xy, arcs, &n, cx, cy, r0, t1, t0, ns, a0);
            xy[2 * n] = cx + r0 * cos(t0); xy[2 * n + 1] = cy + r0 * sin(t0); arcs[n] = -1; n++;
            ok = fem_geom_add_loop(g, xy, arcs, n, false);
            free(xy); free(arcs);
        }
    }
done:
    expr_free(reg);
    if (!ok) fem_geom_free(g);
    return ok;
}

/* ------------------------------------------------------------------ */
/*  Options                                                            */
/* ------------------------------------------------------------------ */

/* Search an option value recursively for key -> number (key a string or a
 * symbol of that name). */
static bool fem_find_opt(const Expr* e, const char* key, double* out) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if ((head_is((Expr*)e, SYM_Rule) || head_is((Expr*)e, SYM_RuleDelayed)) && e->data.function.arg_count == 2) {
        const Expr* l = e->data.function.args[0];
        bool hit = (l->type == EXPR_STRING && strcmp(l->data.string, key) == 0)
                || (l->type == EXPR_SYMBOL && strcmp(l->data.symbol.name, key) == 0);
        if (hit) return fem_num(e->data.function.args[1], out);
        return fem_find_opt(e->data.function.args[1], key, out);
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (fem_find_opt(e->data.function.args[i], key, out)) return true;
    return false;
}

/* ------------------------------------------------------------------ */
/*  The operator                                                       */
/* ------------------------------------------------------------------ */

#define FEM_MAXCOND 16

typedef struct {
    const char* u; const char* xv; const char* yv;
    double c11, c12, c22;
    Expr* bx; Expr* by; Expr* d; Expr* r;     /* fields (owned) or NULL */
    int nnv;
    double nv_s[FEM_MAXCOND];
    Expr* nv_g0[FEM_MAXCOND]; Expr* nv_g1[FEM_MAXCOND]; Expr* nv_pred[FEM_MAXCOND];
    int ndc;
    Expr* dc_val[FEM_MAXCOND]; Expr* dc_pred[FEM_MAXCOND];
} FemPde;

static void fem_pde_free(FemPde* P) {
    expr_free(P->bx); expr_free(P->by); expr_free(P->d); expr_free(P->r);
    for (int k = 0; k < P->nnv; k++) { expr_free(P->nv_g0[k]); expr_free(P->nv_g1[k]); expr_free(P->nv_pred[k]); }
    for (int k = 0; k < P->ndc; k++) { expr_free(P->dc_val[k]); expr_free(P->dc_pred[k]); }
}

/* u[x, y] or Derivative[i, j][u][x, y] literal */
static Expr* fem_ulit(const FemPde* P, int i, int j) {
    Expr* ar[2] = { expr_new_symbol(P->xv), expr_new_symbol(P->yv) };
    if (i == 0 && j == 0) return expr_new_function(expr_new_symbol(P->u), ar, 2);
    Expr* di[2] = { expr_new_integer(i), expr_new_integer(j) };
    Expr* dh = expr_new_function(expr_new_symbol(SYM_Derivative), di, 2);
    Expr* uu[1] = { expr_new_symbol(P->u) };
    return expr_new_function(expr_new_function(dh, uu, 1), ar, 2);
}

/* Collect the distinct NeumannValue[...] subexpressions of e. */
static void fem_collect_nv(const Expr* e, Expr** out, int* n) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (head_is((Expr*)e, SYM_NeumannValue)) {
        for (int k = 0; k < *n; k++) if (expr_eq(out[k], e)) return;
        if (*n < FEM_MAXCOND) out[(*n)++] = (Expr*)e;
        return;
    }
    fem_collect_nv(e->data.function.head, out, n);
    for (size_t i = 0; i < e->data.function.arg_count; i++) fem_collect_nv(e->data.function.args[i], out, n);
}

static const int FEM_DI[6] = { 0, 1, 0, 2, 1, 0 };
static const int FEM_DJ[6] = { 0, 0, 1, 0, 1, 2 };

/* Read the residual R (borrowed; "R == 0") into P. False with *why. When
 * homogeneous, any term free of u is an error (an eigen operator). */
static bool fem_read_operator(FemPde* P, Expr* R, bool homogeneous, const char** why) {
    Expr* nvs[FEM_MAXCOND]; int nnv = 0;
    fem_collect_nv(R, nvs, &nnv);
    if (homogeneous && nnv > 0) { *why = "NeumannValue is not allowed in an eigenvalue operator"; return false; }
    Expr* lits[6 + FEM_MAXCOND]; Expr* syms[6 + FEM_MAXCOND];
    char nm[48];
    for (int k = 0; k < 6; k++) {
        lits[k] = fem_ulit(P, FEM_DI[k], FEM_DJ[k]);
        snprintf(nm, sizeof nm, "NDSolve`FEM`u%d%d", FEM_DI[k], FEM_DJ[k]);
        syms[k] = expr_new_symbol(intern_symbol(nm));
    }
    for (int k = 0; k < nnv; k++) {
        lits[6 + k] = expr_copy(nvs[k]);
        snprintf(nm, sizeof nm, "NDSolve`FEM`nv%d", k);
        syms[6 + k] = expr_new_symbol(intern_symbol(nm));
    }
    int nl = 6 + nnv;
    Expr* Rs = nd_replace_all(expr_copy(R), lits, syms, (size_t)nl);
    Expr* co[6 + FEM_MAXCOND];
    bool ok = Rs != NULL;
    for (int k = 0; k < nl; k++) co[k] = NULL;
    for (int k = 0; k < nl && ok; k++) {
        Expr* da[2] = { expr_copy(Rs), expr_copy(syms[k]) };
        co[k] = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), da, 2));
        for (int j = 0; j < nl && ok; j++)
            if (fem_contains_symbol(co[k], syms[j]->data.symbol.name)) { ok = false; *why = "the equation is not linear in the unknown"; }
    }
    Expr* zeros[6 + FEM_MAXCOND];
    for (int k = 0; k < nl; k++) zeros[k] = expr_new_integer(0);
    Expr* R0 = ok ? nd_replace_all(expr_copy(Rs), syms, zeros, (size_t)nl) : NULL;
    for (int k = 0; k < nl; k++) expr_free(zeros[k]);
    if (ok && (fem_contains_symbol(R0, P->u))) { ok = false; *why = "the unknown appears in an unsupported form"; }
    for (int k = 0; k < nl && ok; k++)
        if (fem_contains_symbol(co[k], P->u)) { ok = false; *why = "the unknown appears in an unsupported form"; }
    if (ok) {
        /* second-order coefficients: constants */
        double cxx, cxy, cyy;
        if (fem_contains_symbol(co[3], P->xv) || fem_contains_symbol(co[3], P->yv)
            || fem_contains_symbol(co[4], P->xv) || fem_contains_symbol(co[4], P->yv)
            || fem_contains_symbol(co[5], P->xv) || fem_contains_symbol(co[5], P->yv)
            || !fem_num(co[3], &cxx) || !fem_num(co[4], &cxy) || !fem_num(co[5], &cyy)) {
            ok = false; *why = "the second-order coefficients must be constants";
        } else {
            P->c11 = -cxx; P->c22 = -cyy; P->c12 = -0.5 * cxy;
            bool posdef = P->c11 > 0 && P->c22 > 0 && P->c11 * P->c22 - P->c12 * P->c12 > 0;
            bool negdef = P->c11 < 0 && P->c22 < 0 && P->c11 * P->c22 - P->c12 * P->c12 > 0;
            if (negdef && !homogeneous) {
                /* Laplacian[u] == f: the same equation as -Laplacian[u] == -f
                 * (a NeumannValue keeps its meaning, n . (c grad u) for the
                 * c of the equation as written, through the sign of s) */
                P->c11 = -P->c11; P->c22 = -P->c22; P->c12 = -P->c12;
                for (int k = 0; k < nl; k++)
                    if (co[k]) co[k] = eval_and_free(nd_call2(SYM_Times, expr_new_integer(-1), co[k]));
                R0 = eval_and_free(nd_call2(SYM_Times, expr_new_integer(-1), R0));
                posdef = true;
            }
            if (!posdef) {
                ok = false; *why = homogeneous
                    ? "the operator is not elliptic (-div(c grad u) with c positive definite; an eigenvalue problem is posed for -Laplacian, not Laplacian)"
                    : "the operator is not elliptic (-div(c grad u) with c definite)";
            }
        }
    }
    if (ok) {
        P->d = co[0]; co[0] = NULL;
        P->bx = co[1]; co[1] = NULL;
        P->by = co[2]; co[2] = NULL;
        double z;
        if (fem_num(R0, &z) && z == 0.0) { expr_free(R0); R0 = NULL; }
        if (homogeneous && R0) { ok = false; *why = "the eigenvalue operator has a term free of the unknown"; }
        P->r = R0; R0 = NULL;
        for (int k = 0; k < nnv && ok; k++) {
            double s;
            if (!fem_num(co[6 + k], &s)) { ok = false; *why = "a NeumannValue term has a non-constant factor"; break; }
            Expr* nv = nvs[k];
            if (nv->data.function.arg_count != 2) { ok = false; *why = "NeumannValue takes a value and a predicate"; break; }
            /* g = g0 + g1 u */
            Expr* Us = expr_new_symbol(intern_symbol("NDSolve`FEM`ug"));
            Expr* ul = fem_ulit(P, 0, 0);
            Expr* gs = nd_replace_all(expr_copy(nv->data.function.args[0]), &ul, &Us, 1);
            Expr* da[2] = { expr_copy(gs), expr_copy(Us) };
            Expr* g1 = eval_and_free(expr_new_function(expr_new_symbol(SYM_D), da, 2));
            Expr* zz = expr_new_integer(0);
            Expr* g0 = nd_replace_all(expr_copy(gs), &Us, &zz, 1);
            expr_free(zz); expr_free(gs); expr_free(ul);
            if (fem_contains_symbol(g1, Us->data.symbol.name) || fem_contains_symbol(g1, P->u)
                || fem_contains_symbol(g0, P->u)) {
                ok = false; *why = "a NeumannValue must be linear in the unknown";
                expr_free(g0); expr_free(g1);
            } else {
                P->nv_s[P->nnv] = s;
                P->nv_g0[P->nnv] = g0; P->nv_g1[P->nnv] = g1;
                P->nv_pred[P->nnv] = expr_copy(nv->data.function.args[1]);
                P->nnv++;
            }
            expr_free(Us);
        }
    }
    expr_free(R0); expr_free(Rs);
    for (int k = 0; k < nl; k++) { expr_free(co[k]); expr_free(lits[k]); expr_free(syms[k]); }
    return ok;
}

/* DirichletCondition[u[x, y] == h, pred] */
static bool fem_read_dirichlet(FemPde* P, Expr* dc, const char** why) {
    if (dc->data.function.arg_count != 2 || !head_is(dc->data.function.args[0], SYM_Equal)
        || dc->data.function.args[0]->data.function.arg_count != 2) {
        *why = "DirichletCondition must be DirichletCondition[u == value, predicate]";
        return false;
    }
    if (P->ndc >= FEM_MAXCOND) { *why = "too many DirichletConditions"; return false; }
    Expr* L = dc->data.function.args[0]->data.function.args[0];
    Expr* Rr = dc->data.function.args[0]->data.function.args[1];
    Expr* ul = fem_ulit(P, 0, 0);
    Expr* val = NULL;
    if (expr_eq(L, ul) && !fem_contains_symbol(Rr, P->u)) val = Rr;
    else if (expr_eq(Rr, ul) && !fem_contains_symbol(L, P->u)) val = L;
    expr_free(ul);
    if (!val) { *why = "DirichletCondition must prescribe u[x, y] itself"; return false; }
    P->dc_val[P->ndc] = expr_copy(val);
    P->dc_pred[P->ndc] = expr_copy(dc->data.function.args[1]);
    P->ndc++;
    return true;
}

/* ------------------------------------------------------------------ */
/*  P2 element                                                         */
/* ------------------------------------------------------------------ */

/* Dunavant degree-6, 12 points: (l1, l2) barycentric with weights (sum 1). */
static const double DQ[12][3] = {
    { 0.249286745170910, 0.249286745170910, 0.116786275726379 },
    { 0.249286745170910, 0.501426509658179, 0.116786275726379 },
    { 0.501426509658179, 0.249286745170910, 0.116786275726379 },
    { 0.063089014491502, 0.063089014491502, 0.050844906370207 },
    { 0.063089014491502, 0.873821971016996, 0.050844906370207 },
    { 0.873821971016996, 0.063089014491502, 0.050844906370207 },
    { 0.310352451033784, 0.636502499121399, 0.082851075618374 },
    { 0.636502499121399, 0.053145049844817, 0.082851075618374 },
    { 0.053145049844817, 0.310352451033784, 0.082851075618374 },
    { 0.636502499121399, 0.310352451033784, 0.082851075618374 },
    { 0.310352451033784, 0.053145049844817, 0.082851075618374 },
    { 0.053145049844817, 0.636502499121399, 0.082851075618374 },
};
#define NQ 12

/* Shape functions on the reference triangle at (xi, eta): vertices 0, 1, 2
 * at (0,0), (1,0), (0,1); midsides 3: (0,1), 4: (1,2), 5: (2,0). */
static void p2_shape(double xi, double eta, double* N, double* Nx, double* Ne) {
    double l0 = 1.0 - xi - eta, l1 = xi, l2 = eta;
    N[0] = l0 * (2 * l0 - 1); N[1] = l1 * (2 * l1 - 1); N[2] = l2 * (2 * l2 - 1);
    N[3] = 4 * l0 * l1; N[4] = 4 * l1 * l2; N[5] = 4 * l2 * l0;
    if (!Nx) return;
    /* dl0 = (-1, -1), dl1 = (1, 0), dl2 = (0, 1) */
    Nx[0] = -(4 * l0 - 1); Ne[0] = -(4 * l0 - 1);
    Nx[1] = 4 * l1 - 1;    Ne[1] = 0;
    Nx[2] = 0;             Ne[2] = 4 * l2 - 1;
    Nx[3] = 4 * (l0 - l1); Ne[3] = -4 * l1;
    Nx[4] = 4 * l2;        Ne[4] = 4 * l1;
    Nx[5] = -4 * l2;       Ne[5] = 4 * (l0 - l2);
}

/* Jacobian of the isoparametric map at (xi, eta); det and the physical
 * gradients of the six shape functions. */
static double p2_jac(const double* X, const double* Y, const double* Nx, const double* Ne,
                     double* gx, double* gy) {
    double j11 = 0, j12 = 0, j21 = 0, j22 = 0;
    for (int k = 0; k < 6; k++) {
        j11 += Nx[k] * X[k]; j12 += Ne[k] * X[k];
        j21 += Nx[k] * Y[k]; j22 += Ne[k] * Y[k];
    }
    double det = j11 * j22 - j12 * j21;
    if (gx) for (int k = 0; k < 6; k++) {
        gx[k] = ( j22 * Nx[k] - j21 * Ne[k]) / det;
        gy[k] = (-j12 * Nx[k] + j11 * Ne[k]) / det;
    }
    return det;
}

/* ------------------------------------------------------------------ */
/*  Sparse assembly                                                    */
/* ------------------------------------------------------------------ */

typedef struct { int64_t n, nnz, cap; int64_t* r; int64_t* c; double* v; } Coo;

static void coo_add(Coo* A, int64_t i, int64_t j, double v) {
    if (A->nnz == A->cap) {
        A->cap = A->cap ? 2 * A->cap : 4096;
        A->r = realloc(A->r, sizeof(int64_t) * (size_t)A->cap);
        A->c = realloc(A->c, sizeof(int64_t) * (size_t)A->cap);
        A->v = realloc(A->v, sizeof(double) * (size_t)A->cap);
    }
    A->r[A->nnz] = i; A->c[A->nnz] = j; A->v[A->nnz] = v; A->nnz++;
}
static void coo_free(Coo* A) { free(A->r); free(A->c); free(A->v); memset(A, 0, sizeof *A); }

/* CSR (owned arrays) with duplicates summed and columns sorted. */
typedef struct { int64_t n, nnz; int64_t* rp; int64_t* ci; double* x; } Csr;

static void csr_from_coo(const Coo* A, double scale_b, const Coo* B, Csr* out) {
    int64_t n = A->n, tot = A->nnz + (B ? B->nnz : 0);
    int64_t* cnt = calloc((size_t)n + 1, sizeof(int64_t));
    for (int64_t k = 0; k < A->nnz; k++) cnt[A->r[k] + 1]++;
    if (B) for (int64_t k = 0; k < B->nnz; k++) cnt[B->r[k] + 1]++;
    for (int64_t i = 0; i < n; i++) cnt[i + 1] += cnt[i];
    int64_t* ci = malloc(sizeof(int64_t) * (size_t)(tot ? tot : 1));
    double* x = malloc(sizeof(double) * (size_t)(tot ? tot : 1));
    int64_t* pos = malloc(sizeof(int64_t) * (size_t)(n + 1));
    memcpy(pos, cnt, sizeof(int64_t) * (size_t)(n + 1));
    for (int64_t k = 0; k < A->nnz; k++) { int64_t p = pos[A->r[k]]++; ci[p] = A->c[k]; x[p] = A->v[k]; }
    if (B) for (int64_t k = 0; k < B->nnz; k++) { int64_t p = pos[B->r[k]]++; ci[p] = B->c[k]; x[p] = scale_b * B->v[k]; }
    /* sort each row by column (insertion sort: rows are short) and merge */
    int64_t* rp = malloc(sizeof(int64_t) * (size_t)(n + 1));
    int64_t w = 0;
    rp[0] = 0;
    for (int64_t i = 0; i < n; i++) {
        int64_t s = cnt[i], e = cnt[i + 1];
        for (int64_t a = s + 1; a < e; a++) {
            int64_t c = ci[a]; double v = x[a]; int64_t b = a - 1;
            while (b >= s && ci[b] > c) { ci[b + 1] = ci[b]; x[b + 1] = x[b]; b--; }
            ci[b + 1] = c; x[b + 1] = v;
        }
        for (int64_t a = s; a < e; a++) {
            if (w > rp[i] && ci[w - 1] == ci[a]) x[w - 1] += x[a];
            else { ci[w] = ci[a]; x[w] = x[a]; w++; }
        }
        rp[i + 1] = w;
    }
    free(cnt); free(pos);
    out->n = n; out->nnz = w; out->rp = rp; out->ci = ci; out->x = x;
}
static void csr_free(Csr* A) { free(A->rp); free(A->ci); free(A->x); memset(A, 0, sizeof *A); }
static void csr_view(const Csr* A, SpMat* M) {
    memset(M, 0, sizeof *M);
    M->m = A->n; M->n = A->n; M->nnz = A->nnz; M->rp = A->rp; M->ci = A->ci; M->x = A->x; M->cplx = false;
}
static void csr_mv(const Csr* A, const double* x, double* y) {
    for (int64_t i = 0; i < A->n; i++) {
        double s = 0.0;
        for (int64_t p = A->rp[i]; p < A->rp[i + 1]; p++) s += A->x[p] * x[A->ci[p]];
        y[i] = s;
    }
}

/* A factorisation of a free-node matrix: the nested-dissection Cholesky of
 * fem_chol.c for a symmetric matrix (always, unless it is not positive
 * definite), SpFactor (UMFPACK, or the band / dense LU) otherwise. */
typedef struct { FemChol* C; SpFactor* F; } FemSolver;

static bool fem_solver_factor(FemSolver* S, const Csr* A, bool symmetric, const double* xyf) {
    memset(S, 0, sizeof *S);
    if (symmetric) {
        const char* why;
        S->C = fem_chol_factor(A->n, A->rp, A->ci, A->x, xyf, &why);
        if (S->C) return true;
    }
    SpMat Am; csr_view(A, &Am);
    SpfStatus st;
    S->F = spf_factor(&Am, SPF_WANT_AUTO, &st);
    return S->F != NULL;
}
static bool fem_solver_solve(const FemSolver* S, const double* b, double* x) {
    if (S->C) { fem_chol_solve(S->C, b, x); return true; }
    return S->F && spf_solve(S->F, b, x, 1, false);
}
static void fem_solver_free(FemSolver* S) {
    fem_chol_free(S->C);
    if (S->F) spf_free(S->F);
    memset(S, 0, sizeof *S);
}

/* The discrete problem on a mesh. */
typedef struct {
    FemMesh mesh;
    int64_t nfree;
    int64_t* map;        /* node -> free index, or -1 (Dirichlet) */
    double* ufix;        /* Dirichlet values at fixed nodes */
    Coo K, M;            /* over the free nodes */
    double* F;           /* load over the free nodes */
    bool symmetric;
} FemSystem;

static void fem_system_free(FemSystem* S) {
    fem_mesh_free(&S->mesh);
    free(S->map); free(S->ufix); free(S->F);
    coo_free(&S->K); coo_free(&S->M);
}

/* Coordinates of the free nodes, for the nested-dissection ordering. */
static double* fem_free_xy(const FemSystem* S) {
    double* xy = malloc(sizeof(double) * 2 * (size_t)(S->nfree ? S->nfree : 1));
    for (int64_t i = 0; i < S->mesh.nn; i++)
        if (S->map[i] >= 0) {
            xy[2 * S->map[i]] = S->mesh.x[2 * i];
            xy[2 * S->map[i] + 1] = S->mesh.x[2 * i + 1];
        }
    return xy;
}

/* Assemble K (and M when want_mass, F when want_load) for P on S->mesh. */
static bool fem_assemble(const FemPde* P, FemSystem* S, bool want_mass, bool want_load, const char** why) {
    FemMesh* m = &S->mesh;
    int64_t ne = m->ne, nn = m->nn;
    /* Dirichlet nodes */
    S->map = malloc(sizeof(int64_t) * (size_t)nn);
    S->ufix = calloc((size_t)nn, sizeof(double));
    int64_t nf = 0;
    for (int64_t i = 0; i < nn; i++) {
        S->map[i] = -2;
        if (!m->bnode[i]) continue;
        double x = m->x[2 * i], y = m->x[2 * i + 1];
        for (int k = 0; k < P->ndc; k++) {
            int pv = fem_pred(P->dc_pred[k], P->xv, P->yv, x, y);
            if (pv < 0) { *why = "a DirichletCondition predicate did not evaluate to True or False"; return false; }
            if (pv) {
                double v;
                if (!fem_field(P->dc_val[k], P->xv, P->yv, &x, &y, 1, &v)) {
                    *why = "a Dirichlet value is not numeric"; return false;
                }
                S->map[i] = -1; S->ufix[i] = v;
                break;
            }
        }
    }
    for (int64_t i = 0; i < nn; i++) if (S->map[i] == -2) S->map[i] = nf++;
    S->nfree = nf;
    S->K.n = nf; S->M.n = nf;
    if (want_load) S->F = calloc((size_t)(nf ? nf : 1), sizeof(double));

    /* quadrature points of every element, for batch field evaluation */
    int64_t nqp = ne * NQ;
    double* QX = malloc(sizeof(double) * (size_t)nqp);
    double* QY = malloc(sizeof(double) * (size_t)nqp);
    double Nq[NQ][6], Nxq[NQ][6], Neq[NQ][6];
    for (int q = 0; q < NQ; q++) p2_shape(DQ[q][0], DQ[q][1], Nq[q], Nxq[q], Neq[q]);
    for (int64_t e = 0; e < ne; e++) {
        const int64_t* t = m->tri + 6 * e;
        for (int q = 0; q < NQ; q++) {
            double x = 0, y = 0;
            for (int k = 0; k < 6; k++) { x += Nq[q][k] * m->x[2 * t[k]]; y += Nq[q][k] * m->x[2 * t[k] + 1]; }
            QX[e * NQ + q] = x; QY[e * NQ + q] = y;
        }
    }
    double* fd = NULL; double* fbx = NULL; double* fby = NULL; double* fr = NULL;
    bool ok = true;
    if (P->d)  { fd = malloc(sizeof(double) * (size_t)nqp);  ok = ok && fem_field(P->d, P->xv, P->yv, QX, QY, nqp, fd); }
    if (P->bx) { fbx = malloc(sizeof(double) * (size_t)nqp); ok = ok && fem_field(P->bx, P->xv, P->yv, QX, QY, nqp, fbx); }
    if (P->by) { fby = malloc(sizeof(double) * (size_t)nqp); ok = ok && fem_field(P->by, P->xv, P->yv, QX, QY, nqp, fby); }
    if (want_load && P->r) { fr = malloc(sizeof(double) * (size_t)nqp); ok = ok && fem_field(P->r, P->xv, P->yv, QX, QY, nqp, fr); }
    free(QX); free(QY);
    if (!ok) { *why = "a coefficient or the source did not evaluate to a real number"; free(fd); free(fbx); free(fby); free(fr); return false; }
    bool any_b = false, any_d = false;
    for (int64_t i = 0; fbx && i < nqp && !any_b; i++) any_b = fbx[i] != 0.0;
    for (int64_t i = 0; fby && i < nqp && !any_b; i++) any_b = fby[i] != 0.0;
    for (int64_t i = 0; fd && i < nqp && !any_d; i++) any_d = fd[i] != 0.0;
    S->symmetric = !any_b;

    for (int64_t e = 0; e < ne; e++) {
        const int64_t* t = m->tri + 6 * e;
        double X[6], Y[6];
        for (int k = 0; k < 6; k++) { X[k] = m->x[2 * t[k]]; Y[k] = m->x[2 * t[k] + 1]; }
        double Ke[6][6], Me[6][6], Fe[6];
        memset(Ke, 0, sizeof Ke); memset(Me, 0, sizeof Me); memset(Fe, 0, sizeof Fe);
        for (int q = 0; q < NQ; q++) {
            double gx[6], gy[6];
            double det = p2_jac(X, Y, Nxq[q], Neq[q], gx, gy);
            if (!(det > 0)) { ok = false; break; }
            double w = 0.5 * DQ[q][2] * det;
            int64_t qi = e * NQ + q;
            double dv = (any_d && fd) ? fd[qi] : 0.0;
            double bxv = (any_b && fbx) ? fbx[qi] : 0.0, byv = (any_b && fby) ? fby[qi] : 0.0;
            const double* N = Nq[q];
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    double a = P->c11 * gx[j] * gx[i] + P->c12 * (gx[j] * gy[i] + gy[j] * gx[i])
                             + P->c22 * gy[j] * gy[i];
                    a += (bxv * gx[j] + byv * gy[j]) * N[i] + dv * N[j] * N[i];
                    Ke[i][j] += w * a;
                    if (want_mass) Me[i][j] += w * N[j] * N[i];
                }
                if (fr) Fe[i] -= w * fr[qi] * N[i];
            }
        }
        if (!ok) break;
        for (int i = 0; i < 6; i++) {
            int64_t gi = S->map[t[i]];
            if (gi < 0) continue;
            for (int j = 0; j < 6; j++) {
                int64_t gj = S->map[t[j]];
                if (gj >= 0) {
                    coo_add(&S->K, gi, gj, Ke[i][j]);
                    if (want_mass) coo_add(&S->M, gi, gj, Me[i][j]);
                } else if (want_load) {
                    S->F[gi] -= Ke[i][j] * S->ufix[t[j]];
                }
            }
            if (want_load) S->F[gi] += Fe[i];
        }
    }
    free(fd); free(fbx); free(fby); free(fr);
    if (!ok) { *why = "the mesh has an inverted element"; return false; }

    /* NeumannValue terms on boundary edges: s Int g v = s Int (g0 + g1 u) v */
    if (P->nnv > 0) {
        static const double GL[5][2] = {
            { 0.0469100770306680, 0.1184634425280945 }, { 0.2307653449471585, 0.2393143352496832 },
            { 0.5, 0.2844444444444444 }, { 0.7692346550528415, 0.2393143352496832 },
            { 0.9530899229693320, 0.1184634425280945 } };
        for (int64_t b = 0; b < m->nb; b++) {
            int64_t nd3[3] = { m->bedge[3 * b], m->bedge[3 * b + 2], m->bedge[3 * b + 1] };
            double ex[3], ey[3];
            for (int k = 0; k < 3; k++) { ex[k] = m->x[2 * nd3[k]]; ey[k] = m->x[2 * nd3[k] + 1]; }
            for (int k = 0; k < P->nnv; k++) {
                int pv = fem_pred(P->nv_pred[k], P->xv, P->yv, ex[1], ey[1]);
                if (pv < 0) { *why = "a NeumannValue predicate did not evaluate to True or False"; return false; }
                if (!pv) continue;
                double qx[5], qy[5], wq[5], L[5][3];
                for (int q = 0; q < 5; q++) {
                    double s = GL[q][0];
                    L[q][0] = (1 - s) * (1 - 2 * s); L[q][1] = 4 * s * (1 - s); L[q][2] = s * (2 * s - 1);
                    double d0 = 4 * s - 3, d1 = 4 - 8 * s, d2 = 4 * s - 1;
                    double tx = d0 * ex[0] + d1 * ex[1] + d2 * ex[2], ty = d0 * ey[0] + d1 * ey[1] + d2 * ey[2];
                    qx[q] = L[q][0] * ex[0] + L[q][1] * ex[1] + L[q][2] * ex[2];
                    qy[q] = L[q][0] * ey[0] + L[q][1] * ey[1] + L[q][2] * ey[2];
                    wq[q] = GL[q][1] * sqrt(tx * tx + ty * ty);
                }
                double g0[5], g1[5];
                if (!fem_field(P->nv_g0[k], P->xv, P->yv, qx, qy, 5, g0)
                    || !fem_field(P->nv_g1[k], P->xv, P->yv, qx, qy, 5, g1)) {
                    *why = "a NeumannValue did not evaluate to a real number"; return false;
                }
                double s = P->nv_s[k];
                for (int i = 0; i < 3; i++) {
                    int64_t gi = S->map[nd3[i]];
                    if (gi < 0) continue;
                    for (int q = 0; q < 5; q++) {
                        if (want_load) S->F[gi] -= s * wq[q] * g0[q] * L[q][i];
                        for (int j = 0; j < 3; j++) {
                            double a = s * wq[q] * g1[q] * L[q][j] * L[q][i];
                            if (a == 0.0) continue;
                            int64_t gj = S->map[nd3[j]];
                            if (gj >= 0) coo_add(&S->K, gi, gj, a);
                            else if (want_load) S->F[gi] -= a * S->ufix[nd3[j]];
                        }
                    }
                }
                if (g1[0] != 0 || g1[2] != 0 || g1[4] != 0) S->symmetric = S->symmetric && true;
            }
        }
    }
    return true;
}

/* ------------------------------------------------------------------ */
/*  Mesh functions (InterpolatingFunction over the mesh)              */
/* ------------------------------------------------------------------ */

/* NDSolve`FEM`Mesh[coords (n x 2), elements (ne x 6, 0-based, as reals),
 * curved (ne), region] -- one expression shared by every function on the
 * mesh; the region (as given, evaluated) lets NIntegrate over that region
 * integrate element by element. */
static Expr* fem_mesh_expr(const FemMesh* m, Expr* region) {
    double* c = malloc(sizeof(double) * 2 * (size_t)m->nn);
    memcpy(c, m->x, sizeof(double) * 2 * (size_t)m->nn);
    double* t = malloc(sizeof(double) * 6 * (size_t)m->ne);
    for (int64_t i = 0; i < 6 * m->ne; i++) t[i] = (double)m->tri[i];
    double* cv = malloc(sizeof(double) * (size_t)m->ne);
    for (int64_t i = 0; i < m->ne; i++) cv[i] = m->curved[i] ? 1.0 : 0.0;
    Expr* a[4] = { fem_packed(c, 2, m->nn, 2), fem_packed(t, 2, m->ne, 6), fem_packed(cv, 1, m->ne, 0),
                   eval_and_free(expr_copy(region)) };
    return expr_new_function(expr_new_symbol(SYM_FEMMesh), a, 4);
}

static Expr* fem_function_expr(Expr* meshE, const FemMesh* m, const double* vals) {
    double* v = malloc(sizeof(double) * (size_t)m->nn);
    memcpy(v, vals, sizeof(double) * (size_t)m->nn);
    Expr* table_args[2] = { expr_copy(meshE), fem_packed(v, 1, m->nn, 0) };
    Expr* table = expr_new_function(expr_new_symbol(SYM_FEMFunction), table_args, 2);
    Expr* r0[2] = { expr_new_real(m->bbox[0]), expr_new_real(m->bbox[1]) };
    Expr* r1[2] = { expr_new_real(m->bbox[2]), expr_new_real(m->bbox[3]) };
    Expr* dom[2] = { fem_list(r0, 2), fem_list(r1, 2) };
    Expr* ia[2] = { fem_list(dom, 2), table };
    return expr_new_function(expr_new_symbol(SYM_InterpolatingFunction), ia, 2);
}

bool fem_is_mesh_function(const Expr* table) {
    return head_is((Expr*)table, SYM_FEMFunction) && table->data.function.arg_count == 2
        && head_is(table->data.function.args[0], SYM_FEMMesh)
        && table->data.function.args[0]->data.function.arg_count == 4;
}

/* Point location: a uniform bucket grid over the elements' bounding boxes. */
typedef struct {
    const Expr* key; Expr* hold;
    int64_t nn, ne;
    const double* X;      /* coords (borrowed from hold) */
    double* T;            /* elements as doubles (borrowed) */
    const double* C;      /* curved flags */
    double x0, y0, dx, dy; int nbx, nby;
    int64_t* start; int64_t* list;
} FemLocator;

#define FEM_CACHE 6
static FemLocator* g_loc[FEM_CACHE];
static int g_loc_next = 0;

static const double* fem_buffer(const Expr* e, int64_t want) {
    if (!e || e->type != EXPR_NDARRAY || e->data.ndarray.dtype != NDT_FLOAT64) return NULL;
    int64_t tot = 1;
    for (int k = 0; k < e->data.ndarray.rank; k++) tot *= e->data.ndarray.dims[k];
    return tot == want ? (const double*)e->data.ndarray.data : NULL;
}

static void fem_loc_free(FemLocator* L) {
    if (!L) return;
    expr_free(L->hold); free(L->start); free(L->list); free(L);
}

static FemLocator* fem_locator(const Expr* meshE) {
    for (int i = 0; i < FEM_CACHE; i++)
        if (g_loc[i] && g_loc[i]->key == meshE) return g_loc[i];
    const Expr* ce = meshE->data.function.args[0];
    const Expr* te = meshE->data.function.args[1];
    const Expr* ve = meshE->data.function.args[2];
    if (ce->type != EXPR_NDARRAY || ce->data.ndarray.rank != 2 || te->type != EXPR_NDARRAY
        || te->data.ndarray.rank != 2 || te->data.ndarray.dims[1] != 6) return NULL;
    int64_t nn = ce->data.ndarray.dims[0], ne = te->data.ndarray.dims[0];
    const double* X = fem_buffer(ce, 2 * nn);
    const double* T = fem_buffer(te, 6 * ne);
    const double* C = fem_buffer(ve, ne);
    if (!X || !T || !C || ne < 1) return NULL;
    FemLocator* L = calloc(1, sizeof *L);
    L->key = meshE; L->hold = expr_copy((Expr*)meshE);
    L->nn = nn; L->ne = ne; L->X = X; L->T = (double*)T; L->C = C;
    double xmin = INFINITY, xmax = -INFINITY, ymin = INFINITY, ymax = -INFINITY;
    for (int64_t i = 0; i < nn; i++) {
        double x = X[2 * i], y = X[2 * i + 1];
        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
    }
    double w = xmax - xmin, h = ymax - ymin;
    double cells = (double)ne;
    int nbx = (int)ceil(sqrt(cells * w / (h > 0 ? h : 1))), nby;
    if (nbx < 1) nbx = 1;
    if (nbx > 4096) nbx = 4096;
    nby = (int)ceil(cells / nbx);
    if (nby < 1) nby = 1;
    if (nby > 4096) nby = 4096;
    L->nbx = nbx; L->nby = nby;
    L->x0 = xmin; L->y0 = ymin;
    L->dx = (w > 0 ? w : 1) / nbx * (1 + 1e-12); L->dy = (h > 0 ? h : 1) / nby * (1 + 1e-12);
    int64_t ncell = (int64_t)nbx * nby;
    L->start = calloc((size_t)ncell + 1, sizeof(int64_t));
    for (int pass = 0; pass < 2; pass++) {
        if (pass == 1) {
            for (int64_t c = 0; c < ncell; c++) L->start[c + 1] += L->start[c];
            L->list = malloc(sizeof(int64_t) * (size_t)(L->start[ncell] ? L->start[ncell] : 1));
        }
        int64_t* fill = pass == 1 ? calloc((size_t)ncell, sizeof(int64_t)) : NULL;
        for (int64_t e = 0; e < ne; e++) {
            double bx0 = INFINITY, bx1 = -INFINITY, by0 = INFINITY, by1 = -INFINITY;
            for (int k = 0; k < 6; k++) {
                int64_t v = (int64_t)T[6 * e + k];
                double x = X[2 * v], y = X[2 * v + 1];
                if (x < bx0) bx0 = x;
                if (x > bx1) bx1 = x;
                if (y < by0) by0 = y;
                if (y > by1) by1 = y;
            }
            int i0 = (int)floor((bx0 - L->x0) / L->dx), i1 = (int)floor((bx1 - L->x0) / L->dx);
            int j0 = (int)floor((by0 - L->y0) / L->dy), j1 = (int)floor((by1 - L->y0) / L->dy);
            if (i0 < 0) i0 = 0;
            if (j0 < 0) j0 = 0;
            if (i1 >= nbx) i1 = nbx - 1;
            if (j1 >= nby) j1 = nby - 1;
            for (int j = j0; j <= j1; j++)
                for (int i = i0; i <= i1; i++) {
                    int64_t c = (int64_t)j * nbx + i;
                    if (pass == 0) L->start[c + 1]++;
                    else L->list[L->start[c] + fill[c]++] = e;
                }
        }
        free(fill);
    }
    fem_loc_free(g_loc[g_loc_next]);
    g_loc[g_loc_next] = L;
    g_loc_next = (g_loc_next + 1) % FEM_CACHE;
    return L;
}

/* Reference coordinates of (px, py) in element e: affine from the vertices,
 * then Newton on the isoparametric map for a curved element. Returns the
 * smallest barycentric coordinate (>= 0 inside). */
static double fem_local(const FemLocator* L, int64_t e, double px, double py, double* xi, double* eta) {
    const double* t = L->T + 6 * e;
    double X[6], Y[6];
    for (int k = 0; k < 6; k++) { int64_t v = (int64_t)t[k]; X[k] = L->X[2 * v]; Y[k] = L->X[2 * v + 1]; }
    double a11 = X[1] - X[0], a12 = X[2] - X[0], a21 = Y[1] - Y[0], a22 = Y[2] - Y[0];
    double det = a11 * a22 - a12 * a21;
    double rx = px - X[0], ry = py - Y[0];
    double s = ( a22 * rx - a12 * ry) / det;
    double r = (-a21 * rx + a11 * ry) / det;
    if (L->C[e] != 0.0) {
        for (int it = 0; it < 12; it++) {
            double N[6], Nx[6], Ne[6];
            p2_shape(s, r, N, Nx, Ne);
            double fx = -px, fy = -py, j11 = 0, j12 = 0, j21 = 0, j22 = 0;
            for (int k = 0; k < 6; k++) {
                fx += N[k] * X[k]; fy += N[k] * Y[k];
                j11 += Nx[k] * X[k]; j12 += Ne[k] * X[k]; j21 += Nx[k] * Y[k]; j22 += Ne[k] * Y[k];
            }
            double dj = j11 * j22 - j12 * j21;
            if (dj == 0) break;
            double ds = ( j22 * fx - j12 * fy) / dj, dr = (-j21 * fx + j11 * fy) / dj;
            s -= ds; r -= dr;
            if (fabs(ds) + fabs(dr) < 1e-15) break;
        }
    }
    *xi = s; *eta = r;
    double l0 = 1 - s - r, mn = l0;
    if (s < mn) mn = s;
    if (r < mn) mn = r;
    return mn;
}

/* Value (dx, dy derivative orders) of the mesh function at (px, py). */
static bool fem_eval_point(const FemLocator* L, const double* V, int dx, int dy,
                           double px, double py, double* out, bool* outside) {
    int i = (int)floor((px - L->x0) / L->dx), j = (int)floor((py - L->y0) / L->dy);
    int64_t best = -1; double bestm = -INFINITY, bxi = 0, beta = 0;
    for (int ring = 0; ring <= 2 && bestm < -1e-10; ring++) {
        for (int jj = j - ring; jj <= j + ring; jj++)
            for (int ii = i - ring; ii <= i + ring; ii++) {
                if (ring > 0 && jj > j - ring && jj < j + ring && ii > i - ring && ii < i + ring) continue;
                if (ii < 0 || jj < 0 || ii >= L->nbx || jj >= L->nby) continue;
                int64_t c = (int64_t)jj * L->nbx + ii;
                for (int64_t p = L->start[c]; p < L->start[c + 1]; p++) {
                    int64_t e = L->list[p];
                    double xi, eta;
                    double mn = fem_local(L, e, px, py, &xi, &eta);
                    if (mn > bestm) { bestm = mn; best = e; bxi = xi; beta = eta; }
                    if (mn >= -1e-12) goto found;
                }
            }
    }
found:
    if (best < 0) {
        /* far outside: the nearest element by a full scan */
        for (int64_t e = 0; e < L->ne; e++) {
            double xi, eta, mn = fem_local(L, e, px, py, &xi, &eta);
            if (mn > bestm) { bestm = mn; best = e; bxi = xi; beta = eta; }
        }
        if (best < 0) return false;
    }
    *outside = bestm < -1e-8;
    const double* t = L->T + 6 * best;
    double N[6], Nx[6], Ne[6];
    p2_shape(bxi, beta, N, Nx, Ne);
    double U[6], X[6], Y[6];
    for (int k = 0; k < 6; k++) {
        int64_t v = (int64_t)t[k];
        U[k] = V[v]; X[k] = L->X[2 * v]; Y[k] = L->X[2 * v + 1];
    }
    if (dx == 0 && dy == 0) {
        double s = 0; for (int k = 0; k < 6; k++) s += N[k] * U[k];
        *out = s; return true;
    }
    double gx[6], gy[6];
    p2_jac(X, Y, Nx, Ne, gx, gy);
    if (dx + dy == 1) {
        double s = 0;
        for (int k = 0; k < 6; k++) s += (dx ? gx[k] : gy[k]) * U[k];
        *out = s; return true;
    }
    if (dx + dy == 2 && L->C[best] == 0.0) {
        /* affine element: constant second derivatives. N_k'' in reference
         * coordinates, mapped by the constant inverse Jacobian. */
        static const double Hss[6] = { 4, 4, 0, -8, 0, 0 };
        static const double Hrr[6] = { 4, 0, 4, 0, 0, -8 };
        static const double Hsr[6] = { 4, 0, 0, -4, 4, -4 };
        double a11 = X[1] - X[0], a12 = X[2] - X[0], a21 = Y[1] - Y[0], a22 = Y[2] - Y[0];
        double det = a11 * a22 - a12 * a21;
        /* d(s,r)/d(x,y) */
        double sx = a22 / det, sy = -a12 / det, rx = -a21 / det, ry = a11 / det;
        double s = 0;
        for (int k = 0; k < 6; k++) {
            double h;
            if (dx == 2) h = Hss[k] * sx * sx + 2 * Hsr[k] * sx * rx + Hrr[k] * rx * rx;
            else if (dy == 2) h = Hss[k] * sy * sy + 2 * Hsr[k] * sy * ry + Hrr[k] * ry * ry;
            else h = Hss[k] * sx * sy + Hsr[k] * (sx * ry + sy * rx) + Hrr[k] * rx * ry;
            s += h * U[k];
        }
        *out = s; return true;
    }
    if (dx + dy >= 2) { *out = 0.0; return dx + dy > 2 || true; }
    return false;
}

Expr* fem_interp_apply(Expr* ifun, Expr** call_args, size_t argc) {
    Expr* table = ifun->data.function.args[1];
    Expr* meshE = table->data.function.args[0];
    if (argc == 1 && call_args[0]->type == EXPR_STRING) {
        const char* p = call_args[0]->data.string;
        if (strcmp(p, "Domain") == 0) return expr_copy(ifun->data.function.args[0]);
        if (strcmp(p, "Coordinates") == 0 || strcmp(p, "Grid") == 0)
            return expr_copy(meshE->data.function.args[0]);
        if (strcmp(p, "ValuesOnGrid") == 0) return expr_copy(table->data.function.args[1]);
        if (strcmp(p, "ElementMesh") == 0) return expr_copy(meshE);
        if (strcmp(p, "Properties") == 0) {
            const char* ps[5] = { "Coordinates", "Domain", "ElementMesh", "Grid", "ValuesOnGrid" };
            Expr* it[5];
            for (int k = 0; k < 5; k++) it[k] = expr_new_string(ps[k]);
            return fem_list(it, 5);
        }
        return NULL;
    }
    if (argc != 2) return NULL;
    int dx = 0, dy = 0;
    if (ifun->data.function.arg_count >= 3) {
        Expr* ds = ifun->data.function.args[2];
        if (!head_is(ds, SYM_List) || ds->data.function.arg_count != 2
            || ds->data.function.args[0]->type != EXPR_INTEGER || ds->data.function.args[1]->type != EXPR_INTEGER)
            return NULL;
        dx = (int)ds->data.function.args[0]->data.integer;
        dy = (int)ds->data.function.args[1]->data.integer;
        if (dx < 0 || dy < 0 || dx + dy > 2) return NULL;
    }
    FemLocator* L = fem_locator(meshE);
    if (!L) return NULL;
    const double* V = fem_buffer(table->data.function.args[1], L->nn);
    if (!V) return NULL;
    double px, py;
    bool outside = false;
    if (nd_to_double(call_args[0], &px) && nd_to_double(call_args[1], &py)) {
        if (!isfinite(px) || !isfinite(py)) return NULL;
        double v;
        if (!fem_eval_point(L, V, dx, dy, px, py, &v, &outside)) return NULL;
        if (outside)
            mth_message("InterpolatingFunction", "femdmval",
                        "Input value {%.6g, %.6g} lies outside the range of data in the interpolating "
                        "function. Extrapolation will be used.", px, py);
        return expr_new_real(v);
    }
    /* vectors of points (packed or plain numeric Lists of equal length) */
    Expr* a = call_args[0]; Expr* b = call_args[1];
    int64_t n = -1;
    if (a->type == EXPR_NDARRAY && a->data.ndarray.rank == 1) n = a->data.ndarray.dims[0];
    else if (head_is(a, SYM_List)) n = (int64_t)a->data.function.arg_count;
    if (n < 1) return NULL;
    double* xs = malloc(sizeof(double) * (size_t)n);
    double* ys = malloc(sizeof(double) * (size_t)n);
    if (!fem_read_reals(a, n, xs) || !fem_read_reals(b, n, ys)) { free(xs); free(ys); return NULL; }
    double* out = malloc(sizeof(double) * (size_t)n);
    bool any_out = false;
    for (int64_t i = 0; i < n; i++) {
        if (!fem_eval_point(L, V, dx, dy, xs[i], ys[i], &out[i], &outside)) { free(xs); free(ys); free(out); return NULL; }
        any_out = any_out || outside;
    }
    if (any_out)
        mth_message("InterpolatingFunction", "femdmval",
                    "Input values lie outside the range of data in the interpolating function. "
                    "Extrapolation will be used.");
    free(xs); free(ys);
    return fem_packed(out, 1, n, 0);
}

/* ------------------------------------------------------------------ */
/*  Problem setup shared by NDSolve and NDEigensystem                  */
/* ------------------------------------------------------------------ */

typedef struct {
    const char* u; bool applied;
    const char* xv; const char* yv;
    Expr* region;          /* borrowed */
    double max_area, max_bedge;
} FemSpec;

/* u or u[x, y]; {x, y} in reg */
static bool fem_read_spec(Expr* uspec, Expr* rspec, FemSpec* F, const char** why) {
    memset(F, 0, sizeof *F);
    if (!fem_is_region_spec(rspec) || rspec->data.function.args[0]->data.function.arg_count != 2) {
        *why = "the spatial specification must be {x, y} ∈ region"; return false;
    }
    Expr* vs = rspec->data.function.args[0];
    if (vs->data.function.args[0]->type != EXPR_SYMBOL || vs->data.function.args[1]->type != EXPR_SYMBOL) {
        *why = "the spatial variables must be symbols"; return false;
    }
    F->xv = vs->data.function.args[0]->data.symbol.name;
    F->yv = vs->data.function.args[1]->data.symbol.name;
    F->region = rspec->data.function.args[1];
    Expr* u = uspec;
    if (head_is(u, SYM_List) && u->data.function.arg_count == 1) u = u->data.function.args[0];
    if (u->type == EXPR_SYMBOL) { F->u = u->data.symbol.name; F->applied = false; }
    else if (u->type == EXPR_FUNCTION && u->data.function.head->type == EXPR_SYMBOL
             && u->data.function.arg_count == 2) {
        F->u = u->data.function.head->data.symbol.name; F->applied = true;
    } else { *why = "the unknown must be u or u[x, y]"; return false; }
    return true;
}

static void fem_mesh_options(Expr** opts, size_t nopt, FemSpec* F) {
    F->max_area = -1; F->max_bedge = -1;
    for (size_t i = 0; i < nopt; i++) {
        double v;
        if (fem_find_opt(opts[i], "MaxCellMeasure", &v) && v > 0) F->max_area = v;
        if (fem_find_opt(opts[i], "MaxBoundaryCellMeasure", &v) && v > 0) F->max_bedge = v;
    }
}

static bool fem_build_mesh(const FemSpec* F, const FemPde* P, FemMesh* mesh, const char** why) {
    FemGeom g;
    if (!fem_region(F->region, &g, why)) return false;
    /* boundary-condition type of each input edge (at its midpoint), so the
     * mesher grades the corners where a Dirichlet edge meets a natural one */
    g.edge_dir = malloc((size_t)g.nv);
    for (int k = 0; k < g.nloop; k++) {
        int s0 = g.loop_start[k], n = g.loop_n[k];
        for (int i = 0; i < n; i++) {
            int a = s0 + i, b = s0 + (i + 1) % n;
            double mx = 0.5 * (g.xy[2 * a] + g.xy[2 * b]), my = 0.5 * (g.xy[2 * a + 1] + g.xy[2 * b + 1]);
            if (g.edge_arc[a] >= 0) {
                const FemArc* A = &g.arc[g.edge_arc[a]];
                double ta = atan2(g.xy[2 * a + 1] - A->cy, g.xy[2 * a] - A->cx);
                double tb = atan2(g.xy[2 * b + 1] - A->cy, g.xy[2 * b] - A->cx);
                double d = tb - ta;
                while (d > M_PI) d -= 2 * M_PI;
                while (d <= -M_PI) d += 2 * M_PI;
                mx = A->cx + A->r * cos(ta + 0.5 * d); my = A->cy + A->r * sin(ta + 0.5 * d);
            }
            unsigned char dir = 0;
            for (int c = 0; c < P->ndc && !dir; c++)
                if (fem_pred(P->dc_pred[c], P->xv, P->yv, mx, my) == 1) dir = 1;
            g.edge_dir[a] = dir;
        }
    }
    double area = fabs(fem_geom_area(&g));
    FemMeshOpts mo;
    mo.max_area = F->max_area > 0 ? F->max_area : area / 1500.0;
    mo.max_bedge = F->max_bedge;
    mo.max_vertices = 4000000;
    mo.grade_corners = true;
    bool ok = fem_mesh_build(&g, &mo, mesh, why);
    fem_geom_free(&g);
    return ok;
}

static bool fem_is_option(const Expr* e) {
    return (head_is((Expr*)e, SYM_Rule) || head_is((Expr*)e, SYM_RuleDelayed)) && e->data.function.arg_count == 2;
}

/* ------------------------------------------------------------------ */
/*  NDSolve                                                            */
/* ------------------------------------------------------------------ */

Expr* fem_ndsolve(Expr* res) {
    Expr** A = res->data.function.args;
    size_t argc = res->data.function.arg_count;
    const char* why = NULL;
    FemSpec F;
    if (argc < 3 || !fem_read_spec(A[1], A[2], &F, &why)) {
        if (why) mth_message("NDSolve", "femspec", "%s.", why);
        return NULL;
    }
    fem_mesh_options(A + 3, argc - 3, &F);
    FemPde P; memset(&P, 0, sizeof P);
    P.u = F.u; P.xv = F.xv; P.yv = F.yv;
    Expr** eqs; size_t neq;
    Expr* eqlist = A[0];
    if (head_is(eqlist, SYM_List)) { eqs = eqlist->data.function.args; neq = eqlist->data.function.arg_count; }
    else { eqs = &A[0]; neq = 1; }
    Expr* pde = NULL;
    bool ok = true;
    for (size_t i = 0; i < neq && ok; i++) {
        Expr* e = eqs[i];
        if (head_is(e, SYM_DirichletCondition)) ok = fem_read_dirichlet(&P, e, &why);
        else if (head_is(e, SYM_Equal) && e->data.function.arg_count == 2) {
            if (pde) { ok = false; why = "only one partial differential equation is supported"; }
            else pde = e;
        } else { ok = false; why = "each item must be an equation or a DirichletCondition"; }
    }
    if (ok && !pde) { ok = false; why = "there is no partial differential equation"; }
    if (ok) {
        Expr* R = eval_and_free(nd_call2(SYM_Subtract, expr_copy(pde->data.function.args[0]),
                                         expr_copy(pde->data.function.args[1])));
        ok = fem_read_operator(&P, R, false, &why);
        expr_free(R);
    }
    FemSystem S; memset(&S, 0, sizeof S);
    if (ok) ok = fem_build_mesh(&F, &P, &S.mesh, &why);
    if (ok) ok = fem_assemble(&P, &S, false, true, &why);
    Expr* out = NULL;
    if (ok) {
        double* U = malloc(sizeof(double) * (size_t)S.mesh.nn);
        double* xf = malloc(sizeof(double) * (size_t)(S.nfree ? S.nfree : 1));
        bool solved = true;
        if (S.nfree > 0) {
            Csr K; csr_from_coo(&S.K, 0.0, NULL, &K);
            double* xyf = fem_free_xy(&S);
            FemSolver fs;
            solved = fem_solver_factor(&fs, &K, S.symmetric, xyf) && fem_solver_solve(&fs, S.F, xf);
            fem_solver_free(&fs);
            free(xyf);
            csr_free(&K);
            if (!solved) why = "the finite element system is singular (a pure Neumann problem has no unique solution)";
        }
        if (solved) {
            for (int64_t i = 0; i < S.mesh.nn; i++) U[i] = S.map[i] >= 0 ? xf[S.map[i]] : S.ufix[i];
            Expr* meshE = fem_mesh_expr(&S.mesh, F.region);
            Expr* ifun = fem_function_expr(meshE, &S.mesh, U);
            expr_free(meshE);
            Expr* lhs;
            if (F.applied) {
                Expr* ar[2] = { expr_new_symbol(F.xv), expr_new_symbol(F.yv) };
                lhs = expr_new_function(expr_new_symbol(F.u), ar, 2);
                Expr* ar2[2] = { expr_new_symbol(F.xv), expr_new_symbol(F.yv) };
                ifun = expr_new_function(ifun, ar2, 2);
            } else lhs = expr_new_symbol(F.u);
            Expr* rule = nd_call2(SYM_Rule, lhs, ifun);
            Expr* inner = fem_list(&rule, 1);
            out = fem_list(&inner, 1);
        } else ok = false;
        free(U); free(xf);
    }
    if (!ok && why) mth_message("NDSolve", "femfail", "%s.", why);
    fem_system_free(&S);
    fem_pde_free(&P);
    return out;
}

/* ------------------------------------------------------------------ */
/*  NDEigenvalues / NDEigensystem                                      */
/* ------------------------------------------------------------------ */

typedef struct { const FemSolver* F; const Csr* M; double* tmp; } FemOp;

static bool fem_op(void* ctx, const double* x, double* y) {
    FemOp* o = ctx;
    csr_mv(o->M, x, o->tmp);
    return fem_solver_solve(o->F, o->tmp, y);
}

static Expr* fem_eigen(Expr* res, bool want_vectors) {
    const char* head = want_vectors ? "NDEigensystem" : "NDEigenvalues";
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 4) return NULL;
    Expr** A = res->data.function.args;
    size_t argc = res->data.function.arg_count;
    const char* why = NULL;
    FemSpec F;
    if (!fem_read_spec(A[1], A[2], &F, &why)) { mth_message(head, "femspec", "%s.", why); return NULL; }
    if (A[3]->type != EXPR_INTEGER || A[3]->data.integer < 1) {
        mth_message(head, "intp", "Positive integer expected at position 4 in %s.", head);
        return NULL;
    }
    int nev = (int)A[3]->data.integer;
    for (size_t i = 4; i < argc; i++) {
        if (!fem_is_option(A[i])) return NULL;
        const Expr* l = A[i]->data.function.args[0];
        const char* nm = l->type == EXPR_SYMBOL ? l->data.symbol.name : NULL;
        if (nm != SYM_Method && nm != intern_symbol("DependentVariables"))
            mth_message(head, "woptx", "Unknown option %s is not one of {DependentVariables -> Automatic, "
                        "Method -> Automatic}; it is ignored.", nm ? nm : "given");
    }
    fem_mesh_options(A + 4, argc - 4, &F);
    FemPde P; memset(&P, 0, sizeof P);
    P.u = F.u; P.xv = F.xv; P.yv = F.yv;
    Expr* op = NULL;
    bool ok = true;
    Expr* spec = A[0];
    if (head_is(spec, SYM_List)) {
        for (size_t i = 0; i < spec->data.function.arg_count && ok; i++) {
            Expr* e = spec->data.function.args[i];
            if (head_is(e, SYM_DirichletCondition)) ok = fem_read_dirichlet(&P, e, &why);
            else if (!op) op = e;
            else { ok = false; why = "only one operator is supported"; }
        }
    } else op = spec;
    if (ok && !op) { ok = false; why = "there is no operator"; }
    if (ok) {
        Expr* R = eval_and_free(expr_copy(op));
        ok = fem_read_operator(&P, R, true, &why);
        expr_free(R);
    }
    for (int k = 0; k < P.ndc && ok; k++) {
        double z;
        if (!(fem_num(P.dc_val[k], &z) && z == 0.0)) { ok = false; why = "an eigenvalue problem needs homogeneous DirichletConditions (u == 0)"; }
    }
    FemSystem S; memset(&S, 0, sizeof S);
    if (ok) ok = fem_build_mesh(&F, &P, &S.mesh, &why);
    if (ok) ok = fem_assemble(&P, &S, true, false, &why);
    if (ok && nev > S.nfree) { ok = false; why = "more eigenvalues were requested than the discretization has"; }
    Expr* out = NULL;
    if (ok) {
        int64_t n = S.nfree;
        Csr K, M, Asig;
        csr_from_coo(&S.K, 0.0, NULL, &K);
        csr_from_coo(&S.M, 0.0, NULL, &M);
        /* shift: 0 when K is definite (Dirichlet somewhere), else a small
         * negative multiple of the spectrum's scale */
        double sigma = 0.0;
        FemSolver fsol; memset(&fsol, 0, sizeof fsol);
        bool fac = false;
        double* xyf = fem_free_xy(&S);
        /* no Dirichlet node: K is singular (the constants), so shift away
         * from 0 at once -- a factor of the singular K can succeed in
         * rounding and swamps the Krylov tolerance with a huge 1/theta */
        bool any_fixed = S.nfree < S.mesh.nn;
        for (int attempt = any_fixed ? 0 : 1; attempt < 2 && !fac; attempt++) {
            if (attempt == 1) {
                double dk = 0, dm = 0;
                for (int64_t i = 0; i < n; i++)
                    for (int64_t p = K.rp[i]; p < K.rp[i + 1]; p++)
                        if (K.ci[p] == i) { dk += fabs(K.x[p]); }
                for (int64_t i = 0; i < n; i++)
                    for (int64_t p = M.rp[i]; p < M.rp[i + 1]; p++)
                        if (M.ci[p] == i) { dm += fabs(M.x[p]); }
                /* trace(K)/trace(M)/n: of the order of the smallest
                 * nonzero eigenvalue */
                sigma = -(dk / (dm > 0 ? dm : 1)) / (double)n;
                if (sigma == 0.0) sigma = -1.0;
            }
            csr_from_coo(&S.K, -sigma, &S.M, &Asig);
            fac = fem_solver_factor(&fsol, &Asig, S.symmetric, xyf);
            if (fac && fsol.F && spf_badly_conditioned(fsol.F) && attempt == 0) { fem_solver_free(&fsol); fac = false; }
            csr_free(&Asig);
        }
        free(xyf);
        if (!fac) { ok = false; why = "the shifted operator could not be factored"; }
        KsResult kr; memset(&kr, 0, sizeof kr);
        if (ok) {
            FemOp o = { &fsol, &M, malloc(sizeof(double) * (size_t)n) };
            KsProblem kp; memset(&kp, 0, sizeof kp);
            kp.n = n; kp.cplx = false; kp.hermitian = false; kp.op = fem_op; kp.ctx = &o;
            int want = nev + (nev < 4 ? 2 : 3);       /* guard pairs for clusters */
            if (want > n) want = (int)n;
            kp.nev = want; kp.which = KS_LM; kp.ncv = 0; kp.maxit = 0; kp.tol = 1e-13;
            kp.v0 = NULL; kp.want_vectors = true;
            KsStatus ks = ks_solve(&kp, &kr);
            free(o.tmp);
            if ((ks != KS_OK && ks != KS_MAXIT) || kr.nconv < nev) {
                ok = false; why = "the eigensolver did not converge";
            }
        }
        fem_solver_free(&fsol);
        if (ok) {
            int nc = kr.nconv;
            double* lam = malloc(sizeof(double) * (size_t)nc);
            double* vec = malloc(sizeof(double) * (size_t)nc * (size_t)n);
            double* Kx = malloc(sizeof(double) * (size_t)n);
            double* Mx = malloc(sizeof(double) * (size_t)n);
            for (int k = 0; k < nc; k++) {
                double* v = vec + (size_t)k * (size_t)n;
                for (int64_t i = 0; i < n; i++)
                    v[i] = kr.vec_cplx ? kr.vec[2 * ((size_t)k * (size_t)n + (size_t)i)] : kr.vec[(size_t)k * (size_t)n + (size_t)i];
                csr_mv(&K, v, Kx); csr_mv(&M, v, Mx);
                double num = 0, den = 0, sum = 0;
                for (int64_t i = 0; i < n; i++) { num += v[i] * Kx[i]; den += v[i] * Mx[i]; sum += Mx[i]; }
                lam[k] = den > 0 ? num / den : sigma + 1.0 / kr.re[k];
                double sc = den > 0 ? 1.0 / sqrt(den) : 1.0;
                if (sum < 0) sc = -sc;
                for (int64_t i = 0; i < n; i++) v[i] *= sc;
            }
            /* ascending order */
            int* ord = malloc(sizeof(int) * (size_t)nc);
            for (int k = 0; k < nc; k++) ord[k] = k;
            for (int a = 1; a < nc; a++) {
                int t = ord[a], b = a - 1;
                while (b >= 0 && lam[ord[b]] > lam[t]) { ord[b + 1] = ord[b]; b--; }
                ord[b + 1] = t;
            }
            Expr** vals = malloc(sizeof(Expr*) * (size_t)nev);
            for (int k = 0; k < nev; k++) vals[k] = expr_new_real(lam[ord[k]]);
            Expr* vlist = fem_list(vals, (size_t)nev);
            free(vals);
            if (!want_vectors) out = vlist;
            else {
                Expr* meshE = fem_mesh_expr(&S.mesh, F.region);
                Expr** fs = malloc(sizeof(Expr*) * (size_t)nev);
                double* U = malloc(sizeof(double) * (size_t)S.mesh.nn);
                for (int k = 0; k < nev; k++) {
                    const double* v = vec + (size_t)ord[k] * (size_t)n;
                    for (int64_t i = 0; i < S.mesh.nn; i++) U[i] = S.map[i] >= 0 ? v[S.map[i]] : 0.0;
                    Expr* ifun = fem_function_expr(meshE, &S.mesh, U);
                    if (F.applied) {
                        Expr* ar[2] = { expr_new_symbol(F.xv), expr_new_symbol(F.yv) };
                        ifun = expr_new_function(ifun, ar, 2);
                    }
                    fs[k] = ifun;
                }
                free(U);
                expr_free(meshE);
                Expr* pair[2] = { vlist, fem_list(fs, (size_t)nev) };
                free(fs);
                out = fem_list(pair, 2);
            }
            free(ord); free(lam); free(vec); free(Kx); free(Mx);
        }
        ks_result_free(&kr);
        csr_free(&K); csr_free(&M);
    }
    if (!ok && why) mth_message(head, "femfail", "%s.", why);
    fem_system_free(&S);
    fem_pde_free(&P);
    return out;
}

static Expr* builtin_ndeigenvalues(Expr* res) { return fem_eigen(res, false); }
static Expr* builtin_ndeigensystem(Expr* res) { return fem_eigen(res, true); }

/* ------------------------------------------------------------------ */
/*  Registration                                                       */
/* ------------------------------------------------------------------ */

void fem_init(void) {
    SYM_DirichletCondition = intern_symbol("DirichletCondition");
    SYM_NeumannValue = intern_symbol("NeumannValue");
    SYM_FEMMesh = intern_symbol("NDSolve`FEM`Mesh");
    SYM_FEMFunction = intern_symbol("NDSolve`FEM`MeshFunction");

    symtab_add_builtin("NDEigenvalues", builtin_ndeigenvalues);
    symtab_add_builtin("NDEigensystem", builtin_ndeigensystem);
    symtab_get_def("NDEigenvalues")->attributes |= ATTR_PROTECTED;
    symtab_get_def("NDEigensystem")->attributes |= ATTR_PROTECTED;
    symtab_get_def("DirichletCondition")->attributes |= ATTR_PROTECTED;
    symtab_get_def("NeumannValue")->attributes |= ATTR_PROTECTED;
    {
        Expr* o[2] = { nd_call2(SYM_Rule, expr_new_symbol(intern_symbol("DependentVariables")),
                                expr_new_symbol(SYM_Automatic)),
                       nd_call2(SYM_Rule, expr_new_symbol(SYM_Method), expr_new_symbol(SYM_Automatic)) };
        symtab_set_options("NDEigenvalues", fem_list(o, 2));
        Expr* o2[2] = { nd_call2(SYM_Rule, expr_new_symbol(intern_symbol("DependentVariables")),
                                 expr_new_symbol(SYM_Automatic)),
                        nd_call2(SYM_Rule, expr_new_symbol(SYM_Method), expr_new_symbol(SYM_Automatic)) };
        symtab_set_options("NDEigensystem", fem_list(o2, 2));
    }
    symtab_set_docstring("NDEigenvalues",
        "NDEigenvalues[{op, DirichletCondition[u[x, y] == 0, pred], ...}, u[x, y], {x, y} ∈ reg, n]\n"
        "\tgives the n smallest eigenvalues of the linear operator op (e.g. -Laplacian[u[x, y], {x, y}])\n"
        "\ton the region reg (Rectangle, Polygon, Triangle, Disk or a sector, Annulus), by quadratic\n"
        "\tfinite elements. Edges without a DirichletCondition have the natural (zero Neumann)\n"
        "\tcondition. Method -> {\"PDEDiscretization\" -> {\"FiniteElement\", {\"MeshOptions\" ->\n"
        "\t{\"MaxCellMeasure\" -> a}}}} sets the largest triangle area.");
    symtab_set_docstring("NDEigensystem",
        "NDEigensystem[{op, DirichletCondition[...]}, u[x, y], {x, y} ∈ reg, n] gives\n"
        "\t{{lambda1, ..., lambdan}, {f1, ..., fn}}: the n smallest eigenvalues of op on reg and their\n"
        "\teigenfunctions (InterpolatingFunction[...][x, y], normalized to unit L2 norm), by\n"
        "\tquadratic finite elements; options as for NDEigenvalues.");
    symtab_set_docstring("DirichletCondition",
        "DirichletCondition[u[x, y] == h, pred] represents the boundary condition u = h on the parts\n"
        "\tof the boundary where pred is True, in NDSolve, NDEigenvalues and NDEigensystem over a region.");
    symtab_set_docstring("NeumannValue",
        "NeumannValue[g, pred] represents the flux n . (c grad u) = g on the parts of the boundary\n"
        "\twhere pred is True: op u == f + NeumannValue[g, pred] in NDSolve over a region. g may be\n"
        "\tlinear in u (a Robin condition). Boundaries with neither condition have zero flux.");
}
