/* vectoranal.c -- Vector-analysis differential operators.
 *
 *   Grad[f, {x1,...,xn}]            gradient / Jacobian (Cartesian)
 *   Grad[f, {x1,...,xn}, chart]     gradient in an orthogonal coordinate chart
 *   Div[f, {x1,...,xn}[, chart]]    divergence (contracts the innermost slot)
 *   Curl[f, {x1,...,xn}[, chart]]   curl (generalized Levi-Civita form)
 *   Laplacian[f, {x1,...,xn}[, chart]]  Laplacian (Laplace-Beltrami)
 *
 * DESIGN.  Nothing here re-implements differentiation: every operator is
 * assembled as a Mathilda expression built out of `D[...]` and reduced with a
 * single evaluate().  The interpreter's array-derivative already does the heavy
 * lifting -- in particular `D[f, {{x1,...,xn}}]` is the gradient/Jacobian (it
 * appends a new innermost tensor slot), so Grad is a direct passthrough.
 *
 * Cartesian (2-arg) forms are fully general in the rank of f:
 *   Grad       D[f, {{vars}}]
 *   Laplacian  Sum_i D[f, {x_i, 2}]           (D threads over an array f)
 *   Div        contract the innermost slot of f with vars:
 *                vector  -> Sum_i D[f_i, x_i]
 *                tensor  -> map Div over the outer structure
 *   Curl       (1/k!) Sum eps_{a.. i j..} d_{x_i} f_{j..}, where k = depth(f);
 *              result depth n-k-1.  Covers 2-D vector -> scalar, 3-D vector ->
 *              vector, and rank-2 tensor -> scalar.
 *
 * Curvilinear (3-arg) forms use the orthogonal scale factors h_i (Lame
 * coefficients) of the chart, in the orthonormal (physical) basis:
 *   Grad[scalar]     { (1/h_i) D[f, x_i] }
 *   Div[vector]      (1/J) Sum_i D[(J/h_i) f_i, x_i],   J = Prod h_i
 *   Laplacian[scalar](1/J) Sum_i D[(J/h_i^2) D[f, x_i], x_i]
 *   Curl[vector]     3-D: (1/(h_j h_k))[D[h_k f_k, x_j] - D[h_j f_j, x_k]]
 *                    2-D: scalar (1/(h_1 h_2))[D[h_2 f_2, x_1] - D[h_1 f_1, x_2]]
 * Supported charts: "Cartesian", "Polar" (2-D), "Cylindrical" (3-D),
 * "Spherical" (3-D).  Any 3-arg call whose field rank does not match a clean
 * orthonormal formula (tensor fields, the vector Laplacian, Grad of a vector in
 * a chart -- all of which need Christoffel symbols) is left unevaluated, as is
 * an unrecognized chart.  This is Mathilda's honest "can't evaluate" contract:
 * a builtin returns NULL and the expression flows on unchanged.
 */

#include "vectoranal.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "sym_names.h"
#include "pack.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------- */
/* Tiny expression builders (per-module convention; cf. deriv.c)           */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }

static Expr* mk_fn1(const char* name, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(mk_sym(name), args, 1);
}

static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(mk_sym(name), args, 2);
}

/* Take ownership of `items` and wrap as name[items...]; frees the array. */
static Expr* mk_fnN_adopt(const char* name, Expr** items, size_t n) {
    Expr* r = expr_new_function(mk_sym(name), items, n);
    free(items);
    return r;
}

/* Times[-1, a]. */
static Expr* mk_neg(Expr* a) { return mk_fn2("Times", mk_int(-1), a); }

/* D[f, x] (both operands adopted). */
static Expr* mk_d(Expr* f, Expr* x) { return mk_fn2("D", f, x); }

/* ---------------------------------------------------------------------- */
/* Predicates and shape helpers                                            */
/* ---------------------------------------------------------------------- */

static bool is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* True iff `f` is a flat List of exactly `n` non-list elements (a vector). */
static bool is_vector_len(const Expr* f, int n) {
    if (!is_list(f) || (int)f->data.function.arg_count != n) return false;
    for (size_t i = 0; i < f->data.function.arg_count; i++)
        if (is_list(f->data.function.args[i])) return false;
    return true;
}

/* Depth of a uniform n-cube: 0 for a scalar (non-list) leaf, else the number
 * of nested List levels provided every level has exactly `n` elements and all
 * branches share the same depth.  Returns -1 for a ragged shape or a level
 * whose length is not n. */
static int ncube_depth(const Expr* f, int n) {
    if (!is_list(f)) return 0;
    if ((int)f->data.function.arg_count != n) return -1;
    int d0 = ncube_depth(f->data.function.args[0], n);
    if (d0 < 0) return -1;
    for (size_t i = 1; i < f->data.function.arg_count; i++)
        if (ncube_depth(f->data.function.args[i], n) != d0) return -1;
    return d0 + 1;
}

/* Levi-Civita sign of a permutation p[0..n-1] of {1..n} = (-1)^inversions. */
static int perm_sign(const int* p, int n) {
    int inv = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (p[i] > p[j]) inv++;
    return (inv % 2 == 0) ? 1 : -1;
}

static int64_t ifact(int k) {
    int64_t r = 1;
    for (int i = 2; i <= k; i++) r *= i;
    return r;
}

/* Return an owned copy of `a`, first materialising a packed NDArray to a plain
 * nested List (a small numeric constant field can arrive packed). */
static Expr* normalized_copy(Expr* a) {
    if (a->type == EXPR_NDARRAY) {
        Expr* u = pack_unpack(a);
        if (u) return u;
    }
    return expr_copy(a);
}

/* ---------------------------------------------------------------------- */
/* Cartesian Grad / Laplacian / Div                                        */
/* ---------------------------------------------------------------------- */

/* Grad[f, {{vars}}]: build the double-list array-derivative spec and evaluate.
 * Handles scalar -> vector, vector -> Jacobian, tensor -> +1 rank uniformly. */
static Expr* grad_cartesian(Expr* f, Expr** vars, int n) {
    Expr** vcopy = malloc((size_t)n * sizeof(Expr*));
    for (int i = 0; i < n; i++) vcopy[i] = expr_copy(vars[i]);
    Expr* varlist = mk_fnN_adopt("List", vcopy, (size_t)n);   /* {v1,...,vn} */
    Expr* spec    = mk_fn1("List", varlist);                  /* {{v1,...,vn}} */
    Expr* dcall   = mk_d(expr_copy(f), spec);                 /* D[f, {{...}}] */
    return eval_and_free(dcall);
}

/* Laplacian[f, vars] = Sum_i D[f, {x_i, 2}].  D threads over an explicit array
 * f, so the result carries f's own dimensions (element-wise Laplacian). */
static Expr* laplacian_cartesian(Expr* f, Expr** vars, int n) {
    if (n < 1) return NULL;
    Expr** terms = malloc((size_t)n * sizeof(Expr*));
    for (int i = 0; i < n; i++) {
        Expr* spec = mk_fn2("List", expr_copy(vars[i]), mk_int(2)); /* {x_i, 2} */
        terms[i] = mk_d(expr_copy(f), spec);
    }
    Expr* sum = (n == 1) ? terms[0] : mk_fnN_adopt("Plus", terms, (size_t)n);
    if (n == 1) free(terms);
    return eval_and_free(sum);
}

/* Build the (unevaluated) Cartesian divergence expression, contracting the
 * innermost slot of `f` with `vars`.  Returns NULL for a scalar or a
 * shape/length mismatch. */
static Expr* build_div(const Expr* f, Expr** vars, int n) {
    if (!is_list(f)) return NULL;                 /* scalar has no divergence */
    size_t m = f->data.function.arg_count;

    /* Tensor (list of lists): map Div over the outer structure. */
    if (m > 0 && is_list(f->data.function.args[0])) {
        Expr** kids = malloc(m * sizeof(Expr*));
        for (size_t i = 0; i < m; i++) {
            Expr* d = is_list(f->data.function.args[i])
                        ? build_div(f->data.function.args[i], vars, n) : NULL;
            if (!d) {
                for (size_t t = 0; t < i; t++) expr_free(kids[t]);
                free(kids);
                return NULL;
            }
            kids[i] = d;
        }
        return mk_fnN_adopt("List", kids, m);
    }

    /* Vector: Sum_i D[f_i, x_i]; require length n and no nested lists. */
    if ((int)m != n) return NULL;
    for (size_t i = 0; i < m; i++)
        if (is_list(f->data.function.args[i])) return NULL;
    Expr** terms = malloc(m * sizeof(Expr*));
    for (size_t i = 0; i < m; i++)
        terms[i] = mk_d(expr_copy(f->data.function.args[i]), expr_copy(vars[i]));
    return mk_fnN_adopt("Plus", terms, m);
}

static Expr* div_cartesian(Expr* f, Expr** vars, int n) {
    Expr* pre = build_div(f, vars, n);
    return pre ? eval_and_free(pre) : NULL;
}

/* ---------------------------------------------------------------------- */
/* Cartesian Curl -- generalized Levi-Civita contraction                   */
/* ---------------------------------------------------------------------- */

typedef struct {
    Expr* f;          /* the uniform n-cube field (owned by the caller) */
    Expr** vars;      /* n coordinate expressions (borrowed) */
    int n, k, dep;    /* dim, depth of f, result depth = n - k - 1 */
    int64_t kf;       /* k! */
    size_t ncells;    /* n^dep */
    Expr*** cterms;   /* [ncells] growable arrays of term expressions */
    size_t* ccount;
    size_t* ccap;
} curl_ctx;

static void curl_cell_append(curl_ctx* c, size_t idx, Expr* term) {
    if (c->ccount[idx] == c->ccap[idx]) {
        size_t nc = c->ccap[idx] ? c->ccap[idx] * 2 : 4;
        c->cterms[idx] = realloc(c->cterms[idx], nc * sizeof(Expr*));
        c->ccap[idx] = nc;
    }
    c->cterms[idx][c->ccount[idx]++] = term;
}

/* For one signed permutation split as (a_part | i | j_part), accumulate
 * sign * D[f_{j_part}, x_i] into the cell indexed by a_part. */
static void curl_emit_perm(curl_ctx* c, const int* perm) {
    int s = perm_sign(perm, c->n);
    int i = perm[c->dep];
    const Expr* leaf = c->f;
    for (int t = 0; t < c->k; t++)
        leaf = leaf->data.function.args[perm[c->dep + 1 + t] - 1];
    Expr* d = mk_d(expr_copy((Expr*)leaf), expr_copy(c->vars[i - 1]));
    Expr* term = (s == 1) ? d : mk_neg(d);
    size_t idx = 0;
    for (int t = 0; t < c->dep; t++)
        idx = idx * (size_t)c->n + (size_t)(perm[t] - 1);
    curl_cell_append(c, idx, term);
}

static void curl_perm_recur(curl_ctx* c, int* perm, bool* used, int depth) {
    if (depth == c->n) { curl_emit_perm(c, perm); return; }
    for (int v = 1; v <= c->n; v++) {
        if (used[v - 1]) continue;
        used[v - 1] = true;
        perm[depth] = v;
        curl_perm_recur(c, perm, used, depth + 1);
        used[v - 1] = false;
    }
}

/* Collapse one cell's accumulated terms into Plus[...] / k! (adopts terms). */
static Expr* curl_cell_expr(curl_ctx* c, size_t idx) {
    size_t m = c->ccount[idx];
    Expr* sum;
    if (m == 0) {
        sum = mk_int(0);
    } else if (m == 1) {
        sum = c->cterms[idx][0];
        c->cterms[idx][0] = NULL;
    } else {
        Expr** items = malloc(m * sizeof(Expr*));
        for (size_t t = 0; t < m; t++) {
            items[t] = c->cterms[idx][t];
            c->cterms[idx][t] = NULL;
        }
        sum = mk_fnN_adopt("Plus", items, m);
    }
    if (c->k > 1) sum = mk_fn2("Divide", sum, mk_int(c->kf));
    return sum;
}

/* Assemble the depth-`dep` nested-list result (row-major over the cells). */
static Expr* curl_build_nested(curl_ctx* c, int level, size_t base) {
    if (level == c->dep) return curl_cell_expr(c, base);
    Expr** kids = malloc((size_t)c->n * sizeof(Expr*));
    for (int t = 0; t < c->n; t++)
        kids[t] = curl_build_nested(c, level + 1, base * (size_t)c->n + (size_t)t);
    return mk_fnN_adopt("List", kids, (size_t)c->n);
}

static Expr* curl_cartesian(Expr* f, Expr** vars, int n) {
    if (n < 2 || n > 6) return NULL;            /* bound the n! enumeration */
    int k = ncube_depth(f, n);
    if (k < 1) return NULL;                     /* need a vector/tensor field */
    int dep = n - k - 1;
    if (dep < 0) return NULL;                   /* k > n-1: undefined */

    size_t ncells = 1;
    for (int t = 0; t < dep; t++) ncells *= (size_t)n;

    curl_ctx c;
    c.f = f; c.vars = vars; c.n = n; c.k = k; c.dep = dep; c.kf = ifact(k);
    c.ncells = ncells;
    c.cterms = calloc(ncells, sizeof(Expr**));
    c.ccount = calloc(ncells, sizeof(size_t));
    c.ccap   = calloc(ncells, sizeof(size_t));

    int* perm = malloc((size_t)n * sizeof(int));
    bool* used = calloc((size_t)n, sizeof(bool));
    curl_perm_recur(&c, perm, used, 0);
    free(perm);
    free(used);

    Expr* pre = curl_build_nested(&c, 0, 0);    /* consumes every term */

    for (size_t i = 0; i < ncells; i++) free(c.cterms[i]);
    free(c.cterms);
    free(c.ccount);
    free(c.ccap);

    return eval_and_free(pre);
}

/* ---------------------------------------------------------------------- */
/* Coordinate charts (orthogonal scale factors / Lame coefficients)        */
/* ---------------------------------------------------------------------- */

/* Fill h[0..n-1] with freshly-owned scale-factor expressions for `chart` in
 * the positional coordinates `vars`.  Returns true on success; sets *known to
 * whether the chart name is recognized at all (false + known==true means the
 * name is known but its dimension does not match n). */
static bool chart_scale_factors(const char* chart, Expr** vars, int n,
                                Expr** h, bool* known) {
    *known = true;
    if (strcmp(chart, "Cartesian") == 0) {
        for (int i = 0; i < n; i++) h[i] = mk_int(1);
        return true;
    }
    if (strcmp(chart, "Polar") == 0) {
        if (n != 2) return false;
        h[0] = mk_int(1);
        h[1] = expr_copy(vars[0]);                     /* r */
        return true;
    }
    if (strcmp(chart, "Cylindrical") == 0) {
        if (n != 3) return false;
        h[0] = mk_int(1);
        h[1] = expr_copy(vars[0]);                     /* r */
        h[2] = mk_int(1);
        return true;
    }
    if (strcmp(chart, "Spherical") == 0) {
        if (n != 3) return false;
        h[0] = mk_int(1);
        h[1] = expr_copy(vars[0]);                     /* r */
        h[2] = mk_fn2("Times", expr_copy(vars[0]),
                      mk_fn1("Sin", expr_copy(vars[1]))); /* r Sin[theta] */
        return true;
    }
    *known = false;
    return false;
}

/* J = Prod_i h_i, built fresh from copies of h. */
static Expr* mk_jacobian(Expr** h, int n) {
    if (n == 1) return expr_copy(h[0]);
    Expr** items = malloc((size_t)n * sizeof(Expr*));
    for (int i = 0; i < n; i++) items[i] = expr_copy(h[i]);
    return mk_fnN_adopt("Times", items, (size_t)n);
}

/* Resolve the chart string and its scale factors; emits Head::chart for an
 * unrecognized name.  On success the caller owns h[0..n-1] and *h_out. */
static bool resolve_chart(Expr* chart, Expr** vars, int n, const char* head,
                          Expr*** h_out) {
    if (chart->type != EXPR_STRING) return false;    /* nested/metric spec */
    Expr** h = malloc((size_t)n * sizeof(Expr*));
    bool known = false;
    if (!chart_scale_factors(chart->data.string, vars, n, h, &known)) {
        free(h);
        if (!known)
            fprintf(stderr,
                    "%s::chart: \"%s\" is not a supported coordinate chart.\n",
                    head, chart->data.string);
        return false;
    }
    *h_out = h;
    return true;
}

static void free_h(Expr** h, int n) {
    for (int i = 0; i < n; i++) expr_free(h[i]);
    free(h);
}

/* Grad[scalar f, vars, chart] = { (1/h_i) D[f, x_i] }. */
static Expr* grad_chart(Expr* f, Expr** vars, int n, Expr* chart,
                        const char* head) {
    if (is_list(f)) return NULL;                 /* vector grad needs a metric */
    Expr** h;
    if (!resolve_chart(chart, vars, n, head, &h)) return NULL;
    Expr** comps = malloc((size_t)n * sizeof(Expr*));
    for (int i = 0; i < n; i++)
        comps[i] = mk_fn2("Divide", mk_d(expr_copy(f), expr_copy(vars[i])),
                          expr_copy(h[i]));
    Expr* pre = mk_fnN_adopt("List", comps, (size_t)n);
    free_h(h, n);
    return eval_and_free(pre);
}

/* Div[vector f, vars, chart] = (1/J) Sum_i D[(J/h_i) f_i, x_i]. */
static Expr* div_chart(Expr* f, Expr** vars, int n, Expr* chart,
                       const char* head) {
    if (!is_vector_len(f, n)) return NULL;       /* tensor div needs a metric */
    Expr** h;
    if (!resolve_chart(chart, vars, n, head, &h)) return NULL;
    Expr** terms = malloc((size_t)n * sizeof(Expr*));
    for (int i = 0; i < n; i++) {
        Expr* coeff = mk_fn2("Divide", mk_jacobian(h, n), expr_copy(h[i]));
        Expr* inner = mk_fn2("Times", coeff,
                             expr_copy(f->data.function.args[i]));
        terms[i] = mk_d(inner, expr_copy(vars[i]));
    }
    Expr* sum = (n == 1) ? terms[0] : mk_fnN_adopt("Plus", terms, (size_t)n);
    if (n == 1) free(terms);
    Expr* pre = mk_fn2("Divide", sum, mk_jacobian(h, n));
    free_h(h, n);
    return eval_and_free(pre);
}

/* Laplacian[scalar f, vars, chart] = (1/J) Sum_i D[(J/h_i^2) D[f, x_i], x_i]. */
static Expr* laplacian_chart(Expr* f, Expr** vars, int n, Expr* chart,
                             const char* head) {
    if (is_list(f)) return NULL;                 /* vector Laplacian: Christoffel */
    Expr** h;
    if (!resolve_chart(chart, vars, n, head, &h)) return NULL;
    Expr** terms = malloc((size_t)n * sizeof(Expr*));
    for (int i = 0; i < n; i++) {
        Expr* h2 = mk_fn2("Power", expr_copy(h[i]), mk_int(2));
        Expr* coeff = mk_fn2("Divide", mk_jacobian(h, n), h2);
        Expr* inner = mk_fn2("Times", coeff,
                             mk_d(expr_copy(f), expr_copy(vars[i])));
        terms[i] = mk_d(inner, expr_copy(vars[i]));
    }
    Expr* sum = (n == 1) ? terms[0] : mk_fnN_adopt("Plus", terms, (size_t)n);
    if (n == 1) free(terms);
    Expr* pre = mk_fn2("Divide", sum, mk_jacobian(h, n));
    free_h(h, n);
    return eval_and_free(pre);
}

/* One orthonormal curl component (1/(h_j h_k))[D[h_k f_k, x_j] - D[h_j f_j, x_k]]. */
static Expr* curl_component(Expr* f, Expr** vars, Expr** h, int j, int k) {
    Expr* fj = f->data.function.args[j];
    Expr* fk = f->data.function.args[k];
    Expr* a = mk_d(mk_fn2("Times", expr_copy(h[k]), expr_copy(fk)),
                   expr_copy(vars[j]));
    Expr* b = mk_d(mk_fn2("Times", expr_copy(h[j]), expr_copy(fj)),
                   expr_copy(vars[k]));
    Expr* num = mk_fn2("Plus", a, mk_neg(b));
    Expr* den = mk_fn2("Times", expr_copy(h[j]), expr_copy(h[k]));
    return mk_fn2("Divide", num, den);
}

/* Curl[vector f, vars, chart]: 3-D vector or 2-D scalar in the orthonormal basis. */
static Expr* curl_chart(Expr* f, Expr** vars, int n, Expr* chart,
                        const char* head) {
    if (n != 2 && n != 3) return NULL;
    if (!is_vector_len(f, n)) return NULL;       /* tensor curl needs a metric */
    Expr** h;
    if (!resolve_chart(chart, vars, n, head, &h)) return NULL;

    Expr* pre;
    if (n == 2) {
        pre = curl_component(f, vars, h, 0, 1);  /* scalar */
    } else {
        Expr** comps = malloc(3 * sizeof(Expr*));
        for (int i = 0; i < 3; i++)              /* cyclic (i, i+1, i+2) */
            comps[i] = curl_component(f, vars, h, (i + 1) % 3, (i + 2) % 3);
        pre = mk_fnN_adopt("List", comps, 3);
    }
    free_h(h, n);
    return eval_and_free(pre);
}

/* ---------------------------------------------------------------------- */
/* Builtin dispatch                                                        */
/* ---------------------------------------------------------------------- */

/* Common front end: validate arity, normalize the coordinate list, and hand
 * off to the Cartesian (argc==2) or chart (argc==3) implementation.  The two
 * function pointers build a fully-evaluated result or return NULL. */
typedef Expr* (*cart_fn)(Expr* f, Expr** vars, int n);
typedef Expr* (*chart_fn)(Expr* f, Expr** vars, int n, Expr* chart,
                          const char* head);

static Expr* vecop(Expr* res, const char* head, cart_fn cart, chart_fn chart) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;

    Expr* vlist = normalized_copy(res->data.function.args[1]);
    if (!is_list(vlist)) { expr_free(vlist); return NULL; }
    int n = (int)vlist->data.function.arg_count;
    Expr** vars = vlist->data.function.args;

    Expr* f = normalized_copy(res->data.function.args[0]);
    Expr* result = (argc == 2)
        ? cart(f, vars, n)
        : chart(f, vars, n, res->data.function.args[2], head);

    expr_free(f);
    expr_free(vlist);
    return result;
}

Expr* builtin_grad(Expr* res) {
    return vecop(res, "Grad", grad_cartesian, grad_chart);
}

Expr* builtin_div(Expr* res) {
    return vecop(res, "Div", div_cartesian, div_chart);
}

Expr* builtin_curl(Expr* res) {
    return vecop(res, "Curl", curl_cartesian, curl_chart);
}

Expr* builtin_laplacian(Expr* res) {
    return vecop(res, "Laplacian", laplacian_cartesian, laplacian_chart);
}

void vectoranal_init(void) {
    symtab_add_builtin("Grad", builtin_grad);
    symtab_add_builtin("Div", builtin_div);
    symtab_add_builtin("Curl", builtin_curl);
    symtab_add_builtin("Laplacian", builtin_laplacian);
    symtab_get_def("Grad")->attributes      |= ATTR_PROTECTED;
    symtab_get_def("Div")->attributes        |= ATTR_PROTECTED;
    symtab_get_def("Curl")->attributes       |= ATTR_PROTECTED;
    symtab_get_def("Laplacian")->attributes  |= ATTR_PROTECTED;
}
