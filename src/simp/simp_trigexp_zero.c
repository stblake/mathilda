/*
 * simp_trigexp_zero.c — see simp_trigexp_zero.h.
 *
 * ALGORITHM (exact, fast, no numeric sampling)
 * --------------------------------------------
 * An expression that is a rational function of a single trig/exp kernel
 * t = E^(i x) (times possibly some opaque transcendental subterms and free
 * parameters, which enter as independent generators) is decided to be
 * identically 0 by EXACT polynomial identity testing: evaluate it at a tensor
 * grid of distinct RATIONAL points whose size exceeds the numerator degree in
 * every generator. A polynomial of degree ≤ d_v in generator v that vanishes
 * on d_v+1 distinct values per axis is identically zero (combinatorial
 * Nullstellensatz) — so a pole-free grid of all-zero evaluations is a PROOF.
 *
 * Generators (each an independent grid axis):
 *   • the kernel t = E^(i x): every Sin[k x], Cos[k x], Sec[k x], …, E^(k i x)
 *     is replaced by its EXACT value at t = t0 (a rational):
 *        Cos[k x] = (t0^k + t0^-k)/2       E^(k i x) = t0^k
 *        Sin[k x] = (t0^k - t0^-k)/(2 i)   Sec/Csc/Tan/Cot accordingly
 *   • opaque subterms that depend on x but are not kernels (Log[...], Sqrt[...],
 *     …) — a Risch diff-back D[G]-f keeps Log terms from the product rule that
 *     cancel only when treated as independent transcendentals, exactly as
 *     Together treats them. Each is substituted by a fresh independent rational.
 *   • free parameters (symbols other than x) — coefficients like a, b.
 *
 * Evaluation is done on the ORIGINAL compact expression (a diff-back is only a
 * few hundred leaves), collapsing each point to an exact (Gaussian) rational in
 * sub-millisecond time — we never build the exponentially-larger TrigToExp form
 * nor call the generic (seconds-slow) Together. Points are rational, so there is
 * no floating-point sampling and no branch-cut hazard: an algebraic decision.
 *
 * Soundness of declining: treating opaque subterms as independent can only make
 * the test MORE conservative (an identity holding via a transcendental relation
 * among the opaques may read as non-zero → we decline, never a false zero).
 *
 * Reused by the Simplify fast path `transform_trigexp_vanish`
 * (src/simp/simp_builtins.c) and by `decide_trigexp_rational` (src/zero_test.c).
 */
#include "simp_trigexp_zero.h"
#include "simp_internal.h"

#include "eval.h"
#include "expr.h"
#include "sym_names.h"

#include <stdlib.h>
#include <string.h>

/* Guardrails: decline rather than grind on pathological sizes. */
#define TEZ_MAX_AXIS_DEGREE 2000  /* per-generator numerator-degree cap        */
#define TEZ_MAX_GRID_POINTS 1500  /* total tensor-grid evaluation cap (~1-2 s) */
#define TEZ_MAX_AXES        8     /* kernel + opaque generators + params       */
#define TEZ_POLE_RETRIES    8     /* grid-offset retries when a pole is hit    */

/* ------------------------------------------------------------------ */
/*  small constructors / predicates                                   */
/* ------------------------------------------------------------------ */

static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }
static Expr* mk_int(int64_t n)     { return expr_new_integer(n); }

static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(mk_sym(h), args, 2);
}

static bool sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL &&
           (e->data.symbol.name == name ||
            strcmp(e->data.symbol.name, name) == 0);
}

static bool is_known_const(const Expr* e) {
    return sym_is(e, SYM_Pi) || sym_is(e, SYM_E) || sym_is(e, SYM_I) ||
           sym_is(e, "EulerGamma") || sym_is(e, "Degree") ||
           sym_is(e, "GoldenRatio") || sym_is(e, "Catalan") ||
           sym_is(e, "Infinity") || sym_is(e, "Indeterminate") ||
           sym_is(e, "ComplexInfinity");
}

static bool contains_named_symbol(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return sym_is(e, name);
    if (e->type == EXPR_FUNCTION) {
        if (contains_named_symbol(e->data.function.head, name)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (contains_named_symbol(e->data.function.args[i], name)) return true;
    }
    return false;
}

/* ------------------------------------------------------------------ */
/*  kernel classification                                             */
/* ------------------------------------------------------------------ */

typedef enum { K_SIN, K_COS, K_TAN, K_COT, K_SEC, K_CSC, K_EXP } KernelKind;

/* Evaluate q = num / (mul * var); Integer when q is a genuine integer. */
static Expr* ratio_over_var(const Expr* num, const Expr* var, const char* mul) {
    Expr* denom = mul ? mk_fn2(SYM_Times, mk_sym(mul), expr_copy((Expr*)var))
                      : expr_copy((Expr*)var);
    Expr* prod = mk_fn2(SYM_Times, expr_copy((Expr*)num),
                        mk_fn2(SYM_Power, denom, mk_int(-1)));
    return eval_and_free(prod);
}

/* If `e` is a supported single trig/exp kernel of `var`, set *kind and the
 * integer multiple *m and return true. */
static bool classify_kernel(const Expr* e, const Expr* var,
                            KernelKind* kind, long* m) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    size_t n = e->data.function.arg_count;
    if (!h || h->type != EXPR_SYMBOL) return false;
    const char* hn = h->data.symbol.name;

    const Expr* exponent = NULL;
    if (sym_is(h, SYM_Exp) && n == 1) exponent = e->data.function.args[0];
    else if (sym_is(h, SYM_Power) && n == 2 &&
             sym_is(e->data.function.args[0], SYM_E))
        exponent = e->data.function.args[1];
    if (exponent) {
        Expr* q = ratio_over_var(exponent, var, SYM_I);
        bool okint = q && q->type == EXPR_INTEGER;
        if (okint) { *kind = K_EXP; *m = q->data.integer; }
        if (q) expr_free(q);
        return okint;
    }

    KernelKind k;
    if (sym_is(h, SYM_Sin))          k = K_SIN;
    else if (sym_is(h, SYM_Cos))     k = K_COS;
    else if (strcmp(hn, "Tan") == 0) k = K_TAN;
    else if (strcmp(hn, "Cot") == 0) k = K_COT;
    else if (strcmp(hn, "Sec") == 0) k = K_SEC;
    else if (strcmp(hn, "Csc") == 0) k = K_CSC;
    else return false;
    if (n != 1) return false;
    Expr* q = ratio_over_var(e->data.function.args[0], var, NULL);
    bool okint = q && q->type == EXPR_INTEGER && q->data.integer != 0;
    if (okint) { *kind = k; *m = q->data.integer; }
    if (q) expr_free(q);
    return okint;
}

/* ------------------------------------------------------------------ */
/*  generator collection                                              */
/* ------------------------------------------------------------------ */

typedef struct { Expr* node; KernelKind kind; long m; } Kernel;
typedef struct { Kernel* items; size_t count, cap; } KernelSet;

/* Non-kernel independent generators: opaque x-dependent subterms and, appended
 * later, the free-parameter symbols. Each becomes one grid axis. */
typedef struct { Expr** items; size_t count, cap; } AtomSet;

static void kset_add(KernelSet* ks, const Expr* e, KernelKind kind, long m) {
    for (size_t i = 0; i < ks->count; i++)
        if (expr_eq(ks->items[i].node, e)) return;
    if (ks->count == ks->cap) {
        ks->cap = ks->cap ? ks->cap * 2 : 8;
        ks->items = realloc(ks->items, ks->cap * sizeof(Kernel));
    }
    ks->items[ks->count].node = expr_copy((Expr*)e);
    ks->items[ks->count].kind = kind;
    ks->items[ks->count].m    = m;
    ks->count++;
}
static void kset_free(KernelSet* ks) {
    for (size_t i = 0; i < ks->count; i++) expr_free(ks->items[i].node);
    free(ks->items);
}

static void aset_add(AtomSet* as, Expr* owned) {
    for (size_t i = 0; i < as->count; i++)
        if (expr_eq(as->items[i], owned)) { expr_free(owned); return; }
    if (as->count == as->cap) {
        as->cap = as->cap ? as->cap * 2 : 8;
        as->items = realloc(as->items, as->cap * sizeof(Expr*));
    }
    as->items[as->count++] = owned;
}
static void aset_free(AtomSet* as) {
    for (size_t i = 0; i < as->count; i++) expr_free(as->items[i]);
    free(as->items);
}

/* Walk `e`, recording trig/exp kernels and opaque x-dependent generators.
 * *supported becomes false on a bare `var` (non-kernel polynomial dependence)
 * or any structure that cannot be treated as a rational function of the
 * generators. `vn` is the kernel-variable name. */
static void collect(const Expr* e, const Expr* var, const char* vn,
                    KernelSet* ks, AtomSet* opaque, bool* supported) {
    if (!e || !*supported) return;

    KernelKind kind; long m;
    if (classify_kernel(e, var, &kind, &m)) { kset_add(ks, e, kind, m); return; }

    if (e->type == EXPR_SYMBOL) {
        if (sym_is(e, vn)) *supported = false;   /* bare var, not in a kernel */
        return;                                  /* param/const → coefficient */
    }
    if (e->type != EXPR_FUNCTION) return;        /* number */

    const Expr* h = e->data.function.head;
    size_t n = e->data.function.arg_count;

    if (sym_is(h, SYM_Plus) || sym_is(h, SYM_Times)) {
        for (size_t i = 0; i < n; i++)
            collect(e->data.function.args[i], var, vn, ks, opaque, supported);
        return;
    }
    if (sym_is(h, SYM_Power) && n == 2 &&
        e->data.function.args[1]->type == EXPR_INTEGER) {
        collect(e->data.function.args[0], var, vn, ks, opaque, supported);
        return;
    }
    /* Any other head. If it involves the variable it is an opaque independent
     * generator (Log[...], Sqrt[...], non-integer power, unknown f[...]); if
     * not, it is a constant coefficient. */
    if (contains_named_symbol(e, vn)) aset_add(opaque, expr_copy((Expr*)e));
}

/* ------------------------------------------------------------------ */
/*  degree bounds (numerator degree in one generator axis)            */
/* ------------------------------------------------------------------ */

typedef struct { long n, d; bool ok; } Deg;

/* Numerator/denominator degree bounds of `e` in one axis.
 * axis == -1 → the kernel t; kernels contribute per head/|m|, atoms are t-free.
 * axis >= 0  → atoms[axis]; that atom has degree 1, kernels and other atoms 0. */
static Deg poly_deg(const Expr* e, const Expr* var,
                    Expr** atoms, size_t natoms, long axis) {
    Deg r = { 0, 0, true };
    if (!e) return r;

    KernelKind kind; long m;
    if (classify_kernel(e, var, &kind, &m)) {
        if (axis != -1) return r;                /* kernel is t-only */
        long a = m < 0 ? -m : m;
        switch (kind) {
            case K_COS: case K_SIN: r.n = 2 * a; r.d = a;         return r;
            case K_SEC: case K_CSC: r.n = a;     r.d = 2 * a;     return r;
            case K_TAN: case K_COT: r.n = 3 * a; r.d = 3 * a;     return r;
            case K_EXP: r.n = m > 0 ? m : 0; r.d = m < 0 ? -m : 0; return r;
        }
    }
    for (size_t i = 0; i < natoms; i++)
        if (expr_eq(e, atoms[i])) { if (axis == (long)i) r.n = 1; return r; }

    if (e->type != EXPR_FUNCTION) return r;      /* number / const / param */
    const Expr* h = e->data.function.head;
    size_t n = e->data.function.arg_count;

    if (sym_is(h, SYM_Plus)) {
        long sumd = 0;
        Deg* cs = malloc(n * sizeof(Deg));
        for (size_t i = 0; i < n; i++) {
            cs[i] = poly_deg(e->data.function.args[i], var, atoms, natoms, axis);
            if (!cs[i].ok) { r.ok = false; free(cs); return r; }
            sumd += cs[i].d;
        }
        long best = 0;
        for (size_t i = 0; i < n; i++) {
            long cand = cs[i].n + (sumd - cs[i].d);
            if (cand > best) best = cand;
        }
        r.n = best; r.d = sumd; free(cs); return r;
    }
    if (sym_is(h, SYM_Times)) {
        for (size_t i = 0; i < n; i++) {
            Deg c = poly_deg(e->data.function.args[i], var, atoms, natoms, axis);
            if (!c.ok) { r.ok = false; return r; }
            r.n += c.n; r.d += c.d;
        }
        return r;
    }
    if (sym_is(h, SYM_Power) && n == 2 &&
        e->data.function.args[1]->type == EXPR_INTEGER) {
        long p = e->data.function.args[1]->data.integer;
        Deg base = poly_deg(e->data.function.args[0], var, atoms, natoms, axis);
        if (!base.ok) { r.ok = false; return r; }
        long ap = p < 0 ? -p : p;
        if (p >= 0) { r.n = ap * base.n; r.d = ap * base.d; }
        else        { r.n = ap * base.d; r.d = ap * base.n; }
        return r;
    }
    r.ok = false;                                /* unbounded head */
    return r;
}

/* ------------------------------------------------------------------ */
/*  exact kernel value at t = t0                                      */
/* ------------------------------------------------------------------ */

static Expr* pw(long t0, long m) { return mk_fn2(SYM_Power, mk_int(t0), mk_int(m)); }

static Expr* kernel_value(KernelKind kind, long m, long t0) {
    Expr* I = mk_sym(SYM_I);
    switch (kind) {
        case K_EXP:
            expr_free(I); return pw(t0, m);
        case K_COS:
            expr_free(I);
            return mk_fn2(SYM_Times, mk_fn2(SYM_Power, mk_int(2), mk_int(-1)),
                          mk_fn2(SYM_Plus, pw(t0, m), pw(t0, -m)));
        case K_SIN:
            return mk_fn2(SYM_Times,
                mk_fn2(SYM_Power, mk_fn2(SYM_Times, mk_int(2), I), mk_int(-1)),
                mk_fn2(SYM_Plus, pw(t0, m), mk_fn2(SYM_Times, mk_int(-1), pw(t0, -m))));
        case K_SEC:
            expr_free(I);
            return mk_fn2(SYM_Times, mk_int(2),
                mk_fn2(SYM_Power, mk_fn2(SYM_Plus, pw(t0, m), pw(t0, -m)), mk_int(-1)));
        case K_CSC:
            return mk_fn2(SYM_Times, mk_fn2(SYM_Times, mk_int(2), I),
                mk_fn2(SYM_Power,
                    mk_fn2(SYM_Plus, pw(t0, m), mk_fn2(SYM_Times, mk_int(-1), pw(t0, -m))),
                    mk_int(-1)));
        case K_TAN:
            return mk_fn2(SYM_Times, mk_fn2(SYM_Power, I, mk_int(-1)),
                mk_fn2(SYM_Times,
                    mk_fn2(SYM_Plus, pw(t0, m), mk_fn2(SYM_Times, mk_int(-1), pw(t0, -m))),
                    mk_fn2(SYM_Power, mk_fn2(SYM_Plus, pw(t0, m), pw(t0, -m)), mk_int(-1))));
        case K_COT:
            return mk_fn2(SYM_Times, I,
                mk_fn2(SYM_Times, mk_fn2(SYM_Plus, pw(t0, m), pw(t0, -m)),
                    mk_fn2(SYM_Power,
                        mk_fn2(SYM_Plus, pw(t0, m), mk_fn2(SYM_Times, mk_int(-1), pw(t0, -m))),
                        mk_int(-1))));
    }
    expr_free(I);
    return mk_int(0);
}

/* ------------------------------------------------------------------ */
/*  point evaluation                                                  */
/* ------------------------------------------------------------------ */

static bool num_is_zero(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: return e->data.integer == 0;
        case EXPR_REAL:    return e->data.real == 0.0;
        case EXPR_BIGINT:  return mpz_sgn(e->data.bigint) == 0;
        default: break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        sym_is(e->data.function.head, SYM_Rational))
        return num_is_zero(e->data.function.args[0]);
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        sym_is(e->data.function.head, SYM_Complex))
        return num_is_zero(e->data.function.args[0]) &&
               num_is_zero(e->data.function.args[1]);
    return false;
}

typedef enum { PT_ZERO, PT_NONZERO, PT_POLE, PT_UNDECIDED } PointResult;

static PointResult classify_point(const Expr* e) {
    if (!e) return PT_UNDECIDED;
    if (sym_is(e, "ComplexInfinity") || sym_is(e, "Indeterminate")) return PT_POLE;
    if (e->type == EXPR_FUNCTION &&
        (sym_is(e->data.function.head, "DirectedInfinity") ||
         sym_is(e->data.function.head, "Infinity")))
        return PT_POLE;
    if (num_is_zero(e)) return PT_ZERO;
    if (expr_is_numeric_like(e)) return PT_NONZERO;
    return PT_UNDECIDED;
}

/* Evaluate `df` with kernels → value at t0 and each atom → its axis value. */
static PointResult eval_point(const Expr* df, const KernelSet* ks, long t0,
                              Expr** atoms, const long* avals, size_t natoms) {
    size_t nrules = ks->count + natoms;
    Expr** rules = malloc(nrules * sizeof(Expr*));
    for (size_t i = 0; i < ks->count; i++)
        rules[i] = mk_fn2(SYM_Rule, expr_copy(ks->items[i].node),
                          kernel_value(ks->items[i].kind, ks->items[i].m, t0));
    for (size_t j = 0; j < natoms; j++)
        rules[ks->count + j] = mk_fn2(SYM_Rule, expr_copy(atoms[j]), mk_int(avals[j]));

    Expr* rulelist = expr_new_function(mk_sym(SYM_List), rules, nrules);
    free(rules);
    Expr* val = eval_and_free(mk_fn2(SYM_ReplaceAll, expr_copy((Expr*)df), rulelist));
    PointResult r = classify_point(val);
    if (val) expr_free(val);
    return r;
}

/* ------------------------------------------------------------------ */
/*  tensor-grid certification                                         */
/* ------------------------------------------------------------------ */

static bool grid_next(long* idx, const long* deg, size_t nv) {
    for (size_t v = 0; v < nv; v++) {
        if (idx[v] < deg[v]) { idx[v]++; return true; }
        idx[v] = 0;
    }
    return false;
}

/* Append distinct free-parameter symbols (not the kernel var, not constants)
 * as further Atom axes; aset_add dedups against opaque generators. */
static void gather_params(const Expr* e, const char* vn, AtomSet* atoms) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (!sym_is(e, vn) && !is_known_const(e))
            aset_add(atoms, expr_copy((Expr*)e));
        return;
    }
    if (e->type == EXPR_FUNCTION)
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            gather_params(e->data.function.args[i], vn, atoms);
}

static TrigExpZeroResult decide_with_var(const Expr* e, const Expr* var) {
    const char* vn = var->data.symbol.name;

    KernelSet ks = {0};
    AtomSet atoms = {0};                          /* opaque generators, then params */
    bool supported = true;
    collect(e, var, vn, &ks, &atoms, &supported);
    if (!supported || ks.count == 0) { kset_free(&ks); aset_free(&atoms); return TRIGEXP_ZERO_UNKNOWN; }

    /* Append free-parameter symbols (symbols other than the kernel variable)
     * as further independent axes; aset_add dedups against the opaques. */
    gather_params(e, vn, &atoms);

    size_t naxes = 1 + atoms.count;               /* axis 0 = t, then atoms */
    if (naxes > TEZ_MAX_AXES) { kset_free(&ks); aset_free(&atoms); return TRIGEXP_ZERO_UNKNOWN; }

    long* deg = malloc(naxes * sizeof(long));
    Deg dt = poly_deg(e, var, atoms.items, atoms.count, -1);
    if (!dt.ok || dt.n > TEZ_MAX_AXIS_DEGREE) {
        free(deg); kset_free(&ks); aset_free(&atoms); return TRIGEXP_ZERO_UNKNOWN;
    }
    deg[0] = dt.n;
    for (size_t j = 0; j < atoms.count; j++) {
        Deg da = poly_deg(e, var, atoms.items, atoms.count, (long)j);
        if (!da.ok || da.n > TEZ_MAX_AXIS_DEGREE) {
            free(deg); kset_free(&ks); aset_free(&atoms); return TRIGEXP_ZERO_UNKNOWN;
        }
        deg[j + 1] = da.n;
    }

    double grid = 1.0;
    for (size_t v = 0; v < naxes; v++) grid *= (double)(deg[v] + 1);
    if (grid > (double)TEZ_MAX_GRID_POINTS) {
        free(deg); kset_free(&ks); aset_free(&atoms); return TRIGEXP_ZERO_UNKNOWN;
    }

    TrigExpZeroResult verdict = TRIGEXP_ZERO_UNKNOWN;
    for (int attempt = 0; attempt < TEZ_POLE_RETRIES; attempt++) {
        long base = 2 + attempt;
        long* idx = calloc(naxes, sizeof(long));
        long* avals = atoms.count ? malloc(atoms.count * sizeof(long)) : NULL;
        bool pole = false, undecided = false, nonzero = false;
        do {
            long t0 = base + idx[0];
            for (size_t j = 0; j < atoms.count; j++)
                avals[j] = base + idx[j + 1] + (long)(2 * j);   /* distinct axes */
            PointResult pr = eval_point(e, &ks, t0, atoms.items, avals, atoms.count);
            if (pr == PT_POLE)      { pole = true; break; }
            if (pr == PT_UNDECIDED) { undecided = true; break; }
            if (pr == PT_NONZERO)   { nonzero = true; break; }
        } while (grid_next(idx, deg, naxes));
        free(idx); free(avals);

        if (nonzero)   { verdict = TRIGEXP_ZERO_FALSE;   break; }
        if (undecided) { verdict = TRIGEXP_ZERO_UNKNOWN; break; }
        if (pole)      continue;
        verdict = TRIGEXP_ZERO_TRUE;
        break;
    }

    free(deg); kset_free(&ks); aset_free(&atoms);
    return verdict;
}

/* ------------------------------------------------------------------ */
/*  kernel-variable detection + public entry points                   */
/* ------------------------------------------------------------------ */

typedef struct { char** names; size_t count, cap; } NameSet;
static void nset_add(NameSet* s, const char* name) {
    for (size_t i = 0; i < s->count; i++) if (strcmp(s->names[i], name) == 0) return;
    if (s->count == s->cap) { s->cap = s->cap ? s->cap * 2 : 4;
        s->names = realloc(s->names, s->cap * sizeof(char*)); }
    s->names[s->count++] = mathilda_strdup(name);
}
static void nset_free(NameSet* s) {
    for (size_t i = 0; i < s->count; i++) free(s->names[i]);
    free(s->names);
}
static void collect_leaf_syms(const Expr* e, NameSet* s) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) { if (!is_known_const(e)) nset_add(s, e->data.symbol.name); return; }
    if (e->type == EXPR_FUNCTION)
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            collect_leaf_syms(e->data.function.args[i], s);
}
/* Kernel-variable candidates: free symbols inside a trig/hyp arg or exp exponent. */
static void collect_kernel_vars(const Expr* e, NameSet* s) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const Expr* h = e->data.function.head;
    size_t n = e->data.function.arg_count;
    bool trig = h && h->type == EXPR_SYMBOL &&
        (sym_is(h, SYM_Sin) || sym_is(h, SYM_Cos) ||
         strcmp(h->data.symbol.name, "Tan") == 0 || strcmp(h->data.symbol.name, "Cot") == 0 ||
         strcmp(h->data.symbol.name, "Sec") == 0 || strcmp(h->data.symbol.name, "Csc") == 0 ||
         strcmp(h->data.symbol.name, "Sinh") == 0 || strcmp(h->data.symbol.name, "Cosh") == 0 ||
         strcmp(h->data.symbol.name, "Tanh") == 0 || strcmp(h->data.symbol.name, "Coth") == 0 ||
         strcmp(h->data.symbol.name, "Sech") == 0 || strcmp(h->data.symbol.name, "Csch") == 0);
    if (trig && n == 1) collect_leaf_syms(e->data.function.args[0], s);
    else if (sym_is(h, SYM_Exp) && n == 1) collect_leaf_syms(e->data.function.args[0], s);
    else if (sym_is(h, SYM_Power) && n == 2 && sym_is(e->data.function.args[0], SYM_E))
        collect_leaf_syms(e->data.function.args[1], s);
    for (size_t i = 0; i < n; i++) collect_kernel_vars(e->data.function.args[i], s);
}

TrigExpZeroResult trigexp_rational_is_zero(const Expr* e) {
    if (!e) return TRIGEXP_ZERO_UNKNOWN;
    NameSet kv = {0};
    collect_kernel_vars(e, &kv);
    TrigExpZeroResult r = TRIGEXP_ZERO_UNKNOWN;
    if (kv.count == 1) {
        Expr* var = mk_sym(kv.names[0]);
        r = decide_with_var(e, var);
        expr_free(var);
    }
    nset_free(&kv);
    return r;
}

Expr* transform_trigexp_vanish(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (!sym_is(e->data.function.head, SYM_Plus)) return NULL;
    if (!contains_trig_or_hyperbolic(e) && !contains_exp_form(e)) return NULL;
    if (trigexp_rational_is_zero(e) == TRIGEXP_ZERO_TRUE) return mk_int(0);
    return NULL;
}
