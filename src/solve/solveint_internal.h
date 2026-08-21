/*
 * solveint_internal.h
 *
 * Shared internal substrate for the Solve[..., Integers] engine, split across
 * src/solve/solveint*.c and src/solve/solve_common.c.  Holds the constraint
 * store (SICtx), the leaf-search state (SearchState), the SI_* limits, the
 * tiny Expr-construction / predicate helpers (inline), and the cross-file
 * prototypes for every shared helper and per-method entry point.
 *
 * NOT a public header -- the public surface is solveint.h.  This file is
 * included only by the solveint*.c / solve_common.c translation units.
 */
#ifndef SOLVEINT_INTERNAL_H
#define SOLVEINT_INTERNAL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <gmp.h>
#include "expr.h"
#include "poly/mpoly.h"

/* ------------------------------------------------------------------ *
 *  Limits.                                                            *
 * ------------------------------------------------------------------ */
#define SI_MAX_VARS   10
/* Refuse (return unevaluated) rather than enumerate a search box with more
 * than this many raw nodes -- prevents a hang without ever truncating a
 * result. */
#define SI_MAX_NODES  200000000LL
#define SI_EQ 0
#define SI_LE 1
#define SI_LEAF_MAXDEG 8

/* ------------------------------------------------------------------ *
 *  Constraint store + search state.                                  *
 * ------------------------------------------------------------------ */
typedef struct {
    Expr*  poly_src;   /* borrowed: the (in)equation as an Expr (unused after conv) */
    MPoly* Q;          /* residual polynomial; interpret with `kind` */
    int    kind;       /* SI_EQ (Q == 0, both orientations) or SI_LE (Q <= 0) */
} BoundConstraint;

typedef struct {
    Expr** var;                 /* n borrowed symbols */
    int    n;

    int64_t lo[SI_MAX_VARS];
    int64_t hi[SI_MAX_VARS];
    bool    has_lo[SI_MAX_VARS];
    bool    has_hi[SI_MAX_VARS];

    /* orderings: var[a] < var[b] (strict) or <= (non-strict) */
    int   ord_a[SI_MAX_VARS * SI_MAX_VARS];
    int   ord_b[SI_MAX_VARS * SI_MAX_VARS];
    bool  ord_strict[SI_MAX_VARS * SI_MAX_VARS];
    int   n_ord;

    /* abs-orderings: |var[a]| < |var[b]| (strict) or <= (non-strict).  A
     * separate store because it constrains magnitudes, not signed values:
     * bound propagation and candidate filtering treat it via |.|, and it never
     * reduces to a signed ordering (both signs of each variable survive). */
    int   abs_ord_a[SI_MAX_VARS * SI_MAX_VARS];
    int   abs_ord_b[SI_MAX_VARS * SI_MAX_VARS];
    bool  abs_ord_strict[SI_MAX_VARS * SI_MAX_VARS];
    int   n_abs_ord;

    /* disequations: var[a] != var[b] */
    int   neq_a[SI_MAX_VARS * SI_MAX_VARS];
    int   neq_b[SI_MAX_VARS * SI_MAX_VARS];
    int   n_neq;

    /* equation MPolys used both for solving and for bounding */
    MPoly* eq[SI_MAX_VARS * 2];
    int    neq;

    /* every relation normalised for the bounder */
    BoundConstraint bc[SI_MAX_VARS * 4];
    int    nbc;

    /* True iff every constraint is fully represented by the bound / ordering /
     * disequation store, so a candidate can be checked numerically without a
     * symbolic re-evaluation of the original conjunction. */
    bool   all_captured;

    Expr*  original;            /* borrowed: full conjunction, for verification */
} SICtx;

typedef struct {
    SICtx*  ctx;
    int     leaf;               /* leaf variable index */
    int     order[SI_MAX_VARS]; /* search-var indices, outer-first */
    int     n_search;
    int64_t val[SI_MAX_VARS];   /* current assignment */
    int64_t* sols;              /* flattened n per solution */
    int      nsol, cap;
    int64_t  visits;            /* leaf nodes visited */
    int64_t  max_visits;        /* runtime backstop */
    bool     overflow;

    /* Per-equation int64 coefficient cache (leaf-independent): the int64 forms
     * of eq[q]->coefs, built once so the hot leaf loop never touches GMP.
     * eqc[q] is NULL / eqc_ok[q] false when a coefficient does not fit int64. */
    int64_t* eqc[SI_MAX_VARS * 2];
    bool     eqc_ok[SI_MAX_VARS * 2];
    bool     eqc_built;

    /* Multi-leaf staged elimination (A2): `n_peel` determined variables, each
     * resolved from a single equation given the enumerated free variables.
     * peel[k] is the k-th determined var; peel_eq[k] its resolving equation. */
    bool     multileaf;
    int      peel[SI_MAX_VARS];
    int      peel_eq[SI_MAX_VARS];
    int      n_peel;
} SearchState;

/* ------------------------------------------------------------------ *
 *  Tiny construction helpers (shared, inlined per TU).               *
 * ------------------------------------------------------------------ */
static inline Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static inline Expr* mk_sym(const char* s) { return expr_new_symbol(s); }
static inline Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static inline Expr* mk_fn1(const char* head, Expr* a) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a }, 1);
}
static inline Expr* mk_rule(Expr* lhs, Expr* rhs) { return mk_fn2("Rule", lhs, rhs); }
static inline Expr* mk_list(Expr** args, size_t n) {
    return expr_new_function(mk_sym("List"), args, n);
}
static inline Expr* mk_mpz(const mpz_t z) {
    return expr_bigint_normalize(expr_new_bigint_from_mpz(z));
}

/* ------------------------------------------------------------------ *
 *  Small predicates / accessors (shared, inlined per TU).            *
 * ------------------------------------------------------------------ */
static inline int find_var_index(const SICtx* c, const char* name) {
    for (int i = 0; i < c->n; i++)
        if (c->var[i]->data.symbol.name == name) return i;
    return -1;
}
static inline bool is_sym(const Expr* e, const char* head) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == head;
}
static inline bool is_fun(const Expr* e, const char* head, size_t argc) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == head
        && e->data.function.arg_count == argc;
}
static inline bool expr_as_i64(const Expr* e, int64_t* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = e->data.integer; return true; }
    if (e->type == EXPR_BIGINT) {
        if (mpz_fits_slong_p(e->data.bigint)) { *out = mpz_get_si(e->data.bigint); return true; }
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Cross-file shared helpers and per-method entry points.            *
 * ------------------------------------------------------------------ */
void si_eval_mpoly(const MPoly* p, const int64_t* vals, mpz_t out);
MPoly* relation_to_mpoly(const SICtx* c, Expr* a, Expr* b);
bool classify_conjunct(SICtx* c, Expr* e);
void derive_bounds(SICtx* c);
void derive_even_only_bounds(SICtx* c);
int univariate_roots(const mpz_t* a, int d, int64_t lo, int64_t hi, int64_t* roots);
void si_build_leaf_cache(SearchState* st);
void si_free_leaf_cache(SearchState* st);
int64_t si_isqrt_i64(int64_t n);
void search_rec(SearchState* st, int depth);
int si_longest_chain(const SICtx* c, const int* vars, int nv);
bool si_solve_box_modsieve(SICtx* c, SearchState* st);
Expr* build_result(SearchState* st);
bool si_verify(SICtx* c, const int64_t* vals);
void emit_full(SearchState* st, const int64_t* vals);
bool si_solve_multileaf(SICtx* c, SearchState* st);
bool mitm_solve(SearchState* st);
void ctx_free(SICtx* c);
void flatten_conjuncts(Expr* e, Expr*** out, int* n);
void si_warn_free_symbols(const SICtx* c, Expr** conj, int ncj);
MPoly* si_resid_to_mpoly(Expr* resid, Expr** vars, int n);
bool si_build_total_order(const SICtx* c, const int* active, int ka, int* order);
bool si_solve_reciprocal(SICtx* c, SearchState* st);
bool si_solve_linelim_bilinear(SICtx* c, SearchState* st);
bool si_pell_cf(const mpz_t D, const mpz_t N, mpz_t u, mpz_t v, mpz_t bx, mpz_t by);
bool si_solve_pell(SICtx* c, SearchState* st);
Expr* si_solve_pell_parametric(SICtx* c);
Expr* si_solve_genpell_parametric(SICtx* c);
Expr* si_solve_ternary_quadratic(SICtx* c);
Expr* si_solve_linear_system_ray(SICtx* c);
Expr* si_solve_power_sum_equal(SICtx* c);
Expr* si_solve_mordell(SICtx* c);
Expr* si_solve_linear_parametric(SICtx* c);
Expr* si_solve_linear_system_hnf(SICtx* c);
bool si_solve_linear_bounded(SICtx* c, SearchState* st);
void si_two_power_solve(SICtx* c, SearchState* st, int ip, int jp, int e, const mpz_t m, int64_t* vals);
bool si_solve_powersum_divisor(SICtx* c, SearchState* st);
bool si_solve_conic(SICtx* c, SearchState* st);
bool si_solve_factorable_conic(SICtx* c, SearchState* st);
bool si_solve_elliptic_bqf(SICtx* c, SearchState* st);
int64_t si_isqrt_i128(__int128 n);
Expr* builtin_cube_roots_mod(Expr* res);
bool si_solve_three_cubes_booker(SICtx* c, SearchState* st);
Expr* si_solve_three_cubes_mod9(SICtx* c);
bool si_solve_biquadrate_frye(SICtx* c, SearchState* st);
bool si_solve_separable_mitm(SICtx* c, SearchState* st);
/* Additive power term  coef * (base)^(exp), used by the exponential-Diophantine
 * shape parser (solveint_exp.c) and the Ramanujan-Nagell solver (solveint_rn.c).
 * bvar/evar are variable indices (or -1); bconst/econst the constant base/exp. */
typedef struct { int coef; int bvar; int64_t bconst; int evar; int64_t econst; } SIExpTerm;
bool si_exp_collect(Expr* e, int sign, Expr** var, int n, SIExpTerm* t, int* nt, mpz_t K);
Expr* si_exp_one_tuple(Expr** var, const int64_t* vals, int n);
bool si_verify_symbolic(Expr* expr, Expr** var, const int64_t* vals, int n);
Expr* si_solve_exponential(Expr* expr, Expr** var, int n);
Expr* si_solve_ramanujan_nagell(Expr* expr, Expr** var, int n);
Expr* si_solve_bounded_powerleaf(Expr* expr, Expr** var, int n);
Expr* si_solve_thue(SICtx* c);

#endif /* SOLVEINT_INTERNAL_H */
