/*
 * solve_explog.c -- exact log-linear algebra for exponential equations.
 *
 * See solve_explog.h for the mathematics.  Layout:
 *
 *   1. ExplogBasis   -- a coprime, primitive multiplicative basis of a set of
 *                       positive integers (gcd refinement + perfect-power
 *                       reduction; no factoring).
 *   2. Shape parse   -- split the residual into exponential factors
 *                       b_i^g_i(x) (with sign s_i = +-1 for the side they sit
 *                       on) and the constant K of  prod b_i^(s_i g_i) == K.
 *   3. Generators    -- the rational basis elements p_k, then E (Log = 1),
 *                       then opaque numeric bases (Log[Pi], Log[-2], ...),
 *                       then the residual constant R of K (Log[R]).
 *   4. Matrix        -- M[gen][atom]: the coefficient of Log[gen] * atom(x)
 *                       in  sum_i s_i g_i(x) Log[b_i] - Log[K],  where the
 *                       atoms are the var-dependent monomials of the expanded
 *                       exponents (and 1).  All entries are exact rationals.
 *   5. Reading       -- rank 1: M = w P^T, the equation is the RATIONAL
 *                       polynomial P(x) == 2 Pi I C/Log[prod p^w]; rank >= 2:
 *                       sum_m atom_m * (column m . Log) == 2 Pi I C.  The
 *                       rational decision polynomial is the gcd of the rows.
 *
 * Everything is exact (GMP mpz/mpq).  No numeric approximation is used except
 * to certify the sign of an opaque constant (Pi > 0), which N does exactly
 * enough for a nonzero number.
 */

#include "solve_explog.h"

#include <stdlib.h>
#include <string.h>

#include "eval.h"
#include "expr.h"
#include "reduce_eq.h"   /* reduce_value_reality(): parameters real over Reals */
#include "sym_intern.h"
#include "sym_names.h"

/* ------------------------------------------------------------------ *
 *  Small builders (each consumes its Expr* arguments).                *
 * ------------------------------------------------------------------ */

static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }
static Expr* fn1(const char* h, Expr* a) {
    return expr_new_function(mk_sym(h), (Expr*[]){ a }, 1);
}
static Expr* fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(h), (Expr*[]){ a, b }, 2);
}
static Expr* ev(Expr* e) { return eval_and_free(e); }

static bool head_is(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

static bool is_sym(const Expr* e, const char* sym) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == sym;
}

static bool contains_var(const Expr* e, const Expr* var) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == var->data.symbol.name;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_var(e->data.function.head, var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_var(e->data.function.args[i], var)) return true;
    return false;
}

/* An inexact number anywhere: the exact algebra here does not apply (Solve's
 * own rationalisation pre-pass normally removes these first). */
static bool has_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type != EXPR_FUNCTION) return false;
    if (has_inexact(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_inexact(e->data.function.args[i])) return true;
    return false;
}

static Expr* mpz_expr(const mpz_t z) {
    return expr_bigint_normalize(expr_new_bigint_from_mpz(z));
}

/* Exact rational q as an evaluated Integer/BigInt/Rational. */
static Expr* mpq_expr(const mpq_t q) {
    Expr* n = mpz_expr(mpq_numref(q));
    if (mpz_cmp_ui(mpq_denref(q), 1) == 0) return n;
    Expr* d = mpz_expr(mpq_denref(q));
    return ev(fn2("Times", n, fn2("Power", d, expr_new_integer(-1))));
}

/* n/d (d > 0, coprime) of an exact rational Integer / BigInt / Rational. */
static bool rat_parts(const Expr* e, mpz_t n, mpz_t d) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) {
        expr_to_mpz(e, n);
        mpz_set_ui(d, 1);
        return true;
    }
    if (head_is(e, SYM_Rational) && e->data.function.arg_count == 2
        && expr_is_integer_like(e->data.function.args[0])
        && expr_is_integer_like(e->data.function.args[1])) {
        expr_to_mpz(e->data.function.args[0], n);
        expr_to_mpz(e->data.function.args[1], d);
        if (mpz_sgn(d) == 0) return false;
        if (mpz_sgn(d) < 0) { mpz_neg(n, n); mpz_neg(d, d); }
        return true;
    }
    return false;
}

/* z = v for any int64_t (mpz_set_si takes a long, 32-bit on LLP64). */
static void mpz_set_i64(mpz_t z, int64_t v) {
    uint64_t a = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
    mpz_set_ui(z, (unsigned long)(a >> 32));
    mpz_mul_2exp(z, z, 32);
    mpz_add_ui(z, z, (unsigned long)(a & UINT64_C(0xffffffff)));
    if (v < 0) mpz_neg(z, z);
}

static bool rat_mpq(const Expr* e, mpq_t q) {
    mpz_t n, d;
    mpz_inits(n, d, NULL);
    bool ok = rat_parts(e, n, d);
    if (ok) {
        mpz_set(mpq_numref(q), n);
        mpz_set(mpq_denref(q), d);
        mpq_canonicalize(q);
    }
    mpz_clears(n, d, NULL);
    return ok;
}

static int const_sign(const Expr* c);

/* ================================================================== *
 *  1. Multiplicative basis.                                           *
 * ================================================================== */

/* Stack-safety bound for mpz_perfect_power_p (it allocates on the stack in
 * proportion to the operand); matches power.c's find_min_perfect_base. */
#define EXPLOG_PERFECT_POWER_MAX_BITS 200000u

/* Reduce n > 1 in place to its primitive root r (n = r^m, r not a perfect
 * power).  Only prime root orders need trying, none above the bit length. */
static void int_primitive_root(mpz_t n) {
    if (mpz_cmp_ui(n, 1) <= 0) return;
    mpz_t r;
    mpz_init(r);
    for (;;) {
        size_t bits = mpz_sizeinbase(n, 2);
        if (bits > EXPLOG_PERFECT_POWER_MAX_BITS || !mpz_perfect_power_p(n)) break;
        bool found = false;
        for (unsigned long p = 2; p <= bits && !found; p++) {
            bool prime = true;
            for (unsigned long q = 2; q * q <= p; q++)
                if (p % q == 0) { prime = false; break; }
            if (prime && mpz_root(r, n, p)) { mpz_set(n, r); found = true; }
        }
        if (!found) break;
    }
    mpz_clear(r);
}

typedef struct { mpz_t* v; size_t n, cap; } MpzList;

static void ml_push(MpzList* L, const mpz_t x) {
    if (mpz_cmp_ui(x, 1) <= 0) return;
    for (size_t i = 0; i < L->n; i++) if (mpz_cmp(L->v[i], x) == 0) return;
    if (L->n == L->cap) {
        L->cap = L->cap ? 2 * L->cap : 8;
        L->v = (mpz_t*)realloc(L->v, L->cap * sizeof(mpz_t));
    }
    mpz_init_set(L->v[L->n++], x);
}

static int cmp_mpz(const void* a, const void* b) {
    return mpz_cmp(*(const mpz_t*)a, *(const mpz_t*)b);
}

void explog_basis_build(const mpz_t* vals, size_t nv, ExplogBasis* B) {
    MpzList L = { NULL, 0, 0 };
    for (size_t i = 0; i < nv; i++) ml_push(&L, vals[i]);

    /* gcd refinement: while two elements share a factor g, replace them by
     * their g-free parts and g.  The product of the list strictly decreases
     * (a b -> (a/g^i)(b/g^j) g with i, j >= 1), so this terminates; the fixed
     * point is pairwise coprime and still generates every input. */
    mpz_t g, a, b;
    mpz_inits(g, a, b, NULL);
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < L.n && !changed; i++) {
            for (size_t j = i + 1; j < L.n && !changed; j++) {
                mpz_gcd(g, L.v[i], L.v[j]);
                if (mpz_cmp_ui(g, 1) == 0) continue;
                mpz_remove(a, L.v[i], g);
                mpz_remove(b, L.v[j], g);
                /* Rebuild the list without i, j; then add a, b, g. */
                MpzList N = { NULL, 0, 0 };
                for (size_t k = 0; k < L.n; k++)
                    if (k != i && k != j) ml_push(&N, L.v[k]);
                ml_push(&N, a); ml_push(&N, b); ml_push(&N, g);
                for (size_t k = 0; k < L.n; k++) mpz_clear(L.v[k]);
                free(L.v);
                L = N;
                changed = true;
            }
        }
    }
    mpz_clears(g, a, b, NULL);

    /* Primitive roots: coprime elements stay coprime (and distinct). */
    for (size_t i = 0; i < L.n; i++) int_primitive_root(L.v[i]);
    if (L.n > 1) qsort(L.v, L.n, sizeof(mpz_t), cmp_mpz);

    B->n = L.n;
    B->p = L.v;
}

bool explog_basis_exponents(const ExplogBasis* B, const mpz_t v, int64_t* e) {
    mpz_t t;
    mpz_init_set(t, v);
    for (size_t k = 0; k < B->n; k++)
        e[k] = (int64_t)mpz_remove(t, t, B->p[k]);
    bool ok = mpz_cmp_ui(t, 1) == 0;
    mpz_clear(t);
    return ok;
}

void explog_basis_clear(ExplogBasis* B) {
    if (!B) return;
    for (size_t k = 0; k < B->n; k++) mpz_clear(B->p[k]);
    free(B->p);
    B->p = NULL;
    B->n = 0;
}

/* ================================================================== *
 *  2. Shape: exponential factors and the constant side.               *
 * ================================================================== */

/* Classification of a variable-free base (or constant factor). */
typedef enum { BK_RAT, BK_RATPOW, BK_E, BK_OPAQUE } BaseKind;

typedef struct {
    BaseKind kind;
    mpz_t n, d;      /* BK_RAT / BK_RATPOW: the positive rational n/d      */
    mpq_t r;         /* BK_RATPOW: the rational power (n/d)^r              */
    const Expr* b;   /* borrowed: the base itself                          */
} BaseInfo;

typedef struct {
    BaseInfo base;
    Expr* g;         /* owned exponent g_i(x)                              */
    int sign;        /* +1: left side, -1: moved from the right            */
} ExpFac;

typedef struct { ExpFac* f; size_t n, cap; } FacList;

static void fl_free(FacList* F) {
    for (size_t i = 0; i < F->n; i++) {
        mpz_clears(F->f[i].base.n, F->f[i].base.d, NULL);
        mpq_clear(F->f[i].base.r);
        expr_free(F->f[i].g);
    }
    free(F->f);
    F->f = NULL; F->n = F->cap = 0;
}

/* Classify a variable-free base.  false: not numeric (a symbolic parameter
 * base is left to the generic inverse-function path) or zero. */
static bool classify_base(const Expr* b, BaseInfo* bi) {
    mpz_inits(bi->n, bi->d, NULL);
    mpq_init(bi->r);
    bi->b = b;
    if (rat_parts(b, bi->n, bi->d)) {
        if (mpz_sgn(bi->n) > 0) { bi->kind = BK_RAT; return true; }
        if (mpz_sgn(bi->n) == 0) goto fail;
        bi->kind = BK_OPAQUE;                       /* negative rational */
        return true;
    }
    if (is_sym(b, SYM_E)) { bi->kind = BK_E; return true; }
    if (head_is(b, SYM_Power) && b->data.function.arg_count == 2
        && rat_parts(b->data.function.args[0], bi->n, bi->d)
        && mpz_sgn(bi->n) > 0
        && rat_mpq(b->data.function.args[1], bi->r)) {
        bi->kind = BK_RATPOW;
        return true;
    }
    {
        Expr* nq = ev(fn1("NumericQ", expr_copy((Expr*)b)));
        bool numeric = is_sym(nq, SYM_True);
        expr_free(nq);
        if (numeric) { bi->kind = BK_OPAQUE; return true; }
    }
fail:
    mpz_clears(bi->n, bi->d, NULL);
    mpq_clear(bi->r);
    return false;
}

static void fl_push(FacList* F, const BaseInfo* bi, Expr* g, int sign) {
    if (F->n == F->cap) {
        F->cap = F->cap ? 2 * F->cap : 8;
        F->f = (ExpFac*)realloc(F->f, F->cap * sizeof(ExpFac));
    }
    F->f[F->n].base = *bi;   /* takes the mpz/mpq handles */
    F->f[F->n].g = g;
    F->f[F->n].sign = sign;
    F->n++;
}

/* Split one variable-carrying term into exponential factors (pushed onto F
 * with `sign`) and its constant coefficient (*cst, owned, evaluated).  false
 * when some factor carrying the variable is not b^g(x) with a numeric b. */
static bool split_term(const Expr* t, const Expr* var, int sign,
                       FacList* F, Expr** cst) {
    size_t nf = head_is(t, SYM_Times) ? t->data.function.arg_count : 1;
    Expr** consts = (Expr**)malloc(sizeof(Expr*) * (nf ? nf : 1));
    size_t nc = 0;
    bool ok = true;
    for (size_t i = 0; i < nf && ok; i++) {
        const Expr* f = head_is(t, SYM_Times) ? t->data.function.args[i] : t;
        if (!contains_var(f, var)) { consts[nc++] = expr_copy((Expr*)f); continue; }
        const Expr* base = NULL;
        Expr* g = NULL;
        if (head_is(f, SYM_Power) && f->data.function.arg_count == 2
            && !contains_var(f->data.function.args[0], var)) {
            base = f->data.function.args[0];
            g = expr_copy(f->data.function.args[1]);
        } else if (head_is(f, SYM_Power) && f->data.function.arg_count == 2
                   && f->data.function.args[1]->type == EXPR_INTEGER
                   && head_is(f->data.function.args[0], SYM_Power)
                   && f->data.function.args[0]->data.function.arg_count == 2
                   && !contains_var(f->data.function.args[0]->data.function.args[0], var)) {
            /* (b^g)^n == b^(g n) for integer n. */
            const Expr* in = f->data.function.args[0];
            base = in->data.function.args[0];
            g = ev(fn2("Times", expr_copy(in->data.function.args[1]),
                       expr_copy(f->data.function.args[1])));
        }
        BaseInfo bi;
        if (!base || !classify_base(base, &bi)) { if (g) expr_free(g); ok = false; break; }
        fl_push(F, &bi, g, sign);
    }
    if (!ok) {
        for (size_t i = 0; i < nc; i++) expr_free(consts[i]);
        free(consts);
        return false;
    }
    *cst = nc == 0 ? expr_new_integer(1)
         : nc == 1 ? consts[0]
         : ev(expr_new_function(mk_sym("Times"), consts, nc));
    free(consts);   /* elements consumed / moved; the array is ours */
    return true;
}

/* ================================================================== *
 *  3./4. Generators and the coefficient matrix.                       *
 * ================================================================== */

typedef enum { GK_RAT, GK_E, GK_OPAQUE } GenKind;

typedef struct {
    GenKind kind;
    Expr* key;    /* GK_OPAQUE: the base (or the constant R), owned        */
    Expr* logv;   /* GK_OPAQUE: its evaluated Log, owned                   */
} Gen;

typedef struct {
    ExplogBasis B;          /* rational generators are 0 .. B.n-1          */
    Gen* gen;               /* all generators (rational ones have no key)  */
    size_t ng, cap;
    long ge;                /* index of the E generator, -1 if none        */
} GenTable;

static size_t gt_add_opaque(GenTable* T, const Expr* key) {
    for (size_t i = T->B.n; i < T->ng; i++)
        if (T->gen[i].kind == GK_OPAQUE && expr_eq(T->gen[i].key, key)) return i;
    if (T->ng == T->cap) {
        T->cap = T->cap ? 2 * T->cap : 8;
        T->gen = (Gen*)realloc(T->gen, T->cap * sizeof(Gen));
    }
    T->gen[T->ng].kind = GK_OPAQUE;
    T->gen[T->ng].key = expr_copy((Expr*)key);
    T->gen[T->ng].logv = ev(fn1("Log", expr_copy((Expr*)key)));
    return T->ng++;
}

static size_t gt_add_e(GenTable* T) {
    if (T->ge >= 0) return (size_t)T->ge;
    if (T->ng == T->cap) {
        T->cap = T->cap ? 2 * T->cap : 8;
        T->gen = (Gen*)realloc(T->gen, T->cap * sizeof(Gen));
    }
    T->gen[T->ng].kind = GK_E;
    T->gen[T->ng].key = NULL;
    T->gen[T->ng].logv = NULL;
    T->ge = (long)T->ng;
    return T->ng++;
}

static void gt_free(GenTable* T) {
    for (size_t i = 0; i < T->ng; i++) {
        if (T->gen[i].key) expr_free(T->gen[i].key);
        if (T->gen[i].logv) expr_free(T->gen[i].logv);
    }
    free(T->gen);
    explog_basis_clear(&T->B);
}

/* A dense rational vector over the generators. */
typedef struct { mpq_t* q; size_t n; } QVec;

static void qv_init(QVec* v, size_t n) {
    v->n = n;
    v->q = (mpq_t*)malloc(sizeof(mpq_t) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) mpq_init(v->q[i]);
}
static void qv_clear(QVec* v) {
    for (size_t i = 0; i < v->n; i++) mpq_clear(v->q[i]);
    free(v->q);
    v->q = NULL; v->n = 0;
}
static bool qv_is_zero(const QVec* v) {
    for (size_t i = 0; i < v->n; i++) if (mpq_sgn(v->q[i]) != 0) return false;
    return true;
}

/* v += scale * (vector of the positive rational n/d over the basis). */
static bool qv_add_rat(QVec* v, const ExplogBasis* B, const mpz_t n,
                       const mpz_t d, const mpq_t scale) {
    int64_t* en = (int64_t*)calloc(B->n ? B->n : 1, sizeof(int64_t));
    int64_t* ed = (int64_t*)calloc(B->n ? B->n : 1, sizeof(int64_t));
    bool ok = explog_basis_exponents(B, n, en) && explog_basis_exponents(B, d, ed);
    if (ok) {
        mpq_t t;
        mpq_init(t);
        for (size_t k = 0; k < B->n; k++) {
            int64_t e = en[k] - ed[k];
            if (e == 0) continue;
            mpz_set_i64(mpq_numref(t), e);
            mpz_set_ui(mpq_denref(t), 1);
            mpq_mul(t, t, scale);
            mpq_add(v->q[k], v->q[k], t);
        }
        mpq_clear(t);
    }
    free(en); free(ed);
    return ok;
}

/* v += scale * Log-vector of a classified base. */
static bool qv_add_base(QVec* v, GenTable* T, const BaseInfo* bi, const mpq_t scale) {
    switch (bi->kind) {
    case BK_RAT:
        return qv_add_rat(v, &T->B, bi->n, bi->d, scale);
    case BK_RATPOW: {
        mpq_t s;
        mpq_init(s);
        mpq_mul(s, scale, bi->r);
        bool ok = qv_add_rat(v, &T->B, bi->n, bi->d, s);
        mpq_clear(s);
        return ok;
    }
    case BK_E: {
        size_t k = gt_add_e(T);
        if (k >= v->n) return false;
        mpq_add(v->q[k], v->q[k], scale);
        return true;
    }
    case BK_OPAQUE: {
        size_t k = gt_add_opaque(T, bi->b);
        if (k >= v->n) return false;
        mpq_add(v->q[k], v->q[k], scale);
        return true;
    }
    }
    return false;
}

/* Decomposition of the constant side K = q * prod (n/d)^r * E^e * R. */
typedef struct {
    FacList pw;      /* BK_RAT / BK_RATPOW / BK_E factors, g = NULL        */
    Expr* R;         /* the rest (sign, I, Pi, parameters, ...), owned     */
} ConstDecomp;

static void cd_free(ConstDecomp* C) {
    for (size_t i = 0; i < C->pw.n; i++) {
        mpz_clears(C->pw.f[i].base.n, C->pw.f[i].base.d, NULL);
        mpq_clear(C->pw.f[i].base.r);
    }
    free(C->pw.f);
    if (C->R) expr_free(C->R);
}

static void cd_build(const Expr* K, ConstDecomp* C) {
    C->pw.f = NULL; C->pw.n = C->pw.cap = 0;
    size_t nf = head_is(K, SYM_Times) ? K->data.function.arg_count : 1;
    Expr** rest = (Expr**)malloc(sizeof(Expr*) * (nf + 1));
    size_t nr = 0;
    for (size_t i = 0; i < nf; i++) {
        const Expr* f = head_is(K, SYM_Times) ? K->data.function.args[i] : K;
        BaseInfo bi;
        mpz_inits(bi.n, bi.d, NULL);
        mpq_init(bi.r);
        bi.b = f;
        if (rat_parts(f, bi.n, bi.d) && mpz_sgn(bi.n) != 0) {
            if (mpz_sgn(bi.n) < 0) { mpz_neg(bi.n, bi.n); rest[nr++] = expr_new_integer(-1); }
            if (mpz_cmp(bi.n, bi.d) != 0) {
                bi.kind = BK_RAT;
                fl_push(&C->pw, &bi, NULL, 1);
                continue;
            }
        } else if (head_is(f, SYM_Power) && f->data.function.arg_count == 2
                   && rat_parts(f->data.function.args[0], bi.n, bi.d)
                   && mpz_sgn(bi.n) > 0
                   && rat_mpq(f->data.function.args[1], bi.r)) {
            bi.kind = BK_RATPOW;
            fl_push(&C->pw, &bi, NULL, 1);
            continue;
        } else if (head_is(f, SYM_Power) && f->data.function.arg_count == 2
                   && is_sym(f->data.function.args[0], SYM_E)
                   && rat_mpq(f->data.function.args[1], bi.r)) {
            bi.kind = BK_E;
            fl_push(&C->pw, &bi, NULL, 1);
            continue;
        } else if (is_sym(f, SYM_E)) {
            mpq_set_ui(bi.r, 1, 1);
            bi.kind = BK_E;
            fl_push(&C->pw, &bi, NULL, 1);
            continue;
        } else {
            /* A positive real numeric constant b, or a rational power b^r of
             * one, is its own (opaque) generator: Log[b^r] = r Log[b] holds
             * exactly for b > 0, so Pi^x == Pi^2 reads x - 2 == 0. */
            const Expr* b = f;
            mpq_set_ui(bi.r, 1, 1);
            if (head_is(f, SYM_Power) && f->data.function.arg_count == 2
                && rat_mpq(f->data.function.args[1], bi.r))
                b = f->data.function.args[0];
            else
                mpq_set_ui(bi.r, 1, 1);
            Expr* nq = ev(fn1("NumericQ", expr_copy((Expr*)b)));
            bool numeric = is_sym(nq, SYM_True);
            expr_free(nq);
            if (numeric && const_sign(b) == 1) {
                bi.kind = BK_OPAQUE;
                bi.b = b;
                fl_push(&C->pw, &bi, NULL, 1);
                continue;
            }
            rest[nr++] = expr_copy((Expr*)f);
        }
        mpz_clears(bi.n, bi.d, NULL);
        mpq_clear(bi.r);
    }
    C->R = nr == 0 ? expr_new_integer(1)
         : ev(expr_new_function(mk_sym("Times"), rest, nr));
    free(rest);
}

/* Sign of a variable-free constant: 1 positive, -1 negative, 2 non-real,
 * 0 undecided (symbolic). */
static int const_sign(const Expr* c) {
    Expr* nv = ev(fn1("N", expr_copy((Expr*)c)));
    int s = 0;
    if (nv->type == EXPR_INTEGER) s = nv->data.integer > 0 ? 1 : nv->data.integer < 0 ? -1 : 0;
    else if (nv->type == EXPR_REAL) s = nv->data.real > 0 ? 1 : nv->data.real < 0 ? -1 : 0;
    else if (nv->type == EXPR_BIGINT) s = mpz_sgn(nv->data.bigint);
#ifdef USE_MPFR
    else if (nv->type == EXPR_MPFR) s = mpfr_sgn(nv->data.mpfr);
#endif
    else if (head_is(nv, SYM_Complex) && nv->data.function.arg_count == 2) {
        const Expr* im = nv->data.function.args[1];
        bool nz = (im->type == EXPR_REAL && im->data.real != 0.0)
               || (im->type == EXPR_INTEGER && im->data.integer != 0);
        s = nz ? 2 : 0;
    }
    expr_free(nv);
    return s;
}

/* Term of an expanded exponent: coefficient (rational) x atom. */
typedef struct { Expr** a; size_t n, cap; } AtomList;   /* a[0] == NULL: 1 */

static size_t al_index(AtomList* A, Expr* atom /* consumed, NULL = 1 */) {
    for (size_t i = 0; i < A->n; i++) {
        if ((!atom && !A->a[i]) || (atom && A->a[i] && expr_eq(atom, A->a[i]))) {
            if (atom) expr_free(atom);
            return i;
        }
    }
    if (A->n == A->cap) {
        A->cap = A->cap ? 2 * A->cap : 8;
        A->a = (Expr**)realloc(A->a, A->cap * sizeof(Expr*));
    }
    A->a[A->n] = atom;
    return A->n++;
}

static void al_free(AtomList* A) {
    for (size_t i = 0; i < A->n; i++) if (A->a[i]) expr_free(A->a[i]);
    free(A->a);
}

/* Polynomial degree of a monomial atom in var: 0 for the constant atom,
 * k >= 1 for var^k, -1 otherwise (non-polynomial or parametric). */
static int atom_degree(const Expr* a, const Expr* var) {
    if (!a) return 0;
    if (a->type == EXPR_SYMBOL && a->data.symbol.name == var->data.symbol.name) return 1;
    if (head_is(a, SYM_Power) && a->data.function.arg_count == 2
        && a->data.function.args[0]->type == EXPR_SYMBOL
        && a->data.function.args[0]->data.symbol.name == var->data.symbol.name
        && a->data.function.args[1]->type == EXPR_INTEGER
        && a->data.function.args[1]->data.integer >= 1
        && a->data.function.args[1]->data.integer <= 1000000)
        return (int)a->data.function.args[1]->data.integer;
    return -1;
}

/* Is the atom real-valued for real var (where defined)?  var^k for any
 * integer k, times variable-free factors that are real -- parameters are real
 * over the Reals, as for every Reals-domain solve (reduce_value_reality). */
static bool atom_real(const Expr* a, const Expr* var) {
    if (!a) return true;
    size_t n = head_is(a, SYM_Times) ? a->data.function.arg_count : 1;
    for (size_t i = 0; i < n; i++) {
        const Expr* f = head_is(a, SYM_Times) ? a->data.function.args[i] : a;
        if (f->type == EXPR_SYMBOL && f->data.symbol.name == var->data.symbol.name) continue;
        if (head_is(f, SYM_Power) && f->data.function.arg_count == 2
            && f->data.function.args[0]->type == EXPR_SYMBOL
            && f->data.function.args[0]->data.symbol.name == var->data.symbol.name
            && f->data.function.args[1]->type == EXPR_INTEGER) continue;
        if (contains_var(f, var)) return false;
        Expr* cond = NULL;
        int r = reduce_value_reality(f, &cond);
        if (cond) expr_free(cond);
        if (r != 1) return false;
    }
    return true;
}

/* Expand g and decompose it into rational coefficients over atoms:
 * coef[idx] += scale * c for each term c * atom.  false: some coefficient is
 * not an exact rational (then the general path is used). */
typedef struct { size_t* idx; mpq_t* c; size_t n; } AtomTerms;

static bool atom_terms(const Expr* g, const Expr* var, AtomList* A, AtomTerms* out) {
    Expr* x = ev(fn1("Expand", expr_copy((Expr*)g)));
    size_t nt = head_is(x, SYM_Plus) ? x->data.function.arg_count : 1;
    out->idx = (size_t*)malloc(sizeof(size_t) * nt);
    out->c = (mpq_t*)malloc(sizeof(mpq_t) * nt);
    out->n = 0;
    bool ok = true;
    for (size_t i = 0; i < nt && ok; i++) {
        const Expr* t = head_is(x, SYM_Plus) ? x->data.function.args[i] : x;
        mpq_init(out->c[out->n]);
        if (!contains_var(t, var)) {
            if (!rat_mpq(t, out->c[out->n])) { mpq_clear(out->c[out->n]); ok = false; break; }
            out->idx[out->n++] = al_index(A, NULL);
            continue;
        }
        Expr* atom = NULL;
        mpq_set_ui(out->c[out->n], 1, 1);
        if (head_is(t, SYM_Times) && t->data.function.arg_count >= 2
            && rat_mpq(t->data.function.args[0], out->c[out->n])) {
            size_t m = t->data.function.arg_count - 1;
            if (m == 1) atom = expr_copy(t->data.function.args[1]);
            else {
                Expr** fs = (Expr**)malloc(sizeof(Expr*) * m);
                for (size_t j = 0; j < m; j++) fs[j] = expr_copy(t->data.function.args[j + 1]);
                atom = expr_new_function(mk_sym("Times"), fs, m);
                free(fs);
            }
        } else {
            atom = expr_copy((Expr*)t);
        }
        out->idx[out->n++] = al_index(A, atom);
    }
    expr_free(x);
    if (!ok) {
        for (size_t i = 0; i < out->n; i++) mpq_clear(out->c[i]);
        free(out->idx); free(out->c);
        out->idx = NULL; out->c = NULL; out->n = 0;
    }
    return ok;
}

static void atom_terms_free(AtomTerms* t) {
    for (size_t i = 0; i < t->n; i++) mpq_clear(t->c[i]);
    free(t->idx); free(t->c);
}

/* ================================================================== *
 *  5. Reading the matrix back as expressions.                         *
 * ================================================================== */

/* The logarithm  sum_k v_k Log[gen_k]  as an expression.  Over the rational
 * generators an integer vector is written as ONE logarithm, Log[prod p^v]
 * (Log[8] rather than 3 Log[2]), with a value below 1 written as -Log[1/q]. */
static Expr* log_combo(const GenTable* T, const QVec* v) {
    size_t cap = T->ng + 2, nt = 0;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * cap);
    bool all_int = true;
    for (size_t k = 0; k < T->B.n; k++)
        if (mpz_cmp_ui(mpq_denref(v->q[k]), 1) != 0) { all_int = false; break; }
    if (all_int) {
        mpz_t num, den, pw;
        mpz_inits(num, den, pw, NULL);
        mpz_set_ui(num, 1); mpz_set_ui(den, 1);
        for (size_t k = 0; k < T->B.n; k++) {
            int sg = mpq_sgn(v->q[k]);
            if (sg == 0) continue;
            mpz_abs(pw, mpq_numref(v->q[k]));
            unsigned long e = mpz_get_ui(pw);
            mpz_pow_ui(pw, T->B.p[k], e);
            if (sg > 0) mpz_mul(num, num, pw); else mpz_mul(den, den, pw);
        }
        int c = mpz_cmp(num, den);
        if (c != 0) {
            mpq_t q;
            mpq_init(q);
            if (c > 0) { mpz_set(mpq_numref(q), num); mpz_set(mpq_denref(q), den); }
            else       { mpz_set(mpq_numref(q), den); mpz_set(mpq_denref(q), num); }
            mpq_canonicalize(q);
            Expr* lg = fn1("Log", mpq_expr(q));
            terms[nt++] = c > 0 ? lg : fn2("Times", expr_new_integer(-1), lg);
            mpq_clear(q);
        }
        mpz_clears(num, den, pw, NULL);
    } else {
        for (size_t k = 0; k < T->B.n; k++) {
            if (mpq_sgn(v->q[k]) == 0) continue;
            terms[nt++] = fn2("Times", mpq_expr(v->q[k]), fn1("Log", mpz_expr(T->B.p[k])));
        }
    }
    for (size_t k = T->B.n; k < T->ng; k++) {
        if (mpq_sgn(v->q[k]) == 0) continue;
        if (T->gen[k].kind == GK_E) terms[nt++] = mpq_expr(v->q[k]);
        else terms[nt++] = fn2("Times", mpq_expr(v->q[k]), expr_copy(T->gen[k].logv));
    }
    Expr* out = nt == 0 ? expr_new_integer(0)
              : ev(expr_new_function(mk_sym("Plus"), terms, nt));
    free(terms);
    return out;
}

/* True when every nonzero entry of v is on an independent generator
 * (rational basis element or E) -- then v . Log == 0 iff v == 0. */
static bool qv_independent(const GenTable* T, const QVec* v) {
    for (size_t k = T->B.n; k < T->ng; k++)
        if (T->gen[k].kind == GK_OPAQUE && mpq_sgn(v->q[k]) != 0) return false;
    return true;
}

/* sum_m coef_m atom_m. */
static Expr* atoms_poly(const AtomList* A, mpq_t* coef) {
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (A->n ? A->n : 1));
    size_t nt = 0;
    for (size_t m = 0; m < A->n; m++) {
        if (mpq_sgn(coef[m]) == 0) continue;
        Expr* c = mpq_expr(coef[m]);
        terms[nt++] = A->a[m] ? fn2("Times", c, expr_copy(A->a[m])) : c;
    }
    Expr* out = nt == 0 ? expr_new_integer(0)
              : ev(expr_new_function(mk_sym("Plus"), terms, nt));
    free(terms);
    return out;
}

static Expr* two_pi_i(void) {
    return expr_new_function(mk_sym("Times"),
        (Expr*[]){ expr_new_integer(2), mk_sym("I"), mk_sym("Pi") }, 3);
}

/* ================================================================== *
 *  explog_linearize                                                    *
 * ================================================================== */

void explog_linear_free(ExplogLinear* L) {
    if (!L) return;
    if (L->lhs) expr_free(L->lhs);
    if (L->period) expr_free(L->period);
    if (L->real_cond) expr_free(L->real_cond);
    if (L->rational_poly) expr_free(L->rational_poly);
    memset(L, 0, sizeof(*L));
}

/* Exponent real-valued on the reals: a polynomial in var whose coefficients
 * are real numbers (the general path's test; the rational path reads it off
 * the atoms). */
static bool exponent_real_poly(const Expr* g, const Expr* var) {
    Expr* pq = ev(fn2("PolynomialQ", expr_copy((Expr*)g), expr_copy((Expr*)var)));
    bool poly = is_sym(pq, SYM_True);
    expr_free(pq);
    if (!poly) return false;
    Expr* cl = ev(fn2("CoefficientList", expr_copy((Expr*)g), expr_copy((Expr*)var)));
    bool ok = head_is(cl, SYM_List);
    for (size_t i = 0; ok && i < cl->data.function.arg_count; i++) {
        const Expr* c = cl->data.function.args[i];
        if (c->type == EXPR_INTEGER) continue;
        int s = const_sign(c);
        if (s == 1 || s == -1) continue;
        Expr* cond = NULL;
        int r = reduce_value_reality(c, &cond);   /* parameters are real */
        if (cond) expr_free(cond);
        if (r != 1) ok = false;
    }
    expr_free(cl);
    return ok;
}

bool explog_linearize(const Expr* residual, const Expr* var, ExplogLinear* out) {
    memset(out, 0, sizeof(*out));
    out->kind = EXPLOG_NO_SHAPE;
    if (!residual || !var || var->type != EXPR_SYMBOL) return false;
    if (!contains_var(residual, var) || has_inexact(residual)) return false;

    /* ---- shape: var terms and the free part ---------------------------- */
    size_t nt = head_is(residual, SYM_Plus) ? residual->data.function.arg_count : 1;
    const Expr* vt[2] = { NULL, NULL };
    size_t nv = 0;
    Expr** fr = (Expr**)malloc(sizeof(Expr*) * nt);
    size_t nfr = 0;
    for (size_t i = 0; i < nt; i++) {
        const Expr* t = head_is(residual, SYM_Plus) ? residual->data.function.args[i] : residual;
        if (contains_var(t, var)) {
            if (nv == 2) { nv = 3; break; }
            vt[nv++] = t;
        } else {
            fr[nfr++] = expr_copy((Expr*)t);
        }
    }
    if (nv == 0 || nv > 2 || (nv == 2 && nfr > 0)) {
        for (size_t i = 0; i < nfr; i++) expr_free(fr[i]);
        free(fr);
        return false;
    }
    Expr* freep = nfr == 0 ? expr_new_integer(0)
                : ev(expr_new_function(mk_sym("Plus"), fr, nfr));
    free(fr);

    FacList F = { NULL, 0, 0 };
    Expr* c1 = NULL;
    Expr* c2 = NULL;
    bool ok = split_term(vt[0], var, 1, &F, &c1);
    if (ok && nv == 2) ok = split_term(vt[1], var, -1, &F, &c2);
    if (!ok || F.n == 0) {
        if (c1) expr_free(c1);
        if (c2) expr_free(c2);
        expr_free(freep);
        fl_free(&F);
        return false;
    }

    /* prod b_i^(s_i g_i) == K. */
    Expr* K;
    if (nv == 2) {
        K = ev(fn2("Times", fn2("Times", expr_new_integer(-1), c2),
                   fn2("Power", c1, expr_new_integer(-1))));
    } else {
        bool zero = freep->type == EXPR_INTEGER && freep->data.integer == 0;
        if (zero) {
            /* c prod b^g == 0 with nonzero numeric bases: no solution. */
            expr_free(c1); expr_free(freep); fl_free(&F);
            out->kind = EXPLOG_NO_SOLUTION;
            return true;
        }
        K = ev(fn2("Times", fn2("Times", expr_new_integer(-1), expr_copy(freep)),
                   fn2("Power", c1, expr_new_integer(-1))));
    }
    expr_free(freep);
    if (contains_var(K, var)) { expr_free(K); fl_free(&F); return false; }

    ConstDecomp CD;
    cd_build(K, &CD);          /* borrows factors of K: K lives until cd_free */

    /* ---- generators -------------------------------------------------- */
    GenTable T;
    memset(&T, 0, sizeof(T));
    T.ge = -1;
    {
        size_t cap = 2 * (F.n + CD.pw.n) + 1, nvals = 0;
        mpz_t* vals = (mpz_t*)malloc(sizeof(mpz_t) * cap);
        for (size_t i = 0; i < F.n; i++) {
            const BaseInfo* bi = &F.f[i].base;
            if (bi->kind == BK_RAT || bi->kind == BK_RATPOW) {
                mpz_init_set(vals[nvals++], bi->n);
                mpz_init_set(vals[nvals++], bi->d);
            }
        }
        for (size_t i = 0; i < CD.pw.n; i++) {
            const BaseInfo* bi = &CD.pw.f[i].base;
            if (bi->kind == BK_RAT || bi->kind == BK_RATPOW) {
                mpz_init_set(vals[nvals++], bi->n);
                mpz_init_set(vals[nvals++], bi->d);
            }
        }
        explog_basis_build((const mpz_t*)vals, nvals, &T.B);
        for (size_t i = 0; i < nvals; i++) mpz_clear(vals[i]);
        free(vals);
    }
    T.cap = T.B.n + F.n + 4;
    T.gen = (Gen*)calloc(T.cap, sizeof(Gen));
    for (size_t k = 0; k < T.B.n; k++) T.gen[k].kind = GK_RAT;
    T.ng = T.B.n;
    /* Pre-register the E / opaque generators so every vector has full width. */
    for (size_t i = 0; i < F.n; i++) {
        if (F.f[i].base.kind == BK_E) gt_add_e(&T);
        else if (F.f[i].base.kind == BK_OPAQUE) gt_add_opaque(&T, F.f[i].base.b);
    }
    for (size_t i = 0; i < CD.pw.n; i++) {
        if (CD.pw.f[i].base.kind == BK_E) gt_add_e(&T);
        else if (CD.pw.f[i].base.kind == BK_OPAQUE) gt_add_opaque(&T, CD.pw.f[i].base.b);
    }
    long gR = -1;
    if (!(CD.R->type == EXPR_INTEGER && CD.R->data.integer == 1))
        gR = (long)gt_add_opaque(&T, CD.R);
    size_t ng = T.ng;

    /* Base vectors and the constant vector  -Log[K]. */
    QVec* vb = (QVec*)malloc(sizeof(QVec) * F.n);
    QVec vk;
    qv_init(&vk, ng);
    mpq_t one, mone;
    mpq_inits(one, mone, NULL);
    mpq_set_si(one, 1, 1);
    mpq_set_si(mone, -1, 1);
    ok = true;
    for (size_t i = 0; i < F.n; i++) {
        qv_init(&vb[i], ng);
        mpq_t s;
        mpq_init(s);
        mpq_set_si(s, F.f[i].sign, 1);
        if (!qv_add_base(&vb[i], &T, &F.f[i].base, s)) ok = false;
        mpq_clear(s);
    }
    for (size_t i = 0; i < CD.pw.n && ok; i++) {
        const BaseInfo* bi = &CD.pw.f[i].base;
        if (bi->kind == BK_E || bi->kind == BK_OPAQUE) {
            size_t k = bi->kind == BK_E ? gt_add_e(&T) : gt_add_opaque(&T, bi->b);
            mpq_sub(vk.q[k], vk.q[k], bi->r);
        } else if (!qv_add_base(&vk, &T, bi, mone)) {
            ok = false;
        }
    }
    if (gR >= 0) mpq_sub(vk.q[gR], vk.q[gR], one);
    if (T.ng != ng) ok = false;   /* a generator appeared late: cannot happen */

    /* ---- real branch admissibility ----------------------------------- */
    bool bases_pos = true;
    for (size_t i = 0; i < F.n && bases_pos; i++) {
        const BaseInfo* bi = &F.f[i].base;
        if (bi->kind == BK_OPAQUE && const_sign(bi->b) != 1) bases_pos = false;
    }
    ExplogRhs rhs = EXPLOG_RHS_POSITIVE;
    Expr* rcond = NULL;
    if (gR >= 0) {
        int s = const_sign(CD.R);
        if (s == 1) rhs = EXPLOG_RHS_POSITIVE;
        else if (s == -1 || s == 2) rhs = EXPLOG_RHS_NONE;
        else {
            rhs = EXPLOG_RHS_COND;
            rcond = ev(fn2("Greater", expr_copy(CD.R), expr_new_integer(0)));
            if (is_sym(rcond, SYM_True)) { expr_free(rcond); rcond = NULL; rhs = EXPLOG_RHS_POSITIVE; }
            else if (is_sym(rcond, SYM_False)) { expr_free(rcond); rcond = NULL; rhs = EXPLOG_RHS_NONE; }
        }
    }

    /* ---- atoms: the rational path ------------------------------------ */
    AtomList A = { NULL, 0, 0 };
    al_index(&A, NULL);                              /* atom 0: the constant */
    AtomTerms* at = (AtomTerms*)calloc(F.n ? F.n : 1, sizeof(AtomTerms));
    bool rational_path = ok;
    for (size_t i = 0; i < F.n && rational_path; i++)
        if (!atom_terms(F.f[i].g, var, &A, &at[i])) rational_path = false;

    if (ok && rational_path) {
        size_t na = A.n;
        /* M[k][m], row-major. */
        mpq_t* M = (mpq_t*)malloc(sizeof(mpq_t) * ng * na);
        for (size_t z = 0; z < ng * na; z++) mpq_init(M[z]);
        mpq_t t;
        mpq_init(t);
        for (size_t i = 0; i < F.n; i++)
            for (size_t j = 0; j < at[i].n; j++)
                for (size_t k = 0; k < ng; k++) {
                    if (mpq_sgn(vb[i].q[k]) == 0) continue;
                    mpq_mul(t, vb[i].q[k], at[i].c[j]);
                    mpq_add(M[k * na + at[i].idx[j]], M[k * na + at[i].idx[j]], t);
                }
        for (size_t k = 0; k < ng; k++) mpq_add(M[k * na + 0], M[k * na + 0], vk.q[k]);

        bool poly_atoms = true, real_atoms = true, has_var_atom = false;
        for (size_t m = 1; m < na; m++) {
            if (atom_degree(A.a[m], var) < 1) poly_atoms = false;
            if (!atom_real(A.a[m], var)) real_atoms = false;
        }

        /* Rank-1 test. */
        long r0 = -1;
        for (size_t k = 0; k < ng && r0 < 0; k++)
            for (size_t m = 0; m < na; m++)
                if (mpq_sgn(M[k * na + m]) != 0) { r0 = (long)k; break; }
        if (r0 < 0) {
            /* Identically satisfied (2^x 3^x == 6^x): left to the evaluator. */
            ok = false;
        }
        bool rank1 = ok;
        size_t m0 = 0;
        QVec w;
        qv_init(&w, ng);
        mpq_t* P = (mpq_t*)malloc(sizeof(mpq_t) * na);
        for (size_t m = 0; m < na; m++) mpq_init(P[m]);
        if (ok) {
            while (mpq_sgn(M[r0 * na + m0]) == 0) m0++;
            for (size_t k = 0; k < ng && rank1; k++) {
                mpq_div(w.q[k], M[k * na + m0], M[r0 * na + m0]);
                for (size_t m = 0; m < na; m++) {
                    mpq_mul(t, w.q[k], M[r0 * na + m]);
                    if (!mpq_equal(t, M[k * na + m])) { rank1 = false; break; }
                }
            }
            for (size_t m = 0; m < na; m++) mpq_set(P[m], M[r0 * na + m]);
            for (size_t m = 1; m < na; m++) if (mpq_sgn(P[m]) != 0) has_var_atom = true;
        }

        out->real_ok = bases_pos && real_atoms;
        out->rhs_real = rhs;
        out->real_cond = rcond; rcond = NULL;

        if (ok && rank1) {
            /* w -> primitive integer vector, P scaled inversely. */
            mpz_t L, G, tmp;
            mpz_inits(L, G, tmp, NULL);
            mpz_set_ui(L, 1);
            mpz_set_ui(G, 0);
            for (size_t k = 0; k < ng; k++) {
                if (mpq_sgn(w.q[k]) == 0) continue;
                mpz_lcm(L, L, mpq_denref(w.q[k]));
            }
            for (size_t k = 0; k < ng; k++) {
                if (mpq_sgn(w.q[k]) == 0) continue;
                mpz_divexact(tmp, L, mpq_denref(w.q[k]));
                mpz_mul(tmp, tmp, mpq_numref(w.q[k]));
                mpz_gcd(G, G, tmp);
            }
            mpq_t s;
            mpq_init(s);
            mpz_set(mpq_numref(s), L);
            mpz_set(mpq_denref(s), G);
            mpq_canonicalize(s);                       /* s = L / G */
            for (size_t k = 0; k < ng; k++) mpq_mul(w.q[k], w.q[k], s);
            for (size_t m = 0; m < na; m++) mpq_div(P[m], P[m], s);

            /* Orientation: prod p^w > 1 when w lies on the rational basis,
             * else the first nonzero component positive. */
            int flip = 0;
            bool rat_only = true;
            for (size_t k = T.B.n; k < ng; k++) if (mpq_sgn(w.q[k]) != 0) rat_only = false;
            if (rat_only) {
                /* log(prod p^w) = sum w_k log p_k: compare prod over w>0
                 * against prod over w<0 exactly. */
                mpz_t num, den, pw;
                mpz_inits(num, den, pw, NULL);
                mpz_set_ui(num, 1); mpz_set_ui(den, 1);
                for (size_t k = 0; k < T.B.n; k++) {
                    int sg = mpq_sgn(w.q[k]);
                    if (!sg) continue;
                    mpz_abs(pw, mpq_numref(w.q[k]));
                    mpz_pow_ui(pw, T.B.p[k], mpz_get_ui(pw));
                    if (sg > 0) mpz_mul(num, num, pw); else mpz_mul(den, den, pw);
                }
                flip = mpz_cmp(num, den) < 0;
                mpz_clears(num, den, pw, NULL);
            } else {
                for (size_t k = 0; k < ng; k++)
                    if (mpq_sgn(w.q[k]) != 0) { flip = mpq_sgn(w.q[k]) < 0; break; }
            }
            if (flip) {
                for (size_t k = 0; k < ng; k++) mpq_neg(w.q[k], w.q[k]);
                for (size_t m = 0; m < na; m++) mpq_neg(P[m], P[m]);
            }

            if (!has_var_atom) {
                /* const * Log[...] == 2 Pi I C: with independent (real)
                 * logarithms the left side is a nonzero real, never a
                 * nonzero imaginary multiple of 2 Pi. */
                if (qv_independent(&T, &w)) out->kind = EXPLOG_NO_SOLUTION;
                else ok = false;
            } else {
                /* Content of the variable part into the unit: 8^x == 2^201
                 * reads x - 67 == 2 Pi I C/Log[8], not 3x - 201 == .../Log[2]. */
                mpz_set_ui(G, 0);
                for (size_t m = 1; m < na; m++)
                    if (mpq_sgn(P[m]) != 0) mpz_gcd(G, G, mpq_numref(P[m]));
                if (mpz_cmp_ui(G, 1) > 0) {
                    mpq_set_z(s, G);
                    for (size_t k = 0; k < ng; k++) mpq_mul(w.q[k], w.q[k], s);
                    for (size_t m = 0; m < na; m++) mpq_div(P[m], P[m], s);
                }
                /* Leading (highest-degree) coefficient positive: C -> -C. */
                size_t lead = 0;
                int ldeg = -2;
                for (size_t m = 1; m < na; m++) {
                    if (mpq_sgn(P[m]) == 0) continue;
                    int dg = atom_degree(A.a[m], var);
                    if (dg > ldeg) { ldeg = dg; lead = m; }
                }
                if (lead && mpq_sgn(P[lead]) < 0)
                    for (size_t m = 0; m < na; m++) mpq_neg(P[m], P[m]);

                out->kind = EXPLOG_LINEAR;
                out->lhs = atoms_poly(&A, P);
                Expr* unit = log_combo(&T, &w);
                out->period = ev(fn2("Times", two_pi_i(),
                                     fn2("Power", unit, expr_new_integer(-1))));
                if (poly_atoms && qv_independent(&T, &w))
                    out->rational_poly = expr_copy(out->lhs);
            }
            mpq_clear(s);
            mpz_clears(L, G, tmp, NULL);
        } else if (ok) {
            /* rank >= 2: sum_m atom_m (column_m . Log) == 2 Pi I C. */
            QVec col;
            qv_init(&col, ng);
            Expr** terms = (Expr**)malloc(sizeof(Expr*) * na);
            size_t ntm = 0;
            for (size_t m = 0; m < na; m++) {
                for (size_t k = 0; k < ng; k++) mpq_set(col.q[k], M[k * na + m]);
                if (qv_is_zero(&col)) continue;
                Expr* c = log_combo(&T, &col);
                terms[ntm++] = A.a[m] ? fn2("Times", c, expr_copy(A.a[m])) : c;
            }
            qv_clear(&col);
            out->kind = EXPLOG_LINEAR;
            out->lhs = ev(expr_new_function(mk_sym("Plus"), terms, ntm));
            free(terms);
            out->period = two_pi_i();
            /* Rational solutions: common roots of the rows, when every row
             * sits on an independent generator. */
            bool indep = true;
            for (size_t k = T.B.n; k < ng && indep; k++) {
                if (T.gen[k].kind != GK_OPAQUE) continue;
                for (size_t m = 0; m < na; m++)
                    if (mpq_sgn(M[k * na + m]) != 0) { indep = false; break; }
            }
            if (indep && poly_atoms) {
                Expr* gcd = NULL;
                for (size_t k = 0; k < ng; k++) {
                    Expr* row = atoms_poly(&A, &M[k * na]);
                    if (row->type == EXPR_INTEGER && row->data.integer == 0) { expr_free(row); continue; }
                    gcd = gcd ? ev(fn2("PolynomialGCD", gcd, row)) : row;
                }
                out->rational_poly = gcd ? gcd : expr_new_integer(0);
            }
        }

        qv_clear(&w);
        for (size_t m = 0; m < na; m++) mpq_clear(P[m]);
        free(P);
        mpq_clear(t);
        for (size_t z = 0; z < ng * na; z++) mpq_clear(M[z]);
        free(M);
    } else if (ok) {
        /* General path: some exponent has a non-rational coefficient.
         * sum_i s_i g_i Log[b_i] - Log[K] == 2 Pi I C, each Log exact. */
        Expr** terms = (Expr**)malloc(sizeof(Expr*) * (F.n + 1));
        bool real_exps = true;
        for (size_t i = 0; i < F.n; i++) {
            terms[i] = fn2("Times", expr_copy(F.f[i].g), log_combo(&T, &vb[i]));
            if (real_exps && !exponent_real_poly(F.f[i].g, var)) real_exps = false;
        }
        terms[F.n] = log_combo(&T, &vk);
        out->kind = EXPLOG_LINEAR;
        out->lhs = ev(expr_new_function(mk_sym("Plus"), terms, F.n + 1));
        free(terms);
        out->period = two_pi_i();
        out->real_ok = bases_pos && real_exps;
        out->rhs_real = rhs;
        out->real_cond = rcond; rcond = NULL;
    }

    for (size_t i = 0; i < F.n; i++) {
        if (at[i].idx) atom_terms_free(&at[i]);
        qv_clear(&vb[i]);
    }
    free(at);
    free(vb);
    qv_clear(&vk);
    mpq_clears(one, mone, NULL);
    al_free(&A);
    if (rcond) expr_free(rcond);
    gt_free(&T);
    cd_free(&CD);
    expr_free(K);
    fl_free(&F);

    if (!ok || out->kind == EXPLOG_NO_SHAPE) {
        explog_linear_free(out);
        out->kind = EXPLOG_NO_SHAPE;
        return false;
    }
    return true;
}

/* ================================================================== *
 *  explog_kernel                                                       *
 * ================================================================== */

void explog_kernel_free(ExplogKernel* K) {
    if (!K) return;
    if (K->poly) expr_free(K->poly);
    if (K->u) expr_free(K->u);
    if (K->kernel) expr_free(K->kernel);
    memset(K, 0, sizeof(*K));
}

bool explog_kernel(const Expr* residual, const Expr* var, ExplogKernel* out) {
    memset(out, 0, sizeof(*out));
    if (!head_is(residual, SYM_Plus) || !var || var->type != EXPR_SYMBOL) return false;
    if (has_inexact(residual)) return false;
    size_t nt = residual->data.function.arg_count;

    /* Per var term: its exponential factors and coefficient. */
    FacList* Fs = (FacList*)calloc(nt, sizeof(FacList));
    Expr** cs = (Expr**)calloc(nt, sizeof(Expr*));
    Expr** fr = (Expr**)malloc(sizeof(Expr*) * nt);
    size_t nfr = 0, nvt = 0;
    bool ok = true;
    for (size_t i = 0; i < nt && ok; i++) {
        const Expr* t = residual->data.function.args[i];
        if (!contains_var(t, var)) { fr[nfr++] = expr_copy((Expr*)t); continue; }
        if (!split_term(t, var, 1, &Fs[i], &cs[i]) || Fs[i].n == 0) { ok = false; break; }
        for (size_t j = 0; j < Fs[i].n; j++) {
            BaseKind bk = Fs[i].f[j].base.kind;
            if (bk != BK_RAT && bk != BK_RATPOW) ok = false;
        }
        nvt++;
    }
    if (ok && (nvt < 2 || (nvt == 2 && nfr == 0))) ok = false;   /* linearize's job */

    ExplogBasis B = { 0, NULL };
    if (ok) {
        size_t nvals = 0, cap = 1;
        for (size_t i = 0; i < nt; i++) cap += 2 * Fs[i].n;
        mpz_t* vals = (mpz_t*)malloc(sizeof(mpz_t) * cap);
        for (size_t i = 0; i < nt; i++)
            for (size_t j = 0; j < Fs[i].n; j++) {
                mpz_init_set(vals[nvals++], Fs[i].f[j].base.n);
                mpz_init_set(vals[nvals++], Fs[i].f[j].base.d);
            }
        explog_basis_build((const mpz_t*)vals, nvals, &B);
        for (size_t i = 0; i < nvals; i++) mpz_clear(vals[i]);
        free(vals);
    }

    /* alpha_t (slope vector) and the coefficient c_t prod b^(const part). */
    size_t nb = B.n;
    mpq_t* alpha = (mpq_t*)malloc(sizeof(mpq_t) * (nt * (nb ? nb : 1)));
    for (size_t z = 0; z < nt * (nb ? nb : 1); z++) mpq_init(alpha[z]);
    Expr** coef = (Expr**)calloc(nt, sizeof(Expr*));
    GenTable T;
    memset(&T, 0, sizeof(T));
    T.B = B;
    T.ng = nb;
    T.ge = -1;
    for (size_t i = 0; i < nt && ok; i++) {
        if (!cs[i]) continue;
        Expr* c = expr_copy(cs[i]);
        for (size_t j = 0; j < Fs[i].n && ok; j++) {
            Expr* x = ev(fn1("Expand", expr_copy(Fs[i].f[j].g)));
            Expr* a = ev(fn2("Coefficient", expr_copy(x), expr_copy((Expr*)var)));
            Expr* b0 = ev(fn2("ReplaceAll", expr_copy(x),
                              fn2("Rule", expr_copy((Expr*)var), expr_new_integer(0))));
            Expr* chk = ev(fn2("Plus", expr_copy(x),
                fn2("Times", expr_new_integer(-1),
                    fn2("Plus", fn2("Times", expr_copy(a), expr_copy((Expr*)var)), expr_copy(b0)))));
            mpq_t aq;
            mpq_init(aq);
            bool lin = chk->type == EXPR_INTEGER && chk->data.integer == 0
                    && rat_mpq(a, aq);
            if (lin) {
                QVec tmp;
                qv_init(&tmp, nb);
                lin = qv_add_base(&tmp, &T, &Fs[i].f[j].base, aq);
                for (size_t k = 0; k < nb; k++)
                    mpq_add(alpha[i * nb + k], alpha[i * nb + k], tmp.q[k]);
                qv_clear(&tmp);
                c = ev(fn2("Times", c, fn2("Power", expr_copy((Expr*)Fs[i].f[j].base.b), expr_copy(b0))));
            }
            if (!lin) ok = false;
            mpq_clear(aq);
            expr_free(x); expr_free(a); expr_free(b0); expr_free(chk);
        }
        coef[i] = c;
    }
    T.B.n = 0; T.B.p = NULL;   /* B is owned below, not by T */

    /* All slopes on one ray: alpha_t = lambda_t w. */
    long r0 = -1;
    for (size_t i = 0; i < nt && ok && r0 < 0; i++) {
        if (!cs[i]) continue;
        for (size_t k = 0; k < nb; k++) if (mpq_sgn(alpha[i * nb + k]) != 0) { r0 = (long)i; break; }
        if (r0 < 0) ok = false;            /* a var term with zero slope */
    }
    mpq_t* lam = (mpq_t*)malloc(sizeof(mpq_t) * nt);
    for (size_t i = 0; i < nt; i++) mpq_init(lam[i]);
    size_t k0 = 0;
    if (ok && r0 >= 0) {
        while (mpq_sgn(alpha[r0 * nb + k0]) == 0) k0++;
        mpq_t t;
        mpq_init(t);
        for (size_t i = 0; i < nt && ok; i++) {
            if (!cs[i]) continue;
            mpq_div(lam[i], alpha[i * nb + k0], alpha[r0 * nb + k0]);
            for (size_t k = 0; k < nb; k++) {
                mpq_mul(t, lam[i], alpha[r0 * nb + k]);
                if (!mpq_equal(t, alpha[i * nb + k])) { ok = false; break; }
            }
            if (mpq_sgn(lam[i]) == 0) ok = false;
        }
        mpq_clear(t);
    } else ok = false;

    if (ok) {
        /* Base B0 = prod p^w with w = alpha_r0 made a primitive integer
         * vector (lambda rescaled), oriented so that B0 > 1. */
        mpz_t L, G, tmp, num, den, pw;
        mpz_inits(L, G, tmp, num, den, pw, NULL);
        mpz_set_ui(L, 1); mpz_set_ui(G, 0);
        for (size_t k = 0; k < nb; k++)
            if (mpq_sgn(alpha[r0 * nb + k])) mpz_lcm(L, L, mpq_denref(alpha[r0 * nb + k]));
        mpq_t* w = (mpq_t*)malloc(sizeof(mpq_t) * (nb ? nb : 1));
        for (size_t k = 0; k < nb; k++) {
            mpq_init(w[k]);
            mpq_set(w[k], alpha[r0 * nb + k]);
            if (!mpq_sgn(w[k])) continue;
            mpz_divexact(tmp, L, mpq_denref(w[k]));
            mpz_mul(tmp, tmp, mpq_numref(w[k]));
            mpz_gcd(G, G, tmp);
        }
        mpq_t s;
        mpq_init(s);
        mpz_set(mpq_numref(s), L); mpz_set(mpq_denref(s), G);
        mpq_canonicalize(s);
        for (size_t k = 0; k < nb; k++) mpq_mul(w[k], w[k], s);
        for (size_t i = 0; i < nt; i++) mpq_div(lam[i], lam[i], s);
        mpz_set_ui(num, 1); mpz_set_ui(den, 1);
        for (size_t k = 0; k < nb; k++) {
            int sg = mpq_sgn(w[k]);
            if (!sg) continue;
            mpz_abs(pw, mpq_numref(w[k]));
            mpz_pow_ui(pw, B.p[k], mpz_get_ui(pw));
            if (sg > 0) mpz_mul(num, num, pw); else mpz_mul(den, den, pw);
        }
        if (mpz_cmp(num, den) < 0) {
            mpz_swap(num, den);
            for (size_t i = 0; i < nt; i++) mpq_neg(lam[i], lam[i]);
        }
        /* delta = gcd of the lambdas; n_t = lambda_t / delta are integers. */
        mpz_t gn, ld;
        mpz_inits(gn, ld, NULL);
        mpz_set_ui(ld, 1);
        for (size_t i = 0; i < nt; i++) {
            if (!cs[i]) continue;
            mpz_gcd(gn, gn, mpq_numref(lam[i]));
            mpz_lcm(ld, ld, mpq_denref(lam[i]));
        }
        mpq_t delta;
        mpq_init(delta);
        mpz_set(mpq_numref(delta), gn); mpz_set(mpq_denref(delta), ld);
        mpq_canonicalize(delta);
        int64_t nmin = 0;
        int64_t* ex = (int64_t*)calloc(nt, sizeof(int64_t));
        mpq_t q;
        mpq_init(q);
        for (size_t i = 0; i < nt && ok; i++) {
            if (!cs[i]) continue;
            mpq_div(q, lam[i], delta);
            if (mpz_cmp_ui(mpq_denref(q), 1) != 0 || !mpz_fits_slong_p(mpq_numref(q))) { ok = false; break; }
            ex[i] = (int64_t)mpz_get_si(mpq_numref(q));
            if (ex[i] < -64 || ex[i] > 64) ok = false;   /* keep the polynomial small */
            if (ex[i] < nmin) nmin = ex[i];
        }
        if (ok) {
            Expr* u = mk_sym(intern_symbol("Solve`Private`explogU"));
            Expr** terms = (Expr**)malloc(sizeof(Expr*) * (nt + 1));
            size_t ntm = 0;
            for (size_t i = 0; i < nt; i++) {
                if (!cs[i]) continue;
                terms[ntm++] = fn2("Times", expr_copy(coef[i]),
                    fn2("Power", expr_copy(u), expr_new_integer(ex[i] - nmin)));
            }
            if (nfr) {
                Expr* f = ev(expr_new_function(mk_sym("Plus"), fr, nfr));
                nfr = 0;
                terms[ntm++] = fn2("Times", f,
                    fn2("Power", expr_copy(u), expr_new_integer(-nmin)));
            }
            out->poly = ev(fn1("Expand", expr_new_function(mk_sym("Plus"), terms, ntm)));
            free(terms);
            mpq_t base;
            mpq_init(base);
            mpz_set(mpq_numref(base), num); mpz_set(mpq_denref(base), den);
            mpq_canonicalize(base);
            out->kernel = ev(fn2("Power", mpq_expr(base),
                                 fn2("Times", mpq_expr(delta), expr_copy((Expr*)var))));
            mpq_clear(base);
            out->u = u;
        }
        mpq_clear(q);
        free(ex);
        mpq_clear(delta);
        mpz_clears(gn, ld, NULL);
        mpq_clear(s);
        for (size_t k = 0; k < nb; k++) mpq_clear(w[k]);
        free(w);
        mpz_clears(L, G, tmp, num, den, pw, NULL);
    }

    for (size_t i = 0; i < nt; i++) mpq_clear(lam[i]);
    free(lam);
    for (size_t z = 0; z < nt * (nb ? nb : 1); z++) mpq_clear(alpha[z]);
    free(alpha);
    for (size_t i = 0; i < nt; i++) {
        if (coef[i]) expr_free(coef[i]);
        if (cs[i]) expr_free(cs[i]);
        fl_free(&Fs[i]);
    }
    free(coef); free(cs); free(Fs);
    for (size_t i = 0; i < nfr; i++) expr_free(fr[i]);
    free(fr);
    free(T.gen);
    explog_basis_clear(&B);

    if (!ok || !out->poly) { explog_kernel_free(out); return false; }
    return true;
}
