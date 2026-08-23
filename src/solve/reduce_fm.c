/*
 * reduce_fm.c
 *
 * Fourier-Motzkin elimination for linear real systems (REDUCE_PLAN.md, Phase 3).
 * See reduce_fm.h.  Coefficient vectors are exact rational Exprs; all arithmetic
 * reuses the evaluator (Plus/Times/Coefficient) so no separate rational library
 * plumbing is needed.
 */
#include "reduce_fm.h"

#include "eval.h"
#include "sym_names.h"

#include <stdlib.h>

#define FM_LE 0
#define FM_LT 1
#define FM_MAX_CON 4000     /* elimination blow-up guard */

/* A constraint  a[0] + a[1] x_0 + ... + a[nv] x_{nv-1}  (rel)  0,
 * where rel is FM_LE or FM_LT and each a[i] is an exact rational Expr. */
typedef struct { Expr** a; int rel; } FMCon;
typedef struct { FMCon* c; int n, cap; } FMSet;

/* ------------------------------------------------------------------ *
 *  Rational-scalar arithmetic (reusing the evaluator)                 *
 * ------------------------------------------------------------------ */

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

static bool is_num_scalar(const Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_REAL:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_FUNCTION: return is_head(e, SYM_Rational);
        default: return false;
    }
}

static int q_sign(const Expr* e) {
    if (e->type == EXPR_INTEGER) return (e->data.integer > 0) - (e->data.integer < 0);
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint);
    if (e->type == EXPR_REAL)    return (e->data.real > 0) - (e->data.real < 0);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    return mpfr_sgn(e->data.mpfr);
#endif
    if (is_head(e, SYM_Rational) && e->data.function.arg_count == 2) {
        const Expr* nu = e->data.function.args[0];
        if (nu->type == EXPR_INTEGER) return (nu->data.integer > 0) - (nu->data.integer < 0);
        if (nu->type == EXPR_BIGINT)  return mpz_sgn(nu->data.bigint);
    }
    return 0;
}

static bool q_is_zero(const Expr* e) {
    return e->type == EXPR_INTEGER && e->data.integer == 0;
}

static Expr* q_add(const Expr* a, const Expr* b) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy((Expr*)a), expr_copy((Expr*)b) }, 2));
}
static Expr* q_mul(const Expr* a, const Expr* b) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_copy((Expr*)a), expr_copy((Expr*)b) }, 2));
}
static Expr* q_neg(const Expr* a) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)a) }, 2));
}

/* ------------------------------------------------------------------ *
 *  Constraint sets                                                    *
 * ------------------------------------------------------------------ */

static void fmcon_free(FMCon* c, int nv) {
    if (!c->a) return;
    for (int i = 0; i <= nv; i++) if (c->a[i]) expr_free(c->a[i]);
    free(c->a); c->a = NULL;
}
static FMCon fmcon_copy(const FMCon* c, int nv) {
    FMCon r; r.rel = c->rel;
    r.a = malloc((size_t)(nv + 1) * sizeof(Expr*));
    for (int i = 0; i <= nv; i++) r.a[i] = expr_copy(c->a[i]);
    return r;
}
static FMSet* fmset_new(void) { return calloc(1, sizeof(FMSet)); }
static void fmset_free(FMSet* s, int nv) {
    if (!s) return;
    for (int i = 0; i < s->n; i++) fmcon_free(&s->c[i], nv);
    free(s->c); free(s);
}
/* Structural equality of two constraints (same rel and identical vectors). */
static bool fmcon_eq(const FMCon* x, const FMCon* y, int nv) {
    if (x->rel != y->rel) return false;
    for (int i = 0; i <= nv; i++) if (!expr_eq(x->a[i], y->a[i])) return false;
    return true;
}
/* Push, taking ownership of `c.a`; drops a structural duplicate. */
static void fmset_push(FMSet* s, FMCon c, int nv) {
    for (int i = 0; i < s->n; i++)
        if (fmcon_eq(&s->c[i], &c, nv)) { fmcon_free(&c, nv); return; }
    if (s->n == s->cap) { s->cap = s->cap ? s->cap * 2 : 8;
                          s->c = realloc(s->c, (size_t)s->cap * sizeof(FMCon)); }
    s->c[s->n++] = c;
}

/* ------------------------------------------------------------------ *
 *  Parse one conjunction into a constraint set                        *
 * ------------------------------------------------------------------ */

/* Coefficient of vars[j] in poly, as an exact rational; *ok=false if symbolic. */
static Expr* lin_coeff(const Expr* poly, const Expr* var, bool* ok) {
    Expr* c = eval_and_free(expr_new_function(expr_new_symbol(SYM_Coefficient),
        (Expr*[]){ expr_copy((Expr*)poly), expr_copy((Expr*)var) }, 2));
    if (!is_num_scalar(c)) { *ok = false; }
    return c;
}
/* Constant term: poly with every var set to 0. */
static Expr* lin_const(const Expr* poly, Expr** vars, int nv, bool* ok) {
    Expr** rr = malloc((size_t)nv * sizeof(Expr*));
    for (int j = 0; j < nv; j++)
        rr[j] = expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_copy(vars[j]), expr_new_integer(0) }, 2);
    Expr* rlist = expr_new_function(expr_new_symbol(SYM_List), rr, (size_t)nv);
    free(rr);
    Expr* c = eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll),
        (Expr*[]){ expr_copy((Expr*)poly), rlist }, 2));
    if (!is_num_scalar(c)) { *ok = false; }
    return c;
}

/* Build the coefficient vector of an atom's poly; verify it is truly linear
 * (poly == const + sum coeff*var).  Returns NULL on decline. */
static Expr** atom_vector(const Expr* poly, Expr** vars, int nv, bool* ok) {
    if (nv > 15) { *ok = false; return NULL; }   /* lin_const cap */
    Expr** a = malloc((size_t)(nv + 1) * sizeof(Expr*));
    a[0] = lin_const(poly, vars, nv, ok);
    for (int j = 0; j < nv; j++) a[j + 1] = lin_coeff(poly, vars[j], ok);
    if (!*ok) { for (int i = 0; i <= nv; i++) expr_free(a[i]); free(a); return NULL; }
    /* linearity check: poly - (a0 + sum a_{j+1} vars[j]) == 0 */
    Expr* recon = expr_copy(a[0]);
    for (int j = 0; j < nv; j++) {
        Expr* term = q_mul(a[j + 1], vars[j]);
        recon = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ recon, term }, 2));
    }
    Expr* diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Subtract),
        (Expr*[]){ expr_copy((Expr*)poly), recon }, 2));
    if (!q_is_zero(diff)) { *ok = false; }
    expr_free(diff);
    if (!*ok) { for (int i = 0; i <= nv; i++) expr_free(a[i]); free(a); return NULL; }
    return a;
}

static FMSet* parse_conjunction(const RConj* conj, Expr** vars, int nv, bool* ok) {
    FMSet* S = fmset_new();
    for (int k = 0; k < conj->n; k++) {
        const RAtom* at = &conj->a[k];
        if (at->rel == R_ELEM || at->rel == R_NE || at->nonconst_denom) { *ok = false; break; }
        Expr** a = atom_vector(at->poly, vars, nv, ok);
        if (!*ok) break;
        if (at->rel == R_EQ) {
            /* poly == 0  ==>  poly <= 0 && -poly <= 0 */
            FMCon lo; lo.rel = FM_LE; lo.a = a;
            FMCon hi; hi.rel = FM_LE;
            hi.a = malloc((size_t)(nv + 1) * sizeof(Expr*));
            for (int i = 0; i <= nv; i++) hi.a[i] = q_neg(a[i]);
            fmset_push(S, lo, nv);
            fmset_push(S, hi, nv);
        } else {
            FMCon c; c.rel = (at->rel == R_LT) ? FM_LT : FM_LE; c.a = a;
            fmset_push(S, c, nv);
        }
    }
    if (!*ok) { fmset_free(S, nv); return NULL; }
    return S;
}

/* ------------------------------------------------------------------ *
 *  Elimination                                                        *
 * ------------------------------------------------------------------ */

/* combined = (-a_q[j+1]) * p + (a_p[j+1]) * q  -- cancels variable j. */
static FMCon combine(const FMCon* p, const FMCon* q, int j, int nv) {
    Expr* mp = q_neg(q->a[j + 1]);     /* > 0 */
    const Expr* mq = p->a[j + 1];      /* > 0 */
    FMCon r; r.rel = (p->rel == FM_LT || q->rel == FM_LT) ? FM_LT : FM_LE;
    r.a = malloc((size_t)(nv + 1) * sizeof(Expr*));
    for (int i = 0; i <= nv; i++) {
        Expr* t1 = q_mul(mp, p->a[i]);
        Expr* t2 = q_mul(mq, q->a[i]);
        r.a[i] = q_add(t1, t2);
        expr_free(t1); expr_free(t2);
    }
    expr_free(mp);
    return r;
}

/* Project variable j out of S.  Returns a new set (or NULL if the blow-up guard
 * trips). */
static FMSet* eliminate(const FMSet* S, int j, int nv, bool* ok) {
    FMSet* out = fmset_new();
    for (int i = 0; i < S->n; i++)                 /* constraints without var j */
        if (q_sign(S->c[i].a[j + 1]) == 0) fmset_push(out, fmcon_copy(&S->c[i], nv), nv);
    for (int p = 0; p < S->n; p++) {
        if (q_sign(S->c[p].a[j + 1]) <= 0) continue;
        for (int q = 0; q < S->n; q++) {
            if (q_sign(S->c[q].a[j + 1]) >= 0) continue;
            fmset_push(out, combine(&S->c[p], &S->c[q], j, nv), nv);
            if (out->n > FM_MAX_CON) { *ok = false; return out; }
        }
    }
    return out;
}

/* ------------------------------------------------------------------ *
 *  Bound expression / emission                                        *
 * ------------------------------------------------------------------ */

/* B = -(a0 + sum_{i<j} a_{i+1} vars[i]) / a_{j+1}, for a level-j constraint. */
static Expr* bound_expr(const FMCon* c, int j, Expr** vars) {
    Expr* R = expr_copy(c->a[0]);
    for (int i = 0; i < j; i++) {
        Expr* term = q_mul(c->a[i + 1], vars[i]);
        R = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ R, term }, 2));
    }
    Expr* negR = q_neg(R); expr_free(R);
    Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy(c->a[j + 1]), expr_new_integer(-1) }, 2));
    Expr* B = q_mul(negR, inv);
    expr_free(negR); expr_free(inv);
    /* Expand so `-(-1 + x)` prints as `1 - x` and equal lower/upper bounds
     * become structurally identical (letting an equation re-emit as `==`). */
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand),
        (Expr*[]){ B }, 1));
}

typedef struct { Expr* b; bool strict; } Bound;

/* Emit the atom(s) constraining vars[j], given the level-(j+1) projection. */
static Expr* emit_level(const FMSet* proj_hi, int j, Expr** vars) {
    Bound* los = NULL; int nlo = 0;
    Bound* his = NULL; int nhi = 0;
    for (int i = 0; i < proj_hi->n; i++) {
        const FMCon* c = &proj_hi->c[i];
        int s = q_sign(c->a[j + 1]);
        if (s == 0) continue;
        Expr* B = bound_expr(c, j, vars);
        bool strict = (c->rel == FM_LT);
        Bound** dst = (s > 0) ? &his : &los;
        int* n = (s > 0) ? &nhi : &nlo;
        bool dup = false;
        for (int t = 0; t < *n; t++)
            if ((*dst)[t].strict == strict && expr_eq((*dst)[t].b, B)) { dup = true; break; }
        if (dup) { expr_free(B); continue; }
        *dst = realloc(*dst, (size_t)(*n + 1) * sizeof(Bound));
        (*dst)[*n].b = B; (*dst)[*n].strict = strict; (*n)++;
    }

    Expr* result = NULL;   /* NULL == no constraint on this (free) variable */

    /* Any non-strict lower bound equal to a non-strict upper bound pins the
     * variable to an equation, even when there are further bounds: those are
     * redundant, because Fourier-Motzkin projection has already emitted their
     * cross-implications (lo' <= hi') as constraints on the earlier variables
     * (that is exactly how `x + y == 1 && x - y == 3` yields x == 2 at the
     * x-level, leaving the y-level free to collapse to y == 1 - x). */
    const Expr* eqval = NULL;
    for (int t = 0; t < nlo && !eqval; t++) {
        if (los[t].strict) continue;
        for (int u = 0; u < nhi; u++) {
            if (his[u].strict) continue;
            if (expr_eq(los[t].b, his[u].b)) { eqval = los[t].b; break; }
        }
    }
    if (eqval) {
        result = expr_new_function(expr_new_symbol(SYM_Equal),
            (Expr*[]){ expr_copy(vars[j]), expr_copy((Expr*)eqval) }, 2);
    } else if (nlo == 1 && nhi == 1) {
        /* chained  lo (op) x (op) hi */
        const char* op1 = los[0].strict ? SYM_Less : SYM_LessEqual;
        const char* op2 = his[0].strict ? SYM_Less : SYM_LessEqual;
        result = expr_new_function(expr_new_symbol(SYM_Inequality),
            (Expr*[]){ expr_copy(los[0].b), expr_new_symbol(op1),
                       expr_copy(vars[j]),  expr_new_symbol(op2),
                       expr_copy(his[0].b) }, 5);
    } else {
        /* general: AND of each lower/upper bound atom */
        Expr** parts = NULL; int npar = 0;
        for (int t = 0; t < nlo; t++) {
            const char* op = los[t].strict ? SYM_Greater : SYM_GreaterEqual;
            Expr* atom = expr_new_function(expr_new_symbol(op),
                (Expr*[]){ expr_copy(vars[j]), expr_copy(los[t].b) }, 2);
            parts = realloc(parts, (size_t)(npar + 1) * sizeof(Expr*)); parts[npar++] = atom;
        }
        for (int t = 0; t < nhi; t++) {
            const char* op = his[t].strict ? SYM_Less : SYM_LessEqual;
            Expr* atom = expr_new_function(expr_new_symbol(op),
                (Expr*[]){ expr_copy(vars[j]), expr_copy(his[t].b) }, 2);
            parts = realloc(parts, (size_t)(npar + 1) * sizeof(Expr*)); parts[npar++] = atom;
        }
        if (npar == 1) { result = parts[0]; }
        else if (npar > 1) {
            result = expr_new_function(expr_new_symbol(SYM_And), parts, (size_t)npar);
        }
        free(parts);
    }

    for (int t = 0; t < nlo; t++) expr_free(los[t].b);
    for (int t = 0; t < nhi; t++) expr_free(his[t].b);
    free(los); free(his);
    return result;
}

/* ------------------------------------------------------------------ *
 *  One conjunction -> formula                                         *
 * ------------------------------------------------------------------ */

static Expr* fm_conjunction(const RConj* conj, Expr** vars, int nv) {
    bool ok = true;
    FMSet* full = parse_conjunction(conj, vars, nv, &ok);
    if (!ok) return NULL;

    /* proj[k] = projection onto vars[0..k-1]; proj[nv] = full. */
    FMSet** proj = calloc((size_t)(nv + 1), sizeof(FMSet*));
    proj[nv] = full;
    for (int j = nv - 1; j >= 0; j--) {
        proj[j] = eliminate(proj[j + 1], j, nv, &ok);
        if (!ok) break;
    }

    Expr* result = NULL;
    if (ok) {
        /* Feasibility: proj[0] holds only constant constraints. */
        bool infeasible = false;
        for (int i = 0; i < proj[0]->n && !infeasible; i++) {
            const FMCon* c = &proj[0]->c[i];
            int s = q_sign(c->a[0]);
            if (c->rel == FM_LE && s > 0) infeasible = true;   /* c0 <= 0 with c0 > 0 */
            if (c->rel == FM_LT && s >= 0) infeasible = true;  /* c0 < 0  with c0 >= 0 */
        }
        if (infeasible) {
            result = expr_new_symbol(SYM_False);
        } else {
            Expr** parts = NULL; int npar = 0;
            for (int j = 0; j < nv; j++) {
                Expr* lv = emit_level(proj[j + 1], j, vars);
                if (lv) { parts = realloc(parts, (size_t)(npar + 1) * sizeof(Expr*));
                          parts[npar++] = lv; }
            }

            /* Back-substitution: forward-propagate any level that pins its
             * variable to a numeric constant (x == 2) into the later levels, so
             * a fully-determined system reads x == 2 && y == -1 rather than the
             * triangular x == 2 && y == 1 - x.  Bounds/ranges are untouched
             * (their RHS is non-numeric), matching x > 0 && y == 1 - x. */
            {
                Expr** subrules = NULL; int nsub = 0;
                for (int t = 0; t < npar; t++) {
                    if (nsub > 0) {
                        Expr** rc = malloc((size_t)nsub * sizeof(Expr*));
                        for (int s = 0; s < nsub; s++) rc[s] = expr_copy(subrules[s]);
                        Expr* rlist = expr_new_function(expr_new_symbol(SYM_List),
                                                        rc, (size_t)nsub);
                        free(rc);
                        parts[t] = eval_and_free(expr_new_function(
                            expr_new_symbol(SYM_ReplaceAll),
                            (Expr*[]){ parts[t], rlist }, 2));
                    }
                    if (is_head(parts[t], SYM_Equal)
                        && parts[t]->data.function.arg_count == 2
                        && is_num_scalar(parts[t]->data.function.args[1])) {
                        subrules = realloc(subrules,
                                           (size_t)(nsub + 1) * sizeof(Expr*));
                        subrules[nsub++] = expr_new_function(expr_new_symbol(SYM_Rule),
                            (Expr*[]){ expr_copy(parts[t]->data.function.args[0]),
                                       expr_copy(parts[t]->data.function.args[1]) }, 2);
                    }
                }
                for (int s = 0; s < nsub; s++) expr_free(subrules[s]);
                free(subrules);
            }

            if (npar == 0) result = expr_new_symbol(SYM_True);
            else if (npar == 1) result = parts[0];
            else result = expr_new_function(expr_new_symbol(SYM_And), parts, (size_t)npar);
            free(parts);
        }
    }

    for (int j = 0; j <= nv; j++) if (proj[j]) fmset_free(proj[j], nv);
    free(proj);
    return result;   /* NULL if !ok */
}

/* ------------------------------------------------------------------ *
 *  Driver                                                             *
 * ------------------------------------------------------------------ */

Expr* reduce_fm(const RForm* F, Expr** vars, int nv) {
    if (F->is_true) return expr_new_symbol(SYM_True);
    if (F->n == 0)  return expr_new_symbol(SYM_False);

    Expr** disj = NULL; int nd = 0;
    for (int i = 0; i < F->n; i++) {
        Expr* sub = fm_conjunction(F->c[i], vars, nv);
        if (!sub) {                                   /* decline whole formula */
            for (int t = 0; t < nd; t++) expr_free(disj[t]);
            free(disj);
            return NULL;
        }
        disj = realloc(disj, (size_t)(nd + 1) * sizeof(Expr*));
        disj[nd++] = sub;
    }

    Expr* out;
    if (nd == 1) out = disj[0];
    else out = expr_new_function(expr_new_symbol(SYM_Or), disj, (size_t)nd);
    free(disj);
    return eval_and_free(out);
}
