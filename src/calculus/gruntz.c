/*
 * gruntz.c -- Gruntz's mrv algorithm for computing limits of exp-log
 * functions. Faithful C port of the method in Dominik Gruntz's 1996 ETH PhD
 * thesis (Chapter 3 + Appendix A), following the proven structure of SymPy's
 * gruntz.py.
 *
 * The engine works at x -> +oo. gruntz_limit() reduces every other request
 * (finite point, -oo, one-sided) to that case by a substitution.
 *
 * Key pieces (see the thesis / SymPy for the theory):
 *   mrv(e,x)         -- the most-rapidly-varying subexpression set Omega,
 *                       plus e rewritten in terms of fresh dummies for each
 *                       Omega element (Algorithm 3.12 + SubsSet bookkeeping).
 *   compare(a,b,x)   -- comparability-class order via limitinf(ln a/ln b).
 *   rewrite(...)     -- express the mrv elements as A*w^c, w->0+ (rule 3.2).
 *   mrv_leadterm     -- leading (coeff, exponent) of the series in w.
 *   limitinf(e,x)    -- lim_{x->oo} e, driving the leading-term recursion.
 *   sign(e,x)        -- sign of e as x->oo (an oracle, per Gruntz).
 *
 * Memory: inputs to gruntz_limit are borrowed. Internally every helper
 * returns a freshly-owned Expr* (or a SubsSet the caller must ss_free). On
 * any undecidable oracle / cap / series failure a file-static g_fail flag is
 * raised; the driver then returns NULL so Limit is left unevaluated.
 */

#include "gruntz.h"
#include "expr.h"
#include "eval.h"
#include "common.h"        /* head_is */
#include "sym_names.h"
#include "zero_test.h"
#include "numeric.h"
#include "arithmetic.h"    /* expr_numeric_sign */

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Direction codes must match limit.c's LIMIT_DIR_* for the finite-point
 * reduction. Only the sign matters here. */
#define GZ_DIR_TWOSIDED   0
#define GZ_DIR_FROMABOVE  1
#define GZ_DIR_FROMBELOW (-1)

#define GRUNTZ_MAX_DEPTH  80
#define GRUNTZ_SERIES_START 2
#define GRUNTZ_SERIES_MAX   16

/* ---------------------------------------------------------------------- */
/* Failure / depth state (poor man's exception).                          */
/* ---------------------------------------------------------------------- */
static int  g_fail;
static int  g_depth;
static long g_work;              /* total mrv_leadterm invocations this call  */
#define GRUNTZ_MAX_WORK 6000     /* hard budget: abort (unevaluated) if hit    */

#define GZFAIL() do { g_fail = 1; } while (0)


/* ---------------------------------------------------------------------- */
/* Tiny builders                                                           */
/* ---------------------------------------------------------------------- */
static Expr* mk_int(int64_t v)     { return expr_new_integer(v); }
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }

static Expr* mk_fn1(const char* n, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(mk_sym(n), args, 1);
}
static Expr* mk_fn2(const char* n, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(mk_sym(n), args, 2);
}
static Expr* mk_times(Expr* a, Expr* b) { return mk_fn2("Times", a, b); }
static Expr* mk_plus(Expr* a, Expr* b)  { return mk_fn2("Plus", a, b); }
static Expr* mk_pow(Expr* a, Expr* b)   { return mk_fn2("Power", a, b); }
static Expr* mk_neg(Expr* a)            { return mk_times(mk_int(-1), a); }
static Expr* mk_exp(Expr* a)            { return mk_pow(mk_sym("E"), a); }
static Expr* mk_log(Expr* a)            { return mk_fn1("Log", a); }

/* Evaluate + free the source. evaluate() copies its input internally. */
static Expr* simp(Expr* e) {
    if (!e) return NULL;
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* ---------------------------------------------------------------------- */
/* Memoization caches. Without them the mutual recursion re-derives the     */
/* same sub-limits exponentially (SymPy caches limitinf/sign/mrv_leadterm). */
/* Keyed on the expression alone: the limit variable x is fixed per call.   */
/* All caches are cleared at the start of each gruntz_limit().              */
/* ---------------------------------------------------------------------- */
#define GZ_CACHE_BUCKETS 8192
typedef struct CNode {
    Expr* key;
    Expr* val;   /* limitinf result / mrv_leadterm c0        */
    Expr* val2;  /* mrv_leadterm e0                          */
    int   ival;  /* sign result                             */
    struct CNode* next;
} CNode;

static CNode* g_cache_limit[GZ_CACHE_BUCKETS];
static CNode* g_cache_sign[GZ_CACHE_BUCKETS];
static CNode* g_cache_lt[GZ_CACHE_BUCKETS];

static void cache_clear(CNode** tbl) {
    for (int i = 0; i < GZ_CACHE_BUCKETS; i++) {
        CNode* n = tbl[i];
        while (n) {
            CNode* nx = n->next;
            if (n->key)  expr_free(n->key);
            if (n->val)  expr_free(n->val);
            if (n->val2) expr_free(n->val2);
            free(n);
            n = nx;
        }
        tbl[i] = NULL;
    }
}
static void cache_clear_all(void) {
    cache_clear(g_cache_limit);
    cache_clear(g_cache_sign);
    cache_clear(g_cache_lt);
}
static CNode* cache_find(CNode** tbl, const Expr* key) {
    unsigned idx = (unsigned)(expr_hash(key) & (GZ_CACHE_BUCKETS - 1));
    for (CNode* n = tbl[idx]; n; n = n->next)
        if (expr_eq(n->key, (Expr*)key)) return n;
    return NULL;
}
static CNode* cache_insert(CNode** tbl, const Expr* key) {
    unsigned idx = (unsigned)(expr_hash(key) & (GZ_CACHE_BUCKETS - 1));
    CNode* n = (CNode*)calloc(1, sizeof(CNode));
    n->key = expr_copy((Expr*)key);
    n->next = tbl[idx];
    tbl[idx] = n;
    return n;
}

/* ---------------------------------------------------------------------- */
/* Predicates                                                              */
/* ---------------------------------------------------------------------- */
static bool is_sym_name(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL &&
           strcmp(e->data.symbol.name, name) == 0;
}
static bool is_E(const Expr* e) { return is_sym_name(e, "E"); }

static bool has_x(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (expr_eq((Expr*)e, (Expr*)x)) return true;
    if (e->type == EXPR_FUNCTION) {
        if (has_x(e->data.function.head, x)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (has_x(e->data.function.args[i], x)) return true;
    }
    return false;
}

static bool is_mul(const Expr* e) { return head_is(e, SYM_Times); }
static bool is_add(const Expr* e) { return head_is(e, SYM_Plus); }
static bool is_pow(const Expr* e) { return head_is(e, SYM_Power); }
static bool is_log(const Expr* e) { return head_is(e, SYM_Log) &&
                                           e->data.function.arg_count == 1; }

/* Exponential detection. Returns the exponent argument (borrowed) if e is
 * Power[E, arg] or Exp[arg]; NULL otherwise. */
static Expr* exp_arg(const Expr* e) {
    if (is_pow(e) && e->data.function.arg_count == 2 &&
        is_E(e->data.function.args[0]))
        return e->data.function.args[1];
    if (head_is(e, SYM_Exp) && e->data.function.arg_count == 1)
        return e->data.function.args[0];
    return NULL;
}
static bool is_exp(const Expr* e) { return exp_arg(e) != NULL; }

static bool is_infinity_val(const Expr* e) {
    if (!e) return false;
    if (is_sym_name(e, "Infinity") || is_sym_name(e, "ComplexInfinity"))
        return true;
    if (head_is(e, SYM_DirectedInfinity)) return true;
    if (is_mul(e)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (is_infinity_val(e->data.function.args[i])) return true;
    }
    return false;
}

static bool is_lit_zero(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL)    return e->data.real == 0.0;
    return false;
}

/* ---------------------------------------------------------------------- */
/* Structural replace (xreplace): replace exact matches of `from` with a    */
/* copy of `to`, without descending into replacements, without evaluating.  */
/* ---------------------------------------------------------------------- */
static Expr* xreplace(const Expr* e, const Expr* from, const Expr* to) {
    if (expr_eq((Expr*)e, (Expr*)from)) return expr_copy((Expr*)to);
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    Expr* head = xreplace(e->data.function.head, from, to);
    size_t n = e->data.function.arg_count;
    Expr** args = (n ? (Expr**)malloc(n * sizeof(Expr*)) : NULL);
    for (size_t i = 0; i < n; i++)
        args[i] = xreplace(e->data.function.args[i], from, to);
    Expr* r = expr_new_function(head, args, n);
    if (args) free(args);
    return r;
}

/* ---------------------------------------------------------------------- */
/* Numeric sign of an x-free expression: +1 / -1 / 0, or 2 = "unknown".    */
/* ---------------------------------------------------------------------- */
static int numeric_sign_of(const Expr* e) {
    NumericSpec spec = numeric_machine_spec();
    Expr* v = numericalize(e, spec);
    int r = 2;
    if (v) {
        if (v->type == EXPR_INTEGER || v->type == EXPR_REAL ||
            v->type == EXPR_BIGINT) {
            int s = expr_numeric_sign(v);
            r = s; /* -1/0/1 */
        }
#ifdef USE_MPFR
        else if (v->type == EXPR_MPFR) {
            /* mpfr sign via the numeric helper is not exposed here; fall
             * back to machine double comparison through re-eval. */
        }
#endif
        expr_free(v);
    }
    /* Distinguish a genuine zero from a near-zero cancellation using the
     * zero-test oracle. */
    if (r == 0) {
        if (zero_test_decide(e) == ZERO_TEST_FALSE) r = 2; /* nonzero, sign? */
    }
    return r;
}

/* ---------------------------------------------------------------------- */
/* SubsSet -- {expr -> dummy} plus rewrites {dummy -> expr}.               */
/* ---------------------------------------------------------------------- */
typedef struct {
    Expr** key;   /* the mrv exprs (owned)          */
    Expr** var;   /* fresh dummy symbols (owned)    */
    size_t n;
    Expr** rw_var;   /* dummy (owned)  */
    Expr** rw_expr;  /* expr  (owned)  */
    size_t rw_n;
} SubsSet;

static long g_dummy_ctr;

static Expr* fresh_dummy(void) {
    char buf[32];
    snprintf(buf, sizeof(buf), "$gz%ld$", ++g_dummy_ctr);
    return mk_sym(buf);
}

static void ss_init(SubsSet* s) { memset(s, 0, sizeof(*s)); }

static void ss_free(SubsSet* s) {
    for (size_t i = 0; i < s->n; i++) { expr_free(s->key[i]); expr_free(s->var[i]); }
    for (size_t i = 0; i < s->rw_n; i++) { expr_free(s->rw_var[i]); expr_free(s->rw_expr[i]); }
    free(s->key); free(s->var); free(s->rw_var); free(s->rw_expr);
    ss_init(s);
}

static bool ss_empty(const SubsSet* s) { return s->n == 0; }

/* Index of key (by structural equality), or -1. */
static int ss_find(const SubsSet* s, const Expr* key) {
    for (size_t i = 0; i < s->n; i++)
        if (expr_eq(s->key[i], (Expr*)key)) return (int)i;
    return -1;
}

/* Get the dummy for `key` (borrowed); create a fresh one if absent. */
static Expr* ss_get(SubsSet* s, const Expr* key) {
    int i = ss_find(s, key);
    if (i >= 0) return s->var[i];
    s->key = (Expr**)realloc(s->key, (s->n + 1) * sizeof(Expr*));
    s->var = (Expr**)realloc(s->var, (s->n + 1) * sizeof(Expr*));
    s->key[s->n] = expr_copy((Expr*)key);
    s->var[s->n] = fresh_dummy();
    s->n++;
    return s->var[s->n - 1];
}

/* Append a rewrite mapping dummy->expr (copies). */
static void ss_add_rewrite(SubsSet* s, const Expr* var, const Expr* expr) {
    s->rw_var  = (Expr**)realloc(s->rw_var,  (s->rw_n + 1) * sizeof(Expr*));
    s->rw_expr = (Expr**)realloc(s->rw_expr, (s->rw_n + 1) * sizeof(Expr*));
    s->rw_var[s->rw_n]  = expr_copy((Expr*)var);
    s->rw_expr[s->rw_n] = expr_copy((Expr*)expr);
    s->rw_n++;
}
static Expr* ss_rewrite_of(const SubsSet* s, const Expr* var) {
    for (size_t i = 0; i < s->rw_n; i++)
        if (expr_eq(s->rw_var[i], (Expr*)var)) return s->rw_expr[i]; /* borrowed */
    return NULL;
}

static void ss_copy(SubsSet* dst, const SubsSet* src) {
    ss_init(dst);
    for (size_t i = 0; i < src->n; i++) { (void)ss_get(dst, src->key[i]);
        /* replace fresh dummy with the source's dummy to keep identity */
        expr_free(dst->var[dst->n - 1]);
        dst->var[dst->n - 1] = expr_copy(src->var[i]);
    }
    for (size_t i = 0; i < src->rw_n; i++) ss_add_rewrite(dst, src->rw_var[i], src->rw_expr[i]);
}

/* do_subs: replace every dummy in `e` by its key expr. Returns owned. */
static Expr* ss_do_subs(const SubsSet* s, const Expr* e) {
    Expr* cur = expr_copy((Expr*)e);
    for (size_t i = 0; i < s->n; i++) {
        Expr* nx = xreplace(cur, s->var[i], s->key[i]);
        expr_free(cur); cur = nx;
    }
    return cur;
}

static bool ss_meets(const SubsSet* a, const SubsSet* b) {
    for (size_t i = 0; i < a->n; i++)
        if (ss_find(b, a->key[i]) >= 0) return true;
    return false;
}

/* union: merge s2 into a copy of s1. If exps_io != NULL, translate s2's
 * dummies (for keys already in s1) to s1's dummies within *exps_io. */
static void ss_union(const SubsSet* s1, const SubsSet* s2,
                     SubsSet* out, Expr** exps_io) {
    ss_copy(out, s1);
    /* translation map for colliding keys: s2->var  ==> s1->var */
    for (size_t i = 0; i < s2->n; i++) {
        int j = ss_find(s1, s2->key[i]);
        if (j >= 0) {
            if (exps_io && *exps_io) {
                Expr* nx = xreplace(*exps_io, s2->var[i], s1->var[j]);
                expr_free(*exps_io); *exps_io = nx;
            }
        } else {
            /* add key with s2's dummy preserved */
            out->key = (Expr**)realloc(out->key, (out->n + 1) * sizeof(Expr*));
            out->var = (Expr**)realloc(out->var, (out->n + 1) * sizeof(Expr*));
            out->key[out->n] = expr_copy(s2->key[i]);
            out->var[out->n] = expr_copy(s2->var[i]);
            out->n++;
        }
    }
    /* carry over s2 rewrites, translating colliding dummies */
    for (size_t i = 0; i < s2->rw_n; i++) {
        Expr* rexpr = expr_copy(s2->rw_expr[i]);
        for (size_t k = 0; k < s2->n; k++) {
            int j = ss_find(s1, s2->key[k]);
            if (j >= 0) {
                Expr* nx = xreplace(rexpr, s2->var[k], s1->var[j]);
                expr_free(rexpr); rexpr = nx;
            }
        }
        ss_add_rewrite(out, s2->rw_var[i], rexpr);
        expr_free(rexpr);
    }
}

/* ---------------------------------------------------------------------- */
/* Forward declarations of the mutually-recursive core.                    */
/* ---------------------------------------------------------------------- */
static void  mrv(const Expr* e, const Expr* x, SubsSet* out, Expr** expr_out);
static char  compare(const Expr* a, const Expr* b, const Expr* x);
static int   gz_sign(const Expr* e, const Expr* x);
static Expr* limitinf(const Expr* e, const Expr* x);
static void  mrv_leadterm(const Expr* e, const Expr* x,
                          Expr** c0_out, Expr** e0_out);
static int   gz_sign_core(const Expr* e, const Expr* x);
static Expr* limitinf_core(const Expr* e, const Expr* x);
static void  mrv_leadterm_core(const Expr* e, const Expr* x,
                               Expr** c0_out, Expr** e0_out);
static Expr* make_exponent(int64_t num, int64_t den);
static bool  contains_complex_head(const Expr* e);
static int64_t gz_leaf_count(const Expr* e);

/* ---------------------------------------------------------------------- */
/* powsimp(combine='exp'): merge exponentials in a product / power tower.   */
/* We rely mostly on the evaluator, but explicitly combine E^a * E^b and    */
/* (E^a)^b so mrv sees a single exponential.                                */
/* ---------------------------------------------------------------------- */
static Expr* exp_combine(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    /* recurse first */
    size_t n = e->data.function.arg_count;
    Expr* head = exp_combine(e->data.function.head);
    Expr** args = (n ? (Expr**)malloc(n * sizeof(Expr*)) : NULL);
    for (size_t i = 0; i < n; i++) args[i] = exp_combine(e->data.function.args[i]);
    Expr* r = expr_new_function(head, args, n);
    if (args) free(args);

    /* (E^a)^b -> E^(a*b) */
    if (is_pow(r) && r->data.function.arg_count == 2) {
        Expr* base = r->data.function.args[0];
        Expr* ex   = r->data.function.args[1];
        Expr* ba   = exp_arg(base);
        if (ba) {
            Expr* merged = mk_exp(simp(mk_times(expr_copy(ba), expr_copy(ex))));
            expr_free(r);
            return merged;
        }
    }
    /* E^a * E^b * ... -> E^(a+b+...) ; keep non-exp factors */
    if (is_mul(r)) {
        Expr* sumarg = NULL;
        Expr* rest = NULL; /* product of non-exp factors */
        for (size_t i = 0; i < r->data.function.arg_count; i++) {
            Expr* f = r->data.function.args[i];
            Expr* a = exp_arg(f);
            if (a) sumarg = sumarg ? simp(mk_plus(sumarg, expr_copy(a))) : expr_copy(a);
            else   rest   = rest   ? simp(mk_times(rest, expr_copy(f)))  : expr_copy(f);
        }
        if (sumarg) {
            Expr* ex = mk_exp(sumarg);
            Expr* out = rest ? simp(mk_times(rest, ex)) : ex;
            expr_free(r);
            return out;
        }
        if (rest) expr_free(rest);
    }
    return r;
}

/* ---------------------------------------------------------------------- */
/* as_two_terms / independent split helpers for Mul & Add.                 */
/* Collect the x-dependent operands (dfac) and the x-free product/sum (i).  */
/* ---------------------------------------------------------------------- */
static Expr* combine_head(const char* head, Expr** parts, size_t n) {
    if (n == 0) return NULL;
    if (n == 1) return parts[0];
    Expr* acc = parts[0];
    for (size_t i = 1; i < n; i++) acc = simp(expr_new_function(mk_sym(head),
                                     (Expr*[]){ acc, parts[i] }, 2));
    return acc;
}

/* ---------------------------------------------------------------------- */
/* mrv_max: choose the higher comparability class (or the union).          */
/* mrv_max3 mirrors SymPy: f,expsf / g,expsg / union,expsboth.             */
/* Takes ownership of expsf,expsg,expsboth (frees the two not returned).    */
/* Copies whichever SubsSet wins into `out`; caller owns f,g,u separately.  */
/* ---------------------------------------------------------------------- */
static void mrv_max3(const SubsSet* f, Expr* expsf,
                     const SubsSet* g, Expr* expsg,
                     const SubsSet* u, Expr* expsboth,
                     const Expr* x, SubsSet* out, Expr** expr_out) {
    if (g_fail) { *expr_out = expsf; expr_free(expsg); expr_free(expsboth);
                  ss_copy(out, f); return; }
    if (ss_empty(f)) { ss_copy(out, g); *expr_out = expsg;
                       expr_free(expsf); expr_free(expsboth); return; }
    if (ss_empty(g)) { ss_copy(out, f); *expr_out = expsf;
                       expr_free(expsg); expr_free(expsboth); return; }
    if (ss_meets(f, g)) { ss_copy(out, u); *expr_out = expsboth;
                          expr_free(expsf); expr_free(expsg); return; }
    char c = compare(f->key[0], g->key[0], x);
    if (c == '>') { ss_copy(out, f); *expr_out = expsf;
                    expr_free(expsg); expr_free(expsboth); }
    else if (c == '<') { ss_copy(out, g); *expr_out = expsg;
                         expr_free(expsf); expr_free(expsboth); }
    else { /* '=' or fail */ ss_copy(out, u); *expr_out = expsboth;
           expr_free(expsf); expr_free(expsg); }
}

/* mrv_max1(s1,e1, s2,e2): combine two mrv results whose combined expression
 * (already built from e1,e2) is `exps`. Takes ownership of exps. */
static void mrv_max1(const SubsSet* s1, const SubsSet* s2, Expr* exps,
                     const Expr* x, SubsSet* out, Expr** expr_out) {
    SubsSet u; ss_init(&u);
    Expr* b = expr_copy(exps);
    ss_union(s1, s2, &u, &b);             /* b = exps with s2->s1 dummy remap */
    Expr* expsf = ss_do_subs(s2, exps);   /* if f(=s1) wins, expand s2 dummies */
    Expr* expsg = ss_do_subs(s1, exps);   /* if g(=s2) wins, expand s1 dummies */
    expr_free(exps);
    mrv_max3(s1, expsf, s2, expsg, &u, b, x, out, expr_out);
    ss_free(&u);
}

/* ---------------------------------------------------------------------- */
/* mrv(e,x) -- Algorithm 3.12 with SubsSet bookkeeping.                    */
/* Returns Omega in *out and e-rewritten-in-dummies in *expr_out (owned).   */
/* ---------------------------------------------------------------------- */
static void mrv(const Expr* e_in, const Expr* x, SubsSet* out, Expr** expr_out) {
    ss_init(out);
    if (g_fail || g_depth > GRUNTZ_MAX_DEPTH) { GZFAIL(); *expr_out = mk_int(0); return; }
    g_depth++;

    Expr* e = exp_combine(e_in);

    if (!has_x(e, x)) { *expr_out = e; g_depth--; return; }

    if (expr_eq(e, (Expr*)x)) {
        Expr* d = expr_copy(ss_get(out, x));
        *expr_out = d;
        expr_free(e); g_depth--; return;
    }

    if (is_mul(e) || is_add(e)) {
        const char* hd = is_mul(e) ? "Times" : "Plus";
        size_t n = e->data.function.arg_count;
        /* split into x-free (ifac) and x-dependent (dfac) operands */
        Expr** ifac = (Expr**)malloc(n * sizeof(Expr*)); size_t ni = 0;
        Expr** dfac = (Expr**)malloc(n * sizeof(Expr*)); size_t nd = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* a = e->data.function.args[i];
            if (has_x(a, x)) dfac[nd++] = expr_copy(a);
            else             ifac[ni++] = expr_copy(a);
        }
        Expr* ipart = combine_head(hd, ifac, ni); /* may be NULL; consumes ifac elems */
        free(ifac);
        if (nd == 0) { /* unreachable: e has x but no operand does */
            free(dfac); if (ipart) expr_free(ipart);
            GZFAIL(); *expr_out = mk_int(0); expr_free(e); g_depth--; return;
        }
        if (nd == 1) {
            SubsSet s; Expr* ex; mrv(dfac[0], x, &s, &ex);
            expr_free(dfac[0]); free(dfac);
            Expr* combined;
            if (ipart) { Expr* parts[2] = { ipart, ex };
                         combined = combine_head(hd, parts, 2); }
            else combined = ex;
            ss_copy(out, &s); ss_free(&s);
            *expr_out = combined;
        } else {
            /* a = dfac[0], b = combine(rest) */
            Expr* a = dfac[0];
            Expr* b = combine_head(hd, dfac + 1, nd - 1);
            free(dfac);
            SubsSet s1, s2; Expr *e1, *e2;
            mrv(a, x, &s1, &e1);
            mrv(b, x, &s2, &e2);
            expr_free(a); expr_free(b);
            /* combined = hd(ipart, e1, e2) */
            Expr* parts[3]; size_t np = 0;
            if (ipart) parts[np++] = ipart;
            parts[np++] = e1; parts[np++] = e2;
            Expr* combined = combine_head(hd, parts, np);
            mrv_max1(&s1, &s2, combined, x, out, expr_out);
            ss_free(&s1); ss_free(&s2);
        }
        expr_free(e); g_depth--; return;
    }

    if (is_pow(e) && !is_E(e->data.function.args[0])) {
        /* collapse a power tower b1^(e1) */
        Expr* b1 = expr_copy(e->data.function.args[0]);
        Expr* e1 = expr_copy(e->data.function.args[1]);
        while (is_pow(b1) && b1->data.function.arg_count == 2 &&
               !is_E(b1->data.function.args[0])) {
            Expr* nb = expr_copy(b1->data.function.args[0]);
            Expr* ne = simp(mk_times(e1, expr_copy(b1->data.function.args[1])));
            expr_free(b1); b1 = nb; e1 = ne;
        }
        if (has_x(e1, x)) {
            /* b1^e1 = exp(e1 * log(b1)) */
            Expr* rewritten = mk_exp(simp(mk_times(e1, mk_log(b1))));
            mrv(rewritten, x, out, expr_out);
            expr_free(rewritten);
        } else {
            SubsSet s; Expr* ex; mrv(b1, x, &s, &ex);
            expr_free(b1);
            *expr_out = simp(mk_pow(ex, e1));
            ss_copy(out, &s); ss_free(&s);
        }
        expr_free(e); g_depth--; return;
    }

    if (is_log(e)) {
        SubsSet s; Expr* ex; mrv(e->data.function.args[0], x, &s, &ex);
        *expr_out = simp(mk_log(ex));
        ss_copy(out, &s); ss_free(&s);
        expr_free(e); g_depth--; return;
    }

    {
        Expr* ea = exp_arg(e);
        if (ea) {
            if (is_log(ea)) { /* exp(log(u)) -> mrv(u) */
                mrv(ea->data.function.args[0], x, out, expr_out);
                expr_free(e); g_depth--; return;
            }
            Expr* li = limitinf(ea, x);
            bool inf = is_infinity_val(li);
            expr_free(li);
            if (g_fail) { *expr_out = mk_int(0); expr_free(e); g_depth--; return; }
            if (inf) {
                SubsSet s1; ss_init(&s1);
                Expr* e1 = expr_copy(ss_get(&s1, e));    /* dummy for exp(ea) */
                SubsSet s2; Expr* e2; mrv(ea, x, &s2, &e2);
                SubsSet su; ss_init(&su);
                ss_union(&s1, &s2, &su, NULL);
                Expr* expea = mk_exp(expr_copy(e2));
                ss_add_rewrite(&su, e1, expea);          /* rewrites[e1]=exp(e2) */
                /* mrv_max3(s1,e1, s2,exp(e2), su, e1) */
                mrv_max3(&s1, expr_copy(e1), &s2, expea, &su, expr_copy(e1),
                         x, out, expr_out);
                expr_free(e1); expr_free(e2);
                ss_free(&s1); ss_free(&s2); ss_free(&su);
            } else {
                SubsSet s; Expr* ex; mrv(ea, x, &s, &ex);
                *expr_out = mk_exp(ex);
                ss_copy(out, &s); ss_free(&s);
            }
            expr_free(e); g_depth--; return;
        }
    }

    /* Other functions of a single x-bearing argument: the mrv set is that of
     * the argument and the function is applied to the rewritten argument.
     * This is only sound for functions that are TRACTABLE at the limit of
     * their argument -- i.e. analytic there, so Mathilda's Series can expand
     * them (the circular/hyperbolic trig family and their inverses; Abs). It
     * is NOT sound for functions with an essential singularity at infinity
     * (Erf, ExpIntegralEi, Gamma, PolyGamma, Zeta, Bessel*, ...): those
     * require Gruntz's semi-tractable rewrites (thesis 5.2), not yet
     * implemented, and treating them generically yields WRONG answers. We
     * therefore whitelist the safe heads and refuse the rest, so an
     * unsupported special function becomes an honest "unevaluated" rather
     * than a silently wrong limit. */
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL) {
        static const char* SAFE[] = {
            "Sin","Cos","Tan","Cot","Sec","Csc",
            "ArcSin","ArcCos","ArcTan","ArcCot","ArcSec","ArcCsc",
            "Sinh","Cosh","Tanh","Coth","Sech","Csch",
            "ArcSinh","ArcCosh","ArcTanh","ArcCoth","Abs", NULL };
        const char* hn = e->data.function.head->data.symbol.name;
        bool safe = false;
        for (int i = 0; SAFE[i]; i++) if (strcmp(hn, SAFE[i]) == 0) { safe = true; break; }
        int xarg = -1; size_t n = e->data.function.arg_count; int cnt = 0;
        for (size_t i = 0; i < n; i++)
            if (has_x(e->data.function.args[i], x)) { xarg = (int)i; cnt++; }
        if (safe && cnt == 1) {
            SubsSet s; Expr* ex; mrv(e->data.function.args[(size_t)xarg], x, &s, &ex);
            Expr** args = (Expr**)malloc(n * sizeof(Expr*));
            for (size_t i = 0; i < n; i++)
                args[i] = (i == (size_t)xarg) ? ex : expr_copy(e->data.function.args[i]);
            Expr* r = expr_new_function(expr_copy(e->data.function.head), args, n);
            free(args);
            *expr_out = r;
            ss_copy(out, &s); ss_free(&s);
            expr_free(e); g_depth--; return;
        }
    }

    /* Unhandled shape. */
    GZFAIL();
    *expr_out = mk_int(0);
    expr_free(e); g_depth--;
}

/* ---------------------------------------------------------------------- */
/* compare(a,b,x): comparability-class order via limitinf(ln a / ln b).    */
/*   '<' if a < b, '=' if same class, '>' if a > b, '?' on failure.        */
/* ---------------------------------------------------------------------- */
static char compare(const Expr* a, const Expr* b, const Expr* x) {
    Expr* la = is_exp(a) ? expr_copy(exp_arg(a)) : mk_log(expr_copy((Expr*)a));
    Expr* lb = is_exp(b) ? expr_copy(exp_arg(b)) : mk_log(expr_copy((Expr*)b));
    Expr* ratio = simp(mk_times(la, mk_pow(lb, mk_int(-1))));
    Expr* c = limitinf(ratio, x);
    expr_free(ratio);
    char r;
    if (g_fail)                 r = '?';
    else if (is_lit_zero(c))    r = '<';
    else if (is_infinity_val(c))r = '>';
    else                        r = '=';
    expr_free(c);
    return r;
}

/* ---------------------------------------------------------------------- */
/* as_two_terms of a Mul/Add: a = first operand, b = combine(rest).        */
/* ---------------------------------------------------------------------- */

/* quick positivity test for x-bearing e (used by the log branch of sign). */
static bool is_positive_quick(const Expr* e, const Expr* x) {
    if (expr_eq((Expr*)e, (Expr*)x)) return true;
    if (is_exp(e)) return true;
    if (!has_x(e, x)) return numeric_sign_of(e) == 1;
    return false;
}

/* ---------------------------------------------------------------------- */
/* sign(e,x): sign of e as x->oo, per Gruntz. 0 with g_fail set = unknown. */
/* ---------------------------------------------------------------------- */
static int gz_sign_core(const Expr* e, const Expr* x) {
    if (g_fail || g_depth > GRUNTZ_MAX_DEPTH) { GZFAIL(); return 0; }
    g_depth++;
    int ret;

    if (!has_x(e, x)) {
        int s = numeric_sign_of(e);
        if (s == 2) { GZFAIL(); ret = 0; } else ret = s;
        g_depth--; return ret;
    }
    if (expr_eq((Expr*)e, (Expr*)x)) { g_depth--; return 1; }

    if (is_mul(e)) {
        int s = 1;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            int si = gz_sign(e->data.function.args[i], x);
            if (g_fail) { g_depth--; return 0; }
            if (si == 0) { g_depth--; return 0; }
            s *= si;
        }
        g_depth--; return s;
    }
    if (is_exp(e)) { g_depth--; return 1; }
    if (is_pow(e) && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        if (is_E(base)) { g_depth--; return 1; }
        int s = gz_sign(base, x);
        if (g_fail) { g_depth--; return 0; }
        if (s == 1) { g_depth--; return 1; }
        if (ex->type == EXPR_INTEGER) {
            int r = (ex->data.integer % 2 == 0) ? (s == 0 ? 0 : 1) : s;
            g_depth--; return r;
        }
    }
    if (is_log(e) && is_positive_quick(e->data.function.args[0], x)) {
        Expr* am1 = simp(mk_plus(expr_copy(e->data.function.args[0]), mk_int(-1)));
        int r = gz_sign(am1, x);
        expr_free(am1);
        g_depth--; return r;
    }

    /* fallback: sign of the leading coefficient of the series. */
    {
        Expr* c0; Expr* e0;
        mrv_leadterm(e, x, &c0, &e0);
        expr_free(e0);
        if (g_fail) { if (c0) expr_free(c0); g_depth--; return 0; }
        int r = gz_sign(c0, x);
        expr_free(c0);
        g_depth--; return r;
    }
}

/* ---------------------------------------------------------------------- */
/* Expression-tree height for choosing omega (least-nested exponential).   */
/* ---------------------------------------------------------------------- */
static int node_height(const SubsSet* Om, size_t idx, int* visiting) {
    if (visiting[idx]) return 1; /* cycle guard */
    visiting[idx] = 1;
    Expr* rexpr = ss_rewrite_of(Om, Om->var[idx]);
    int h = 1;
    if (rexpr) {
        for (size_t j = 0; j < Om->n; j++) {
            if (j == idx) continue;
            if (has_x(rexpr, Om->var[j])) h += node_height(Om, j, visiting);
        }
    }
    visiting[idx] = 0;
    return h;
}

/* ---------------------------------------------------------------------- */
/* rewrite: express the mrv elements as A*w^c, w->0+ (rule 3.2).           */
/* Returns the rewritten expr (owned) in *f_out and log(w) in *logw_out.   */
/* ---------------------------------------------------------------------- */
static void rewrite_omega(const Expr* exps, const SubsSet* Omega,
                          const Expr* x, const Expr* wsym,
                          Expr** f_out, Expr** logw_out) {
    size_t n = Omega->n;
    /* order indices by height desc; g = last (min height). */
    int* order = (int*)malloc(n * sizeof(int));
    int* heights = (int*)malloc(n * sizeof(int));
    int* visiting = (int*)calloc(n, sizeof(int));
    for (size_t i = 0; i < n; i++) { order[i] = (int)i; heights[i] = node_height(Omega, i, visiting); }
    /* simple insertion sort desc by height */
    for (size_t i = 1; i < n; i++) {
        int oi = order[i]; int hi = heights[oi]; size_t j = i;
        while (j > 0 && heights[order[j-1]] < hi) { order[j] = order[j-1]; j--; }
        order[j] = oi;
    }
    size_t gidx = (size_t)order[n - 1];
    Expr* g_exp = exp_arg(Omega->key[gidx]);   /* borrowed */
    int sig = gz_sign(g_exp, x);
    if (g_fail || (sig != 1 && sig != -1)) { GZFAIL();
        free(order); free(heights); free(visiting);
        *f_out = expr_copy((Expr*)exps); *logw_out = mk_int(0); return; }
    int wsign = (sig == 1) ? -1 : 1; /* wsym exponent sign factor */

    /* Build O2 substitutions and apply them in the sorted order. */
    Expr* f = expr_copy((Expr*)exps);
    for (size_t t = 0; t < n; t++) {
        size_t i = (size_t)order[t];
        Expr* f_exp = exp_arg(Omega->key[i]);      /* borrowed */
        /* c = limitinf(f_exp / g_exp) */
        Expr* ratio = simp(mk_times(expr_copy(f_exp), mk_pow(expr_copy(g_exp), mk_int(-1))));
        Expr* c = limitinf(ratio, x);
        expr_free(ratio);
        if (g_fail) { expr_free(c); expr_free(f);
            free(order); free(heights); free(visiting);
            *f_out = expr_copy((Expr*)exps); *logw_out = mk_int(0); return; }
        /* arg = rewrites[var] present ? its exp-arg : f_exp */
        Expr* rw = ss_rewrite_of(Omega, Omega->var[i]);
        Expr* arg = rw ? exp_arg(rw) : f_exp;      /* borrowed */
        /* repl = exp(arg - c*g_exp) * wsym^(wsign*c) */
        Expr* cg = simp(mk_times(expr_copy(c), expr_copy(g_exp)));
        Expr* aminus = simp(mk_plus(expr_copy(arg), mk_neg(cg)));
        Expr* wexp = (wsign == 1) ? c : simp(mk_neg(c)); /* consumes c if wsign 1 */
        Expr* repl = simp(mk_times(mk_exp(aminus), mk_pow(expr_copy((Expr*)wsym), wexp)));
        /* xreplace var -> repl */
        Expr* nf = xreplace(f, Omega->var[i], repl);
        expr_free(f); f = nf;
        expr_free(repl);
    }
    Expr* logw = expr_copy(g_exp);
    if (sig == 1) { Expr* nl = simp(mk_neg(logw)); logw = nl; }
    free(order); free(heights); free(visiting);
    *f_out = simp(f);
    *logw_out = logw;
}

/* ---------------------------------------------------------------------- */
/* Read the leading (coef, exponent) of a Series result in `w`.            */
/* Returns true and fills c0 (owned), num/den on success.                  */
/* ---------------------------------------------------------------------- */
static bool read_series_leading(Expr* sd, const Expr* w,
                                Expr** c0, int64_t* num, int64_t* den) {
    if (!sd) return false;
    /* constant in w */
    if (!has_x(sd, w)) { *c0 = expr_copy(sd); *num = 0; *den = 1; return true; }

    Expr* prefactor = NULL;
    Expr* s = sd;
    if (head_is(sd, SYM_Times)) {
        Expr* found = NULL; Expr* pf = mk_int(1);
        for (size_t i = 0; i < sd->data.function.arg_count; i++) {
            Expr* a = sd->data.function.args[i];
            if (head_is(a, SYM_SeriesData) && !found) found = a;
            else pf = simp(mk_times(pf, expr_copy(a)));
        }
        if (found) { s = found; prefactor = pf; } else expr_free(pf);
    }
    if (!head_is(s, SYM_SeriesData) || s->data.function.arg_count != 6) {
        if (prefactor) expr_free(prefactor);
        return false;
    }
    Expr* coefs  = s->data.function.args[2];
    Expr* nmin_e = s->data.function.args[3];
    Expr* den_e  = s->data.function.args[5];
    if (!head_is(coefs, SYM_List) || nmin_e->type != EXPR_INTEGER ||
        den_e->type != EXPR_INTEGER) { if (prefactor) expr_free(prefactor); return false; }
    int64_t nmin = nmin_e->data.integer;
    int64_t dd   = den_e->data.integer;
    size_t k = 0, kn = coefs->data.function.arg_count;
    while (k < kn && is_lit_zero(coefs->data.function.args[k])) k++;
    if (k == kn) { if (prefactor) expr_free(prefactor); return false; } /* escalate */
    Expr* lead = expr_copy(coefs->data.function.args[k]);
    if (prefactor) lead = simp(mk_times(prefactor, lead));
    *c0  = lead;
    *num = nmin + (int64_t)k;
    *den = dd;
    return true;
}

/* Expand every Log[G] with G depending on w (and G != w) into its truncated
 * power series in w (which surfaces the -k*Log[w] singular part), post-order
 * so inner logs are handled first, freezing already-produced Log[w] as the
 * constant `lw_sym`. This realises SymPy's logx semantics: after this pass
 * (plus a final Log[w] -> lw_sym replacement) the expression is a proper
 * Puiseux series in w whose coefficients may contain lw_sym, x, etc., which
 * Mathilda's Series can then expand even through an enclosing Exp[...]. */
static Expr* expand_logs(const Expr* e, const Expr* w, const Expr* x, const Expr* lw_sym,
                         int64_t order) {
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr* head = expand_logs(e->data.function.head, w, x, lw_sym, order);
    Expr** args = (n ? (Expr**)malloc(n * sizeof(Expr*)) : NULL);
    for (size_t i = 0; i < n; i++)
        args[i] = expand_logs(e->data.function.args[i], w, x, lw_sym, order);
    Expr* node = expr_new_function(head, args, n);
    if (args) free(args);

    if (is_log(node)) {
        Expr* G = node->data.function.args[0];
        if (has_x(G, w) && !expr_eq(G, (Expr*)w)) {
            /* freeze any inner Log[w] (from a just-expanded child) as lw_sym,
             * then expand this log via Series + Normal. */
            Expr* logw_node = mk_log(expr_copy((Expr*)w));
            Expr* G2 = xreplace(G, logw_node, lw_sym);
            expr_free(logw_node);
            /* Guard: expanding a log whose (frozen) argument mixes both w and
             * the limit variable x, or is large, makes the inner Series blow
             * up / loop (its coefficients are big x-expressions). Leave it
             * unexpanded; series_leadterm then bails to the leading-order
             * fallback. Pure-w log arguments (the well-behaved case) still
             * expand. */
            bool mixed = has_x(G2, w) && has_x(G2, (Expr*)x);
            if (has_x(G2, w) && !mixed && gz_leaf_count(G2) <= 120 &&
                g_work <= GRUNTZ_MAX_WORK) {
                /* Factor out the leading w-power first: Series would otherwise
                 * split a w-pole  Log[1/w - a]  into the branch-contaminated
                 * constants  Log[-a] + Log[-1/a]  (= 2 Pi I, not 0, because the
                 * frozen lw_sym has lost its sign). Writing
                 *   Log[G2] = k Log[w] + Log[H],   H = G2 w^{-k} -> nonzero const,
                 * makes Series[Log[H]] a clean pole-free Taylor series and keeps
                 * the singular part as the explicit  k Log[w]  the freeze step
                 * turns into  k lw_sym.  Only taken when H's w=0 value is a clean
                 * constant free of w and lw_sym; otherwise the direct expansion
                 * (and its guard) still applies. */
                Expr* factored = NULL;
                Expr* vspec = expr_new_function(mk_sym("List"),
                    (Expr*[]){ expr_copy((Expr*)w), mk_int(0), mk_int(order) }, 3);
                Expr* gsd = simp(mk_fn2("Series", expr_copy(G2), vspec));
                Expr* clead; int64_t vnum = 0, vden = 1;
                bool rok = read_series_leading(gsd, w, &clead, &vnum, &vden);
                /* H(0) = clead (the leading coefficient); usable only if it is a
                 * clean nonzero constant free of w and of the frozen lw_sym --
                 * otherwise Log[H] itself would carry a Log of the log-scale. */
                if (rok && vnum != 0 && !has_x(clead, w) && !has_x(clead, lw_sym) &&
                    zero_test_decide(clead) != ZERO_TEST_TRUE) {
                    /* Hh = Expand[G2 w^{-k}]: distributes so the w-pole cancels
                     * (w (1/w - a) -> 1 - a w) and Series[Log[Hh]] sees no pole. */
                    Expr* Hh = simp(mk_fn1("Expand",
                                    mk_times(expr_copy(G2),
                                        mk_pow(expr_copy((Expr*)w),
                                               make_exponent(-vnum, vden)))));
                    Expr* hspec = expr_new_function(mk_sym("List"),
                        (Expr*[]){ expr_copy((Expr*)w), mk_int(0), mk_int(order) }, 3);
                    Expr* logH = simp(mk_fn1("Normal",
                        mk_fn2("Series", mk_log(Hh), hspec)));
                    Expr* klogw = simp(mk_times(make_exponent(vnum, vden),
                                                mk_log(expr_copy((Expr*)w))));
                    factored = simp(mk_plus(klogw, logH));
                }
                if (rok) expr_free(clead);
                expr_free(gsd);   /* vspec consumed by mk_fn2 */
                /* A surviving Complex[] means the sign-tracking was still beaten
                 * at a deeper level; drop the contaminated expansion and leave
                 * the log unexpanded (series_leadterm's wsym/complex guards then
                 * abstain FAST rather than feeding Series a blow-up). */
                if (factored && contains_complex_head(factored)) {
                    expr_free(factored); factored = NULL;
                }
                if (factored) { expr_free(G2); expr_free(node); return factored; }

                Expr* spec = expr_new_function(mk_sym("List"),
                    (Expr*[]){ expr_copy((Expr*)w), mk_int(0), mk_int(order) }, 3);
                Expr* ser = simp(mk_fn1("Normal",
                    mk_fn2("Series", mk_log(G2), spec)));  /* mk_log consumes G2 */
                if (contains_complex_head(ser)) { expr_free(ser); return node; }
                expr_free(node);
                return ser;
            }
            expr_free(G2);
        }
    }
    return node;
}

/* series leading term of f in wsym; log(w) treated as the finite logw. */
/* leaf count (cheap complexity measure). */
static int64_t gz_leaf_count(const Expr* e) {
    if (!e) return 0;
    if (e->type != EXPR_FUNCTION) return 1;
    int64_t c = gz_leaf_count(e->data.function.head);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        c += gz_leaf_count(e->data.function.args[i]);
    return c;
}
/* true if some Log[g] in e has g depending on `sym` (the frozen log(w)
 * placeholder). Such a nested Log[log(w)] is a semi-tractable log that
 * Mathilda's Series mis-handles (it does not know log(w) < 0), so we route
 * those to the leading-order fallback instead. */
static bool contains_log_of_sym(const Expr* e, const Expr* sym) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (is_log(e) && has_x(e->data.function.args[0], sym)) return true;
    if (contains_log_of_sym(e->data.function.head, sym)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_log_of_sym(e->data.function.args[i], sym)) return true;
    return false;
}

/* True if e contains an imaginary unit (Complex[..]) anywhere. Used to reject a
 * frozen-scale Series input contaminated by a branch artefact: the mrv engine
 * works with real, positive log-scales, so once the P = -Log[w] substitution
 * has run, any surviving Complex[] means the sign-tracking was defeated at a
 * deeper log level (a 3+-level tower re-introduces Log[-P] -> Log[P] + I Pi
 * inside expand_logs) and feeding it to Series would blow up / return a wrong
 * complex value. Bail to the fallback / abstention instead. */
static bool contains_complex_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Complex)) return true;
    if (contains_complex_head(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_complex_head(e->data.function.args[i])) return true;
    return false;
}

static bool series_leadterm(const Expr* f, const Expr* wsym, const Expr* x, const Expr* logw,
                            Expr** c0_out, int64_t* num_out, int64_t* den_out) {
    Expr* lw_sym = fresh_dummy();
    Expr* logw_node = mk_log(expr_copy((Expr*)wsym));
    /* The frozen Log[w] is NEGATIVE (w -> 0+, so Log[w] -> -oo), hence the
     * log-scale  -Log[w] = P > 0.  Mathilda's Series does not know that, so a
     * Log[-lw_sym] (= Log[P], real) gets mis-branched to Log[lw_sym] + I Pi.
     * We therefore run Series in terms of the positive placeholder P = -lw_sym:
     * substitute lw_sym -> -P before Series (Log[-lw_sym] -> Log[P], clean) and
     * restore P -> -logw afterwards. */
    Expr* pos_sym = fresh_dummy();                       /* P = -Log[w] > 0 */
    Expr* neg_logw = simp(mk_neg(expr_copy((Expr*)logw)));/* P restores to -logw */

    bool ok = false;
    for (int64_t ord = GRUNTZ_SERIES_START; ord <= GRUNTZ_SERIES_MAX && !ok; ord *= 2) {
        /* expand inner logs and freeze Log[w] -> lw_sym (logx semantics) */
        Expr* fe = expand_logs(f, wsym, x, lw_sym, ord + 4);
        Expr* f3 = xreplace(fe, logw_node, lw_sym);
        expr_free(fe);

        /* Guard: an unexpanded Log[..w..] (expand_logs gave up on a big
         * argument) or an oversized expanded form makes Mathilda's Series
         * diverge/blow up; bail to the fallback. A Log[lw_sym] is NOT a bail
         * reason here -- the P-substitution below makes Series handle it. */
        if (contains_log_of_sym(f3, wsym) || gz_leaf_count(f3) > 300) {
            expr_free(f3); break;
        }

        /* Series in the positive log-scale P (= -lw_sym). */
        Expr* negP = mk_neg(expr_copy(pos_sym));
        Expr* f4 = simp(xreplace(f3, lw_sym, negP));
        expr_free(negP); expr_free(f3);

        /* A surviving Complex[] (I Pi) means the sign-tracking was defeated at a
         * deeper log level; Series would hang / go complex. Abstain instead. */
        if (contains_complex_head(f4)) { expr_free(f4); break; }

        Expr* spec = expr_new_function(mk_sym("List"),
            (Expr*[]){ expr_copy((Expr*)wsym), mk_int(0), mk_int(ord) }, 3);
        Expr* sd = simp(mk_fn2("Series", f4, spec));
        Expr* c0; int64_t num, den;
        if (read_series_leading(sd, wsym, &c0, &num, &den)) {
            Expr* c0r = xreplace(c0, pos_sym, neg_logw);   /* restore -Log[w] */
            expr_free(c0);
            if (has_x(c0r, wsym)) { expr_free(c0r); }  /* not a proper series */
            else { *c0_out = c0r; *num_out = num; *den_out = den; ok = true; }
        }
        expr_free(sd);
    }
    expr_free(logw_node); expr_free(lw_sym);
    expr_free(pos_sym); expr_free(neg_logw);
    return ok;
}

/* ---------------------------------------------------------------------- */
/* moveup: substitute x -> exp(x), collapsing Log[E^u] -> u.               */
/* ---------------------------------------------------------------------- */
static Expr* collapse_log_exp(const Expr* e) {
    if (is_log(e)) {
        Expr* a = exp_arg(e->data.function.args[0]);
        if (a) return collapse_log_exp(a);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr* head = collapse_log_exp(e->data.function.head);
    Expr** args = (n ? (Expr**)malloc(n * sizeof(Expr*)) : NULL);
    for (size_t i = 0; i < n; i++) args[i] = collapse_log_exp(e->data.function.args[i]);
    Expr* r = expr_new_function(head, args, n);
    if (args) free(args);
    return r;
}
static Expr* moveup_expr(const Expr* e, const Expr* x) {
    Expr* ex = mk_exp(expr_copy((Expr*)x));
    Expr* sub = xreplace(e, x, ex);
    expr_free(ex);
    Expr* col = collapse_log_exp(sub);
    expr_free(sub);
    Expr* r = simp(col);
    /* re-collapse after evaluation in case simp reintroduced Log[E^..] */
    Expr* r2 = collapse_log_exp(r);
    expr_free(r);
    return r2;
}

/* Build the exponent Expr num/den (den > 0). */
static Expr* make_exponent(int64_t num, int64_t den) {
    if (num == 0) return mk_int(0);
    if (den == 1) return mk_int(num);
    return simp(mk_times(mk_int(num), mk_pow(mk_int(den), mk_int(-1))));
}

/* sign of an exponent Expr: -1/0/+1, or 2 = "unknown". */
static int exponent_sign(const Expr* e0) {
    if (e0->type == EXPR_INTEGER)
        return e0->data.integer > 0 ? 1 : (e0->data.integer < 0 ? -1 : 0);
    if (e0->type == EXPR_REAL)
        return e0->data.real > 0.0 ? 1 : (e0->data.real < 0.0 ? -1 : 0);
    if (zero_test_decide(e0) == ZERO_TEST_TRUE) return 0;
    return numeric_sign_of(e0);
}
/* compare two exponents: -1 if a<b, 0 if a==b, +1 if a>b, 2 unknown. */
static int exponent_cmp(const Expr* a, const Expr* b) {
    Expr* d = simp(mk_plus(expr_copy((Expr*)a), mk_neg(expr_copy((Expr*)b))));
    int s = exponent_sign(d);
    expr_free(d);
    return s;
}

/* ---------------------------------------------------------------------- */
/* Leading-order calculus (fallback for when Mathilda's SeriesData cannot  */
/* represent the (irrational/symbolic) w-exponents). Computes (c0, e0) of   */
/* the leading term of e in w, with log(w) carried as the symbol LW.        */
/* Cannot resolve leading-order cancellation in a sum: returns false then    */
/* (the Series path handles those, always with rational exponents).         */
/* ---------------------------------------------------------------------- */
static bool leadterm_lo(const Expr* e, const Expr* w, const Expr* x, const Expr* LW,
                        Expr** c_out, Expr** e_out) {
    if (g_fail || g_depth > GRUNTZ_MAX_DEPTH) { GZFAIL(); return false; }
    if (!has_x(e, w)) { *c_out = expr_copy((Expr*)e); *e_out = mk_int(0); return true; }
    if (expr_eq((Expr*)e, (Expr*)w)) { *c_out = mk_int(1); *e_out = mk_int(1); return true; }

    if (is_pow(e) && e->data.function.arg_count == 2 &&
        !is_E(e->data.function.args[0])) {
        Expr* base = e->data.function.args[0];
        Expr* ex   = e->data.function.args[1];
        /* A constant exponent (free of both w and x) multiplies the base's
         * leading exponent directly. If the exponent depends on w or x, the
         * result would be an invalid (w-/x-dependent) leading exponent, so
         * rewrite base^ex = exp(ex*log(base)) and let the exp branch handle
         * it (log(w) becomes the constant symbol LW). */
        if (!has_x(ex, w) && !has_x(ex, x)) {
            Expr* cg; Expr* eg;
            if (!leadterm_lo(base, w, x, LW, &cg, &eg)) return false;
            *c_out = simp(mk_pow(cg, expr_copy(ex)));
            *e_out = simp(mk_times(eg, expr_copy(ex)));
            return true;
        } else {
            Expr* re = mk_exp(simp(mk_times(expr_copy(ex), mk_log(expr_copy(base)))));
            bool ok = leadterm_lo(re, w, x, LW, c_out, e_out);
            expr_free(re); return ok;
        }
    }
    if (is_exp(e)) {
        Expr* H = exp_arg(e);
        Expr* ch; Expr* eh;
        if (!leadterm_lo(H, w, x, LW, &ch, &eh)) return false;
        int s = exponent_sign(eh);
        expr_free(eh);
        if (s == 1) { expr_free(ch); *c_out = mk_int(1); *e_out = mk_int(0); return true; }
        if (s == 0) { *c_out = simp(mk_exp(ch)); *e_out = mk_int(0); return true; }
        expr_free(ch); return false; /* exp of a blow-up: fallback can't */
    }
    if (is_log(e)) {
        Expr* cg; Expr* eg;
        if (!leadterm_lo(e->data.function.args[0], w, x, LW, &cg, &eg)) return false;
        *c_out = simp(mk_plus(mk_log(cg), mk_times(eg, expr_copy((Expr*)LW))));
        *e_out = mk_int(0);
        return true;
    }
    if (is_mul(e)) {
        Expr* c = mk_int(1); Expr* ex = mk_int(0);
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            Expr* cf; Expr* ef;
            if (!leadterm_lo(e->data.function.args[i], w, x, LW, &cf, &ef)) { expr_free(c); expr_free(ex); return false; }
            c = simp(mk_times(c, cf)); ex = simp(mk_plus(ex, ef));
        }
        *c_out = c; *e_out = ex; return true;
    }
    if (is_add(e)) {
        size_t n = e->data.function.arg_count;
        Expr** cs = (Expr**)malloc(n * sizeof(Expr*));
        Expr** es = (Expr**)malloc(n * sizeof(Expr*));
        for (size_t i = 0; i < n; i++) {
            if (!leadterm_lo(e->data.function.args[i], w, x, LW, &cs[i], &es[i])) {
                for (size_t j = 0; j < i; j++) { expr_free(cs[j]); expr_free(es[j]); }
                free(cs); free(es); return false;
            }
        }
        /* find min exponent */
        size_t mi = 0;
        for (size_t i = 1; i < n; i++) {
            int c = exponent_cmp(es[i], es[mi]);
            if (c == 2) { for (size_t j = 0; j < n; j++) { expr_free(cs[j]); expr_free(es[j]); } free(cs); free(es); return false; }
            if (c < 0) mi = i;
        }
        /* sum coefficients sharing the min exponent */
        Expr* csum = mk_int(0); Expr* emin = expr_copy(es[mi]);
        for (size_t i = 0; i < n; i++) {
            if (exponent_cmp(es[i], emin) == 0) csum = simp(mk_plus(csum, expr_copy(cs[i])));
        }
        for (size_t i = 0; i < n; i++) { expr_free(cs[i]); expr_free(es[i]); }
        free(cs); free(es);
        if (zero_test_decide(csum) == ZERO_TEST_TRUE) { expr_free(csum); expr_free(emin); return false; /* cancellation */ }
        *c_out = csum; *e_out = emin; return true;
    }
    /* other function f[..g(w)..]: only the leading constant argument. */
    if (e->type == EXPR_FUNCTION) {
        int xa = -1; int cnt = 0; size_t n = e->data.function.arg_count;
        for (size_t i = 0; i < n; i++)
            if (has_x(e->data.function.args[i], w)) { xa = (int)i; cnt++; }
        if (cnt == 1) {
            Expr* cg; Expr* eg;
            if (!leadterm_lo(e->data.function.args[(size_t)xa], w, x, LW, &cg, &eg)) return false;
            int s = exponent_sign(eg); expr_free(eg);
            if (s == 0) {
                Expr** args = (Expr**)malloc(n * sizeof(Expr*));
                for (size_t i = 0; i < n; i++)
                    args[i] = (i == (size_t)xa) ? cg : expr_copy(e->data.function.args[i]);
                *c_out = simp(expr_new_function(expr_copy(e->data.function.head), args, n));
                free(args); *e_out = mk_int(0); return true;
            }
            expr_free(cg); return false;
        }
    }
    return false;
}

static bool leadterm_fallback(const Expr* f, const Expr* w, const Expr* x, const Expr* logw,
                              Expr** c0, Expr** e0) {
    Expr* LW = fresh_dummy();
    Expr* c; Expr* ex;
    bool ok = leadterm_lo(f, w, x, LW, &c, &ex);
    if (ok) {
        *c0 = xreplace(c, LW, logw); expr_free(c);
        *e0 = xreplace(ex, LW, logw); expr_free(ex);
    }
    expr_free(LW);
    return ok && !g_fail;
}

/* ---------------------------------------------------------------------- */
/* mrv_leadterm: (c0, e0) of e's series in its mrv variable (e0 an Expr).  */
/* ---------------------------------------------------------------------- */
static void mrv_leadterm_core(const Expr* e, const Expr* x,
                         Expr** c0_out, Expr** e0_out) {
    if (g_fail || g_depth > GRUNTZ_MAX_DEPTH || ++g_work > GRUNTZ_MAX_WORK) {
        GZFAIL(); *c0_out = mk_int(0); *e0_out = mk_int(0); return; }
    g_depth++;
    if (!has_x(e, x)) { *c0_out = expr_copy((Expr*)e); *e0_out = mk_int(0); g_depth--; return; }

    SubsSet Omega; Expr* exps;
    mrv(e, x, &Omega, &exps);
    if (g_fail) { ss_free(&Omega); expr_free(exps); *c0_out = mk_int(0); *e0_out = mk_int(0); g_depth--; return; }
    if (ss_empty(&Omega)) { *c0_out = exps; *e0_out = mk_int(0); ss_free(&Omega); g_depth--; return; }

    /* x in Omega? move up one level in the scale. */
    if (ss_find(&Omega, x) >= 0) {
        SubsSet up; ss_init(&up);
        for (size_t i = 0; i < Omega.n; i++) {
            Expr* nk = moveup_expr(Omega.key[i], x);
            (void)ss_get(&up, nk);
            expr_free(up.var[up.n - 1]); up.var[up.n - 1] = expr_copy(Omega.var[i]);
            expr_free(nk);
        }
        for (size_t i = 0; i < Omega.rw_n; i++) {
            Expr* ne = moveup_expr(Omega.rw_expr[i], x);
            ss_add_rewrite(&up, Omega.rw_var[i], ne);
            expr_free(ne);
        }
        Expr* eu = moveup_expr(exps, x);
        expr_free(exps); exps = eu;
        ss_free(&Omega); Omega = up;
    }

    Expr* wsym = fresh_dummy();
    Expr* f; Expr* logw;
    rewrite_omega(exps, &Omega, x, wsym, &f, &logw);
    expr_free(exps);
    if (g_fail) { expr_free(f); expr_free(logw); expr_free(wsym); ss_free(&Omega);
                  *c0_out = mk_int(0); *e0_out = mk_int(0); g_depth--; return; }

    Expr* c0; int64_t num, den;
    if (series_leadterm(f, wsym, x, logw, &c0, &num, &den)) {
        *c0_out = c0; *e0_out = make_exponent(num, den);
    } else {
        Expr* c0f; Expr* e0f;
        if (leadterm_fallback(f, wsym, x, logw, &c0f, &e0f)) {
            *c0_out = c0f; *e0_out = e0f;
        } else { GZFAIL(); *c0_out = mk_int(0); *e0_out = mk_int(0); }
    }
    expr_free(f); expr_free(logw); expr_free(wsym); ss_free(&Omega);
    g_depth--;
}

/* ---------------------------------------------------------------------- */
/* limitinf(e,x): lim_{x->oo} e.                                           */
/* ---------------------------------------------------------------------- */
static Expr* signed_infinity(int s) {
    if (s > 0) return mk_sym("Infinity");
    if (s < 0) return simp(mk_neg(mk_sym("Infinity")));
    return mk_sym("Indeterminate");
}

static Expr* limitinf_core(const Expr* e, const Expr* x) {
    if (g_fail) return mk_int(0);
    if (g_depth > GRUNTZ_MAX_DEPTH) { GZFAIL(); return mk_int(0); }
    if (!has_x(e, x)) return expr_copy((Expr*)e);
    g_depth++;

    Expr* c0; Expr* e0;
    mrv_leadterm(e, x, &c0, &e0);
    if (g_fail) { expr_free(c0); expr_free(e0); g_depth--; return mk_int(0); }

    int s = exponent_sign(e0);
    expr_free(e0);
    Expr* ret;
    if (s == 2) { expr_free(c0); GZFAIL(); ret = mk_int(0); }
    else if (s == 1) { expr_free(c0); ret = mk_int(0); }
    else if (s == -1) {
        int sc = gz_sign(c0, x);
        expr_free(c0);
        if (g_fail || sc == 0) { GZFAIL(); ret = mk_int(0); }
        else ret = signed_infinity(sc);
    } else {
        /* e0 == 0: limit is limit of the leading coefficient. */
        if (expr_eq(c0, (Expr*)e)) {
            Expr* cc = simp(mk_fn1("Cancel", expr_copy(c0)));
            if (expr_eq(cc, c0)) { expr_free(cc); expr_free(c0); GZFAIL(); ret = mk_int(0); }
            else { expr_free(c0); ret = limitinf(cc, x); expr_free(cc); }
        } else {
            ret = limitinf(c0, x);
            expr_free(c0);
        }
    }
    g_depth--;
    return ret;
}

/* ---------------------------------------------------------------------- */
/* Caching wrappers around the recursive cores.                            */
/* ---------------------------------------------------------------------- */
static int gz_sign(const Expr* e, const Expr* x) {
    if (has_x(e, x)) { CNode* c = cache_find(g_cache_sign, e); if (c) return c->ival; }
    int r = gz_sign_core(e, x);
    if (!g_fail && has_x(e, x)) { CNode* n = cache_insert(g_cache_sign, e); n->ival = r; }
    return r;
}
static void mrv_leadterm(const Expr* e, const Expr* x, Expr** c0_out, Expr** e0_out) {
    if (has_x(e, x)) {
        CNode* c = cache_find(g_cache_lt, e);
        if (c) { *c0_out = expr_copy(c->val); *e0_out = expr_copy(c->val2); return; }
    }
    mrv_leadterm_core(e, x, c0_out, e0_out);
    if (!g_fail && has_x(e, x)) {
        CNode* n = cache_insert(g_cache_lt, e);
        n->val = expr_copy(*c0_out); n->val2 = expr_copy(*e0_out);
    }
}
static Expr* limitinf(const Expr* e, const Expr* x) {
    if (has_x(e, x)) { CNode* c = cache_find(g_cache_limit, e); if (c) return expr_copy(c->val); }
    Expr* r = limitinf_core(e, x);
    if (!g_fail && r && has_x(e, x)) { CNode* n = cache_insert(g_cache_limit, e); n->val = expr_copy(r); }
    return r;
}

/* ---------------------------------------------------------------------- */
/* Driver: reduce lim_{x->point,dir} f to limitinf at +oo.                 */
/* ---------------------------------------------------------------------- */
static Expr* subst_eval(const Expr* e, const Expr* from, Expr* to /*owned*/) {
    Expr* rule = mk_fn2("Rule", expr_copy((Expr*)from), to);
    Expr* ra   = mk_fn2("ReplaceAll", expr_copy((Expr*)e), rule);
    return simp(ra);
}

/* Flatten Log[x^a] -> a Log[x], (x y)^c -> x^c y^c under x -> +oo (x > 0),
 * which shrinks log/exp towers so the mrv recursion stays shallow. Restricted
 * to {x}: coefficients and other variables are untouched. Consumes `e`. */
static Expr* gz_powerexpand(Expr* e, const Expr* x) {
    Expr* wlist = expr_new_function(mk_sym("List"),
                    (Expr*[]){ expr_copy((Expr*)x) }, 1);
    return simp(mk_fn2("PowerExpand", e, wlist));
}

/* Reject a "result" that is not a clean limit value: exactly +/-Infinity is
 * fine, but a ComplexInfinity / Indeterminate / DirectedInfinity, or an
 * Infinity buried inside a function (e.g. Sin[ComplexInfinity] from an
 * oscillating-at-a-pole limit), means the engine did not actually resolve the
 * limit -- return NULL so Limit is left for other layers / left unevaluated. */
static bool contains_bad_inf(const Expr* e) {
    if (!e) return false;
    if (is_sym_name(e, "Infinity") || is_sym_name(e, "ComplexInfinity") ||
        is_sym_name(e, "Indeterminate") || head_is(e, SYM_DirectedInfinity))
        return true;
    if (e->type == EXPR_FUNCTION) {
        if (contains_bad_inf(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (contains_bad_inf(e->data.function.args[i])) return true;
    }
    return false;
}
static bool gz_result_ok(const Expr* r) {
    if (!r) return false;
    if (is_sym_name(r, "Infinity")) return true;              /* +Infinity     */
    if (head_is(r, SYM_Times) && r->data.function.arg_count == 2) {
        Expr* a = r->data.function.args[0];
        Expr* b = r->data.function.args[1];
        if ((is_sym_name(a, "Infinity") && b->type == EXPR_INTEGER && b->data.integer < 0) ||
            (is_sym_name(b, "Infinity") && a->type == EXPR_INTEGER && a->data.integer < 0))
            return true;                                      /* -Infinity     */
    }
    return !contains_bad_inf(r);                              /* else: finite  */
}

/* ---------------------------------------------------------------------- */
/* Phase 2: essential-singularity isolation (thesis 5.2).                  */
/*                                                                         */
/* The mrv engine only understands exp-log expressions. A function with an */
/* essential singularity at +oo (Erf, Erfc, ExpIntegralEi) is invisible to */
/* it -- and Mathilda's Series cannot expand F[1/w] about w=0 either. But  */
/* Series CAN expand these functions at Infinity directly, surfacing the   */
/* singular part as an explicit exp factor (thesis transforms 5.6-5.7):    */
/*     Erf[y]  ->  1 - E^(-y^2) erf_s(y)                                    */
/*     Ei[y]   ->  E^y Ei_s(y)                                             */
/* where the residual erf_s/Ei_s are ordinary (Poincare) power series in   */
/* 1/y. We obtain that expansion from  Normal[Series[F[y],{y,Infinity,N}]] */
/* and substitute y -> g (the real argument), turning F[g] into a pure     */
/* exp-log expression the Phase-1 engine already handles. The number of    */
/* asymptotic terms N is escalated and two orders must agree, so a limit   */
/* whose value hides behind cancellation deeper than N is reported as an   */
/* honest gap rather than a truncation artefact.                           */
/* ---------------------------------------------------------------------- */
static bool is_semitractable_head(const Expr* e, const char** hn_out) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || e->data.function.arg_count != 1) return false;
    const char* hn = h->data.symbol.name;
    /* Erf/Erfc/Ei carry the singularity as an explicit Exp[-x^2]/Exp[x] prefix
     * times a Laurent tail; LogGamma's Stirling expansion is an *additive*
     * exp-log head plus a Laurent tail. Both keep the engine's series
     * expansions cheap. Gamma[g] = Exp[LogGamma[g]] wraps the additive Stirling
     * tower in an Exp, yielding an x^x-scale tower; the isolation branch below
     * builds it from the (cheap) LogGamma series and hands the mrv engine an
     * explicit Exp head. It was historically excluded because the deep
     * cancellations it appears in drove the exact multivariate zero test
     * (is_zero_poly) to exponential cost; with the FLINT-backed zero test that
     * bottleneck is gone, so Gamma is admitted again. */
    if (strcmp(hn, "Erf") == 0 || strcmp(hn, "Erfc") == 0 ||
        strcmp(hn, "ExpIntegralEi") == 0 || strcmp(hn, "LogGamma") == 0 ||
        strcmp(hn, "Gamma") == 0) {
        *hn_out = hn; return true;
    }
    return false;
}

static bool contains_semitractable(const Expr* e) {
    const char* hn;
    if (is_semitractable_head(e, &hn)) return true;
    if (e && e->type == EXPR_FUNCTION) {
        if (contains_semitractable(e->data.function.head)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (contains_semitractable(e->data.function.args[i])) return true;
    }
    return false;
}

static bool is_neg_infinity(const Expr* r) {
    if (!r || !head_is(r, SYM_Times) || r->data.function.arg_count != 2) return false;
    Expr* a = r->data.function.args[0];
    Expr* b = r->data.function.args[1];
    return (is_sym_name(a, "Infinity") && b->type == EXPR_INTEGER && b->data.integer < 0) ||
           (is_sym_name(b, "Infinity") && a->type == EXPR_INTEGER && a->data.integer < 0);
}

/* Normal[Series[head[y], {y, Infinity, nterms}]] with y -> arg. Returns a
 * fresh owned expr; on a Series failure it degrades to head[arg] (which the
 * engine then leaves unevaluated -- never wrong). */
static Expr* asymptotic_expansion(const char* head, const Expr* arg, int nterms) {
    Expr* d = fresh_dummy();
    Expr* fd = mk_fn1(head, expr_copy(d));
    Expr* spec = expr_new_function(mk_sym("List"),
        (Expr*[]){ expr_copy(d), mk_sym("Infinity"), mk_int(nterms) }, 3);
    Expr* norm = simp(mk_fn1("Normal", mk_fn2("Series", fd, spec)));
    Expr* sub  = xreplace(norm, d, arg);
    expr_free(norm);
    expr_free(d);
    return simp(sub);
}

/* Rewrite every semi-tractable node F[g] (g -> +/-oo) as its exp-log
 * asymptotic expansion; recurse bottom-up so nested arguments are isolated
 * first. Leaves the node untouched when g's limit is finite/undecidable. */
static Expr* isolate_semitractable(const Expr* e, const Expr* x, int nterms) {
    const char* hn;
    if (is_semitractable_head(e, &hn)) {
        Expr* g2 = isolate_semitractable(e->data.function.args[0], x, nterms);

        int saved = g_fail; g_fail = 0;
        Expr* L = limitinf(g2, x);
        bool lfail = g_fail;
        g_fail = saved;
        int cls = 0; /* +1: g->+oo, -1: g->-oo, 0: finite/unknown */
        if (!lfail && L) {
            if (is_sym_name(L, "Infinity")) cls = 1;
            else if (is_neg_infinity(L))    cls = -1;
        }
        if (L) expr_free(L);

        bool is_loggamma = (strcmp(hn, "LogGamma") == 0);

        if (cls == 1) {
            if (strcmp(hn, "Gamma") == 0) {
                /* Gamma[g] = Exp[LogGamma[g]]: build the additive Stirling tower
                 * from LogGamma's cheap series and wrap it in Exp, giving the
                 * mrv engine an explicit exp-log head to expand and cancel. */
                Expr* lg = asymptotic_expansion("LogGamma", g2, nterms);
                expr_free(g2);
                return simp(mk_exp(lg));
            }
            Expr* r = asymptotic_expansion(hn, g2, nterms);
            expr_free(g2);
            return simp(r);
        }
        if (cls == -1 && !is_loggamma) {
            if (strcmp(hn, "ExpIntegralEi") == 0) {
                /* Ei's asymptotic series E^y(1/y+1/y^2+...) is valid for
                 * y -> -oo as well (E^y -> 0). */
                Expr* r = asymptotic_expansion(hn, g2, nterms);
                expr_free(g2);
                return r;
            }
            /* Erf/Erfc: reflect to a +oo argument (thesis 5.16).
             *   Erf[g]  = -Erf[-g],   Erfc[g] = 2 - Erfc[-g].  */
            Expr* ng = simp(mk_neg(expr_copy(g2)));
            Expr* r  = asymptotic_expansion(hn, ng, nterms);
            expr_free(ng);
            expr_free(g2);
            return (strcmp(hn, "Erf") == 0)
                 ? simp(mk_neg(r))
                 : simp(mk_plus(mk_int(2), mk_neg(r)));
        }
        /* finite / undecidable / LogGamma at -oo (poles): keep the node
         * (honest gap). */
        return mk_fn1(hn, g2);
    }

    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    Expr* head = isolate_semitractable(e->data.function.head, x, nterms);
    size_t n = e->data.function.arg_count;
    Expr** args = n ? (Expr**)malloc(n * sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < n; i++)
        args[i] = isolate_semitractable(e->data.function.args[i], x, nterms);
    Expr* r = expr_new_function(head, args, n);
    if (args) free(args);
    return r;
}

/* Two computed limits agree if structurally equal or their difference is a
 * proven zero (guards against cosmetically-different equal constants). */
static bool same_limit_value(const Expr* a, const Expr* b) {
    if (expr_eq((Expr*)a, (Expr*)b)) return true;
    if (is_infinity_val(a) || is_infinity_val(b)) return false;
    Expr* diff = simp(mk_plus(expr_copy((Expr*)a), mk_neg(expr_copy((Expr*)b))));
    bool z = is_lit_zero(diff) || (zero_test_decide(diff) == ZERO_TEST_TRUE);
    expr_free(diff);
    return z;
}

/* Compute the limit at several asymptotic orders and accept a value only once
 * two orders concur (a lone order could be a truncation artefact). The needed
 * cancellation depth for these singularities is small, so we sweep a dense low
 * range: higher orders put ever-longer truncated series into denominators,
 * which the engine may fail to expand -- but as long as any two orders agree
 * on a clean value we trust it, and otherwise report an honest gap (NULL). */
static Expr* gruntz_semitractable_limit(const Expr* e0, const Expr* x) {
    static const int NS[] = { 3, 4, 5 };
    enum { NORD = (int)(sizeof(NS) / sizeof(NS[0])) };
    Expr* vals[NORD];
    int nv = 0;
    for (int i = 0; i < NORD; i++) {
        /* Each order gets a full work budget: the isolation-classification and
         * the main limit must not share (or exhaust) a single call's budget,
         * else a later order spuriously fails. */
        g_fail = 0; g_work = 0; cache_clear_all();
        Expr* iso = isolate_semitractable(e0, x, NS[i]);
        iso = gz_powerexpand(iso, x);
        /* Distribute so the constant parts of the asymptotic expansions (e.g.
         * the leading 1 in Erf -> 1 + ...) cancel across a difference before
         * the engine sees it; otherwise the engine must resolve a spurious
         * leading cancellation by deep series expansion and can give up. */
        iso = simp(mk_fn1("Expand", iso));
        g_fail = 0; g_work = 0; cache_clear_all();
        Expr* r = limitinf(iso, x);
        expr_free(iso);
        if (g_fail || !r || has_x(r, x) || !gz_result_ok(r)) {
            if (r) expr_free(r);
            continue;
        }
        for (int k = 0; k < nv; k++) {
            if (same_limit_value(vals[k], r)) {
                for (int j = 0; j < nv; j++) expr_free(vals[j]);
                return r;   /* two orders concur -> trust it */
            }
        }
        if (nv < NORD) vals[nv++] = r; else expr_free(r);
    }
    for (int j = 0; j < nv; j++) expr_free(vals[j]);
    return NULL;
}

Expr* gruntz_limit(Expr* f, Expr* x, Expr* point, int dir, int depth) {
    if (!f || !x || !point) return NULL;
    if (x->type != EXPR_SYMBOL) return NULL;

    g_fail = 0;
    g_depth = depth;
    g_work = 0;
    if (g_depth < 1) g_depth = 1;
    cache_clear_all();

    Expr* e0 = NULL;
    bool pos_inf = is_sym_name(point, "Infinity");
    bool neg_inf = (head_is(point, SYM_Times) && point->data.function.arg_count == 2 &&
                    (is_sym_name(point->data.function.args[0], "Infinity") ||
                     is_sym_name(point->data.function.args[1], "Infinity")));
    /* recompute neg_inf robustly: -Infinity is Times[-1, Infinity] */
    neg_inf = false;
    if (head_is(point, SYM_Times) && point->data.function.arg_count == 2) {
        Expr* a = point->data.function.args[0];
        Expr* b = point->data.function.args[1];
        Expr* inf = is_sym_name(a, "Infinity") ? a : (is_sym_name(b, "Infinity") ? b : NULL);
        Expr* oth = (inf == a) ? b : a;
        if (inf && oth->type == EXPR_INTEGER && oth->data.integer < 0) neg_inf = true;
    }

    if (pos_inf) {
        e0 = expr_copy(f);
    } else if (neg_inf) {
        e0 = subst_eval(f, x, mk_neg(expr_copy(x)));
    } else if (point->type == EXPR_SYMBOL || point->type == EXPR_INTEGER ||
               point->type == EXPR_REAL || head_is(point, SYM_Rational) ||
               head_is(point, SYM_Plus) || head_is(point, SYM_Times) ||
               point->type == EXPR_BIGINT) {
        /* finite point: x -> point +/- 1/x */
        Expr* recip = mk_pow(expr_copy(x), mk_int(-1));
        Expr* shifted = (dir == GZ_DIR_FROMBELOW)
            ? simp(mk_plus(expr_copy(point), mk_neg(recip)))
            : simp(mk_plus(expr_copy(point), recip));
        e0 = subst_eval(f, x, shifted);
        if (dir == GZ_DIR_TWOSIDED) {
            /* also compute from below and require agreement */
            Expr* recip2 = mk_pow(expr_copy(x), mk_int(-1));
            Expr* below = simp(mk_plus(expr_copy(point), mk_neg(recip2)));
            Expr* e0b = subst_eval(f, x, below);
            Expr* rb = limitinf(e0b, x);
            expr_free(e0b);
            Expr* ra = limitinf(e0, x);
            expr_free(e0);
            if (g_fail || has_x(ra, x) || has_x(rb, x) ||
                !gz_result_ok(ra) || !gz_result_ok(rb)) {
                expr_free(ra); expr_free(rb); return NULL; }
            Expr* result = expr_eq(ra, rb) ? ra : NULL;
            if (result) { expr_free(rb); return result; }
            expr_free(ra); expr_free(rb);
            return NULL;
        }
    } else {
        return NULL; /* unsupported point */
    }

    e0 = gz_powerexpand(e0, x);

    /* Phase 2: a special function with an essential singularity at infinity
     * needs the isolation pre-pass (thesis 5.2) before the mrv engine. */
    if (contains_semitractable(e0)) {
        Expr* r = gruntz_semitractable_limit(e0, x);
        expr_free(e0);
        return r; /* already validated, or NULL for an honest gap */
    }

    Expr* r = limitinf(e0, x);
    expr_free(e0);
    if (g_fail || !r) { if (r) expr_free(r); return NULL; }
    if (has_x(r, x) || !gz_result_ok(r)) { expr_free(r); return NULL; }
    return r;
}

void gruntz_init(void) { /* engine is reached via limit.c; nothing to register. */ }
