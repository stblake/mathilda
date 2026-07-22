/* ratcanon.c — unified rational normalization (see RATCANON_REWRITE_PLAN.md).
 *
 * PHASE 1 prototype only.  Proves the one-front-end + one-reduction pipeline on
 * the four representative regimes.  Throwaway — replaced by rat_canon_build /
 * rat_canon_reduce in Phases 2-3.
 */
#include "ratcanon.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "flint_bridge.h"
#include "risch_tower.h"    /* rt_expand_logs, rt_expand_exp_sums */
#include "rat_internal.h"   /* extract_num_den */
#include "expand.h"         /* expr_expand */
#include "core.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* ---- generator substitution map ---------------------------------------- */

typedef struct { Expr* kernel; char* name; } RcpGen;   /* kernel <-> $rcpN$ */
typedef struct {
    RcpGen* g; size_t n, cap;
    int ctr;
    const char* iname;   /* shared name for the imaginary unit, or NULL */
} RcpMap;

static void rcp_map_free(RcpMap* m) {
    for (size_t i = 0; i < m->n; i++) { expr_free(m->g[i].kernel); free(m->g[i].name); }
    free(m->g);
}

/* Return the fresh symbol name assigned to kernel `k` (dedup by expr_eq),
 * creating it on first sight.  `k` is copied. */
static const char* rcp_bind(RcpMap* m, const Expr* k) {
    for (size_t i = 0; i < m->n; i++)
        if (expr_eq(m->g[i].kernel, (Expr*)k)) return m->g[i].name;
    if (m->n == m->cap) {
        m->cap = m->cap ? m->cap * 2 : 8;
        m->g = realloc(m->g, m->cap * sizeof(RcpGen));
    }
    char buf[32];
    snprintf(buf, sizeof buf, "$rcp%d$", m->ctr++);
    char* nm = malloc(strlen(buf) + 1);
    strcpy(nm, buf);
    m->g[m->n].kernel = expr_copy((Expr*)k);
    m->g[m->n].name = nm;
    return m->g[m->n++].name;
}

/* Shared fresh symbol for the imaginary unit (kernel = I, so map-back + eval
 * applies I^2 -> -1). */
static const char* rcp_i_gen(RcpMap* m) {
    if (m->iname) return m->iname;
    Expr* i = expr_new_symbol("I");
    m->iname = rcp_bind(m, i);
    expr_free(i);
    return m->iname;
}

static int rcp_is_number(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return 1;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational) return 1;
    return 0;
}

static int rcp_is_indep_kernel_head(const char* h) {
    static const char* const ks[] = {
        "Log", "ArcSin","ArcCos","ArcTan","ArcCot","ArcSec","ArcCsc",
        "ArcSinh","ArcCosh","ArcTanh","ArcCoth","ArcSech","ArcCsch", NULL };
    for (int i = 0; ks[i]; i++) if (strcmp(h, ks[i]) == 0) return 1;
    return 0;
}

/* Substitute every kernel / algebraic constant / radical to a fresh free symbol.
 * Kernels are treated as opaque atoms (arguments not descended into — Phase 2
 * adds recursive argument normalization). */
static Expr* rcp_forward(const Expr* e, RcpMap* m) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const Expr* h = e->data.function.head;
    size_t ac = e->data.function.arg_count;
    Expr** av = e->data.function.args;
    if (h->type == EXPR_SYMBOL) {
        const char* hn = h->data.symbol.name;
        /* Complex[a,b] -> a + b*Igen */
        if (hn == SYM_Complex && ac == 2) {
            const char* ig = rcp_i_gen(m);
            return expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){
                expr_copy(av[0]),
                expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
                    expr_copy(av[1]), expr_new_symbol(ig) }, 2) }, 2);
        }
        /* Sqrt[symbolic] / symbolic^(p/q) -> fresh algebraic symbol */
        if (hn == SYM_Sqrt && ac == 1 && !rcp_is_number(av[0]))
            return expr_new_symbol(rcp_bind(m, e));
        if (hn == SYM_Power && ac == 2) {
            int64_t p, q;
            if (is_rational(av[1], &p, &q) && q != 1 && !rcp_is_number(av[0]))
                return expr_new_symbol(rcp_bind(m, e));
        }
        /* Log / inverse-trig -> fresh transcendental symbol */
        if (rcp_is_indep_kernel_head(hn) && ac == 1)
            return expr_new_symbol(rcp_bind(m, e));
    }
    /* recurse */
    Expr* nh = rcp_forward(h, m);
    Expr** na = malloc(sizeof(Expr*) * (ac ? ac : 1));
    for (size_t i = 0; i < ac; i++) na[i] = rcp_forward(av[i], m);
    Expr* r = expr_new_function(nh, na, ac);
    free(na);
    return r;
}

/* Substitute the fresh $rcpN$ symbols back to their kernels. */
static Expr* rcp_backward(const Expr* e, RcpMap* m) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < m->n; i++)
            if (strcmp(e->data.symbol.name, m->g[i].name) == 0)
                return expr_copy(m->g[i].kernel);
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t ac = e->data.function.arg_count;
    Expr* nh = rcp_backward(e->data.function.head, m);
    Expr** na = malloc(sizeof(Expr*) * (ac ? ac : 1));
    for (size_t i = 0; i < ac; i++) na[i] = rcp_backward(e->data.function.args[i], m);
    Expr* r = expr_new_function(nh, na, ac);
    free(na);
    return r;
}

/* RatCanonPrototype[expr]: substitute all kernels -> fresh free symbols, reduce
 * over Q with the existing plain-Q engine, map back, eval to apply the algebraic
 * relations.  Returns NULL (unevaluated) if the reduction declines. */
static Expr* builtin_ratcanon_prototype(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    RcpMap m; memset(&m, 0, sizeof m);
    Expr* sub = rcp_forward(arg, &m);

    Expr* reduced = flint_rational_together(sub);   /* one reduction over Q */
    expr_free(sub);
    if (!reduced) { rcp_map_free(&m); return NULL; }

    Expr* back = rcp_backward(reduced, &m);
    expr_free(reduced);
    rcp_map_free(&m);
    return eval_and_free(back);   /* applies I^2->-1, Sqrt[k]^2->k, ... */
}

/* ======================================================================== *
 *  Phase 2: the tower IR + builder (rat_canon_build)
 * ======================================================================== */

static int64_t rcp_igcd(int64_t a, int64_t b) {
    if (a < 0) a = -a; if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a ? a : 1;
}
static int64_t rcp_ilcm(int64_t a, int64_t b) {
    if (a < 0) a = -a; if (b < 0) b = -b;
    if (!a || !b) return 1;
    return a / rcp_igcd(a, b) * b;
}

/* rational gcd: gcd(an/ad, bn/bd) = gcd(an,bn)/lcm(ad,bd), reduced. */
static void rcp_qgcd(int64_t an, int64_t ad, int64_t bn, int64_t bd,
                     int64_t* gn, int64_t* gd) {
    int64_t n = rcp_igcd(an, bn), d = rcp_ilcm(ad, bd);
    int64_t g = rcp_igcd(n, d);
    *gn = n / g; *gd = d / g;
}

/* Split an exponent w into rational coeff (cn/cd) * base.  base is a fresh
 * owned Expr (Integer 1 when w is a pure number). */
static void rcp_exp_split(const Expr* w, int64_t* cn, int64_t* cd, Expr** base) {
    int64_t p, q;
    if (w->type == EXPR_INTEGER) { *cn = w->data.integer; *cd = 1; *base = expr_new_integer(1); return; }
    if (is_rational(w, &p, &q)) { *cn = p; *cd = q; *base = expr_new_integer(1); return; }
    if (w->type == EXPR_FUNCTION && w->data.function.head->type == EXPR_SYMBOL &&
        w->data.function.head->data.symbol.name == SYM_Times &&
        w->data.function.arg_count >= 2) {
        Expr* f0 = w->data.function.args[0];
        int64_t fp, fq;
        if (f0->type == EXPR_INTEGER || is_rational(f0, &fp, &fq)) {
            if (f0->type == EXPR_INTEGER) { *cn = f0->data.integer; *cd = 1; }
            else { *cn = fp; *cd = fq; }
            size_t rest = w->data.function.arg_count - 1;
            if (rest == 1) { *base = expr_copy(w->data.function.args[1]); return; }
            Expr** ra = malloc(sizeof(Expr*) * rest);
            for (size_t i = 0; i < rest; i++) ra[i] = expr_copy(w->data.function.args[i + 1]);
            *base = expr_new_function(expr_new_symbol(SYM_Times), ra, rest);
            free(ra);
            return;
        }
    }
    *cn = 1; *cd = 1; *base = expr_copy((Expr*)w);
}

static bool rcp_is_indep_head(const char* h) {
    if (h == SYM_Log || rcp_is_indep_kernel_head(h)) return true;
    /* Tangent-family: a single such kernel is a valid free monomial (Risch
     * RT_TAN).  Sin/Cos/Sec/Csc/Sinh/Cosh are algebraically dependent and are
     * left un-substituted (out of tower scope; declined downstream). */
    static const char* const tg[] = { "Tan","Cot","Tanh","Coth", NULL };
    for (int i = 0; tg[i]; i++) if (strcmp(h, tg[i]) == 0) return true;
    return false;
}

static bool rcp_is_const_name(const char* n) {
    static const char* const ks[] = { "Pi","E","I","EulerGamma","GoldenRatio",
        "Catalan","Degree","Glaisher","Khinchin","Infinity", NULL };
    for (int i = 0; ks[i]; i++) if (strcmp(n, ks[i]) == 0) return true;
    return false;
}

/* Exp-fundamental bookkeeping during a build. */
typedef struct { Expr* base; int64_t gn, gd; char* sym; } ExpFund;
/* Radical fundamental: all radicals of a common radicand r (Sqrt[r]=r^(1/2),
 * r^(p/q)) are powers of one generator r^(1/Q), Q = lcm of the roots.  The
 * radical analogue of ExpFund — r^(p/q) = (r^(1/Q))^(p*Q/q). */
typedef struct { Expr* radicand; int64_t Q; char* sym; } RadFund;
typedef struct {
    RatCanonForm* f;
    int ctr;
    ExpFund* ef; size_t nef, capef;
    RadFund* rf; size_t nrf, caprf;
} Builder;

static char* rcp_fresh(Builder* b) {
    char buf[32]; snprintf(buf, sizeof buf, "$rcg%d$", b->ctr++);
    char* s = malloc(strlen(buf) + 1); strcpy(s, buf); return s;
}

/* Pass 1: collect exp fundamentals.  For each E^w with an x-dependent exponent,
 * fold w's rational coeff into the gcd of its base's fundamental. */
static void rcp_collect_exps(const Expr* e, Builder* b) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power &&
        e->data.function.arg_count == 2 &&
        e->data.function.args[0]->type == EXPR_SYMBOL &&
        e->data.function.args[0]->data.symbol.name == SYM_E) {
        int64_t cn, cd; Expr* base;
        rcp_exp_split(e->data.function.args[1], &cn, &cd, &base);
        if (!(base->type == EXPR_INTEGER)) {   /* skip constant exponents */
            size_t i = 0;
            for (; i < b->nef; i++) if (expr_eq(b->ef[i].base, base)) break;
            if (i == b->nef) {
                if (b->nef == b->capef) { b->capef = b->capef ? b->capef * 2 : 4;
                    b->ef = realloc(b->ef, b->capef * sizeof(ExpFund)); }
                b->ef[i].base = expr_copy(base); b->ef[i].gn = cn < 0 ? -cn : cn;
                b->ef[i].gd = cd; b->ef[i].sym = NULL; b->nef++;
            } else {
                int64_t gn, gd; rcp_qgcd(b->ef[i].gn, b->ef[i].gd, cn, cd, &gn, &gd);
                b->ef[i].gn = gn; b->ef[i].gd = gd;
            }
        }
        expr_free(base);
        return;   /* do not descend into the exponent */
    }
    rcp_collect_exps(h, b);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rcp_collect_exps(e->data.function.args[i], b);
}

/* Note a radical of radicand `r` with root `q` (q>=2): fold q into the lcm Q of
 * r's fundamental group. */
static void rcp_rad_note(Builder* b, const Expr* r, int64_t q) {
    size_t i = 0;
    for (; i < b->nrf; i++) if (expr_eq(b->rf[i].radicand, (Expr*)r)) break;
    if (i == b->nrf) {
        if (b->nrf == b->caprf) { b->caprf = b->caprf ? b->caprf * 2 : 4;
            b->rf = realloc(b->rf, b->caprf * sizeof(RadFund)); }
        b->rf[i].radicand = expr_copy((Expr*)r); b->rf[i].Q = q; b->rf[i].sym = NULL; b->nrf++;
    } else b->rf[i].Q = rcp_ilcm(b->rf[i].Q, q);
}

/* Pass 1 (radicals): collect the fundamental root Q = lcm(roots) per radicand. */
static void rcp_collect_rads(const Expr* e, Builder* b) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        const char* hn = h->data.symbol.name;
        if (hn == SYM_Sqrt && e->data.function.arg_count == 1) {
            rcp_rad_note(b, e->data.function.args[0], 2);
            rcp_collect_rads(e->data.function.args[0], b);
            return;
        }
        if (hn == SYM_Power && e->data.function.arg_count == 2) {
            int64_t p, q;
            if (is_rational(e->data.function.args[1], &p, &q) && q != 1) {
                rcp_rad_note(b, e->data.function.args[0], q);
                rcp_collect_rads(e->data.function.args[0], b);
                return;
            }
        }
    }
    rcp_collect_rads(h, b);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rcp_collect_rads(e->data.function.args[i], b);
}

static Expr* rcp_radical_relation(const char* sym, int64_t deg, const Expr* radicand);
static const char* rcf_gen(RatCanonForm* f, Builder* b, Expr* kernel,
                           RcGenKind kind, Expr* relation);

/* Fundamental generator symbol for radicand `r`; creates the gen r^(1/Q) with
 * relation sym^Q - r on first use.  *Q_out receives the group's Q. */
static const char* rcp_rad_sym(Builder* b, RatCanonForm* f, const Expr* r, int64_t* Q_out) {
    size_t i = 0;
    for (; i < b->nrf; i++) if (expr_eq(b->rf[i].radicand, (Expr*)r)) break;
    if (i == b->nrf) { *Q_out = 0; return NULL; }
    RadFund* rd = &b->rf[i];
    *Q_out = rd->Q;
    if (!rd->sym) {
        Expr* onQ = expr_new_function(expr_new_symbol(SYM_Rational),
                        (Expr*[]){ expr_new_integer(1), expr_new_integer(rd->Q) }, 2);
        Expr* fk = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_copy((Expr*)r), onQ }, 2);   /* r^(1/Q) */
        const char* s = rcf_gen(f, b, fk, RCG_ALGEBRAIC,
                                rcp_radical_relation("$tmp$", rd->Q, r));
        for (size_t j = 0; j < f->n; j++)
            if (strcmp(f->gens[j].sym, s) == 0 && f->gens[j].relation) {
                expr_free(f->gens[j].relation);
                f->gens[j].relation = rcp_radical_relation(s, rd->Q, r);
                break;
            }
        rd->sym = (char*)s;
    }
    return rd->sym;
}

/* Add (or find) an algebraic/transcendental generator by its surface kernel;
 * returns its fresh symbol name (owned by the form). `relation` and `kernel` are
 * adopted on first insert, freed on a dedup hit. */
static const char* rcf_gen(RatCanonForm* f, Builder* b, Expr* kernel,
                           RcGenKind kind, Expr* relation) {
    for (size_t i = 0; i < f->n; i++)
        if (expr_eq(f->gens[i].kernel, kernel)) {
            expr_free(kernel); if (relation) expr_free(relation);
            return f->gens[i].sym;
        }
    if (f->n == f->cap) { f->cap = f->cap ? f->cap * 2 : 8;
        f->gens = realloc(f->gens, f->cap * sizeof(RcGen)); }
    f->gens[f->n].kernel = kernel;
    f->gens[f->n].sym = rcp_fresh(b);
    f->gens[f->n].kind = kind;
    f->gens[f->n].relation = relation;
    return f->gens[f->n++].sym;
}

/* Build relation `sym^deg - radicand` (deg>=2). */
static Expr* rcp_radical_relation(const char* sym, int64_t deg, const Expr* radicand) {
    Expr* pw = expr_new_function(expr_new_symbol(SYM_Power),
                 (Expr*[]){ expr_new_symbol(sym), expr_new_integer(deg) }, 2);
    Expr* neg = expr_new_function(expr_new_symbol(SYM_Times),
                 (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)radicand) }, 2);
    return expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ pw, neg }, 2);
}

/* Pass 2: substitute every kernel -> its generator symbol (E^(c*base) ->
 * sym^(c/g); Complex[a,b] -> a + b*symI; radical Power[r,p/q] -> sym^p). */
static Expr* rcp_build_forward(const Expr* e, Builder* b) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    RatCanonForm* f = b->f;
    const Expr* h = e->data.function.head;
    size_t ac = e->data.function.arg_count;
    Expr** av = e->data.function.args;
    if (h->type == EXPR_SYMBOL) {
        const char* hn = h->data.symbol.name;
        /* exponential */
        if (hn == SYM_Power && ac == 2 && av[0]->type == EXPR_SYMBOL &&
            av[0]->data.symbol.name == SYM_E) {
            int64_t cn, cd; Expr* base;
            rcp_exp_split(av[1], &cn, &cd, &base);
            if (base->type == EXPR_INTEGER) { expr_free(base); }
            else {
                size_t i = 0; for (; i < b->nef; i++) if (expr_eq(b->ef[i].base, base)) break;
                if (i < b->nef) {
                    ExpFund* ef = &b->ef[i];
                    if (!ef->sym) {
                        /* fundamental kernel E^((gn/gd)*base) */
                        Expr* fexp;
                        if (ef->gn == 1 && ef->gd == 1) fexp = expr_copy(ef->base);
                        else {
                            Expr* c = (ef->gd == 1) ? expr_new_integer(ef->gn)
                                : expr_new_function(expr_new_symbol(SYM_Rational),
                                    (Expr*[]){ expr_new_integer(ef->gn), expr_new_integer(ef->gd) }, 2);
                            fexp = expr_new_function(expr_new_symbol(SYM_Times),
                                    (Expr*[]){ c, expr_copy(ef->base) }, 2);
                        }
                        Expr* fk = expr_new_function(expr_new_symbol(SYM_Power),
                                    (Expr*[]){ expr_new_symbol(SYM_E), fexp }, 2);
                        ef->sym = (char*)rcf_gen(f, b, fk, RCG_TRANSCENDENTAL, NULL);
                    }
                    /* m = (cn/cd) / (gn/gd) = (cn*gd)/(cd*gn), integer */
                    int64_t m = (cn * ef->gd) / (cd * ef->gn);
                    expr_free(base);
                    if (m == 1) return expr_new_symbol(ef->sym);
                    return expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_new_symbol(ef->sym), expr_new_integer(m) }, 2);
                }
                expr_free(base);
            }
        }
        /* imaginary unit */
        if (hn == SYM_Complex && ac == 2) {
            Expr* rel = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){
                expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_new_symbol("I"), expr_new_integer(2) }, 2),
                expr_new_integer(1) }, 2);
            const char* s = rcf_gen(f, b, expr_new_symbol("I"), RCG_ALGEBRAIC, rel);
            return expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){
                expr_copy(av[0]),
                expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_copy(av[1]), expr_new_symbol(s) }, 2) }, 2);
        }
        /* Sqrt[r] = r^(1/2) = (r^(1/Q))^(Q/2) with Q = the group's lcm root. */
        if (hn == SYM_Sqrt && ac == 1) {
            int64_t Q; const char* s = rcp_rad_sym(b, f, av[0], &Q);
            if (s) {
                int64_t m = Q / 2;
                if (m == 1) return expr_new_symbol(s);
                return expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_new_symbol(s), expr_new_integer(m) }, 2);
            }
        }
        /* radical Power[r, p/q], q>1 -> sym^(p*Q/q) with sym = r^(1/Q). */
        if (hn == SYM_Power && ac == 2) {
            int64_t p, q;
            if (is_rational(av[1], &p, &q) && q != 1) {
                int64_t Q; const char* s = rcp_rad_sym(b, f, av[0], &Q);
                if (s) {
                    int64_t m = p * (Q / q);
                    if (m == 1) return expr_new_symbol(s);
                    return expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_new_symbol(s), expr_new_integer(m) }, 2);
                }
            }
        }
        /* Log / inverse-trig / Tan (transcendental, opaque) */
        if (rcp_is_indep_head(hn) && ac == 1)
            return expr_new_symbol(rcf_gen(f, b, expr_copy((Expr*)e), RCG_TRANSCENDENTAL, NULL));
    }
    /* recurse */
    Expr* nh = rcp_build_forward(h, b);
    Expr** na = malloc(sizeof(Expr*) * (ac ? ac : 1));
    for (size_t i = 0; i < ac; i++) na[i] = rcp_build_forward(av[i], b);
    Expr* r = expr_new_function(nh, na, ac);
    free(na);
    return r;
}

/* Pick a heuristic main variable: the first free symbol (not a generator, not a
 * numeric constant) encountered in num/den. */
static void rcp_find_var(const Expr* e, RatCanonForm* f, Expr** out) {
    if (*out || !e) return;
    if (e->type == EXPR_SYMBOL) {
        const char* n = e->data.symbol.name;
        if (n[0] == '$') return;                       /* generator symbol */
        if (rcp_is_const_name(n)) return;              /* Pi, E, I, ... */
        for (size_t i = 0; i < f->n; i++)
            if (f->gens[i].kernel->type == EXPR_SYMBOL &&
                strcmp(f->gens[i].kernel->data.symbol.name, n) == 0) return;
        *out = expr_new_symbol(n);
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    rcp_find_var(e->data.function.head, f, out);
    for (size_t i = 0; i < e->data.function.arg_count && !*out; i++)
        rcp_find_var(e->data.function.args[i], f, out);
}

/* stable sort: algebraic generators before transcendental (LEX-leading). */
static void rcp_order_gens(RatCanonForm* f) {
    for (size_t i = 1; i < f->n; i++) {
        RcGen key = f->gens[i]; size_t j = i;
        while (j > 0 && f->gens[j - 1].kind == RCG_TRANSCENDENTAL &&
               key.kind == RCG_ALGEBRAIC) { f->gens[j] = f->gens[j - 1]; j--; }
        f->gens[j] = key;
    }
}

RatCanonForm* rat_canon_build(const Expr* e) {
    RatCanonForm* f = calloc(1, sizeof(RatCanonForm));
    if (!f) return NULL;
    Builder b; memset(&b, 0, sizeof b); b.f = f;

    /* Pre-normalize: expose atomic kernels (Log[ab]->Log a+Log b, split E^(a+b)).
     * rt_expand_logs / rt_expand_exp_sums BORROW their argument (return a fresh
     * tree, do not free the input), so free each intermediate explicitly. */
    Expr* c0 = expr_copy((Expr*)e);
    Expr* mid = rt_expand_logs(c0);
    expr_free(c0);
    Expr* pre = rt_expand_exp_sums(mid);
    expr_free(mid);

    rcp_collect_exps(pre, &b);
    rcp_collect_rads(pre, &b);
    Expr* sub = rcp_build_forward(pre, &b);
    expr_free(pre);

    extract_num_den(sub, &f->num, &f->den);
    expr_free(sub);

    rcp_order_gens(f);
    rcp_find_var(f->num, f, &f->var);
    if (!f->var) rcp_find_var(f->den, f, &f->var);

    for (size_t i = 0; i < b.nef; i++) expr_free(b.ef[i].base);
    free(b.ef);
    for (size_t i = 0; i < b.nrf; i++) expr_free(b.rf[i].radicand);
    free(b.rf);
    return f;
}

void rat_canon_free(RatCanonForm* f) {
    if (!f) return;
    for (size_t i = 0; i < f->n; i++) {
        expr_free(f->gens[i].kernel); free(f->gens[i].sym);
        if (f->gens[i].relation) expr_free(f->gens[i].relation);
    }
    free(f->gens);
    if (f->num) expr_free(f->num);
    if (f->den) expr_free(f->den);
    if (f->var) expr_free(f->var);
    free(f);
}

/* Substitute generator symbols back to their kernels in `e`. */
static Expr* rcp_subst_back(const Expr* e, const RatCanonForm* f) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < f->n; i++)
            if (strcmp(e->data.symbol.name, f->gens[i].sym) == 0)
                return expr_copy(f->gens[i].kernel);
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t ac = e->data.function.arg_count;
    Expr* nh = rcp_subst_back(e->data.function.head, f);
    Expr** na = malloc(sizeof(Expr*) * (ac ? ac : 1));
    for (size_t i = 0; i < ac; i++) na[i] = rcp_subst_back(e->data.function.args[i], f);
    Expr* r = expr_new_function(nh, na, ac);
    free(na);
    return r;
}

Expr* rat_canon_subst_back(const RatCanonForm* f, const Expr* e) {
    return rcp_subst_back(e, f);
}

Expr* rat_canon_roundtrip(const RatCanonForm* f) {
    Expr* deninv = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_copy(f->den), expr_new_integer(-1) }, 2);
    Expr* frac = expr_new_function(expr_new_symbol(SYM_Times),
                     (Expr*[]){ expr_copy(f->num), deninv }, 2);
    Expr* back = rcp_subst_back(frac, f);
    expr_free(frac);
    return eval_and_free(back);
}

/* ======================================================================== *
 *  Phase 3: the reduction (rat_canon_reduce / rat_canon_normalize)
 * ======================================================================== */

/* §0.4 denominator-sign convention: flip only the FLINT all-negative-denominator
 * -(N)/(-D) artifact (never a mixed-sign polynomial).  Mirrors rat.c's
 * rc_sign_normalize; the one place the output convention is applied. */
/* Total degree of a monomial term: sum of integer powers of its variable factors
 * (radicals / numbers / kernels count 0).  Used to pick a Plus's leading terms. */
static int rco_term_degree(const Expr* t) {
    if (!t || t->type != EXPR_FUNCTION) return t && t->type == EXPR_SYMBOL ? 1 : 0;
    const char* h = t->data.function.head->type == EXPR_SYMBOL
                    ? t->data.function.head->data.symbol.name : NULL;
    if (h == SYM_Power && t->data.function.arg_count == 2 &&
        t->data.function.args[1]->type == EXPR_INTEGER) {
        return t->data.function.args[0]->type == EXPR_SYMBOL
               ? (int)t->data.function.args[1]->data.integer : 0;
    }
    if (h == SYM_Times) {
        int s = 0;
        for (size_t i = 0; i < t->data.function.arg_count; i++)
            s += rco_term_degree(t->data.function.args[i]);
        return s;
    }
    return 0;
}

static Expr* rc_expand(Expr* e);   /* fwd: Expand[e], defined below */

static Expr* rco_sign_normalize(Expr* r) {
    if (!r) return r;
    Expr* n; Expr* d;
    extract_num_den(r, &n, &d);
    bool flip;
    if (d->type == EXPR_FUNCTION && d->data.function.head->type == EXPR_SYMBOL &&
        d->data.function.head->data.symbol.name == SYM_Plus &&
        d->data.function.arg_count > 0) {
        /* Flip iff every MAXIMUM-total-degree term is negative — i.e. the leading
         * coefficient is negative (5 - x^2 -> x^2 - 5), never a mixed-degree
         * denominator whose leading term is positive (a^2 - a b, x - 1). */
        int maxdeg = 0;
        for (size_t i = 0; i < d->data.function.arg_count; i++) {
            int dg = rco_term_degree(d->data.function.args[i]);
            if (dg > maxdeg) maxdeg = dg;
        }
        flip = true;
        for (size_t i = 0; i < d->data.function.arg_count; i++)
            if (rco_term_degree(d->data.function.args[i]) == maxdeg &&
                !is_superficially_negative(d->data.function.args[i])) { flip = false; break; }
    } else flip = is_superficially_negative(d);
    if (!flip) {
        /* Even without a flip, Expand the numerator so an undistributed
         * Times[-1, Plus[...]] artifact of extract_num_den collapses. */
        Expr* nn = rc_expand(n);
        if (d->type == EXPR_INTEGER && d->data.integer == 1) { expr_free(d); expr_free(r); return nn; }
        expr_free(r);
        Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ d, expr_new_integer(-1) }, 2));
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ nn, inv }, 2));
    }
    Expr* nn = rc_expand(negate_expr(n));   /* distribute the sign */
    Expr* nd;
    if (d->type == EXPR_FUNCTION && d->data.function.head->type == EXPR_SYMBOL &&
        d->data.function.head->data.symbol.name == SYM_Plus) {
        size_t c = d->data.function.arg_count;
        Expr** ar = malloc(sizeof(Expr*) * c);
        for (size_t i = 0; i < c; i++) ar[i] = negate_expr(d->data.function.args[i]);
        nd = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), ar, c));
        free(ar);
    } else nd = negate_expr(d);
    expr_free(n); expr_free(d); expr_free(r);
    if (nd->type == EXPR_INTEGER && nd->data.integer == 1) { expr_free(nd); return nn; }
    Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ nd, expr_new_integer(-1) }, 2));
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ nn, inv }, 2));
}

/* True if `e` carries an algebraic atom (radical / Complex / root of unity) that
 * remains in a denominator — the reduction did not eliminate it, so either it is
 * WL-faithfully kept (1/(x-Sqrt2)) or the pre-formed cancellation needs a field
 * GCD this engine does not do (cube roots).  Either way, decline to classical. */
static bool rco_has_radical(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) {
        if (h->data.symbol.name == SYM_Sqrt || h->data.symbol.name == SYM_Complex) return true;
        if (h->data.symbol.name == SYM_Power && e->data.function.arg_count == 2) {
            int64_t p, q;
            if (is_rational(e->data.function.args[1], &p, &q) && q != 1) return true;
        }
    }
    if (rco_has_radical(h)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rco_has_radical(e->data.function.args[i])) return true;
    return false;
}

/* Phase 3c: complete a CONSTANT-radicand pre-formed cancellation via the
 * number-field GCD.  `num`/`den` are borrowed (den carries a radical the ideal
 * reduction could not eliminate).  Returns the cancelled, Expand-cleaned result
 * when the GCD is non-trivial AND the cancelled denominator is radical-free
 * ((x^2-2)/(x-Sqrt2) -> x+Sqrt2), else NULL (coprime / WL-kept / partial — the
 * classical path gives the canonical form). */
static Expr* rc_expand(Expr* e) {  /* owns e; returns Expand[e] */
    return eval_and_free(expr_new_function(expr_new_symbol("Expand"),
                                           (Expr*[]){ e }, 1));
}
static Expr* rat_canon_nf_complete(const Expr* num, const Expr* den) {
    Expr* g = flint_extension_gcd(num, den);
    if (!g) return NULL;
    Expr* ge = eval_and_free(g);
    if (ge->type == EXPR_INTEGER && ge->data.integer == 1) { expr_free(ge); return NULL; }
    Expr* nn = flint_extension_divexact(num, ge);
    Expr* nd = flint_extension_divexact(den, ge);
    expr_free(ge);
    if (!nn || !nd) { if (nn) expr_free(nn); if (nd) expr_free(nd); return NULL; }
    nn = rc_expand(nn);
    nd = rc_expand(nd);
    if (rco_has_radical(nd)) { expr_free(nn); expr_free(nd); return NULL; }
    /* Fold the divexact result nn/nd.  A numeric nd (often the sign unit -1)
     * makes the quotient a polynomial: Expand distributes the scalar so
     * (-Sqrt2-x)/(-1) collapses to Sqrt2+x rather than -(-Sqrt2-x). */
    int nd_numeric = (nd->type == EXPR_INTEGER || nd->type == EXPR_BIGINT ||
                      (nd->type == EXPR_FUNCTION && nd->data.function.head->type == EXPR_SYMBOL &&
                       nd->data.function.head->data.symbol.name == SYM_Rational));
    Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ nd, expr_new_integer(-1) }, 2));
    Expr* out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ nn, inv }, 2));
    return nd_numeric ? rc_expand(out) : out;
}

Expr* rat_canon_reduce(const RatCanonForm* f, RcMode mode) {
    (void)mode;   /* single-fraction reduction is mode-agnostic; the Cancel
                     per-term split happens in rat_canon_normalize */
    /* collect algebraic generators + relations, and the eliminable radicand
     * variables: a gen whose radicand (kernel argument) is a BARE SYMBOL
     * (Sqrt[b], y^(1/3)) — ordering that symbol above the generator lets the
     * ideal reduction eliminate it (enabling the cancellation).  A constant or
     * POLYNOMIAL radicand (Sqrt2, Sqrt[1+x^2]) is NOT eliminable — its variables
     * must stay below the generator so x^2 is never rewritten to g^2-1. */
    const char** asyms = malloc(sizeof(char*) * (f->n ? f->n : 1));
    const Expr** rels = malloc(sizeof(Expr*) * (f->n ? f->n : 1));
    const char** elim = malloc(sizeof(char*) * (f->n ? f->n : 1));
    int n_alg = 0, n_elim = 0;
    for (size_t i = 0; i < f->n; i++)
        if (f->gens[i].kind == RCG_ALGEBRAIC) {
            asyms[n_alg] = f->gens[i].sym; rels[n_alg] = f->gens[i].relation; n_alg++;
            Expr* k = f->gens[i].kernel;   /* Sqrt[r] / Power[r,1/q] / I */
            if (k->type == EXPR_FUNCTION && k->data.function.arg_count >= 1) {
                const char* hn = k->data.function.head->type == EXPR_SYMBOL
                                 ? k->data.function.head->data.symbol.name : NULL;
                Expr* r = k->data.function.args[0];
                if ((hn == SYM_Sqrt || hn == SYM_Power) && r->type == EXPR_SYMBOL &&
                    r->data.symbol.name[0] != '$')
                    elim[n_elim++] = r->data.symbol.name;
            }
        }

    /* frac = num / den in generator symbols */
    Expr* frac;
    if (f->den->type == EXPR_INTEGER && f->den->data.integer == 1) frac = expr_copy(f->num);
    else frac = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){
        expr_copy(f->num), expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ expr_copy(f->den), expr_new_integer(-1) }, 2) }, 2);

    Expr* reduced = flint_tower_reduce(frac, asyms, rels, n_alg, elim, n_elim);
    expr_free(frac); free(asyms); free(rels); free(elim);
    if (!reduced) return NULL;

    Expr* res = eval_and_free(rat_canon_subst_back(f, reduced));
    expr_free(reduced);

    /* Accept when the DENOMINATOR is radical-free — the result is fully reduced.
     * Covers: free/transcendental; sum-of-conjugates (Q(i), Sqrt[k]/Sqrt[d] sums
     * -> radical-free denominator); and, via the radicand-variable ordering in
     * flint_tower_reduce, the pre-formed VARIABLE-radicand cancellations whose
     * radical lands only in the numerator ((a^2-b)/(a-Sqrt b) -> a+Sqrt b;
     * (y-1)/(y^(1/3)-1) -> 1+y^(1/3)+y^(2/3)).  A radical REMAINING in the
     * denominator is either WL-faithfully kept (1/(x-Sqrt2)), a coprime
     * multi-radical sum, or a CONSTANT-radicand pre-formed cancellation needing a
     * number-field GCD (Phase 3c): decline to the classical path, which produces
     * the canonical form for all of those. */
    if (n_alg > 0) {
        Expr* num; Expr* den; extract_num_den(res, &num, &den);
        if (rco_has_radical(den)) {
            /* Phase 3c: constant-radicand pre-formed cancellation via number-field
             * GCD; NULL if coprime / WL-kept (decline to classical). */
            /* Phase 3c: constant-radicand pre-formed cancellation via number-field
             * GCD; NULL if coprime / WL-kept / unhandled -> decline to classical,
             * which gives the canonical form (and fully reduces the cases the
             * builder under-represents, e.g. commensurate radicals
             * y^(1/2)/y^(1/3)/y^(1/6)).  Blanket-keeping those regressed rat/
             * simplify — see plan Phase 3d. */
            Expr* done = rat_canon_nf_complete(num, den);
            expr_free(num); expr_free(den);
            if (!done) { expr_free(res); return NULL; }
            expr_free(res); res = done;
        } else { expr_free(num); expr_free(den); }
    }
    return rco_sign_normalize(res);
}

Expr* rat_canon_normalize(const Expr* e, RcMode mode) {
    /* Cancel maps the reduction over the top-level additive terms. */
    if (mode == RCM_CANCEL && e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Plus) {
        size_t nt = e->data.function.arg_count;
        Expr** parts = malloc(sizeof(Expr*) * nt);
        for (size_t i = 0; i < nt; i++) {
            parts[i] = rat_canon_normalize(e->data.function.args[i], mode);
            if (!parts[i]) { for (size_t j = 0; j < i; j++) expr_free(parts[j]); free(parts); return NULL; }
        }
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), parts, nt);
        free(parts);
        return eval_and_free(sum);
    }
    RatCanonForm* f = rat_canon_build(e);
    if (!f) return NULL;
    Expr* r = rat_canon_reduce(f, mode);
    rat_canon_free(f);
    return r;
}

/* Test builtins for the reduction (Phase 3). */
static Expr* builtin_ratcanon_normalize(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return rat_canon_normalize(res->data.function.args[0], RCM_TOGETHER);
}
static Expr* builtin_ratcanon_cancel(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return rat_canon_normalize(res->data.function.args[0], RCM_CANCEL);
}

void ratcanon_init(void) {
    symtab_add_builtin("RatCanonNormalize", builtin_ratcanon_normalize);
    symtab_get_def("RatCanonNormalize")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("RatCanonCancel", builtin_ratcanon_cancel);
    symtab_get_def("RatCanonCancel")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("RatCanonPrototype", builtin_ratcanon_prototype);
    symtab_get_def("RatCanonPrototype")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("RatCanonPrototype",
        "RatCanonPrototype[expr] (Phase-1 prototype) reduces a rational function "
        "over the differential/algebraic tower of expr via one FLINT reduction.");
}
