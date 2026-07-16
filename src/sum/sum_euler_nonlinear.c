/*
 * sum_euler_nonlinear.c -- Sum`EulerNonlinear: nonlinear Euler sums via MZV
 * reduction (weight <= 7).
 *
 * Evaluates
 *
 *     Sum_{k>=1} ( prod_i HarmonicNumber[k, p_i] ) / k^q          (H-case)
 *     Sum_{k>=1} k^{-q} prod_i PolyGamma[m_i, k]                  (tail-case)
 *
 * as polynomials in Riemann zeta values, whenever the total weight is <= 7 so
 * every multiple zeta value produced reduces to single zetas (the first
 * irreducible MZV is zeta(6,2) at weight 8).
 *
 * Pipeline (all exact / algorithmic):
 *   1. STUFFLE:  prod_i H_k^{(p_i)}  ->  Z-linear combo of multiple harmonic
 *      sums H_k(s_1,...,s_r), via the quasi-shuffle recursion.
 *   2. LIFT:  attach the outer k^{-q} -- H-case  T(q;s) = zeta(q,s)+zeta(q+s_1,..),
 *      tail-case  That(q;s) = zeta(s,q)+zeta(..,s_r+q)  (min-counting dual);
 *      producing a combo of multiple zeta values.
 *   3. REDUCE each MZV to single zetas: depth 1 direct; depth 2 by Euler's
 *      zeta(a,1), the reflection, and the odd-weight double-zeta formula (shared
 *      with sum_euler.c); depth >=3 by a curated weight-<=7 table.  Any word
 *      that does not fully reduce (weight-8 irreducibles) forces the whole sum
 *      to fall through (NULL), so a wrong or non-closed value is never emitted.
 *
 * PolyGamma factors use  PolyGamma[m,k] = (-1)^(m+1) m! (zeta(m+1) - H_{k-1}^(m+1))
 * = (-1)^(m+1) m! * tail_k^(m+1), and the tail (min-counting) lift keeps the
 * individually-divergent pieces together.
 *
 * Runs AFTER the narrow Sum`Euler stage, so the existing linear/quadratic fast
 * paths and their tests are untouched.
 *
 * Memory contract: takes ownership of res but must not free it; returns an
 * owned closed form or NULL to fall through.
 */

#include "sum_internal.h"
#include "sum_euler_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "sym_names.h"
#include "arithmetic.h"   /* is_infinity_sym, make_rational */
#include <stdlib.h>
#include <string.h>

#define NL_MAX_WEIGHT 7
#define NLM 10           /* max MZV depth (weight<=7 => depth<=7) */
#define NL_MAXFAC 16     /* max harmonic/polygamma factors */

/* ---------------- word polynomials (Z-linear combos of MZV words) --------- */

typedef struct { int s[NLM]; int d; } Word;
typedef struct { Word w; int64_t num, den; } WT;
typedef struct { WT* v; int n, cap; } WP;

static void wp_init(WP* p) { p->v = NULL; p->n = 0; p->cap = 0; }
static void wp_free(WP* p) { free(p->v); p->v = NULL; p->n = 0; p->cap = 0; }

static bool word_eq(const Word* a, const Word* b) {
    if (a->d != b->d) return false;
    for (int i = 0; i < a->d; i++) if (a->s[i] != b->s[i]) return false;
    return true;
}
static int64_t nl_gcd(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a ? a : 1;
}
static void wp_add(WP* p, Word w, int64_t num, int64_t den) {
    if (num == 0) return;
    if (den < 0) { den = -den; num = -num; }
    int64_t g = nl_gcd(num, den); num /= g; den /= g;
    for (int i = 0; i < p->n; i++) {
        if (word_eq(&p->v[i].w, &w)) {
            int64_t n2 = p->v[i].num * den + num * p->v[i].den;
            int64_t d2 = p->v[i].den * den;
            int64_t gg = nl_gcd(n2, d2); n2 /= gg; d2 /= gg;
            if (d2 < 0) { d2 = -d2; n2 = -n2; }
            p->v[i].num = n2; p->v[i].den = d2;
            if (n2 == 0) { p->v[i] = p->v[p->n - 1]; p->n--; }
            return;
        }
    }
    if (p->n == p->cap) {
        p->cap = p->cap ? p->cap * 2 : 8;
        p->v = realloc(p->v, (size_t)p->cap * sizeof(WT));
    }
    p->v[p->n].w = w; p->v[p->n].num = num; p->v[p->n].den = den; p->n++;
}

/* stuffle (quasi-shuffle) of two words a, b; accumulate coeff-1 words into out. */
static void nl_stuffle(const int* a, int na, const int* b, int nb, WP* out) {
    if (na == 0) { Word w; w.d = nb; for (int i = 0; i < nb; i++) w.s[i] = b[i]; wp_add(out, w, 1, 1); return; }
    if (nb == 0) { Word w; w.d = na; for (int i = 0; i < na; i++) w.s[i] = a[i]; wp_add(out, w, 1, 1); return; }
    /* recurse then prepend the leading letter */
    struct { int lead; const int* aa; int naa; const int* bb; int nbb; } cases[3] = {
        { a[0],        a + 1, na - 1, b,     nb     },
        { b[0],        a,     na,     b + 1, nb - 1 },
        { a[0] + b[0], a + 1, na - 1, b + 1, nb - 1 },
    };
    for (int c = 0; c < 3; c++) {
        WP t; wp_init(&t);
        nl_stuffle(cases[c].aa, cases[c].naa, cases[c].bb, cases[c].nbb, &t);
        for (int i = 0; i < t.n; i++) {
            Word w; w.d = t.v[i].w.d + 1; w.s[0] = cases[c].lead;
            for (int j = 0; j < t.v[i].w.d; j++) w.s[j + 1] = t.v[i].w.s[j];
            wp_add(out, w, t.v[i].num, t.v[i].den);
        }
        wp_free(&t);
    }
}

/* prod_i H_k^{(orders[i])} as a Z-combo of multiple harmonic sums. */
static WP nl_mhs_product(const int* orders, int d) {
    WP acc; wp_init(&acc);
    { Word w0; w0.d = 1; w0.s[0] = orders[0]; wp_add(&acc, w0, 1, 1); }
    for (int i = 1; i < d; i++) {
        WP na; wp_init(&na);
        int single[1] = { orders[i] };
        for (int j = 0; j < acc.n; j++) {
            WP tmp; wp_init(&tmp);
            nl_stuffle(acc.v[j].w.s, acc.v[j].w.d, single, 1, &tmp);
            for (int k = 0; k < tmp.n; k++)
                wp_add(&na, tmp.v[k].w, tmp.v[k].num * acc.v[j].num, tmp.v[k].den * acc.v[j].den);
            wp_free(&tmp);
        }
        wp_free(&acc); acc = na;
    }
    return acc;
}

/* H-case lift:  T(q; s_1..s_r) = zeta(q,s) + zeta(q+s_1, s_2..). */
static WP nl_lift_T(const WP* mhs, int q) {
    WP out; wp_init(&out);
    for (int i = 0; i < mhs->n; i++) {
        Word w = mhs->v[i].w;
        Word w1; w1.d = w.d + 1; w1.s[0] = q;
        for (int j = 0; j < w.d; j++) w1.s[j + 1] = w.s[j];
        wp_add(&out, w1, mhs->v[i].num, mhs->v[i].den);
        Word w2 = w; w2.s[0] = w.s[0] + q;                  /* merge q into first */
        wp_add(&out, w2, mhs->v[i].num, mhs->v[i].den);
    }
    return out;
}

/* tail-case (PolyGamma) lift with outer k^{-q} on the min index.
 *   q >= 1:  That(q; s) = zeta(s_1..s_r, q) + zeta(s_1..s_{r-1}, s_r+q)
 *   q == 0:  lower the innermost index by 1 (H_min^{(0)} = min). */
static WP nl_lift_tail(const WP* mhs, int q) {
    WP out; wp_init(&out);
    for (int i = 0; i < mhs->n; i++) {
        Word w = mhs->v[i].w;
        if (q == 0) {
            Word w1 = w; w1.s[w1.d - 1] -= 1;                /* lower innermost */
            wp_add(&out, w1, mhs->v[i].num, mhs->v[i].den);
        } else {
            Word w1 = w; w1.s[w1.d] = q; w1.d = w.d + 1;     /* append q */
            wp_add(&out, w1, mhs->v[i].num, mhs->v[i].den);
            Word w2 = w; w2.s[w2.d - 1] += q;                /* merge q into innermost */
            wp_add(&out, w2, mhs->v[i].num, mhs->v[i].den);
        }
    }
    return out;
}

/* ---------------- MZV -> single zetas ------------------------------------- */

static Expr* nl_zeta_pow(int n, int e) {
    return expr_new_function(expr_new_symbol(SYM_Power),
               (Expr*[]){ eu_zeta(n), sum_int(e) }, 2);
}

/* Depth-2 MZV zeta(s1,s2), s1>=2, -> zeta values, or NULL if irreducible. */
static Expr* nl_d2(int s1, int s2) {
    if (s1 < 2) return NULL;
    int w = s1 + s2;
    if (s2 == 1)
        return eu_plus2(euler_order1(s1), eu_times2(sum_int(-1), eu_zeta(s1 + 1)));
    if (s1 == s2)
        return eu_times2(eu_half(),
                   eu_plus2(nl_zeta_pow(s1, 2), eu_times2(sum_int(-1), eu_zeta(2 * s1))));
    if (w % 2 == 1) {
        if (s1 % 2 == 1) return euler_Z_odd(s1, s2);
        /* s2 odd: zeta(s1,s2) = zeta(s1)zeta(s2) - zeta(w) - zeta(s2,s1). */
        Expr* refl = eu_plus2(eu_times2(eu_zeta(s1), eu_zeta(s2)),
                              eu_times2(sum_int(-1), eu_zeta(w)));
        return eu_plus2(refl, eu_times2(sum_int(-1), euler_Z_odd(s2, s1)));
    }
    /* even weight */
    if (w == 6 && s1 == 4 && s2 == 2)   /* zeta(4,2) = zeta(3)^2 - 4/3 zeta(6) */
        return eu_plus2(nl_zeta_pow(3, 2), eu_times2(make_rational(-4, 3), eu_zeta(6)));
    if (w == 6 && s1 == 2 && s2 == 4) { /* reflection off zeta(4,2) */
        Expr* refl = eu_plus2(eu_times2(eu_zeta(2), eu_zeta(4)),
                              eu_times2(sum_int(-1), eu_zeta(6)));
        return eu_plus2(refl, eu_times2(sum_int(-1), nl_d2(4, 2)));
    }
    return NULL;   /* weight >= 8 even non-diagonal: irreducible */
}

/* c1/d1 * zeta(a) + c2/d2 * zeta(b)*zeta(cc) + c3/d3 * zeta(e)*zeta(g).
 * A zeta index of 0 means "that factor is absent" (so a lone zeta or a square
 * are built by repeating an index). Used for the weight-6/7 table entries. */
static Expr* nl_term(int64_t n, int64_t d, int z1, int z2) {
    Expr* zz = (z2 == 0) ? eu_zeta(z1)
                         : (z1 == z2 ? nl_zeta_pow(z1, 2)
                                     : eu_times2(eu_zeta(z1), eu_zeta(z2)));
    return eu_times2(make_rational(n, d), zz);
}

/* Depth->=3 curated table (weight <= 7).  Returns owned Expr or NULL. */
static Expr* nl_d3plus(const int* s, int d) {
    if (d == 3) {
        int a = s[0], b = s[1], c = s[2];
        if (a == 3 && b == 2 && c == 1)   /* zeta(3,2,1) = 3 z3^2 - 203/48 z6 */
            return eu_plus2(nl_term(3, 1, 3, 3), nl_term(-203, 48, 6, 0));
        if (a == 3 && b == 1 && c == 2)   /* zeta(3,1,2) = 53/24 z6 - 3/2 z3^2 */
            return eu_plus2(nl_term(53, 24, 6, 0), nl_term(-3, 2, 3, 3));
        if (a == 5 && b == 1 && c == 1)   /* zeta(5,1,1) = 5 z7 - 2 z2 z5 - 5/4 z3 z4 */
            return eu_plus2(nl_term(5, 1, 7, 0),
                       eu_plus2(nl_term(-2, 1, 2, 5), nl_term(-5, 4, 3, 4)));
        if (a == 4 && b == 1 && c == 2)   /* zeta(4,1,2) = 5/8 z7 + 5/2 z2 z5 - 15/4 z3 z4 */
            return eu_plus2(nl_term(5, 8, 7, 0),
                       eu_plus2(nl_term(5, 2, 2, 5), nl_term(-15, 4, 3, 4)));
        if (a == 4 && b == 2 && c == 1)   /* zeta(4,2,1) = -221/16 z7 + 11/2 z2 z5 + 7/2 z3 z4 */
            return eu_plus2(nl_term(-221, 16, 7, 0),
                       eu_plus2(nl_term(11, 2, 2, 5), nl_term(7, 2, 3, 4)));
        return NULL;
    }
    if (d == 4) {
        if (s[0] == 4 && s[1] == 1 && s[2] == 1 && s[3] == 1)
            /* zeta(4,1,1,1) = zeta(5,1,1) = 5 z7 - 2 z2 z5 - 5/4 z3 z4 */
            return eu_plus2(nl_term(5, 1, 7, 0),
                       eu_plus2(nl_term(-2, 1, 2, 5), nl_term(-5, 4, 3, 4)));
        return NULL;
    }
    return NULL;
}

/* Reduce a single MZV word to zeta values, or NULL if irreducible. */
static Expr* nl_mzv_reduce(const int* s, int d) {
    if (d <= 0) return NULL;
    if (s[0] < 2) return NULL;                 /* non-admissible / divergent */
    if (d == 1) return eu_zeta(s[0]);
    if (d == 2) return nl_d2(s[0], s[1]);
    return nl_d3plus(s, d);
}

/* Sum_i (num/den)_i * reduce(word_i), times prefactor; unevaluated.  NULL if
 * any word is irreducible. */
static Expr* nl_zeta_poly(const WP* mzvs, Expr* prefactor) {
    Expr* acc = NULL;
    for (int i = 0; i < mzvs->n; i++) {
        Expr* red = nl_mzv_reduce(mzvs->v[i].w.s, mzvs->v[i].w.d);
        if (!red) { if (acc) expr_free(acc); expr_free(prefactor); return NULL; }
        Expr* term = eu_times2(make_rational(mzvs->v[i].num, mzvs->v[i].den), red);
        acc = acc ? eu_plus2(acc, term) : term;
    }
    if (!acc) acc = sum_int(0);
    return eu_times2(prefactor, acc);
}

/* ---------------- recognizer --------------------------------------------- */

static bool nl_is_harm(Expr* g, Expr* var, int* p_out) {
    if (g->type != EXPR_FUNCTION) return false;
    Expr* h = g->data.function.head;
    if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_HarmonicNumber) return false;
    size_t ac = g->data.function.arg_count;
    if (ac < 1 || ac > 2) return false;
    if (!expr_eq(g->data.function.args[0], var)) return false;
    if (ac == 1) { *p_out = 1; return true; }
    Expr* pa = g->data.function.args[1];
    if (pa->type != EXPR_INTEGER || pa->data.integer < 1) return false;
    *p_out = (int)pa->data.integer;
    return true;
}
static bool nl_is_polygamma(Expr* g, Expr* var, int* m_out) {
    if (g->type != EXPR_FUNCTION) return false;
    Expr* h = g->data.function.head;
    if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_PolyGamma) return false;
    if (g->data.function.arg_count != 2) return false;
    Expr* ma = g->data.function.args[0];
    if (ma->type != EXPR_INTEGER || ma->data.integer < 1) return false;
    if (!expr_eq(g->data.function.args[1], var)) return false;
    *m_out = (int)ma->data.integer;
    return true;
}

Expr* builtin_sum_euler_nonlinear(Expr* res);

Expr* builtin_sum_euler_nonlinear(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite) return NULL;
    if (!is_infinity_sym(imax)) return NULL;
    if (imin->type != EXPR_INTEGER || imin->data.integer != 1) return NULL;

    Expr* fc = expr_copy(f);
    Expr* fn = evaluate(fc);
    expr_free(fc);   /* evaluate() borrows; free the copy we made */

    bool is_times = (fn->type == EXPR_FUNCTION
                     && fn->data.function.head->type == EXPR_SYMBOL
                     && fn->data.function.head->data.symbol.name == SYM_Times);
    size_t n = is_times ? fn->data.function.arg_count : 1;

    int Hord[NL_MAXFAC]; int nH = 0;
    int Tord[NL_MAXFAC]; int nT = 0;
    int64_t pf = 1;                     /* PolyGamma prefactor (integer) */
    Expr* rest = sum_int(1);
    bool bail = false;

    for (size_t i = 0; i < n && !bail; i++) {
        Expr* g = is_times ? fn->data.function.args[i] : fn;
        int pp, mm, e = 1;
        Expr* base = g;
        if (g->type == EXPR_FUNCTION && g->data.function.head->type == EXPR_SYMBOL
            && g->data.function.head->data.symbol.name == SYM_Power
            && g->data.function.arg_count == 2
            && g->data.function.args[1]->type == EXPR_INTEGER
            && g->data.function.args[1]->data.integer >= 1) {
            base = g->data.function.args[0];
            e = (int)g->data.function.args[1]->data.integer;
        }
        if (nl_is_harm(base, var, &pp)) {
            for (int j = 0; j < e && nH < NL_MAXFAC; j++) Hord[nH++] = pp;
        } else if (nl_is_polygamma(base, var, &mm)) {
            /* PolyGamma[m,k] = (-1)^(m+1) m! * tail^(m+1). */
            int64_t fac = 1; for (int t = 2; t <= mm; t++) fac *= t;   /* m! */
            int64_t sgn = (mm % 2 == 0) ? -1 : 1;                      /* (-1)^(m+1) */
            for (int j = 0; j < e; j++) { pf *= sgn * fac; if (nT < NL_MAXFAC) Tord[nT++] = mm + 1; }
        } else {
            Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ rest, expr_copy(g) }, 2);
            rest = evaluate(t); expr_free(t);
        }
    }
    expr_free(fn);

    /* Need at least two harmonic factors (or a squared/higher one), or a
     * PolyGamma product -- but not a mix (mixed H*PolyGamma deferred). */
    if ((nH > 0 && nT > 0) || (nH == 0 && nT == 0)) { expr_free(rest); return NULL; }
    if (nH == 1 && nT == 0) { expr_free(rest); return NULL; }  /* single H: Sum`Euler's job */

    /* rest must be c * var^(-q), c free of var, q >= 0. */
    Expr* one = sum_int(1);
    Expr* c = sum_subst(rest, var, one);
    expr_free(one);
    if (!sum_free_of(c, var)) { expr_free(rest); expr_free(c); return NULL; }
    Expr* cinv = expr_new_function(expr_new_symbol(SYM_Power),
                     (Expr*[]){ expr_copy(c), sum_int(-1) }, 2);
    Expr* mono = sum_eval("Times", (Expr*[]){ rest, cinv }, 2);   /* consumes rest,cinv */

    int q = 0; bool okmono = false;
    if (mono->type == EXPR_INTEGER && mono->data.integer == 1) { q = 0; okmono = true; }
    else if (mono->type == EXPR_FUNCTION && mono->data.function.head->type == EXPR_SYMBOL
             && mono->data.function.head->data.symbol.name == SYM_Power
             && mono->data.function.arg_count == 2
             && expr_eq(mono->data.function.args[0], var)
             && mono->data.function.args[1]->type == EXPR_INTEGER
             && mono->data.function.args[1]->data.integer < 0) {
        q = (int)(-mono->data.function.args[1]->data.integer); okmono = true;
    }
    expr_free(mono);
    if (!okmono) { expr_free(c); return NULL; }

    /* Convergence + weight gates. */
    int weight = q, conv;
    if (nH > 0) {
        for (int i = 0; i < nH; i++) weight += Hord[i];
        conv = q;                                   /* H-case converges for q >= 2 */
    } else {
        conv = q;
        for (int i = 0; i < nT; i++) { weight += Tord[i]; conv += Tord[i] - 1; }
    }
    if (conv < 2 || weight > NL_MAX_WEIGHT) { expr_free(c); return NULL; }

    /* Build MZV combination and reduce. */
    WP mhs = (nH > 0) ? nl_mhs_product(Hord, nH) : nl_mhs_product(Tord, nT);
    WP mzvs = (nH > 0) ? nl_lift_T(&mhs, q) : nl_lift_tail(&mhs, q);
    wp_free(&mhs);

    Expr* prefactor = eu_times2(c, make_rational(pf, 1));   /* consumes c */
    Expr* raw = nl_zeta_poly(&mzvs, prefactor);             /* consumes prefactor */
    wp_free(&mzvs);
    if (!raw) return NULL;

    /* Expand distributes the nested Times/Plus and collects like zeta monomials
     * (e.g. all the Zeta[7], Pi^2 Zeta[5], Pi^4 Zeta[3] terms) into one form. */
    Expr* out = sum_eval("Expand", (Expr*[]){ raw }, 1);   /* consumes raw */
    return out;
}

void sum_euler_nonlinear_init(void) {
    symtab_add_builtin("Sum`EulerNonlinear", builtin_sum_euler_nonlinear);
    symtab_get_def("Sum`EulerNonlinear")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Sum`EulerNonlinear",
        "Sum`EulerNonlinear[f, k, 1, Infinity] evaluates nonlinear Euler sums -- "
        "products of HarmonicNumber[k, p] (or PolyGamma[m, k]) over k^q -- as "
        "polynomials in Riemann zeta values, for total weight <= 7 (via multiple "
        "zeta value reduction). Returns unevaluated otherwise.");
}
