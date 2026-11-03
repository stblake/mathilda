/* zero_test_distrib.c -- PossibleZeroQ for generalized functions.
 *
 * A residual that carries DiracDelta / Derivative[k][DiracDelta] /
 * HeavisideTheta is a DISTRIBUTION, and pointwise sampling cannot decide it:
 * DiracDelta vanishes at every sample point, so a wrong delta coefficient is
 * invisible, while HeavisideTheta does not evaluate at a complex sample, so the
 * sampler cannot see the regular part either.  Both happen in the DSolve corpus:
 * the impulse response y'' + 2 y' + 2 y == DiracDelta[t - Pi] back-substitutes to
 *     -DiracDelta[t-Pi] - 2 E^(Pi-t) Cos[t] DiracDelta[t-Pi]
 *       - E^(Pi-t) Sin[t] DiracDelta'[t-Pi]
 * which is zero (sift: f(t) d(t-a) = f(a) d(t-a), and
 * g(t) d'(t-a) = g(a) d'(t-a) - g'(a) d(t-a)), and PossibleZeroQ said False.
 *
 * The decision is EXACT in structure and delegates every scalar question back to
 * the ordinary zero test:
 *   1. Every distribution atom must have an argument linear in ONE variable t,
 *      c (t - a) with c a non-zero real number and a a real number (the
 *      breakpoints are ordered numerically).  Otherwise decline (UNKNOWN).
 *   2. Atoms are normalised to the point a and replaced by placeholder symbols:
 *      d^(k)(c u) = d^(k)(u) / (c^k |c|),  Theta(c u) = Theta(u) (c > 0) or
 *      1 - Theta(u) (c < 0; equal almost everywhere, which is all a distribution
 *      sees).  The result is expanded in the placeholders.
 *   3. A term with one delta placeholder is sifted onto the basis d^(m)(t - a):
 *      c(t) d^(k)(t-a) = sum_j (-1)^j Binomial[k, j] c^(j)(a) d^(k-j)(t-a), after
 *      replacing any Theta(t - b), b != a, in c by its constant value near a.
 *      Each basis coefficient must be zero.  A product of two deltas, or a delta
 *      times a Theta at the same point, is not a distribution we can name:
 *      decline.
 *   4. The remaining regular part is piecewise analytic: on each interval between
 *      consecutive breakpoints every Theta is 0 or 1, and the piece must be zero.
 *      An analytic piece that vanishes on an open interval vanishes identically,
 *      so the piece is decided by the ordinary (complex-domain) zero test.
 * Verdict: FALSE if any coefficient or piece is non-zero (a non-zero regular
 * piece on an interval, or a non-zero point mass, is a non-zero distribution);
 * TRUE if all are zero; otherwise UNKNOWN. */
#include "zero_test_distrib.h"

#include "core.h"           /* leaf_count_internal */
#include "eval.h"
#include "expand.h"
#include "expr.h"
#include "numeric.h"
#include "sym_names.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ZD_MAX_ATOMS   16       /* distinct (kind, point, order) placeholders  */
#define ZD_MAX_LEAVES  20000    /* beyond this the expansion is not worth it   */
#define ZD_MAX_ORDER   8        /* highest delta derivative sifted             */

typedef enum { ZD_DELTA, ZD_THETA } ZdKind;

typedef struct {
    ZdKind kind;
    int    order;          /* derivative order of a delta; 0 for Theta          */
    Expr*  point;          /* exact breakpoint a (owned)                        */
    double pointd;         /* its numeric value, for ordering / identity        */
    Expr*  sym;            /* placeholder symbol (owned)                        */
} ZdAtom;

typedef struct {
    const char* var;       /* the single variable t (interned name), or NULL    */
    ZdAtom atoms[ZD_MAX_ATOMS];
    int    n;
    bool   fail;           /* an atom we cannot normalise                       */
} ZdCtx;

static bool is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}
static bool head_is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION && is_sym(e->data.function.head, name);
}

/* Derivative[k][DiracDelta][arg]: returns k (>= 1) and sets *arg, else 0.
The head of Derivative[k][DiracDelta][x] is Derivative[k][DiracDelta], whose
 * head is Derivative[k] -- i.e. e->head->head->head == Derivative. */
static int delta_deriv_order(const Expr* e, const Expr** arg) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 1) return 0;
    const Expr* h = e->data.function.head;               /* Derivative[k][DiracDelta] */
    if (!h || h->type != EXPR_FUNCTION || h->data.function.arg_count != 1) return 0;
    if (!is_sym(h->data.function.args[0], SYM_DiracDelta)) return 0;
    const Expr* dk = h->data.function.head;              /* Derivative[k]            */
    if (!head_is_sym(dk, SYM_Derivative) || dk->data.function.arg_count != 1) return 0;
    const Expr* k = dk->data.function.args[0];
    if (k->type != EXPR_INTEGER || k->data.integer < 1) return 0;
    *arg = e->data.function.args[0];
    return (int)k->data.integer;
}

/* Is `e` a distribution atom?  Sets kind/order/arg. */
static bool atom_of(const Expr* e, ZdKind* kind, int* order, const Expr** arg) {
    if (head_is_sym(e, SYM_DiracDelta) && e->data.function.arg_count == 1) {
        *kind = ZD_DELTA; *order = 0; *arg = e->data.function.args[0]; return true;
    }
    if (head_is_sym(e, SYM_HeavisideTheta) && e->data.function.arg_count == 1) {
        *kind = ZD_THETA; *order = 0; *arg = e->data.function.args[0]; return true;
    }
    int k = delta_deriv_order(e, arg);
    if (k > 0) { *kind = ZD_DELTA; *order = k; return true; }
    return false;
}

bool zt_has_distribution(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    ZdKind kd; int ord; const Expr* a;
    if (atom_of(e, &kd, &ord, &a)) return true;
    if (head_is_sym(e, SYM_DiracDelta) || head_is_sym(e, SYM_HeavisideTheta)) return true;
    if (zt_has_distribution(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (zt_has_distribution(e->data.function.args[i])) return true;
    return false;
}

/* The single non-constant symbol of `e` (Pi, E, I, ... are constants), or
 * NULL when there is none; sets *many when there is more than one. */
static void collect_var(const Expr* e, const char** var, bool* many) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (expr_is_numeric_quantity(e)) return;
        const char* n = e->data.symbol.name;
        if (!*var) *var = n;
        else if (*var != n) *many = true;
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect_var(e->data.function.args[i], var, many);
}

static Expr* call1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}
static Expr* call2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* rule(Expr* a, Expr* b) { return call2(SYM_Rule, a, b); }

/* e /. var -> val, evaluated.  e borrowed. */
static Expr* subst_var(const Expr* e, const char* var, const Expr* val) {
    return eval_and_free(call2(SYM_ReplaceAll, expr_copy((Expr*)e),
                               rule(expr_new_symbol(var), expr_copy((Expr*)val))));
}

/* Machine value of a real numeric quantity (Im == 0); false otherwise. */
static bool real_value(const Expr* e, double* out) {
    if (!expr_is_numeric_quantity(e)) return false;
    Expr* n = numericalize(e, numeric_machine_spec());
    bool ok = false;
    if (n) {
        if (n->type == EXPR_REAL) { *out = n->data.real; ok = isfinite(*out); }
        else if (n->type == EXPR_INTEGER) { *out = (double)n->data.integer; ok = true; }
    }
    if (n) expr_free(n);
    return ok;
}

/* Placeholder for (kind, order, point); interns a new one when unseen. */
static const Expr* placeholder(ZdCtx* c, ZdKind kind, int order, const Expr* point,
                               double pd) {
    for (int i = 0; i < c->n; i++) {
        ZdAtom* a = &c->atoms[i];
        if (a->kind == kind && a->order == order &&
            fabs(a->pointd - pd) <= 1e-12 * (1.0 + fabs(pd)))
            return a->sym;
    }
    if (c->n >= ZD_MAX_ATOMS) { c->fail = true; return NULL; }
    char name[64];
    snprintf(name, sizeof name, "ZeroTest`Private`dist%d", c->n);
    ZdAtom* a = &c->atoms[c->n++];
    a->kind = kind; a->order = order; a->point = expr_copy((Expr*)point); a->pointd = pd;
    a->sym = expr_new_symbol(name);
    return a->sym;
}

/* `e` with every distribution atom replaced by its normalised placeholder form
 * (step 2).  Returns NULL (and sets c->fail) on an atom that cannot be
 * normalised.  Unevaluated. */
static Expr* normalise(const Expr* e, ZdCtx* c) {
    if (c->fail || !e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    ZdKind kind; int order; const Expr* arg;
    if (atom_of(e, &kind, &order, &arg)) {
        if (order > ZD_MAX_ORDER) { c->fail = true; return NULL; }
        const char* v = NULL; bool many = false;
        collect_var(arg, &v, &many);
        if (!v || many || (c->var && c->var != v)) { c->fail = true; return NULL; }
        c->var = v;
        /* arg = s t + r with s, r free of t:  s = D[arg, t], linear iff D2 == 0 */
        Expr* s  = eval_and_free(call2(SYM_D, expr_copy((Expr*)arg), expr_new_symbol(v)));
        Expr* s2 = eval_and_free(call2(SYM_D, expr_copy(s), expr_new_symbol(v)));
        double sd = 0.0;
        bool lin = s2 && s2->type == EXPR_INTEGER && s2->data.integer == 0;
        if (s2) expr_free(s2);
        if (!lin || !real_value(s, &sd) || sd == 0.0) {
            expr_free(s); c->fail = true; return NULL;
        }
        Expr* zero = expr_new_integer(0);
        Expr* r = subst_var(arg, v, zero);
        expr_free(zero);
        /* point a = -r / s */
        Expr* pt = eval_and_free(call2(SYM_Times, expr_new_integer(-1),
                                       call2(SYM_Times, r, call2(SYM_Power, expr_copy(s),
                                                                 expr_new_integer(-1)))));
        double pd = 0.0;
        if (!real_value(pt, &pd)) { expr_free(s); expr_free(pt); c->fail = true; return NULL; }
        const Expr* ph = placeholder(c, kind, order, pt, pd);
        expr_free(pt);
        if (!ph) { expr_free(s); return NULL; }
        Expr* out;
        if (kind == ZD_DELTA) {
            /* d^(k)(s u) = d^(k)(u) / (s^k |s|) */
            Expr* den = call2(SYM_Times,
                              call2(SYM_Power, expr_copy(s), expr_new_integer(order)),
                              call1(SYM_Abs, expr_copy(s)));
            out = call2(SYM_Times, expr_copy((Expr*)ph),
                        call2(SYM_Power, den, expr_new_integer(-1)));
        } else {
            out = sd > 0.0 ? expr_copy((Expr*)ph)
                           : call2(SYM_Plus, expr_new_integer(1),
                                   call2(SYM_Times, expr_new_integer(-1), expr_copy((Expr*)ph)));
        }
        expr_free(s);
        return out;
    }
    size_t n = e->data.function.arg_count;
    Expr* h = normalise(e->data.function.head, c);
    if (!h) return NULL;
    Expr** args = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        args[i] = normalise(e->data.function.args[i], c);
        if (!args[i]) {
            for (size_t j = 0; j < i; j++) expr_free(args[j]);
            free(args); expr_free(h);
            return NULL;
        }
    }
    Expr* r = expr_new_function(h, args, n);
    free(args);
    return r;
}

/* Index of the placeholder symbol `e`, or -1. */
static int atom_index(const ZdCtx* c, const Expr* e) {
    if (!e || e->type != EXPR_SYMBOL) return -1;
    for (int i = 0; i < c->n; i++)
        if (c->atoms[i].sym->data.symbol.name == e->data.symbol.name) return i;
    return -1;
}

/* Does `e` mention any placeholder? */
static bool mentions_atom(const ZdCtx* c, const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return atom_index(c, e) >= 0;
    if (e->type != EXPR_FUNCTION) return false;
    if (mentions_atom(c, e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (mentions_atom(c, e->data.function.args[i])) return true;
    return false;
}

/* Split one expanded term into (delta index or -1, coefficient).  A Theta
 * power collapses to Theta (idempotent).  Returns false when the term is not a
 * distribution we can name (two deltas, a delta power, a placeholder inside a
 * non-polynomial position).  *coef is the term without its delta factor. */
static bool split_term(const ZdCtx* c, const Expr* term, int* delta, Expr** coef) {
    *delta = -1;
    const Expr* const* fs; size_t nf; const Expr* one[1];
    if (head_is_sym(term, SYM_Times)) {
        fs = (const Expr* const*)term->data.function.args; nf = term->data.function.arg_count;
    } else { one[0] = term; fs = one; nf = 1; }
    Expr** out = malloc((nf ? nf : 1) * sizeof(Expr*)); size_t no = 0;
    bool ok = true;
    for (size_t i = 0; i < nf && ok; i++) {
        const Expr* f = fs[i];
        int ai = atom_index(c, f);
        int pw = 1;
        if (ai < 0 && head_is_sym(f, SYM_Power) && f->data.function.arg_count == 2) {
            ai = atom_index(c, f->data.function.args[0]);
            const Expr* ex = f->data.function.args[1];
            if (ai >= 0) {
                if (ex->type != EXPR_INTEGER || ex->data.integer < 1) { ok = false; break; }
                pw = (int)ex->data.integer;
            }
        }
        if (ai >= 0) {
            if (c->atoms[ai].kind == ZD_DELTA) {
                if (pw != 1 || *delta >= 0) { ok = false; break; }
                *delta = ai;
            } else {
                out[no++] = expr_copy(c->atoms[ai].sym);   /* Theta^n == Theta */
            }
            continue;
        }
        if (mentions_atom(c, f)) { ok = false; break; }
        out[no++] = expr_copy((Expr*)f);
    }
    if (!ok) {
        for (size_t i = 0; i < no; i++) expr_free(out[i]);
        free(out);
        return false;
    }
    *coef = no == 0 ? expr_new_integer(1)
          : no == 1 ? out[0]
          : expr_new_function(expr_new_symbol(SYM_Times), out, no);
    if (no <= 1) free(out);           /* >1: the array went to expr_new_function */
    return true;
}

/* coef with each Theta placeholder replaced by its value on the side of the
 * point x0: Theta(t - b) = (x0 > b).  Callers never ask at a breakpoint of a
 * Theta that occurs in coef.  Evaluated result; coef borrowed. */
static Expr* theta_at(const ZdCtx* c, const Expr* coef, double x0) {
    int nr = 0;
    Expr* rules[ZD_MAX_ATOMS];
    for (int i = 0; i < c->n; i++) {
        if (c->atoms[i].kind != ZD_THETA) continue;
        rules[nr++] = rule(expr_copy(c->atoms[i].sym),
                           expr_new_integer(x0 > c->atoms[i].pointd ? 1 : 0));
    }
    if (nr == 0) return eval_and_free(expr_copy((Expr*)coef));
    Expr* lst = expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)nr);
    return eval_and_free(call2(SYM_ReplaceAll, expr_copy((Expr*)coef), lst));
}

/* Does `e` mention the placeholder of atom i? */
static bool mentions_one(const ZdCtx* c, int i, const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return atom_index(c, e) == i;
    if (e->type != EXPR_FUNCTION) return false;
    if (mentions_one(c, i, e->data.function.head)) return true;
    for (size_t k = 0; k < e->data.function.arg_count; k++)
        if (mentions_one(c, i, e->data.function.args[k])) return true;
    return false;
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

/* Fold a verdict into the running one: FALSE dominates, then UNKNOWN. */
static void fold(ZeroTestResult* acc, ZeroTestResult r) {
    if (*acc == ZERO_TEST_FALSE) return;
    if (r == ZERO_TEST_FALSE) *acc = ZERO_TEST_FALSE;
    else if (r == ZERO_TEST_UNKNOWN) *acc = ZERO_TEST_UNKNOWN;
}

ZeroTestResult zt_decide_distribution(const Expr* e,
                                      ZeroTestResult (*decide)(const Expr*)) {
    if (!zt_has_distribution(e)) return ZERO_TEST_UNKNOWN;
    if (leaf_count_internal((Expr*)e, true) > ZD_MAX_LEAVES) return ZERO_TEST_UNKNOWN;
    ZdCtx c; memset(&c, 0, sizeof c);
    ZeroTestResult verdict = ZERO_TEST_UNKNOWN;
    Expr* norm = normalise(e, &c);
    Expr* ex = NULL;
    Expr** coefs = NULL;      /* accumulated delta-basis coefficient, per atom */
    Expr*  regular = NULL;
    if (!norm || c.fail || !c.var) goto out;
    norm = eval_and_free(norm);
    ex = expr_expand(norm);
    if (!ex) goto out;

    coefs = calloc((size_t)c.n, sizeof(Expr*));
    /* The basis d^(m)(t - a) shares a placeholder only with the same (a, m);
     * sifting can create a lower order at the same point, so make sure every
     * order 0..k exists for each delta point first. */
    for (int i = 0, n0 = c.n; i < n0; i++) {
        if (c.atoms[i].kind != ZD_DELTA) continue;
        for (int m = 0; m < c.atoms[i].order; m++)
            if (!placeholder(&c, ZD_DELTA, m, c.atoms[i].point, c.atoms[i].pointd)) goto out;
    }
    free(coefs);
    coefs = calloc((size_t)c.n, sizeof(Expr*));

    {
        const Expr* const* ts; size_t nt; const Expr* one[1];
        if (head_is_sym(ex, SYM_Plus)) {
            ts = (const Expr* const*)ex->data.function.args; nt = ex->data.function.arg_count;
        } else { one[0] = ex; ts = one; nt = 1; }
        Expr** reg = malloc((nt ? nt : 1) * sizeof(Expr*)); size_t nreg = 0;
        Expr* tvar = expr_new_symbol(c.var);
        bool ok = true;
        for (size_t i = 0; i < nt && ok; i++) {
            int di; Expr* cf = NULL;
            if (!split_term(&c, ts[i], &di, &cf)) { ok = false; break; }
            if (di < 0) { reg[nreg++] = cf; continue; }
            ZdAtom* d = &c.atoms[di];
            /* a Theta at the delta's own point cannot be sifted */
            bool bad = false;
            for (int k = 0; k < c.n; k++)
                if (c.atoms[k].kind == ZD_THETA &&
                    fabs(c.atoms[k].pointd - d->pointd) <= 1e-12 * (1.0 + fabs(d->pointd)) &&
                    mentions_one(&c, k, cf)) bad = true;
            if (bad) { expr_free(cf); ok = false; break; }
            Expr* cf0 = theta_at(&c, cf, d->pointd);
            expr_free(cf);
            /* sift: c(t) d^(k)(t-a) = sum_j (-1)^j C(k,j) c^(j)(a) d^(k-j)(t-a) */
            Expr* dj = cf0;   /* c^(j)(t), j = 0 */
            int64_t binom = 1;
            for (int j = 0; j <= d->order && ok; j++) {
                if (j > 0) {
                    Expr* nx = eval_and_free(call2(SYM_D, dj, expr_copy(tvar)));
                    dj = nx;
                    binom = binom * (d->order - j + 1) / j;
                }
                Expr* at = subst_var(dj, c.var, d->point);
                int64_t sgn = (j % 2) ? -1 : 1;
                Expr* contrib = call2(SYM_Times, expr_new_integer(sgn * binom), at);
                /* placeholder of (a, k - j) */
                const Expr* ph = placeholder(&c, ZD_DELTA, d->order - j, d->point, d->pointd);
                int pi = atom_index(&c, ph);
                if (pi < 0 || pi >= c.n) { expr_free(contrib); ok = false; break; }
                coefs[pi] = coefs[pi] ? call2(SYM_Plus, coefs[pi], contrib) : contrib;
            }
            expr_free(dj);
        }
        expr_free(tvar);
        if (!ok) {
            for (size_t i = 0; i < nreg; i++) expr_free(reg[i]);
            free(reg);
            goto out;
        }
        regular = nreg == 0 ? expr_new_integer(0)
                : expr_new_function(expr_new_symbol(SYM_Plus), reg, nreg);
        free(reg);
    }

    verdict = ZERO_TEST_TRUE;
    /* step 3: every delta-basis coefficient is zero */
    for (int i = 0; i < c.n && verdict != ZERO_TEST_FALSE; i++) {
        if (!coefs[i]) continue;
        Expr* v = eval_and_free(coefs[i]);
        coefs[i] = NULL;
        fold(&verdict, decide(v));
        expr_free(v);
    }
    /* step 4: the regular part is zero on every interval between breakpoints */
    if (verdict != ZERO_TEST_FALSE) {
        double bp[ZD_MAX_ATOMS]; int nb = 0;
        for (int i = 0; i < c.n; i++)
            if (c.atoms[i].kind == ZD_THETA) bp[nb++] = c.atoms[i].pointd;
        qsort(bp, (size_t)nb, sizeof(double), cmp_double);
        /* one representative point per interval: below the first breakpoint,
         * between consecutive ones, above the last */
        for (int j = 0; j <= nb && verdict != ZERO_TEST_FALSE; j++) {
            double x0 = nb == 0 ? 0.0
                      : j == 0  ? bp[0] - 1.0
                      : j == nb ? bp[nb - 1] + 1.0
                      : 0.5 * (bp[j - 1] + bp[j]);
            if (j > 0 && j < nb && bp[j] - bp[j - 1] <= 1e-12 * (1.0 + fabs(bp[j]))) continue;
            Expr* piece = theta_at(&c, regular, x0);
            fold(&verdict, decide(piece));
            expr_free(piece);
        }
    }

out:
    if (coefs) {
        for (int i = 0; i < c.n; i++) if (coefs[i]) expr_free(coefs[i]);
        free(coefs);
    }
    if (regular) expr_free(regular);
    if (ex) expr_free(ex);
    if (norm) expr_free(norm);
    for (int i = 0; i < c.n; i++) {
        expr_free(c.atoms[i].point);
        expr_free(c.atoms[i].sym);
    }
    return verdict;
}
