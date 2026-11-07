/* Mathilda — NDSolve events: WhenEvent[event, action] (see ndsolve_events.h). */
#include "ndsolve_events.h"
#include "../sym_names.h"
#include "../sym_intern.h"
#include "../eval.h"
#include "../common.h"
#include "../arithmetic.h"
#include "../numeric.h"
#include "../compile/autocompile.h"
#include "../message.h"
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

/* ------------------------------------------------------------------ *
 *  Representation                                                     *
 * ------------------------------------------------------------------ */
typedef enum {
    EVK_ZERO,   /* lhs == rhs: g = lhs - rhs changes sign (either direction)   */
    EVK_RISE,   /* inequality: g > 0 (strict) or g >= 0 becomes true           */
    EVK_BOOL    /* any other condition: becomes True (located by bisection)    */
} EvKind;

typedef struct {
    EvKind  kind;
    bool    strict;        /* EVK_RISE: true iff g > 0 (else g >= 0)           */
    Expr*   g;             /* event function in the reduced state (owned)      */
    Expr*   cond;          /* EVK_BOOL condition (owned)                       */
    bool    stop, remove;  /* "StopIntegration" / "RemoveEvent"                */
    size_t  nreset;
    size_t* ridx;          /* reduced index of each reset target               */
    Expr**  rrhs;          /* reset right-hand sides (owned)                   */
    size_t  nact;
    Expr**  act;           /* expression actions (owned, held)                 */
    bool    active;        /* armed in the current direction                   */
    bool    known;         /* the stored node value below is valid             */
    double  gp;            /* g at the current node (ZERO/RISE)                */
    int     tp;            /* truth at the current node (RISE/BOOL)            */
    double  gscale;        /* |g| across the last located step (ZERO)          */
    double  groot;         /* g at the located root (ZERO)                     */
    AutoCompiled* ac;      /* machine compiled g, or NULL                      */
    bool    ac_tried;
#ifdef USE_MPFR
    AutoCompiled* acp;     /* MPFR compiled g at acp_bits, or NULL             */
    bool    acp_tried;
    long    acp_bits;
    mpfr_t  gp_mp, groot_mp, gscale_mp;
    bool    mp_init;
#endif
} NdEvent;

struct NdEvents {
    size_t   n;
    NdEvent* ev;
};

bool nd_is_when_event(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol.name, "WhenEvent") == 0;
}

void nd_events_free(struct NdEvents* E) {
    if (!E) return;
    for (size_t i = 0; i < E->n; i++) {
        NdEvent* v = &E->ev[i];
        expr_free(v->g); expr_free(v->cond);
        for (size_t j = 0; j < v->nreset; j++) expr_free(v->rrhs[j]);
        free(v->rrhs); free(v->ridx);
        for (size_t j = 0; j < v->nact; j++) expr_free(v->act[j]);
        free(v->act);
        if (v->ac) autocompiled_free(v->ac);
#ifdef USE_MPFR
        if (v->acp) autocompiled_free(v->acp);
        if (v->mp_init) { mpfr_clear(v->gp_mp); mpfr_clear(v->groot_mp); mpfr_clear(v->gscale_mp); }
#endif
    }
    free(E->ev);
    free(E);
}

/* ------------------------------------------------------------------ *
 *  Construction                                                       *
 * ------------------------------------------------------------------ */

/* Structural (non-evaluating) substitution: every subexpression equal to
 * lits[i] is replaced by a copy of subs[i].  WhenEvent holds its arguments, and
 * an action such as Sow[t] must not run while it is being rewritten. */
static Expr* ev_subst(const Expr* e, Expr** lits, Expr** subs, size_t n) {
    for (size_t i = 0; i < n; i++) if (expr_eq(e, lits[i])) return expr_copy(subs[i]);
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t na = e->data.function.arg_count;
    Expr* head = ev_subst(e->data.function.head, lits, subs, n);
    Expr** args = malloc(sizeof(Expr*) * (na ? na : 1));
    for (size_t i = 0; i < na; i++) args[i] = ev_subst(e->data.function.args[i], lits, subs, n);
    Expr* r = expr_new_function(head, args, na);
    free(args);
    return r;
}

static bool ev_head(const Expr* e, const char* name, size_t nargs) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name && e->data.function.arg_count == nargs;
}

/* Match  u[t]  or  Derivative[m][u][t]  with t the independent variable. */
static bool ev_match_state(const Expr* e, const char* tvar, const char** fname, int* order) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 1) return false;
    const Expr* a = e->data.function.args[0];
    if (a->type != EXPR_SYMBOL || a->data.symbol.name != tvar) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL) { *fname = h->data.symbol.name; *order = 0; return true; }
    if (h->type == EXPR_FUNCTION && h->data.function.arg_count == 1
        && h->data.function.args[0]->type == EXPR_SYMBOL
        && ev_head(h->data.function.head, SYM_Derivative, 1)
        && h->data.function.head->data.function.args[0]->type == EXPR_INTEGER) {
        *fname = h->data.function.args[0]->data.symbol.name;
        *order = (int)h->data.function.head->data.function.args[0]->data.integer;
        return true;
    }
    return false;
}

typedef struct {
    NdProblem* P;
    Expr** lits; Expr** subs; size_t nsub;
    const char** fnames; const size_t* base; const int* maxorder; size_t nfun;
} EvBuildCtx;

static void ev_warn(const char* tag, const char* msg) { mth_message("NDSolve", tag, "%s", msg); }

/* Parse one action (recursively for a List) into v. */
static bool ev_parse_action(const EvBuildCtx* B, NdEvent* v, const Expr* a) {
    if (a->type == EXPR_STRING) {
        const char* s = a->data.string;
        if (strcmp(s, "StopIntegration") == 0) { v->stop = true; return true; }
        if (strcmp(s, "RemoveEvent") == 0) { v->remove = true; return true; }
        if (strcmp(s, "CrossDiscontinuity") == 0 || strcmp(s, "RestartIntegration") == 0)
            return true;   /* every event already restarts the step at the event */
        mth_message("NDSolve", "wevact", "unknown WhenEvent action \"%s\"", s);
        return false;
    }
    if (ev_head(a, SYM_List, a->type == EXPR_FUNCTION ? a->data.function.arg_count : 0)) {
        for (size_t i = 0; i < a->data.function.arg_count; i++)
            if (!ev_parse_action(B, v, a->data.function.args[i])) return false;
        return true;
    }
    if (ev_head(a, SYM_Rule, 2) || ev_head(a, SYM_RuleDelayed, 2)) {
        const char* fn; int ord;
        const Expr* lhs = a->data.function.args[0];
        if (ev_match_state(lhs, B->P->tvar, &fn, &ord)) {
            for (size_t k = 0; k < B->nfun; k++) if (B->fnames[k] == fn) {
                if (ord < 0 || ord >= B->maxorder[k]) break;
                v->ridx = realloc(v->ridx, sizeof(size_t) * (v->nreset + 1));
                v->rrhs = realloc(v->rrhs, sizeof(Expr*) * (v->nreset + 1));
                v->ridx[v->nreset] = B->base[k] + (size_t)ord;
                v->rrhs[v->nreset] = ev_subst(a->data.function.args[1], B->lits, B->subs, B->nsub);
                v->nreset++;
                return true;
            }
        }
        ev_warn("wevrs", "a WhenEvent state reset must have the form u[t] -> value or "
                         "u'[t] -> value for a dependent function u below its highest order");
        return false;
    }
    v->act = realloc(v->act, sizeof(Expr*) * (v->nact + 1));
    v->act[v->nact++] = ev_subst(a, B->lits, B->subs, B->nsub);
    return true;
}

/* Add one event condition (a List of conditions adds one event per element,
 * sharing the actions). */
static bool ev_add(const EvBuildCtx* B, struct NdEvents* E, const Expr* cond, const Expr* action) {
    if (ev_head(cond, SYM_List, cond->type == EXPR_FUNCTION ? cond->data.function.arg_count : 0)) {
        for (size_t i = 0; i < cond->data.function.arg_count; i++)
            if (!ev_add(B, E, cond->data.function.args[i], action)) return false;
        return true;
    }
    E->ev = realloc(E->ev, sizeof(NdEvent) * (E->n + 1));
    NdEvent* v = &E->ev[E->n++];
    memset(v, 0, sizeof(*v));
    Expr* c = ev_subst(cond, B->lits, B->subs, B->nsub);
    const char* h = (c->type == EXPR_FUNCTION && c->data.function.head->type == EXPR_SYMBOL)
                    ? c->data.function.head->data.symbol.name : NULL;
    bool two = c->type == EXPR_FUNCTION && c->data.function.arg_count == 2;
    if (two && (h == SYM_Equal || h == SYM_Greater || h == SYM_GreaterEqual
                || h == SYM_Less || h == SYM_LessEqual)) {
        bool flip = (h == SYM_Less || h == SYM_LessEqual);
        Expr* l = expr_copy(c->data.function.args[flip ? 1 : 0]);
        Expr* r = expr_copy(c->data.function.args[flip ? 0 : 1]);
        v->g = eval_and_free(nd_call2(SYM_Subtract, l, r));
        v->kind = (h == SYM_Equal) ? EVK_ZERO : EVK_RISE;
        v->strict = (h == SYM_Greater || h == SYM_Less);
        expr_free(c);
    } else {
        v->kind = EVK_BOOL;
        v->cond = c;
    }
    return ev_parse_action(B, v, action);
}

bool nd_events_build(NdProblem* P, Expr** items, size_t nitems,
                     Expr** lits, Expr** subs, size_t nsub,
                     const char** fnames, const size_t* base,
                     const int* maxorder, size_t nfun) {
    if (nitems == 0) return true;
    struct NdEvents* E = calloc(1, sizeof(*E));
    EvBuildCtx B = { P, lits, subs, nsub, fnames, base, maxorder, nfun };
    bool ok = true;
    for (size_t i = 0; i < nitems && ok; i++) {
        const Expr* w = items[i];
        if (w->data.function.arg_count < 2) {
            ev_warn("wevarg", "WhenEvent needs an event and an action: WhenEvent[event, action]");
            ok = false; break;
        }
        ok = ev_add(&B, E, w->data.function.args[0], w->data.function.args[1]);
    }
    if (!ok) { nd_events_free(E); return false; }
    P->events = E;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Machine precision                                                  *
 * ------------------------------------------------------------------ */

/* g (ZERO/RISE) or the truth (RISE/BOOL) at (t, Y).  Returns false when the
 * condition does not evaluate to a real number / a Boolean there. */
static bool ev_value(NdProblem* P, NdEvent* v, double t, const double* Y, double* g, int* truth) {
    size_t d = P->d;
    if (v->kind == EVK_BOOL) {
        nd_bind_point(P, t, Y);
        arith_warnings_mute_push();
        Expr* r = eval_and_free(expr_copy(v->cond));
        arith_warnings_mute_pop();
        bool ok = r && r->type == EXPR_SYMBOL
                  && (r->data.symbol.name == SYM_True || r->data.symbol.name == SYM_False);
        if (ok) *truth = r->data.symbol.name == SYM_True;
        expr_free(r);
        return ok;
    }
    bool ok = false;
    if (!v->ac_tried) {
        v->ac_tried = true;
        if (!P->eval_monitor) {
            const Expr** vars = malloc(sizeof(Expr*) * (d + 1));
            for (size_t i = 0; i < d; i++) vars[i] = P->ysym[i];
            Expr* ts = expr_new_symbol(P->tvar);
            vars[d] = ts;
            v->ac = autocompile_new(v->g, vars, d + 1);
            expr_free(ts); free(vars);
        }
    }
    if (v->ac) {
        double* xs = malloc(sizeof(double) * (d + 1));
        memcpy(xs, Y, sizeof(double) * d);
        xs[d] = t;
        ok = autocompiled_eval_real(v->ac, xs, g);
        free(xs);
    }
    if (!ok) {
        nd_bind_point(P, t, Y);
        ok = nd_eval_bound(P, v->g, g);
    }
    if (ok) *truth = v->strict ? (*g > 0.0) : (*g >= 0.0);
    return ok;
}

/* Has the event fired at a point with value (g, truth), relative to the
 * stored node value?  `sb` is the sign of g on the crossed side (ZERO). */
static bool ev_fired_side(const NdEvent* v, double g, int truth, int sb) {
    if (v->kind == EVK_ZERO) return g == 0.0 || (g > 0.0 ? 1 : -1) == sb;
    return truth != 0;
}

void nd_events_begin(NdProblem* P, double t0, const double* Y0) {
    struct NdEvents* E = P->events;
    if (!E) return;
    for (size_t i = 0; i < E->n; i++) {
        NdEvent* v = &E->ev[i];
        v->active = true;
        v->known = ev_value(P, v, t0, Y0, &v->gp, &v->tp);
        if (!v->known)
            mth_message("NDSolve", "evnum",
                        "the event condition of WhenEvent did not evaluate to %s at the initial point",
                        v->kind == EVK_BOOL ? "True or False" : "a real number");
    }
}

/* Locate the crossing of event v in the step from a (not fired, value fa) to
 * b (fired, value fb) by a safeguarded Illinois iteration; returns the fired
 * end of the final bracket.  Ys is scratch for the state. */
static double ev_solve(NdProblem* P, NdEvent* v, double a, double fa, double b, double fb,
                       int sb, NdEvAdvance adv, void* ctx, double* Ys, int* nadv) {
    int last = 0;               /* +1: b moved last, -1: a moved last         */
    double w4 = fabs(b - a);    /* bracket width four iterations ago          */
    bool bisect = (v->kind == EVK_BOOL);
    bool force_mid = false;     /* slow progress: take one bisection         */
    for (int it = 0; it < 400; it++) {
        double tol = 2.0 * DBL_EPSILON * fmax(fabs(a), fabs(b));
        if (!(fabs(b - a) > tol)) break;
        double c;
        if (bisect || force_mid || fb == fa) c = 0.5 * (a + b);
        else {
            c = b - fb * (b - a) / (fb - fa);
            if (!isfinite(c) || !((c - a) * (c - b) < 0.0)) c = 0.5 * (a + b);
        }
        if (c == a || c == b) { c = 0.5 * (a + b); if (c == a || c == b) break; }
        double gc = 0.0; int tc = 0;
        (*nadv)++;
        if (!adv(ctx, c, Ys) || !ev_value(P, v, c, Ys, &gc, &tc)) {
            /* unevaluable inside the step: fall back to bisection on b */
            b = c; fb = 0.5 * fb; bisect = true; continue;
        }
        if (v->kind == EVK_ZERO && gc == 0.0) { b = c; break; }
        if (ev_fired_side(v, gc, tc, sb)) { b = c; fb = gc; if (last == 1) fa *= 0.5; last = 1; }
        else                               { a = c; fa = gc; if (last == -1) fb *= 0.5; last = -1; }
        force_mid = false;
        if (!bisect && (it & 3) == 3) {
            if (fabs(b - a) > 0.5 * w4) force_mid = true;
            w4 = fabs(b - a);
        }
    }
    return b;
}

int nd_events_locate(NdProblem* P, double ta, const double* Ya,
                     double tb, const double* Yb,
                     NdEvAdvance adv, void* ctx, double* te, double* Ye) {
    struct NdEvents* E = P->events;
    if (!E) return -1;
    (void)Ya;
    size_t d = P->d;
    double* gb = malloc(sizeof(double) * E->n);
    int* tbv = malloc(sizeof(int) * E->n);
    bool* okb = malloc(sizeof(bool) * E->n);
    double* Ys = malloc(sizeof(double) * d);
    int best = -1; double tbest = tb; int nadv = 0;
    double dir = tb > ta ? 1.0 : -1.0;
    for (size_t i = 0; i < E->n; i++) {
        NdEvent* v = &E->ev[i];
        okb[i] = v->active && ev_value(P, v, tb, Yb, &gb[i], &tbv[i]);
    }
    for (size_t i = 0; i < E->n; i++) {
        NdEvent* v = &E->ev[i];
        if (!okb[i] || !v->known) continue;
        bool fire;
        int sb = 0;
        if (v->kind == EVK_ZERO) {
            sb = gb[i] > 0.0 ? 1 : (gb[i] < 0.0 ? -1 : 0);
            fire = v->gp != 0.0 && (sb == 0 || (v->gp > 0.0) != (sb > 0));
            if (sb == 0) sb = v->gp > 0.0 ? -1 : 1;
        } else fire = !v->tp && tbv[i];
        if (!fire) continue;
        double t;
        if (v->kind == EVK_ZERO && gb[i] == 0.0) t = tb;
        else {
            v->gscale = fabs(v->gp) + fabs(gb[i]);
            t = ev_solve(P, v, ta, v->gp, tb, gb[i], sb, adv, ctx, Ys, &nadv);
        }
        if (best < 0 || (t - tbest) * dir < 0.0) { best = (int)i; tbest = t; }
    }
    if (best >= 0) {
        *te = tbest;
        if (tbest == tb && nadv == 0) memcpy(Ye, Yb, sizeof(double) * d);
        else if (!adv(ctx, tbest, Ye)) { best = -1; }
        if (best >= 0) {
            NdEvent* v = &E->ev[best];
            int tr;
            if (v->kind == EVK_ZERO && !ev_value(P, v, tbest, Ye, &v->groot, &tr)) v->groot = 0.0;
        }
    }
    if (best < 0) {
        for (size_t i = 0; i < E->n; i++) {
            NdEvent* v = &E->ev[i];
            if (!v->active) continue;
            v->known = okb[i];
            if (okb[i]) { v->gp = gb[i]; v->tp = tbv[i]; }
        }
    }
    free(gb); free(tbv); free(okb); free(Ys);
    return best;
}

/* Run the expression actions with t and the state bound to the event point. */
static void ev_run_actions(NdEvent* v) {
    for (size_t j = 0; j < v->nact; j++) {
        Expr* r = eval_and_free(expr_copy(v->act[j]));
        expr_free(r);
    }
}

unsigned nd_events_fire(NdProblem* P, int k, double te, const double* Ye, double* Ynew) {
    struct NdEvents* E = P->events;
    size_t d = P->d;
    memcpy(Ynew, Ye, sizeof(double) * d);
    if (!E || k < 0 || (size_t)k >= E->n) return 0;
    NdEvent* v = &E->ev[k];
    unsigned flags = 0;
    nd_bind_point(P, te, Ye);
    ev_run_actions(v);
    if (v->nreset) {
        /* every right-hand side sees the pre-event state */
        double* val = malloc(sizeof(double) * v->nreset);
        bool ok = true;
        nd_bind_point(P, te, Ye);
        for (size_t j = 0; j < v->nreset && ok; j++) ok = nd_eval_bound(P, v->rrhs[j], &val[j]);
        if (ok) { for (size_t j = 0; j < v->nreset; j++) Ynew[v->ridx[j]] = val[j]; flags |= ND_EVF_RESET; }
        else ev_warn("evrs", "a WhenEvent state reset did not evaluate to a real number; the state is unchanged");
        free(val);
    }
    if (v->stop) flags |= ND_EVF_STOP;
    if (v->remove) v->active = false;
    /* new node values at the post-event state */
    for (size_t i = 0; i < E->n; i++) {
        NdEvent* w = &E->ev[i];
        if (!w->active) continue;
        double g; int tr;
        w->known = ev_value(P, w, te, Ynew, &g, &tr);
        if (!w->known) continue;
        w->gp = g; w->tp = tr;
        /* The fired zero-crossing event sits ON its surface: unless a reset
         * moved it off, its value there is rounding noise whose sign must not
         * register as a second crossing on the next step. */
        if ((int)i == k && w->kind == EVK_ZERO
            && fabs(g - w->groot) <= 1e-6 * w->gscale + DBL_MIN)
            w->gp = 0.0;
    }
    return flags;
}

/* ------------------------------------------------------------------ *
 *  Arbitrary precision (MPFR)                                         *
 * ------------------------------------------------------------------ */
#ifdef USE_MPFR

static void ev_bind_mp(NdProblem* P, const mpfr_t t, mpfr_t* Y) {
    nd_bind_set(&P->bind_t, expr_new_mpfr_copy(t));
    for (size_t i = 0; i < P->d; i++) nd_bind_set(&P->bind_y[i], expr_new_mpfr_copy(Y[i]));
    eval_clock_bump();
}

/* Evaluate a bound expression to a real mpfr value at `bits`. */
static bool ev_eval_mp(Expr* e, mpfr_t out, long bits) {
    arith_warnings_mute_push();
    Expr* raw = eval_and_free(expr_copy(e));
    arith_warnings_mute_pop();
    if (!raw) return false;
    NumericSpec spec; spec.mode = NUMERIC_MODE_MPFR; spec.bits = bits;
    Expr* num = numericalize(raw, spec);
    expr_free(raw);
    mpfr_t im; mpfr_init2(im, bits);
    bool inexact;
    bool ok = num && get_approx_mpfr(num, out, im, &inexact) && mpfr_number_p(out) && mpfr_zero_p(im);
    mpfr_clear(im);
    expr_free(num);
    return ok;
}

static bool ev_value_mp(NdProblem* P, NdEvent* v, const mpfr_t t, mpfr_t* Y, long bits,
                        mpfr_t g, int* truth) {
    size_t d = P->d;
    if (v->kind == EVK_BOOL) {
        ev_bind_mp(P, t, Y);
        arith_warnings_mute_push();
        Expr* r = eval_and_free(expr_copy(v->cond));
        arith_warnings_mute_pop();
        bool ok = r && r->type == EXPR_SYMBOL
                  && (r->data.symbol.name == SYM_True || r->data.symbol.name == SYM_False);
        if (ok) *truth = r->data.symbol.name == SYM_True;
        expr_free(r);
        return ok;
    }
    if (!v->acp_tried || v->acp_bits != bits) {
        if (v->acp) { autocompiled_free(v->acp); v->acp = NULL; }
        v->acp_tried = true; v->acp_bits = bits;
        if (!P->eval_monitor) {
            const Expr** vars = malloc(sizeof(Expr*) * (d + 1));
            for (size_t i = 0; i < d; i++) vars[i] = P->ysym[i];
            Expr* ts = expr_new_symbol(P->tvar);
            vars[d] = ts;
            v->acp = autocompile_new_prec(v->g, vars, d + 1, bits);
            expr_free(ts); free(vars);
        }
    }
    bool ok = false;
    if (v->acp) {
        Expr** xs = malloc(sizeof(Expr*) * (d + 1));
        for (size_t i = 0; i < d; i++) xs[i] = expr_new_mpfr_copy(Y[i]);
        xs[d] = expr_new_mpfr_copy(t);
        Expr* r = autocompiled_eval_mpfr(v->acp, (const Expr* const*)xs);
        mpfr_t im; mpfr_init2(im, bits);
        bool inexact;
        ok = r && get_approx_mpfr(r, g, im, &inexact) && mpfr_number_p(g) && mpfr_zero_p(im);
        mpfr_clear(im);
        expr_free(r);
        for (size_t i = 0; i <= d; i++) expr_free(xs[i]);
        free(xs);
    }
    if (!ok) { ev_bind_mp(P, t, Y); ok = ev_eval_mp(v->g, g, bits); }
    if (ok) *truth = v->strict ? (mpfr_sgn(g) > 0) : (mpfr_sgn(g) >= 0);
    return ok;
}

static void ev_mp_ensure(NdEvent* v, long bits) {
    if (v->mp_init && mpfr_get_prec(v->gp_mp) == (mpfr_prec_t)bits) return;
    if (v->mp_init) { mpfr_clear(v->gp_mp); mpfr_clear(v->groot_mp); mpfr_clear(v->gscale_mp); }
    mpfr_init2(v->gp_mp, bits); mpfr_init2(v->groot_mp, bits); mpfr_init2(v->gscale_mp, bits);
    mpfr_set_zero(v->gp_mp, 1); mpfr_set_zero(v->groot_mp, 1); mpfr_set_zero(v->gscale_mp, 1);
    v->mp_init = true;
}

void nd_events_begin_mp(NdProblem* P, const mpfr_t t0, mpfr_t* Y0, long bits) {
    struct NdEvents* E = P->events;
    if (!E) return;
    for (size_t i = 0; i < E->n; i++) {
        NdEvent* v = &E->ev[i];
        ev_mp_ensure(v, bits);
        v->active = true;
        v->known = ev_value_mp(P, v, t0, Y0, bits, v->gp_mp, &v->tp);
        if (!v->known)
            mth_message("NDSolve", "evnum",
                        "the event condition of WhenEvent did not evaluate to %s at the initial point",
                        v->kind == EVK_BOOL ? "True or False" : "a real number");
    }
}

static bool ev_fired_side_mp(const NdEvent* v, const mpfr_t g, int truth, int sb) {
    if (v->kind == EVK_ZERO) return mpfr_zero_p(g) || (mpfr_sgn(g) > 0 ? 1 : -1) == sb;
    return truth != 0;
}

/* MPFR twin of ev_solve; the located (fired-side) time is written to b. */
static void ev_solve_mp(NdProblem* P, NdEvent* v, mpfr_t a, mpfr_t fa, mpfr_t b, mpfr_t fb,
                        int sb, NdEvAdvanceMp adv, void* ctx, mpfr_t* Ys, long bits, long ubits,
                        int* nadv) {
    int last = 0;
    bool bisect = (v->kind == EVK_BOOL);
    mpfr_t c, gc, w, tol, num, den, w4;
    mpfr_inits2(bits, c, gc, w, tol, num, den, w4, (mpfr_ptr)0);
    mpfr_sub(w4, b, a, MPFR_RNDN); mpfr_abs(w4, w4, MPFR_RNDN);
    int maxit = 8 * (int)ubits + 64;
    for (int it = 0; it < maxit; it++) {
        mpfr_sub(w, b, a, MPFR_RNDN); mpfr_abs(w, w, MPFR_RNDN);
        if (mpfr_cmpabs(a, b) > 0) mpfr_abs(tol, a, MPFR_RNDN); else mpfr_abs(tol, b, MPFR_RNDN);
        mpfr_mul_2si(tol, tol, -(long)ubits + 1, MPFR_RNDN);
        if (mpfr_cmp(w, tol) <= 0) break;
        bool mid = bisect || mpfr_equal_p(fb, fa);
        if (!mid) {
            mpfr_sub(den, fb, fa, MPFR_RNDN);
            mpfr_sub(num, b, a, MPFR_RNDN);
            mpfr_mul(num, num, fb, MPFR_RNDN);
            mpfr_div(num, num, den, MPFR_RNDN);
            mpfr_sub(c, b, num, MPFR_RNDN);
            /* strictly inside (a, b)? */
            bool inside = mpfr_number_p(c)
                && ((mpfr_cmp(c, a) > 0 && mpfr_cmp(c, b) < 0) || (mpfr_cmp(c, a) < 0 && mpfr_cmp(c, b) > 0));
            if (!inside) mid = true;
        }
        if (mid) { mpfr_add(c, a, b, MPFR_RNDN); mpfr_div_2ui(c, c, 1, MPFR_RNDN); }
        if (mpfr_equal_p(c, a) || mpfr_equal_p(c, b)) break;
        int tc = 0;
        (*nadv)++;
        if (!adv(ctx, c, Ys) || !ev_value_mp(P, v, c, Ys, bits, gc, &tc)) {
            mpfr_set(b, c, MPFR_RNDN); mpfr_div_2ui(fb, fb, 1, MPFR_RNDN); bisect = true; continue;
        }
        if (v->kind == EVK_ZERO && mpfr_zero_p(gc)) { mpfr_set(b, c, MPFR_RNDN); break; }
        if (ev_fired_side_mp(v, gc, tc, sb)) {
            mpfr_set(b, c, MPFR_RNDN); mpfr_set(fb, gc, MPFR_RNDN);
            if (last == 1) mpfr_div_2ui(fa, fa, 1, MPFR_RNDN);
            last = 1;
        } else {
            mpfr_set(a, c, MPFR_RNDN); mpfr_set(fa, gc, MPFR_RNDN);
            if (last == -1) mpfr_div_2ui(fb, fb, 1, MPFR_RNDN);
            last = -1;
        }
        if (!bisect && (it & 3) == 3) {
            mpfr_sub(w, b, a, MPFR_RNDN); mpfr_abs(w, w, MPFR_RNDN);
            mpfr_div_2ui(w4, w4, 1, MPFR_RNDN);
            if (mpfr_cmp(w, w4) > 0) { /* slow progress: one bisection */
                mpfr_add(c, a, b, MPFR_RNDN); mpfr_div_2ui(c, c, 1, MPFR_RNDN);
                (*nadv)++;
                if (adv(ctx, c, Ys) && ev_value_mp(P, v, c, Ys, bits, gc, &tc)) {
                    if (ev_fired_side_mp(v, gc, tc, sb)) { mpfr_set(b, c, MPFR_RNDN); mpfr_set(fb, gc, MPFR_RNDN); }
                    else { mpfr_set(a, c, MPFR_RNDN); mpfr_set(fa, gc, MPFR_RNDN); }
                }
                mpfr_sub(w, b, a, MPFR_RNDN); mpfr_abs(w, w, MPFR_RNDN);
            }
            mpfr_set(w4, w, MPFR_RNDN);
        }
    }
    mpfr_clears(c, gc, w, tol, num, den, w4, (mpfr_ptr)0);
}

int nd_events_locate_mp(NdProblem* P, const mpfr_t ta, mpfr_t* Ya,
                        const mpfr_t tb, mpfr_t* Yb,
                        NdEvAdvanceMp adv, void* ctx, mpfr_t te, mpfr_t* Ye,
                        long bits, long ubits) {
    struct NdEvents* E = P->events;
    if (!E) return -1;
    (void)Ya;
    size_t d = P->d;
    mpfr_t* gb = malloc(sizeof(mpfr_t) * E->n);
    int* tbv = malloc(sizeof(int) * E->n);
    bool* okb = malloc(sizeof(bool) * E->n);
    mpfr_t* Ys = malloc(sizeof(mpfr_t) * d);
    for (size_t i = 0; i < d; i++) mpfr_init2(Ys[i], bits);
    mpfr_t a, fa, b, fb, tbest;
    mpfr_inits2(bits, a, fa, b, fb, tbest, (mpfr_ptr)0);
    mpfr_set(tbest, tb, MPFR_RNDN);
    int best = -1, nadv = 0;
    int dir = mpfr_cmp(tb, ta) > 0 ? 1 : -1;
    for (size_t i = 0; i < E->n; i++) {
        NdEvent* v = &E->ev[i];
        ev_mp_ensure(v, bits);
        mpfr_init2(gb[i], bits);
        okb[i] = v->active && ev_value_mp(P, v, tb, Yb, bits, gb[i], &tbv[i]);
    }
    for (size_t i = 0; i < E->n; i++) {
        NdEvent* v = &E->ev[i];
        if (!okb[i] || !v->known) continue;
        bool fire;
        int sb = 0;
        if (v->kind == EVK_ZERO) {
            sb = mpfr_sgn(gb[i]);
            int sp = mpfr_sgn(v->gp_mp);
            fire = sp != 0 && (sb == 0 || (sp > 0) != (sb > 0));
            if (sb == 0) sb = sp > 0 ? -1 : 1;
        } else fire = !v->tp && tbv[i];
        if (!fire) continue;
        if (v->kind == EVK_ZERO && mpfr_zero_p(gb[i])) mpfr_set(b, tb, MPFR_RNDN);
        else {
            mpfr_abs(v->gscale_mp, v->gp_mp, MPFR_RNDN);
            mpfr_abs(fa, gb[i], MPFR_RNDN);
            mpfr_add(v->gscale_mp, v->gscale_mp, fa, MPFR_RNDN);
            mpfr_set(a, ta, MPFR_RNDN); mpfr_set(fa, v->gp_mp, MPFR_RNDN);
            mpfr_set(b, tb, MPFR_RNDN); mpfr_set(fb, gb[i], MPFR_RNDN);
            ev_solve_mp(P, v, a, fa, b, fb, sb, adv, ctx, Ys, bits, ubits, &nadv);
        }
        mpfr_sub(fa, b, tbest, MPFR_RNDN);
        if (best < 0 || mpfr_sgn(fa) * dir < 0) { best = (int)i; mpfr_set(tbest, b, MPFR_RNDN); }
    }
    if (best >= 0) {
        mpfr_set(te, tbest, MPFR_RNDN);
        if (mpfr_equal_p(tbest, tb) && nadv == 0) for (size_t i = 0; i < d; i++) mpfr_set(Ye[i], Yb[i], MPFR_RNDN);
        else if (!adv(ctx, tbest, Ye)) best = -1;
        if (best >= 0) {
            NdEvent* v = &E->ev[best];
            int tr;
            if (v->kind == EVK_ZERO && !ev_value_mp(P, v, tbest, Ye, bits, v->groot_mp, &tr))
                mpfr_set_zero(v->groot_mp, 1);
        }
    }
    if (best < 0) {
        for (size_t i = 0; i < E->n; i++) {
            NdEvent* v = &E->ev[i];
            if (!v->active) continue;
            v->known = okb[i];
            if (okb[i]) { mpfr_set(v->gp_mp, gb[i], MPFR_RNDN); v->tp = tbv[i]; }
        }
    }
    for (size_t i = 0; i < E->n; i++) mpfr_clear(gb[i]);
    for (size_t i = 0; i < d; i++) mpfr_clear(Ys[i]);
    free(gb); free(tbv); free(okb); free(Ys);
    mpfr_clears(a, fa, b, fb, tbest, (mpfr_ptr)0);
    return best;
}

unsigned nd_events_fire_mp(NdProblem* P, int k, const mpfr_t te, mpfr_t* Ye,
                           mpfr_t* Ynew, long bits) {
    struct NdEvents* E = P->events;
    size_t d = P->d;
    for (size_t i = 0; i < d; i++) mpfr_set(Ynew[i], Ye[i], MPFR_RNDN);
    if (!E || k < 0 || (size_t)k >= E->n) return 0;
    NdEvent* v = &E->ev[k];
    unsigned flags = 0;
    ev_bind_mp(P, te, Ye);
    ev_run_actions(v);
    if (v->nreset) {
        mpfr_t* val = malloc(sizeof(mpfr_t) * v->nreset);
        for (size_t j = 0; j < v->nreset; j++) mpfr_init2(val[j], bits);
        bool ok = true;
        ev_bind_mp(P, te, Ye);
        for (size_t j = 0; j < v->nreset && ok; j++) ok = ev_eval_mp(v->rrhs[j], val[j], bits);
        if (ok) { for (size_t j = 0; j < v->nreset; j++) mpfr_set(Ynew[v->ridx[j]], val[j], MPFR_RNDN); flags |= ND_EVF_RESET; }
        else ev_warn("evrs", "a WhenEvent state reset did not evaluate to a real number; the state is unchanged");
        for (size_t j = 0; j < v->nreset; j++) mpfr_clear(val[j]);
        free(val);
    }
    if (v->stop) flags |= ND_EVF_STOP;
    if (v->remove) v->active = false;
    mpfr_t g, lim; mpfr_inits2(bits, g, lim, (mpfr_ptr)0);
    for (size_t i = 0; i < E->n; i++) {
        NdEvent* w = &E->ev[i];
        if (!w->active) continue;
        ev_mp_ensure(w, bits);
        int tr;
        w->known = ev_value_mp(P, w, te, Ynew, bits, g, &tr);
        if (!w->known) continue;
        mpfr_set(w->gp_mp, g, MPFR_RNDN); w->tp = tr;
        if ((int)i == k && w->kind == EVK_ZERO) {
            mpfr_sub(lim, g, w->groot_mp, MPFR_RNDN); mpfr_abs(lim, lim, MPFR_RNDN);
            mpfr_t sc; mpfr_init2(sc, bits);
            mpfr_mul_d(sc, w->gscale_mp, 1e-6, MPFR_RNDN);
            if (mpfr_cmp(lim, sc) <= 0) mpfr_set_zero(w->gp_mp, 1);
            mpfr_clear(sc);
        }
    }
    mpfr_clears(g, lim, (mpfr_ptr)0);
    return flags;
}

#endif /* USE_MPFR */
