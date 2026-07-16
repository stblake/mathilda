/* integrate_beta.c
 *
 * Definite integration by Euler-Beta reduction.  See integrate_beta.h.
 *
 *   [0,1]:      x^(k-1) (1-x)^(l-1) -> Beta[k,l]; Log[x]^i Log[1-x]^j weights
 *               -> the mixed parameter derivative of Beta (Beta already carries
 *               a PolyGamma derivative rule, so plain D suffices).
 *   [0,Pi/2]:   Sin[x]^m Cos[x]^n -> Beta[(m+1)/2,(n+1)/2]/2.
 *   [0,Pi]:     n odd -> 0; n even -> 2 x the quarter-period value.
 *   [0,2 Pi]:   m or n odd -> 0; both even -> 4 x the quarter-period value.
 *
 * Every value is gated on the Beta convergence strip (Re arg > 0), proved by
 * Simplify under the caller's Assumptions; an undecided strip becomes a
 * ConditionalExpression, a refuted one a NULL (fall through).
 */

#include "integrate_beta.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"

#include <string.h>
#include <stdbool.h>

/* ---- construction / evaluation helpers ---------------------------------- */

static Expr* cp(const Expr* e) { return expr_copy((Expr*)e); }
static Expr* mk_int(long v) { return expr_new_integer((int64_t)v); }

static Expr* mk_fn(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
static Expr* mk_fn1(const char* h, Expr* a) { Expr* v[1]={a}; return mk_fn(h,v,1); }
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) { Expr* v[2]={a,b}; return mk_fn(h,v,2); }

static Expr* Times_(Expr* a, Expr* b) { return mk_fn2("Times", a, b); }
static Expr* Plus_(Expr* a, Expr* b)  { return mk_fn2("Plus", a, b); }
static Expr* Beta_(Expr* a, Expr* b)  { return mk_fn2("Beta", a, b); }

static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);   /* evaluate does not free its input */
    expr_free(call);
    return r;
}
static Expr* ev1(const char* name, Expr* a) { return eval_take(mk_fn1(name, a)); }
static Expr* simp(Expr* e) { return ev1("Simplify", e); }
static Expr* simp2(Expr* e, Expr* as) {
    if (!as) return simp(e);
    return eval_take(mk_fn2("Simplify", e, cp(as)));
}

static bool head_name_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}
static bool sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}
static bool is_symbol(const Expr* e, const Expr* x) {
    return e->type == EXPR_SYMBOL && x->type == EXPR_SYMBOL &&
           e->data.symbol.name == x->data.symbol.name;
}

static bool contains_symbol(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return is_symbol(e, x);
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol(e->data.function.args[i], x)) return true;
    return false;
}
static bool free_of_x(const Expr* e, const Expr* x) { return !contains_symbol(e, x); }

/* True iff `e` mentions any divergence / indeterminacy symbol or an unevaluated
 * Integrate -- i.e. it is NOT a usable finite closed form. */
static bool mentions_nonfinite(const Expr* e) {
    if (!e) return true;
    if (e->type == EXPR_SYMBOL)
        return sym_is(e, "Infinity") || sym_is(e, "ComplexInfinity") ||
               sym_is(e, "Indeterminate") || sym_is(e, "Undefined");
    if (head_name_is(e, "DirectedInfinity") || head_name_is(e, "Integrate"))
        return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (mentions_nonfinite(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (mentions_nonfinite(e->data.function.args[i])) return true;
    return false;
}

/* Prove a predicate is True under `as` via Simplify.  Returns +1 (True),
 * -1 (False), 0 (undecided). */
static int prove(Expr* pred, Expr* as) {
    Expr* s = simp2(pred, as);
    int r = 0;
    if (sym_is(s, "True")) r = 1;
    else if (sym_is(s, "False")) r = -1;
    expr_free(s);
    return r;
}

/* Structural test for the factor (1 - x). */
static bool is_one_minus_x(const Expr* g, const Expr* x) {
    if (!head_name_is(g, "Plus") || g->data.function.arg_count != 2) return false;
    Expr* p = g->data.function.args[0];
    Expr* q = g->data.function.args[1];
    /* one summand is the integer 1, the other is -x = Times[-1, x]. */
    Expr* one = NULL; Expr* negx = NULL;
    if (p->type == EXPR_INTEGER && p->data.integer == 1) { one = p; negx = q; }
    else if (q->type == EXPR_INTEGER && q->data.integer == 1) { one = q; negx = p; }
    if (!one) return false;
    return head_name_is(negx, "Times") && negx->data.function.arg_count == 2 &&
           negx->data.function.args[0]->type == EXPR_INTEGER &&
           negx->data.function.args[0]->data.integer == -1 &&
           is_symbol(negx->data.function.args[1], x);
}

/* Iterate the top-level product factors of f, calling visit() on each.  Returns
 * false as soon as visit() returns false. */
typedef bool (*factor_visitor)(Expr* factor, void* ctx);
static bool for_each_factor(Expr* f, factor_visitor visit, void* ctx) {
    if (head_name_is(f, "Times")) {
        for (size_t i = 0; i < f->data.function.arg_count; i++)
            if (!visit(f->data.function.args[i], ctx)) return false;
        return true;
    }
    return visit(f, ctx);
}

/* Apply the Beta convergence-strip gate (defined below). */
static Expr* beta_gate(Expr* value, Expr* A, Expr* B, Expr* assumptions);

/* ---- Beta on [0,1] ------------------------------------------------------ */

typedef struct {
    const Expr* x;
    Expr* aexp;   /* accumulated exponent of x        (owned) */
    Expr* bexp;   /* accumulated exponent of (1-x)     (owned) */
    long  logi;   /* power of Log[x]                            */
    long  logj;   /* power of Log[1-x]                          */
    Expr* C;      /* accumulated x-free constant       (owned) */
    bool  ok;
} BetaCtx;

/* Is factor Log[x] (returns 1) / Log[1-x] (returns 2) / neither (0)? */
static int log_kind(const Expr* g, const Expr* x) {
    if (!head_name_is(g, "Log") || g->data.function.arg_count != 1) return 0;
    Expr* a = g->data.function.args[0];
    if (is_symbol(a, x)) return 1;
    if (is_one_minus_x(a, x)) return 2;
    return 0;
}

static bool beta_visit(Expr* g, void* vctx) {
    BetaCtx* c = (BetaCtx*)vctx;
    const Expr* x = c->x;
    if (free_of_x(g, x)) { c->C = Times_(c->C, cp(g)); return true; }
    if (is_symbol(g, x)) { c->aexp = Plus_(c->aexp, mk_int(1)); return true; }
    if (is_one_minus_x(g, x)) { c->bexp = Plus_(c->bexp, mk_int(1)); return true; }
    { int lk = log_kind(g, x);
      if (lk == 1) { c->logi += 1; return true; }
      if (lk == 2) { c->logj += 1; return true; } }
    if (head_name_is(g, "Power") && g->data.function.arg_count == 2) {
        Expr* base = g->data.function.args[0];
        Expr* e    = g->data.function.args[1];
        if (is_symbol(base, x) && free_of_x(e, x)) {
            c->aexp = Plus_(c->aexp, cp(e)); return true;
        }
        if (is_one_minus_x(base, x) && free_of_x(e, x)) {
            c->bexp = Plus_(c->bexp, cp(e)); return true;
        }
        int lk = log_kind(base, x);
        if ((lk == 1 || lk == 2) && e->type == EXPR_INTEGER && e->data.integer > 0) {
            if (lk == 1) c->logi += e->data.integer; else c->logj += e->data.integer;
            return true;
        }
    }
    c->ok = false;               /* unrecognised factor: not a Beta integrand */
    return false;
}

/* d^i/dP^i d^j/dQ^j Beta[P,Q] with P,Q placeholder symbols, then P->A, Q->B. */
static Expr* beta_param_deriv(Expr* A, Expr* B, long i, long j) {
    const char* P = "Integrate`Beta`P";
    const char* Q = "Integrate`Beta`Q";
    Expr* d = Beta_(expr_new_symbol(P), expr_new_symbol(Q));
    for (long k = 0; k < i; k++) d = eval_take(mk_fn2("D", d, expr_new_symbol(P)));
    for (long k = 0; k < j; k++) d = eval_take(mk_fn2("D", d, expr_new_symbol(Q)));
    Expr* rules = mk_fn2("List", mk_fn2("Rule", expr_new_symbol(P), A),
                                  mk_fn2("Rule", expr_new_symbol(Q), B));
    return eval_take(mk_fn2("ReplaceAll", d, rules));
}

Expr* integrate_beta_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!contains_symbol(f, x)) return NULL;
    /* Interval must be exactly [0,1]. */
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    if (!(b->type == EXPR_INTEGER && b->data.integer == 1)) return NULL;

    BetaCtx c = { x, mk_int(0), mk_int(0), 0, 0, mk_int(1), true };
    for_each_factor(f, beta_visit, &c);
    if (!c.ok) {
        expr_free(c.aexp); expr_free(c.bexp); expr_free(c.C);
        return NULL;
    }
    /* Beta arguments k = aexp + 1, l = bexp + 1. */
    Expr* A = simp(Plus_(c.aexp, mk_int(1)));
    Expr* B = simp(Plus_(c.bexp, mk_int(1)));

    /* Value: constant * (Beta or its mixed parameter derivative). */
    Expr* core = (c.logi == 0 && c.logj == 0)
        ? Beta_(cp(A), cp(B))
        : beta_param_deriv(cp(A), cp(B), c.logi, c.logj);
    Expr* value = simp2(Times_(c.C, core), assumptions);
    if (!value || mentions_nonfinite(value)) {
        if (value) expr_free(value);
        expr_free(A); expr_free(B);
        return NULL;
    }

    /* Convergence strip: Re[k] > 0 && Re[l] > 0.  Prove on the plain positivity
     * k > 0 && l > 0 -- a sound *sufficient* condition (real positivity implies
     * Re > 0, and a complex argument leaves Greater undecided, never wrongly
     * True), which Simplify can actually discharge under a `k > 0` assumption.
     * Display the WL-faithful Re form when the strip stays undecided. */
    Expr* value2 = beta_gate(value, A, B, assumptions);
    expr_free(A); expr_free(B);
    return value2;
}

/* Apply the Beta convergence-strip gate to an already-evaluated `value` whose
 * Beta arguments are A, B (all borrowed except `value`, which is consumed).
 * Returns `value` when the strip is proved, NULL (freeing `value`) when it is
 * refuted, or a ConditionalExpression when it is undecided. */
static Expr* beta_gate(Expr* value, Expr* A, Expr* B, Expr* assumptions) {
    Expr* cond_prove = mk_fn2("And",
        mk_fn2("Greater", cp(A), mk_int(0)),
        mk_fn2("Greater", cp(B), mk_int(0)));
    int pv = prove(cond_prove, assumptions);
    if (pv == 1)  return value;
    if (pv == -1) { expr_free(value); return NULL; }
    Expr* cond_show = mk_fn2("And",
        mk_fn2("Greater", mk_fn1("Re", cp(A)), mk_int(0)),
        mk_fn2("Greater", mk_fn1("Re", cp(B)), mk_int(0)));
    Expr* cs = simp2(cond_show, assumptions);
    return eval_take(mk_fn2("ConditionalExpression", value, cs));
}

/* ---- Sin^m Cos^n over [0, c] -------------------------------------------- */

typedef struct {
    const Expr* x;
    Expr* mexp;   /* exponent of Sin[x] (owned) */
    Expr* nexp;   /* exponent of Cos[x] (owned) */
    Expr* C;      /* x-free constant     (owned) */
    bool  ok;
} TrigCtx;

/* Sin[x] -> 1, Cos[x] -> 2, neither -> 0. */
static int trig_kind(const Expr* g, const Expr* x) {
    if (head_name_is(g, "Sin") && g->data.function.arg_count == 1 &&
        is_symbol(g->data.function.args[0], x)) return 1;
    if (head_name_is(g, "Cos") && g->data.function.arg_count == 1 &&
        is_symbol(g->data.function.args[0], x)) return 2;
    return 0;
}

static bool trig_visit(Expr* g, void* vctx) {
    TrigCtx* c = (TrigCtx*)vctx;
    const Expr* x = c->x;
    if (free_of_x(g, x)) { c->C = Times_(c->C, cp(g)); return true; }
    int tk = trig_kind(g, x);
    if (tk == 1) { c->mexp = Plus_(c->mexp, mk_int(1)); return true; }
    if (tk == 2) { c->nexp = Plus_(c->nexp, mk_int(1)); return true; }
    if (head_name_is(g, "Power") && g->data.function.arg_count == 2) {
        int bk = trig_kind(g->data.function.args[0], x);
        Expr* e = g->data.function.args[1];
        if (bk == 1 && free_of_x(e, x)) { c->mexp = Plus_(c->mexp, cp(e)); return true; }
        if (bk == 2 && free_of_x(e, x)) { c->nexp = Plus_(c->nexp, cp(e)); return true; }
    }
    c->ok = false;
    return false;
}

/* Classify the upper limit as one of the canonical trig intervals.
 * Returns 1 for Pi/2, 2 for Pi, 4 for 2 Pi, 0 otherwise. */
static int trig_interval(Expr* b) {
    Expr* d = simp(mk_fn2("Times", mk_fn2("Power", expr_new_symbol("Pi"), mk_int(-1)), cp(b)));
    int r = 0;
    if (d->type == EXPR_INTEGER && d->data.integer == 2) r = 4;    /* 2 Pi */
    else if (d->type == EXPR_INTEGER && d->data.integer == 1) r = 2; /* Pi */
    else if (head_name_is(d, "Rational") && d->data.function.arg_count == 2 &&
             d->data.function.args[0]->type == EXPR_INTEGER &&
             d->data.function.args[0]->data.integer == 1 &&
             d->data.function.args[1]->type == EXPR_INTEGER &&
             d->data.function.args[1]->data.integer == 2) r = 1;   /* Pi/2 */
    expr_free(d);
    return r;
}

/* Non-negative even integer? */
static bool is_even_nonneg_int(const Expr* e) {
    return e->type == EXPR_INTEGER && e->data.integer >= 0 && (e->data.integer % 2) == 0;
}
static bool is_nonneg_int(const Expr* e) {
    return e->type == EXPR_INTEGER && e->data.integer >= 0;
}

Expr* integrate_trigpower_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!contains_symbol(f, x)) return NULL;
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    int kind = trig_interval(b);
    if (kind == 0) return NULL;

    TrigCtx c = { x, mk_int(0), mk_int(0), mk_int(1), true };
    for_each_factor(f, trig_visit, &c);
    if (!c.ok) { expr_free(c.mexp); expr_free(c.nexp); expr_free(c.C); return NULL; }

    Expr* m = simp(c.mexp);
    Expr* n = simp(c.nexp);

    /* Quarter-period value: (1/2) Beta[(m+1)/2, (n+1)/2]. */
    Expr* Barg1 = mk_fn2("Times", mk_fn2("Rational", mk_int(1), mk_int(2)),
                                   Plus_(cp(m), mk_int(1)));
    Expr* Barg2 = mk_fn2("Times", mk_fn2("Rational", mk_int(1), mk_int(2)),
                                   Plus_(cp(n), mk_int(1)));
    Expr* Barg1s = simp(Barg1), *Barg2s = simp(Barg2);
    Expr* quarter = Times_(mk_fn2("Rational", mk_int(1), mk_int(2)),
                           Beta_(cp(Barg1s), cp(Barg2s)));

    /* Interval multiplier / parity vanishing. */
    long mult = 1;
    bool vanishes = false, need_int = (kind != 1);
    if (kind == 2) {                       /* [0, Pi]: n odd -> 0, n even -> 2x */
        if (!is_nonneg_int(n)) { /* symbolic n: cannot decide parity */
            expr_free(m); expr_free(n); expr_free(Barg1s); expr_free(Barg2s);
            expr_free(quarter); expr_free(c.C); return NULL;
        }
        if (is_even_nonneg_int(n)) mult = 2; else vanishes = true;
    } else if (kind == 4) {                /* [0, 2Pi]: m|n odd -> 0, both even -> 4x */
        if (!is_nonneg_int(m) || !is_nonneg_int(n)) {
            expr_free(m); expr_free(n); expr_free(Barg1s); expr_free(Barg2s);
            expr_free(quarter); expr_free(c.C); return NULL;
        }
        if (is_even_nonneg_int(m) && is_even_nonneg_int(n)) mult = 4; else vanishes = true;
    }
    (void)need_int;

    Expr* value;
    if (vanishes) {
        expr_free(quarter);
        value = simp2(Times_(c.C, mk_int(0)), assumptions);   /* C * 0 = 0 */
    } else {
        value = simp2(Times_(c.C, Times_(mk_int(mult), quarter)), assumptions);
    }
    expr_free(m); expr_free(n);

    if (!value || mentions_nonfinite(value)) {
        if (value) expr_free(value);
        expr_free(Barg1s); expr_free(Barg2s);
        return NULL;
    }
    if (vanishes) { expr_free(Barg1s); expr_free(Barg2s); return value; }

    /* Convergence strip: Re[(m+1)/2] > 0 && Re[(n+1)/2] > 0. */
    Expr* value2 = beta_gate(value, Barg1s, Barg2s, assumptions);
    expr_free(Barg1s); expr_free(Barg2s);
    return value2;
}

/* ---- builtins ----------------------------------------------------------- */

/* Shared parse of `head[f, {x,a,b}, opts]`, extracting the spec + Assumptions. */
static bool parse_definite(Expr* res, Expr** f, Expr** x, Expr** a, Expr** b,
                           Expr** assumptions) {
    if (res->type != EXPR_FUNCTION) return false;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return false;
    *f = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!head_name_is(spec, "List") || spec->data.function.arg_count != 3) return false;
    *x = spec->data.function.args[0];
    *a = spec->data.function.args[1];
    *b = spec->data.function.args[2];
    if ((*x)->type != EXPR_SYMBOL) return false;
    *assumptions = NULL;
    for (size_t t = 2; t < argc; t++) {
        Expr* opt = res->data.function.args[t];
        if (opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2 &&
            opt->data.function.head->type == EXPR_SYMBOL &&
            (opt->data.function.head->data.symbol.name == SYM_Rule ||
             opt->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
            Expr* lhs = opt->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL && strcmp(lhs->data.symbol.name, "Assumptions") == 0) {
                *assumptions = opt->data.function.args[1];
                continue;
            }
        }
        return false;
    }
    return true;
}

Expr* builtin_integrate_beta(Expr* res) {
    Expr *f, *x, *a, *b, *as;
    if (!parse_definite(res, &f, &x, &a, &b, &as)) return NULL;
    return integrate_beta_try(f, x, a, b, as);
}
Expr* builtin_integrate_trigpower(Expr* res) {
    Expr *f, *x, *a, *b, *as;
    if (!parse_definite(res, &f, &x, &a, &b, &as)) return NULL;
    return integrate_trigpower_try(f, x, a, b, as);
}

void integrate_beta_init(void) {
    symtab_add_builtin("Integrate`Beta", builtin_integrate_beta);
    symtab_get_def("Integrate`Beta")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`Beta",
        "Integrate`Beta[f, {x, 0, 1}] evaluates a definite integral of the form "
        "x^(k-1) (1-x)^(l-1) over [0,1] as the Euler Beta function Beta[k, l], "
        "and an integrand additionally weighted by Log[x]^i Log[1-x]^j as the "
        "mixed parameter derivative of Beta.  Gated on the convergence strip "
        "Re[k] > 0 && Re[l] > 0 (an undecided strip yields a "
        "ConditionalExpression).  Returns unevaluated when the integrand is not "
        "of this form or the interval is not [0,1].");

    symtab_add_builtin("Integrate`TrigPower", builtin_integrate_trigpower);
    symtab_get_def("Integrate`TrigPower")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`TrigPower",
        "Integrate`TrigPower[f, {x, 0, c}] evaluates Sin[x]^m Cos[x]^n over a "
        "canonical interval: over [0, Pi/2] it is Beta[(m+1)/2, (n+1)/2]/2; over "
        "[0, Pi] it is 0 when n is an odd integer and twice the quarter-period "
        "value when n is even; over [0, 2 Pi] it is 0 unless both m and n are "
        "even integers, and four times the quarter-period value otherwise.  "
        "Gated on Re[(m+1)/2] > 0 && Re[(n+1)/2] > 0.  Returns unevaluated when "
        "the integrand or interval is not of this form.");
}
