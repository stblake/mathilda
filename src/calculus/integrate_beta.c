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
#include "arithmetic.h"   /* arith_warnings_mute_push/pop */
#include "parse.h"

#include <string.h>
#include <stdbool.h>
#include <stdarg.h>

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
/* True iff some subexpression of e has head `name`. */
static bool contains_head(const Expr* e, const char* name) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is(e, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_head(e->data.function.args[i], name)) return true;
    return false;
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

/* Prove a predicate True under `as` via Refine (the assumption-aware sign
 * prover): Simplify leaves an implication like `c > 0 => c+1 > 0` unproved
 * (it returns `c > -1`), but Refine discharges it.  Consumes pred. */
static bool prove_ref(Expr* pred, Expr* as) {
    Expr* r = as ? eval_take(mk_fn2("Refine", pred, cp(as))) : ev1("Refine", pred);
    bool t = sym_is(r, "True");
    expr_free(r);
    return t;
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

/* d^i/dP^i d^j/dQ^j Beta[P,Q] with P,Q placeholder symbols, then P->A, Q->B.
 * Consumes A and B.
 *
 * At a boundary point of Beta's own strip the substituted derivative is a
 * 0 * ComplexInfinity (Beta[0, 1] PolyGamma[...]) although the integral
 * converges: a Log[1-x]^j weight vanishes like x^j at 0 (a Log[x]^i weight like
 * (1-x)^i at 1), which widens the convergence strip to A > -j, B > -i.  On that
 * wider strip the integral is analytic in A and B, so its value is the
 * (removable) limit of the meromorphic derivative.  The substitution is
 * therefore an internal probe -- muted, so its Infinity::indet neither prints
 * nor trips Check[] -- and a non-finite value falls back to the iterated limit.
 * The caller still gates the result on the widened strip. */
static Expr* beta_param_deriv(Expr* A, Expr* B, long i, long j) {
    const char* P = "Integrate`Beta`P";
    const char* Q = "Integrate`Beta`Q";
    Expr* d = Beta_(expr_new_symbol(P), expr_new_symbol(Q));
    for (long k = 0; k < i; k++) d = eval_take(mk_fn2("D", d, expr_new_symbol(P)));
    for (long k = 0; k < j; k++) d = eval_take(mk_fn2("D", d, expr_new_symbol(Q)));
    arith_warnings_mute_push();
    Expr* rules = mk_fn2("List", mk_fn2("Rule", expr_new_symbol(P), cp(A)),
                                  mk_fn2("Rule", expr_new_symbol(Q), cp(B)));
    Expr* v = eval_take(mk_fn2("ReplaceAll", cp(d), rules));
    if ((!v || mentions_nonfinite(v)) && (i > 0 || j > 0)) {
        /* Iterated limit, the regular variable substituted first. */
        for (int order = 0; order < 2; order++) {
            const char* sv = order ? P : Q;      /* substituted */
            const char* lv = order ? Q : P;      /* limit variable */
            Expr* sval = order ? A : B;
            Expr* lval = order ? B : A;
            Expr* part = eval_take(mk_fn2("ReplaceAll", cp(d),
                             mk_fn2("Rule", expr_new_symbol(sv), cp(sval))));
            if (!part || mentions_nonfinite(part)) { if (part) expr_free(part); continue; }
            Expr* lim = eval_take(mk_fn2("Limit", part,
                            mk_fn2("Rule", expr_new_symbol(lv), cp(lval))));
            if (lim && !mentions_nonfinite(lim) && !contains_head(lim, "Limit")) {
                if (v) expr_free(v);
                v = lim;
                break;
            }
            if (lim) expr_free(lim);
        }
    }
    arith_warnings_mute_pop();
    expr_free(d); expr_free(A); expr_free(B);
    return v;
}

/* PolyGamma[n, k] at a positive integer k (n >= 1) in closed form:
 *   PolyGamma[n, k] = (-1)^(n+1) n! (Zeta[n+1] - Sum_{m=1}^{k-1} m^-(n+1)),
 * which the Beta parameter derivatives produce at integer arguments (the
 * widened-strip limit of Log[x] Log[1-x]/x gives -PolyGamma[2,1]/2 = Zeta[3]).
 * Bounded (k, n <= 64), exact. */
static Expr* beta_polygamma_to_zeta(Expr* v) {
    if (!v || !contains_head(v, "PolyGamma")) return v;
    /* Half-integer (and other rational) arguments go through FunctionExpand:
     * PolyGamma[0, 3/2] stays symbolic since v0.6130, as in the reference,
     * so the Beta value would otherwise keep it. */
    static const char* RULE =
        "{PolyGamma[n_Integer, k_Integer] /; 1 <= n <= 64 && 1 <= k <= 64 :> "
        "(-1)^(n + 1) n! (Zeta[n + 1] - Sum[1/m^(n + 1), {m, 1, k - 1}]), "
        "PolyGamma[n_Integer, q_Rational] /; 0 <= n <= 64 && q > 0 :> "
        "FunctionExpand[PolyGamma[n, q]]}";
    Expr* rule = parse_expression(RULE);
    if (!rule) return v;
    return eval_take(mk_fn2("ReplaceAll", v, rule));
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
        : beta_polygamma_to_zeta(beta_param_deriv(cp(A), cp(B), c.logi, c.logj));
    Expr* value = simp2(Times_(c.C, core), assumptions);
    if (!value || mentions_nonfinite(value)) {
        if (value) expr_free(value);
        expr_free(A); expr_free(B);
        return NULL;
    }

    /* Convergence strip: Re[k] > -logj && Re[l] > -logi (a Log[1-x]^j weight
     * vanishes like x^j at 0, a Log[x]^i weight like (1-x)^i at 1).  Prove on
     * the plain real inequalities -- a sound *sufficient* condition (real
     * positivity implies Re > 0, and a complex argument leaves Greater
     * undecided, never wrongly True), which Simplify can actually discharge
     * under a `k > 0` assumption.  Display the WL-faithful Re form when the
     * strip stays undecided. */
    Expr* Ash = simp(Plus_(cp(A), mk_int(c.logj)));
    Expr* Bsh = simp(Plus_(cp(B), mk_int(c.logi)));
    Expr* value2 = beta_gate(value, Ash, Bsh, assumptions);
    expr_free(Ash); expr_free(Bsh);
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
    /* Simplify leaves some plain implications unreduced (1/2 + n > 0 under
     * n > 0 came back as itself, so the value was needlessly conditioned);
     * the Refine sign prover discharges them. */
    Expr* cond_ref = assumptions ? cp(cond_prove) : NULL;
    int pv = prove(cond_prove, assumptions);
    if (pv == 0 && cond_ref) { if (prove_ref(cond_ref, assumptions)) pv = 1; }
    else if (cond_ref) expr_free(cond_ref);
    if (pv == 1)  return value;
    if (pv == -1) { expr_free(value); return NULL; }
    Expr* cond_show = mk_fn2("And",
        mk_fn2("Greater", mk_fn1("Re", cp(A)), mk_int(0)),
        mk_fn2("Greater", mk_fn1("Re", cp(B)), mk_int(0)));
    Expr* cs = simp2(cond_show, assumptions);
    return eval_take(mk_fn2("ConditionalExpression", value, cs));
}

/* Structural test for the factor (1 - x^2). */
static bool is_one_minus_xsq(const Expr* g, const Expr* x) {
    if (!head_name_is(g, "Plus") || g->data.function.arg_count != 2) return false;
    Expr* p = g->data.function.args[0];
    Expr* q = g->data.function.args[1];
    Expr* one = NULL; Expr* negxsq = NULL;
    if (p->type == EXPR_INTEGER && p->data.integer == 1) { one = p; negxsq = q; }
    else if (q->type == EXPR_INTEGER && q->data.integer == 1) { one = q; negxsq = p; }
    if (!one) return false;
    if (!head_name_is(negxsq, "Times") || negxsq->data.function.arg_count != 2) return false;
    if (!(negxsq->data.function.args[0]->type == EXPR_INTEGER &&
          negxsq->data.function.args[0]->data.integer == -1)) return false;
    Expr* pw = negxsq->data.function.args[1];
    return head_name_is(pw, "Power") && pw->data.function.arg_count == 2 &&
           is_symbol(pw->data.function.args[0], x) &&
           pw->data.function.args[1]->type == EXPR_INTEGER &&
           pw->data.function.args[1]->data.integer == 2;
}

/* Half-Beta on [0,1]:  Integrate[C x^(c-1) (1-x^2)^e Log[x]^i, {x,0,1}]
 *   = C (1/2)^(i+1) d^i/dA^i Beta[A, e+1]|_{A = c/2},
 * the substitution t = x^2 (dt = 2x dx, Log[x] = (1/2)Log[t]).  Handles the
 * quadratic radical/power (1-x^2)^e the ordinary Beta stage (which models only
 * the linear 1-x) cannot.  Requires the (1-x^2) factor to be present (else it is
 * a plain Beta, left to integrate_beta_try).  No Log[1-x^2] weight (that would
 * add a B-parameter derivative).  Gated on the Beta strip A>0 && e+1>0. */
Expr* integrate_halfbeta_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    if (!(b->type == EXPR_INTEGER && b->data.integer == 1)) return NULL;

    Expr* aexp = mk_int(0); Expr* e = mk_int(0); Expr* C = mk_int(1);
    long logi = 0; bool has_xsq = false; bool ok = true;
    Expr** fac; size_t nf; Expr* one[1];
    if (head_name_is(f, "Times")) { fac = f->data.function.args; nf = f->data.function.arg_count; }
    else { one[0] = f; fac = one; nf = 1; }
    for (size_t i = 0; i < nf && ok; i++) {
        Expr* g = fac[i];
        if (free_of_x(g, x)) { C = Times_(C, cp(g)); continue; }
        if (is_symbol(g, x)) { aexp = Plus_(aexp, mk_int(1)); continue; }
        if (is_one_minus_xsq(g, x)) { e = Plus_(e, mk_int(1)); has_xsq = true; continue; }
        if (log_kind(g, x) == 1) { logi += 1; continue; }
        if (head_name_is(g, "Power") && g->data.function.arg_count == 2) {
            Expr* base = g->data.function.args[0];
            Expr* ex   = g->data.function.args[1];
            if (is_symbol(base, x) && free_of_x(ex, x)) { aexp = Plus_(aexp, cp(ex)); continue; }
            if (is_one_minus_xsq(base, x) && free_of_x(ex, x)) { e = Plus_(e, cp(ex)); has_xsq = true; continue; }
            if (log_kind(base, x) == 1 && ex->type == EXPR_INTEGER && ex->data.integer > 0) {
                logi += (long) ex->data.integer; continue;
            }
        }
        ok = false;
    }
    if (!ok || !has_xsq) { expr_free(aexp); expr_free(e); expr_free(C); return NULL; }

    Expr* A = simp(Times_(Plus_(aexp, mk_int(1)), mk_fn2("Rational", mk_int(1), mk_int(2))));
    Expr* B = simp(Plus_(e, mk_int(1)));
    Expr* core = (logi == 0) ? Beta_(cp(A), cp(B)) : beta_param_deriv(cp(A), cp(B), logi, 0);
    Expr* half = mk_fn2("Power", mk_int(2), mk_int(-(logi + 1)));
    Expr* value = simp2(Times_(C, Times_(half, core)), assumptions);
    if (!value || mentions_nonfinite(value)) {
        if (value) expr_free(value);
        expr_free(A); expr_free(B);
        return NULL;
    }
    Expr* value2 = beta_gate(value, A, B, assumptions);
    expr_free(A); expr_free(B);
    return value2;
}

/* ---- Reductions INTO the Beta family ------------------------------------ */
/* Two shape-gated pre-FTC stages that integrate nothing themselves: each
 * rewrites the integrand into a shape the [0,1] Beta stages own and hands it to
 * them.  They run ahead of Newton-Leibniz for the same reason those stages do:
 * with a SYMBOLIC exponent the antiderivative is a non-elementary incomplete
 * Beta, and the indefinite cascade only finds that out after its full search
 * -- 55-66 s on the two integrals below, 45 s of it ParallelMixedTower running
 * into its own backstop -- before the definite cascade reaches a stage that
 * closes them.  Integer and numeric exponents are left to Newton-Leibniz,
 * whose elementary antiderivative is fast there. */

/* Evaluate the pure-function template `t` (#1..#n) on borrowed arguments.
 * Reparsed per call: each stage calls it a few times, and only once its
 * structural gate has passed.  Template pattern variables live in
 * Integrate`BetaR`, never in the user's context. */
static Expr* bt_x(const char* t, int n, ...) {
    Expr* body = parse_expression(t);
    if (!body) return NULL;
    Expr* args[6];
    va_list ap;
    va_start(ap, n);
    for (int i = 0; i < n && i < 6; i++) args[i] = cp(va_arg(ap, Expr*));
    va_end(ap);
    Expr* fn = mk_fn1("Function", body);
    return eval_take(expr_new_function(fn, n ? args : NULL, (size_t)n));
}
static bool bt_true(Expr* e) {            /* consumes e */
    bool r = e && sym_is(e, "True");
    if (e) expr_free(e);
    return r;
}

/* True iff some Power[base, e] in f has a base depending on x and an exponent
 * free of x that is not a number (a, 2 n - 1, -2 a - 2, ...). */
static bool has_symbolic_exponent(const Expr* f, const Expr* x) {
    if (!f || f->type != EXPR_FUNCTION) return false;
    if (head_name_is(f, "Power") && f->data.function.arg_count == 2 &&
        contains_symbol(f->data.function.args[0], x) &&
        free_of_x(f->data.function.args[1], x) &&
        !bt_true(ev1("NumericQ", cp(f->data.function.args[1]))))
        return true;
    for (size_t i = 0; i < f->data.function.arg_count; i++)
        if (has_symbolic_exponent(f->data.function.args[i], x)) return true;
    return false;
}

/* The Beta-family closure of g on [0,1] in t: Beta, Half-Beta, Euler 2F1. */
static Expr* beta_family_01(Expr* g, Expr* t, Expr* assumptions) {
    Expr* zero = mk_int(0); Expr* one = mk_int(1);
    Expr* r = integrate_beta_try(g, t, zero, one, assumptions);
    if (!r) r = integrate_halfbeta_try(g, t, zero, one, assumptions);
    if (!r) r = integrate_euler_2f1_try(g, t, zero, one, assumptions);
    expr_free(zero); expr_free(one);
    return r;
}

/* Beta[p, q] -> Gamma[p] Gamma[q]/Gamma[p + q], then Simplify -- inside the
 * ConditionalExpression when there is one.  The user-facing form of a value
 * the Beta stages build as Beta (FunctionExpand has no Beta rule).  Consumes v. */
static Expr* beta_to_gamma(Expr* v, Expr* assumptions) {
    if (!v) return NULL;
    Expr* as = assumptions ? cp(assumptions) : expr_new_symbol("True");
    Expr* r = bt_x(
        "With[{g = Function[Integrate`BetaR`v, Simplify[Integrate`BetaR`v /. "
        "Beta[Integrate`BetaR`p_, Integrate`BetaR`q_] :> Gamma[Integrate`BetaR`p] "
        "Gamma[Integrate`BetaR`q]/Gamma[Integrate`BetaR`p + Integrate`BetaR`q], #2]]}, "
        "If[Head[#1] === ConditionalExpression, "
        "ConditionalExpression[g[#1[[1]]], #1[[2]]], g[#1]]]", 2, v, as);
    expr_free(v); expr_free(as);
    return r;
}

/* t = Tan[x] on [0, b], 0 < b <= Pi/4:
 *   Integrate[f(x), {x, 0, b}] = Integrate[f(ArcTan t)/(1 + t^2), {t, 0, Tan[b]}],
 * with Cos[x] = (1 + t^2)^(-1/2), Sin[x] = t (1 + t^2)^(-1/2),
 * Cos[2x] = (1 - t^2)/(1 + t^2), Sin[2x] = 2t/(1 + t^2); then t -> Tan[b] t
 * maps the interval onto [0, 1].  Roach's 2008 #9, Cos[x]^(-2a-2) Cos[2x]^a on
 * [0, Pi/4], becomes (1 - t^2)^a, a Half-Beta.
 *
 * Why b <= Pi/4: on 0 < t < Tan[b] <= 1 every base the substitution produces
 * (t, 1 + t^2, 1 - t^2) is positive, so PowerExpand is exact; past Pi/4 the
 * factor 1 - t^2 changes sign and splitting a power of it would be wrong.  The
 * exponents are then Expanded, so that (1 + t^2)^(-1 - a + (2 + 2a)/2)
 * collapses to 1: left as it is, the Beta stages see a spurious (1 + t^2)
 * factor and decline.  Declines unless one of the Beta stages closes. */
Expr* integrate_tanbeta_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    if (!contains_symbol(f, x) || contains_symbol(b, x)) return NULL;
    if (!contains_head(f, "Cos") && !contains_head(f, "Sin") && !contains_head(f, "Tan") &&
        !contains_head(f, "Sec") && !contains_head(f, "Csc") && !contains_head(f, "Cot"))
        return NULL;
    if (!has_symbolic_exponent(f, x)) return NULL;
    if (!bt_true(bt_x("TrueQ[0 < #1 <= Pi/4]", 1, b))) return NULL;

    Expr* t = expr_new_symbol("Integrate`BetaR`t");
    Expr* g = bt_x(
        "With[{g = (#1 /. {Cos[2 #2] -> (1 - #3^2)/(1 + #3^2), "
        "Sin[2 #2] -> 2 #3/(1 + #3^2), Cos[#2] -> 1/Sqrt[1 + #3^2], "
        "Sin[#2] -> #3/Sqrt[1 + #3^2], Tan[#2] -> #3, Sec[#2] -> Sqrt[1 + #3^2], "
        "Csc[#2] -> Sqrt[1 + #3^2]/#3, Cot[#2] -> 1/#3})/(1 + #3^2)}, "
        "If[!FreeQ[g, #2], $Failed, "
        "PowerExpand[Tan[#4] (g /. #3 -> Tan[#4] #3)] /. "
        "Power[Integrate`BetaR`u_, Integrate`BetaR`e_] :> "
        "Power[Integrate`BetaR`u, Expand[Integrate`BetaR`e]]]]",
        4, f, x, t, b);
    Expr* r = NULL;
    if (g && contains_symbol(g, t)) r = beta_family_01(g, t, assumptions);
    if (g) expr_free(g);
    expr_free(t);
    return beta_to_gamma(r, assumptions);
}

/* Integration by parts on [0, 1] for C x^p u(x), u = ArcSin[x] or ArcCos[x]:
 *   = C u(1)/(p + 1) - C/(p + 1) Integrate[x^(p+1) u'(x), {x, 0, 1}],
 * u' = +-(1 - x^2)^(-1/2), so the remainder is a Half-Beta.  The boundary term
 * at 0 vanishes, and the integral converges, exactly when Re[p] > -1: proved
 * under the assumptions (bare value), refuted (decline), or undecided (a
 * ConditionalExpression on Re[p + 1] > 0).  Roach's 2008 #7, ArcSin[x]
 * x^(2n-1), closes in milliseconds instead of after a 55 s indefinite search
 * and a 125 s DiffUnderInt descent. */
Expr* integrate_invtrig_beta_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    if (!(b->type == EXPR_INTEGER && b->data.integer == 1)) return NULL;
    if (!contains_head(f, "ArcSin") && !contains_head(f, "ArcCos")) return NULL;

    Expr** fac; size_t nf; Expr* one[1];
    if (head_name_is(f, "Times")) { fac = f->data.function.args; nf = f->data.function.arg_count; }
    else { one[0] = f; fac = one; nf = 1; }
    int kind = 0;                      /* 1 ArcSin, 2 ArcCos */
    Expr* p = mk_int(0); Expr* C = mk_int(1);
    bool ok = true;
    for (size_t i = 0; i < nf && ok; i++) {
        Expr* g = fac[i];
        if (free_of_x(g, x)) { C = Times_(C, cp(g)); continue; }
        if (is_symbol(g, x)) { p = Plus_(p, mk_int(1)); continue; }
        if (kind == 0 && (head_name_is(g, "ArcSin") || head_name_is(g, "ArcCos")) &&
            g->data.function.arg_count == 1 && is_symbol(g->data.function.args[0], x)) {
            kind = head_name_is(g, "ArcSin") ? 1 : 2;
            continue;
        }
        if (head_name_is(g, "Power") && g->data.function.arg_count == 2 &&
            is_symbol(g->data.function.args[0], x) && free_of_x(g->data.function.args[1], x)) {
            p = Plus_(p, cp(g->data.function.args[1]));
            continue;
        }
        ok = false;
    }
    if (!ok || !kind) { expr_free(p); expr_free(C); return NULL; }
    Expr* p1 = simp(Plus_(p, mk_int(1)));
    /* a numeric exponent is Newton-Leibniz's: its antiderivative is elementary */
    if (!p1 || bt_true(ev1("NumericQ", cp(p1)))) {
        if (p1) expr_free(p1);
        expr_free(C);
        return NULL;
    }

    Expr* conv = mk_fn2("Greater", cp(p1), mk_int(0));
    int pv = prove(cp(conv), assumptions);
    if (pv == 0 && assumptions && prove_ref(cp(conv), assumptions)) pv = 1;
    if (pv == 0 && assumptions &&
        prove_ref(mk_fn2("LessEqual", cp(p1), mk_int(0)), assumptions)) pv = -1;
    if (pv == -1) { expr_free(conv); expr_free(p1); expr_free(C); return NULL; }
    Expr* as2 = assumptions ? mk_fn2("And", cp(assumptions), cp(conv)) : cp(conv);

    /* remainder  s C x^(p+1) (1 - x^2)^(-1/2)/(p + 1):  s = -1 (ArcSin), +1 (ArcCos) */
    Expr* rem = eval_take(mk_fn("Times", (Expr*[]){
        mk_int(kind == 1 ? -1 : 1), cp(C), mk_fn2("Power", cp(p1), mk_int(-1)),
        mk_fn2("Power", cp(x), cp(p1)),
        mk_fn2("Power", Plus_(mk_int(1), Times_(mk_int(-1), mk_fn2("Power", cp(x), mk_int(2)))),
               mk_fn2("Rational", mk_int(-1), mk_int(2))) }, 5));
    Expr* zero = mk_int(0); Expr* onei = mk_int(1);
    Expr* R = rem ? integrate_halfbeta_try(rem, x, zero, onei, as2) : NULL;
    expr_free(zero); expr_free(onei);
    if (rem) expr_free(rem);
    if (R && (head_name_is(R, "ConditionalExpression") || mentions_nonfinite(R))) {
        expr_free(R); R = NULL;
    }
    Expr* out = NULL;
    if (R) {
        /* boundary term C u(1)/(p + 1): u(1) = Pi/2 (ArcSin), 0 (ArcCos) */
        Expr* bterm = kind == 1
            ? Times_(cp(C), Times_(mk_fn2("Rational", mk_int(1), mk_int(2)),
                     Times_(expr_new_symbol("Pi"), mk_fn2("Power", cp(p1), mk_int(-1)))))
            : mk_int(0);
        Expr* v = beta_to_gamma(eval_take(Plus_(bterm, R)), as2);
        if (v && !mentions_nonfinite(v)) {
            if (pv == 1) out = v;
            else {
                Expr* cs = simp2(mk_fn2("Greater", mk_fn1("Re", cp(p1)), mk_int(0)), assumptions);
                out = eval_take(mk_fn2("ConditionalExpression", v, cs));
            }
        } else if (v) expr_free(v);
    }
    expr_free(as2); expr_free(conv); expr_free(p1); expr_free(C);
    return out;
}

/* ---- Euler finite-interval -> Beta * Hypergeometric2F1 ------------------ */
/* Integrate[x^(a-1) (1-x)^(b-1) (alpha + beta x)^e, {x,0,1}]
 *   = alpha^e Beta[a,b] Hypergeometric2F1[-e, a, a+b, -beta/alpha],
 * the Euler integral representation of 2F1.  Requires the extra linear factor
 * (alpha + beta x); a pure Beta integrand is left to integrate_beta_try.  Gated
 * on Re a>0, Re b>0 (convergence) and alpha>0 && alpha+beta>0 (the base stays
 * positive on [0,1], which also puts the 2F1 argument off the branch cut). */

static Expr* mk_fn3b(const char* h, Expr* a, Expr* b, Expr* c) { Expr* v[3]={a,b,c}; return mk_fn(h,v,3); }
static Expr* mk_fn4b(const char* h, Expr* a, Expr* b, Expr* c, Expr* d) { Expr* v[4]={a,b,c,d}; return mk_fn(h,v,4); }
static Expr* Pow_(Expr* a, Expr* b) { return mk_fn2("Power", a, b); }
static Expr* coeff_(const Expr* e, const Expr* x, long k) {
    return eval_take(mk_fn3b("Coefficient", cp(e), cp((Expr*)x), mk_int(k)));
}
static bool is_zero_simp_(Expr* e) {   /* consumes e */
    Expr* s = simp(e);
    bool z = (s->type == EXPR_INTEGER && s->data.integer == 0);
    expr_free(s);
    return z;
}

typedef struct {
    const Expr* x;
    Expr* aexp;        /* exponent of x        (owned) */
    Expr* bexp;        /* exponent of (1-x)    (owned) */
    Expr* C;           /* x-free constant      (owned) */
    Expr* lin_base;    /* captured extra-factor base (owned), or NULL */
    Expr* lin_exp;     /* captured extra-factor exponent (owned), or NULL */
    bool  ok;
} EulerCtx;

static bool euler_visit(Expr* g, void* vctx) {
    EulerCtx* c = (EulerCtx*)vctx;
    const Expr* x = c->x;
    if (free_of_x(g, x)) { c->C = Times_(c->C, cp(g)); return true; }
    if (is_symbol(g, x)) { c->aexp = Plus_(c->aexp, mk_int(1)); return true; }
    if (is_one_minus_x(g, x)) { c->bexp = Plus_(c->bexp, mk_int(1)); return true; }
    if (head_name_is(g, "Power") && g->data.function.arg_count == 2) {
        Expr* base = g->data.function.args[0];
        Expr* e    = g->data.function.args[1];
        if (!free_of_x(e, x)) { c->ok = false; return false; }
        if (is_symbol(base, x))      { c->aexp = Plus_(c->aexp, cp(e)); return true; }
        if (is_one_minus_x(base, x)) { c->bexp = Plus_(c->bexp, cp(e)); return true; }
        if (!c->lin_base) { c->lin_base = cp(base); c->lin_exp = cp(e); return true; }
    }
    c->ok = false;
    return false;
}

Expr* integrate_euler_2f1_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!contains_symbol(f, x)) return NULL;
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    if (!(b->type == EXPR_INTEGER && b->data.integer == 1)) return NULL;

    EulerCtx c = { x, mk_int(0), mk_int(0), mk_int(1), NULL, NULL, true };
    for_each_factor(f, euler_visit, &c);
    if (!c.ok || !c.lin_base) {
        expr_free(c.aexp); expr_free(c.bexp); expr_free(c.C);
        if (c.lin_base) expr_free(c.lin_base);
        if (c.lin_exp) expr_free(c.lin_exp);
        return NULL;
    }
    /* Validate the extra factor's base is genuinely linear alpha + beta x, beta != 0. */
    Expr* beta  = coeff_(c.lin_base, x, 1);
    Expr* alpha = coeff_(c.lin_base, x, 0);
    bool linear = free_of_x(alpha, x) && free_of_x(beta, x) &&
                  !(beta->type == EXPR_INTEGER && beta->data.integer == 0) &&
                  is_zero_simp_(Plus_(cp(c.lin_base),
                      Times_(mk_int(-1), Plus_(cp(alpha), Times_(cp(beta), cp((Expr*)x))))));
    Expr* e = cp(c.lin_exp);
    expr_free(c.lin_base); expr_free(c.lin_exp);
    if (!linear) {
        expr_free(c.aexp); expr_free(c.bexp); expr_free(c.C);
        expr_free(alpha); expr_free(beta); expr_free(e);
        return NULL;
    }
    Expr* A = simp(Plus_(c.aexp, mk_int(1)));   /* a = aexp + 1 */
    Expr* B = simp(Plus_(c.bexp, mk_int(1)));   /* b = bexp + 1 */

    /* Gates: a>0, b>0, alpha>0, alpha+beta>0 (base positive on [0,1]). */
    bool gate = prove_ref(mk_fn2("Greater", cp(A), mk_int(0)), assumptions) &&
                prove_ref(mk_fn2("Greater", cp(B), mk_int(0)), assumptions) &&
                prove_ref(mk_fn2("Greater", cp(alpha), mk_int(0)), assumptions) &&
                prove_ref(mk_fn2("Greater", Plus_(cp(alpha), cp(beta)), mk_int(0)), assumptions);
    if (!gate) {
        expr_free(c.C); expr_free(A); expr_free(B); expr_free(alpha); expr_free(beta); expr_free(e);
        return NULL;
    }
    /* value = C alpha^e Beta[a,b] Hypergeometric2F1[-e, a, a+b, -beta/alpha]. */
    Expr* two_f1 = mk_fn4b("Hypergeometric2F1",
        Times_(mk_int(-1), cp(e)), cp(A), Plus_(cp(A), cp(B)),
        Times_(mk_int(-1), Times_(cp(beta), Pow_(cp(alpha), mk_int(-1)))));
    Expr* val = Times_(c.C, Times_(Pow_(cp(alpha), cp(e)),
                    Times_(Beta_(cp(A), cp(B)), two_f1)));
    expr_free(A); expr_free(B); expr_free(alpha); expr_free(beta); expr_free(e);
    val = simp2(val, assumptions);
    if (!val || mentions_nonfinite(val)) { if (val) expr_free(val); return NULL; }
    return val;
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
            if (lhs->type == EXPR_SYMBOL && lhs->data.symbol.name == SYM_Assumptions) {
                *assumptions = opt->data.function.args[1];
                continue;
            }
        }
        return false;
    }
    return true;
}

/* Log[Sin[x]] / Log[Cos[x]] over a canonical interval starting at 0.  The
 * classical log-sine values:
 *   Integrate[Log[Sin[x]], {x,0,Pi/2}] = Integrate[Log[Cos[x]], {x,0,Pi/2}]
 *                                      = -(Pi/2) Log[2],
 *   Integrate[Log[Sin[x]], {x,0,Pi}]   = -Pi Log[2].
 * Log[Cos[x]] on [0,Pi] is complex (Cos < 0 on (Pi/2,Pi)), so it is declined.
 * The trig argument must be exactly x (coefficient 1). */
Expr* integrate_logtrig_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    (void)assumptions;
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    int kind = trig_interval(b);                 /* 1 = Pi/2, 2 = Pi, 4 = 2 Pi */
    if (kind != 1 && kind != 2) return NULL;
    if (!head_name_is(f, "Log") || f->data.function.arg_count != 1) return NULL;
    Expr* g = f->data.function.args[0];
    bool is_sin = head_name_is(g, "Sin");
    bool is_cos = head_name_is(g, "Cos");
    if ((!is_sin && !is_cos) || g->data.function.arg_count != 1) return NULL;
    if (!expr_eq(g->data.function.args[0], x)) return NULL;   /* argument exactly x */

    /* -(Pi/2) Log[2] on [0, Pi/2] for both Sin and Cos. */
    if (kind == 1)
        return simp(mk_fn2("Times", mk_fn2("Rational", mk_int(-1), mk_int(2)),
                           mk_fn2("Times", expr_new_symbol("Pi"), mk_fn1("Log", mk_int(2)))));
    /* -Pi Log[2] on [0, Pi] for Sin; Cos is complex there -> decline. */
    if (kind == 2 && is_sin)
        return simp(mk_fn2("Times", mk_int(-1),
                           mk_fn2("Times", expr_new_symbol("Pi"), mk_fn1("Log", mk_int(2)))));
    return NULL;
}

/* Poisson log integral over [0, Pi]:
 *   Integrate[Log[A + B Cos[x]], {x,0,Pi}] = Pi Log[(A + Sqrt[A^2 - B^2])/2],
 * valid for A > 0 and A >= |B| (otherwise Log[A + B Cos x] goes negative).  A, B
 * free of x.  E.g. Log[1 - 2a Cos x + a^2] = Log[(1+a^2) + (-2a) Cos x] gives 0
 * for 0<a<1 and 2 Pi Log[a] for a>1. */
Expr* integrate_poissonlog_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    if (trig_interval(b) != 2) return NULL;                 /* [0, Pi] only */
    if (!head_name_is(f, "Log") || f->data.function.arg_count != 1) return NULL;
    Expr* g = f->data.function.args[0];
    if (!contains_symbol(g, x)) return NULL;

    /* g must be exactly A + B Cos[x], with A, B free of x. */
    Expr* cosx = mk_fn1("Cos", cp(x));
    Expr* A = coeff_(g, cosx, 0);
    Expr* B = coeff_(g, cosx, 1);
    Expr* recon = mk_fn2("Plus", cp(A), mk_fn2("Times", cp(B), cp(cosx)));
    Expr* diff  = mk_fn2("Plus", cp(g), mk_fn2("Times", mk_int(-1), recon));
    bool ok_form = is_zero_simp_(diff) &&                   /* consumes diff */
                   !contains_symbol(A, x) && !contains_symbol(B, x);
    expr_free(cosx);
    if (!ok_form) { expr_free(A); expr_free(B); return NULL; }

    /* Soundness: A > 0 and A^2 >= B^2 (A >= |B|) must be proved. */
    bool gate = prove_ref(mk_fn2("Greater", cp(A), mk_int(0)), assumptions) &&
                prove_ref(mk_fn2("GreaterEqual", Pow_(cp(A), mk_int(2)),
                                  Pow_(cp(B), mk_int(2))), assumptions);
    if (!gate) { expr_free(A); expr_free(B); return NULL; }

    Expr* disc = mk_fn2("Plus", Pow_(cp(A), mk_int(2)),
                        mk_fn2("Times", mk_int(-1), Pow_(cp(B), mk_int(2))));
    Expr* root = Pow_(disc, mk_fn2("Rational", mk_int(1), mk_int(2)));   /* consumes disc */
    Expr* inner = mk_fn2("Times", mk_fn2("Plus", cp(A), root),
                         mk_fn2("Rational", mk_int(1), mk_int(2)));
    Expr* res = mk_fn2("Times", expr_new_symbol("Pi"), mk_fn1("Log", inner));
    expr_free(A); expr_free(B);
    return assumptions ? simp2(res, assumptions) : simp(res);
}

/* King's-property reflection: for f = x g(x) on [0,c] with g invariant under
 * x -> c - x, f(x) + f(c-x) = c g(x), so Integrate[x g, {x,0,c}] =
 * (c/2) Integrate[g, {x,0,c}].  E.g. Integrate[x Sin[x]/(1+Cos[x]^2), {x,0,Pi}]
 * = (Pi/2) Integrate[Sin[x]/(1+Cos[x]^2), {x,0,Pi}] = (Pi/2)(Pi/2) = Pi^2/4. */
Expr* integrate_kingreflect_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    /* The reflection x -> c - x is an involution of [0, c] only for a FINITE c;
     * on a half-line there is no such map.  Without this gate the invariance
     * test substituted x -> Infinity - x and Simplify'd the result, which raised
     * Infinity::indet (Infinity - Infinity) -- a spurious diagnostic that made
     * Check[] reject a half-line integral another stage then closed. */
    if (mentions_nonfinite(b)) return NULL;
    /* Cheap gate: f must be a product with an explicit bare x factor (f = x g). */
    if (!head_name_is(f, "Times")) return NULL;
    bool has_x = false;
    for (size_t i = 0; i < f->data.function.arg_count; i++)
        if (expr_eq(f->data.function.args[i], (Expr*)x)) { has_x = true; break; }
    if (!has_x) return NULL;

    Expr* g = eval_take(mk_fn1("Cancel", mk_fn2("Times", cp(f), Pow_(cp(x), mk_int(-1)))));
    if (!g || !contains_symbol(g, x)) { if (g) expr_free(g); return NULL; }

    /* g must be invariant under x -> b - x. */
    Expr* refl = eval_take(mk_fn2("ReplaceAll", cp(g),
                     mk_fn2("Rule", cp(x),
                            mk_fn2("Plus", cp(b), mk_fn2("Times", mk_int(-1), cp(x))))));
    bool invariant = is_zero_simp_(mk_fn2("Plus", refl,
                                   mk_fn2("Times", mk_int(-1), cp(g))));   /* consumes */
    if (!invariant) { expr_free(g); return NULL; }

    /* I = (b/2) Integrate[g, {x,0,b}]. */
    Expr* spec = mk_fn3b("List", cp(x), mk_int(0), cp(b));
    Expr* inner = assumptions
        ? eval_take(mk_fn3b("Integrate", cp(g), spec,
                            mk_fn2("Rule", expr_new_symbol("Assumptions"), cp(assumptions))))
        : eval_take(mk_fn2("Integrate", cp(g), spec));
    expr_free(g);
    if (!inner || head_name_is(inner, "Integrate") || contains_symbol(inner, x)) {
        if (inner) expr_free(inner);
        return NULL;
    }
    Expr* res = mk_fn2("Times",
                       mk_fn2("Times", cp(b), mk_fn2("Rational", mk_int(1), mk_int(2))), inner);
    return assumptions ? simp2(res, assumptions) : simp(res);
}

/* Log-weighted rational on [0,1]:
 *   Integrate[x^(c-1) Log[x]^n / (1 - x^m), {x,0,1}]
 *     = (-1)^n n! / m^(n+1) * HurwitzZeta[n+1, c/m]   (c > 0, m >= 1, n >= 1).
 * This is the l = 0 limit the Beta stage cannot represent (the n = 0 base
 * integral diverges, so the value is not a k-derivative of a finite Beta;
 * Beta[c, 0] = ComplexInfinity makes beta_param_deriv Indeterminate).  Running
 * before the Beta stage both closes the family and avoids that spurious
 * Infinity::indet. */
Expr* integrate_ratlog01_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!(a->type == EXPR_INTEGER && a->data.integer == 0)) return NULL;
    if (!(b->type == EXPR_INTEGER && b->data.integer == 1)) return NULL;

    /* Strip Log[x]^n (n >= 1); collect the remaining factors. */
    int n = 0;
    Expr* restparts[32]; size_t nrest = 0;
    Expr** fac; size_t nf; Expr* one[1];
    if (head_name_is(f, "Times")) { fac = f->data.function.args; nf = f->data.function.arg_count; }
    else { one[0] = f; fac = one; nf = 1; }
    for (size_t i = 0; i < nf; i++) {
        Expr* fe = fac[i];
        if (head_name_is(fe, "Log") && fe->data.function.arg_count == 1 &&
            expr_eq(fe->data.function.args[0], (Expr*)x)) { n += 1; continue; }
        if (head_name_is(fe, "Power") && fe->data.function.arg_count == 2 &&
            head_name_is(fe->data.function.args[0], "Log") &&
            fe->data.function.args[0]->data.function.arg_count == 1 &&
            expr_eq(fe->data.function.args[0]->data.function.args[0], (Expr*)x) &&
            fe->data.function.args[1]->type == EXPR_INTEGER &&
            fe->data.function.args[1]->data.integer >= 1) {
            n += (int) fe->data.function.args[1]->data.integer; continue;
        }
        if (nrest < 32) restparts[nrest++] = fe;
    }
    if (n < 1) return NULL;

    Expr* rest;
    if (nrest == 0) rest = mk_int(1);
    else if (nrest == 1) rest = cp(restparts[0]);
    else {
        Expr** c2 = malloc(nrest * sizeof(Expr*));
        for (size_t i = 0; i < nrest; i++) c2[i] = cp(restparts[i]);
        rest = mk_fn("Times", c2, nrest);
        free(c2);
    }

    /* rest should be C x^(c-1) / (1 - x^m).  Together normalizes the sign (e.g.
     * 1/(1-x) -> -1/(-1+x)), so match den == D0 (1 - x^m) with D0 the (nonzero,
     * x-free) constant term, folding D0 into C. */
    Expr* tog = eval_take(mk_fn1("Together", rest));
    Expr* num = eval_take(mk_fn1("Numerator", cp(tog)));
    Expr* den = eval_take(mk_fn1("Denominator", tog));

    Expr* mE = eval_take(mk_fn2("Exponent", cp(den), cp(x)));
    long m = (mE && mE->type == EXPR_INTEGER) ? (long) mE->data.integer : 0;
    if (mE) expr_free(mE);
    Expr* D0 = coeff_(den, x, 0);                     /* constant term of den */
    bool den_ok = m >= 1 && D0 && !contains_symbol(D0, x) &&
        !(D0->type == EXPR_INTEGER && D0->data.integer == 0) &&
        is_zero_simp_(mk_fn2("Plus", cp(den),
            mk_fn2("Times", mk_int(-1),
                mk_fn2("Times", cp(D0),
                    mk_fn2("Plus", mk_int(1), mk_fn2("Times", mk_int(-1), Pow_(cp(x), mk_int(m))))))));

    /* numC = num / D0 == C x^p, C free of x; c = p + 1. */
    Expr* numC = den_ok ? eval_take(mk_fn2("Times", cp(num), Pow_(cp(D0), mk_int(-1)))) : NULL;
    Expr* pE = numC ? eval_take(mk_fn2("Exponent", cp(numC), cp(x))) : NULL;
    long p = (pE && pE->type == EXPR_INTEGER) ? (long) pE->data.integer : -1;
    if (pE) expr_free(pE);
    Expr* C = (p >= 0) ? eval_take(mk_fn2("Times", cp(numC), Pow_(cp(x), mk_int(-p)))) : NULL;
    bool num_ok = p >= 0 && C && !contains_symbol(C, x) &&
        is_zero_simp_(mk_fn2("Plus", cp(numC),
            mk_fn2("Times", mk_int(-1), mk_fn2("Times", cp(C), Pow_(cp(x), mk_int(p))))));

    expr_free(num); expr_free(den); if (D0) expr_free(D0); if (numC) expr_free(numC);
    if (!den_ok || !num_ok) { if (C) expr_free(C); return NULL; }

    /* value = C (-1)^n n! / m^(n+1) HurwitzZeta[n+1, (p+1)/m]. */
    Expr* hz = mk_fn2("HurwitzZeta", mk_int(n + 1),
                   mk_fn2("Times", mk_int(p + 1), Pow_(mk_int(m), mk_int(-1))));
    Expr* val = mk_fn2("Times", C,
                    mk_fn2("Times", (n % 2 == 0) ? mk_int(1) : mk_int(-1),
                        mk_fn2("Times", mk_fn1("Factorial", mk_int(n)),
                            mk_fn2("Times", Pow_(mk_int(m), mk_int(-(n + 1))), hz))));
    return assumptions ? simp2(val, assumptions) : simp(val);
}

/* Raabe's integral: Integrate[Log[Gamma[x]], {x, t, t+1}] = (1/2)Log[2 Pi] +
 * t Log[t] - t  for t >= 0 (the canonical [0,1] instance is (1/2)Log[2 Pi]).
 * The Gamma argument must be exactly the integration variable and the interval a
 * unit interval [t, t+1] with t >= 0; t = 0 is returned directly (the t Log t
 * term is the removable 0 limit). */
Expr* integrate_raabe_try(Expr* f, Expr* x, Expr* a, Expr* b, Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;
    if (!head_name_is(f, "Log") || f->data.function.arg_count != 1) return NULL;
    Expr* g = f->data.function.args[0];
    if (!head_name_is(g, "Gamma") || g->data.function.arg_count != 1) return NULL;
    if (!expr_eq(g->data.function.args[0], (Expr*)x)) return NULL;
    /* Unit interval: b - a == 1. */
    if (!is_zero_simp_(mk_fn2("Plus", cp(b),
             mk_fn2("Plus", mk_fn2("Times", mk_int(-1), cp(a)), mk_int(-1)))))
        return NULL;
    Expr* half_log2pi = mk_fn2("Times", mk_fn2("Rational", mk_int(1), mk_int(2)),
                            mk_fn1("Log", mk_fn2("Times", mk_int(2), expr_new_symbol("Pi"))));
    /* t = 0: the canonical value. */
    if (a->type == EXPR_INTEGER && a->data.integer == 0) return simp(half_log2pi);
    /* t > 0: require it proved (numeric or under assumptions). */
    bool t_pos = (a->type == EXPR_INTEGER && a->data.integer > 0) ||
                 (a->type == EXPR_REAL && a->data.real > 0.0);
    if (!t_pos && assumptions) {
        Expr* q = eval_take(mk_fn2("Refine", mk_fn2("Greater", cp(a), mk_int(0)), cp(assumptions)));
        t_pos = q && q->type == EXPR_SYMBOL && q->data.symbol.name == SYM_True;
        if (q) expr_free(q);
    }
    if (!t_pos) { expr_free(half_log2pi); return NULL; }
    Expr* val = mk_fn2("Plus", half_log2pi,
                    mk_fn2("Plus", mk_fn2("Times", cp(a), mk_fn1("Log", cp(a))),
                           mk_fn2("Times", mk_int(-1), cp(a))));
    return assumptions ? simp2(val, assumptions) : simp(val);
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
    symtab_get_def("Integrate`Beta")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`Beta",
        "Integrate`Beta[f, {x, 0, 1}] evaluates a definite integral of the form "
        "x^(k-1) (1-x)^(l-1) over [0,1] as the Euler Beta function Beta[k, l], "
        "and an integrand additionally weighted by Log[x]^i Log[1-x]^j as the "
        "mixed parameter derivative of Beta.  Gated on the convergence strip "
        "Re[k] > 0 && Re[l] > 0 (an undecided strip yields a "
        "ConditionalExpression).  Returns unevaluated when the integrand is not "
        "of this form or the interval is not [0,1].");

    symtab_add_builtin("Integrate`TrigPower", builtin_integrate_trigpower);
    symtab_get_def("Integrate`TrigPower")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`TrigPower",
        "Integrate`TrigPower[f, {x, 0, c}] evaluates Sin[x]^m Cos[x]^n over a "
        "canonical interval: over [0, Pi/2] it is Beta[(m+1)/2, (n+1)/2]/2; over "
        "[0, Pi] it is 0 when n is an odd integer and twice the quarter-period "
        "value when n is even; over [0, 2 Pi] it is 0 unless both m and n are "
        "even integers, and four times the quarter-period value otherwise.  "
        "Gated on Re[(m+1)/2] > 0 && Re[(n+1)/2] > 0.  Returns unevaluated when "
        "the integrand or interval is not of this form.");
}
