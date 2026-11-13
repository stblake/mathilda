/* integrate_piecewise.c
 *
 * Definite integration of piecewise integrands.  See integrate_piecewise.h for
 * the contract.  The integrand is a combination of the branch-defining heads
 * Abs / Sign / HeavisideTheta / UnitStep / Min / Max / Piecewise whose
 * arguments depend on the integration variable.  We find the breakpoints where
 * a branch switches, split [a, b] at the interior ones, reduce each sub-interval
 * to a single smooth branch (chosen by a *proved* sign of each governing
 * argument at an interior test point), integrate that with the ordinary
 * cascade, and sum.  Anything undecidable yields NULL (left unevaluated).
 */

#include "integrate_piecewise.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"

#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* ---- small construction / evaluation helpers ---------------------------- */

static Expr* cp(const Expr* e) { return expr_copy((Expr*)e); }

static Expr* mk_fn(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}
static Expr* mk_fn1(const char* head, Expr* a) {
    Expr* args[1] = { a }; return mk_fn(head, args, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    Expr* args[2] = { a, b }; return mk_fn(head, args, 2);
}
static Expr* mk_fn3(const char* head, Expr* a, Expr* b, Expr* c) {
    Expr* args[3] = { a, b, c }; return mk_fn(head, args, 3);
}
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }
static Expr* neg_of(Expr* g) { return mk_fn2("Times", expr_new_integer(-1), g); }

/* Evaluate `call` (taking ownership) and return the result. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);   /* evaluate does not free its input */
    expr_free(call);
    return r;
}

static bool head_name_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}
static bool sym_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

static bool contains_symbol(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return x->type == EXPR_SYMBOL && e->data.symbol.name == x->data.symbol.name;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol(e->data.function.args[i], x)) return true;
    return false;
}
static bool contains_symbol_name(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return strcmp(e->data.symbol.name, name) == 0;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_symbol_name(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_symbol_name(e->data.function.args[i], name)) return true;
    return false;
}
static bool contains_head(const Expr* e, const char* name) {
    if (!e) return false;
    if (head_name_is(e, name)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_head(e->data.function.head, name)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_head(e->data.function.args[i], name)) return true;
    return false;
}

/* ---- infinities (bounds may be +/-Infinity) ----------------------------- */

static bool is_pos_infinity(const Expr* b) {
    if (sym_is(b, "Infinity")) return true;
    return head_name_is(b, "DirectedInfinity") && b->data.function.arg_count == 1 &&
           b->data.function.args[0]->type == EXPR_INTEGER &&
           b->data.function.args[0]->data.integer == 1;
}
static bool is_neg_infinity(const Expr* a) {
    if (head_name_is(a, "DirectedInfinity") && a->data.function.arg_count == 1 &&
        a->data.function.args[0]->type == EXPR_INTEGER &&
        a->data.function.args[0]->data.integer == -1)
        return true;
    if (head_name_is(a, "Times") && a->data.function.arg_count == 2) {
        Expr* p = a->data.function.args[0];
        Expr* q = a->data.function.args[1];
        if (p->type == EXPR_INTEGER && p->data.integer == -1 && is_pos_infinity(q))
            return true;
    }
    return false;
}

/* ---- assumption-aware predicate proving --------------------------------- */

/* True iff Refine[pred, as] proves `pred` to the symbol True.  Consumes pred. */
static bool refine_true(Expr* pred, const Expr* as) {
    Expr* call = as ? mk_fn2("Refine", pred, cp(as)) : mk_fn1("Refine", pred);
    Expr* r = eval_take(call);
    bool t = sym_is(r, "True");
    if (r) expr_free(r);
    return t;
}

/* ReplaceAll[e, x -> val], evaluated. */
static Expr* subst(const Expr* e, const Expr* x, const Expr* val) {
    Expr* rule = mk_fn2("Rule", cp(x), cp(val));
    return eval_take(mk_fn2("ReplaceAll", cp(e), rule));
}

/* Is `lo` < `c` provable (lo, c may be +/-Infinity)? */
static bool lt_bound(const Expr* lo, const Expr* c, const Expr* as) {
    if (is_neg_infinity(lo)) return true;
    if (is_pos_infinity(lo)) return false;
    if (is_pos_infinity(c))  return true;
    if (is_neg_infinity(c))  return false;
    return refine_true(mk_fn2("Less", cp(lo), cp(c)), as);
}

/* ---- breakpoint set ----------------------------------------------------- */

/* `incomplete`: some breakpoint could not be placed -- a periodic root family
 * on an infinite (or non-numeric) range, a Solve that did not answer, a root
 * that still mentions x.  The pieces between the known breakpoints are then
 * not branch-constant, and resolving one at its test point would silently pick
 * one branch for the whole piece: Integrate[Abs[Sin[x]] E^-x, {x, 0, Infinity}]
 * was 1/2, the integral of Sin[x] E^-x (G&R 17.13.34-35). */
typedef struct { Expr** pts; int n, cap; bool incomplete;
                 /* a periodic family's numeric period, when one was seen */
                 double period; bool has_period; Expr* period_exact; } BPSet;

static void bp_init(BPSet* s) { s->pts = NULL; s->n = 0; s->cap = 0; s->incomplete = false;
                                 s->period = 0.0; s->has_period = false;
                                 s->period_exact = NULL; }
static void bp_free(BPSet* s) {
    if (s->period_exact) { expr_free(s->period_exact); s->period_exact = NULL; }
    for (int i = 0; i < s->n; i++) expr_free(s->pts[i]);
    free(s->pts);
    s->pts = NULL; s->n = 0; s->cap = 0;
}
/* Add `p` (takes ownership); drop an exact structural duplicate. */
static void bp_add(BPSet* s, Expr* p) {
    for (int i = 0; i < s->n; i++)
        if (expr_eq(s->pts[i], p)) { expr_free(p); return; }
    if (s->n == s->cap) {
        s->cap = s->cap ? s->cap * 2 : 4;
        s->pts = realloc(s->pts, (size_t)s->cap * sizeof(Expr*));
    }
    s->pts[s->n++] = p;
}

/* A concrete non-integer exponent: a Rational p/q with |q| != 1, or a
 * non-integer real.  This is the fingerprint of a radical Power[g, p].  A
 * SYMBOLIC exponent (a, -b, ...) is deliberately NOT fractional here, so
 * Abs[x]^a / Abs[x]^(-b) stay on the branch-head path and are never treated as
 * radicals. */
static bool is_frac_exponent(const Expr* e) {
    if (head_name_is(e, "Rational") && e->data.function.arg_count == 2) {
        Expr* d = e->data.function.args[1];
        if (d->type == EXPR_INTEGER) return d->data.integer != 1 && d->data.integer != -1;
    }
    if (e->type == EXPR_REAL) {
        double v = e->data.real;
        return v != floor(v);
    }
    return false;
}

/* N[e] as a finite machine real; false if not a concrete real. */
static bool num_double(const Expr* e, double* out) {
    Expr* nv = eval_take(mk_fn1("N", cp(e)));
    bool ok = false;
    if (nv) {
        if (nv->type == EXPR_REAL) { *out = nv->data.real; ok = true; }
        else if (nv->type == EXPR_INTEGER) { *out = (double)nv->data.integer; ok = true; }
        expr_free(nv);
    }
    return ok;
}

/* The integer-constrained parameter of a Solve condition, e.g. C[1] in
 * Element[C[1], Integers] (recursing through And/Or wrappers).  NULL if none. */
static const Expr* integer_param_of(const Expr* cond) {
    if (!cond || cond->type != EXPR_FUNCTION) return NULL;
    if (head_name_is(cond, "Element") && cond->data.function.arg_count == 2 &&
        sym_is(cond->data.function.args[1], "Integers"))
        return cond->data.function.args[0];
    for (size_t i = 0; i < cond->data.function.arg_count; i++) {
        const Expr* p = integer_param_of(cond->data.function.args[i]);
        if (p) return p;
    }
    return NULL;
}

/* Enumerate the concrete roots of a periodic family body(param), param an
 * integer, that fall strictly inside (a, b), appending each as an exact Expr.
 * Needs numeric a, b and a numeric period; a no-op otherwise.  This is what
 * turns Solve[Sin[x] == 0] -> {Pi + 2 C[1] Pi, 2 C[1] Pi} into the actual
 * breakpoints {Pi} on (0, 2 Pi) -- the periodic family Solve cannot localize by
 * itself. */
static void enum_periodic_root(const Expr* body, const Expr* param,
                               const Expr* a, const Expr* b, BPSet* s) {
    double av, bv;
    bool afin = num_double(a, &av), bfin = num_double(b, &bv);
    Expr* z = expr_new_integer(0);
    Expr* one = expr_new_integer(1);
    Expr* base = subst(body, param, z);
    Expr* b1   = subst(body, param, one);
    expr_free(z); expr_free(one);
    double basev, b1v;
    bool okb = base && num_double(base, &basev);
    bool ok1 = b1 && num_double(b1, &b1v);
    if (!okb || !ok1) {
        /* A symbolic family (x = Pi/a + 2 Pi k/a).  The period is recorded for
         * the half-line rule, which proves its sign under the assumptions.  On a
         * finite range the members are still placeable by RATIO: base + k P lies
         * strictly inside (a, b) iff k lies strictly between (a - base)/P and
         * (b - base)/P, and those ratios are numbers whenever the range is
         * measured in periods ([0, 2 Pi/a] for Sin[a x]). */
        Expr* P = (base && b1) ? eval_take(mk_fn2("Plus", cp(b1), neg_of(cp(base)))) : NULL;
        if (P && !s->period_exact) s->period_exact = cp(P);
        bool placed = false;
        if (P && !is_pos_infinity(a) && !is_neg_infinity(a) &&
            !is_pos_infinity(b) && !is_neg_infinity(b)) {
            Expr* r0 = eval_take(mk_fn1("Simplify", mk_fn2("Times", mk_fn2("Plus", cp(a), neg_of(cp(base))),
                                          mk_fn2("Power", cp(P), expr_new_integer(-1)))));
            Expr* r1 = eval_take(mk_fn1("Simplify", mk_fn2("Times", mk_fn2("Plus", cp(b), neg_of(cp(base))),
                                          mk_fn2("Power", cp(P), expr_new_integer(-1)))));
            double k0, k1;
            if (r0 && r1 && num_double(r0, &k0) && num_double(r1, &k1) &&
                fabs(k1 - k0) < 100000.0) {
                double lo = fmin(k0, k1), hi = fmax(k0, k1);
                for (long k = (long)floor(lo) - 1; k <= (long)ceil(hi) + 1; k++) {
                    if ((double)k <= lo + 1e-9 || (double)k >= hi - 1e-9) continue;
                    Expr* kk = expr_new_integer((int64_t)k);
                    Expr* root = subst(body, param, kk);
                    expr_free(kk);
                    if (root) bp_add(s, root);
                }
                placed = true;
            }
            if (r0) expr_free(r0);
            if (r1) expr_free(r1);
        }
        if (P) expr_free(P);
        if (base) expr_free(base);
        if (b1) expr_free(b1);
        if (!placed) s->incomplete = true;
        return;
    }
    if (base) expr_free(base);
    if (b1) expr_free(b1);

    double period = b1v - basev;
    if (fabs(period) < 1e-12) { s->incomplete = true; return; }
    if (!s->has_period || fabs(fabs(period) - s->period) > 1e-12) {
        /* exact period |body(1) - body(0)|; two families of different
         * periods leave the one seen last (the caller checks periodicity) */
        Expr* z0 = expr_new_integer(0), *o1 = expr_new_integer(1);
        Expr* e0 = subst(body, param, z0), *e1 = subst(body, param, o1);
        expr_free(z0); expr_free(o1);
        if (e0 && e1) {
            if (s->period_exact) expr_free(s->period_exact);
            s->period_exact = eval_take(mk_fn1("Abs", mk_fn2("Plus", e1, neg_of(e0))));
        } else { if (e0) expr_free(e0); if (e1) expr_free(e1); }
    }
    s->period = fabs(period); s->has_period = true;
    if (!afin || !bfin) { s->incomplete = true; return; }   /* infinitely many */
    if (bv < av) { double t = av; av = bv; bv = t; }
    double k0 = (av - basev) / period, k1 = (bv - basev) / period;
    double klo = floor(fmin(k0, k1)) - 1.0, khi = ceil(fmax(k0, k1)) + 1.0;
    if (khi - klo > 100000.0) { s->incomplete = true; return; }   /* absurd range */
    for (long k = (long)klo; k <= (long)khi; k++) {
        double pt = basev + period * (double)k;
        if (pt <= av + 1e-9 || pt >= bv - 1e-9) continue;   /* strict interior */
        Expr* kk = expr_new_integer((int64_t)k);
        Expr* root = subst(body, param, kk);                /* exact symbolic root */
        expr_free(kk);
        if (root) bp_add(s, root);
    }
}

/* Real roots of g == 0 in x, appended to s.  Complex roots and roots that still
 * mention x are skipped; a periodic ConditionalExpression family is expanded to
 * its concrete members inside (a, b) via enum_periodic_root. */
static void solve_roots(const Expr* g, const Expr* x, const Expr* a,
                        const Expr* b, BPSet* s) {
    Expr* eq  = mk_fn2("Equal", cp(g), expr_new_integer(0));
    /* Quiet: Solve is an internal probe here; its Solve::ifun note on a
     * transcendental breakpoint argument must not surface from Integrate. */
    Expr* sol = eval_take(mk_fn1("Quiet", mk_fn2("Solve", eq, cp(x))));
    if (head_name_is(sol, "List")) {
        for (size_t i = 0; i < sol->data.function.arg_count; i++) {
            Expr* si = sol->data.function.args[i];
            if (!head_name_is(si, "List")) continue;
            for (size_t j = 0; j < si->data.function.arg_count; j++) {
                Expr* rule = si->data.function.args[j];
                if (!head_name_is(rule, "Rule") ||
                    rule->data.function.arg_count != 2) continue;
                if (!expr_eq(rule->data.function.args[0], (Expr*)x)) continue;
                Expr* root = rule->data.function.args[1];
                if (head_name_is(root, "ConditionalExpression") &&
                    root->data.function.arg_count == 2) {
                    Expr* body = root->data.function.args[0];
                    const Expr* param = integer_param_of(root->data.function.args[1]);
                    if (param && !contains_symbol(body, x) &&
                        !contains_head(body, "Complex"))
                        enum_periodic_root(body, param, a, b, s);
                    else if (!contains_head(body, "Complex"))
                        s->incomplete = true;
                    continue;
                }
                if (contains_head(root, "Complex")) continue;   /* not real */
                if (contains_symbol(root, x)) { s->incomplete = true; continue; }
                bp_add(s, cp(root));
            }
        }
    } else {
        s->incomplete = true;                       /* Solve did not answer */
    }
    if (sol) expr_free(sol);
}

/* Breakpoints contributed by a condition `cond` (relational boundaries). */
static void collect_condition_boundaries(const Expr* cond, const Expr* x,
                                         const Expr* a, const Expr* b, BPSet* s) {
    if (!cond || cond->type != EXPR_FUNCTION) return;
    const char* h = cond->data.function.head->type == EXPR_SYMBOL
                      ? cond->data.function.head->data.symbol.name : "";
    if ((h == SYM_Less || h == SYM_LessEqual ||
         h == SYM_Greater || h == SYM_GreaterEqual ||
         h == SYM_Equal || h == SYM_Unequal) &&
        cond->data.function.arg_count == 2) {
        Expr* diff = mk_fn2("Plus", cp(cond->data.function.args[0]),
                            neg_of(cp(cond->data.function.args[1])));
        if (contains_symbol(diff, x)) solve_roots(diff, x, a, b, s);
        expr_free(diff);
        return;
    }
    /* And / Or / Inequality / Not: recurse into the pieces. */
    for (size_t i = 0; i < cond->data.function.arg_count; i++)
        collect_condition_boundaries(cond->data.function.args[i], x, a, b, s);
}

/* Walk the integrand, collecting every branch-switch breakpoint in x.  A
 * fractional-power radical Power[g, p] contributes the real zeros of its
 * radicand g: the radical is smooth except where g vanishes, so those zeros are
 * kink points exactly like an Abs argument's zeros. */
static void collect_breakpoints(const Expr* e, const Expr* x, const Expr* a,
                                const Expr* b, BPSet* s) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const char* h = e->data.function.head->type == EXPR_SYMBOL
                      ? e->data.function.head->data.symbol.name : "";
    size_t n = e->data.function.arg_count;
    if ((h == SYM_Abs || h == SYM_Sign ||
         h == SYM_HeavisideTheta || h == SYM_UnitStep) &&
        n == 1 && contains_symbol(e->data.function.args[0], x)) {
        solve_roots(e->data.function.args[0], x, a, b, s);
    } else if ((h == SYM_Min || h == SYM_Max) && n >= 2) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = i + 1; j < n; j++) {
                Expr* d = mk_fn2("Plus", cp(e->data.function.args[i]),
                                 neg_of(cp(e->data.function.args[j])));
                if (contains_symbol(d, x)) solve_roots(d, x, a, b, s);
                expr_free(d);
            }
    } else if (h == SYM_Piecewise && n >= 1 &&
               head_name_is(e->data.function.args[0], "List")) {
        Expr* clauses = e->data.function.args[0];
        for (size_t i = 0; i < clauses->data.function.arg_count; i++) {
            Expr* cl = clauses->data.function.args[i];
            if (head_name_is(cl, "List") && cl->data.function.arg_count == 2)
                collect_condition_boundaries(cl->data.function.args[1], x, a, b, s);
        }
    }
    /* Recurse into head + args to catch nested / product forms. */
    collect_breakpoints(e->data.function.head, x, a, b, s);
    for (size_t i = 0; i < n; i++)
        collect_breakpoints(e->data.function.args[i], x, a, b, s);
}

/* Does `e` contain a genuine branch-defining head (Abs/Sign/step/Min/Max/
 * Piecewise) with an x-dependent argument? */
static bool has_branch_head(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const char* h = e->data.function.head->type == EXPR_SYMBOL
                      ? e->data.function.head->data.symbol.name : "";
    size_t n = e->data.function.arg_count;
    if ((h == SYM_Abs || h == SYM_Sign ||
         h == SYM_HeavisideTheta || h == SYM_UnitStep ||
         h == SYM_Min || h == SYM_Max ||
         h == SYM_Piecewise) && contains_symbol(e, x))
        return true;
    if (has_branch_head(e->data.function.head, x)) return true;
    for (size_t i = 0; i < n; i++)
        if (has_branch_head(e->data.function.args[i], x)) return true;
    return false;
}

/* Does `e` contain a fractional-power radical Power[g, p] of an x-dependent
 * radicand?  Cheap structural gate for whether radical normalization is worth
 * attempting. */
static bool has_sqrt_of_x(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2 &&
        is_frac_exponent(e->data.function.args[1]) &&
        contains_symbol(e->data.function.args[0], x))
        return true;
    if (has_sqrt_of_x(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_sqrt_of_x(e->data.function.args[i], x)) return true;
    return false;
}

/* Rewrite each fractional-power radical Power[g, p] of an x-dependent radicand
 * to its Abs normal form when g is a perfect square: a real integration
 * variable makes Simplify turn Sqrt[1 - Cos[2x]] into Sqrt[2] Abs[Sin x], which
 * the Abs branch machinery below integrates exactly.  A non-square radical
 * (Sqrt[a^2 - x^2], ...) gains no Abs and is left untouched, so it falls
 * through to the ordinary cascade.  Only a result that actually introduces Abs
 * replaces the node, keeping this sound (Sqrt[g] = Abs[u] holds only for a real
 * u with u^2 = g). */
static Expr* normalize_radicals(const Expr* e, const Expr* x, const Expr* as) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return cp(e);
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2 &&
        is_frac_exponent(e->data.function.args[1]) &&
        contains_symbol(e->data.function.args[0], x)) {
        Expr* real_x = mk_fn2("Element", cp(x), expr_new_symbol("Reals"));
        Expr* real_as = as ? mk_fn2("And", real_x, cp(as)) : real_x;
        Expr* s = eval_take(mk_fn2("Simplify", cp(e), real_as));
        if (s && contains_head(s, "Abs")) return s;
        if (s) expr_free(s);
        /* fall through: keep the node, still normalize its children */
    }
    Expr* rh = normalize_radicals(e->data.function.head, x, as);
    size_t n = e->data.function.arg_count;
    Expr** ra = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++)
        ra[i] = normalize_radicals(e->data.function.args[i], x, as);
    Expr* out = expr_new_function(rh, ra, n);
    free(ra);
    return out;
}

/* Does `e` contain any branch-defining head with an x-dependent argument? */
static bool has_target(const Expr* e, const Expr* x) {
    return has_branch_head(e, x);
}

/* ---- branch resolution at an interior test point ------------------------ */

/* Compare ei vs ej at x = m under `as`: -1 (<), 0 (==), 1 (>), 2 undecidable. */
static int cmp_at(const Expr* ei, const Expr* ej, const Expr* x,
                  const Expr* m, const Expr* as) {
    Expr* im = subst(ei, x, m);
    Expr* jm = subst(ej, x, m);
    int r;
    if (refine_true(mk_fn2("Less", cp(im), cp(jm)), as)) r = -1;
    else if (refine_true(mk_fn2("Greater", cp(im), cp(jm)), as)) r = 1;
    else if (refine_true(mk_fn2("Equal", cp(im), cp(jm)), as)) r = 0;
    else r = 2;
    expr_free(im); expr_free(jm);
    return r;
}

static Expr* resolve_at(const Expr* e, const Expr* x, const Expr* m,
                        const Expr* as, bool* ok);

/* Abs/Sign/HeavisideTheta/UnitStep of `g` (already resolved; consumed). */
static Expr* resolve_unary(const char* h, Expr* g, const Expr* x,
                           const Expr* m, const Expr* as, bool* ok) {
    Expr* gm = subst(g, x, m);
    bool pos = refine_true(mk_fn2("Greater", cp(gm), expr_new_integer(0)), as);
    bool neg = pos ? false
                   : refine_true(mk_fn2("Less", cp(gm), expr_new_integer(0)), as);
    expr_free(gm);
    if (!pos && !neg) { expr_free(g); *ok = false; return NULL; }
    if (h == SYM_Abs) return pos ? g : neg_of(g);
    expr_free(g);
    if (h == SYM_Sign) return expr_new_integer(pos ? 1 : -1);
    /* HeavisideTheta / UnitStep */
    return expr_new_integer(pos ? 1 : 0);
}

/* Piecewise[clauses, default?] active branch at x = m. */
static Expr* resolve_piecewise(const Expr* e, const Expr* x, const Expr* m,
                               const Expr* as, bool* ok) {
    Expr* clauses = e->data.function.args[0];
    for (size_t i = 0; i < clauses->data.function.arg_count; i++) {
        Expr* cl = clauses->data.function.args[i];
        if (!head_name_is(cl, "List") || cl->data.function.arg_count != 2) {
            *ok = false; return NULL;
        }
        Expr* cond_m = subst(cl->data.function.args[1], x, m);
        bool t = refine_true(cp(cond_m), as);
        bool f = t ? false : refine_true(mk_fn1("Not", cp(cond_m)), as);
        expr_free(cond_m);
        if (t) return resolve_at(cl->data.function.args[0], x, m, as, ok);
        if (!f) { *ok = false; return NULL; }   /* condition undecidable here */
    }
    if (e->data.function.arg_count >= 2)
        return resolve_at(e->data.function.args[1], x, m, as, ok);
    return expr_new_integer(0);   /* Mathematica's Piecewise default */
}

/* Rewrite `e` to the single smooth branch active on the sub-interval whose
 * interior test point is `m`.  *ok is cleared when a branch cannot be proved. */
static Expr* resolve_at(const Expr* e, const Expr* x, const Expr* m,
                        const Expr* as, bool* ok) {
    if (!*ok) return NULL;
    if (!e || e->type != EXPR_FUNCTION) return cp(e);
    const char* h = e->data.function.head->type == EXPR_SYMBOL
                      ? e->data.function.head->data.symbol.name : "";
    size_t n = e->data.function.arg_count;

    if ((h == SYM_Abs || h == SYM_Sign ||
         h == SYM_HeavisideTheta || h == SYM_UnitStep) &&
        n == 1 && contains_symbol(e->data.function.args[0], x)) {
        Expr* g = resolve_at(e->data.function.args[0], x, m, as, ok);
        if (!*ok) { if (g) expr_free(g); return NULL; }
        return resolve_unary(h, g, x, m, as, ok);
    }
    if ((h == SYM_Min || h == SYM_Max) && n >= 1 &&
        contains_symbol(e, x)) {
        bool want_max = (h == SYM_Max);
        size_t best = 0;
        for (size_t j = 1; j < n; j++) {
            int c = cmp_at(e->data.function.args[j], e->data.function.args[best],
                           x, m, as);
            if (c == 2) { *ok = false; return NULL; }
            if (( want_max && c > 0) || (!want_max && c < 0)) best = j;
        }
        return resolve_at(e->data.function.args[best], x, m, as, ok);
    }
    if (h == SYM_Piecewise && n >= 1 &&
        head_name_is(e->data.function.args[0], "List")) {
        return resolve_piecewise(e, x, m, as, ok);
    }

    /* Generic node: resolve head + args and rebuild. */
    Expr* rhead = resolve_at(e->data.function.head, x, m, as, ok);
    Expr** rargs = malloc(n * sizeof(Expr*));
    size_t built = 0;
    for (size_t i = 0; i < n && *ok; i++)
        rargs[built++] = resolve_at(e->data.function.args[i], x, m, as, ok);
    if (!*ok) {
        if (rhead) expr_free(rhead);
        for (size_t i = 0; i < built; i++) if (rargs[i]) expr_free(rargs[i]);
        free(rargs);
        return NULL;
    }
    Expr* out = expr_new_function(rhead, rargs, n);
    free(rargs);
    return out;
}

/* ---- the integrator ----------------------------------------------------- */

/* Integrate[piece, {x, p, q}, Assumptions -> as] -> value, or NULL if it does
 * not close to an x-free, Integrate-free result. */
static Expr* integrate_piece(const Expr* piece, const Expr* x,
                             const Expr* p, const Expr* q, const Expr* as) {
    /* On a sub-interval provably in x <= 0, a resolved Abs branch leaves a
     * (-x)^a factor that the indefinite power rule rewrites to the principal
     * branch (-1)^a x^a = E^(I Pi a) x^a -- complex for x < 0, which corrupts a
     * real integral (e.g. Integrate[Sign[x] Abs[x]^a, {x,-1,1}] should be 0).
     * Reflect x -> -x so the radicand stays positive and the branch stays real:
     * Integrate[f, {x,p,q}] == Integrate[f /. x->-x, {x,-q,-p}].  The reflected
     * interval is in x >= 0. */
    const Expr* ip = piece; const Expr* ix = x; const Expr* ip_lo = p; const Expr* ip_hi = q;
    Expr* rpiece = NULL; Expr* rlo = NULL; Expr* rhi = NULL;
    if (refine_true(mk_fn2("LessEqual", cp(q), expr_new_integer(0)), as)) {
        rpiece = subst(piece, x, neg_of(cp(x)));
        rlo = eval_take(neg_of(cp(q)));
        rhi = eval_take(neg_of(cp(p)));
        ip = rpiece; ip_lo = rlo; ip_hi = rhi;
    }

    Expr* spec = mk_fn("List", (Expr*[]){ cp(ix), cp(ip_lo), cp(ip_hi) }, 3);
    Expr* call;
    if (as) {
        Expr* opt = mk_fn2("Rule", expr_new_symbol("Assumptions"), cp(as));
        call = mk_fn("Integrate", (Expr*[]){ cp(ip), spec, opt }, 3);
    } else {
        call = mk_fn("Integrate", (Expr*[]){ cp(ip), spec }, 2);
    }
    if (rpiece) expr_free(rpiece);
    if (rlo) expr_free(rlo);
    if (rhi) expr_free(rhi);
    Expr* v = eval_take(call);
    if (!v) return NULL;
    if (contains_head(v, "Integrate") || contains_symbol(v, x)) {
        expr_free(v); return NULL;
    }
    return v;
}

/* Integrate[f, {x, a, Infinity}] for f = E^(lam x) P(x) with P periodic of
 * period T and Re[lam] < 0:
 *     Integral_a^Infinity f = Integral_a^(a+T) f / (1 - E^(lam T)),
 * the sum over periods of a geometric series (each period is the first one
 * times E^(lam k T)).  P may carry Abs/Sign/step heads whose breakpoints form
 * an infinite periodic family (Abs[Sin[a x]]), which the piece-by-piece path
 * cannot enumerate; the single period [a, a + T] it can.  NULL unless the
 * shape is certified: lam from the exponential factors (x-linear exponents),
 * P(x + T) - P(x) == 0 by PossibleZeroQ, Re[lam] < 0 under the assumptions,
 * and the one-period integral closes. */
static void collect_exp_rate(const Expr* e, const Expr* x, Expr** acc, bool* ok) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (head_name_is(e, "Times")) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            collect_exp_rate(e->data.function.args[i], x, acc, ok);
        return;
    }
    if (head_name_is(e, "Power") && e->data.function.arg_count == 2 &&
        sym_is(e->data.function.args[0], "E") &&
        contains_symbol(e->data.function.args[1], x)) {
        Expr* c = eval_take(mk_fn2("Coefficient", cp(e->data.function.args[1]), cp((Expr*)x)));
        Expr* rest = eval_take(mk_fn1("Expand", mk_fn2("Plus", cp(e->data.function.args[1]),
                                     neg_of(mk_fn2("Times", cp(c), cp((Expr*)x))))));
        if (!c || !rest || contains_symbol(c, x) || contains_symbol(rest, x)) *ok = false;
        else *acc = eval_take(mk_fn2("Plus", *acc, cp(c)));
        if (c) expr_free(c);
        if (rest) expr_free(rest);
    }
}

static bool has_free_param(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return !(sym_is(e, "Pi") || sym_is(e, "E") || sym_is(e, "I"));
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_free_param(e->data.function.args[i])) return true;
    return false;
}
static size_t leaf_count(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return 1;
    size_t n = 1;
    for (size_t i = 0; i < e->data.function.arg_count; i++) n += leaf_count(e->data.function.args[i]);
    return n;
}

static Expr* periodic_exponential_halfline(const Expr* f, const Expr* x, const Expr* a,
                                           const Expr* b, const Expr* T_in,
                                           const Expr* assumptions) {
    if (!is_pos_infinity((Expr*)b) || is_neg_infinity((Expr*)a) || contains_symbol(a, x))
        return NULL;
    /* T > 0 (a symbolic period |2 Pi/a| is refined under the assumptions) */
    Expr* Tp = eval_take(assumptions ? mk_fn2("Refine", mk_fn1("Abs", cp(T_in)), cp((Expr*)assumptions))
                                     : mk_fn1("Abs", cp(T_in)));
    if (!Tp || !refine_true(mk_fn2("Greater", cp(Tp), expr_new_integer(0)), assumptions)) {
        if (Tp) expr_free(Tp);
        return NULL;
    }
    Expr* T = Tp;   /* owned from here: freed on every exit below */
    Expr* lam = expr_new_integer(0);
    bool ok = true;
    collect_exp_rate(f, x, &lam, &ok);
    if (!ok || (lam->type == EXPR_INTEGER && lam->data.integer == 0)) {
        expr_free(lam); expr_free(T); return NULL;
    }
    /* Re[lam] < 0 over the whole region */
    Expr* neg = mk_fn2("Less", mk_fn1("Re", cp(lam)), expr_new_integer(0));
    if (!refine_true(neg, assumptions)) { expr_free(lam); expr_free(T); return NULL; }
    Expr* P = eval_take(mk_fn2("Times", cp((Expr*)f),
                  mk_fn2("Power", mk_sym("E"), neg_of(mk_fn2("Times", cp(lam), cp((Expr*)x))))));
    Expr* xT = eval_take(mk_fn2("Plus", cp((Expr*)x), cp(T)));
    Expr* shifted = subst(P, x, xT);
    expr_free(xT);
    Expr* diff = eval_take(mk_fn1("PossibleZeroQ", mk_fn2("Plus", shifted, neg_of(cp(P)))));
    bool periodic = diff && sym_is(diff, "True");
    if (diff) expr_free(diff);
    expr_free(P);
    if (!periodic) { expr_free(lam); expr_free(T); return NULL; }
    Expr* top = eval_take(mk_fn2("Plus", cp((Expr*)a), cp(T)));
    /* the one-period integral through the full Integrate: its breakpoints on
     * the finite [a, a + T] enumerate, so the branch machinery closes it */
    Expr* spec = mk_fn3("List", cp((Expr*)x), cp((Expr*)a), top);
    Expr* call = assumptions
        ? expr_new_function(mk_sym("Integrate"), (Expr*[]){ cp((Expr*)f), spec,
              mk_fn2("Rule", mk_sym("Assumptions"), cp((Expr*)assumptions)) }, 3)
        : mk_fn2("Integrate", cp((Expr*)f), spec);
    Expr* one = eval_take(mk_fn1("Quiet", call));
    if (!one || contains_head(one, "Integrate") || contains_symbol(one, x)) {
        if (one) expr_free(one);
        expr_free(lam); expr_free(T); return NULL;
    }
    Expr* den = mk_fn2("Plus", expr_new_integer(1),
                       neg_of(mk_fn2("Power", mk_sym("E"), mk_fn2("Times", lam, cp(T)))));
    Expr* v = eval_take(mk_fn2("Times", one, mk_fn2("Power", den, expr_new_integer(-1))));
    Expr* sv = v ? eval_take(mk_fn1("Quiet", assumptions
                       ? mk_fn2("Simplify", v, cp((Expr*)assumptions))
                       : mk_fn1("Simplify", v))) : NULL;
    /* A parameter-free value carries E^((c - I) Pi)-type factors from the
     * period boundary; ComplexExpand collapses them.  Kept only if it is the
     * same number and no larger. */
    if (sv && contains_head(sv, "Complex") && !has_free_param(sv)) {
        Expr* ce = eval_take(mk_fn1("Quiet", mk_fn1("Simplify", mk_fn1("ComplexExpand", cp(sv)))));
        Expr* dz = ce ? eval_take(mk_fn1("PossibleZeroQ", mk_fn2("Plus", cp(ce), neg_of(cp(sv))))) : NULL;
        if (ce && dz && sym_is(dz, "True") && leaf_count(ce) <= leaf_count(sv)) {
            expr_free(sv); sv = ce; ce = NULL;
        }
        if (ce) expr_free(ce);
        if (dz) expr_free(dz);
    }
    expr_free(T);
    return sv;
}

Expr* integrate_piecewise_try(Expr* f, Expr* x, Expr* a, Expr* b,
                              Expr* assumptions) {
    if (!f || !x || !a || !b || x->type != EXPR_SYMBOL) return NULL;

    /* Expose any perfect-square radical as its Abs normal form (Sqrt[1-Cos[2x]]
     * -> Sqrt[2] Abs[Sin x]) so the Abs branch machinery can integrate it
     * exactly; a non-square radical is left untouched.  fn owns the rewrite. */
    Expr* fn = has_sqrt_of_x(f, x) ? normalize_radicals(f, x, assumptions) : NULL;
    const Expr* ff = fn ? fn : f;

    if (!has_target(ff, x)) { if (fn) expr_free(fn); return NULL; }

    /* 1. Collect breakpoints, keep the strictly-interior ones. */
    BPSet all; bp_init(&all);
    collect_breakpoints(ff, x, a, b, &all);
    if (all.incomplete) {
        /* The pieces are not branch-constant.  The one closed case: a periodic
         * branch structure on a half-line, damped by an exponential. */
        Expr* v = all.period_exact
                ? periodic_exponential_halfline(ff, x, a, b, all.period_exact, assumptions)
                : NULL;
        bp_free(&all);
        if (fn) expr_free(fn);
        return v;
    }
    BPSet inside; bp_init(&inside);
    for (int i = 0; i < all.n; i++) {
        if (lt_bound(a, all.pts[i], assumptions) &&
            lt_bound(all.pts[i], b, assumptions))
            bp_add(&inside, cp(all.pts[i]));
    }
    bp_free(&all);

    /* 2. Sort interior breakpoints (insertion sort, proved order). */
    for (int i = 1; i < inside.n; i++) {
        Expr* key = inside.pts[i];
        int j = i - 1;
        while (j >= 0) {
            if (refine_true(mk_fn2("Less", cp(key), cp(inside.pts[j])),
                            assumptions)) {
                inside.pts[j + 1] = inside.pts[j]; j--;
            } else if (refine_true(mk_fn2("Less", cp(inside.pts[j]), cp(key)),
                                   assumptions) ||
                       refine_true(mk_fn2("Equal", cp(inside.pts[j]), cp(key)),
                                   assumptions)) {
                break;
            } else {
                bp_free(&inside); if (fn) expr_free(fn);
                return NULL;   /* order undecidable */
            }
        }
        inside.pts[j + 1] = key;
    }

    /* 3. Endpoints [a, c1, ..., cn, b]; sum the resolved pieces. */
    int npts = inside.n + 2;
    Expr** ep = malloc((size_t)npts * sizeof(Expr*));
    ep[0] = cp(a);
    for (int i = 0; i < inside.n; i++) ep[i + 1] = cp(inside.pts[i]);
    ep[npts - 1] = cp(b);
    bp_free(&inside);

    Expr** terms = malloc((size_t)(npts - 1) * sizeof(Expr*));
    int nterms = 0;
    bool ok = true;
    for (int i = 0; i < npts - 1 && ok; i++) {
        Expr* p = ep[i];
        Expr* q = ep[i + 1];
        /* interior test point m of [p, q] (p, q may be +/-Infinity) */
        Expr* m;
        bool pinf = is_pos_infinity(p) || is_neg_infinity(p);
        bool qinf = is_pos_infinity(q) || is_neg_infinity(q);
        if (!pinf && !qinf)
            m = eval_take(mk_fn2("Divide", mk_fn2("Plus", cp(p), cp(q)),
                                 expr_new_integer(2)));           /* (p+q)/2 */
        else if (!pinf && qinf)
            m = eval_take(mk_fn2("Plus", cp(p), expr_new_integer(1)));   /* p+1 */
        else if (pinf && !qinf)
            m = eval_take(mk_fn2("Plus", cp(q), expr_new_integer(-1)));  /* q-1 */
        else
            m = expr_new_integer(0);   /* whole line: test at 0 */
        Expr* piece = resolve_at(ff, x, m, assumptions, &ok);
        expr_free(m);
        if (!ok || !piece) { if (piece) expr_free(piece); ok = false; break; }
        Expr* val = integrate_piece(piece, x, p, q, assumptions);
        expr_free(piece);
        if (!val) { ok = false; break; }
        terms[nterms++] = val;
    }
    for (int i = 0; i < npts; i++) expr_free(ep[i]);
    free(ep);

    if (!ok) {
        for (int i = 0; i < nterms; i++) expr_free(terms[i]);
        free(terms);
        if (fn) expr_free(fn);
        return NULL;
    }

    Expr* sum;
    if (nterms == 1) { sum = terms[0]; free(terms); }
    else {
        sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)nterms);
        free(terms);
    }
    if (fn) expr_free(fn);
    /* Clean up under the assumptions (a>0 etc. tidy the closed form).  Quiet: a
     * piece boundary can carry a transient 0^(1-b)-style term (an endpoint value
     * the Limit engine cannot resolve for a symbolic exponent) that raises
     * Power::infy as Simplify folds it away; that is internal, not user-facing. */
    Expr* simp;
    if (assumptions)
        simp = eval_take(mk_fn1("Quiet", mk_fn2("Simplify", sum, cp(assumptions))));
    else
        simp = eval_take(mk_fn1("Quiet", mk_fn1("Simplify", sum)));
    return simp;
}

/* ---- indefinite: one step factor with a linear argument ----------------- */

/* Is `e` an explicit real number (Integer / Rational / Real)? */
static bool is_real_literal(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT ||
                 e->type == EXPR_REAL || head_name_is(e, "Rational"));
}

/* Does `e` mention a divergent / undefined value (the antiderivative is
 * singular at the breakpoint, so no continuous antiderivative exists)? */
static bool has_singular_value(const Expr* e) {
    return contains_head(e, "DirectedInfinity") ||
           contains_symbol_name(e, "Infinity") ||
           contains_symbol_name(e, "ComplexInfinity") ||
           contains_symbol_name(e, "Indeterminate");
}

/* Is `e` a one-argument jump/impulse head of `x`: UnitStep / HeavisideTheta
 * (want_delta false) or DiracDelta (want_delta true)? */
static bool is_x_jump(const Expr* e, const Expr* x, bool want_delta) {
    bool head_ok = want_delta
        ? head_name_is(e, "DiracDelta")
        : (head_name_is(e, "UnitStep") || head_name_is(e, "HeavisideTheta"));
    return head_ok && e->data.function.arg_count == 1 && contains_symbol(e, x);
}

/* The x-free coefficient f / f[[idx]] (or 1 when f is the factor itself),
 * evaluated.  f borrowed. */
static Expr* cofactor_without(const Expr* f, bool is_times, size_t idx) {
    if (!is_times) return expr_new_integer(1);
    size_t n = f->data.function.arg_count;
    Expr** rest = (Expr**)malloc((n > 1 ? n - 1 : 1) * sizeof(Expr*));
    size_t k = 0;
    for (size_t i = 0; i < n; i++)
        if (i != idx) rest[k++] = cp(f->data.function.args[i]);
    Expr* g = eval_take(mk_fn("Times", rest, k));
    free(rest);
    return g;
}

/* Does `e` carry a jump / impulse head whose argument mentions x? */
static bool has_x_jump(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if ((head_name_is(e, "UnitStep") || head_name_is(e, "HeavisideTheta") ||
         head_name_is(e, "DiracDelta")) && contains_symbol(e, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_x_jump(e->data.function.args[i], x)) return true;
    return false;
}

/* A product with a SUM factor carrying jumps -- g(x) (UnitStep[x - a] -
 * UnitStep[x - b]), the variation-of-parameters integrand of a step-pulse
 * forcing -- has no top-level jump factor.  Expand it and integrate term by
 * term: each term then carries its jump as a factor and closes through the
 * rule below.  (The general Linearity stage does not get there for a trig
 * cofactor: it TrigReduces first and declines when that is not a sum.)  NULL
 * unless every term closes. */
static Expr* integrate_jump_sum_try(Expr* f, Expr* x, bool is_times) {
    if (!is_times) return NULL;
    bool sum_factor = false;
    for (size_t i = 0; i < f->data.function.arg_count && !sum_factor; i++) {
        const Expr* fi = f->data.function.args[i];
        if (head_name_is(fi, "Plus") && has_x_jump(fi, x)) sum_factor = true;
    }
    if (!sum_factor) return NULL;
    Expr* ex = eval_take(mk_fn1("Expand", cp(f)));
    if (!head_name_is(ex, "Plus") || expr_eq(ex, f)) { expr_free(ex); return NULL; }
    size_t n = ex->data.function.arg_count;
    Expr** parts = (Expr**)malloc(n * sizeof(Expr*));
    size_t got = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* gi = eval_take(mk_fn2("Integrate", cp(ex->data.function.args[i]), cp(x)));
        if (!gi || contains_head(gi, "Integrate")) {
            if (gi) expr_free(gi);
            for (size_t j = 0; j < got; j++) expr_free(parts[j]);
            free(parts); expr_free(ex);
            return NULL;
        }
        parts[got++] = gi;
    }
    expr_free(ex);
    Expr* sum = eval_take(mk_fn("Plus", parts, got));
    free(parts);
    return sum;
}

Expr* integrate_step_indefinite_try(Expr* f, Expr* x) {
    if (!f || !x || x->type != EXPR_SYMBOL || f->type != EXPR_FUNCTION) return NULL;

    /* Locate the jump factor J: f itself, or one factor of a Times.  An impulse
     * DiracDelta[u] takes precedence (a step cofactor is then just part of the
     * sifted value); otherwise a step UnitStep[u] / HeavisideTheta[u]. */
    bool is_times = head_name_is(f, "Times");
    const Expr* jump = NULL;
    size_t jump_idx = 0;
    bool delta = false;
    for (int pass = 0; pass < 2 && !jump; pass++) {
        bool want_delta = (pass == 0);
        if (is_times) {
            for (size_t i = 0; i < f->data.function.arg_count; i++)
                if (is_x_jump(f->data.function.args[i], x, want_delta)) {
                    jump = f->data.function.args[i]; jump_idx = i; break;
                }
        } else if (is_x_jump(f, x, want_delta)) {
            jump = f;
        }
        if (jump) delta = want_delta;
    }
    if (!jump) return integrate_jump_sum_try(f, x, is_times);

    /* u must be linear in x with a nonzero real-literal slope c: then the jump
     * sits at exactly one point, x0 = -(u - c x)/c. */
    const Expr* u = jump->data.function.args[0];
    Expr* c = eval_take(mk_fn2("D", cp(u), cp(x)));
    if (!is_real_literal(c) || contains_symbol(c, x)) { expr_free(c); return NULL; }
    Expr* czero = eval_take(mk_fn2("Equal", cp(c), expr_new_integer(0)));
    bool c_is_zero = sym_is(czero, "True");
    expr_free(czero);
    if (c_is_zero) { expr_free(c); return NULL; }
    Expr* d = eval_take(mk_fn2("Plus", cp(u),
                  neg_of(mk_fn2("Times", cp(c), cp(x)))));
    if (contains_symbol(d, x)) { expr_free(c); expr_free(d); return NULL; }
    Expr* x0 = eval_take(mk_fn2("Times", neg_of(d),
                   mk_fn2("Power", cp(c), expr_new_integer(-1))));   /* consumes d */

    Expr* g = cofactor_without(f, is_times, jump_idx);

    if (delta) {
        /* Sifting: g(x) DiracDelta[c (x - x0)] = g(x0)/|c| DiracDelta[x - x0],
         * whose antiderivative is g(x0)/|c| HeavisideTheta[x - x0].  A second
         * impulse in g is a product of distributions -- not this rule's. */
        if (contains_head(g, "DiracDelta")) {
            expr_free(g); expr_free(c); expr_free(x0); return NULL;
        }
        Expr* gx0 = subst(g, x, x0);
        expr_free(g);
        if (!gx0 || contains_symbol(gx0, x) || has_singular_value(gx0)) {
            if (gx0) expr_free(gx0);
            expr_free(c); expr_free(x0);
            return NULL;
        }
        Expr* xm = mk_fn2("Plus", cp(x), neg_of(x0));          /* consumes x0 */
        Expr* theta = mk_fn1("HeavisideTheta", xm);
        Expr* scale = mk_fn2("Power", mk_fn1("Abs", c), expr_new_integer(-1));
        Expr* args[3] = { gx0, scale, theta };
        return eval_take(mk_fn("Times", args, 3));
    }
    expr_free(c);

    /* Step.  The derivative of a step is an impulse; a DiracDelta cofactor was
     * claimed by the sifting arm above, so none is left here. */
    Expr* G = eval_take(mk_fn2("Integrate", g, cp(x)));      /* consumes g */
    if (!G || contains_head(G, "Integrate") || has_singular_value(G)) {
        if (G) expr_free(G);
        expr_free(x0);
        return NULL;
    }
    /* G(x0): must be finite, else G has a singularity at the breakpoint and the
     * glued antiderivative below would not be continuous. */
    Expr* Gx0 = subst(G, x, x0);
    expr_free(x0);
    if (!Gx0 || contains_symbol(Gx0, x) || has_singular_value(Gx0)) {
        if (Gx0) expr_free(Gx0);
        expr_free(G);
        return NULL;
    }

    /* The CONTINUOUS antiderivative (G - G(x0)) S.  For an increasing step
     * (c > 0) it is 0 left of x0 and G - G(x0) right of it; for a decreasing
     * one (c < 0) it is G - G(x0) left of x0 and 0 right of it, which is the
     * same expression.  Either way it tends to 0 at the breakpoint from both
     * sides, so it is a valid antiderivative ACROSS x0 -- the property the
     * definite integral (Newton-Leibniz) and a DSolve IC-fit rely on. */
    Expr* diff = mk_fn2("Plus", G, neg_of(Gx0));
    return eval_take(mk_fn2("Times", diff, cp(jump)));
}

/* ---- builtin ------------------------------------------------------------ */

Expr* builtin_integrate_piecewise(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;
    Expr* f    = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];
    if (!head_name_is(spec, "List") || spec->data.function.arg_count != 3)
        return NULL;
    Expr* x = spec->data.function.args[0];
    Expr* a = spec->data.function.args[1];
    Expr* b = spec->data.function.args[2];
    if (x->type != EXPR_SYMBOL) return NULL;

    Expr* assumptions = NULL;
    for (size_t t = 2; t < argc; t++) {
        Expr* opt = res->data.function.args[t];
        if (opt->type == EXPR_FUNCTION && opt->data.function.arg_count == 2 &&
            opt->data.function.head->type == EXPR_SYMBOL &&
            (opt->data.function.head->data.symbol.name == SYM_Rule ||
             opt->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
            Expr* lhs = opt->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL &&
                lhs->data.symbol.name == SYM_Assumptions) {
                assumptions = opt->data.function.args[1];
                continue;
            }
        }
        return NULL;
    }
    return integrate_piecewise_try(f, x, a, b, assumptions);
}

void integrate_piecewise_init(void) {
    symtab_add_builtin("Integrate`Piecewise", builtin_integrate_piecewise);
    symtab_get_def("Integrate`Piecewise")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Integrate`Piecewise",
        "Integrate`Piecewise[f, {x, a, b}] evaluates a definite integral whose "
        "integrand is built from Abs, Sign, HeavisideTheta, UnitStep, Min, Max, "
        "or Piecewise.  It splits [a, b] at the interior breakpoints where a "
        "branch switches (decided assumption-aware with Refine), reduces each "
        "sub-interval to its single active branch, integrates that, and sums.  "
        "Returns unevaluated when a breakpoint position or branch sign cannot be "
        "proved, or a piece does not close.");
}
