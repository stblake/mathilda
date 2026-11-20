/* risch_util.c — shared low-level helpers for the transcendental Risch modules.
 *
 * Small host-grounded eval / predicate / polynomial primitives (the `rt_` core)
 * shared across integrate_risch_transcendental.c and its tower / RDE / case
 * modules.  See risch_util.h for the public surface and the ownership notes on
 * each wrapper.
 */

#include "risch_util.h"
#include "mtime.h"          /* mth_mono_seconds (RT_STAGE) */

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "sym_intern.h"
#include "simp_trigexp_zero.h"
#include "message.h"           /* rt_solve_always_quiet: a probe */
#include "symtab.h"            /* rt_collect_params: symtab_lookup */
#include "attr.h"              /* ATTR_CONSTANT */
#include "context.h"           /* context_home */
#include "print.h"
#include <stdio.h>
#include <stdlib.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================== */
/* Small evaluation helpers.                                          */
/* ================================================================== */

/* True iff `e` is the unevaluated call `name[...]`. */
bool rt_head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == intern_symbol(name);
}

/* True iff `e` is a nonzero integer constant; writes its value to *out. */
static bool rt_is_int_const(const Expr* e, long* out) {
    if (e && e->type == EXPR_INTEGER && e->data.integer != 0) {
        *out = (long)e->data.integer; return true;
    }
    return false;
}

/* True iff `e` is a nonzero rational number constant (Integer, BigInt, or
 * Rational[p, q]) — used to group multiplicatively commensurate exponents. */
bool rt_is_rat_const(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer != 0;
    if (e->type == EXPR_BIGINT) return true;
    return e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == intern_symbol("Rational");
}

/* Numerator / denominator (den > 0) of a rational-constant Expr (Integer or
 * Rational[p, q]).  Returns false if `e` is not a small rational constant
 * (BigInt rationals are out of range here and decline).  Used to synthesize a
 * commensurability-class primitive exponent from non-integer rational ratios. */
static bool rt_rat_num_den(const Expr* e, long* num, long* den) {
    long i;
    if (rt_is_int_const(e, &i)) { *num = i; *den = 1; return true; }
    if (e && e->type == EXPR_FUNCTION && rt_head_is(e, "Rational") &&
        e->data.function.arg_count == 2) {
        const Expr* p = e->data.function.args[0];
        const Expr* q = e->data.function.args[1];
        if (p && p->type == EXPR_INTEGER && q && q->type == EXPR_INTEGER &&
            q->data.integer != 0) {
            long a = (long)p->data.integer, b = (long)q->data.integer;
            if (b < 0) { a = -a; b = -b; }
            *num = a; *den = b; return true;
        }
    }
    return false;
}

static long rt_gcd_l(long a, long b) {
    a = labs(a); b = labs(b);
    while (b) { long t = a % b; a = b; b = t; }
    return a;
}
static long rt_lcm_l(long a, long b) {
    if (a == 0 || b == 0) return 0;
    return labs(a / rt_gcd_l(a, b) * b);
}

/* Build `head[args...]` (adopting the owned `args` element pointers) and
 * evaluate it, freeing the constructed call.  Returns evaluate()'s result. */
/* Count of polynomial-algebra primitives that came back UNEVALUATED.
 *
 * The Risch decision procedure is authoritative only while its polynomial
 * arithmetic over the constant field is: a PolynomialGCD / PolynomialQuotient /
 * PolynomialRemainder / PolynomialExtendedGCD / Resultant left unevaluated (a
 * coefficient ring an engine does not cover) is then computed with as if it were
 * a polynomial, and the "no solution" the algorithm reaches from that garbage is
 * no theorem.  That is how PolynomialGCD[1, Sqrt[2] Sqrt[a]] (once unevaluated)
 * turned the elementary
 *     Cos[x] Sin[x] Cosh[a x] (-1)^(3/4) Sqrt[Pi]/(Sqrt[2] Sqrt[a])
 * into a "proved" Integrate::nonelem.  rt_decide_field snapshots this counter
 * around the field integration and demotes a NONELEMENTARY verdict reached with
 * an unevaluated primitive to UNKNOWN. */
unsigned long g_rt_prim_unevaluated = 0;

static bool rt_is_poly_primitive(const char* head) {
    return strcmp(head, "PolynomialGCD") == 0 || strcmp(head, "PolynomialQuotient") == 0
        || strcmp(head, "PolynomialRemainder") == 0
        || strcmp(head, "PolynomialExtendedGCD") == 0
        || strcmp(head, "PolynomialQuotientRemainder") == 0
        || strcmp(head, "Resultant") == 0;
}

Expr* rt_eval_call(const char* head, Expr** args, size_t n) {
    Expr* call = expr_new_function(expr_new_symbol(head), args, n);
    Expr* r = evaluate(call);
    expr_free(call);
    if (r && r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL
        && strcmp(r->data.function.head->data.symbol.name, head) == 0
        && rt_is_poly_primitive(head)) {
        g_rt_prim_unevaluated++;
        if (getenv("MATHILDA_RDE_CERT_DEBUG")) {
            char* sdbg = expr_to_string(r);
            fprintf(stderr, "[rt-prim-uneval] %s\n", sdbg); free(sdbg);
        }
    }
    return r;
}

/* Evaluate `e`, taking ownership: frees `e` and returns evaluate()'s
 * result (evaluate itself does not consume its argument). */
Expr* rt_eval_own(Expr* e) {
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* ================================================================== */
/* Generic builtin-call and predicate helpers (grounded in the host).  */
/* ================================================================== */

/* NOTE ON CORRECTNESS.  This is a decision procedure: every branch is
 * correct by construction and is NOT checked by differentiating the
 * result.  Each case fires only behind an exact structural certificate
 * (Cancel of a logarithmic derivative to a polynomial, an exactly linear
 * denominator, a genuinely x-free cofactor, ...) that already proves the
 * closed form it emits.                                                */

Expr* rt_eval1(const char* head, Expr* a) {
    return rt_eval_call(head, (Expr*[]){ a }, 1);
}
Expr* rt_eval2(const char* head, Expr* a, Expr* b) {
    return rt_eval_call(head, (Expr*[]){ a, b }, 2);
}
Expr* rt_eval3(const char* head, Expr* a, Expr* b, Expr* c) {
    return rt_eval_call(head, (Expr*[]){ a, b, c }, 3);
}

/* SolveAlways[eqn, vars] for an undetermined-coefficient ANSATZ.  Every Risch /
 * Cherry / Knowles caller uses SolveAlways as a speculative search -- an
 * inconsistent or underdetermined coefficient system just means "this ansatz does
 * not fit" and the caller declines -- so the Solve::nsdim / Solve::svars it may
 * raise is probe noise, not a user diagnostic.  Run it in a suppressed region and
 * roll the fired-count back, so the noise neither prints nor flips an enclosing
 * Check[] (Simplify[1/(E^Cos[x](c + Integrate[E^-Cos[x], x]))] leaked
 * Solve::nsdim from an exp-tower ansatz).  Arguments consumed. */
Expr* rt_solve_always_quiet(Expr* eqn, Expr* vars) {
    unsigned long fired = mth_msg_fired_count();
    mth_msg_suppress_push();
    Expr* r = rt_eval_call("SolveAlways", (Expr*[]){ eqn, vars }, 2);
    mth_msg_suppress_pop();
    mth_msg_fired_restore(fired);
    return r;
}

/* ------------------------------------------------------------------ */
/* Symbolic parameters of a SolveAlways identity (v0.6586 / v0.6588).  */
/* ------------------------------------------------------------------ */

static bool rt_name_in_list(const char* nm, const Expr* vars) {
    if (!vars || vars->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < vars->data.function.arg_count; i++) {
        const Expr* v = vars->data.function.args[i];
        if (v->type == EXPR_SYMBOL && v->data.symbol.name == nm) return true;
    }
    return false;
}

/* Collect (deduplicated interned names) the symbolic PARAMETERS of e: leaf
 * symbols (function heads are not walked) that are not in `vars`, not
 * Constant-attributed (Pi, E, EulerGamma, ...) and not kernel symbols
 * (System`/Mathilda`, e.g. Infinity).  This is the Risch layer's own reading of
 * a free symbol -- a constant of the derivation (risch_field.c: an unlisted
 * symbol has derivative 0), such as the c in E^(c x).
 *
 * `rat` is true while e sits in a rational-function position (under Plus,
 * Times, List, an integer Power's base); a parameter met anywhere else
 * (Sqrt[c], Log[c], E^c, ...) sets *nonrat.  SolveAlways reads such an atom as
 * one more unknown and can then report a solution or a contradiction that is
 * not one for the identity in c (SolveAlways[Sqrt[c] - r c == 0, {c}] ->
 * {{r -> 0}}). */
void rt_collect_params(const Expr* e, const Expr* vars, bool rat, bool* nonrat,
                       const char*** out, size_t* n, size_t* cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        if (rt_name_in_list(nm, vars)) return;
        bool seen = false;
        for (size_t i = 0; i < *n; i++) if ((*out)[i] == nm) seen = true;
        if (seen) { if (!rat) *nonrat = true; return; }
        SymbolDef* def = symtab_lookup(nm);
        if (def && (def->attributes & ATTR_CONSTANT)) return;
        /* A kernel symbol is not a parameter.  `is_system` alone cannot tell:
         * every name the kernel interns at start-up (a, b, c, x, ... among
         * them) carries it, so ask for the symbol's home context instead --
         * a bare user symbol lives in Global`. */
        if (strchr(nm, '`') ? strncmp(nm, "System`", 7) == 0
                            : strcmp(context_home(nm, def), "Global`") != 0)
            return;
        if (!rat) *nonrat = true;
        if (*n == *cap) {
            *cap = *cap ? 2 * *cap : 4;
            *out = realloc(*out, *cap * sizeof(const char*));
        }
        (*out)[(*n)++] = nm;
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    bool sum_prod = rt_head_is(e, "Plus") || rt_head_is(e, "Times") || rt_head_is(e, "List");
    bool int_pow = rt_head_is(e, "Power") && e->data.function.arg_count == 2 &&
                   (e->data.function.args[1]->type == EXPR_INTEGER ||
                    e->data.function.args[1]->type == EXPR_BIGINT);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rt_collect_params(e->data.function.args[i], vars,
                          rat && (sum_prod || (int_pow && i == 0)), nonrat, out, n, cap);
}

/* Solve the polynomial identity of an ansatz for its CONSTANT unknowns
 * (review F1/F3, 2026-11-16; unknowns made explicit 2026-11-17).
 *
 * The unknowns (residues, log coefficients, the fold's c_k) are constants of
 * the derivation, and so is every other symbol that is not a variable: a
 * symbolic parameter.  A solution may legitimately DEPEND on the parameters
 * (c = 1/a beside Log[1 + a u]).  Plain SolveAlways[eqn, {t, x}] solves for
 * every non-variable symbol, the parameters included, so it set a -> 0:
 * Integrate[(1 + a + a E^u + u)/((1 + E^u)(1 + u)), u] gave u - Log[1 + E^u].
 * Putting the parameters among the VARIABLES is wrong the other way: it
 * demands one solution independent of a, so {} no longer means "no solution"
 * (the v0.6590 false certificate).  The identity is linear in the unknowns
 * over C = Q(params), so it is solved there:
 *
 *     Solve[Thread[Flatten[CoefficientList[lhs, vars]] == 0], unknowns]
 *
 * The unknowns are NAMED by the caller (borrowed List of symbols).  v0.6589-91
 * inferred the parameters from a caller-chosen expression instead, and a
 * parameter that lived only in a tower kernel's derivative was missed:
 * Log[a + u]/u is rmR0/u in tower form, t' = 1/(a + u), and the log fold
 * "solved" a -> 0, c -> 1, so Int Log[u] was folded into t = Log[a + u]
 * (a wrong Q, caught by the answer check, but the verdict was lost).
 * Naming the unknowns leaves nothing to infer.
 *
 * Solve's answers are generic in the parameters: a solution holds outside a
 * Zariski-closed set, and {} means inconsistent for generic parameter values
 * -- the Phase 4.4 (E4.4) semantics of every parametric verdict.
 *
 * Positive evidence: lhs must be a polynomial in vars, and every rule must bind
 * one of the unknowns to an expression free of vars.  An underdetermined
 * system leaves unknowns free in the values; the callers pin them (to 0) or
 * reject them, as they did SolveAlways's.
 * Anything else returns NULL (undecided), never {}.  eqn and vars are
 * consumed; unknowns is borrowed. */
Expr* rt_solve_ansatz(Expr* eqn, Expr* vars, const Expr* unknowns) {
    Expr* sol = NULL;
    Expr* lhs = (rt_head_is(eqn, "Equal") && eqn->data.function.arg_count == 2)
        ? rt_eval2("Subtract", expr_copy(eqn->data.function.args[0]), expr_copy(eqn->data.function.args[1]))
        : NULL;
    expr_free(eqn);
    bool poly = false;
    if (lhs && vars && rt_head_is(vars, "List") && unknowns && rt_head_is(unknowns, "List")) {
        Expr* pq = rt_eval2("PolynomialQ", expr_copy(lhs), expr_copy(vars));
        poly = pq && rt_is_true(pq);
        if (pq) expr_free(pq);
        if (!poly && getenv("MATHILDA_RDE_CERT_DEBUG")) {
            char* a1 = expr_to_string(lhs); char* a2 = expr_to_string(vars);
            fprintf(stderr, "[rt-ansatz] not polynomial in %s: %.300s\n", a2, a1); free(a1); free(a2);
        }
    }
    if (poly) {
        Expr* coeffs = rt_eval1("Flatten", rt_eval2("CoefficientList", expr_copy(lhs), expr_copy(vars)));
        if (coeffs && rt_head_is(coeffs, "List")) {
            /* the unknowns that actually occur; none: the identity holds ({{}})
             * or it does not ({}) */
            bool any = false;
            for (size_t i = 0; i < unknowns->data.function.arg_count && !any; i++)
                if (!rt_free_of_x(coeffs, unknowns->data.function.args[i])) any = true;
            if (!any) {
                bool allzero = true;
                for (size_t i = 0; i < coeffs->data.function.arg_count && allzero; i++) {
                    Expr* z = rt_eval1("Together", expr_copy(coeffs->data.function.args[i]));
                    allzero = z && rt_is_zero(z);
                    if (z) expr_free(z);
                }
                sol = allzero
                    ? expr_new_function(expr_new_symbol("List"),
                          (Expr*[]){ expr_new_function(expr_new_symbol("List"), NULL, 0) }, 1)
                    : expr_new_function(expr_new_symbol("List"), NULL, 0);
                expr_free(coeffs);
            } else {
                /* Thread[coeffs == 0]; coeffs is adopted */
                Expr* thr = rt_eval1("Thread", expr_new_function(expr_new_symbol("Equal"),
                    (Expr*[]){ coeffs, expr_new_integer(0) }, 2));
                unsigned long fired = mth_msg_fired_count();
                mth_msg_suppress_push();
                sol = thr ? rt_eval_call("Solve", (Expr*[]){ thr, expr_copy((Expr*)unknowns) }, 2) : NULL;
                mth_msg_suppress_pop();
                mth_msg_fired_restore(fired);
                if (sol && rt_head_is(sol, "List")) {
                    bool ok = true;
                    for (size_t k = 0; k < sol->data.function.arg_count && ok; k++) {
                        const Expr* rl = sol->data.function.args[k];
                        if (!rt_head_is(rl, "List")) { ok = false; break; }
                        for (size_t j = 0; j < rl->data.function.arg_count && ok; j++) {
                            Expr* r = rl->data.function.args[j];
                            if (!rt_head_is(r, "Rule") || r->data.function.arg_count != 2
                                || r->data.function.args[0]->type != EXPR_SYMBOL
                                || !rt_name_in_list(r->data.function.args[0]->data.symbol.name, unknowns)) {
                                ok = false; break;
                            }
                            Expr* v = r->data.function.args[1];
                            for (size_t i = 0; i < vars->data.function.arg_count && ok; i++)
                                if (!rt_free_of_x(v, vars->data.function.args[i])) ok = false;
                        }
                    }
                    if (!ok) { expr_free(sol); sol = NULL; }
                } else if (sol) { expr_free(sol); sol = NULL; }
            }
        } else if (coeffs) expr_free(coeffs);
    }
    if (lhs) expr_free(lhs);
    if (vars) expr_free(vars);
    return sol;
}

/* Given a commensurability class of x-dependent exponents `ws[0..nw-1]` (every
 * ratio ws[j]/ws[0] a rational constant), synthesize the class *primitive*
 * exponent p = ws[0] / M (M = lcm of the ratio denominators) so that every
 * E^(ws[j]) = (E^p)^(kof[j]) with kof[j] a nonzero integer.  This generalizes
 * the earlier "some member is an integer multiple of every other" test — which
 * only handled integer ratios (E^(2u) = (E^u)^2) — to genuine rational ratios
 * (E^(x/2), E^(x/3) → primitive E^(x/6), k = 3, 2).  Returns an owned,
 * Cancel-normalized primitive Expr* and fills kof[], or NULL if any ratio is
 * not a small rational constant (genuinely independent kernels — a tower /
 * sum, handled elsewhere) or on arithmetic overflow. */
Expr* rt_class_primitive(Expr** ws, size_t nw, long* kof) {
    if (nw == 0) return NULL;
    long* nums = malloc(nw * sizeof(long));
    long* dens = malloc(nw * sizeof(long));
    long M = 1;
    bool ok = true;
    for (size_t j = 0; j < nw && ok; j++) {
        Expr* ratio = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(ws[j]), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_copy(ws[0]), expr_new_integer(-1) }, 2) }, 2));
        long a, b;
        if (ratio && rt_rat_num_den(ratio, &a, &b) && a != 0) {
            nums[j] = a; dens[j] = b; M = rt_lcm_l(M, b);
        } else ok = false;
        if (ratio) expr_free(ratio);
    }
    Expr* prim = NULL;
    if (ok && M > 0) {
        for (size_t j = 0; j < nw; j++) kof[j] = nums[j] * (M / dens[j]);
        prim = rt_eval1("Cancel", expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(ws[0]), expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_integer(M), expr_new_integer(-1) }, 2) }, 2));
    }
    free(nums); free(dens);
    return prim;
}

bool rt_is_true(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == intern_symbol("True");
}
/* FreeQ[e, x] */
bool rt_free_of_x(Expr* e, Expr* x) {
    Expr* r = rt_eval2("FreeQ", expr_copy(e), expr_copy(x));
    bool t = rt_is_true(r);
    if (r) expr_free(r);
    return t;
}
/* FreeQ[e, head] — true iff the symbol `head` occurs nowhere in e. */
bool rt_free_of_head(Expr* e, const char* head) {
    Expr* r = rt_eval2("FreeQ", expr_copy(e), expr_new_symbol(head));
    bool t = rt_is_true(r);
    if (r) expr_free(r);
    return t;
}
/* PolynomialQ[e, x] */
bool rt_is_poly(Expr* e, Expr* x) {
    Expr* r = rt_eval2("PolynomialQ", expr_copy(e), expr_copy(x));
    bool t = rt_is_true(r);
    if (r) expr_free(r);
    return t;
}
/* Degree of a polynomial in x via Length[CoefficientList[e, x]] - 1
 * (Mathilda has no Exponent builtin).  Returns -1 if CoefficientList does
 * not reduce to a List.  Call only after rt_is_poly. */
long rt_degree(Expr* e, Expr* x) {
    Expr* cl = rt_eval2("CoefficientList", expr_copy(e), expr_copy(x));
    long d = -1;
    if (cl && cl->type == EXPR_FUNCTION
        && cl->data.function.head->type == EXPR_SYMBOL
        && cl->data.function.head->data.symbol.name == intern_symbol("List"))
        d = (long)cl->data.function.arg_count - 1;
    if (cl) expr_free(cl);
    return d;
}
/* Coefficient[e, x, k] */
Expr* rt_coeff(Expr* e, Expr* x, long k) {
    return rt_eval3("Coefficient", expr_copy(e), expr_copy(x),
                    expr_new_integer(k));
}
/* Together[e] === 0 — an exact zero test for rational functions of the
 * field kernels (Together is exact and, unlike Simplify, cheap; it is only
 * ever applied here to constants and small rational cofactors). */
/* Is `e` identically zero?  A structural recursion handles the shapes that
 * dominate the Risch path — a single fraction Times[num, Power[den, -1]] or a
 * polynomial — WITHOUT combining over a common denominator: a product is zero
 * iff some factor is zero, and a reciprocal factor Power[g, negative] = 1/g is
 * never zero, so its (possibly degree-2n) denominator is never examined. Only a
 * Plus (whose terms might cancel) or an opaque head falls back to the exact
 * Together test, and then only on the smaller numerator polynomials. This keeps
 * the zero-test off the big rational denominators that made it cost ~0.5 s. */
bool rt_is_zero(Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: return e->data.integer == 0;
        case EXPR_REAL:    return e->data.real == 0.0;
        case EXPR_BIGINT:  return false;   /* normalized: |value| >= 2^63 != 0 */
        case EXPR_SYMBOL:
        case EXPR_STRING:  return false;
        default: break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL) {
        if (rt_head_is(e, "Times")) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                Expr* fac = e->data.function.args[i];
                if (rt_head_is(fac, "Power") && fac->data.function.arg_count == 2) {
                    Expr* ex = fac->data.function.args[1];
                    if (ex->type == EXPR_INTEGER && ex->data.integer < 0) continue;
                    if (rt_head_is(ex, "Rational") && ex->data.function.arg_count == 2 &&
                        ex->data.function.args[0]->type == EXPR_INTEGER &&
                        ex->data.function.args[0]->data.integer < 0) continue;
                }
                if (rt_is_zero(fac)) return true;
            }
            return false;
        }
        if (rt_head_is(e, "Power") && e->data.function.arg_count == 2) {
            Expr* ex = e->data.function.args[1];
            if (ex->type == EXPR_INTEGER && ex->data.integer < 0) return false;
            if (rt_head_is(ex, "Rational") && ex->data.function.arg_count == 2 &&
                ex->data.function.args[0]->type == EXPR_INTEGER &&
                ex->data.function.args[0]->data.integer < 0) return false;
            return rt_is_zero(e->data.function.args[0]);
        }
    }
    Expr* s = rt_eval1("Together", expr_copy(e));
    bool z = s && s->type == EXPR_INTEGER && s->data.integer == 0;
    if (s) expr_free(s);
    return z;
}

/* Parse `tmpl`, substitute the named placeholder symbols with `vals`
 * (borrowed; copied in), evaluate, and return the result (or NULL). */
Expr* rt_template(const char* tmpl, const char** names,
                         Expr** vals, size_t n) {
    Expr* t = parse_expression(tmpl);
    if (!t) return NULL;
    Expr** rules = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        rules[i] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_symbol(names[i]), expr_copy(vals[i]) }, 2);
    }
    Expr* rl = expr_new_function(expr_new_symbol("List"), rules, n);
    free(rules);
    Expr* ra = expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ t, rl }, 2);
    Expr* out = evaluate(ra);
    expr_free(ra);
    return out;
}
/* Find (borrowed) the argument u of the first Log[u] subexpression of `e`
 * whose argument depends on x. */
Expr* rt_find_log_of_x(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == intern_symbol("Log")
        && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x))
        return e->data.function.args[0];
    Expr* r = rt_find_log_of_x(e->data.function.head, x);
    if (r) return r;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        r = rt_find_log_of_x(e->data.function.args[i], x);
        if (r) return r;
    }
    return NULL;
}

/* Find (borrowed) the exponent u of the first exponential kernel of `e`
 * whose exponent depends on x, matching either Exp[u] or E^u (Power[E,u]). */
Expr* rt_find_exp_of_x(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol.name : NULL;
    if (h == intern_symbol("Exp") && e->data.function.arg_count == 1
        && !rt_free_of_x(e->data.function.args[0], x))
        return e->data.function.args[0];
    if (h == intern_symbol("Power") && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* ex = e->data.function.args[1];
        if (base->type == EXPR_SYMBOL && base->data.symbol.name == intern_symbol("E")
            && !rt_free_of_x(ex, x))
            return ex;
    }
    Expr* r = rt_find_exp_of_x(e->data.function.head, x);
    if (r) return r;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        r = rt_find_exp_of_x(e->data.function.args[i], x);
        if (r) return r;
    }
    return NULL;
}

/* Diff-back verifier shared across the transcendental Risch modules: is
 * Simplify[D[result,x] - f] identically 0?  Fast exact path first
 * (Together[TrigToExp[...]] and the rational-trig point-grid), Simplify last.
 * Correct-by-construction cases skip this; search-based ones gate on it. */
/* Is the residual `diff` (consumed) of a candidate `answer` (borrowed) zero?
 * Exactly by Together (Expand first when `expand`); when the answer carries a
 * RootSum log part (algebraic residues, review F7), which no rational
 * normalisation can zero, by PossibleZeroQ -- the verifier of record. */
bool rt_residual_zero(Expr* diff, const Expr* answer, bool expand) {
    if (!diff) return false;
    bool rs = answer && !rt_free_of_head((Expr*)answer, "RootSum");
    Expr* keep = rs ? expr_copy(diff) : NULL;
    Expr* chk = expand ? rt_eval1("Together", rt_eval1("Expand", diff))
                       : rt_eval1("Together", diff);
    bool ok = chk && rt_is_zero(chk);
    if (chk) expr_free(chk);
    if (!ok && keep) {
        Expr* pz = rt_eval1("PossibleZeroQ", keep);   /* adopts keep */
        ok = pz && rt_is_true(pz);
        if (pz) expr_free(pz);
    } else if (keep) expr_free(keep);
    return ok;
}

/* Log[a] -> Log[Together[a]] throughout e (owned in, owned out).  A rational
 * zero test treats each Log as an atom, so the same logarithm spelled two ways
 * -- Log[u + (E^u - E^-u)/(E^u + E^-u)] in TrigToExp[f] and
 * Log[(-1 + E^(2 u) + u + E^(2 u) u)/(1 + E^(2 u))] in an antiderivative --
 * never cancels.  Together is a normal form for the rational argument, so this
 * is an exact rewrite (Tanh hunt, 2026-11-19: a correct antiderivative failed
 * verification after 58 s in Simplify). */
static Expr* rt_together_log_args(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return e;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        e->data.function.args[i] = rt_together_log_args(e->data.function.args[i]);
    if (rt_head_is(e, "Log") && e->data.function.arg_count == 1) {
        Expr* a = rt_eval1("Together", expr_copy(e->data.function.args[0]));
        if (a) { expr_free(e->data.function.args[0]); e->data.function.args[0] = a; }
    }
    return e;
}

/* ---- Exponential normal form for the exact zero test ----------------------
 * The same exponential can reach a residual under several spellings that a
 * rational zero test treats as unrelated atoms: E^(E^(2u)/(1 + E^(2u))) and
 * E^(1/(1 + E^(2u))) (their exponents sum to 1), E^(a + 2u) beside E^a E^(2u),
 * E^(-2u) beside E^(2u).  The trig front end's Tanh answers carry all three
 * (Tanh hunt, 2026-11-20: PossibleZeroQ ran out of a 30 s budget on a residual
 * numerically 3e-15).  rt_expnorm rewrites, bottom-up, every x-dependent
 * E^a as E^c * prod z_i^(n_i): the exponent, its own exponentials already
 * rewritten, is expanded and split into partial fractions in each z (Apart),
 * the x-free terms collect into c, and each other term n_i m_i (n_i an integer
 * or I times one) becomes z_i^(n_i) for a fresh symbol z_i keyed on m_i.  A
 * term with any other coefficient stays an atom E^(term).  Each step is an
 * identity and z_i -> E^(m_i) maps the result back, so a rational function
 * that is 0 in the z_i proves the residual is 0; a non-zero proves nothing. */
typedef struct { Expr** key; Expr** sym; size_t n, cap; } RtExpTab;

static Expr* rt_expnorm_sym(RtExpTab* T, Expr* m) {   /* borrows m */
    for (size_t i = 0; i < T->n; i++)
        if (expr_eq(T->key[i], m)) return expr_copy(T->sym[i]);
    if (T->n == T->cap) {
        T->cap = T->cap ? 2 * T->cap : 8;
        T->key = realloc(T->key, T->cap * sizeof(Expr*));
        T->sym = realloc(T->sym, T->cap * sizeof(Expr*));
    }
    char nm[64];
    snprintf(nm, sizeof nm, "Risch`Private`rtZ%zu", T->n + 1);
    T->key[T->n] = expr_copy(m);
    T->sym[T->n] = expr_new_symbol(nm);
    return expr_copy(T->sym[T->n++]);
}

/* the integer n of a coefficient n or n I (Complex[0, n]); *imag says which */
static bool rt_expnorm_unit(const Expr* c, int64_t* n, bool* imag) {
    if (c->type == EXPR_INTEGER) { *n = c->data.integer; *imag = false; return true; }
    if (rt_head_is(c, "Complex") && c->data.function.arg_count == 2
        && c->data.function.args[0]->type == EXPR_INTEGER && c->data.function.args[0]->data.integer == 0
        && c->data.function.args[1]->type == EXPR_INTEGER) {
        *n = c->data.function.args[1]->data.integer; *imag = true; return true;
    }
    return false;
}

static bool rt_is_number_node(const Expr* e) {
    return e->type == EXPR_INTEGER || e->type == EXPR_BIGINT || e->type == EXPR_REAL
        || rt_head_is(e, "Rational") || rt_head_is(e, "Complex");
}

static Expr* rt_expnorm(Expr* e, Expr* x, RtExpTab* T);

/* the key of E^I (a process-lifetime constant, never freed) */
static const Expr* rt_expnorm_i(void) {
    static Expr* i = NULL;
    if (!i) i = expr_new_function(expr_new_symbol("Complex"),
                                  (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2);
    return i;
}

/* one term of a normalised exponent -> its factor (owned); term borrowed */
static Expr* rt_expnorm_term(Expr* term, RtExpTab* T) {
    Expr* c = NULL; Expr* m = NULL;
    int64_t n0; bool im0;
    if (rt_expnorm_unit(term, &n0, &im0) && im0) {   /* n I: (E^I)^n */
        Expr* z = rt_expnorm_sym(T, (Expr*)rt_expnorm_i());
        return expr_new_function(expr_new_symbol("Power"), (Expr*[]){ z, expr_new_integer(n0) }, 2);
    }
    /* a + b I, Gaussian integer: E^a (E^I)^b -- one constant, however the
     * evaluator merged it (E^(1 + I + x) beside E^x and E^(I + I x)) */
    if (rt_head_is(term, "Complex") && term->data.function.arg_count == 2
        && term->data.function.args[0]->type == EXPR_INTEGER
        && term->data.function.args[1]->type == EXPR_INTEGER) {
        Expr* z = rt_expnorm_sym(T, (Expr*)rt_expnorm_i());
        return expr_new_function(expr_new_symbol("Times"), (Expr*[]){
            expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_new_symbol("E"),
                expr_copy(term->data.function.args[0]) }, 2),
            expr_new_function(expr_new_symbol("Power"), (Expr*[]){ z,
                expr_copy(term->data.function.args[1]) }, 2) }, 2);
    }
    if (rt_head_is(term, "Times") && term->data.function.arg_count >= 2
        && rt_is_number_node(term->data.function.args[0])) {
        size_t k = term->data.function.arg_count;
        c = term->data.function.args[0];
        if (k == 2) m = expr_copy(term->data.function.args[1]);
        else {
            Expr** a = malloc((k - 1) * sizeof(Expr*));
            for (size_t i = 1; i < k; i++) a[i - 1] = expr_copy(term->data.function.args[i]);
            m = expr_new_function(expr_new_symbol("Times"), a, k - 1);
            free(a);
        }
    } else {
        m = expr_copy(term);
    }
    int64_t n = 1; bool imag = false;
    if (!c || rt_expnorm_unit(c, &n, &imag)) {
        if (imag) m = rt_eval_own(expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_function(expr_new_symbol("Complex"),
                            (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2), m }, 2));
        Expr* z = rt_expnorm_sym(T, m);
        expr_free(m);
        return expr_new_function(expr_new_symbol("Power"), (Expr*[]){ z, expr_new_integer(n) }, 2);
    }
    /* a Gaussian-integer coefficient (a + b I) m: E^(a m) E^(b I m), the
     * product of the two normal forms.  Left raw, E^((1 + I) x) was an atom
     * unrelated to E^x E^(I x), and a residual that is 0 did not Together to
     * 0 (gauss-elem-23: the answer check then fell to PossibleZeroQ, which
     * ran past two minutes) */
    if (rt_head_is(c, "Complex") && c->data.function.arg_count == 2
        && c->data.function.args[0]->type == EXPR_INTEGER
        && c->data.function.args[1]->type == EXPR_INTEGER) {
        Expr* mi = rt_eval_own(expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Complex"),
                (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2), expr_copy(m) }, 2));
        Expr* zr = rt_expnorm_sym(T, m);
        Expr* zi = mi ? rt_expnorm_sym(T, mi) : NULL;
        if (mi) expr_free(mi);
        int64_t a = c->data.function.args[0]->data.integer;
        int64_t b = c->data.function.args[1]->data.integer;
        expr_free(m);
        if (!zi) { expr_free(zr);
            return expr_new_function(expr_new_symbol("Power"),
                (Expr*[]){ expr_new_symbol("E"), expr_copy(term) }, 2); }
        return expr_new_function(expr_new_symbol("Times"), (Expr*[]){
            expr_new_function(expr_new_symbol("Power"), (Expr*[]){ zr, expr_new_integer(a) }, 2),
            expr_new_function(expr_new_symbol("Power"), (Expr*[]){ zi, expr_new_integer(b) }, 2) }, 2);
    }
    expr_free(m);
    return expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ expr_new_symbol("E"), expr_copy(term) }, 2);
}

/* x-free e (borrowed) with every constant E^(a + b I), b != 0 and a, b
 * integers, written E^a (E^I)^b in the table's symbol for E^I: a constant
 * coefficient E^I and the E^I inside a normalised E^(I + I x) are then the
 * same symbol (gauss-elem-23's residual).  Owned. */
static Expr* rt_expnorm_const(Expr* e, RtExpTab* T) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy(e);
    if (rt_head_is(e, "Power") && e->data.function.arg_count == 2
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == intern_symbol("E")) {
        Expr* w = e->data.function.args[1];
        if (rt_head_is(w, "Complex") && w->data.function.arg_count == 2
            && w->data.function.args[0]->type == EXPR_INTEGER
            && w->data.function.args[1]->type == EXPR_INTEGER)
            return rt_expnorm_term(w, T);
        return expr_copy(e);
    }
    size_t k = e->data.function.arg_count;
    Expr** a = malloc((k ? k : 1) * sizeof(Expr*));
    bool changed = false;
    for (size_t i = 0; i < k; i++) {
        a[i] = rt_expnorm_const(e->data.function.args[i], T);
        if (!expr_eq(a[i], e->data.function.args[i])) changed = true;
    }
    if (!changed) {
        for (size_t i = 0; i < k; i++) expr_free(a[i]);
        free(a);
        return expr_copy(e);
    }
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, k);
    free(a);
    return rt_eval_own(r);
}

/* e borrowed; returns the normalised form (owned) */
static Expr* rt_expnorm(Expr* e, Expr* x, RtExpTab* T) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy(e);
    if (rt_free_of_x(e, x)) return rt_expnorm_const(e, T);
    size_t k = e->data.function.arg_count;
    if (rt_head_is(e, "Power") && k == 2 && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == intern_symbol("E")) {
        Expr* a = rt_eval1("Expand", rt_expnorm(e->data.function.args[1], x, T));
        if (!a) return expr_copy(e);
        for (size_t i = 0; i < T->n && a; i++)
            if (!rt_free_of_x(a, T->sym[i]))
                a = rt_eval1("Expand", rt_eval2("Apart", a, expr_copy(T->sym[i])));
        if (!a) return expr_copy(e);
        size_t nt = rt_head_is(a, "Plus") ? a->data.function.arg_count : 1;
        Expr** fac = malloc((nt + 1) * sizeof(Expr*));
        size_t nf = 0;
        int64_t ecount = 0;   /* the integer part of the constant term: E^ecount */
        for (size_t i = 0; i < nt; i++) {
            Expr* t = (nt == 1 && !rt_head_is(a, "Plus")) ? a : a->data.function.args[i];
            /* an integer constant stays a power of E; every other term, x-free
             * ones included (E^(-2 I) of E^(I (1 + u))), gets its symbol: a
             * transcendental constant left as a coefficient made the Together
             * over Q(E^I)(z) grind past 60 s on a Tan[1 + u] residual */
            if (t->type == EXPR_INTEGER) { ecount += t->data.integer; continue; }
            fac[nf++] = rt_expnorm_term(t, T);
        }
        if (ecount || nf == 0) fac[nf++] = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_new_symbol("E"), expr_new_integer(ecount) }, 2);
        Expr* r = rt_eval_own(expr_new_function(expr_new_symbol("Times"), fac, nf));
        free(fac);
        expr_free(a);
        return r;
    }
    Expr** a = malloc((k ? k : 1) * sizeof(Expr*));
    for (size_t i = 0; i < k; i++) a[i] = rt_expnorm(e->data.function.args[i], x, T);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, k);
    free(a);
    return rt_eval_own(r);
}

/* Exact rational zero test of a residual d (borrowed) in x: Together of
 * TrigToExp[d] with the exponentials in the normal form of rt_expnorm and
 * every logarithm's argument in Together normal form (rt_together_log_args)
 * is 0.  One Together on the normalised form subsumes one on the raw form:
 * a rational function that is identically 0 in the raw atoms stays 0 under
 * any substitution for them.  And it is the cheaper one: the raw Together
 * ground 30 s on a Tanh residual whose normal form is 0 in 0.05 s (2026-11-20).
 * TrigToExp must come first -- Together on the raw hyperbolic form does GCDs
 * over the opaque Cosh/Sinh kernels and hangs (the Sech[x] / Sec[x]^3 hang).
 * A true is a proof; a false proves nothing. */
bool rt_exact_zero_trigexp(Expr* d, Expr* x) {
    Expr* te = rt_eval1("TrigToExp", expr_copy(d));
    if (!te) return false;
    if (rt_find_exp_of_x(te, x)) {
        RtExpTab T = { NULL, NULL, 0, 0 };
        Expr* n = rt_expnorm(te, x, &T);
        for (size_t i = 0; i < T.n; i++) { expr_free(T.key[i]); expr_free(T.sym[i]); }
        free(T.key); free(T.sym);
        expr_free(te);
        te = n;
        if (!te) return false;
    }
    if (!rt_free_of_head(te, "Log")) te = rt_eval_own(rt_together_log_args(te));
    Expr* t = te ? rt_eval1("Together", te) : NULL;   /* adopts te */
    bool z = t && rt_is_zero(t);
    if (t) expr_free(t);
    return z;
}

bool rt_verify_antideriv(Expr* result, Expr* f, Expr* x) {
    Expr* d = rt_eval2("D", expr_copy(result), expr_copy(x));
    if (!d) return false;
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ d, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(f) }, 2) }, 2);
    /* Fast EXACT zero-test before the (occasionally hang-prone) Simplify: TrigToExp
     * collapses the circular/hyperbolic kernels to exponentials, and Together over
     * that single rational function of E^(...) reduces a genuinely-zero difference
     * to 0 at once — a sound certificate.  This makes the verification of a
     * trig/hyperbolic antiderivative decidable in O(rational).  NB: TrigToExp must
     * come FIRST — a plain Together on the raw hyperbolic form itself blows up doing
     * GCD over the opaque Cosh/Sinh kernels (the Sech[x] / Sec[x]^3 hang):
     * Simplify[D[2 ArcTan[Cosh x + Sinh x]] - Sech x] and Together of it both hang,
     * yet Together[TrigToExp[...]] = 0 immediately. */
    if (rt_exact_zero_trigexp(diff, x)) { expr_free(diff); return true; }
    /* Exact single-kernel rational-trig decision (Together over the raw Laurent
     * TrigToExp form above leaves E^(-i x) poles uncombined for a Log/ArcTan-of-Tan
     * diff-back — e.g. the §5.10 hypertangent answer for a rational function of
     * Tan[x] — so it reads non-zero; and Simplify then mis-reduces the Sqrt-wrapped
     * form to a spurious residual).  trigexp_rational_is_zero clears the kernel by
     * an EXACT rational point-grid (a proof, never a numeric sample), so a genuine
     * rational-in-Tan antiderivative certifies here in sub-ms.  Rigorous both ways;
     * declines (UNKNOWN) to Simplify on multi-kernel / opaque-log forms. */
    TrigExpZeroResult tez = trigexp_rational_is_zero(diff);
    if (tez == TRIGEXP_ZERO_TRUE) { expr_free(diff); return true; }
    /* ComplexExpand+Together exact test, tried BEFORE the hang-prone Simplify.
     * The Cherry Q(i sqrt d) ei conjugate pairs differentiate to E^x/(x+a) with a in
     * Q(i sqrt d); the difference is zero but a plain Simplify over Q(i sqrt d) + the
     * E^x kernel blows up the generic multivariate GCD (hangs).  ComplexExpand splits
     * the I*Sqrt[d] algebra into real/imaginary parts over R, where Together decides
     * exactly — certifying the pair without ever reaching Simplify.  Sound:
     * ComplexExpand is an exact identity, so Together of it is 0 only for a genuine
     * zero.  On a real-radical diff (sqrt2) it returns nonzero and we fall through to
     * Simplify (which handles the reals), so this only ADDS certifications. */
    Expr* ce = rt_eval1("Together",
        rt_eval1("ComplexExpand", rt_eval1("TrigToExp", expr_copy(diff))));
    bool zce = ce && rt_is_zero(ce);
    if (ce) expr_free(ce);
    if (zce) { expr_free(diff); return true; }
    /* A RootSum answer (algebraic residues, review F7): no rational
     * normalisation sees D[RootSum[Q, z Log[S]]] - f as zero, and Simplify does
     * not either.  PossibleZeroQ is the verifier of record (rt_answer_verified). */
    if (!rt_free_of_head(result, "RootSum")) {
        Expr* pz = rt_eval1("PossibleZeroQ", diff);   /* adopts diff */
        bool zz = pz && rt_is_true(pz);
        if (pz) expr_free(pz);
        return zz;
    }
    /* PossibleZeroQ, the verifier of record (rt_answer_verified), before the
     * unbounded Simplify.  The trig front end hands back an answer in its own
     * TrigToExp spelling of the kernels: E^(E^(1+u)/(E^(-1-u) + E^(1+u))) where
     * f has E^Tanh[1+u], which equals E^(1/2) times its square root's square,
     * not the same atom.  No rational normalisation cancels that; Simplify ground
     * for 120 s, PossibleZeroQ decides in 0.6 s (Tanh hunt, 2026-11-20).  Its False
     * is not a proof of non-zero, so Simplify still runs after it. */
    {
        /* on the TrigToExp form (an exact rewrite): the answer's own spelling,
         * 0.56 s against 5 s on the raw residual */
        Expr* pz = rt_eval1("PossibleZeroQ", rt_eval1("TrigToExp", expr_copy(diff)));
        bool zz = pz && rt_is_true(pz);
        if (pz) expr_free(pz);
        if (zz) { expr_free(diff); return true; }
    }
    /* Residual (real algebraic / log) forms fall back to Simplify (original gate). */
    Expr* s = rt_eval1("Simplify", diff);
    bool z = s && s->type == EXPR_INTEGER && s->data.integer == 0;
    if (s) expr_free(s);
    return z;
}

/* True iff `e` carries an ALGEBRAIC function of x anywhere: a radical (an
 * x-dependent base raised to a non-integer Rational power — Sqrt[x] = x^(1/2),
 * x^(1/3), Sqrt[Sin[x]] = Sin[x]^(1/2)) or a Surd / Root / AlgebraicNumber of x.
 * The single-transcendental-extension procedures (and the exponential tower
 * builder) are decision procedures ONLY over a purely transcendental tower over
 * C(x); an algebraic extension (C(x, Sqrt[x])) is out of scope.  An x-FREE
 * radical (Sqrt[2], an algebraic CONSTANT coefficient) is deliberately NOT
 * flagged — those are legitimate transcendental-integrand coefficients. */
/* True iff `e` carries a TRANSCENDENTAL CONSTANT: an x-free subexpression with a
 * transcendental head (Log[2], E^2, ArcTan[1/3], Sin[1], Gamma[1/3], ...) or the
 * symbols Pi, E, EulerGamma, Catalan.  The field decision works over the
 * constant field Q(i) (algebraic constants are handled); with such a constant
 * the true constant field is larger, and a "no solution" over Q(i) is not a
 * proof: Int D[ArcTan[u] Log[2 u]] needs the coefficient (I/2)(Log[2] + Log[u])
 * and was certified NOT elementary (adversarial review, 2026-11-14).  Symbolic
 * parameters (x-free symbols) are not flagged: rational arithmetic carries them
 * exactly. */
static bool rt_transcendental_head(const char* h) {
    static const char* const heads[] = { "Log", "Exp", "Sin", "Cos", "Tan", "Cot", "Sec", "Csc",
        "ArcSin", "ArcCos", "ArcTan", "ArcCot", "ArcSec", "ArcCsc", "Sinh", "Cosh", "Tanh", "Coth",
        "Sech", "Csch", "ArcSinh", "ArcCosh", "ArcTanh", "ArcCoth", "ArcSech", "ArcCsch", "Gamma",
        "Zeta", "PolyLog", "Erf", "ProductLog", "LogGamma", NULL };
    for (int i = 0; heads[i]; i++) if (h == intern_symbol(heads[i])) return true;
    return false;
}
bool rt_has_transcendental_constant(Expr* e, Expr* x) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        const char* n = e->data.symbol.name;
        return n == intern_symbol("Pi") || n == intern_symbol("E") ||
               n == intern_symbol("EulerGamma") || n == intern_symbol("Catalan");
    }
    if (e->type != EXPR_FUNCTION) return false;
    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol.name : NULL;
    if (h && rt_transcendental_head(h) && rt_free_of_x(e, x)) return true;
    /* E^c with c free of x (E^x itself is a monomial, not a constant) */
    if (h == intern_symbol("Power") && e->data.function.arg_count == 2 && rt_free_of_x(e, x)) {
        Expr* b = e->data.function.args[0];
        if (b->type == EXPR_SYMBOL && b->data.symbol.name == intern_symbol("E")) return true;
    }
    /* E as the base of a monomial E^(f(x)) is not a constant: skip that head */
    if (h == intern_symbol("Power") && e->data.function.arg_count == 2) {
        Expr* b = e->data.function.args[0];
        if (b->type == EXPR_SYMBOL && b->data.symbol.name == intern_symbol("E"))
            return rt_has_transcendental_constant(e->data.function.args[1], x);
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rt_has_transcendental_constant(e->data.function.args[i], x)) return true;
    return false;
}

bool rt_has_algebraic_of_x(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const char* h = (e->data.function.head->type == EXPR_SYMBOL)
        ? e->data.function.head->data.symbol.name : NULL;
    if (h == intern_symbol("Power") && e->data.function.arg_count == 2) {
        Expr* b = e->data.function.args[0];
        Expr* p = e->data.function.args[1];
        /* radical: x-dependent base raised to a non-integer rational power */
        if (p && p->type == EXPR_FUNCTION && rt_head_is(p, "Rational")
            && !rt_free_of_x(b, x))
            return true;
    }
    if ((h == intern_symbol("Surd") || h == intern_symbol("Root")
         || h == intern_symbol("AlgebraicNumber")) && !rt_free_of_x(e, x))
        return true;
    if (rt_has_algebraic_of_x(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rt_has_algebraic_of_x(e->data.function.args[i], x)) return true;
    return false;
}

/* A single-extension case is valid only when its kernel's defining function `u`
 * (a Log argument or an exp exponent) is a RATIONAL function of x ALONE — i.e.
 * free of any other exp/log of x AND of any algebraic (radical) function of x.  A
 * NESTED kernel (e.g. u = E^x for the outer kernel of E^(E^x)) is a two-extension
 * tower: the single-kernel derivation Dt = u' theta would carry the unsubstituted
 * inner kernel, and SolveAlways would treat it as a free parameter and certify a
 * WRONG residue.  Such integrands must be left to the tower cases
 * (rt_log_tower_case / rt_exp_tower_case).  An ALGEBRAIC exponent (u = I Sqrt[x]
 * in E^(I Sqrt[x]), from TrigToExp of Cos[Sqrt[x]]) is likewise out of scope: the
 * derivation u' = I/(2 Sqrt[x]) makes the RDE / Hermite / SolveAlways solve run
 * over C(x, Sqrt[x]) — outside the pure-transcendental single-extension field,
 * where the coefficient arithmetic does not terminate (the Cos[Sqrt[x]] /
 * Sin[x^(1/3)] hang).  Declining here is sound: the transcendental Risch path
 * bows out and DerivativeDivides / the algebraic-substitution routes handle it. */
bool rt_kernel_simple(Expr* u, Expr* x) {
    return rt_find_exp_of_x(u, x) == NULL && rt_find_log_of_x(u, x) == NULL
        && !rt_has_algebraic_of_x(u, x);
}

/* MATHILDA_RT_PROFILE=1: per-stage wall time of the integrator cascades on
 * stderr (a stage that declines slowly is otherwise invisible).  Used through
 * RT_STAGE (risch_util.h). */
int rt_profile_on(void) {
    static int on = -1;
    if (on < 0) on = getenv("MATHILDA_RT_PROFILE") != NULL;
    return on;
}
Expr* rt_stage_timed(const char* name, Expr* (*fn)(Expr*, Expr*), Expr* f, Expr* x) {
    double t0 = mth_mono_seconds();
    Expr* r = fn(f, x);
    double dt = mth_mono_seconds() - t0;
    if (dt > 0.05) fprintf(stderr, "[rt-stage] %-32s %8.3f s  %s\n", name, dt, r ? "answer" : "declined");
    return r;
}
