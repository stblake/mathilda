/* integrate_risch_transcendental.c — recursive transcendental Risch integrator.
 *
 * The recursive Risch decision procedure for transcendental elementary
 * functions (Bronstein/Roach lineage), with arithmetic grounded in
 * Mathilda's existing Expr/poly/rat primitives (see the header for the
 * full contract).
 *
 * Structure: a differential transcendental tower over a single integration
 * variable, dispatched through the rational base case (delegated to
 * Integrate`BronsteinRational), the logarithmic / exponential / coupled /
 * tower cases, a trig-hyperbolic front-end, and (flag-gated) special-
 * function outputs.  Every branch is correct by construction behind an
 * exact structural certificate.
 */

#include "integrate_risch_transcendental.h"
#include "integrate_risch_rde.h"
#include "risch_util.h"
#include "risch_tower.h"
#include "risch_field_integrate.h"
#include "risch_singleext.h"
#include "risch_trig_frontend.h"
#include "risch_special.h"
#include "cherry_driver.h"

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "flint_bridge.h"
#include "risch_field.h"
#include "risch_hermite.h"
#include "risch_canonical.h"
#include "simp_trigexp_zero.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>





/* Dispatch the transcendental cases: the primitive (logarithmic) polynomial
 * reduction, the exponential (hyperexponential / Risch-DE) reduction, the
 * fractional (Rothstein-Trager) log-part, the real hypertangent case, and the
 * trig/hyperbolic (TrigToExp) front-end.  The general Hermite reduction for
 * repeated poles lands in a subsequent increment. */
static Expr* rt_transcendental_case(Expr* f, Expr* x) {
    Expr* r = rt_log_poly_case(f, x);
    if (r) return r;
    r = rt_exp_poly_case(f, x);
    if (r) return r;
    r = rt_frac_case(f, x);
    if (r) return r;
    /* Fast rational-of-a-single-exp path (kernelize -> pure rational integral in
     * t, which the FLINT-accelerated rational integrator closes) BEFORE the
     * repeated-pole Hermite ansatz: when the integrand is a rational function of
     * one exponential with a linear exponent (F free of x, u' free of x), this
     * closes it in O(deg) instead of the O(mult)-variable SolveAlways Hermite
     * solve.  It declines (F still carries x, nonlinear exponent, multi-kernel)
     * exactly when the genuine coupled/tower cases below are needed, so their
     * cleaner forms still win; rt_frac_case above keeps precedence for the
     * squarefree ArcTan/Log parts.  Correctness rests on the same diff-back gate. */
    r = rt_exp_ratreduce_case(f, x);
    if (r) return r;
    r = rt_hermite_case(f, x);
    if (r) return r;
    r = rt_hyperexp_case(f, x);
    if (r) return r;
    r = rt_expsum_case(f, x);   /* direct multi-kernel exponential sums */
    if (r) return r;
    /* The flat-tower cases run before the recursive one because they close a
     * PRIMITIVE-POLYNOMIAL class the recursive path still declines — e.g.
     * Log[Log[x]]/(x Log[x]) -> Log[Log[x]]^2/2 needs the new-logarithm fold-back
     * (§5.8 IntegratePrimitivePolynomial / LimitedIntegrate) that the deterministic
     * recursion does not yet do; their SolveAlways ansatz finds it directly.  So
     * they are NOT redundant with rt_recursive_tower_case and cannot be retired
     * (Gap 4) until the recursive primitive-polynomial path gains LimitedIntegrate. */
    /* Gap 4: the genuine one-extension recursion is now the PRIMARY tower path — with
     * the §5.8 LimitedIntegrate new-logarithm fold-back (rt_limited_field_integrate)
     * it subsumes the primitive-polynomial class the flat SolveAlways ansätze were
     * kept for (Log[Log[x]]/(x Log[x]) -> Log[Log[x]]^2/2, ...).  The flat cases
     * (rt_log_tower_case / rt_exp_tower_case) remain only as a FALLBACK for anything
     * the deterministic recursion declines (they also carry the rt_tower_solve Cherry
     * P5 substrate). */
    r = rt_recursive_tower_case(f, x);   /* one-extension recursion (mixed / rational coeff) */
    if (r) return r;
    r = rt_log_tower_case(f, x);   /* flat ansatz fallback (nested log tower) */
    if (r) return r;
    r = rt_exp_tower_case(f, x);    /* flat ansatz fallback (nested exp tower) */
    if (r) return r;
    /* Real tangent monomial (Bronstein §5.10): integrate rational functions of a
     * single Tan[u] directly and real, retiring the TrigToExp route below for
     * real-tangent integrands (which it strands at an I-laden complex-log form). */
    r = rt_hypertangent_case(f, x);
    if (r) return r;
    r = rt_trig_frontend(f, x);
    if (r) return r;
    /* Dependent-logarithm fallback.  When every case above declined, the
     * integrand may carry Q-linearly DEPENDENT logarithmic generators that no
     * single case models — e.g. Log[x/(1+x)] + Log[1+x], two distinct Log
     * kernels that collapse to the single generator Log[x].  Normalize with the
     * log-of-product / log-of-power expansion (Log[a b] -> Log a + Log b,
     * Log[b^p] -> p Log b); the evaluator then recombines the dependent pair
     * (Log[x/(1+x)] + Log[1+x] -> Log[x]), and rt_log_poly_case integrates the
     * collapsed Log[x]/x -> Log[x]^2/2.  Only runs after the original form failed
     * (so it never changes an already-handled case's cleaner output), value-
     * preserving (rt_expand_logs is an identity on the function's value), and
     * guarded against re-entry so the retry itself cannot loop. */
    static int rt_in_logcombine = 0;
    if (!rt_in_logcombine) {
        Expr* fn = rt_eval_own(rt_expand_logs(f));
        if (fn && !expr_eq(fn, f)) {
            rt_in_logcombine = 1;
            r = rt_transcendental_case(fn, x);
            rt_in_logcombine = 0;
            expr_free(fn);
            if (r) return r;
        } else if (fn) {
            expr_free(fn);
        }
    }
    return NULL;
}

/* ================================================================== */
/* Constant-base exponential support: a^u -> E^(u Log a).             */
/* ================================================================== */
/* The evaluator canonicalizes E^(c Log[a]) back to a^c (power.c make_power),
 * so a base-a exponential (a != e) is stored as Power[a, u] and NONE of the
 * exponential recognizers above (all keyed to Exp / Power[E, .]) fire — the
 * integrand strands unintegrated even though a^u = E^(u Log a) is an ordinary
 * hyperexponential monomial.  We bridge that by DEBASING: rewrite every
 * x-dependent constant-base power a^u to E^(u * (sum_j e_j L_{p_j})), where
 * a = prod_j p_j^{e_j} is the prime factorisation and each L_{p_j} is a fresh
 * OPAQUE constant symbol standing for Log[p_j].  Because the exponent carries
 * L_{p_j} (not the literal Log[p_j]), make_power does NOT collapse it, so the
 * base-e machinery sees a genuine E^(...) tower with a constant log-derivative.
 * Sharing L by prime makes commensurate bases commensurate (4^x = E^(2 x L_2),
 * 2^x = E^(x L_2) -> one kernel), which the naive per-base symbol would miss.
 * After integration the map L_{p_j} -> Log[p_j] is substituted back (and the
 * evaluator re-collapses E^(x Log a) -> a^x), so the result is in the original
 * base.  The rewrite is an exact identity, so correctness is preserved. */

/* When false, the a^u debasing pre-pass is suppressed (see the header): the
 * Integrate inexact path disables it so a rationalised float base is not given a
 * spuriously exact Log-based closed form. */
static bool rt_debase_enabled = true;
bool rt_transcendental_set_debase(bool enabled) {
    bool prev = rt_debase_enabled;
    rt_debase_enabled = enabled;
    return prev;
}

/* Opaque interned symbol name L_p standing for Log[p] (p a prime literal). */
static const char* rt_baselog_name(Expr* p) {
    char* ps = expr_to_string(p);
    char buf[96];
    snprintf(buf, sizeof buf, "Integrate`bLog$%s", ps ? ps : "0");
    if (ps) free(ps);
    return intern_symbol(buf);
}

/* a is a positive rational constant (Integer >= 2, positive Bigint, or a
 * positive Rational) other than 1 — a usable exponential base. */
static bool rt_is_pos_rational_base(Expr* a) {
    if (!a) return false;
    if (a->type == EXPR_INTEGER) return a->data.integer >= 2;
    if (a->type == EXPR_BIGINT)  return mpz_cmp_ui(a->data.bigint, 1) > 0;
    if (rt_head_is(a, "Rational") && a->data.function.arg_count == 2) {
        Expr* num = a->data.function.args[0];
        return (num->type == EXPR_INTEGER && num->data.integer > 0)
            || (num->type == EXPR_BIGINT && mpz_sgn(num->data.bigint) > 0);
    }
    return false;
}

/* Collect (owned, deduplicated) every distinct constant base a of an
 * x-dependent power a^u in `e`. */
static void rt_collect_rational_bases(Expr* e, Expr* x,
                                      Expr*** arr, size_t* n, size_t* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (rt_head_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* a = e->data.function.args[0];
        Expr* u = e->data.function.args[1];
        if (rt_is_pos_rational_base(a) && rt_free_of_x(a, x) && !rt_free_of_x(u, x)) {
            bool dup = false;
            for (size_t i = 0; i < *n; i++) if (expr_eq((*arr)[i], a)) { dup = true; break; }
            if (!dup) {
                if (*n == *cap) { *cap = *cap ? *cap * 2 : 4;
                                  *arr = realloc(*arr, *cap * sizeof(Expr*)); }
                (*arr)[(*n)++] = expr_copy(a);
            }
        }
    }
    rt_collect_rational_bases(e->data.function.head, x, arr, n, cap);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        rt_collect_rational_bases(e->data.function.args[i], x, arr, n, cap);
}

/* Log a = sum_j e_j L_{p_j} for a = prod_j p_j^{e_j}, via FactorInteger (which
 * returns signed exponents for a rational a).  Records each prime p_j (owned,
 * deduplicated) into *primes for the back-substitution map.  Owned coeff, or
 * NULL if a factors to a unit. */
static Expr* rt_base_logcoeff(Expr* a, Expr*** primes, size_t* np, size_t* pcap) {
    Expr* fi = rt_eval1("FactorInteger", expr_copy(a));
    if (!fi || !rt_head_is(fi, "List")) { if (fi) expr_free(fi); return NULL; }
    Expr** terms = NULL; size_t nt = 0, tcap = 0;
    for (size_t i = 0; i < fi->data.function.arg_count; i++) {
        Expr* pe = fi->data.function.args[i];
        if (!rt_head_is(pe, "List") || pe->data.function.arg_count != 2) continue;
        Expr* p = pe->data.function.args[0];
        Expr* ex = pe->data.function.args[1];
        if (p->type == EXPR_INTEGER && p->data.integer <= 1) continue;
        bool dup = false;
        for (size_t k = 0; k < *np; k++) if (expr_eq((*primes)[k], p)) { dup = true; break; }
        if (!dup) {
            if (*np == *pcap) { *pcap = *pcap ? *pcap * 2 : 4;
                                *primes = realloc(*primes, *pcap * sizeof(Expr*)); }
            (*primes)[(*np)++] = expr_copy(p);
        }
        Expr* term = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(ex), expr_new_symbol(rt_baselog_name(p)) }, 2);
        if (nt == tcap) { tcap = tcap ? tcap * 2 : 4; terms = realloc(terms, tcap * sizeof(Expr*)); }
        terms[nt++] = term;
    }
    expr_free(fi);
    if (nt == 0) { free(terms); return NULL; }
    Expr* coeff = rt_eval_own(expr_new_function(expr_new_symbol("Plus"), terms, nt));
    free(terms);
    return coeff;
}

/* Debase every x-dependent constant-base power of `f` (see block comment).
 * Returns the rewritten integrand and sets *backsub_out to the owned rule list
 * {L_{p_j} -> Log[p_j]} — or NULL (leaving *backsub_out NULL) when `f` carries
 * no such base, so the caller integrates it directly. */
static Expr* rt_debase_exponentials(Expr* f, Expr* x, Expr** backsub_out) {
    *backsub_out = NULL;
    if (!rt_debase_enabled) return NULL;
    Expr** bases = NULL; size_t nb = 0, bcap = 0;
    rt_collect_rational_bases(f, x, &bases, &nb, &bcap);
    if (nb == 0) { free(bases); return NULL; }

    Expr** primes = NULL; size_t np = 0, pcap = 0;
    Expr** rules = malloc(nb * sizeof(Expr*)); size_t nr = 0;
    for (size_t i = 0; i < nb; i++) {
        Expr* coeff = rt_base_logcoeff(bases[i], &primes, &np, &pcap);
        if (!coeff) continue;
        Expr* upat = expr_new_function(expr_new_symbol("Pattern"),
            (Expr*[]){ expr_new_symbol("Integrate`bU"),
                       expr_new_function(expr_new_symbol("Blank"), NULL, 0) }, 2);
        Expr* lhs = expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ expr_copy(bases[i]), upat }, 2);
        Expr* rhs = expr_new_function(expr_new_symbol("Exp"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ expr_new_symbol("Integrate`bU"), coeff }, 2) }, 1);
        rules[nr++] = expr_new_function(expr_new_symbol("RuleDelayed"),
            (Expr*[]){ lhs, rhs }, 2);
    }
    for (size_t i = 0; i < nb; i++) expr_free(bases[i]);
    free(bases);
    if (nr == 0) {
        free(rules);
        for (size_t i = 0; i < np; i++) expr_free(primes[i]);
        free(primes);
        return NULL;
    }
    Expr* rl = expr_new_function(expr_new_symbol("List"), rules, nr);
    free(rules);
    Expr* fd = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceRepeated"),
        (Expr*[]){ expr_copy(f), rl }, 2));

    Expr** brules = malloc((np ? np : 1) * sizeof(Expr*));
    for (size_t i = 0; i < np; i++)
        brules[i] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_new_symbol(rt_baselog_name(primes[i])),
                       expr_new_function(expr_new_symbol("Log"),
                           (Expr*[]){ expr_copy(primes[i]) }, 1) }, 2);
    *backsub_out = expr_new_function(expr_new_symbol("List"), brules, np);
    free(brules);
    for (size_t i = 0; i < np; i++) expr_free(primes[i]);
    free(primes);
    return fd;
}

/* ================================================================== */
/* Top-level integration dispatch.                                    */
/* ================================================================== */

/* Returns a fresh antiderivative (also self-verified by the recognizers,
 * and re-verified by the caller's diff-back gate) or NULL if no case
 * applies.  Dispatch order of the recursive Risch algorithm:
 *   1. rational base case (delegated to the recursive rational Risch,
 *      Integrate`BronsteinRational);
 *   2. transcendental case over a single logarithmic / exponential
 *      monomial extension (rt_transcendental_case — the recursive Risch
 *      proper: Hermite reduction, residue log-part, and the polynomial /
 *      Risch-differential-equation reductions);
 *   3. special-function outputs (Erf / dilog / Ei / li forms).
 * Every branch is verified by differentiation, so a mis-reduction can
 * only decline, never emit a wrong closed form.  NB: this must NOT fall
 * back on the parallel-Risch (pmint) engine Integrate`RischNorman —
 * that is a different algorithm; RischTranscendental is the recursive Risch. */
static Expr* rt_integrate_core(Expr* f, Expr* x) {
    Expr* r = rt_rational_case(f, x);
    if (r) return r;
    r = rt_transcendental_case(f, x);
    if (r) return r;
    r = extended_liouville_solve(f, x, RT_SF_TOP_ANY);   /* Cherry dispatch (C0 seam) */
    if (r) return r;
    return NULL;
}

/* Wrapper: debase any constant-base exponential a^u (a != e) to the base-e form
 * the recursive Risch machinery handles, integrate, then map the opaque Log
 * symbols back so the antiderivative is expressed in the original base. */
static Expr* rt_integrate(Expr* f, Expr* x) {
    Expr* backsub = NULL;
    Expr* fd = rt_debase_exponentials(f, x, &backsub);
    if (!fd) return rt_integrate_core(f, x);      /* no a^u base present */
    Expr* r0 = rt_integrate_core(fd, x);
    expr_free(fd);
    if (!r0) { expr_free(backsub); return NULL; }
    Expr* r = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ r0, backsub }, 2));             /* adopts r0, backsub; L_p -> Log[p] */
    /* Fuse the split prime-logs a COMPOSITE base leaves behind.  A single-prime
     * base collapses cleanly (make_power: E^(k Log p) -> p^k), but 6 = 2·3 debases
     * to Log2 + Log3 sums (in both exponents and the 1/Log(a) constant), e.g.
     * Log[1 + E^((Log2+Log3) x)]/(Log2+Log3).  n Log a + m Log b = Log[a^n b^m] is
     * exact for positive integers a, b, so recombining pairs restores Log[6]
     * (whereupon E^(x Log6) re-collapses to 6^x).  It fires only when two
     * integer-logs actually pair, so it is a no-op on the clean single-base forms
     * (a lone Log[3] denominator is untouched). */
    Expr* logfuse = parse_expression(
        "{ (nF_. Log[aF_Integer] + mF_. Log[bF_Integer] + rF_. /; aF > 0 && bF > 0)"
        "    :> Log[aF^nF bF^mF] + rF }");
    if (r && logfuse) {
        Expr* r2 = rt_eval_own(expr_new_function(expr_new_symbol("ReplaceRepeated"),
            (Expr*[]){ r, logfuse }, 2));            /* adopts r, logfuse */
        if (r2) r = r2;
    } else if (logfuse) {
        expr_free(logfuse);
    }
    return r;
}

/* ================================================================== */
/* Elementary-integrability decision procedure (P3).                  */
/* ================================================================== */

/* True iff the head names a function that is NOT elementary — the special
 * functions the integrator emits for a non-elementary antiderivative — or an
 * unevaluated Integrate (a partial / declined result).  RootSum and Root are
 * elementary (finite sums of logs / algebraic numbers) and are NOT listed. */
static bool rt_head_nonelementary(const char* h) {
    if (!h) return false;
    static const char* NE[] = {
        "ExpIntegralEi", "ExpIntegralE", "LogIntegral",
        "SinIntegral", "CosIntegral", "SinhIntegral", "CoshIntegral",
        "Erf", "Erfc", "Erfi", "FresnelS", "FresnelC",
        "PolyLog", "Integrate", NULL };
    for (int i = 0; NE[i]; i++)
        if (h == intern_symbol(NE[i])) return true;
    return false;
}

/* An antiderivative expressed WITHOUT any non-elementary special function (and
 * without an unevaluated Integrate) is an elementary closed form: exhibiting it
 * (correct by construction from the integrator) proves f is elementary-integrable. */
static bool rt_expr_is_elementary(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    if (e->data.function.head->type == EXPR_SYMBOL
        && rt_head_nonelementary(e->data.function.head->data.symbol.name))
        return false;
    if (!rt_expr_is_elementary(e->data.function.head)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!rt_expr_is_elementary(e->data.function.args[i])) return false;
    return true;
}

/* rt_has_algebraic_of_x (radical / Surd / Root of x anywhere) is shared from
 * risch_util.c — the single-extension kernel gate (rt_kernel_simple) and the
 * trig front-end scope guard use the same predicate.  See risch_util.h. */

/* Field-path elementary-integrability decision.  Routes f through the AUTHORITATIVE
 * recursive field integrator (rt_field_integrate) in decision mode — the same
 * tower build / substitution / gate as rt_recursive_tower_case, but reading the
 * VERDICT rather than assembling an answer.  Returns:
 *   RT_DEC_ELEMENTARY     — a full antiderivative (no unintegrated remainder);
 *   RT_DEC_NONELEMENTARY  — an authoritative certificate fired: a non-constant
 *                           residue (Thm 5.6.1(ii)), a Risch DE with no rational
 *                           solution (Ch.6, rde_base / exact-bound ansatz), or the
 *                           §5.8 Dc!=0 primitive certificate;
 *   RT_DEC_UNKNOWN        — out of the single-tower field scope (algebraic
 *                           extensions, tower-builder rejects, or a genuine gap).
 * Sound by construction: an ELEMENTARY / NONELEMENTARY verdict is emitted only
 * behind an exact certificate — otherwise UNKNOWN. */
static RtDecision rt_decide_field(Expr* f, Expr* x) {
    /* Scope guard: an algebraic function of x (radical / Surd / Root) puts the
     * integrand outside the purely-transcendental tower the field decision is a
     * decision procedure over — decline to UNKNOWN rather than emit an unsound
     * certificate.  Checked on the ORIGINAL f (Sqrt[x] survives inside Cos[Sqrt[x]]
     * and E^Sqrt[x] alike), before exponentialization. */
    if (rt_has_algebraic_of_x(f, x)) return RT_DEC_UNKNOWN;
    /* Mirror the integrator's trig/hyperbolic front-end (rt_trig_frontend) so the
     * DECISION routes trig integrands through the SAME Gaussian exponential tower
     * the True side uses.  Without this the tower builder cannot form a Risch tower
     * from a bare Sin/Cos/Tan kernel and the field decision falls to UNKNOWN — so
     * genuinely non-elementary trig integrands (Sin[x]/x -> SinIntegral, Sin[x^2] ->
     * FresnelS, Cos[x Log x], Cos[E^x]) returned `undec` instead of the authoritative
     * `False`.  TrigToExp exponentializes the trig kernels; rt_powers_to_exp then
     * re-exposes the E^(c Log b) kernels that TrigToExp leaves collapsed as general
     * powers b^e (Cos[x Log x] -> x^(±I x)).  Both rewrites are exact, so the tower —
     * and hence the residue / Risch-DE certificate read off it — is authoritative. */
    Expr* fe = rt_eval1("TrigToExp", expr_copy(f));
    Expr* src;
    if (fe && !expr_eq(fe, f)) {
        src = fe;                       /* trig/hyperbolic present: use exp form */
    } else {
        if (fe) expr_free(fe);
        src = expr_copy(f);             /* no trig: decide f directly */
    }
    Expr* pw = rt_powers_to_exp(src, x);   /* re-expose b^e exponential kernels */
    expr_free(src);
    Expr* fx = rt_expand_exp_sums(pw);
    expr_free(pw);
    RtTower T;
    /* min_n = 1: the decision routes SINGLE-kernel integrands (E^x/x, E^(x^2),
     * 1/Log[x]) through the same authoritative field integrator too — the flat-tower
     * integrate cases that normally serve n=1 do not carry the decision certificates. */
    if (!rt_tower_build_min(fx, x, &T, 1)) { rt_tower_free(&T); expr_free(fx); return RT_DEC_UNKNOWN; }
    Expr* F = rt_eval1("Together", rt_subst_kernels(fx, &T));
    RtDecision verdict = RT_DEC_UNKNOWN;
    if (F && rt_find_exp_of_x(F, x) == NULL && rt_find_log_of_x(F, x) == NULL) {
        Expr* num = rt_eval1("Numerator", expr_copy(F));
        Expr* den = rt_eval1("Denominator", expr_copy(F));
        Expr** vv = malloc((T.n + 1) * sizeof(Expr*));
        vv[0] = expr_copy(x);
        for (size_t i = 0; i < T.n; i++) vv[i + 1] = expr_copy(T.t[i]);
        Expr* vlist = expr_new_function(expr_new_symbol("List"), vv, T.n + 1);
        free(vv);
        Expr* pqn = num ? rt_eval2("PolynomialQ", expr_copy(num), expr_copy(vlist)) : NULL;
        Expr* pqd = den ? rt_eval2("PolynomialQ", expr_copy(den), expr_copy(vlist)) : NULL;
        bool gate = num && den && rt_is_true(pqn) && rt_is_true(pqd);
        if (pqn) expr_free(pqn);
        if (pqd) expr_free(pqd);
        expr_free(vlist);
        if (num) expr_free(num);
        if (den) expr_free(den);
        if (gate) {
            bool       save_mode = g_rt_decide_mode;
            RtDecision save_dec  = g_rt_decision;
            g_rt_decide_mode = true;
            g_rt_decision    = RT_DEC_UNKNOWN;
            Expr* rem = NULL;
            Expr* Q = rt_field_integrate(F, &T, (long)T.n - 1, x, &rem);
            if (Q && !rem)                               verdict = RT_DEC_ELEMENTARY;
            else if (Q && rem)                           verdict = RT_DEC_NONELEMENTARY;
            else if (g_rt_decision == RT_DEC_NONELEMENTARY) verdict = RT_DEC_NONELEMENTARY;
            else                                         verdict = RT_DEC_UNKNOWN;
            if (Q) expr_free(Q);
            if (rem) expr_free(rem);
            g_rt_decide_mode = save_mode;
            g_rt_decision    = save_dec;
        }
    }
    if (F) expr_free(F);
    rt_tower_free(&T);
    expr_free(fx);
    return verdict;
}

/* Complete elementary-integrability verdict for f w.r.t. x.  True side: exhibit an
 * elementary antiderivative via the full integrator (an existence proof — correct
 * by construction).  False side: the field decision's authoritative non-elementary
 * certificate.  Otherwise UNKNOWN.  Checking True first makes an exhibited
 * antiderivative dominate any certificate, so the verdict is robust. */
static RtDecision rt_decide(Expr* f, Expr* x) {
    Expr* anti = rt_integrate(f, x);
    if (anti) {
        bool elem = rt_expr_is_elementary(anti);
        expr_free(anti);
        if (elem) return RT_DEC_ELEMENTARY;   /* elementary closed form exhibited */
    }
    RtDecision fld = rt_decide_field(f, x);
    if (fld == RT_DEC_ELEMENTARY)    return RT_DEC_ELEMENTARY;
    if (fld == RT_DEC_NONELEMENTARY) return RT_DEC_NONELEMENTARY;
    return RT_DEC_UNKNOWN;
}

/* ================================================================== */
/* Public builtin.                                                    */
/* ================================================================== */

Expr* builtin_rischtranscendental(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 2) return NULL;

    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];

    /* The integration variable must be a single symbol. */
    if (x->type != EXPR_SYMBOL) return NULL;

    /* Scope gate: the recursive transcendental Risch algorithm is a decision
     * procedure ONLY over a purely transcendental tower over C(x).  An algebraic
     * function of x — a radical (Sqrt[1+Sin[x]] = (1+Sin[x])^(1/2)), Surd, or
     * Root of x — puts the integrand in an algebraic extension it does not
     * handle, so bail immediately rather than churn through rt_integrate (which
     * would decline anyway, after needless work).  An x-free algebraic constant
     * (Sqrt[2]) is a legitimate transcendental coefficient and is NOT flagged. */
    if (rt_has_algebraic_of_x(f, x)) return NULL;

    /* Correct by construction: rt_integrate returns a result only behind an
     * exact certificate, so no differentiation check is applied (a Risch
     * integrator is a decision procedure, not a guess-and-verify search).
     *
     * Mute arithmetic warnings across the decision procedure: its internal
     * algebraic manipulations (kernel substitutions, degree/coefficient probes,
     * Together/Cancel on tower normal forms) legitimately form transient
     * singular expressions — e.g. E^(x/E^(1/x))/E^(2/x) provokes a 1/0 while a
     * candidate exponent is reduced — none of which are part of the returned
     * antiderivative.  Surfacing Power::infy from them would be spurious (as it
     * is in Mathematica's Integrate). */
    arith_warnings_mute_push();
    Expr* result = rt_integrate(f, x);
    /* Decision half (P3): when the recursive algorithm produced no elementary
     * antiderivative, consult the authoritative field decision.  If it PROVES the
     * integrand has no elementary integral (a non-constant residue, a Risch DE with
     * no rational solution, or the Dc!=0 primitive certificate — never a bounded-
     * ansatz decline), report it Mathematica-style.  A stray certificate from an
     * abandoned attempt cannot fire here: rt_decide_field re-derives the verdict
     * with elementary-result precedence. */
    if (!result && rt_decide_field(f, x) == RT_DEC_NONELEMENTARY) {
        char* fs = expr_to_string(f);
        char* xs = expr_to_string(x);
        fprintf(stderr,
            "Integrate::nonelem: The integrand %s has no antiderivative "
            "elementary in %s.\n", fs ? fs : "?", xs ? xs : "?");
        free(fs); free(xs);
    }
    arith_warnings_mute_pop();
    return result;
}

/* Risch`ElementaryIntegralQ[f, x] — the Bronstein elementary-integrability decision
 * predicate.  Returns True when f has an antiderivative expressible in elementary
 * terms (exhibited by the integrator — an existence proof), False when the recursive
 * Risch decision procedure PROVES none exists (§5.6 residue criterion Thm 5.6.1(ii),
 * or a Ch.6 Risch DE with no rational solution, or the §5.8 Dc!=0 certificate), and
 * stays UNEVALUATED (with an ElementaryIntegralQ::undec message) when the verdict is
 * outside the single-tower field scope (algebraic extensions, deeper structures).
 * Sound by construction: a Boolean is returned only behind an exact certificate. */
static Expr* builtin_elementaryintegralq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    arith_warnings_mute_push();
    RtDecision d = rt_decide(f, x);
    arith_warnings_mute_pop();
    if (d == RT_DEC_ELEMENTARY)    return expr_new_symbol("True");
    if (d == RT_DEC_NONELEMENTARY) return expr_new_symbol("False");
    /* Undecided within scope: emit a diagnostic and leave the call unevaluated
     * (the *Q-returns-Boolean convention allows NULL only with an emitted tag). */
    {
        char* fs = expr_to_string(f);
        char* xs = expr_to_string(x);
        fprintf(stderr,
            "Risch`ElementaryIntegralQ::undec: Cannot decide elementary "
            "integrability of %s with respect to %s.\n", fs ? fs : "?", xs ? xs : "?");
        free(fs); free(xs);
    }
    return NULL;
}

/* Risch`RischDE[f, g, x] — solve the base-field Risch differential equation
 * D[y] + f y = g for y in C(x) (Bronstein Ch.6), exposing the internal rde_base
 * driver so the solver (and the coupled-system layer built on it) is directly
 * unit-testable.  Returns y (the unique rational solution, or 0 when g == 0), or
 * leaves the call unevaluated when no rational solution exists ("no solution" —
 * the term is non-elementary in C(x)).  Coefficients may lie in C(i)(x): the
 * Gaussian case falls through the FLINT fast path to the Expr pipeline. */
static Expr* builtin_risch_rischde(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* g = res->data.function.args[1];
    Expr* x = res->data.function.args[2];
    if (x->type != EXPR_SYMBOL) return NULL;
    /* If f or g carries a Log/Exp kernel of x, solve over the transcendental tower
     * (Gap 1, rde_tower); otherwise the base field C(x) (rde_base).  NULL bubbles
     * back unevaluated (no rational solution). */
    bool has_kernel = rt_find_exp_of_x(f, x) || rt_find_log_of_x(f, x)
                   || rt_find_exp_of_x(g, x) || rt_find_log_of_x(g, x);
    if (has_kernel) return rde_solve_tower(f, g, x);
    return rde_base(f, g, x);
}

/* Risch`SPDE[a, b, c, x, n] — Rothstein's SPDE box (Bronstein §6.4, base field
 * C(x)): degree-reducing reduction of a Dq + b q = c to a bounded polynomial
 * equation.  Returns {b', c', m, alpha, beta} (any solution q of degree <= n is
 * q = alpha H + beta with deg(H) <= m and D H + b' H = c'), or $Failed when there
 * is no solution of degree <= n.  Exposed for direct unit testing against the
 * book. */
static Expr* builtin_risch_spde(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 5) return NULL;
    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    Expr* c = res->data.function.args[2];
    Expr* x = res->data.function.args[3];
    Expr* nn = res->data.function.args[4];
    if (x->type != EXPR_SYMBOL || nn->type != EXPR_INTEGER) return NULL;
    RdeCtx C = { NULL, -1, x, x };
    RdeSpde sp;
    if (!rde_spde_field(a, b, c, &C, (long)nn->data.integer, &sp))
        return expr_new_symbol("$Failed");
    Expr* out = expr_new_function(expr_new_symbol("List"),
        (Expr*[]){ expr_copy(sp.b), expr_copy(sp.c), expr_new_integer(sp.m),
                   expr_copy(sp.alpha), expr_copy(sp.beta) }, 5);
    rde_spde_free(&sp);
    return out;
}

/* Risch`PolyRischDENoCancel[b, c, x, n] — the non-cancellation polynomial Risch
 * DE solver (Bronstein §6.5, base field C(x)): solve D q + b q = c for a
 * polynomial q with deg(q) <= n, b != 0.  Returns q or $Failed. */
static Expr* builtin_risch_polyrischde_nocancel(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 4) return NULL;
    Expr* b = res->data.function.args[0];
    Expr* c = res->data.function.args[1];
    Expr* x = res->data.function.args[2];
    Expr* nn = res->data.function.args[3];
    if (x->type != EXPR_SYMBOL || nn->type != EXPR_INTEGER) return NULL;
    RdeCtx C = { NULL, -1, x, x };
    Expr* q = rde_polyrischde_nocancel_field(b, c, &C, (long)nn->data.integer);
    return q ? q : expr_new_symbol("$Failed");
}

/* ================================================================== */
/* Registration.                                                      */
/* ================================================================== */

static void rt_install(const char* name, Expr* (*fn)(Expr*), const char* doc) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    if (doc) symtab_set_docstring(name, doc);
}

void integrate_risch_transcendental_init(void) {
    rt_install("Integrate`RischTranscendental", builtin_rischtranscendental,
        "Integrate`RischTranscendental[f, x] integrates f with respect to x using the\n"
        "recursive transcendental Risch algorithm: a decision\n"
        "procedure over a differential transcendental tower, with rational,\n"
        "logarithmic, exponential, and special-function (Erf, ExpIntegralEi,\n"
        "LogIntegral, PolyLog) cases.  Each case is correct by construction (no\n"
        "differentiation check).  Distinct from Integrate`RischNorman, which is\n"
        "the parallel-Risch (pmint) heuristic.  Out-of-scope integrands\n"
        "(algebraic extensions, non-elementary answers) return unevaluated.");
    rt_install("Risch`RischDE", builtin_risch_rischde,
        "Risch`RischDE[f, g, x] solves the Risch differential equation\n"
        "D[y] + f y == g for a rational y in C(x), or in the transcendental TOWER\n"
        "C(x)(t_1..t_n) inferred from the Log/Exp kernels of f and g (Bronstein\n"
        "Chapter 6, recursive over the tower): normal-denominator reduction, SPDE,\n"
        "and the polynomial non-cancellation solve, with the derivation lifted to\n"
        "the tower.  Returns y (0 when g is 0), or stays unevaluated when no\n"
        "rational solution exists.  Coefficients may be Gaussian (C(i)(x)).");
    rt_install("Risch`SPDE", builtin_risch_spde,
        "Risch`SPDE[a, b, c, x, n] applies Rothstein's SPDE reduction (Bronstein\n"
        "§6.4) to a D[q] + b q == c over C(x): returns {b', c', m, alpha, beta}\n"
        "so any solution q with deg[q] <= n is q = alpha H + beta, deg[H] <= m,\n"
        "D[H] + b' H == c'; $Failed when there is no solution of degree <= n.");
    rt_install("Risch`PolyRischDENoCancel", builtin_risch_polyrischde_nocancel,
        "Risch`PolyRischDENoCancel[b, c, x, n] solves the non-cancellation\n"
        "polynomial Risch differential equation D[q] + b q == c (b != 0) for a\n"
        "polynomial q with deg[q] <= n over C(x) (Bronstein §6.5).  Returns q or\n"
        "$Failed.");
    rt_install("Risch`ElementaryIntegralQ", builtin_elementaryintegralq,
        "Risch`ElementaryIntegralQ[f, x] decides whether f has an antiderivative\n"
        "elementary in x: True if the recursive Risch integrator exhibits an\n"
        "elementary closed form, False if it proves none exists (Bronstein residue\n"
        "criterion / Risch-DE no-solution certificates), or unevaluated when the\n"
        "verdict is outside the transcendental-tower field scope.");
}
