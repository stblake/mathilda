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
#include "integrate.h"           /* g_integrate_depth (nonelem-message gating) */
#include "integrate_risch_rde.h"
#include "risch_util.h"
#include "mtime.h"
#include "risch_tower.h"
#include "risch_field_integrate.h"
#include "risch_singleext.h"
#include "risch_trig_frontend.h"
#include "risch_special.h"
#include "cherry_driver.h"
#include "message.h"           /* mth_message: Quiet / Check */

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "core.h"   /* tc_register_unwind_counter */
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
static bool rt_logcombine_verified(Expr* r, Expr* f, Expr* x);
static RtDecision rt_decide_field(Expr* f, Expr* x);
/* Re-entry guard of the F6 flat-ansatz gate (rt_transcendental_case); an int
 * registered with the TimeConstrained unwind registry, like every flag that
 * brackets a sub-solve: a timeout that longjmps out mid-decision would
 * otherwise leave it, and g_rt_decide_mode, set for the rest of the session. */
static int rt_in_flat_gate = 0;
static int rt_in_logcombine = 0;   /* the log-combine retry's guard (same registry) */
static Expr* rt_transcendental_case(Expr* f, Expr* x) {
    Expr* r = RT_STAGE(rt_log_poly_case, f, x);
    if (r) return r;
    r = RT_STAGE(rt_exp_poly_case, f, x);
    if (r) return r;
    r = RT_STAGE(rt_frac_case, f, x);
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
    r = RT_STAGE(rt_exp_ratreduce_case, f, x);
    if (r) return r;
    r = RT_STAGE(rt_hermite_case, f, x);
    if (r) return r;
    r = RT_STAGE(rt_hyperexp_case, f, x);
    if (r) return r;
    r = RT_STAGE(rt_expsum_case, f, x);   /* direct multi-kernel exponential sums */
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
    r = RT_STAGE(rt_recursive_tower_case, f, x);   /* one-extension recursion (mixed / rational coeff) */
    if (r) return r;
    /* The flat ansätze only find ELEMENTARY antiderivatives, and their size is
     * uncapped (175 unknowns for TrigToExp[u^5 ArcTan[u]^3]; review F6).  When the
     * field decision proves the integrand non-elementary they cannot succeed, so
     * a decline there no longer costs the maximum.  The certificate is the one
     * Risch`ElementaryIntegralQ reports (sound, hunted); the decision runs once,
     * not re-entrantly. */
    bool proved_nonelem = false;
    if (!rt_in_flat_gate) {
        rt_in_flat_gate = 1;
        proved_nonelem = rt_decide_field(f, x) == RT_DEC_NONELEMENTARY;
        rt_in_flat_gate = 0;
    }
    if (!proved_nonelem) {
        r = RT_STAGE(rt_log_tower_case, f, x);   /* flat ansatz fallback (nested log tower) */
        if (r) return r;
        r = RT_STAGE(rt_exp_tower_case, f, x);    /* flat ansatz fallback (nested exp tower) */
        if (r) return r;
    }
    /* Real tangent monomial (Bronstein §5.10): integrate rational functions of a
     * single Tan[u] directly and real, retiring the TrigToExp route below for
     * real-tangent integrands (which it strands at an I-laden complex-log form). */
    r = RT_STAGE(rt_hypertangent_case, f, x);
    if (r) return r;
    r = RT_STAGE(rt_trig_frontend, f, x);
    if (r) return r;
    /* Single-kernel field recursion: a proper part in t plus a k-part whose
     * integral is a new logarithm over k (1/(1 + x) + 1/(x (1 + Log[x]))), which
     * no flat single-kernel case above covers; exactly verified. */
    r = RT_STAGE(rt_recursive_tower_case_single, f, x);
    if (r) return r;
    /* Dependent-logarithm fallback.  When every case above declined, the
     * integrand may carry Q-linearly DEPENDENT logarithmic generators that no
     * single case models — e.g. Log[x/(1+x)] + Log[1+x], two distinct Log
     * kernels that collapse to the single generator Log[x].  Normalize with the
     * log-of-product / log-of-power expansion (Log[a b] -> Log a + Log b,
     * Log[b^p] -> p Log b); the evaluator then recombines the dependent pair
     * (Log[x/(1+x)] + Log[1+x] -> Log[x]), and rt_log_poly_case integrates the
     * collapsed Log[x]/x -> Log[x]^2/2.  Only runs after the original form failed
     * (so it never changes an already-handled case's cleaner output), and
     * guarded against re-entry so the retry itself cannot loop.  The expansion
     * is NOT an identity off the positive axis (Log[x^2] != 2 Log[x] for x < 0,
     * so Int Log[x^2]/x came back as Log[x]^2), so the retry's answer is
     * checked against the ORIGINAL integrand and dropped when it fails. */
    if (!rt_in_logcombine) {
        Expr* fn = rt_eval_own(rt_expand_logs(f));
        if (fn && !expr_eq(fn, f)) {
            rt_in_logcombine = 1;
            r = rt_transcendental_case(fn, x);
            rt_in_logcombine = 0;
            expr_free(fn);
            if (r && !rt_logcombine_verified(r, f, x)) { expr_free(r); r = NULL; }
            if (r) return r;
        } else if (fn) {
            expr_free(fn);
        }
    }
    return NULL;
}

/* PossibleZeroQ[D[r, x] - f]: the check on the dependent-logarithm retry,
 * whose integrand rt_expand_logs changed off the positive axis.  The
 * project's verifier (PossibleZeroQ) also certifies the retries that ARE
 * identities, Log[2 x] = Log[2] + Log[x], which Simplify does not reduce. */
static bool rt_logcombine_verified(Expr* r, Expr* f, Expr* x) {
    Expr* d = rt_eval2("D", expr_copy(r), expr_copy(x));
    if (!d) return false;
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ d, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(f) }, 2) }, 2);
    Expr* z = rt_eval1("PossibleZeroQ", diff);
    bool ok = z && z->type == EXPR_SYMBOL && z->data.symbol.name == intern_symbol("True");
    if (z) expr_free(z);
    return ok;
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
/* Set by rt_decide around its rt_integrate: the elementarity verdict needs an
 * ELEMENTARY antiderivative, and the Cherry stage only ever produces special-
 * function (Ei / li / Erf / dilog) forms, which rt_expr_is_elementary rejects.
 * On a non-elementary integrand Cherry is also where a decline spends its time
 * (Risch`ElementaryIntegralQ[1/(1 + Log[u] + Log[u]^3), u]: 85 s, review F6). */
static int g_rt_elementary_only = 0;   /* tc_register_unwind_counter */

bool rt_decision_only(void) { return g_rt_elementary_only != 0; }

static Expr* rt_integrate_core(Expr* f, Expr* x) {
    Expr* r = rt_rational_case(f, x);
    if (r) return r;
    r = rt_transcendental_case(f, x);
    if (r) return r;
    if (g_rt_elementary_only) return NULL;
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

/* ================================================================== */
/* Quasiquadratic radical admission (Knowles Part I §6).              */
/*                                                                    */
/* The scope gate below rejects every algebraic function of x, which  */
/* is right for Sqrt[x] — an algebraic extension the recursive Risch  */
/* algorithm is not a decision procedure over.  But a fractional      */
/* power of a TRANSCENDENTAL KERNEL, Log[x]^(3/2), is a different     */
/* animal: it is Knowles' quasiquadratic case, and knowles_erf.c      */
/* solves it by working in s = Sqrt[kernel].  That engine is only     */
/* reachable through rt_integrate, so the gate was hiding it — and    */
/* only from the *already-reduced* spelling, since the merged form    */
/* E^(1/2 Log[Log x] - 1/Log x)/(x Log x^2) carries no radical and    */
/* sailed through.  Two spellings of one integrand, two answers.      */
/*                                                                    */
/* Admission works by writing the fractional part of each such power  */
/* back into exp-log form, g^(p/q) = g^n E^(r Log g) with n = the     */
/* integer part and r = p/q - n in (0,1) — the exact inverse of       */
/* knowles' collapse_exp_of_log, and precisely the merged spelling.   */
/*                                                                    */
/* SOUNDNESS.  The rewritten tower hides an algebraic relation: the   */
/* monomial E^(r Log g) is not transcendental over the lower field    */
/* (its q-th power is g^p).  Nothing on this path may therefore be    */
/* trusted by construction — the answer is returned ONLY behind a     */
/* diff-back check against the ORIGINAL integrand, and the            */
/* non-elementary decision half never runs here.  A verified          */
/* antiderivative is correct whatever field it was found in; an       */
/* unverified one is discarded.                                       */
/* ================================================================== */

/* Is `e` a Power[g, Rational[p, q]] whose base is a transcendental kernel of x
 * (it contains an Exp or a Log of x) rather than a rational function of x?
 * Sqrt[x] and Sqrt[1 + x^2] are NOT this; Sqrt[Log[x]] and (1 + E^x)^(1/2) are. */
static bool rt_frac_power_of_kernel(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!rt_head_is(e, "Power") || e->data.function.arg_count != 2) return false;
    Expr* b = e->data.function.args[0];
    Expr* p = e->data.function.args[1];
    if (!p || p->type != EXPR_FUNCTION || !rt_head_is(p, "Rational")) return false;
    if (rt_free_of_x(b, x)) return false;
    return rt_find_exp_of_x(b, x) != NULL || rt_find_log_of_x(b, x) != NULL;
}

/* Every algebraic-of-x site in `e` is a fractional power of a transcendental
 * kernel — i.e. admitting `e` costs us nothing but the Knowles case.  A single
 * Sqrt[x] / Surd / Root anywhere disqualifies the whole integrand. */
static bool rt_algebraic_only_in_kernels(Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return true;

    /* An admissible fractional kernel power.  Accept the node itself, but its
     * base must not hide a bare radical of x underneath (Sqrt[Log[Sqrt[x]]]). */
    if (rt_frac_power_of_kernel(e, x))
        return !rt_has_algebraic_of_x(e->data.function.args[0], x);

    /* Surd / Root / AlgebraicNumber of x is never admissible. */
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if ((h == intern_symbol("Surd") || h == intern_symbol("Root")
             || h == intern_symbol("AlgebraicNumber")) && !rt_free_of_x(e, x))
            return false;
    }
    /* Nor is a fractional power of an x-dependent base that is NOT a kernel:
     * Sqrt[x], Sqrt[1 + x^2] — the genuine algebraic extensions the gate is for. */
    if (rt_head_is(e, "Power") && e->data.function.arg_count == 2) {
        Expr* p = e->data.function.args[1];
        if (p && p->type == EXPR_FUNCTION && rt_head_is(p, "Rational")
            && !rt_free_of_x(e->data.function.args[0], x))
            return false;
    }

    if (!rt_algebraic_only_in_kernels(e->data.function.head, x)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!rt_algebraic_only_in_kernels(e->data.function.args[i], x)) return false;
    return true;
}

/* Rewrite every fractional power of a transcendental kernel as
 * g^n E^(r Log[g]), n = Floor[p/q], r = p/q - n.  Exact for the branch the
 * rest of the machinery already assumes (the same one collapse_exp_of_log
 * inverts), and a no-op elsewhere.  Returns an owned copy. */
static Expr* rt_kernel_radicals_to_explog(Expr* e, Expr* x) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    size_t k = e->data.function.arg_count;
    Expr* nh = rt_kernel_radicals_to_explog(e->data.function.head, x);
    Expr** na = malloc((k ? k : 1) * sizeof(Expr*));
    for (size_t i = 0; i < k; i++)
        na[i] = rt_kernel_radicals_to_explog(e->data.function.args[i], x);
    Expr* r = expr_new_function(nh, na, k);
    free(na);
    if (!rt_frac_power_of_kernel(r, x)) return r;

    Expr* g = r->data.function.args[0];
    Expr* p = r->data.function.args[1];
    Expr* n = rt_eval1("Floor", expr_copy(p));                    /* integer part  */
    Expr* frac = rt_eval1("Together", expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ expr_copy(p), expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(n) }, 2) }, 2));   /* r = p - n */
    if (!n || !frac) { if (n) expr_free(n); if (frac) expr_free(frac); return r; }

    Expr* expo = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ frac, expr_new_function(expr_new_symbol("Log"),
                             (Expr*[]){ expr_copy(g) }, 1) }, 2);          /* r Log[g] */
    Expr* out = expr_new_function(expr_new_symbol("Times"),
        (Expr*[]){ expr_new_function(expr_new_symbol("Power"),
                       (Expr*[]){ expr_copy(g), n }, 2),                   /* g^n      */
                   expr_new_function(expr_new_symbol("Power"),
                       (Expr*[]){ expr_new_symbol("E"), expo }, 2) }, 2);  /* E^(r Lg) */
    expr_free(r);
    return out;
}

/* The admission itself: rewrite, integrate constructively, and return the
 * result ONLY if it differentiates back to the original integrand. */
static Expr* rt_try_kernel_radical(Expr* f, Expr* x) {
    if (!rt_algebraic_only_in_kernels(f, x)) return NULL;
    Expr* g = rt_kernel_radicals_to_explog(f, x);
    if (!g) return NULL;
    if (expr_eq(g, f)) { expr_free(g); return NULL; }   /* nothing to admit */

    Expr* r = rt_integrate(g, x);
    expr_free(g);
    if (!r) return NULL;
    if (!rt_verify_antideriv(r, f, x)) { expr_free(r); return NULL; }
    return r;
}

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
/* Every Log[a] (a depending on x) rewritten as Log[Together[a]], innermost
 * first; returns a new expression (f is borrowed).  An identity (the same
 * value), and NECESSARY for the field decision: with a log argument left as an
 * uncombined sum -- Log[(u^2 - 1)/(2 u) + (u^2 + 1)/(2 u)], which is Log[u],
 * as Bronstein's Euler substitution produced it (integration shootout
 * 2026-11-12, bbottom-3-9) -- the decision returned a FALSE "not elementary"
 * for D[Log[(u^2 - 1)/(2 u)] Log[u], u], while the canonical spelling decides
 * True and integrates. */
static Expr* rt_canon_log_args(const Expr* e, const Expr* x) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = rt_canon_log_args(e->data.function.args[i], x);
    Expr* r = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    if (rt_head_is(r, "Log") && n == 1 && !rt_free_of_x(r->data.function.args[0], (Expr*)x)) {
        Expr* t = rt_eval1("Together", expr_copy(r->data.function.args[0]));
        if (t) {
            Expr* l = expr_new_function(expr_new_symbol("Log"), (Expr*[]){ t }, 1);
            expr_free(r);
            r = l;
        }
    }
    return r;
}

static RtDecision rt_decide_field_canon(Expr* f, Expr* x);
/* Pi as an ordinary transcendental parameter (campaign step 3b): every Erf-type
 * integrand carries 1/Sqrt[Pi] from D[Erf].  Q(Sqrt[Pi]) is isomorphic to
 * Q(a), a symbol, and a verdict over Q(a) is accepted already (Prop. T1); a
 * transcendental Sqrt[Pi] avoids the exceptional set.  What makes Pi unlike a
 * parameter is Log[-1] = I Pi: a logarithm (or a trigonometric kernel, through
 * the complex exponentials) can produce it inside the tower.  So: Pi is the
 * only transcendental constant, and f has no Log, LogIntegral, or trigonometric
 * / hyperbolic / inverse kernel of x. */
static bool rt_has_logtrig_of_x(const Expr* e, Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        static const char* const heads[] = { "Log", "LogIntegral", "Sin", "Cos", "Tan", "Cot", "Sec",
            "Csc", "ArcSin", "ArcCos", "ArcTan", "ArcCot", "ArcSec", "ArcCsc", "Sinh", "Cosh", "Tanh",
            "Coth", "Sech", "Csch", "ArcSinh", "ArcCosh", "ArcTanh", "ArcCoth", "ArcSech", "ArcCsch", NULL };
        const char* h = e->data.function.head->data.symbol.name;
        for (int i = 0; heads[i]; i++)
            if (h == intern_symbol(heads[i]) && !rt_free_of_x((Expr*)e, x)) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rt_has_logtrig_of_x(e->data.function.args[i], x)) return true;
    return false;
}
static bool rt_pi_is_a_parameter(Expr* f, Expr* x) {
    if (rt_has_logtrig_of_x(f, x)) return false;
    Expr* a = expr_new_symbol("rmPiParam");
    Expr* fa = rt_eval2("ReplaceAll", expr_copy(f), expr_new_function(expr_new_symbol("Rule"),
        (Expr*[]){ expr_new_symbol("Pi"), expr_copy(a) }, 2));
    bool only = fa && !rt_has_transcendental_constant(fa, x);
    if (fa) expr_free(fa);
    expr_free(a);
    return only;
}

static RtDecision rt_decide_field(Expr* f, Expr* x) {
    /* Scope guard: an algebraic function of x (radical / Surd / Root) puts the
     * integrand outside the purely-transcendental tower the field decision is a
     * decision procedure over — decline to UNKNOWN rather than emit an unsound
     * certificate.  Checked on the ORIGINAL f (Sqrt[x] survives inside Cos[Sqrt[x]]
     * and E^Sqrt[x] alike), before exponentialization. */
    if (rt_has_algebraic_of_x(f, x)) { rde_cert_mark("decide: algebraic function of x"); return RT_DEC_UNKNOWN; }
    Expr* fc = rt_canon_log_args(f, x);       /* canonical log arguments: see the helper */
    RtDecision vc = rt_decide_field_canon(fc, x);
    /* a transcendental constant (Log[2], also one rt_canon_log_args split off
     * Log[2 u]) enlarges the constant field beyond Q(i): a "no" over Q(i) is
     * not a proof (rt_has_transcendental_constant) */
    if (vc == RT_DEC_NONELEMENTARY &&
        (rt_has_transcendental_constant(f, x) || rt_has_transcendental_constant(fc, x))
        && !rt_pi_is_a_parameter(f, x)) {
        vc = RT_DEC_UNKNOWN;
        rde_cert_mark("decide: transcendental constant");
    }
    expr_free(fc);
    return vc;
}

/* The answer check of the cascade gate and the elementarity verdict:
 * PossibleZeroQ[D[r, x] - f], or failing that the same of its PowerExpand.
 * The second pass forgives the principal-branch convention only (an answer
 * that holds where Log[u^k] == k Log[u], the same convention DSolve's verifier
 * and the GammaPower stage adopt; the li producers no longer need it since
 * v0.6601, which writes x/Log[x]^2 as 2 ExpIntegralEi[2 Log[x]] - x^2/Log[x],
 * valid everywhere, not 2 LogIntegral[x^2] - x^2/Log[x]); a wrong answer (review F1: residual -0.186 on the
 * real line) fails both.  An antiderivative valid on an open set still proves
 * elementarity.  rt_logcombine_verified stays strict: it guards a retry that
 * is valid off the positive axis only. */
bool rt_answer_verified(Expr* r, Expr* f, Expr* x) {
    Expr* d = rt_eval2("D", expr_copy(r), expr_copy(x));
    if (!d) return false;
    Expr* diff = expr_new_function(expr_new_symbol("Plus"),
        (Expr*[]){ d, expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), expr_copy(f) }, 2) }, 2);
    /* Cheapest first; PossibleZeroQ's cost depends on the spelling, and no one
     * spelling is cheap for every answer (Tan/Tanh hunts, 2026-11-20):
     *   1. exact: Together of the TrigToExp residual, logs normalised -- 0.04 s
     *      where PossibleZeroQ on the TrigToExp form ran past 30 s;
     *   2. PossibleZeroQ on D[r] - TrigToExp[f], an exact rewrite in a trig
     *      front-end answer's own exponentials: 0.6 s against 12 s for (3) on a
     *      Tanh answer that spells E^Tanh[1 + u] another way;
     *   3. PossibleZeroQ on the residual itself (what rt_logcombine_verified
     *      checks);
     *   4. PossibleZeroQ on the PowerExpand (principal-branch convention).
     * Only a True is taken from any of them. */
    if (rt_exact_zero_trigexp(diff, x)) { expr_free(diff); return true; }
    {
        Expr* fe = rt_eval1("TrigToExp", expr_copy(f));
        bool okt = false;
        if (fe && !expr_eq(fe, f)) {
            Expr* dt = expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ rt_eval2("D", expr_copy(r), expr_copy(x)),
                    expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-1), fe }, 2) }, 2);   /* adopts fe */
            Expr* zt = rt_eval1("PossibleZeroQ", dt);     /* adopts dt */
            okt = zt && rt_is_true(zt);
            if (zt) expr_free(zt);
        } else if (fe) expr_free(fe);
        if (okt) { expr_free(diff); return true; }
    }
    {
        Expr* z0 = rt_eval1("PossibleZeroQ", expr_copy(diff));
        bool ok0 = z0 && rt_is_true(z0);
        if (z0) expr_free(z0);
        if (ok0) { expr_free(diff); return true; }
    }
    Expr* z = rt_eval1("PossibleZeroQ", rt_eval1("PowerExpand", diff));
    bool ok = z && rt_is_true(z);
    if (z) expr_free(z);
    return ok;
}

/* Q (a field-integrator answer in the tower variables) written back in the
 * kernels t_i -> kernel_i, then checked by PossibleZeroQ[D[Q, x] - f]. */
static bool rt_tower_answer_verified(Expr* Q, RtTower* T, Expr* f, Expr* x) {
    Expr** rs = malloc((T->n ? T->n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < T->n; i++)
        rs[i] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(T->t[i]), expr_copy(T->kernel[i]) }, 2);
    Expr* rules = expr_new_function(expr_new_symbol("List"), rs, T->n);
    free(rs);
    Expr* back = rt_eval2("ReplaceAll", expr_copy(Q), rules);
    bool ok = back && rt_answer_verified(back, f, x);
    if (back) expr_free(back);
    return ok;
}

/* The field decision over the differential tower of fx (borrowed): build it,
 * gate on a rational function of the tower, integrate in decide mode and read
 * the verdict.  tan_route asks for the REAL tangent tower (Bronstein §5.10,
 * the RDE over a hypertangent field, review F8): it is decided only when the
 * tower has exactly one hypertangent level and no Liouvillian primitive, so
 * that tan u is a monomial -- two tangents can be dependent (Tan[2 x] is
 * rational in Tan[x]) and an ArcTan-type primitive can make Tan algebraic. */
static RtDecision rt_decide_tower(Expr* fx, Expr* x, bool tan_route) {
    RtTower T;
    /* min_n = 1: the decision routes SINGLE-kernel integrands (E^x/x, E^(x^2),
     * 1/Log[x]) through the same authoritative field integrator too — the flat-tower
     * integrate cases that normally serve n=1 do not carry the decision certificates. */
    /* every UNKNOWN exit below is marked (the real-tangent attempt excepted: it
     * falls back to the exponential route), so a silent decline is visible to
     * tools/check_bronstein_completeness.py and to the hunts' mark sets */
    if (!rt_tower_build_decide(fx, x, &T, 1)) {
        if (!tan_route) rde_cert_mark("decide: tower could not be built");
        rt_tower_free(&T); return RT_DEC_UNKNOWN;
    }
    if (tan_route) {
        size_t ntan = 0; bool scope = true;
        for (size_t i = 0; i < T.n; i++) {
            if (T.kind[i] == RT_TAN) { ntan++; if (T.tsg[i] != 1) scope = false; }
            if (T.kind[i] == RT_PRIM) scope = false;
        }
        if (ntan != 1 || !scope) { rt_tower_free(&T); return RT_DEC_UNKNOWN; }
    }
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
        if (!gate && !tan_route) {
            if (getenv("MATHILDA_RDE_CERT_DEBUG")) {
                char* a1 = expr_to_string(F);
                fprintf(stderr, "[rt-decide-gate] F=%.600s\n", a1);
                free(a1);
            }
            rde_cert_mark("decide: integrand not rational in the tower");
        }
        if (gate) {
            int        save_mode = g_rt_decide_mode;
            RtDecision save_dec  = g_rt_decision;
            g_rt_decide_mode = 1;
            g_rt_decision    = RT_DEC_UNKNOWN;
            Expr* rem = NULL;
            unsigned long prim0 = g_rt_prim_unevaluated;
            Expr* Q = rt_field_integrate(F, &T, (long)T.n - 1, x, &rem);
            /* Elementary only on positive evidence: Q, written back in the
             * kernels, must differentiate to the integrand.  An unverified Q
             * made Log[u] Log[1 + u] and Log[u]/(u + a) "elementary" (review
             * F2-F4, 2026-11-16). */
            if (Q && !rem) {
                bool okq = rt_tower_answer_verified(Q, &T, fx, x);
                /* over an eliminated dependent logarithm (Log[-1 - u] = Log[1 + u]
                 * + I Pi, a constant of the positive reals; campaign step 8) the
                 * kernel form of Q fails off that region; D_tower[Q] = F, exact in
                 * the tower, is then the evidence -- F is the integrand on that
                 * open set, and an antiderivative there proves elementarity */
                if (!okq && T.nlm > 0) {
                    Expr* dq = rt_tower_deriv(Q, &T, x);
                    okq = dq && rt_residual_zero(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){ dq,
                        expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), expr_copy(F) }, 2) }, 2),
                        Q, false);
                }
                verdict = okq ? RT_DEC_ELEMENTARY : RT_DEC_UNKNOWN;
                if (!okq) rde_cert_mark("decide: antiderivative not verified");
            }
            else if (Q && rem)                           verdict = RT_DEC_NONELEMENTARY;
            else if (g_rt_decision == RT_DEC_NONELEMENTARY) verdict = RT_DEC_NONELEMENTARY;
            else {
                verdict = RT_DEC_UNKNOWN;
                rde_cert_mark("decide: field integrator gave no verdict");
            }
            /* A certificate reached through an unevaluated polynomial primitive
             * (PolynomialGCD over a coefficient ring no engine covered) proves
             * nothing: the algorithm computed with an unevaluated head as if it
             * were a polynomial.  Demote it to UNKNOWN (see risch_util.c). */
            if (getenv("MATHILDA_RDE_CERT_DEBUG"))   /* the verdict path, for tracing */
                fprintf(stderr, "[rt-decide] Q=%d rem=%d dec=%d verdict=%d prim=%lu dep=%d\n",
                        Q != NULL, rem != NULL, (int)g_rt_decision, (int)verdict,
                        g_rt_prim_unevaluated - prim0, (int)T.dependent);
            if (verdict == RT_DEC_NONELEMENTARY && g_rt_prim_unevaluated != prim0) {
                verdict = RT_DEC_UNKNOWN;
                rde_cert_mark("decide: a polynomial primitive did not evaluate");
            }
            /* a certificate presupposes algebraically independent monomials */
            if (verdict == RT_DEC_NONELEMENTARY && T.dependent) {
                verdict = RT_DEC_UNKNOWN;
                rde_cert_mark("decide: dependent monomials");
            }
            /* ... and a Liouvillian primitive level (Erf, Ei, li) is never checked
             * for independence (Erfc = 1 - Erf, LogIntegral[u] = Ei[Log[u]]), and the
             * field integrator models it as an exponential: E^u (Erf[u] + Erfc[u])
             * and 1/u + E^u (Ei[Log[u]] - LogIntegral[u]) were certified
             * non-elementary (review step 0.2, v0.6603).  No certificate over such
             * a tower until the primitive top is implemented. */
            /* Since step 3b/3c (v0.6621) the field integrator has the primitive
             * branch, and a certificate needs every such level proved a
             * monomial with a non-elementary Int eta below it. */
            if (verdict == RT_DEC_NONELEMENTARY) {
                bool anyprim = false;
                for (size_t i = 0; i < T.n; i++) if (T.kind[i] == RT_PRIM) anyprim = true;
                if (anyprim && rt_tower_prims_certified(&T, x) != 1) {
                    verdict = RT_DEC_UNKNOWN;
                    rde_cert_mark("decide: Liouvillian level not proved a monomial");
                }
            }
            if (Q) expr_free(Q);
            if (rem) expr_free(rem);
            g_rt_decide_mode = save_mode;
            g_rt_decision    = save_dec;
        }
    }
    else if (!tan_route) rde_cert_mark("decide: kernels left after substitution");
    if (F) expr_free(F);
    rt_tower_free(&T);
    return verdict;
}

static RtDecision rt_decide_field_canon(Expr* f, Expr* x) {
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
    /* Real tangent route first (review F8): a real integrand rational in one
     * Tan[u] and the exp/log kernels is decided over the hypertangent tower
     * itself, where the Risch DE runs Bronstein's delta = 2 boxes; TrigToExp
     * would turn E^Tan[u] into an exponential of a rational function of
     * E^(I u), outside any tower.  Undecided there: the exponential route. */
    if (rt_free_of_head(f, "Complex")) {
        Expr* pt = rt_powers_to_exp(f, x);
        Expr* ft = pt ? rt_expand_exp_sums(pt) : NULL;
        if (pt) expr_free(pt);
        RtDecision vt = ft ? rt_decide_tower(ft, x, true) : RT_DEC_UNKNOWN;
        if (ft) expr_free(ft);
        if (vt != RT_DEC_UNKNOWN) return vt;
    }
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
    RtDecision verdict = rt_decide_tower(fx, x, false);
    expr_free(fx);
    return verdict;
}

/* Complete elementary-integrability verdict for f w.r.t. x.  True side: exhibit an
 * elementary antiderivative via the full integrator (an existence proof — correct
 * by construction).  False side: the field decision's authoritative non-elementary
 * certificate.  Otherwise UNKNOWN.  Checking True first makes an exhibited
 * antiderivative dominate any certificate, so the verdict is robust. */
static RtDecision rt_decide(Expr* f, Expr* x) {
    int save_eo = g_rt_elementary_only;
    g_rt_elementary_only = 1;
    double tp0 = rt_profile_on() ? mth_mono_seconds() : 0;
    Expr* anti = rt_integrate(f, x);
    g_rt_elementary_only = save_eo;
    if (rt_profile_on()) fprintf(stderr, "[rt-stage] decide: integrator %.3f s -> %s\n",
                                 mth_mono_seconds() - tp0, anti ? "answer" : "none");
    if (anti) {
        /* an exhibited antiderivative proves elementarity only once checked:
         * an unverified wrong one (review F1/F2) is no evidence */
        double tv0 = rt_profile_on() ? mth_mono_seconds() : 0;
        bool elem = rt_expr_is_elementary(anti) && rt_answer_verified(anti, f, x);
        if (rt_profile_on()) fprintf(stderr, "[rt-stage] decide: answer check %.3f s -> %d\n",
                                     mth_mono_seconds() - tv0, (int)elem);
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
     * (Sqrt[2]) is a legitimate transcendental coefficient and is NOT flagged.
     *
     * One exception is admitted: when every algebraic site is a fractional power
     * of a TRANSCENDENTAL kernel (Log[x]^(3/2)), that is Knowles' quasiquadratic
     * case and the erf engine solves it in s = Sqrt[kernel].  See
     * rt_try_kernel_radical — the answer is gated on a diff-back check, and the
     * non-elementary decision half deliberately does not run on that path. */
    if (rt_has_algebraic_of_x(f, x)) {
        arith_warnings_mute_push();
        Expr* kr = rt_try_kernel_radical(f, x);
        arith_warnings_mute_pop();
        return kr;
    }

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
    /* canonical log arguments (rt_canon_log_args): an identity, and the tower
     * and the decision both need it */
    Expr* fcanon = rt_canon_log_args(f, x);
    f = fcanon;
    Expr* result = rt_integrate(f, x);
    /* A real exponential the trig route's ExpToTrig wrote Cosh[a] + Sinh[a]
     * (D[ArcTan[E^u]^3] came out ArcTan[E^u]^2 ArcTan[Cosh[u] + Sinh[u]]): put it
     * back when the integrand has no hyperbolic function.  Exact identities. */
    if (result && rt_free_of_head(f, "Cosh") && rt_free_of_head(f, "Sinh")
        && !rt_free_of_head(result, "Sinh")) {
        Expr* pa = expr_new_symbol("Risch`Private`rtChSh");
        Expr* rl = expr_new_function(expr_new_symbol("List"), (Expr*[]){
            expr_new_function(expr_new_symbol("RuleDelayed"), (Expr*[]){
                expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
                    expr_new_function(expr_new_symbol("Cosh"), (Expr*[]){ expr_new_function(expr_new_symbol("Pattern"),
                        (Expr*[]){ expr_copy(pa), expr_new_function(expr_new_symbol("Blank"), NULL, 0) }, 2) }, 1),
                    expr_new_function(expr_new_symbol("Sinh"), (Expr*[]){ expr_new_function(expr_new_symbol("Pattern"), (Expr*[]){ expr_copy(pa), expr_new_function(expr_new_symbol("Blank"), NULL, 0) }, 2) }, 1) }, 2),
                expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_new_symbol("E"), expr_copy(pa) }, 2) }, 2),
            expr_new_function(expr_new_symbol("RuleDelayed"), (Expr*[]){
                expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
                    expr_new_function(expr_new_symbol("Cosh"), (Expr*[]){ expr_new_function(expr_new_symbol("Pattern"),
                        (Expr*[]){ expr_copy(pa), expr_new_function(expr_new_symbol("Blank"), NULL, 0) }, 2) }, 1),
                    expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1),
                        expr_new_function(expr_new_symbol("Sinh"), (Expr*[]){ expr_new_function(expr_new_symbol("Pattern"), (Expr*[]){ expr_copy(pa), expr_new_function(expr_new_symbol("Blank"), NULL, 0) }, 2) }, 1) }, 2) }, 2),
                expr_new_function(expr_new_symbol("Power"), (Expr*[]){ expr_new_symbol("E"),
                    expr_new_function(expr_new_symbol("Times"), (Expr*[]){ expr_new_integer(-1), expr_copy(pa) }, 2) }, 2) }, 2) }, 2);
        expr_free(pa);
        Expr* clean = rt_eval2("ReplaceRepeated", expr_copy(result), rl);
        if (clean) { expr_free(result); result = clean; }
    }
    /* Decision half (P3): when the recursive algorithm produced no elementary
     * antiderivative, consult the authoritative field decision.  If it PROVES the
     * integrand has no elementary integral (a non-constant residue, a Risch DE with
     * no rational solution, or the Dc!=0 primitive certificate — never a bounded-
     * ansatz decline), report it Mathematica-style.  A stray certificate from an
     * abandoned attempt cannot fire here: rt_decide_field re-derives the verdict
     * with elementary-result precedence. */
    /* The nonelem diagnostic is a user-facing message about the ORIGINAL
     * integrand, so emit it only at the outermost integration (g_integrate_depth
     * 0 for a direct Integrate`RischTranscendental[f, x] call, 1 for the top of
     * the Integrate cascade).  A sub-integral spawned by an internal
     * substitution — e.g. DerivativeDivides' u = Log[x], reaching this decision
     * at depth >= 2 — is not the user's integrand: naming its gensym variable
     * (Integrate`DerivativeDivides`u$N) in a message would be meaningless, so we
     * stay silent there and let the outermost frame speak. */
    if (!result && g_integrate_depth <= 1 && g_integrate_quiet == 0
        && rt_decide_field(f, x) == RT_DEC_NONELEMENTARY) {
        /* Shared with ParallelMixedTower's certificate path so the two methods
         * stay byte-identical and cannot double-print in the Automatic cascade.
         * The cheap depth/quiet gate above stays here to avoid the expensive
         * rt_decide_field call on internal recursions; the helper re-checks it. */
        integrate_announce_nonelementary(res->data.function.args[0], x);
    }
    arith_warnings_mute_pop();
    expr_free(fcanon);
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
        /* through the funnel: a raw fprintf ignored Quiet and Check */
        mth_message("Risch`ElementaryIntegralQ", "undec",
            "Cannot decide elementary integrability of %s with respect to %s.", fs ? fs : "?", xs ? xs : "?");
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

/* Risch`RischDECertified[f, g, x] / [f, g, x, {kernels}] — the tri-state
 * Risch DE used by the Bronstein package (src/internal/bronstein/RDE.m), whose
 * non-elementary certificates need a PROVEN "no solution".  Solves
 * D[y] + f y = g in C(x) or in the tower of the Log/Exp kernels of f, g and the
 * optional kernel list, inside a certification scope (rde_cert_begin), and
 * returns
 *     y             a solution (the same one Risch`RischDE finds);
 *     "NoSolution"  no solution exists in that field -- every NULL path taken
 *                   was an authoritative one (Bronstein Ch.6 theorems);
 *     $Failed       undecided: an out-of-scope input or an unimplemented
 *                   sub-case (Risch`RischDECertifiedReason[] names the first).
 * Inputs must be rational in x and the tower kernels (numeric constants
 * allowed); anything else is $Failed.  The five-argument form
 * [f, g, x, {kernels}, {params}] also admits the symbols params as
 * transcendental constants (D = 0): the verdict is then over Q(params)
 * (Bronstein Phase 4.4). */
static const char* g_rde_cert_last_reason = NULL;

/* Does e contain a circular trig function (Tan, Cot, Sec, Csc, Sin, Cos) of an
 * x-dependent argument?  Such a term is rational in the tangent tower variable
 * (Sec[x]^2 = 1 + Tan[x]^2), which rt_subst_kernels rationalises. */
static bool rt_has_tan_of_x(const Expr* e, Expr* x) {
    if (e->type != EXPR_FUNCTION) return false;
    static const char* heads[] = { "Tan", "Cot", "Sec", "Csc", "Sin", "Cos" };
    for (size_t h = 0; h < sizeof heads / sizeof heads[0]; h++)
        if (rt_head_is(e, heads[h]) && e->data.function.arg_count == 1
            && !rt_free_of_x(e->data.function.args[0], x)) return true;
    if (rt_has_tan_of_x(e->data.function.head, x)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rt_has_tan_of_x(e->data.function.args[i], x)) return true;
    return false;
}

/* Does e contain a Liouvillian primitive kernel (Erf, Erfi, Erfc, Ei, li) of an
 * x-dependent argument? */
static bool rt_has_prim_of_x(Expr* e, Expr* x) {
    Expr** pr = NULL; size_t np = 0, cap = 0;
    rt_collect_primitives(e, x, &pr, &np, &cap);
    for (size_t i = 0; i < np; i++) expr_free(pr[i]);
    free(pr);
    return np > 0;
}

static Expr* builtin_risch_rischde_certified(Expr* res) {
    size_t ac = (res->type == EXPR_FUNCTION) ? res->data.function.arg_count : 0;
    if (ac < 3 || ac > 5) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* g = res->data.function.args[1];
    Expr* x = res->data.function.args[2];
    Expr* extra = (ac >= 4) ? res->data.function.args[3] : NULL;
    Expr* params = (ac == 5) ? res->data.function.args[4] : NULL;
    if (x->type != EXPR_SYMBOL) return NULL;
    if (extra && !rt_head_is(extra, "List")) return NULL;
    if (params) {
        if (!rt_head_is(params, "List")) return NULL;
        for (size_t i = 0; i < params->data.function.arg_count; i++) {
            Expr* p = params->data.function.args[i];
            if (p->type != EXPR_SYMBOL || p->data.symbol.name == x->data.symbol.name) return NULL;
        }
        if (params->data.function.arg_count == 0) params = NULL;
    }
    bool has_kernel = rt_find_exp_of_x(f, x) || rt_find_log_of_x(f, x)
                   || rt_find_exp_of_x(g, x) || rt_find_log_of_x(g, x)
                   || (extra && (rt_find_exp_of_x(extra, x) || rt_find_log_of_x(extra, x)))
                   /* a tangent kernel: the hypertangent tower (review F8) */
                   || rt_has_tan_of_x(f, x) || rt_has_tan_of_x(g, x)
                   || (extra && rt_has_tan_of_x(extra, x))
                   /* a Liouvillian primitive kernel (Erf, Ei, li; step 3b) */
                   || rt_has_prim_of_x(f, x) || rt_has_prim_of_x(g, x)
                   || (extra && rt_has_prim_of_x(extra, x));
    RdeCertState saved = rde_cert_begin();
    rde_cert_set_params(params);
    Expr* y = NULL;
    if (has_kernel) {
        y = rde_solve_tower_ext(f, g, x, (extra && extra->data.function.arg_count) ? extra : NULL);
    } else {
        Expr* xl = expr_new_function(expr_new_symbol("List"), (Expr*[]){ expr_copy(x) }, 1);
        if (rde_rational_in(f, xl) && rde_rational_in(g, xl)) y = rde_base(f, g, x);
        else rde_cert_mark("base: coefficients not rational in x");
        expr_free(xl);
    }
    rde_cert_set_params(NULL);
    bool unc = rde_cert_uncertain();
    g_rde_cert_last_reason = y ? NULL : (unc ? rde_cert_reason() : NULL);
    rde_cert_end(saved);
    if (y) return y;
    return unc ? expr_new_symbol("$Failed") : expr_new_string("NoSolution");
}

/* Risch`ParametricRischDE[f, {g0, g1, ..., gm}, x] / [..., {kernels}] /
 * [..., {kernels}, {params}] — the parametric Risch DE ([SI] ch. 7, campaign
 * step 3a): find y in the field and constants c1..cm with
 *     D[y] + f y == g0 + c1 g1 + ... + cm gm.
 * Returns {y, {c1, ..., cm}} (a free constant is 0; the identity is checked
 * exactly), "NoSolution" when no constants admit a solution (proved), or $Failed
 * when undecided (Risch`RischDECertifiedReason[] names why).  The field is that
 * of Risch`RischDECertified: C(x), or the tower of the Log/Exp/Tan kernels of f,
 * the g_i and the optional kernel list, over Q(params). */
static Expr* builtin_risch_parametric_rischde(Expr* res) {
    size_t ac = (res->type == EXPR_FUNCTION) ? res->data.function.arg_count : 0;
    if (ac < 3 || ac > 5) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* gl = res->data.function.args[1];
    Expr* x = res->data.function.args[2];
    Expr* extra = (ac >= 4) ? res->data.function.args[3] : NULL;
    Expr* params = (ac == 5) ? res->data.function.args[4] : NULL;
    if (x->type != EXPR_SYMBOL || !rt_head_is(gl, "List") || gl->data.function.arg_count < 1)
        return NULL;
    if (extra && !rt_head_is(extra, "List")) return NULL;
    if (params) {
        if (!rt_head_is(params, "List")) return NULL;
        for (size_t i = 0; i < params->data.function.arg_count; i++) {
            Expr* p = params->data.function.args[i];
            if (p->type != EXPR_SYMBOL || p->data.symbol.name == x->data.symbol.name) return NULL;
        }
    }
    size_t m = gl->data.function.arg_count - 1;
    /* the unknown constants rmPc1..rmPcm, which the input must not mention */
    Expr** us = malloc((m ? m : 1) * sizeof(Expr*));
    for (size_t i = 0; i < m; i++) {
        char nm[32]; snprintf(nm, sizeof(nm), "rmPc%zu", i + 1);
        us[i] = expr_new_symbol(nm);
    }
    Expr* unk = expr_new_function(expr_new_symbol("List"), us, m);
    free(us);
    {   /* argument by argument: FreeQ of the call itself would re-enter it */
        Expr* w = expr_new_symbol("rmPw");
        bool clash = false;
        for (size_t a = 0; a < ac && !clash; a++) {
            Expr* arg = res->data.function.args[a];
            if (!rt_free_of_x(arg, w)) clash = true;
            for (size_t i = 0; i < m && !clash; i++)
                if (!rt_free_of_x(arg, unk->data.function.args[i])) clash = true;
        }
        expr_free(w);
        if (clash) { expr_free(unk); return NULL; }
    }
    /* g = g0 + Sum c_i g_i */
    Expr** tm = malloc((m + 1) * sizeof(Expr*));
    tm[0] = expr_copy(gl->data.function.args[0]);
    for (size_t i = 0; i < m; i++)
        tm[i + 1] = expr_new_function(expr_new_symbol("Times"),
            (Expr*[]){ expr_copy(unk->data.function.args[i]), expr_copy(gl->data.function.args[i + 1]) }, 2);
    Expr* g = rt_eval1("Expand", expr_new_function(expr_new_symbol("Plus"), tm, m + 1));
    free(tm);
    /* the unknowns are constants of the field, like the declared params */
    size_t np = params ? params->data.function.arg_count : 0;
    Expr** pa = malloc((np + m + 1) * sizeof(Expr*));
    for (size_t i = 0; i < np; i++) pa[i] = expr_copy(params->data.function.args[i]);
    for (size_t i = 0; i < m; i++) pa[np + i] = expr_copy(unk->data.function.args[i]);
    Expr* allp = expr_new_function(expr_new_symbol("List"), pa, np + m);
    free(pa);
    bool has_kernel = false;
    for (size_t i = 0; i < gl->data.function.arg_count && !has_kernel; i++) {
        Expr* gi = gl->data.function.args[i];
        has_kernel = rt_find_exp_of_x(gi, x) || rt_find_log_of_x(gi, x) || rt_has_tan_of_x(gi, x)
                  || rt_has_prim_of_x(gi, x);
    }
    has_kernel = has_kernel || rt_find_exp_of_x(f, x) || rt_find_log_of_x(f, x) || rt_has_tan_of_x(f, x)
        || rt_has_prim_of_x(f, x)
        || (extra && (rt_find_exp_of_x(extra, x) || rt_find_log_of_x(extra, x) || rt_has_tan_of_x(extra, x)
                      || rt_has_prim_of_x(extra, x)));
    RdeCertState saved = rde_cert_begin();
    rde_cert_set_params(allp);
    Expr *y = NULL, *sol = NULL;
    int st;
    if (has_kernel) {
        st = rde_solve_param_ext(f, g, x, (extra && extra->data.function.arg_count) ? extra : NULL,
                                 unk, &y, &sol);
    } else {
        Expr* xl = expr_new_function(expr_new_symbol("List"), (Expr*[]){ expr_copy(x) }, 1);
        if (rde_rational_in(f, xl) && rde_rational_in(g, xl)) {
            RdeCtx C = { NULL, -1, x, x };
            st = rde_tower_param(f, g, &C, unk, &y, &sol);
        } else { rde_cert_mark("base: coefficients not rational in x"); st = -1; }
        expr_free(xl);
    }
    rde_cert_set_params(NULL);
    g_rde_cert_last_reason = st == 1 ? NULL : (st < 0 ? rde_cert_reason() : NULL);
    rde_cert_end(saved);
    expr_free(g); expr_free(allp); expr_free(unk);
    if (st == 1 && y && sol)
        return expr_new_function(expr_new_symbol("List"), (Expr*[]){ y, sol }, 2);
    if (y) expr_free(y);
    if (sol) expr_free(sol);
    return st == 0 ? expr_new_string("NoSolution") : expr_new_symbol("$Failed");
}

/* Risch`RischDECertifiedReason[] — the first uncertain path of the last
 * Risch`RischDECertified call that returned $Failed (a String), or None. */
static Expr* builtin_risch_rischde_certified_reason(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 0) return NULL;
    return g_rde_cert_last_reason ? expr_new_string(g_rde_cert_last_reason)
                                  : expr_new_symbol("None");
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
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    if (doc) symtab_set_docstring(name, doc);
}

void integrate_risch_transcendental_init(void) {
    /* flags that bracket a sub-solve: restored when a TimeConstrained timeout
     * unwinds out of it (otherwise a stuck decide mode degrades every later
     * integral of the session) */
    tc_register_unwind_counter(&g_rt_decide_mode);
    tc_register_unwind_counter(&g_rt_elementary_only);
    tc_register_unwind_counter(&rt_in_flat_gate);
    tc_register_unwind_counter(&rt_in_logcombine);
    tc_register_unwind_counter(rde_alg_mode_ptr());   /* the base RDE's algebraic mode */
    tc_register_unwind_counter(rde_par_mode_ptr());   /* the parametric RDE (step 3a) */
    rt_install("Integrate`RischTranscendental", builtin_rischtranscendental,
        "Integrate`RischTranscendental[f, x] integrates f with respect to x using the\n"
        "recursive transcendental Risch algorithm: a decision\n"
        "procedure over a differential transcendental tower, with rational,\n"
        "logarithmic, exponential, and special-function (Erf, ExpIntegralEi,\n"
        "LogIntegral, PolyLog) cases.  Each case is correct by construction (no\n"
        "differentiation check).  Out-of-scope integrands (algebraic\n"
        "extensions, non-elementary answers) return unevaluated.");
    rt_install("Risch`RischDE", builtin_risch_rischde,
        "Risch`RischDE[f, g, x] solves the Risch differential equation\n"
        "D[y] + f y == g for a rational y in C(x), or in the transcendental TOWER\n"
        "C(x)(t_1..t_n) inferred from the Log/Exp kernels of f and g (Bronstein\n"
        "Chapter 6, recursive over the tower): normal-denominator reduction, SPDE,\n"
        "and the polynomial non-cancellation solve, with the derivation lifted to\n"
        "the tower.  Returns y (0 when g is 0), or stays unevaluated when no\n"
        "rational solution exists.  Coefficients may be Gaussian (C(i)(x)).");
    rt_install("Risch`RischDECertified", builtin_risch_rischde_certified,
        "Risch`RischDECertified[f, g, x] solves D[y] + f y == g like Risch`RischDE,\n"
        "but distinguishes a proof from a give-up: it returns a solution y, the\n"
        "string \"NoSolution\" when no solution exists in C(x) or in the tower of the\n"
        "Log/Exp kernels of f and g (every decision taken is a theorem of Bronstein\n"
        "Chapter 6), or $Failed when undecided (an unimplemented sub-case or input\n"
        "outside rational functions of x and the kernels).\n"
        "Risch`RischDECertified[f, g, x, {k1, ...}] adds the kernels k_i to the tower,\n"
        "so the solution field is the caller's whole tower.\n"
        "Risch`RischDECertified[f, g, x, {k1, ...}, {a1, ...}] treats the symbols a_i as\n"
        "transcendental constants (D = 0): the verdict is over Q(a1, ...).");
    rt_install("Risch`ParametricRischDE", builtin_risch_parametric_rischde,
        "Risch`ParametricRischDE[f, {g0, g1, ..., gm}, x] solves the parametric Risch\n"
        "differential equation D[y] + f y == g0 + c1 g1 + ... + cm gm for y and\n"
        "constants c_i together (Bronstein, Symbolic Integration I, ch. 7), in C(x)\n"
        "or the tower of the Log/Exp/Tan kernels of the inputs.  Returns\n"
        "{y, {c1, ..., cm}} (a free constant is 0), \"NoSolution\" when no constants\n"
        "admit a solution, or $Failed when undecided (see\n"
        "Risch`RischDECertifiedReason).  Optional 4th and 5th arguments: extra tower\n"
        "kernels {k1, ...} and symbolic constants {a1, ...}, as for\n"
        "Risch`RischDECertified.");
    rt_install("Risch`RischDECertifiedReason", builtin_risch_rischde_certified_reason,
        "Risch`RischDECertifiedReason[] names the first undecided path of the last\n"
        "Risch`RischDECertified call that returned $Failed (a String), or None.");
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
