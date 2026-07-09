/* test_intrischnorman.c — Risch-Norman heuristic integrator tests.
 *
 * The corpus grows phase-by-phase per plans/RISCH_NORMAN_PLAN.md:
 *   Phase 1: scaffolding / dispatcher integration.
 *   Phase 2: convert_to_tan + indet collection unit tests.
 *   Phase 3: vector field + splitFactor + deflation unit tests.
 *   Phase 4: 10 Q-rational integrands (Exp[x], Log[x], x Exp[x], ...).
 *   Phase 5: 15+ log-bearing integrands + K=I retry cases.
 *   Phase 6: full ~60-integrand corpus.
 *
 * Universal correctness predicate (Phase 4+):
 *   ASSERT_INTEGRAL_OK(f, x)  ≡  Cancel[Together[Expand[
 *       D[Integrate`RischNorman[f, x], x] - f]]] === 0
 * Mirrors the helper in tests/test_intrat.c:56-72.
 */

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Test helpers.                                                       */
/* ------------------------------------------------------------------ */

/* Compare evaluated output to a string. */
static void run_eq(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n",
               input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Check that Integrate`RischNorman[f, x] differentiates back to f.
 *
 * Pmint's result lives in Tan[x/2] / Tanh[x/2] form.  We reconcile by:
 *   (1) Applying x-specific Weierstrass rules to BOTH integrand and
 *       result.  The rules are keyed on Tan[x] / Sin[x] / etc. (not
 *       Tan[x/2]) so the result's Tan[x/2] stays intact while the
 *       integrand's Tan[x] is rewritten.
 *   (2) Differentiating the result.
 *   (3) Applying PMPythagoreanRewrite (Sec^2 → 1+Tan^2 etc.).
 *   (4) Cancel + Together + Expand.
 *
 * NOTE: this helper deliberately does NOT call TrigExpand on the
 * integrand — TrigExpand blows up products like Sin[x]^N and would
 * regress the Phase 6 Sin[x]^2 / Cos[x]^2 / Tan[x]^2 tests.  For
 * integrands with compound trig arguments (Sin[2 x], Cos[a x + b],
 * Sin[x^2], etc.) use assert_rischnorman_correct_trigexpand. */
__attribute__((unused))
static void assert_rischnorman_correct(const char* integrand) {
    char buf[4096];
    snprintf(buf, sizeof(buf),
        "Module[{specific, generic, f, r, lhs},"
        "  specific = {Tan[x] -> 2 Tan[x/2]/(1 - Tan[x/2]^2),"
        "              Cot[x] -> (1 - Tan[x/2]^2)/(2 Tan[x/2]),"
        "              Sin[x] -> 2 Tan[x/2]/(1 + Tan[x/2]^2),"
        "              Cos[x] -> (1 - Tan[x/2]^2)/(1 + Tan[x/2]^2),"
        "              Sec[x] -> (1 + Tan[x/2]^2)/(1 - Tan[x/2]^2),"
        "              Csc[x] -> (1 + Tan[x/2]^2)/(2 Tan[x/2]),"
        "              Tanh[x] -> 2 Tanh[x/2]/(1 + Tanh[x/2]^2),"
        "              Coth[x] -> (1 + Tanh[x/2]^2)/(2 Tanh[x/2]),"
        "              Sinh[x] -> 2 Tanh[x/2]/(1 - Tanh[x/2]^2),"
        "              Cosh[x] -> (1 + Tanh[x/2]^2)/(1 - Tanh[x/2]^2),"
        "              Sech[x] -> (1 - Tanh[x/2]^2)/(1 + Tanh[x/2]^2),"
        "              Csch[x] -> (1 - Tanh[x/2]^2)/(2 Tanh[x/2])};"
        "  generic = {Sec[u_] Csc[u_] -> (1 + Tan[u]^2)/Tan[u],"
        "             Sech[u_] Csch[u_] -> (1 - Tanh[u]^2)/Tanh[u],"
        "             Cot[u_] -> 1/Tan[u],"
        "             Coth[u_] -> 1/Tanh[u]};"
        "  f = (%s) /. specific;"
        "  r = (Integrate`RischNorman[%s, x]) /. specific;"
        "  lhs = Integrate`Helpers`PMPythagoreanRewrite[D[r, x]] /. generic;"
        "  Cancel[Together[Expand[lhs - f]]]]",
        integrand, integrand);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, "0") != 0) {
        printf("FAIL: D[Integrate`RischNorman[%s, x], x] - %s != 0\n  Got: %s\n",
               integrand, integrand, got);
        ASSERT_STR_EQ(got, "0");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Variant of assert_rischnorman_correct for integrands containing
 * compound trig arguments (Sin[2 x], Cos[2 x], Tan[a x + b], ...).
 * Pre-applies TrigExpand to break those down into Sin[x] / Cos[x]
 * products that the literal Weierstrass rules below can rewrite.
 *
 * Kept separate from assert_rischnorman_correct because TrigExpand
 * blows up powers like Sin[x]^N and would regress the simpler Phase 6
 * tests on Sin[x]^2 / Cos[x]^2 / Tan[x]^2. */
__attribute__((unused))
static void assert_rischnorman_correct_trigexpand(const char* integrand) {
    char buf[4096];
    snprintf(buf, sizeof(buf),
        "Module[{specific, generic, f, r, lhs},"
        "  specific = {Tan[x] -> 2 Tan[x/2]/(1 - Tan[x/2]^2),"
        "              Cot[x] -> (1 - Tan[x/2]^2)/(2 Tan[x/2]),"
        "              Sin[x] -> 2 Tan[x/2]/(1 + Tan[x/2]^2),"
        "              Cos[x] -> (1 - Tan[x/2]^2)/(1 + Tan[x/2]^2),"
        "              Sec[x] -> (1 + Tan[x/2]^2)/(1 - Tan[x/2]^2),"
        "              Csc[x] -> (1 + Tan[x/2]^2)/(2 Tan[x/2]),"
        "              Tanh[x] -> 2 Tanh[x/2]/(1 + Tanh[x/2]^2),"
        "              Coth[x] -> (1 + Tanh[x/2]^2)/(2 Tanh[x/2]),"
        "              Sinh[x] -> 2 Tanh[x/2]/(1 - Tanh[x/2]^2),"
        "              Cosh[x] -> (1 + Tanh[x/2]^2)/(1 - Tanh[x/2]^2),"
        "              Sech[x] -> (1 - Tanh[x/2]^2)/(1 + Tanh[x/2]^2),"
        "              Csch[x] -> (1 - Tanh[x/2]^2)/(2 Tanh[x/2])};"
        "  generic = {Sec[u_] Csc[u_] -> (1 + Tan[u]^2)/Tan[u],"
        "             Sech[u_] Csch[u_] -> (1 - Tanh[u]^2)/Tanh[u],"
        "             Cot[u_] -> 1/Tan[u],"
        "             Coth[u_] -> 1/Tanh[u]};"
        "  f = TrigExpand[%s] /. specific;"
        "  r = TrigExpand[Integrate`RischNorman[%s, x]] /. specific;"
        "  lhs = TrigExpand[Integrate`Helpers`PMPythagoreanRewrite[D[r, x]]] /. generic;"
        "  Cancel[Together[Expand[lhs - f]]]]",
        integrand, integrand);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, "0") != 0) {
        printf("FAIL (trigexpand): D[Integrate`RischNorman[%s, x], x] - %s != 0\n  Got: %s\n",
               integrand, integrand, got);
        ASSERT_STR_EQ(got, "0");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Numerical correctness check for integrands where the symbolic
 * simplifier (Together / Cancel) can't collapse the residual — typically
 * because Mathilda's Together doesn't factor through transcendental
 * generators (e.g. it treats `x^2 Log[x] - x E^(1+x+x^2)` as opaque
 * rather than `x*(x Log[x] - E^(1+x+x^2))`).
 *
 * We sample D[r, x] - f at four rational points (x = 5/2, 7/2, 11/4, 3)
 * via N[] and require each to be effectively zero (|.| < 1e-9).  Four
 * non-special-value samples make a false-pass on a non-trivial residual
 * astronomically unlikely. */
__attribute__((unused))
static void assert_rischnorman_correct_numeric(const char* integrand) {
    static const char* probes[] = {"5/2", "7/2", "11/4", "3"};
    char buf[3072];
    /* Compute the integral once, then evaluate the residual at each probe.
     * Accept Real residuals with |.| < 1e-9 OR Complex residuals where both
     * components are < 1e-9 — Mathilda's D often produces a Complex form
     * with a tiny imaginary part when Sec/Csc are involved. */
    for (size_t i = 0; i < sizeof(probes) / sizeof(probes[0]); ++i) {
        snprintf(buf, sizeof(buf),
            "Module[{r, diff, n}, r = Integrate`RischNorman[%s, x];"
            "  diff = D[r, x] - (%s);"
            "  n = N[diff /. x -> %s];"
            "  If[Head[n] === Real && Abs[n] < 10^-9, 0,"
            "    If[Head[n] === Complex && Abs[Re[n]] < 10^-9 && Abs[Im[n]] < 10^-9, 0, n]]]",
            integrand, integrand, probes[i]);
        Expr* e = parse_expression(buf);
        Expr* res = evaluate(e);
        char* got = expr_to_string(res);
        if (strcmp(got, "0") != 0) {
            printf("FAIL (numeric): D[Integrate`RischNorman[%s, x], x] - (%s) "
                   "!= 0 at x=%s.\n  Residual: %s\n",
                   integrand, integrand, probes[i], got);
            ASSERT_STR_EQ(got, "0");
        }
        free(got);
        expr_free(e);
        expr_free(res);
    }
}

/* Check that Integrate[f, x] (dispatcher path) differentiates back to f.
 * Used in Phase 6 for the end-to-end dispatcher tests. */
__attribute__((unused))
static void assert_integrate_correct(const char* integrand) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Cancel[Together[Expand[D[Integrate[%s, x], x] - (%s)]]]",
        integrand, integrand);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, "0") != 0) {
        printf("FAIL: D[Integrate[%s, x], x] - %s != 0\n  Got: %s\n",
               integrand, integrand, got);
        ASSERT_STR_EQ(got, "0");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ */
/* Phase 1 — scaffolding tests.                                        */
/* ------------------------------------------------------------------ */

/* The Integrate`RischNorman symbol must be registered after
 * integrate_init() runs (via core_init()). */
static void test_phase1_symbol_registered(void) {
    SymbolDef* def = symtab_lookup("Integrate`RischNorman");
    ASSERT(def != NULL);
    ASSERT(def->builtin_func != NULL);
    /* Must be protected so the user can't redefine it accidentally. */
    ASSERT((def->attributes & ATTR_PROTECTED) != 0);
}

/* Genuinely non-elementary integrands bubble back as the unevaluated
 * package head. */
static void test_phase1_unevaluated_returns_self(void) {
    run_eq("Integrate`RischNorman[1/Log[x], x]",
           "Integrate`RischNorman[1/Log[x], x]");
    /* Exp[x^2] is non-elementary (related to the error function). */
    run_eq("Integrate`RischNorman[Exp[x^2], x]",
           "Integrate`RischNorman[E^x^2, x]");
}

/* Regression: an inexact-base exponential integrand (e.g. 2.71828^(-x), as
 * produced by N[] over an Exp[] integrand) once drove split_factor into
 * unbounded recursion and overflowed the C stack, killing the whole process.
 * The recursion is now depth-guarded (PMINT_MAX_SPLIT_DEPTH): the integrator
 * must give up gracefully and leave the integral unevaluated. */
static void test_phase1_inexact_base_no_stack_overflow(void) {
    run_eq("Head[Integrate[x*2.71828^(-x), x]]",         "Integrate");
    run_eq("Head[Integrate[x^2*2.71828^(-x), x]]",       "Integrate");
    run_eq("Head[N[Integrate[x*Exp[-x], {x, 0, Infinity}]]]", "Integrate");
}

/* Wrong arity must also bubble back unevaluated (the BuiltinFunc
 * returns NULL early). */
static void test_phase1_wrong_arity_unevaluated(void) {
    run_eq("Integrate`RischNorman[x]",
           "Integrate`RischNorman[x]");
    run_eq("Integrate`RischNorman[x, y, z]",
           "Integrate`RischNorman[x, y, z]");
}

/* Non-symbol second argument bubbles back unevaluated. */
static void test_phase1_non_symbol_var(void) {
    run_eq("Integrate`RischNorman[Sin[x], 5]",
           "Integrate`RischNorman[Sin[x], 5]");
}

/* The Integrate dispatcher must still produce the same behaviour
 * for rational integrands (Phase 1 only changes the non-rational
 * fall-through). */
static void test_phase1_rational_unchanged(void) {
    run_eq("Integrate[x, x]", "1/2 x^2");
    run_eq("Integrate[1/(x^2 + 1), x]", "ArcTan[x]");
}

/* Integrate now routes through RischNorman for non-rational
 * integrands.  Genuinely non-elementary integrands (1/Log[x],
 * Exp[x^2]) still bubble back unevaluated. */
static void test_phase1_known_fail_unevaluated(void) {
    run_eq("Integrate[1/Log[x], x]", "Integrate[1/Log[x], x]");
    run_eq("Integrate[Exp[x^2], x]", "Integrate[E^x^2, x]");
}

/* ------------------------------------------------------------------ */
/* Phase 2 — convert_to_tan / collect_indets / subst_map.              */
/* ------------------------------------------------------------------ */

/* PMConvertToTan rewrites Sin / Cos / Sec / Csc to half-angle Tan
 * rationals; Sinh / Cosh / Sech / Csch via Tanh[u/2].  Tan, Cot,
 * Tanh, Coth are reciprocally normalised but the Tan / Tanh atoms
 * remain as field generators. */
static void test_phase2_convert_sin(void) {
    run_eq("Integrate`Helpers`PMConvertToTan[Sin[x], x]",
           "(2 Tan[1/2 x])/(1 + Tan[1/2 x]^2)");
}

static void test_phase2_convert_cos(void) {
    run_eq("Integrate`Helpers`PMConvertToTan[Cos[x], x]",
           "(1 - Tan[1/2 x]^2)/(1 + Tan[1/2 x]^2)");
}

static void test_phase2_convert_sec(void) {
    run_eq("Integrate`Helpers`PMConvertToTan[Sec[x], x]",
           "(1 + Tan[1/2 x]^2)/(1 - Tan[1/2 x]^2)");
}

static void test_phase2_convert_csc(void) {
    /* Csc[u] = (1+T^2)/(2T) — Mathilda's evaluator collapses 1/T to Cot,
     * yielding the equivalent Cot form.  Algebraically the same; pmint
     * uses the decot'd internal tree, so this REPL-visible form is OK. */
    run_eq("Integrate`Helpers`PMConvertToTan[Csc[x], x]",
           "1/2 Cot[1/2 x] (1 + Tan[1/2 x]^2)");
}

static void test_phase2_convert_cot(void) {
    /* Cot[u] = (1 - T^2) / (2T),  T = Tan[u/2].  Mathilda re-collapses
     * the 1/T factor to Cot[u/2]; the externally visible form is the
     * algebraically equivalent 1/2 Cot[1/2 x] (1 - Tan[1/2 x]^2). */
    run_eq("Integrate`Helpers`PMConvertToTan[Cot[x], x]",
           "1/2 Cot[1/2 x] (1 - Tan[1/2 x]^2)");
}

static void test_phase2_convert_tan_full_weierstrass(void) {
    /* Tan[u] = 2T / (1 - T^2),  T = Tan[u/2].  Full Weierstrass keeps
     * Tan[u/2] as the unique trig generator. */
    run_eq("Integrate`Helpers`PMConvertToTan[Tan[x], x]",
           "(2 Tan[1/2 x])/(1 - Tan[1/2 x]^2)");
}

static void test_phase2_convert_sinh(void) {
    /* Sinh[u] = 2 T / (1 - T^2),  T = Tanh[u/2]. */
    run_eq("Integrate`Helpers`PMConvertToTan[Sinh[x], x]",
           "(2 Tanh[1/2 x])/(1 - Tanh[1/2 x]^2)");
}

static void test_phase2_convert_cosh(void) {
    run_eq("Integrate`Helpers`PMConvertToTan[Cosh[x], x]",
           "(1 + Tanh[1/2 x]^2)/(1 - Tanh[1/2 x]^2)");
}

static void test_phase2_convert_skips_free_of_x(void) {
    /* Sin[y] doesn't mention x — should remain unchanged. */
    run_eq("Integrate`Helpers`PMConvertToTan[Sin[y], x]", "Sin[y]");
}

static void test_phase2_convert_exp_unchanged(void) {
    /* Exp / E^x is not in the trig family — left alone. */
    run_eq("Integrate`Helpers`PMConvertToTan[Exp[x], x]", "E^x");
}

/* PMCollectIndets returns the list of transcendental atoms in f
 * (assuming f has already been run through PMConvertToTan) with
 * derivative wrt x non-zero, closed under one round of diff.  The
 * first element is always x itself. */
static void test_phase2_indets_x_and_exp(void) {
    /* Exp[x] has D[E^x, x] = E^x — same atom, closure is a no-op. */
    run_eq("Integrate`Helpers`PMCollectIndets[Exp[x], x]",
           "{x, E^x}");
}

static void test_phase2_indets_x_and_log(void) {
    /* Log[x] derivative is 1/x — no new atoms after closure. */
    run_eq("Integrate`Helpers`PMCollectIndets[Log[x], x]",
           "{x, Log[x]}");
}

static void test_phase2_indets_x_tan(void) {
    /* Tan[x] survives convert_to_tan; D[Tan[x],x] = Sec[x]^2 but we
     * don't recurse into the Sec head (Sec isn't a collected atom). */
    run_eq("Integrate`Helpers`PMCollectIndets[Tan[x], x]",
           "{x, Tan[x]}");
}

static void test_phase2_indets_mixed(void) {
    run_eq("Integrate`Helpers`PMCollectIndets[Tan[x] + Log[x] + Exp[x], x]",
           "{x, E^x, Log[x], Tan[x]}");
}

static void test_phase2_indets_skips_constant_atoms(void) {
    /* Log[y] doesn't mention x — D[Log[y], x] = 0, not collected. */
    run_eq("Integrate`Helpers`PMCollectIndets[Log[y] + Exp[x], x]",
           "{x, E^x}");
}

/* PMSubstMap returns the forward substitution rules
 * `term -> pmint$v_i`.  When an atom is Tan[u] or Tanh[u] the map
 * also includes a companion `Cot[u] -> 1/v_k` / `Coth[u] -> 1/v_k`
 * rule so the reciprocal generator collapses through the same
 * fresh variable. */
static void test_phase2_subst_map_basic(void) {
    run_eq("Integrate`Helpers`PMSubstMap[Tan[x] + Exp[x], x]",
           "{x -> pmint$v1, E^x -> pmint$v2, Tan[x] -> pmint$v3, Cot[x] -> 1/pmint$v3}");
}

/* Round-trip: subs[lout, subs[lin, ff]] == ff.  Approximate the
 * round-trip by applying the forward substitution then a reverse
 * built by inverting it manually with Cases. */
static void test_phase2_subst_round_trip(void) {
    /* The subst map maps Tan[x] -> pmint$v_k.  Replacing back via
     * ReplaceAll should restore the original.  Build inv = Reverse
     * each rule. */
    run_eq("With[{m = Integrate`Helpers`PMSubstMap[Tan[x] + Log[x], x],"
           "      f = Tan[x] + Log[x]},"
           "  Module[{inv}, inv = (#[[2]] -> #[[1]]) & /@ m;"
           "    Expand[((f /. m) /. inv) - f]]]",
           "0");
}

/* ------------------------------------------------------------------ */
/* Phase 3 — vector field / split_factor / deflation / monomials.      */
/* ------------------------------------------------------------------ */

/* PMVectorField returns {vars, l, q}.  For integrand Exp[x] the field
 * has vars = {v1=x, v2=E^x}, l = {1, v2}, q = 1. */
static void test_phase3_vector_field_exp(void) {
    run_eq("Integrate`Helpers`PMVectorField[Exp[x], x]",
           "{{pmint$v1, pmint$v2}, {1, pmint$v2}, 1}");
}

/* For Log[x] the field has D[Log[x], x] = 1/x.  q = lcm of denominators
 * including the x denominator from Log's derivative. */
static void test_phase3_vector_field_log(void) {
    run_eq("Integrate`Helpers`PMVectorField[Log[x], x]",
           "{{pmint$v1, pmint$v2}, {pmint$v1, 1}, pmint$v1}");
}

/* apply_d should give back D[f, x] (up to the q scaling) for any f in
 * the integrand's differential field.  Since q = 1 for Exp[x], apply_d
 * of Exp[x] returns Exp[x] itself. */
static void test_phase3_apply_d_exp(void) {
    run_eq("Integrate`Helpers`PMApplyD[Exp[x], x]", "E^x");
}

/* apply_d on x in Exp[x]'s field: D[x, x] = 1; q = 1; l[x] = 1.
 * Result: 1. */
static void test_phase3_apply_d_x_exp(void) {
    run_eq("Integrate`Helpers`PMApplyD[x, x]", "1");
}

/* splitFactor[V (V+1)^2, dx] should return {V, (V+1)^2} where V = E^x
 * has d(V) = V.  The print form distributes (V+1)^2 = 1 + 2V + V^2. */
static void test_phase3_split_factor_v_times_vp1sq(void) {
    run_eq("Integrate`Helpers`PMSplitFactor[Exp[x] (Exp[x]+1)^2, x]",
           "{E^x, 1 + 2 E^x + E^(2 x)}");
}

/* deflation((V+1)^2 (V-1)) should give (V+1)(V-1)... actually pmint's
 * deflation returns just the d-gcd of the primitive part, which for
 * (V+1)^2 (V-1) is (V+1). */
static void test_phase3_deflation_basic(void) {
    run_eq("Integrate`Helpers`PMDeflation[(Exp[x]+1)^2 (Exp[x]-1), x]",
           "1 + E^x");
}

/* enumerate_monomials({x, V}, 3) should produce C(5, 2) = 10 monomials. */
static void test_phase3_enumerate_monoms_count(void) {
    run_eq("Length[Integrate`Helpers`PMEnumerateMonoms[{x, V}, 3]]", "10");
}

/* enumerate_monomials with zero vars returns {1}. */
static void test_phase3_enumerate_monoms_empty(void) {
    run_eq("Integrate`Helpers`PMEnumerateMonoms[{}, 5]", "{1}");
}

/* enumerate_monomials with a single var and degree 0 returns {1}. */
static void test_phase3_enumerate_monoms_deg_zero(void) {
    run_eq("Integrate`Helpers`PMEnumerateMonoms[{x}, 0]", "{1}");
}

/* enumerate_monomials with a single var and degree d returns d+1 items
 * {1, x, x^2, ..., x^d}. */
static void test_phase3_enumerate_monoms_single_var(void) {
    run_eq("Integrate`Helpers`PMEnumerateMonoms[{x}, 4]",
           "{1, x, x^2, x^3, x^4}");
}

/* ------------------------------------------------------------------ */
/* QMatBuild hash index unit tests.                                    */
/*                                                                      */
/* Exercises the open-addressed exp_vec → row hash that backs          */
/* try_solve_direct_q.  Each test calls the testable surface           */
/* Integrate`Helpers`PMQMBStress[nv, ncols, ntrials, mod, seed], which */
/* runs an internal stress loop and returns True/False.                */
/* ------------------------------------------------------------------ */

/* nv = 0 — degenerate single-vector case.  All inserts collapse to    */
/* one row regardless of ntrials.                                       */
static void test_qmb_hash_nv_zero(void) {
    run_eq("Integrate`Helpers`PMQMBStress[0, 1, 50, 1, 1]", "True");
}

/* nv = 1 — heavy collisions (mod = 5 so only 5 unique vectors but     */
/* 1000 trials).  Tests idempotency under repeated re-insertion.       */
static void test_qmb_hash_high_collision(void) {
    run_eq("Integrate`Helpers`PMQMBStress[1, 1, 1000, 5, 999]", "True");
}

/* Typical pmint shape — small nv, modest unique count, plenty of       */
/* trials to force at least 2-3 bucket grows (initial cap = 32, then   */
/* doubling).                                                            */
static void test_qmb_hash_typical(void) {
    run_eq("Integrate`Helpers`PMQMBStress[3, 2, 200, 4, 17]", "True");
}

/* Larger working set — 5000 trials with mod = 3, nv = 8.  Unique       */
/* vector count is bounded by 3^8 = 6561 so we exercise growth         */
/* through several rehashes.                                            */
static void test_qmb_hash_many_rehashes(void) {
    run_eq("Integrate`Helpers`PMQMBStress[8, 4, 5000, 3, 42]", "True");
}

/* Maximum supported nv (PMINT_MAX_INDETS = 32).  Verifies the hash    */
/* and memcmp handle the longest legal key.                             */
static void test_qmb_hash_max_indets(void) {
    run_eq("Integrate`Helpers`PMQMBStress[32, 1, 200, 2, 7]", "True");
}

/* ntrials = 0 — no inserts.  Verifies qmb_free copes with an empty    */
/* table (both rows and buckets still NULL).                            */
static void test_qmb_hash_empty(void) {
    run_eq("Integrate`Helpers`PMQMBStress[3, 1, 0, 1, 0]", "True");
}

/* Determinism — same seed must reproduce the same result, twice in    */
/* a row.  Guards against any iteration-order reliance.                 */
static void test_qmb_hash_deterministic(void) {
    run_eq("Integrate`Helpers`PMQMBStress[4, 2, 500, 3, 12345] && "
           "Integrate`Helpers`PMQMBStress[4, 2, 500, 3, 12345]",
           "True");
}

/* ------------------------------------------------------------------ */
/* Phase 4 — Q-rational integrals.  Uses the universal correctness    */
/* predicate (differentiate-and-cancel).                                */
/* ------------------------------------------------------------------ */

/* Phase 5 — log-bearing integrals. */
static void test_phase5_inv_x(void)              { assert_rischnorman_correct("1/x"); }
static void test_phase5_sin_x(void)              { assert_rischnorman_correct("Sin[x]"); }
static void test_phase5_cos_x(void)              { assert_rischnorman_correct("Cos[x]"); }
static void test_phase5_tan_x(void)              { assert_rischnorman_correct("Tan[x]"); }
static void test_phase5_cot_x(void)              { assert_rischnorman_correct("Cot[x]"); }
static void test_phase5_inv_1_plus_exp(void)     { assert_rischnorman_correct("1/(1 + Exp[x])"); }
static void test_phase5_exp_over_1_plus_exp(void){ assert_rischnorman_correct("Exp[x]/(1 + Exp[x])"); }
static void test_phase5_inv_x_log_x(void)        { assert_rischnorman_correct("1/(x Log[x])"); }
static void test_phase5_log_x_squared(void)      { assert_rischnorman_correct("Log[x]^2"); }
static void test_phase5_one_plus_log_over_x_log(void){ assert_rischnorman_correct("(1 + Log[x])/(x Log[x])"); }

/* Phase 6 — extended corpus (Geddes Table II + Bronstein + Davenport). */
static void test_phase6_exp_ax(void)             { assert_rischnorman_correct("Exp[a x]"); }
static void test_phase6_x3_exp_x(void)           { assert_rischnorman_correct("x^3 Exp[x]"); }
static void test_phase6_log_x_cubed(void)        { assert_rischnorman_correct("Log[x]^3"); }
static void test_phase6_sinh_x(void)             { assert_rischnorman_correct("Sinh[x]"); }
static void test_phase6_cosh_x(void)             { assert_rischnorman_correct("Cosh[x]"); }
static void test_phase6_tanh_x(void)             { assert_rischnorman_correct("Tanh[x]"); }
static void test_phase6_sin_squared(void)        { assert_rischnorman_correct("Sin[x]^2"); }
static void test_phase6_cos_squared(void)        { assert_rischnorman_correct("Cos[x]^2"); }
static void test_phase6_tan_squared(void)        { assert_rischnorman_correct("Tan[x]^2"); }
static void test_phase6_one_over_x2_plus_one(void) { assert_rischnorman_correct("1/(1+x^2)"); }

/* Known-fail integrals: pmint should bubble back unevaluated cleanly. */
static void assert_rischnorman_unevaluated(const char* integrand) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "Integrate`RischNorman[%s, x]", integrand);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    /* Result must remain a function whose head is
     * Integrate`RischNorman — i.e., the call did not produce a closed
     * form. */
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.head
           && res->data.function.head->type == EXPR_SYMBOL
           && strcmp(res->data.function.head->data.symbol,
                     "Integrate`RischNorman") == 0);
    expr_free(e);
    expr_free(res);
}

static void test_phase6_fail_inv_log(void)       { assert_rischnorman_unevaluated("1/Log[x]"); }
static void test_phase6_fail_exp_x2(void)        { assert_rischnorman_unevaluated("Exp[x^2]"); }
static void test_phase6_fail_log_cos(void)       { assert_rischnorman_unevaluated("Log[x] Cos[x]"); }
static void test_phase6_fail_exp_over_log(void)  { assert_rischnorman_unevaluated("Exp[x]/Log[x]"); }
static void test_phase6_fail_exp_sin(void)       { assert_rischnorman_unevaluated("Exp[Sin[x]]"); }

static void test_phase4_exp_x(void)    { assert_rischnorman_correct("Exp[x]"); }
static void test_phase4_exp_2x(void)   { assert_rischnorman_correct("Exp[2 x]"); }
static void test_phase4_x_exp_x(void)  { assert_rischnorman_correct("x Exp[x]"); }
static void test_phase4_x2_exp_x(void) { assert_rischnorman_correct("x^2 Exp[x]"); }
static void test_phase4_log_x(void)    { assert_rischnorman_correct("Log[x]"); }
static void test_phase4_x_log_x(void)  { assert_rischnorman_correct("x Log[x]"); }
static void test_phase4_sin_exp(void)  { assert_rischnorman_correct("Sin[x] Exp[x]"); }
static void test_phase4_cos_exp(void)  { assert_rischnorman_correct("Cos[x] Exp[x]"); }

/* ------------------------------------------------------------------ */
/* Phase 7 — Bronstein / Geddes pmint reference corpus.                */
/* Source: the canonical Pmint test list from Bronstein's Risch-Norman */
/* paper plus the Davenport/Geddes textbook examples.  All integrals   */
/* admit an elementary closed form; some residuals require the          */
/* numerical-fallback predicate because Mathilda's Together can't       */
/* factor through transcendental generators.                            */
/* ------------------------------------------------------------------ */

/* --- 7a: pure-rational integrands. ---------------------------------- */
static void test_phase7_pmint_paper_1(void) {
    assert_rischnorman_correct(
        "(x^7 - 24 x^4 - 4 x^2 + 8 x - 8)/(x^8 + 6 x^6 + 12 x^4 + 8 x^2)");
}
static void test_phase7_pmint_paper_2(void) {
    assert_rischnorman_correct(
        "(-4 x^2 - 4 x^3 - x^4)/((-1 + x^2) (1 + x + x^2)^2)");
}
static void test_phase7_x3_over_xp1(void) {
    assert_rischnorman_correct("x^3/(1 + x)");
}
static void test_phase7_rat_log_complex(void) {
    assert_rischnorman_correct(
        "(x^4 - 3 x^2 + 6)/(x^6 - 5 x^4 + 5 x^2 + 4)");
}

/* --- 7b: log-bearing integrands. ------------------------------------ */
static void test_phase7_log_cube(void) {
    assert_rischnorman_correct("Log[x]^3");
}
static void test_phase7_log_neg2(void) {
    assert_rischnorman_correct("(-1 + Log[x])/Log[x]^2");
}
static void test_phase7_exp_over_log_shifted(void) {
    assert_rischnorman_correct(
        "(E^x (-1 - Log[-1 + x] + x Log[-1 + x]))/((-1 + x) Log[-1 + x]^2)");
}
static void test_phase7_log_diff_sq(void) {
    assert_rischnorman_correct("(1 - Log[x])/(x^2 - Log[x]^2)");
}
static void test_phase7_sin_over_x_plus_log_cos(void) {
    assert_rischnorman_correct("Sin[x]/x + Log[x] Cos[x]");
}

/* --- 7c: exp-bearing integrands. ------------------------------------ */
static void test_phase7_exp_inv_x(void) {
    assert_rischnorman_correct("((x + 1)/x^4) Exp[1/x]");
}
static void test_phase7_inv_1_plus_exp(void) {
    assert_rischnorman_correct("1/(1 + Exp[x])");
}
static void test_phase7_exp_x_plus_inv_log(void) {
    assert_rischnorman_correct("(1 - 1/(x Log[x]^2)) Exp[1/Log[x] + x]");
}
static void test_phase7_exp_rat(void) {
    assert_rischnorman_correct(
        "((-1 - x - x^2 + x^3)/(1 - 2 x^2 + x^4)) Exp[x]");
}
static void test_phase7_exp_inv_log(void) {
    assert_rischnorman_correct("(Exp[1/Log[x]] (Log[x]^2 - 1))/Log[x]^2");
}
static void test_phase7_exp_quad_denom(void) {
    assert_rischnorman_correct(
        "(E^x (-1 + 1289 x + 278 x^2 + 55 x^3 + 56 x^4))/(-392 + x - 56 x^2)^2");
}
static void test_phase7_exp_x2_diff(void) {
    assert_rischnorman_correct(
        "((4 x^2 + 4 x - 1) (Exp[x^2] + 1) (Exp[x^2] - 1))/(x + 1)^2");
}
static void test_phase7_x3_exp_x2(void) {
    assert_rischnorman_correct("x^3 Exp[x^2]");
}
static void test_phase7_exp_inv_x_x3(void) {
    assert_rischnorman_correct("Exp[1/x]/x^3");
}
static void test_phase7_x2_minus_x_exp3x(void) {
    assert_rischnorman_correct("(x^2 - x) Exp[3 x]");
}

/* --- 7d: trig+polynomial integrands. -------------------------------- */
static void test_phase7_trig_x2_polynomial(void) {
    assert_rischnorman_correct("(x^3 + 2 x) Cos[x^2] + x Sin[x^2]");
}
static void test_phase7_x3_sin_x2_plus_1(void) {
    assert_rischnorman_correct("x^3 Sin[x^2 + 1]");
}
static void test_phase7_arctan_over_x2(void) {
    assert_rischnorman_correct("ArcTan[x]/x^2");
}

/* --- 7e: trig integrands.  Numerical verification keeps the test    */
/* fast — the symbolic-Weierstrass route works (see                    */
/* assert_rischnorman_correct_trigexpand) but TrigExpand applied to    */
/* the full integral result blows up in expression size and is slow.   */
static void test_phase7_tan_div(void) {
    assert_rischnorman_correct_numeric(
        "(x - Tan[x])/Tan[x]^2 + Tan[x]");
}
static void test_phase7_cot_half(void) {
    assert_rischnorman_correct_numeric(
        "(1 - Cos[x] + Sin[x])/(-1 + Cos[x])^2");
}
static void test_phase7_davenport_exp_sin2x(void) {
    assert_rischnorman_correct_numeric("Exp[x] Sin[2 x]");
}
static void test_phase7_sin_sin_2x(void) {
    assert_rischnorman_correct_numeric("Sin[x] Sin[2 x]");
}

/* --- 7f: integrals that pmint integrates correctly but whose         */
/* residual Mathilda's Together can't collapse symbolically (because    */
/* it doesn't factor through transcendental generators).  Verified    */
/* numerically instead. ------------------------------------------------ */
static void test_phase7_num_quadratic_exp(void) {
    /* Result: -Log[x] + Log[-E^(1+x+x^2) + x Log[x]].
     * The residual contains denominators x, (x L - E), and
     * (x^2 L - x E) = x (x L - E), but Mathilda's Together treats the
     * third denominator as opaque and fails to see the common factor. */
    assert_rischnorman_correct_numeric(
        "(x + Exp[x^2 + x + 1] (1 - x - 2 x^2))"
        "/(x (x Log[x] - Exp[x^2 + x + 1]))");
}
static void test_phase7_num_two_exps(void) {
    /* Result: E^(-1+x^2)/(-1+E^x).  Together fails to identify
     * E^(x+x^2) with E^x * E^(x^2) under its current normalisation. */
    assert_rischnorman_correct_numeric(
        "(E^(-1 + x^2) (-2 x + E^x (-1 + 2 x)))/(-1 + E^x)^2");
}

/* --- 7g: constant-of-integration stripping. -------------------------- */
/* An antiderivative is defined only up to an additive constant, so the
 * result must never carry a stray term that is free of x.  pmint's
 * parametric solve used to leave one — e.g. x Sin[x^2] integrated to
 * 1/2 (-1 - Cos[x^2]) with a spurious -1/2.  We assert both that the
 * derivative checks out AND that Expand[result] has no x-free summand. */
static void assert_no_constant_of_integration(const char* integrand) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Module[{r, e, terms},"
        "  r = Integrate`RischNorman[%s, x];"
        "  e = Expand[r];"
        "  terms = If[Head[e] === Plus, List @@ e, {e}];"
        "  Select[terms, FreeQ[#, x] &]]",
        integrand);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, "{}") != 0) {
        printf("FAIL: Integrate`RischNorman[%s, x] carries a constant of "
               "integration.\n  x-free summands: %s\n", integrand, got);
        ASSERT_STR_EQ(got, "{}");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

static void test_phase7_no_const_x_sin_x2(void) {
    assert_rischnorman_correct("x Sin[x^2]");
    assert_no_constant_of_integration("x Sin[x^2]");
}
static void test_phase7_no_const_x_cos_x2(void) {
    assert_rischnorman_correct("x Cos[x^2]");
    assert_no_constant_of_integration("x Cos[x^2]");
}

/* ------------------------------------------------------------------ */
/* main().                                                              */
/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    /* Phase 1. */
    TEST(test_phase1_symbol_registered);
    TEST(test_phase1_unevaluated_returns_self);
    TEST(test_phase1_inexact_base_no_stack_overflow);
    TEST(test_phase1_wrong_arity_unevaluated);
    TEST(test_phase1_non_symbol_var);
    TEST(test_phase1_rational_unchanged);
    TEST(test_phase1_known_fail_unevaluated);

    /* Phase 2. */
    TEST(test_phase2_convert_sin);
    TEST(test_phase2_convert_cos);
    TEST(test_phase2_convert_sec);
    TEST(test_phase2_convert_csc);
    TEST(test_phase2_convert_cot);
    TEST(test_phase2_convert_tan_full_weierstrass);
    TEST(test_phase2_convert_sinh);
    TEST(test_phase2_convert_cosh);
    TEST(test_phase2_convert_skips_free_of_x);
    TEST(test_phase2_convert_exp_unchanged);
    TEST(test_phase2_indets_x_and_exp);
    TEST(test_phase2_indets_x_and_log);
    TEST(test_phase2_indets_x_tan);
    TEST(test_phase2_indets_mixed);
    TEST(test_phase2_indets_skips_constant_atoms);
    TEST(test_phase2_subst_map_basic);
    TEST(test_phase2_subst_round_trip);

    /* Phase 3. */
    TEST(test_phase3_vector_field_exp);
    TEST(test_phase3_vector_field_log);
    TEST(test_phase3_apply_d_exp);
    TEST(test_phase3_apply_d_x_exp);
    TEST(test_phase3_split_factor_v_times_vp1sq);
    TEST(test_phase3_deflation_basic);
    TEST(test_phase3_enumerate_monoms_count);
    TEST(test_phase3_enumerate_monoms_empty);
    TEST(test_phase3_enumerate_monoms_deg_zero);
    TEST(test_phase3_enumerate_monoms_single_var);

    /* QMatBuild hash index (backs try_solve_direct_q). */
    TEST(test_qmb_hash_nv_zero);
    TEST(test_qmb_hash_high_collision);
    TEST(test_qmb_hash_typical);
    TEST(test_qmb_hash_many_rehashes);
    TEST(test_qmb_hash_max_indets);
    TEST(test_qmb_hash_empty);
    TEST(test_qmb_hash_deterministic);

    /* Phase 4 — Q-rational integrals via assert_rischnorman_correct. */
    TEST(test_phase4_exp_x);
    TEST(test_phase4_exp_2x);
    TEST(test_phase4_x_exp_x);
    TEST(test_phase4_x2_exp_x);
    TEST(test_phase4_log_x);
    TEST(test_phase4_x_log_x);
    TEST(test_phase4_sin_exp);
    TEST(test_phase4_cos_exp);

    /* Phase 5 — log candidates + getSpecial. */
    TEST(test_phase5_inv_x);
    TEST(test_phase5_sin_x);
    TEST(test_phase5_cos_x);
    TEST(test_phase5_tan_x);
    TEST(test_phase5_cot_x);
    TEST(test_phase5_inv_1_plus_exp);
    TEST(test_phase5_exp_over_1_plus_exp);
    TEST(test_phase5_inv_x_log_x);
    TEST(test_phase5_log_x_squared);
    TEST(test_phase5_one_plus_log_over_x_log);

    /* Phase 6 — extended corpus. */
    TEST(test_phase6_exp_ax);
    TEST(test_phase6_x3_exp_x);
    TEST(test_phase6_log_x_cubed);
    TEST(test_phase6_sinh_x);
    TEST(test_phase6_cosh_x);
    TEST(test_phase6_tanh_x);
    TEST(test_phase6_sin_squared);
    TEST(test_phase6_cos_squared);
    TEST(test_phase6_tan_squared);
    TEST(test_phase6_one_over_x2_plus_one);
    TEST(test_phase6_fail_inv_log);
    TEST(test_phase6_fail_exp_x2);
    TEST(test_phase6_fail_log_cos);
    TEST(test_phase6_fail_exp_over_log);
    TEST(test_phase6_fail_exp_sin);

    /* Phase 7 — Bronstein / Geddes pmint reference corpus. */
    TEST(test_phase7_pmint_paper_1);
    TEST(test_phase7_pmint_paper_2);
    TEST(test_phase7_x3_over_xp1);
    TEST(test_phase7_rat_log_complex);
    TEST(test_phase7_log_cube);
    TEST(test_phase7_log_neg2);
    TEST(test_phase7_exp_over_log_shifted);
    TEST(test_phase7_log_diff_sq);
    TEST(test_phase7_sin_over_x_plus_log_cos);
    TEST(test_phase7_exp_inv_x);
    TEST(test_phase7_inv_1_plus_exp);
    TEST(test_phase7_exp_x_plus_inv_log);
    TEST(test_phase7_exp_rat);
    TEST(test_phase7_exp_inv_log);
    TEST(test_phase7_exp_quad_denom);
    TEST(test_phase7_exp_x2_diff);
    TEST(test_phase7_x3_exp_x2);
    TEST(test_phase7_exp_inv_x_x3);
    TEST(test_phase7_x2_minus_x_exp3x);
    TEST(test_phase7_trig_x2_polynomial);
    TEST(test_phase7_x3_sin_x2_plus_1);
    TEST(test_phase7_arctan_over_x2);
    TEST(test_phase7_tan_div);
    TEST(test_phase7_cot_half);
    TEST(test_phase7_davenport_exp_sin2x);
    TEST(test_phase7_sin_sin_2x);
    TEST(test_phase7_num_quadratic_exp);
    TEST(test_phase7_num_two_exps);
    TEST(test_phase7_no_const_x_sin_x2);
    TEST(test_phase7_no_const_x_cos_x2);

    printf("All intrischnorman tests passed!\n");
    return 0;
}
