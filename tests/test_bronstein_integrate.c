/* test_bronstein_integrate.c -- the Bronstein pipeline end to end
 * (src/internal/bronstein/Pipeline.m) and its C method glue
 * (Integrate[..., Method -> "Bronstein"], Integrate`Bronstein), BRONSTEIN_PLAN.md
 * Phase 7.
 *
 * Closes: the paper's examples ([B] p.134, p.147 in its theta-form, p.158),
 * the arcsinh family (residues at the special places), Shape A families over
 * Log and Exp, the Moebius flattening, and the algebraic / transcendental
 * hand-offs; every answer is checked by PossibleZeroQ[D[F, x] - f].
 * Certificates: [B] p.139 (resultant criterion), Sqrt[Log x] (algebraic part),
 * lifted Nagell-Lutz (infinite-order residue divisor by specialization).
 * Implementation incomplete: nested radicals and a non-cyclic compositum (outside
 * one simple radical extension); trigonometric kernels and a further monomial
 * over a radical-argument kernel (inside it, not implemented yet).
 * Message routing: Integrate::nonelem fires once, Quiet silences it, Check
 * sees it.  Property: derivatives of random elementary Shape A functions are
 * integrated or declined, never answered wrongly.
 */
#include "bronstein_test_utils.h"
#include "symtab.h"

static const char* const closes[] = {
    "(x + 1)/((x Log[x] + x^2) Sqrt[x + Log[x]])",                                  /* [B] p.134 */
    "(((3 x + 1) Log[x] - x^3 + x^2) Sqrt[x + Log[x]] - (2 x^2 - x - 1) Log[x] - 2 x^3 + x^2 + x)/"
    "(x Log[x]^2 - (x^3 - 2 x^2) Log[x] - x^4 + x^3)",                               /* [B] p.147 */
    "(3 (x + E^x)^(1/3) + (2 x^2 + 3 x) E^x + 5 x^2)/(x (x + E^x)^(1/3))",          /* [B] p.158 */
    "E^x/Sqrt[E^(2 x) + 1]",
    "1/Sqrt[E^x + 1]",
    "Sqrt[E^x + 1]",
    "E^x Sqrt[1 + E^x]",
    "1/(x Sqrt[Log[x]^2 + 1])",
    "Sqrt[Log[x]^2 + 1]/(x Log[x])",
    "1/(x Sqrt[Log[x]])",
    "Log[x]/(x Sqrt[1 + Log[x]])",
    "1/(x (Log[x] + 2) Sqrt[Log[x]^2 + 2])",
    "D[Log[Sqrt[Log[x]^3 + x] + x], x]",
    "D[x Log[x] Sqrt[x + Log[x]], x]",
    "D[(x + E^x)^(2/3) E^x, x]",
    "Log[x]/Sqrt[x + 1]",                                                            /* Moebius flattening */
    "1/(x Sqrt[1 - x^3])",                                                           /* Trager hand-off */
    "x E^x + Log[x]",                                                                /* Risch hand-off */
    /* stress regressions (2026-11-12): a log with a constant factor inside the
     * radical (the presentation once dropped Log[2] there: a wrong answer); the
     * exponential special places at 0 AND infinity; an ArcTan in a special value */
    "(-1 - 4*x - 4*x^2 + 8*(-3 - 2*x + Log[2 + 2*x]) + 8*x*(-3 - 2*x + Log[2 + 2*x]))/"
    "(4*(-3 - 2*x + Log[2 + 2*x])^(3/4) + 4*x*(-3 - 2*x + Log[2 + 2*x])^(3/4))",
    "(2 + Log[3 + 3*x])/(Sqrt[2 + Log[3 + 3*x]^2] + x*Sqrt[2 + Log[3 + 3*x]^2])",
    "(2*E^(4*x))/Sqrt[2 + E^(4*x)] - (2*E^(2*x)*(E^(2*x) + Sqrt[2 + E^(4*x)]))/(2 + E^(4*x) + E^(2*x)*Sqrt[2 + E^(4*x)])",
    "(-2*E^x - 2*Sqrt[-3 + 2*E^x] + 2*(-3 + 2*E^x) + (-3 + 2*E^x)^(3/2) + 2*E^x*x + E^x*Sqrt[-3 + 2*E^x]*x + "
    "(-3 + 2*E^x)*x + E^x*x^2)/(-3 + 2*E^x + 2*Sqrt[-3 + 2*E^x] + Sqrt[-3 + 2*E^x]*x)",
    NULL
};

static const char* const certificates[] = {
    "Log[1 + E^x]^(1/3)/(1 + Log[1 + E^x])",                                        /* [B] p.139 */
    "Sqrt[Log[x]]",
    "1/Sqrt[Log[x]]",
    "5/(x (Log[x] - 3) Sqrt[Log[x]^3 - 2])",                                         /* lifted Nagell-Lutz */
    "Sqrt[x + Log[x]]/x",
    NULL
};

static const char* const declines[] = {
    "Log[x] + Sqrt[Log[x] + Sqrt[Log[x]]]",                                          /* nested */
    "Sqrt[Log[x]] Sqrt[Log[x] + 1]",                                                 /* compositum */
    "Sqrt[Sqrt[2] + Log[x]]",                                                         /* a number-field certificate (Phase 4.3b) */
    "Log[x] E^(x Sqrt[x^3 + 1])",                                                    /* Shape B, curve monomial over a tower */
    NULL
};

static void test_closes(void) {
    /* Tan[x] Sqrt[Log[x]] was a decline until the trig front end (v0.6554) */
    BR_EQ("First[BronsteinIntegrate[Tan[x] Sqrt[Log[x]], x]]", "\"not elementary\"");
    for (int i = 0; closes[i]; i++) BR_ANTIDERIV(closes[i], "Bronstein");
    BR_EQ("Integrate[(x + 1)/((x Log[x] + x^2) Sqrt[x + Log[x]]), x, Method -> \"Bronstein\"]", "-2/Sqrt[x + Log[x]]");
    BR_TRUE("PossibleZeroQ[Integrate[(3 (x + E^x)^(1/3) + (2 x^2 + 3 x) E^x + 5 x^2)/(x (x + E^x)^(1/3)), x,"
            " Method -> \"Bronstein\"] - (3 x (x + E^x)^(2/3) + 3 Log[x])]");
}

static void test_certificates(void) {
    char buf[1024];
    for (int i = 0; certificates[i]; i++) {
        snprintf(buf, sizeof buf, "First[BronsteinIntegrate[%s, x]]", certificates[i]);
        BR_EQ(buf, "\"not elementary\"");
        snprintf(buf, sizeof buf, "Quiet[Head[Integrate[%s, x, Method -> \"Bronstein\"]]]", certificates[i]);
        BR_EQ(buf, "Integrate");
    }
    BR_EQ("BronsteinIntegrate[Log[1 + E^x]^(1/3)/(1 + Log[1 + E^x]), x]",
          "{\"not elementary\", \"resultant criterion ([B] Prop 2.9): R(z) has non-constant coefficients\"}");
    BR_EQ("BronsteinIntegrate[Sqrt[Log[x]], x]",
          "{\"not elementary\", \"no algebraic part: the Risch DE D b + s b = c has no solution ([B] Prop 2.15/2.18)\"}");
}

static void test_declines(void) {
    char buf[1024];
    for (int i = 0; declines[i]; i++) {
        snprintf(buf, sizeof buf, "First[BronsteinIntegrate[%s, x]]", declines[i]);
        BR_EQ(buf, "\"implementation incomplete\"");
        snprintf(buf, sizeof buf, "Quiet[Head[Integrate[%s, x, Method -> \"Bronstein\"]]]", declines[i]);
        BR_EQ(buf, "Integrate");
    }
    BR_EQ("BronsteinIntegrate[Log[x] + Sqrt[Log[x] + Sqrt[Log[x]]], x]", "{\"implementation incomplete\", \"nested radicals\"}");
    BR_EQ("BronsteinIntegrate[Log[x] E^(x Sqrt[x^3 + 1]), x]",
          "{\"implementation incomplete\", \"a further monomial over a kernel with a radical argument (transcendental Risch over an algebraic base)\"}");
    /* a non-elementary input is a silent decline, not "incomplete" */
    BR_EQ("BronsteinIntegrate[Erf[x] Sqrt[Log[x]], x]", "{\"failed\", \"non-elementary function in the integrand\"}");
}

static void test_messages_and_glue(void) {
    /* a certificate is announced once as Integrate::nonelem: Check sees it, Quiet hides it */
    BR_EQ("Check[Integrate[Sqrt[Log[x]], x, Method -> \"Bronstein\"], \"msg\"]", "\"msg\"");
    BR_EQ("Check[Quiet[Integrate[Sqrt[Log[x]], x, Method -> \"Bronstein\"]]; \"ok\", \"msg\"]", "\"ok\"");
    /* "implementation incomplete" is announced as Integrate::incomplete:
       Check sees it, Quiet hides it; a non-elementary input stays silent */
    BR_EQ("Check[Integrate[Tan[x] Sqrt[Log[x]], x, Method -> \"Bronstein\"]; \"ok\", \"msg\"]", "\"msg\"");
    BR_EQ("Check[Quiet[Integrate[Tan[x] Sqrt[Log[x]], x, Method -> \"Bronstein\"]]; \"ok\", \"msg\"]", "\"ok\"");
    BR_EQ("Check[Integrate[Erf[x] Sqrt[Log[x]], x, Method -> \"Bronstein\"]; \"ok\", \"msg\"]", "\"ok\"");
    BR_EQ("Integrate`Bronstein[1/(x Sqrt[Log[x]]), x]", "2 Sqrt[Log[x]]");
    BR_EQ("Quiet[Integrate`Bronstein[Sin[x] Sqrt[Log[x]], x]][[1]]", "\"not elementary\"");   /* trig front end (v0.6554) */
    BR_EQ("Head[Integrate`Bronstein[x]]", "Integrate`Bronstein");
    BR_TRUE("MemberQ[Attributes[Integrate`Bronstein], Protected]");
    br_assertions++;
    if (!symtab_get_docstring("Integrate`Bronstein")) {
        fprintf(stderr, "FAIL: Integrate`Bronstein has no docstring\n");
        exit(1);
    }
    /* the Method name is accepted (no Integrate::method) */
    BR_EQ("Check[Integrate[x, x, Method -> \"Bronstein\"], \"msg\"]", "1/2 x^2");
}

static void test_pipeline_routines(void) {
    /* brPresentable re-verifies the presented form, falling back to F */
    BR_EQ("brPresentable[1/x, Log[x], x]", "Log[x]");
    BR_EQ("brPresentable[Log[2 + 2 x]/(x + 1), Log[2 + 2 x]^2/2, x] // Simplify // FreeQ[#, Log[1 + x]] &", "True");
    /* the presentation is optional: bounded by the option, and by a third of an
     * enclosing budget, so it can never abort a caller holding a verified answer */
    BR_EQ("brPresentationTime[]", "5");
    BR_TRUE("TimeConstrained[0.9 < brPresentationTime[] <= 1., 3]");
    BR_TRUE("TimeConstrained[brPresentationTime[] == 5, 60]");
    BR_TRUE("Internal`TimeConstrainedRemaining[] === Infinity");
    BR_TRUE("TimeConstrained[TimeConstrained[Internal`TimeConstrainedRemaining[], 50], 2] <= 2.");
    BR_EQ("TimeConstrained[brPresentable[1/x, Log[x], x], 0.5]", "Log[x]");
    BR_EQ("brIntegrateTranscendental[x E^x, x]", "E^x (-1 + x)");
    BR_EQ("brIntegrateTranscendental[0, x]", "0");
    BR_THROWS("brIntegrateTranscendental[E^x/x, x]", "NonElementary", "transcendental part has no elementary integral (Risch)");
    BR_EQ("brTermExpression[None, x^2]", "x^2");
    BR_EQ("brTermExpression[fieldFn[{x}], x]", "x");
    BR_EQ("brTermExpression[fieldFn[{x}, {t, t^2 - 2, Sqrt[2]}], t x + t^3]", "2 Sqrt[2] + Sqrt[2] x");
    br_run("T134 = brTower[(x + 1)/((x Log[x] + x^2) Sqrt[x + Log[x]]), x]; True");
    br_run("cv134 = brCurve[T134, 2, brT1 + x]; res134 = brIntegrateComponent[T134, cv134, (1 + x)/(x (brT1 + x)^2)]; True");
    BR_EQ("{res134[\"Algebraic\"], res134[\"Logs\"]}", "{-2/(brT1 + x), {}}");
    BR_EQ("brAssembleComponent[T134, res134, Sqrt[brT1 + x], x]", "-2/Sqrt[x + Log[x]]");
    BR_EQ("brShapeA[T134, (x + 1)/((x Log[x] + x^2) Sqrt[x + Log[x]]), x][\"F\"]", "-2/Sqrt[x + Log[x]]");
    BR_TRUE("brShapeA[T134, (x + 1)/((x Log[x] + x^2) Sqrt[x + Log[x]]), x][\"verified\"]");
    /* the exact per-component identity: accepts the true integral, rejects a wrong one */
    BR_TRUE("brComponentIdentityQ[T134, cv134, (1 + x)/(x (brT1 + x)^2), res134]");
    BR_FALSE("brComponentIdentityQ[T134, cv134, (1 + x)/(x (brT1 + x)^2), Join[res134, <|\"Algebraic\" -> -3/(brT1 + x)|>]]");
    BR_FALSE("brComponentIdentityQ[T134, cv134, 1/x, res134]");
    /* with logarithms over a residue field: the arcsinh component */
    br_run("TL1b = brTower[1/(x Sqrt[Log[x]^2 + 1]), x]; True");
    br_run("cvL1b = brCurve[TL1b, 2, brT1^2 + 1]; resL1b = brIntegrateComponent[TL1b, cvL1b, 1/((1 + brT1^2) x)]; True");
    BR_TRUE("brComponentIdentityQ[TL1b, cvL1b, 1/((1 + brT1^2) x), resL1b]");
    BR_FALSE("brComponentIdentityQ[TL1b, cvL1b, 2/((1 + brT1^2) x), resL1b]");
    BR_THROWS("brIntegrateComponent[T134, cv134, 1/x]", "NonElementary",
              "no algebraic part: the Risch DE D b + s b = c has no solution ([B] Prop 2.15/2.18)");
    BR_THROWS("brIntegrateComponent[T134, <|\"m\" -> 2, \"G\" -> (brT1 + x)^2, \"s\" -> 0|>, 1]", "Internal", "brCurveNotReduced");
    BR_EQ("brIntegrateBody[1/(x Sqrt[Log[x]]), x]", "2 Sqrt[Log[x]]");
    BR_THROWS("brIntegrateBody[Sqrt[Log[x]], x]", "NonElementary",
              "no algebraic part: the Risch DE D b + s b = c has no solution ([B] Prop 2.15/2.18)");
    BR_EQ("BronsteinIntegrate[x^2, x]", "1/3 x^3");
    /* brSplitLogs: a Y-free argument is split into factors (constants dropped),
     * equal arguments merged, zero sums dropped; a Y-argument is kept */
    BR_EQ("brSplitLogs[None, {{1, 6 x^2 (1 + x)}}, brYq, {x}]", "{{2, x}, {1, 1 + x}}");
    BR_EQ("brSplitLogs[None, {{1, x^2}, {-2, 3 x}}, brYq, {x}]", "{}");
    BR_EQ("brSplitLogs[None, {{1/2, (x + 1)/x}, {1/2, x}}, brYq, {x}]", "{{1/2, 1 + x}}");
    BR_EQ("brSplitLogs[None, {{1, x + brYq}, {1, x}}, brYq, {x}]", "{{1, brYq + x}, {1, x}}");
    BR_EQ("brSplitLogs[None, {}, brYq, {x}]", "{}");
    BR_EQ("brSplitLogs[fieldFn[{x}, {t, t^2 + 3, I Sqrt[3]}], {{t, 2 t x}, {-t, x}}, brYq, {x}]", "{}");
    /* property: the derivative of Sum c Log[u] is unchanged */
    BR_TRUE("SeedRandom[33]; And @@ Table[Module[{L = Table[{RandomInteger[{-2, 2}], RandomInteger[{1, 5}] x^RandomInteger[{0, 2}] (x + RandomInteger[{1, 3}])^RandomInteger[{-1, 2}]}, {3}], S},"
            "   S = brSplitLogs[None, L, brYq, {x}];"
            "   Together[D[Total[#[[1]] Log[#[[2]]] & /@ L] - Total[#[[1]] Log[#[[2]]] & /@ S], x]] === 0], {20}]");
    /* brCollectLogTerms: equal logarithms merged at the surface */
    BR_EQ("brCollectLogTerms[2/3 Log[1 + x] + x - 2/3 Log[1 + x], x]", "x");
    BR_EQ("brCollectLogTerms[Log[1 + x] + Log[1 + x]/2 + Sqrt[x], x]", "Sqrt[x] + 3/2 Log[1 + x]");
    BR_EQ("brCollectLogTerms[Log[2] + Log[x], x]", "Log[2] + Log[x]");
    BR_EQ("brCollectLogTerms[Log[x], x]", "Log[x]");
    BR_EQ("brCollectLogTerms[0, x]", "0");
    /* every public name of another package is resolved (an unresolved one becomes an
     * undefined capitalised Bronstein`Private` symbol: an unevaluated call that once
     * produced a false certificate) */
    BR_EQ("Select[Names[\"Bronstein`Private`*\"], UpperCaseQ[StringTake[Last[StringSplit[#, \"`\"]], 1]] &]", "{}");
}

/* Property: never a wrong answer on derivatives of random elementary
 * functions b(theta) Y + c Log[p + q Y] (theta Log or Exp, Y a square or cube
 * root of a radicand linear in theta); most of them close. */
static void test_random_derivatives(void) {
    BR_TRUE("SeedRandom[91]; brOK = 0; brAll = 0;"
            " And @@ Table[Module[{th = RandomChoice[{Log[x], E^x, Log[x + 1]}], n = RandomChoice[{2, 2, 3}], Y, P, Fx, f, r},"
            "   P = th + RandomInteger[{1, 3}] x + RandomInteger[{0, 2}];"
            "   Y = P^(1/n);"
            "   Fx = (RandomInteger[{1, 3}] x + RandomInteger[{-2, 2}] th) Y + RandomChoice[{0, 1}] Log[x + Y];"
            "   f = Together[D[Fx, x]];"
            "   r = BronsteinIntegrate[f, x]; brAll++;"
            "   If[ListQ[r], First[r] =!= \"not elementary\", brOK++; TrueQ[PossibleZeroQ[D[r, x] - f]]]], {16}]");
    BR_TRUE("brOK >= 12");
}

int main(void) {
    br_init();
    test_closes();
    test_certificates();
    test_declines();
    test_messages_and_glue();
    test_pipeline_routines();
    test_random_derivatives();
    BR_DONE("bronstein_integrate_tests");
    return 0;
}
