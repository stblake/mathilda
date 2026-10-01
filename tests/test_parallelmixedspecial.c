/* Tests for Integrate`ParallelMixedSpecial -- the special-function stage of the
 * parallel (Risch-Norman) integrator over a mixed tower.
 *
 * Correctness of an antiderivative is asserted by a NUMERIC RESIDUAL at a
 * rational point rather than by a fixed output string, so the tests survive a
 * change of surface form -- this stage legitimately returns a different (and
 * here simpler) Legendre reduction of the same elliptic pencil than the
 * Mathematica reference does, and both are right.
 *
 * ALWAYS with a head guard alongside the residual. `D[unevaluated Integrate]`
 * returns the integrand, so a residual test on its own lets a DECLINE pass
 * vacuously -- the single most important idiom in this file, inherited from
 * test_parallelmixedtower.c.
 *
 * Assertions use ASSERT_MSG / ASSERT_STR_EQ (hard exit(1)), NOT assert_eval_eq,
 * whose libc assert() is a no-op under -DNDEBUG.
 *
 * The erf / incomplete-Gamma family ANSWERS, as of v0.241. Until then it did not:
 * it computed the right answer and strict mode (Theorem 8.6) withheld it, because
 * Part II's certificate was hard-gated to a tower carrying a curve -- the
 * `q =!= None` T2/T10 false-certificate protection. That gate was belt-and-braces
 * over a residue-realisation bug fixed at its root, so it retired; the boundary
 * this file pins therefore MOVED rather than softened. What guards the soundness
 * now is test_no_false_certificate_curve_free in test_parallelmixedtower.c, which
 * asserts that sixteen ELEMENTARY curve-free integrands still SOLVE -- stronger
 * than "do not certify", since a solution proves the certificate branch (entered
 * only on an inconsistent system) was never reached. A false certificate remains
 * the failure mode that matters most.
 */
#define _POSIX_C_SOURCE 200809L

#include "attr.h"
#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char* eval_str(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "parse failed: %s", input);
    Expr* r = evaluate(e);
    expr_free(e);
    char* s = expr_to_string(r);
    expr_free(r);
    return s;
}

static void assert_eval(const char* input, const char* expected) {
    char* s = eval_str(input);
    ASSERT_STR_EQ(s, expected);
    free(s);
}

static void assert_true_msg(const char* input, const char* what) {
    char* s = eval_str(input);
    ASSERT_MSG(strcmp(s, "True") == 0, "%s: expected True, got %s", what, s);
    free(s);
}

/* The workhorse: the method answers, the head really changed, and the answer
 * differentiates back to the integrand at `pt`. */
static void assert_solves(const char* f, const char* pt) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Module[{r = Integrate`ParallelMixedSpecial[%s, x]},"
        " Head[r] =!= Integrate`ParallelMixedSpecial && Head[r] =!= List &&"
        " Abs[N[(D[r, x] - (%s)) /. x -> %s, 30]] < 10^-20]", f, f, pt);
    char* s = eval_str(buf);
    ASSERT_MSG(strcmp(s, "True") == 0,
               "ParallelMixedSpecial(%s): expected an answer differentiating "
               "back at x = %s, got %s", f, pt, s);
    free(s);
}

/* The method declines, and its status is NOT a certificate. A {"failed", ...}
 * proves nothing and is honest; a {"not elementary", ...} here would be a false
 * certificate. */
static void assert_honest_decline(const char* f) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Module[{r = Integrate`ParallelMixedSpecial[%s, x]},"
        " ListQ[r] && r[[1]] =!= \"not elementary\" && r[[1]] =!= \"not in class\"]", f);
    char* s = eval_str(buf);
    ASSERT_MSG(strcmp(s, "True") == 0,
               "ParallelMixedSpecial(%s): expected an honest {\"failed\", ...} "
               "decline, got %s", f, s);
    free(s);
}

/* ---- the kernels, one case per kind ---------------------------------- */

static void test_kernel_ei(void) {
    assert_solves("Exp[x]/x", "7/5");                 /* ExpIntegralEi        */
    assert_solves("Exp[2 x]/x", "7/5");
}

static void test_kernel_li(void) {
    /* Through Present[]: Ei[Log[a]] is rewritten to LogIntegral[a], which is
     * what the reference emits and a better form than the Cherry engine's. */
    assert_solves("1/Log[x]", "9/4");
    assert_true_msg("FreeQ[Integrate`ParallelMixedSpecial[1/Log[x], x], ExpIntegralEi]",
                    "1/Log[x] is presented as LogIntegral, not Ei[Log[..]]");
}

static void test_kernel_sici(void) {
    assert_solves("Sin[x]/x", "7/5");                 /* SinIntegral          */
    assert_solves("Cos[x]/x", "7/5");                 /* CosIntegral          */
}

static void test_kernel_elliptic(void) {
    /* An elliptic pencil: y^2 = x^3 - x, genus 1.  The reference answers
     * -Sqrt[2] EllipticF[ArcSin[Sqrt[2]/Sqrt[1+x]], 1/2]; this stage finds
     * (2 I) EllipticF[ArcSin[Sqrt[x]], -1], a different Legendre reduction of
     * the same pencil differing by a constant.  Hence the residual, not a
     * string -- and hence the head guard, since a decline would pass a residual
     * test vacuously. */
    assert_solves("1/Sqrt[x^3 - x]", "7/5");
    assert_solves("x/Sqrt[x^3 - x]", "9/4");
    assert_true_msg("! FreeQ[Integrate`ParallelMixedSpecial[1/Sqrt[x^3 - x], x], "
                    "EllipticF | EllipticE | EllipticPi]",
                    "the elliptic pencil answers with an elliptic integral");
}

/* An elementary integrand must come back elementary -- a special function here
 * would mean the stage reached for a kernel it did not need (Section 8.1, E1). */
static void test_elementary_stays_elementary(void) {
    assert_solves("Log[x]", "7/5");
    assert_true_msg("FreeQ[Integrate`ParallelMixedSpecial[Log[x], x], "
                    "ExpIntegralEi | LogIntegral | Erf | Erfc | Erfi | SinIntegral | "
                    "CosIntegral | EllipticF | EllipticE | EllipticPi]",
                    "Log[x] integrates without a special-function kernel");
    assert_eval("Integrate`ParallelMixedSpecial[Log[x], x]", "-x + x Log[x]");
}

/* ---- the certificate boundary ---------------------------------------- */

/* The erf / incomplete-Gamma family: the answer was always computed and always
 * the reference's; what was missing was the certificate strict mode (Theorem
 * 8.6) demands before releasing it. Part II now certifies on a curve-free tower
 * -- proved bounds, a complete logand set, a verified residue-free residual,
 * which is Proposition 9.2(b) itself rather than its second-kind specialisation
 * -- so these ANSWER. Before v0.241 every one of them was an HONEST decline. */
static void test_certificate_now_released(void) {
    assert_solves("Exp[-x^2]", "7/5");
    assert_solves("Exp[x^2]", "7/5");
    assert_solves("Exp[-x^3]", "7/5");
    assert_solves("Sqrt[x] Exp[x]", "7/5");
    /* and the released answer is the special-function one, not an elementary
     * near-miss the stage talked itself into */
    assert_true_msg(
        "! FreeQ[Integrate`ParallelMixedSpecial[Exp[-x^2], x], Erf | Erfc | Erfi]",
        "Exp[-x^2] answers with an error function");
    assert_true_msg(
        "! FreeQ[Integrate`ParallelMixedSpecial[Exp[-x^3], x], Gamma]",
        "Exp[-x^3] answers with an incomplete Gamma");
}

/* The tangent tower, both ways it used to fail (corpus #87/#172, #175, #89, #90).
 *
 * A deep (order-2) pole: the Ei columns were always found -- that is the Ci --
 * but the ELEMENTARY remainder they leave, -Sin[x]/x, was rejected because Can's
 * number-field detour silently failed to canonicalise it (see Can in
 * ParallelMixed.m).
 *
 * A non-linear tangent argument: Sin[x^2] needs the conjugate Erf PAIR, which
 * needs W = I xi^2/x^2 recognised as an exact square -- and its numerator
 * -I (t - I)^2 is irreducible over Q, so ExactRoot must factor over Q(i). The
 * freed constant is carried as the kernel's cw, so Sqrt[I] never enters the
 * constant field. */
static void test_tangent_tower_kernels(void) {
    assert_solves("Sin[x]/x^2", "7/5");
    assert_solves("Sin[3 x]/x^2", "7/5");
    assert_solves("Sin[x^2]", "7/5");
    assert_solves("Cos[x^2]", "7/5");
    assert_true_msg(
        "! FreeQ[Integrate`ParallelMixedSpecial[Sin[x]/x^2, x], CosIntegral]",
        "Sin[x]/x^2 answers with CosIntegral");
    assert_true_msg(
        "! FreeQ[Integrate`ParallelMixedSpecial[Sin[x^2], x], Erf | Erfc | Erfi]",
        "Sin[x^2] answers with an error function");
}

/* Branch fidelity: the second embedding of the constant field. The ansatz is
 * solved over Q(i) and nothing in the linear algebra chooses between i -> +-i;
 * both solve, and only one differentiates back once the surface radicals are
 * read on their principal branches. Corpus #34, the last FAIL of the 312: the
 * surface is 2 Sqrt[x] ArcSin[Sqrt[1+x]] - 2 I Sqrt[1+x] and the antiderivative
 * is the same with +2 I Sqrt[1+x], residual -2 I/Sqrt[1+x] otherwise.
 *
 * Note this integrand is COMPLEX-VALUED on the whole real line (ArcSin of
 * something bigger than 1), which is the other half of the fix: the verify gate
 * used to fall back to genuinely complex sample points, where a branch cut puts
 * surface and integrand on opposite sides and neither embedding verifies. */
static void test_second_embedding_branch(void) {
    assert_solves("ArcSin[Sqrt[x + 1]]/Sqrt[x]", "7/5");
    /* and Part II alone, which is where the gate lives */
    assert_true_msg(
        "Module[{r = Integrate`ParallelMixedTower[ArcSin[Sqrt[x + 1]]/Sqrt[x], x]},"
        " Head[r] =!= List && "
        " Abs[N[(D[r, x] - ArcSin[Sqrt[x + 1]]/Sqrt[x]) /. x -> 9/4, 30]] < 10^-20]",
        "Part II verifies the second-embedding surface");
}

/* The floor: an integrand the stage still cannot close must decline HONESTLY --
 * a {"failed", ...} status, never a certificate. The Cherry family below is
 * declined by all four reference ports too (corpus #274-#277), so it is the
 * honest boundary rather than a Mathilda gap. */
static void test_honest_decline_still_honest(void) {
    assert_honest_decline("Exp[-Log[x]^2]");
    assert_honest_decline("Exp[-Log[x]^2]/x");
    assert_honest_decline("x^2 Exp[-Log[x]^2]");
}

/* An integrand genuinely outside the class must decline, never answer. Here
 * {"not in class", ...} is the RIGHT status -- it is what the corpus's `nic`
 * group scores a PASS on -- so unlike assert_honest_decline this only forbids a
 * claimed ANSWER and a claimed non-elementarity certificate. */
static void test_not_in_class(void) {
    assert_true_msg(
        "Module[{r = Integrate`ParallelMixedSpecial[Exp[Exp[Exp[x]]]/x, x]},"
        " ListQ[r] && r[[1]] =!= \"not elementary\"]",
        "an out-of-class integrand declines without claiming non-elementarity");
}

/* ---- the three surfaces ---------------------------------------------- */

static void test_surfaces(void) {
    /* (a) the qualified symbol */
    assert_eval("Integrate`ParallelMixedSpecial[Exp[x]/x, x]", "ExpIntegralEi[x]");
    /* (b) Method -> */
    assert_eval("Integrate[Exp[x]/x, x, Method -> \"ParallelMixedSpecial\"]",
                "ExpIntegralEi[x]");
    /* (c) the Automatic cascade reaches it for an integrand nothing cheaper
     *     closes -- and the answer differentiates back. */
    assert_true_msg(
        "Module[{r = Integrate[1/Sqrt[x^3 - x], x]},"
        " Head[r] =!= Integrate &&"
        " Abs[N[(D[r, x] - 1/Sqrt[x^3 - x]) /. x -> 7/5, 30]] < 10^-20]",
        "plain Integrate closes the elliptic pencil through the cascade");
    /* The method name is listed in the Integrate::method diagnostic. */
    assert_true_msg("StringContainsQ[Information[\"Integrate\"], "
                    "\"ParallelMixedSpecial\"]",
                    "the method is documented in Integrate's docstring");
}

/* The cascade must NEVER return a partial answer from plain Integrate: the
 * partial mode is reachable only through the package's own entry point. */
static void test_cascade_never_partial(void) {
    assert_true_msg("FreeQ[Integrate[Log[x]/(x + 1) + 1/x, x], Inactive]",
                    "plain Integrate returns no Inactive[Integrate] remainder");
}

/* ---- the partial mode, through the package ------------------------- */

static void test_partial_mode(void) {
    /* IntegrateSurfacePartial answers `elementary part + Inactive[Integrate][r, x]`.
     * Reached by its qualified private name: this module is loaded INTO
     * ParallelMixed`Private`. */
    assert_true_msg(
        "Module[{r}, Quiet[Integrate`ParallelMixedSpecial[x, x]];"      /* lazy-load */
        " r = ParallelMixed`Private`IntegrateSurfacePartial["
        "       Log[x]/(x + 1) + 1/x, x];"
        " ListQ[r] && Length[r] == 2 && r[[2]] === True &&"
        " ! FreeQ[r[[1]], Inactive[Integrate]]]",
        "the partial mode returns answer + Inactive[Integrate][remainder, x]");
}

/* ---- the additive decomposition --------------------------------------- */

/* Cost in this stage is MULTIPLICATIVE in the number of independent generator
 * families appearing additively, because one ansatz spans the product of the
 * families. Splitting the integrand into blocks whose families are pairwise
 * disjoint and integrating each on its own tower turns 20.1 s into 0.29 s on
 * the three-family case below.
 *
 * What these tests pin is the SOUNDNESS GUARD, not the speed. Disjoint
 * families do NOT by themselves make the split complete: for two primitives
 * whose derivatives lie in the shared base,
 *     D(c t_A t_B) = c (D t_A) t_B + c t_A (D t_B)
 * is a two-term sum with DISJOINT support, so a sum can be elementary while
 * each half is not even in the class. The guard declines whenever two or more
 * COUPLING-CAPABLE blocks (carrying a Log or an ArcTan-family head at a
 * rational argument) fail to close, and the joint path then runs unchanged. */
static void test_additive_split_soundness(void) {
    /* The counterexample. Each half alone is "not in class"; the sum is
     * Log[x] Log[1+x]. The families {Log[x]} and {Log[1+x]} are disjoint, so a
     * naive split would turn this PASS into a decline. */
    assert_solves("Log[1 + x]/x + Log[x]/(1 + x)", "7/5");
    assert_true_msg("Module[{r = Integrate`ParallelMixedSpecial["
                    "  Log[1 + x]/x + Log[x]/(1 + x), x]},"
                    " FreeQ[r, Inactive[Integrate]] && FreeQ[r, List]]",
                    "the two-coupling-capable-block sum stays COMPLETE");
    /* Same shape through the inverse-tangent family. Both halves are corpus
     * `nic` cases; the sum is ArcTan[x] Log[x]. */
    assert_solves("ArcTan[x]/x + Log[x]/(1 + x^2)", "7/5");

    /* A block's STATUS must never become the sum's. Each of these is out of
     * class on its own and must still report so, not be composed away. */
    assert_true_msg("Module[{r = Integrate`ParallelMixedSpecial[ArcTan[x]/x, x]},"
                    " ListQ[r] && r[[1]] === \"not in class\"]",
                    "a single out-of-class integrand still says so");

    /* Cross-term cancellation inside ONE family: each half is non-elementary
     * (Ei), the sum is E^x/x. These share the generator Exp[x], so they land
     * in the same block and the split never separates them. */
    assert_solves("Exp[x]/x - Exp[x]/x^2", "7/5");
    assert_true_msg("FreeQ[Integrate`ParallelMixedSpecial[Exp[x]/x - Exp[x]/x^2, x],"
                    " ExpIntegralEi]",
                    "the Ei pair cancels to an elementary answer, as before");
}

/* The split must agree with the joint path. "Additive" -> False is a real
 * kill switch, so the old path stays live and directly comparable.
 *
 * Compared through the PARTIAL entry, which is uniform over both outcomes:
 * #129 and #134 below have no complete answer at all (their Exp[-x^2] Log[x]
 * term is out of class), so the complete entry answers them with a status,
 * not a pair. Two things are asserted, and the second is the one that
 * matters: the answers agree up to a constant, AND the inert remainders are
 * equal -- the remainder is what the corpus judge scores, so a split that
 * handed back a differently-apportioned remainder would be a behaviour
 * change even with a correct total. */
static void test_additive_split_agrees(void) {
    const char* cases[] = {
        "2 Log[x]/x + Sin[x]/x",                      /* complete, 24x        */
        "Log[x]/x + Sin[x]/x + Exp[-x^2]",            /* complete, 70x        */
        "Sin[x]/x + Exp[-x^2]",                       /* complete             */
        "x Exp[x]/(x + 1)^2 + Exp[-x^2] Log[x]",      /* partial, 8.9x        */
        "Exp[-x^2] Log[x] - 3 Log[x]/x^2 + 3/x^2",    /* partial, one block   */
        "1/(x + Exp[x]) + Cos[x]/x",                  /* partial, 3.9x        */
        "(Log[x] + 1/x) Exp[x] + 1/x",                /* elementary, one block*/
        NULL
    };
    char buf[2560];
    for (int i = 0; cases[i]; i++) {
        snprintf(buf, sizeof(buf),
            "Module[{a, b, ra, rb, pa, pb},"
            " Quiet[Integrate`ParallelMixedSpecial[x, x]];"
            " a = ParallelMixed`Private`IntegrateSurfacePartial[%s, x, \"Additive\" -> True];"
            " b = ParallelMixed`Private`IntegrateSurfacePartial[%s, x, \"Additive\" -> False];"
            " If[! (ListQ[a] && ListQ[b] && Length[a] == 2 && Length[b] == 2), False,"
            "  ra = Total[Cases[{a[[1]]}, Inactive[Integrate][q_, x] :> q, Infinity]];"
            "  rb = Total[Cases[{b[[1]]}, Inactive[Integrate][q_, x] :> q, Infinity]];"
            "  pa = a[[1]] /. Inactive[Integrate][_, x] -> 0;"
            "  pb = b[[1]] /. Inactive[Integrate][_, x] -> 0;"
            "  a[[2]] === b[[2]] && Simplify[ra - rb] === 0 &&"
            "  Simplify[D[pa - pb, x]] === 0]]", cases[i], cases[i]);
        char* s = eval_str(buf);
        ASSERT_MSG(strcmp(s, "True") == 0,
                   "additive split disagrees with the joint path on %s: got %s",
                   cases[i], s);
        free(s);
    }
}

/* A complete split answer must carry NO Inactive[Integrate]: the joint path
 * reaches its Inactive[Integrate] only for a genuine partial result, so an
 * Inactive[Integrate][0, x] would make a complete answer read as partial --
 * which the corpus judge scores as a different outcome. */
static void test_additive_split_complete_is_complete(void) {
    assert_true_msg(
        "Module[{r}, Quiet[Integrate`ParallelMixedSpecial[x, x]];"
        " r = ParallelMixed`Private`IntegrateSurfacePartial["
        "       2 Log[x]/x + Sin[x]/x, x];"
        " ListQ[r] && r[[2]] === True && FreeQ[r[[1]], Inactive[Integrate]]]",
        "a complete split answer carries no Inactive[Integrate]");
    /* and a genuinely partial one still does */
    assert_true_msg(
        "Module[{r}, Quiet[Integrate`ParallelMixedSpecial[x, x]];"
        " r = ParallelMixed`Private`IntegrateSurfacePartial["
        "       1/(x + Exp[x]) + Cos[x]/x, x];"
        " ListQ[r] && r[[2]] === True && ! FreeQ[r[[1]], Inactive[Integrate]] &&"
        " ! FreeQ[r[[1]], CosIntegral]]",
        "a partial split answer keeps its remainder inert");
}

/* Two different radicals in two different additive terms: each block builds
 * its own tower. The joint path also answers these, but through one Euler
 * substitution carrying both radicals, which returns a far heavier form --
 * so what is pinned here is that the split answers and is FAITHFUL, i.e. the
 * split answer is exactly the sum of what the blocks give on their own.
 *
 * That faithfulness property is the real invariant of this transform, and it
 * is worth asserting separately from correctness: when it holds, a wrong
 * split answer can only come from a wrong block answer, never from the
 * splitting itself. (It is how the Sqrt[1+x^3] case was cleared during
 * development -- the split reproduced the block sum exactly, and the residual
 * came entirely from that block's own elliptic branch, which is equally wrong
 * through the ordinary single-integrand path.) */
static void test_additive_split_two_radicals(void) {
    assert_solves("Sqrt[1 - x^2] + Sqrt[1 + x^2]", "1/3");
    assert_solves("Sqrt[1 - x^2] + 1/Sqrt[4 + x^2]", "1/3");
    assert_true_msg(
        "Module[{s, a, b}, Quiet[Integrate`ParallelMixedSpecial[x, x]];"
        " s = ParallelMixed`Private`IntegrateSurfaceSpecial[Sqrt[1 - x^2] + Sqrt[1 + x^2], x];"
        " a = ParallelMixed`Private`IntegrateSurfaceSpecial[Sqrt[1 - x^2], x];"
        " b = ParallelMixed`Private`IntegrateSurfaceSpecial[Sqrt[1 + x^2], x];"
        " ListQ[s] && ListQ[a] && ListQ[b] &&"
        " Simplify[s[[1]] - (a[[1]] + b[[1]])] === 0]",
        "the split answer is exactly the sum of the per-block answers");
}

/* ---- Part II is untouched by the merge -------------------------------- */

/* The stage's additive delta to Part II is off by default, so
 * ParallelMixedTower must answer exactly as before.  These duplicate a few
 * test_parallelmixedtower.c cases on purpose: they are the tripwire for the
 * hand-merge, and a failure here localises to Part II rather than to this stage. */
static void test_part2_unchanged(void) {
    assert_eval("Simplify[D[Integrate`ParallelMixedTower[Log[x], x], x] - Log[x]]", "0");
    assert_true_msg("Module[{r = Integrate`ParallelMixedTower[Sqrt[x + Sqrt[x]], x]},"
                    " Head[r] =!= List && FreeQ[r, Power[x, 1/4]] &&"
                    " Abs[N[(D[r, x] - Sqrt[x + Sqrt[x]]) /. x -> 2 + I, 25]] < 10^-15]",
                    "Part II still fuses the nested radical and solves it");
    /* The structure theorem, the delta's visible new capability: Exp[2x] is
     * recognised as t^2 over t = E^x. */
    assert_true_msg("Module[{r = ParallelMixed`ParallelIntegrateMixed["
                    "  Exp[2 x]/(1 + Exp[x]), x, \"StructureTheorem\" -> True]},"
                    " Head[r] =!= List &&"
                    " Abs[N[(D[r, x] - Exp[2 x]/(1 + Exp[x])) /. x -> 7/5, 25]] < 10^-15]",
                    "BuildTower's structure theorem closes Exp[2x]/(1+Exp[x])");
}

/* ExtendedBounds Blocks four of Part II's bound-decision symbols and defines
 * DownValues on them inside. It must restore them -- otherwise the FIRST
 * extended call silently repoints Part II for the rest of the session. */
static void test_extended_bounds_restores(void) {
    assert_true_msg(
        "Module[{before, inside, after},"
        " Quiet[Integrate`ParallelMixedSpecial[x, x]];"
        " before = ! FreeQ[DownValues[ParallelMixed`Private`DecideBound],"
        "                  ParallelMixed`Private`DecideBoundOrig];"
        " inside = ParallelMixed`Private`ExtendedBounds[True,"
        "            ! FreeQ[DownValues[ParallelMixed`Private`DecideBound],"
        "                    ParallelMixed`Private`DecideBoundExt]];"
        " after  = ! FreeQ[DownValues[ParallelMixed`Private`DecideBound],"
        "                  ParallelMixed`Private`DecideBoundOrig];"
        " before && inside && after]",
        "ExtendedBounds installs the extended hooks and restores the originals");
}

int main(void) {
    symtab_init();
    core_init();
    /* Production-like environment + chdir to the repo root so the lazy loader
     * resolves src/internal/ParallelMixedSpecial.m. */
    test_load_init_m();

    TEST(test_kernel_ei);
    TEST(test_kernel_li);
    TEST(test_kernel_sici);
    TEST(test_kernel_elliptic);
    TEST(test_elementary_stays_elementary);
    TEST(test_certificate_now_released);
    TEST(test_tangent_tower_kernels);
    TEST(test_second_embedding_branch);
    TEST(test_honest_decline_still_honest);
    TEST(test_not_in_class);
    TEST(test_surfaces);
    TEST(test_cascade_never_partial);
    TEST(test_partial_mode);
    TEST(test_additive_split_soundness);
    TEST(test_additive_split_agrees);
    TEST(test_additive_split_complete_is_complete);
    TEST(test_additive_split_two_radicals);
    TEST(test_part2_unchanged);
    TEST(test_extended_bounds_restores);

    printf("All ParallelMixedSpecial tests passed!\n");
    return 0;
}
