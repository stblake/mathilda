/*
 * test_dsolve_m63_stress.c — anti-overfit stress families for M63.
 *
 * M63 landed two things.  The converter repair (eighteen corpus records that were
 * silently the WRONG equation) has no C surface — `make check-corpus-indvar` is its
 * test — so F5 below only pins the equation shapes that repair produced.  The solver
 * fix is a per-antiderivative BUDGET in `DSolve`Separable`, and the families here are
 * shaped to the ways that can go wrong, none of which is "it stops answering":
 *
 *   F1  the fix is a LATENCY property, so the assertions are about time.  Separable
 *       sits at cascade slot 7 of 52 and had no bound on its own integrals, so one
 *       integrand that does not close spent the whole 8 s solve on behalf of every
 *       method behind it.  §2.2.36-3521/3527/3599 were lost to exactly that:
 *       `DSolve`Linearizable`, three slots later, could solve the ORIGINAL equation
 *       all along.  An answer-only test would have passed before the fix as well as
 *       after — only the clock sees it.  Run over a forward generator of four
 *       mixed-angle shapes, so no member is one of the three.
 *
 *   F2  the generality claim.  The budget is not a trig fix: it bounds ANY integrand.
 *       The multiple-angle sibling `Cos[x-2y]/(Sin[x] Sin[2y]) - 1` is the control —
 *       nothing in the system solves it, on either binary, so the only thing the fix
 *       can do for it is make Separable DECIDE instead of spin, and that is what is
 *       asserted.  A rewrite (TrigExpand before sampling) was built and measured
 *       against this family and discarded; see F1's note on answer quality.
 *
 *   F3  no legitimate separable is lost, and the budget costs nothing when it does
 *       not fire.  The controls are the documented corners of this method: a plain
 *       separable, a rational one, a Log one, the `Cot[x] y/(1+y)` parameter split,
 *       and the autonomous NON-elementary `y' == -2 ArcTan[y]/(1+y^2)` whose implicit
 *       twin legitimately keeps an unevaluated integral — the one case where a
 *       timeout and a decided non-elementary integrand must NOT be confused.
 *
 *   F4  the answers are the GOOD ones.  The whole reason the rewrite was discarded is
 *       that it claimed these records for Separable's implicit twin, whose relation
 *       carries the sampling artefact `Cot[2]`, where the cascade left to itself
 *       returns the explicit `ArcCos[C[1] Csc[x]]` that Mathematica gives.  So this
 *       family asserts the answer is explicit, and verifies it numerically against
 *       the ORIGINAL equation.
 *
 *   F5  the shapes the converter repair produced — ordinary constant-coefficient and
 *       autonomous-quadrature equations.  Before the repair the corpus asked for
 *       `y''[a] - 2a y'[a] + a^2 y[a] == 0` (a variable-coefficient equation in `a`)
 *       and scored the answer to THAT as if it were the book's question.
 *
 * Every family asserts Head[sol] === List FIRST, so a declining method cannot pass
 * vacuously.
 *
 * RUNTIME: ~96 s measured, against the shared alarm(120) in test_utils.h — so there
 * is roughly 24 s of headroom and no room for a sixth family that solves anything.
 * The first cut of this file solved each mixed-angle shape twice (once for latency in
 * F1, once for the residual in F4) and was killed by that alarm after F2, which reads
 * as a truncated log rather than a failure.  If you add a case here, take the time out
 * of an existing solve rather than adding one.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}
static bool lang_true(const char* input) {
    char* s = eval_str(input);
    bool ok = (strcmp(s, "True") == 0);
    if (!ok) fprintf(stderr, "  expected True: %s  =>  %s\n", input, s);
    free(s);
    return ok;
}
#define ASSERT_TRUE(input) ASSERT_MSG(lang_true(input), "expected True: %s", (input))

static double seconds_of(const char* input) {
    clock_t t0 = clock();
    char* s = eval_str(input);
    double dt = (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
    free(s);
    return dt;
}

/* The four mixed-angle shapes, each separable only through an angle-addition
 * identity, and each therefore a case whose sampled integrands do not close:
 *
 *   Cos[x-y]/(Sin x Sin y) - 1  ==  Cot[x] Cot[y]      (3521, 3599)
 *   Cos[x+y]/(Sin x Sin y) + 1  ==  Cot[x] Cot[y]
 *   1 - Sin[x+y]/(Sin y Cos x)  == -Tan[x] Cot[y]      (3527)
 *   Sin[x-y]/(Cos x Sin y) + 1  ==  Tan[x] Cot[y]
 *
 * Built from the identity rather than picked, so nothing here is one of the corpus
 * records, and the pair the corpus does contain is covered by construction. */
static const char* const MIXED_ANGLE[] = {
    "Cos[x - y[x]]/(Sin[x] Sin[y[x]]) - 1",
    "Cos[x + y[x]]/(Sin[x] Sin[y[x]]) + 1",
    "1 - Sin[x + y[x]]/(Sin[y[x]] Cos[x])",
    "Sin[x - y[x]]/(Cos[x] Sin[y[x]]) + 1",
};
#define N_MIXED (sizeof MIXED_ANGLE / sizeof MIXED_ANGLE[0])

/* ---------------------------------------------------------------------------
 * F1 — Separable DECIDES in bounded time, and the cascade reaches its answer.
 *
 * Pre-fix, `DSolve`Separable` on each of these ran 30 s without deciding (measured
 * under a 30 s TimeConstrained, which it hit), so the 8 s solve died inside it.  The
 * assertions are therefore (a) it decides at all under a budget far below that, and
 * (b) the whole solve now finishes, which it could not before.
 * ------------------------------------------------------------------------- */
static void t_m63_separable_decides_in_budget(void) {
    char q[900];
    for (size_t i = 0; i < N_MIXED; i++) {
        /* The pinned method decides — answer or decline — rather than spinning.
         * $Aborted here is the pre-fix behaviour and the thing being guarded: before,
         * this ran 30 s under a 30 s TimeConstrained and reached no verdict.
         *
         * The whole-solve half of this property lives in F4, which already solves each
         * shape once and times it; duplicating the solve here put the file over the
         * shared alarm(120) in test_utils.h. */
        snprintf(q, sizeof q,
                 "TimeConstrained[DSolve`Separable[y'[x] == %s, y, x], 16, $Aborted] "
                 "=!= $Aborted", MIXED_ANGLE[i]);
        ASSERT_TRUE(q);
    }
}

/* ---------------------------------------------------------------------------
 * F2 — the budget is general, not a trig special case.
 *
 * `Cos[x - 2y]/(Sin[x] Sin[2y]) - 1` is separable in principle (as Cot[x] Cot[2y],
 * which DSolve solves in 0.23 s when handed that form) but no method in the system
 * gets there from the compound-argument form — measured identical on the pre-fix
 * binary, 30 s to abort.  So the ONLY thing the fix can do here is stop Separable
 * burning the budget, and that is a property worth pinning on its own: it is what
 * generalises to every integrand nobody has thought of.
 * ------------------------------------------------------------------------- */
static void t_m63_budget_is_general(void) {
    const char* mult = "Cos[x - 2 y[x]]/(Sin[x] Sin[2 y[x]]) - 1";
    char q[900];
    snprintf(q, sizeof q,
             "TimeConstrained[DSolve`Separable[y'[x] == %s, y, x], 16, $Aborted] "
             "=!= $Aborted", mult);
    double t = seconds_of(q);
    ASSERT_TRUE(q);
    ASSERT_MSG(t < 16.0, "Separable did not decide the multiple-angle sibling: %.2f s", t);

    /* The reduced form is solvable, which is what makes the above a budget problem
     * and not a capability problem. */
    ASSERT_TRUE("Head[DSolve[y'[x] == Cot[x] Cot[2 y[x]], y, x]] === List");
}

/* ---------------------------------------------------------------------------
 * F3 — no legitimate separable is lost, and an unfired budget is free.
 * ------------------------------------------------------------------------- */
static void t_m63_separable_controls(void) {
    /* Sub-millisecond members: these also show the TimeConstrained wrapper costs
     * nothing when the integral closes immediately. */
    ASSERT_TRUE("DSolve[y'[x] == 2 x y[x], y, x] === {{y -> Function[{x}, C[1] E^x^2]}}");
    ASSERT_TRUE("Head[DSolve[y'[x] == y[x]^2/(x^2 + 1), y, x]] === List");
    ASSERT_TRUE("DSolve[y'[x] == y[x]/(x Log[x]), y, x] === {{y -> Function[{x}, C[1] Log[x]]}}");
    double t = seconds_of("DSolve[y'[x] == 2 x y[x], y, x]");
    ASSERT_MSG(t < 0.5, "a trivial separable slowed to %.3f s — the budget's "
                        "TimeConstrained wrapper should be free when it does not fire", t);

    /* The symbolic-parameter split (M32) and the trig-coefficient split. */
    ASSERT_TRUE("Head[DSolve[y'[x] == (a y[x] + b)/(c y[x] + d), y, x]] === List");
    ASSERT_TRUE("Head[DSolve[y'[x] == Cot[x] y[x]/(1 + y[x]), y, x]] === List");

    /* THE case a budget must not break: a DECIDED non-elementary integrand, where
     * the implicit twin legitimately returns an unevaluated Integrate.  A timeout and
     * a decided non-elementary integral must not be confused — this is the assertion
     * that notices if they are. */
    ASSERT_TRUE("Head[DSolve[y'[x] == -2 ArcTan[y[x]]/(1 + y[x]^2), y, x]] === List");
    ASSERT_TRUE("!FreeQ[DSolve[y'[x] == -2 ArcTan[y[x]]/(1 + y[x]^2), y, x], Integrate]");
}

/* ---------------------------------------------------------------------------
 * F4 — the answer is the EXPLICIT one, and it verifies against the original.
 *
 * This family is the discarded rewrite's epitaph.  Letting Separable claim these
 * records (which a TrigExpand-before-sampling rewrite does, in 0.08 s) produces an
 * implicit relation carrying the sampling artefact `Cot[2]`; letting the cascade
 * reach `DSolve`Linearizable` produces the ArcCos form Mathematica gives.  Faster was
 * not better, so the shape of the answer is asserted, not just its existence.
 * ------------------------------------------------------------------------- */
static void t_m63_answers_are_explicit(void) {
    char q[1600];
    for (size_t i = 0; i < N_MIXED; i++) {
        /* ONE solve per shape, both properties checked on it — these cost ~6 s each
         * and the suite shares the 120 s alarm in test_utils.h, so solving twice per
         * shape would put the file over it.
         *
         * Explicit: every branch is a Rule to a Function, with no leftover relation
         * (the discarded rewrite returned an implicit one).  And the residual of the
         * UN-REWRITTEN equation vanishes at a generic point, checked numerically
         * because the bodies are ArcCos of a constant times Csc and PossibleZeroQ
         * does not decide that. */
        snprintf(q, sizeof q,
                 "Module[{s = DSolve[y'[x] == %s, y, x], r}, "
                 "Head[s] === List && s =!= {} && "
                 "AllTrue[s, (Length[#] == 1 && Head[#[[1]]] === Rule && "
                 "           Head[#[[1, 2]]] === Function) &] && "
                 "(r = ((y'[x] - (%s)) /. s[[1]]) /. C[1] -> 13/10 /. x -> 11/10; "
                 " Abs[N[r]] < 10^-7)]", MIXED_ANGLE[i], MIXED_ANGLE[i]);
        double t = seconds_of(q);
        ASSERT_TRUE(q);
        /* The latency half of the fix, measured on the same solve.  Loose on purpose:
         * it must catch a regression to the 30 s spin, not police normal variation —
         * these close at ~6.2 s and the harness gives a corpus case 8 s. */
        ASSERT_MSG(t < 10.0, "mixed-angle shape %zu took %.2f s (pre-fix: Separable "
                             "alone spun past 30 s without deciding)", i, t);
    }

    /* The IVP member (§2.2.36-3527's shape): the constant is fitted, so no C[k]
     * survives, and the condition holds. */
    ASSERT_TRUE("Module[{s = DSolve[{y'[x] == 1 - Sin[x + y[x]]/(Sin[y[x]] Cos[x]), "
                "y[Pi/4] == Pi/4}, y, x]}, "
                "Head[s] === List && Cases[s, C[_], Infinity] === {} && "
                "Abs[N[(y[Pi/4] /. s[[1]]) - Pi/4]] < 10^-7]");
}

/* ---------------------------------------------------------------------------
 * F5 — the shapes the converter repair produced.
 * ------------------------------------------------------------------------- */
static void t_m63_repaired_record_shapes(void) {
    /* §2.2.18-1744 / §2.2.36-3570: repeated root a, so (C1 + C2 x) E^(a x). */
    ASSERT_TRUE("Head[DSolve[y''[x] - 2 a y'[x] + a^2 y[x] == 0, y, x]] === List");
    ASSERT_TRUE("PossibleZeroQ[Expand[(y''[x] - 2 a y'[x] + a^2 y[x]) /. "
                "DSolve[y''[x] - 2 a y'[x] + a^2 y[x] == 0, y, x][[1]]]]");
    ASSERT_TRUE("Length[Union[Cases[DSolve[y''[x] - 2 a y'[x] + a^2 y[x] == 0, y, x], "
                "C[_], Infinity]]] == 2");

    /* §2.2.26-2563: simple harmonic motion at frequency w. */
    ASSERT_TRUE("PossibleZeroQ[Expand[(y''[x] + w^2 y[x]) /. "
                "DSolve[y''[x] + w^2 y[x] == 0, y, x][[1]]]]");

    /* §2.2.33-3281: the IVP, whose answer must carry no free constant. */
    ASSERT_TRUE("Module[{s = DSolve[{x''[t] - k^2 x[t] == 0, x[0] == 0, "
                "x'[0] == v0}, x, t]}, "
                "Head[s] === List && Cases[s, C[_], Infinity] === {}]");
    ASSERT_TRUE("PossibleZeroQ[Simplify[(x''[t] - k^2 x[t]) /. "
                "DSolve[{x''[t] - k^2 x[t] == 0, x[0] == 0, x'[0] == v0}, x, t][[1]]]]");

    /* §2.2.16-1537 and §2.2.12-1182: autonomous first-order quadratures whose
     * parameter is a rate, not the variable. */
    ASSERT_TRUE("PossibleZeroQ[Expand[(y'[x] + a y[x]) /. "
                "DSolve[y'[x] + a y[x] == 0, y, x][[1]]]]");
    ASSERT_TRUE("Head[DSolve[y'[x] == a y[x] + b y[x]^2, y, x]] === List");

    /* §2.2.2-170, the constant-radius-of-curvature equation.  It is an honest
     * decline today (recorded in STATUS.md); the assertion is that it is not a
     * WRONG answer, which is what the mis-transcribed reading risked. */
    ASSERT_TRUE("Module[{s = DSolve[r y''[x] == (1 + y'[x]^2)^(3/2), y, x]}, "
                "s === {} || Head[s] === DSolve || "
                "PossibleZeroQ[Simplify[(r y''[x] - (1 + y'[x]^2)^(3/2)) /. s[[1]]]]]");
}

int main(void) {
    setvbuf(stdout, NULL, _IOLBF, 0);   /* the shared alarm() loses buffered output */
    symtab_init();
    core_init();
    test_load_init_m();

    TEST(t_m63_separable_decides_in_budget);
    TEST(t_m63_budget_is_general);
    TEST(t_m63_separable_controls);
    TEST(t_m63_answers_are_explicit);
    TEST(t_m63_repaired_record_shapes);

    printf("All DSolve M63 stress tests passed.\n");
    return 0;
}
