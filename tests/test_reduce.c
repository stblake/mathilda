/*
 * test_reduce.c
 *
 * Unit + stress tests for `Reduce` (src/solve/reduce.c) and its internal
 * logical normal-form layer (src/solve/reduce_form.{h,c}, reduce_atom.c).
 *
 * Coverage spans the shipped phases 0-5:
 *   0  front-end (arg parsing, True/False short-circuit, bad-variable) + DNF layer
 *   1  complete univariate polynomial equations over Complexes
 *   2  univariate real sign diagram (equations + inequalities over Reals)
 *   3  linear real systems (Fourier-Motzkin)
 *   4  parametric linear systems over Complexes (case analysis)
 *   5  Integers / Rationals
 *
 * The subsystem's cardinal invariant is SOUNDNESS OVER COMPLETENESS: an
 * undecidable sign/ordering, an unsupported construct, or a domain/shape the
 * shipped engines do not cover must leave `Reduce` UNEVALUATED (return itself),
 * never emit a wrong formula.  Whole test groups below (test_*_decline,
 * test_decline_soundness) exist only to pin that invariant.
 *
 * Black-box results are compared against FullForm strings, each captured from
 * the built binary and cross-checked for soundness.  White-box tests drive the
 * DNF layer directly on RAW (parse-only, un-evaluated) input, since the
 * evaluator pre-decides most constant logical combinations before Reduce ever
 * sees them.
 *
 * Three earlier "sound but imperfect" behaviours are now fixed and pinned as
 * regressions: `Reduce[a x == 0, x]` solves (was a decline; test_equations),
 * fully-determined Reals equation systems back-substitute to x == 2 && y == -1
 * (was a 4-bound box; test_linear_systems), and parametric conditions print in
 * minimal solved form a != 1 / a == 1/2 (test_parametric_systems).
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "reduce_form.h"
#include "test_utils.h"

/* ------------------------------------------------------------------ *
 *  Black-box: parse -> evaluate -> FullForm                           *
 * ------------------------------------------------------------------ */

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) { printf("FAIL: parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        free(got); expr_free(res); expr_free(e); ASSERT(0); return;
    }
    printf("PASS: %s -> %s\n", input, got);
    free(got); expr_free(res); expr_free(e);
}

/* Structural assertions for version-sensitive output (e.g. radical solutions,
 * whose exact FullForm is brittle): check a substring is present / absent. */
static void run_contains(const char* input, const char* needle) {
    Expr* e = parse_expression(input);
    if (!e) { printf("FAIL: parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (!strstr(got, needle)) {
        printf("FAIL: %s\n  expected to contain: %s\n  got: %s\n", input, needle, got);
        free(got); expr_free(res); expr_free(e); ASSERT(0); return;
    }
    printf("PASS(contains %s): %s\n", needle, input);
    free(got); expr_free(res); expr_free(e);
}

static void run_not_contains(const char* input, const char* needle) {
    Expr* e = parse_expression(input);
    if (!e) { printf("FAIL: parse: %s\n", input); ASSERT(0); return; }
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strstr(got, needle)) {
        printf("FAIL: %s\n  expected NOT to contain: %s\n  got: %s\n", input, needle, got);
        free(got); expr_free(res); expr_free(e); ASSERT(0); return;
    }
    printf("PASS(not-contains %s): %s\n", needle, input);
    free(got); expr_free(res); expr_free(e);
}

/* ------------------------------------------------------------------ *
 *  Phase 0 - constant / logical folding                              *
 * ------------------------------------------------------------------ */

/* Statements that fully decide collapse to True. */
static void test_decides_true(void) {
    run_test("Reduce[True, x]", "True");
    run_test("Reduce[1 < 2, x]", "True");
    run_test("Reduce[x == x, x]", "True");
    run_test("Reduce[2 <= 2, x]", "True");
    run_test("Reduce[Implies[True, 1 < 2], x]", "True");
    run_test("Reduce[! (3 < 2), x]", "True");
    /* Logical connectives over decidable operands. */
    run_test("Reduce[Xor[True, False], x]", "True");
    run_test("Reduce[Implies[1 < 2, 2 < 3], x]", "True");
    run_test("Reduce[!(3 < 2) && 1 < 2, x]", "True");
    run_test("Reduce[1 < 2 || 3 < 2, x]", "True");
    run_test("Reduce[(1 < 2 && 2 < 3) || 5 < 1, x]", "True");
}

/* Statements that fully decide collapse to False. */
static void test_decides_false(void) {
    run_test("Reduce[False, x]", "False");
    run_test("Reduce[3 < 2, x]", "False");
    run_test("Reduce[2 > 5, x]", "False");
    run_test("Reduce[1 < 2 && 3 < 2, x]", "False");
    run_test("Reduce[Xor[True, True], x]", "False");
    run_test("Reduce[x + 1 == x, x]", "False");
    run_test("Reduce[0 == 1, x]", "False");
}

/* ------------------------------------------------------------------ *
 *  Phase 1 - complete univariate equations over Complexes            *
 * ------------------------------------------------------------------ */

static void test_equations(void) {
    /* Parametric linear: the complete set keeps the degenerate a==0 branch. */
    run_test("Reduce[a x == b, x]",
             "Or[And[Unequal[a, 0], Equal[x, Times[Power[a, -1], b]]], "
             "And[Equal[a, 0], Equal[b, 0]]]");
    run_test("Reduce[a x + b == 0, x]",
             "Or[And[Unequal[a, 0], Equal[x, Times[-1, Power[a, -1], b]]], "
             "And[Equal[a, 0], Equal[b, 0]]]");
    /* a x == 0: the a == 0 branch's residual polynomial vanishes identically
     * (0 == 0, true for all x), so the degenerate branch is simply `a == 0`. */
    run_test("Reduce[a x == 0, x]",
             "Or[And[Unequal[a, 0], Equal[x, 0]], Equal[a, 0]]");
    /* Parametric right-hand side under a numeric-degree quadratic. */
    run_test("Reduce[x^2 == a, x]",
             "Or[Equal[x, Times[-1, Power[a, Rational[1, 2]]]], "
             "Equal[x, Power[a, Rational[1, 2]]]]");
    /* Numeric leading coefficient -> terminal, no case split. */
    run_test("Reduce[2 x == 6, x]", "Equal[x, 3]");
    run_test("Reduce[x/2 == 3, x]", "Equal[x, 6]");
    run_test("Reduce[x - 5 == 0, x]", "Equal[x, 5]");
    run_test("Reduce[x == 2, x]", "Equal[x, 2]");
    run_test("Reduce[x^2 == 4, x]", "Or[Equal[x, -2], Equal[x, 2]]");
    run_test("Reduce[x^2 == 4, x, Complexes]", "Or[Equal[x, -2], Equal[x, 2]]");
    /* A single-element variable list behaves like the bare variable. */
    run_test("Reduce[x^2 == 4, {x}]", "Or[Equal[x, -2], Equal[x, 2]]");
    run_test("Reduce[x^2 - 5 x + 6 == 0, x]", "Or[Equal[x, 2], Equal[x, 3]]");
    run_test("Reduce[(x - 1)(x - 2) == 0, x]", "Or[Equal[x, 1], Equal[x, 2]]");
    /* Default domain is Complexes: complex roots are kept. */
    run_test("Reduce[x^2 == -1, x]",
             "Or[Equal[x, Complex[0, -1]], Equal[x, Complex[0, 1]]]");
    run_test("Reduce[x^2 + 1 == 0, x]",
             "Or[Equal[x, Complex[0, -1]], Equal[x, Complex[0, 1]]]");
    /* Quadratic surd roots. */
    run_test("Reduce[x^2 - 2 == 0, x]",
             "Or[Equal[x, Times[-1, Power[2, Rational[1, 2]]]], "
             "Equal[x, Power[2, Rational[1, 2]]]]");
    /* Cubic: rational-power radical form. */
    run_test("Reduce[x^3 - 1 == 0, x]",
             "Or[Equal[x, 1], Equal[x, Times[-1, Power[-1, Rational[1, 3]]]], "
             "Equal[x, Power[-1, Rational[2, 3]]]]");
    /* Irreducible cubic: Root[] objects (real-first ascending index). */
    run_test("Reduce[x^3 - x - 1 == 0, x]",
             "Or[Equal[x, Root[Function[Plus[-1, Times[-1, Slot[1]], "
             "Power[Slot[1], 3]]], 1]], Equal[x, Root[Function[Plus[-1, "
             "Times[-1, Slot[1]], Power[Slot[1], 3]]], 2]], Equal[x, "
             "Root[Function[Plus[-1, Times[-1, Slot[1]], Power[Slot[1], 3]]], 3]]]");
    /* Roots of unity. */
    run_test("Reduce[x^4 == 1, x]",
             "Or[Equal[x, -1], Equal[x, 1], Equal[x, Complex[0, -1]], "
             "Equal[x, Complex[0, 1]]]");
    run_test("Reduce[x^5 == 1, x]",
             "Or[Equal[x, 1], Equal[x, Times[-1, Power[-1, Rational[1, 5]]]], "
             "Equal[x, Power[-1, Rational[2, 5]]], Equal[x, Times[-1, "
             "Power[-1, Rational[3, 5]]]], Equal[x, Power[-1, Rational[4, 5]]]]");
    /* The flagship 3-level split for a fully parametric quadratic. */
    run_test("Reduce[a x^2 + b x + c == 0, x]",
             "Or[And[Unequal[a, 0], Equal[x, Times[Rational[1, 2], Power[a, -1], "
             "Plus[Times[-1, b], Power[Plus[Power[b, 2], Times[-4, Times[a, c]]], "
             "Rational[1, 2]]]]]], And[Unequal[a, 0], Equal[x, Times[Rational[1, 2], "
             "Power[a, -1], Plus[Times[-1, b], Times[-1, Power[Plus[Power[b, 2], "
             "Times[-4, Times[a, c]]], Rational[1, 2]]]]]]], And[Equal[a, 0], "
             "Unequal[b, 0], Equal[x, Times[-1, Power[b, -1], c]]], "
             "And[Equal[a, 0], Equal[b, 0], Equal[c, 0]]]");
}

/* Phase-1 soundness boundary: Complexes routes only on a single equation atom;
 * anything else (!=, conjunction/disjunction of equations, transcendental,
 * the a x==0 asymmetry, no-variable form) must NOT invent an answer. */
static void test_equations_decline(void) {
    /* !=  over Complexes with no equation to pin the variety -> declines. */
    run_test("Reduce[x != 0, x]", "Reduce[Unequal[x, 0], x]");
    run_test("Reduce[x^2 != 4, x]", "Reduce[Unequal[Power[x, 2], 4], x]");
    /* A nonlinear conjunction / disjunction of equations is a zero-dimensional
     * system: the linear engine declines and reduce_zerodim solves it exactly.
     * (x^2==4 && x^3==8 has the single common root 2; the Or unions its cases.) */
    run_test("Reduce[x^2 == 4 && x^3 == 8, x]", "Equal[x, 2]");
    run_test("Reduce[x^2 == 4 || x == 5, x]",
             "Or[Equal[x, -2], Equal[x, 2], Equal[x, 5]]");
    /* No variable argument. */
    run_test("Reduce[x^2 == 4]", "Reduce[Equal[Power[x, 2], 4]]");
    /* Transcendental equation: Solve returns a bare relation, Reduce passes it
     * through unchanged (a distinct outcome from both solved and Reduce[...]). */
    run_test("Reduce[Sin[x] == 0, x]", "Equal[Sin[x], 0]");
}

/* ------------------------------------------------------------------ *
 *  Phase 2 - univariate real sign diagram                            *
 * ------------------------------------------------------------------ */

static void test_real_inequalities(void) {
    run_test("Reduce[x^2 > 1, x, Reals]",  "Or[Less[x, -1], Greater[x, 1]]");
    run_test("Reduce[x^2 >= 1, x, Reals]", "Or[LessEqual[x, -1], GreaterEqual[x, 1]]");
    run_test("Reduce[x^2 < 1, x, Reals]",  "Inequality[-1, Less, x, Less, 1]");
    run_test("Reduce[x^2 <= 1, x, Reals]", "Inequality[-1, LessEqual, x, LessEqual, 1]");
    /* Cofinite: complement of finitely many points. */
    run_test("Reduce[x^2 != 1, x, Reals]", "And[Unequal[x, -1], Unequal[x, 1]]");
    /* No real breakpoints / degenerate ranges -> the whole line decides. */
    run_test("Reduce[x^2 < 0, x, Reals]", "False");
    run_test("Reduce[x^2 >= 0, x, Reals]", "True");
    run_test("Reduce[x^2 <= 0, x, Reals]", "Equal[x, 0]");
    run_test("Reduce[x^2 + 1 > 0, x, Reals]", "True");
    run_test("Reduce[x^2 + 1 < 0, x, Reals]", "False");
    /* Odd multiplicity: a single sign change. */
    run_test("Reduce[x^3 > 0, x, Reals]", "Greater[x, 0]");
    /* Algebraic (radical) breakpoints, ordered and signed via the qqbar oracle. */
    run_test("Reduce[x^2 - 2 > 0, x, Reals]",
             "Or[Less[x, Times[-1, Power[2, Rational[1, 2]]]], "
             "Greater[x, Power[2, Rational[1, 2]]]]");
    run_test("Reduce[x^3 > 2, x, Reals]",
             "Greater[x, Power[2, Rational[1, 3]]]");
    run_test("Reduce[x^2 < 2, x, Reals]",
             "Inequality[Times[-1, Power[2, Rational[1, 2]]], Less, x, Less, "
             "Power[2, Rational[1, 2]]]");
    /* Multiple simple roots: alternating sign cells. */
    run_test("Reduce[(x - 1)(x - 2)(x - 3) > 0, x, Reals]",
             "Or[Inequality[1, Less, x, Less, 2], Greater[x, 3]]");
    run_test("Reduce[x^2 - 5 x + 6 > 0, x, Reals]",
             "Or[Less[x, 2], Greater[x, 3]]");
    /* Conjunctions of bounds. */
    run_test("Reduce[x > 0 && x < 1, x, Reals]", "Inequality[0, Less, x, Less, 1]");
    run_test("Reduce[x < 1 && x > 0, x, Reals]", "Inequality[0, Less, x, Less, 1]");
    run_test("Reduce[x^2 > 1 && x < 5, x, Reals]",
             "Or[Less[x, -1], Inequality[1, Less, x, Less, 5]]");
    run_test("Reduce[x >= 2 && x >= 5, x, Reals]", "GreaterEqual[x, 5]");
    run_test("Reduce[x^2 > 1 && x^2 < 9, x, Reals]",
             "Or[Inequality[-3, Less, x, Less, -1], Inequality[1, Less, x, Less, 3]]");
    run_test("Reduce[x^2 <= 1 && x >= 0, x, Reals]",
             "Inequality[0, LessEqual, x, LessEqual, 1]");
    run_test("Reduce[x >= 1 && x <= 3, x, Reals]",
             "Inequality[1, LessEqual, x, LessEqual, 3]");
    /* Disjunctions. */
    run_test("Reduce[x^2 > 1 || x < -5, x, Reals]",
             "Or[Less[x, -1], Greater[x, 1]]");
    run_test("Reduce[x^2 <= 1 || x >= 3, x, Reals]",
             "Or[Inequality[-1, LessEqual, x, LessEqual, 1], GreaterEqual[x, 3]]");
    run_test("Reduce[!(x^2 < 1), x, Reals]",
             "Or[LessEqual[x, -1], GreaterEqual[x, 1]]");
    /* Mixed equation + inequality over Reals (works, unlike Complexes). */
    run_test("Reduce[x^2 == 4 && x > 0, x, Reals]", "Equal[x, 2]");
    run_test("Reduce[x^2 == 4 || x == 5, x, Reals]",
             "Or[Equal[x, -2], Equal[x, 2], Equal[x, 5]]");
    /* !=  over Reals is a satisfiable condition, not a decline. */
    run_test("Reduce[x != 0, x, Reals]", "Unequal[x, 0]");
    /* Equations over Reals fall out of the same sign diagram. */
    run_test("Reduce[x^2 == 4, x, Reals]", "Or[Equal[x, -2], Equal[x, 2]]");
    run_test("Reduce[x^2 - 2 == 0, x, Reals]",
             "Or[Equal[x, Times[-1, Power[2, Rational[1, 2]]]], "
             "Equal[x, Power[2, Rational[1, 2]]]]");
}

/* ------------------------------------------------------------------ *
 *  Rational-function inequalities over the Reals: the sign diagram    *
 *  adds the poles (roots of the denominator) as breakpoints and signs *
 *  p/q as sign(p)*sign(q), so clearing the denominator never flips    *
 *  the sense.  Poles are excluded (p/q is undefined there).           *
 * ------------------------------------------------------------------ */

static void test_rational_inequalities(void) {
    /* Explicit-domain forms. */
    run_test("Reduce[1/x < 1, x, Reals]", "Or[Less[x, 0], Greater[x, 1]]");
    run_test("Reduce[7/x < 22, x, Reals]",
             "Or[Less[x, 0], Greater[x, Rational[7, 22]]]");
    /* p/q >= 0 excludes the pole: x > 0, not x >= 0. */
    run_test("Reduce[1/x >= 0, x, Reals]", "Greater[x, 0]");
    /* p/q != 0 still excludes the pole: 1/x != 0 is x != 0. */
    run_test("Reduce[1/x != 0, x, Reals]", "Unequal[x, 0]");
    /* p/q == 0 with a nonzero-constant numerator is unsatisfiable. */
    run_test("Reduce[1/x == 0, x, Reals]", "False");
    /* Two poles from a quadratic denominator. */
    run_test("Reduce[1/(x^2 - 1) < 0, x, Reals]",
             "Inequality[-1, Less, x, Less, 1]");
    /* Denominator with no real root: the pole set is empty. */
    run_test("Reduce[1/(x^2 + 1) < 0, x, Reals]", "False");
    run_test("Reduce[(x - 1)/(x - 2) > 0, x, Reals]",
             "Or[Less[x, 1], Greater[x, 2]]");
    /* The reported case: a chained rational inequality with no explicit
     * domain defaults to the Reals and is fully solved. */
    run_test("Reduce[-5 < 3 x + 7/x <= 22, x]",
             "Inequality[Rational[1, 3], LessEqual, x, LessEqual, 7]");
}

/* ------------------------------------------------------------------ *
 *  Phase 3 - multivariate linear systems over Reals (Fourier-Motzkin)*
 * ------------------------------------------------------------------ */

static void test_linear_systems(void) {
    /* The plan's flagship: a triangular description of the feasible region. */
    run_test("Reduce[x + y < 1 && x > 0 && y > 0, {x, y}, Reals]",
             "And[Inequality[0, Less, x, Less, 1], "
             "Inequality[0, Less, y, Less, Plus[1, Times[-1, x]]]]");
    /* Infeasible system -> False. */
    run_test("Reduce[x > 1 && x < 0 && y > 0, {x, y}, Reals]", "False");
    /* An equation is re-detected as `==` in the triangular output; a variable
     * pinned to a range stays symbolic (y == 1 - x, x not determined). */
    run_test("Reduce[x + y == 1 && x > 0, {x, y}, Reals]",
             "And[Greater[x, 0], Equal[y, Plus[1, Times[-1, x]]]]");
    /* A fully-determined equation system back-substitutes the pinned variable:
     * x == 2 propagates into y == 1 - x, giving y == -1 (not a 4-bound box). */
    run_test("Reduce[x + y == 1 && x - y == 3, {x, y}, Reals]",
             "And[Equal[x, 2], Equal[y, -1]]");
    /* Rational bound coefficients. */
    run_test("Reduce[2 x + 3 y <= 6 && x >= 0 && y >= 0, {x, y}, Reals]",
             "And[Inequality[0, LessEqual, x, LessEqual, 3], "
             "Inequality[0, LessEqual, y, LessEqual, Plus[2, Times[Rational[-2, 3], x]]]]");
    run_test("Reduce[2 x + 3 y <= 6 && x >= 0, {x, y}, Reals]",
             "And[GreaterEqual[x, 0], LessEqual[y, Plus[2, Times[Rational[-2, 3], x]]]]");
    /* A free variable is simply omitted. */
    run_test("Reduce[x > 0, {x, y}, Reals]", "Greater[x, 0]");
    /* Disjunction: each conjunct solved and OR-ed. */
    run_test("Reduce[x < 0 || x > 1, {x, y}, Reals]",
             "Or[Less[x, 0], Greater[x, 1]]");
}

/* ------------------------------------------------------------------ *
 *  Phase 4 - parametric linear systems over Complexes                *
 * ------------------------------------------------------------------ */

static void test_parametric_systems(void) {
    /* Numeric determined system. */
    run_test("Reduce[x + y == 3 && x - y == 1, {x, y}]",
             "And[Equal[y, 1], Equal[x, 2]]");
    /* Underdetermined: a variable stays free. */
    run_test("Reduce[x + y == 1, {x, y}]", "Equal[x, Plus[1, Times[-1, y]]]");
    /* Parametric square system: genericity condition + Cramer solution.  The
     * condition prints in minimal solved form `a != 1` (not `1 - a != 0`). */
    run_test("Reduce[a x + y == 1 && x + y == 0, {x, y}]",
             "And[Unequal[a, 1], Equal[x, Power[Plus[-1, a], -1]], "
             "Equal[y, Times[-1, Power[Plus[-1, a], -1]]]]");
    /* Overdetermined: a consistency condition on the parameter, minimal form
     * `a == 1/2` (not `-1 + 2 a == 0`). */
    run_test("Reduce[a x == 1 && x == 2, {x}]",
             "And[Equal[a, Rational[1, 2]], Equal[x, 2]]");
    run_test("Reduce[a x == 1 && x == 2, x]",
             "And[Equal[a, Rational[1, 2]], Equal[x, 2]]");
    /* Three variables. */
    run_test("Reduce[x + y + z == 6 && x - y == 0 && z == 2, {x, y, z}]",
             "And[Equal[z, 2], Equal[y, 2], Equal[x, 2]]");
}

/* Phase-4 boundary and beyond: a linear system is grafted by the linear engine;
 * a ZERO-DIMENSIONAL nonlinear system is solved exactly by reduce_zerodim; only
 * a genuinely POSITIVE-dimensional system (fewer independent equations than
 * variables) still declines. */
static void test_parametric_systems_decline(void) {
    /* Zero-dimensional: x+y==3, xy==1 -> the two roots (3 +- Sqrt[5])/2. */
    run_contains("Reduce[x y == 1 && x + y == 3, {x, y}]", "Power[5, Rational[1, 2]]");
    /* Zero-dimensional: x==0 then y^2==1 -> the two points (0, -+1). */
    run_test("Reduce[x^2 + y^2 == 1 && x == 0, {x, y}]",
             "Or[And[Equal[x, 0], Equal[y, -1]], And[Equal[x, 0], Equal[y, 1]]]");
    /* A single non-linear equation in two variables is positive-dimensional
     * (a curve): zero-dim declines (parametric), CAD over Complexes is unwired,
     * so Reduce stays unevaluated. */
    run_test("Reduce[x^2 + y^2 == 1, {x, y}]",
             "Reduce[Equal[Plus[Power[x, 2], Power[y, 2]], 1], List[x, y]]");
}

/* ------------------------------------------------------------------ *
 *  Phase 5 - Integers / Rationals                                    *
 * ------------------------------------------------------------------ */

static void test_integer_domain(void) {
    /* Equation -> finite solution set. */
    run_test("Reduce[x^2 == 4, x, Integers]", "Or[Equal[x, -2], Equal[x, 2]]");
    run_test("Reduce[x^2 == 2, x, Integers]", "False");
    run_test("Reduce[2 x == 5, x, Integers]", "False");
    run_test("Reduce[2 x + 1 == 5, x, Integers]", "Equal[x, 2]");
    run_test("Reduce[2 x == 4, x, Integers]", "Equal[x, 2]");
    /* Bounded inequality -> integer enumeration (Solve declines, we fall back). */
    run_test("Reduce[0 <= x <= 5, x, Integers]",
             "Or[Equal[x, 0], Equal[x, 1], Equal[x, 2], Equal[x, 3], "
             "Equal[x, 4], Equal[x, 5]]");
    run_test("Reduce[x^2 <= 4, x, Integers]",
             "Or[Equal[x, -2], Equal[x, -1], Equal[x, 0], Equal[x, 1], Equal[x, 2]]");
    run_test("Reduce[x^2 < 5 && x > -3, x, Integers]",
             "Or[Equal[x, -2], Equal[x, -1], Equal[x, 0], Equal[x, 1], Equal[x, 2]]");
    run_test("Reduce[x >= 1 && x <= 3, x, Integers]",
             "Or[Equal[x, 1], Equal[x, 2], Equal[x, 3]]");
    run_test("Reduce[x^2 < 10 && x > 0, x, Integers]",
             "Or[Equal[x, 1], Equal[x, 2], Equal[x, 3]]");
    run_test("Reduce[1 <= x <= 3, x, Integers]",
             "Or[Equal[x, 1], Equal[x, 2], Equal[x, 3]]");
    /* Unbounded inequalities -> one-sided rays (a satisfied tail beyond the
     * extreme root emits x<=k / x>=k rather than declining). */
    run_test("Reduce[x > 0, x, Integers]", "GreaterEqual[x, 1]");
    run_test("Reduce[x < 5, x, Integers]", "LessEqual[x, 4]");
    run_test("Reduce[x^2 > 1, x, Integers]",
             "Or[LessEqual[x, -2], GreaterEqual[x, 2]]");
    /* Bounded system with an equation. */
    run_test("Reduce[x + y == 5 && x > 0 && y > 0, {x, y}, Integers]",
             "Or[And[Equal[x, 1], Equal[y, 4]], And[Equal[x, 2], Equal[y, 3]], "
             "And[Equal[x, 3], Equal[y, 2]], And[Equal[x, 4], Equal[y, 1]]]");
    /* Parametric linear Diophantine -> a C[k] family with an Element condition. */
    run_test("Reduce[x + y == 3, {x, y}, Integers]",
             "And[Element[C[1], Integers], Equal[x, C[1]], "
             "Equal[y, Plus[3, Times[-1, C[1]]]]]");
    run_test("Reduce[2 x + 3 y == 1, {x, y}, Integers]",
             "And[Element[C[1], Integers], Equal[x, Plus[-1, Times[3, C[1]]]], "
             "Equal[y, Plus[1, Times[-2, C[1]]]]]");
}

static void test_rational_domain(void) {
    run_test("Reduce[x^2 == 4, x, Rationals]", "Or[Equal[x, -2], Equal[x, 2]]");
    run_test("Reduce[x^2 == 2, x, Rationals]", "False");
}

/* ------------------------------------------------------------------ *
 *  Domain defaulting: an ordering inequality (< <= > >=) with no      *
 *  explicit domain is solved over the Reals, since ordering is        *
 *  undefined over the default Complexes (matches Mathematica).        *
 *  Equations and Unequal (!=) stay on the Complexes default.          *
 * ------------------------------------------------------------------ */

static void test_inequality_defaults_reals(void) {
    /* Bare univariate inequalities -> same result as the explicit ,Reals form. */
    run_test("Reduce[x > 0, x]", "Greater[x, 0]");
    run_test("Reduce[x^2 > 1, x]", "Or[Less[x, -1], Greater[x, 1]]");
    run_test("Reduce[-5 < 3 x + 7 <= 22, x]",
             "Inequality[-4, Less, x, LessEqual, 5]");
    /* Mixed equation + inequality: the inequality still triggers the Reals
     * default, so this solves rather than declining. */
    run_test("Reduce[x^2 == 4 && x > 0, x]", "Equal[x, 2]");
    /* Multivariate real system, no domain -> CAD/Fourier-Motzkin over Reals. */
    run_test("Reduce[x > 1 && x < 0 && y > 0, {x, y}]", "False");
    /* An equation alone keeps the Complexes default (no over-broadening). */
    run_test("Reduce[x^2 == 4, x]", "Or[Equal[x, -2], Equal[x, 2]]");
    /* Unequal (!=) is meaningful over Complexes, so it does NOT force Reals:
     * this stays on the Complexes path and is left unevaluated as before. */
    run_test("Reduce[x != 1, x]", "Reduce[Unequal[x, 1], x]");
}

/* ------------------------------------------------------------------ *
 *  Soundness net: inputs the shipped engines do NOT cover must leave *
 *  Reduce unevaluated -- never a wrong (or guessed) formula.         *
 * ------------------------------------------------------------------ */

static void test_decline_soundness(void) {
    /* Unsupported / unknown domains. */
    run_test("Reduce[x^2 == 4, x, Booleans]",
             "Reduce[Equal[Power[x, 2], 4], x, Booleans]");
    run_test("Reduce[x^2 == 4, x, GaussianIntegers]",
             "Reduce[Equal[Power[x, 2], 4], x, GaussianIntegers]");
    /* Element condition without a domain argument. */
    run_test("Reduce[Element[x, Integers] && x^2 < 5, x]",
             "Reduce[And[Element[x, Integers], Less[Power[x, 2], 5]], x]");
    /* (Multivariate nonlinear over Reals and parametric-linear real systems are
     * now solved by the CAD engine -- see test_cad_real.) */
    /* Free parameter in a univariate real inequality/equation. */
    run_test("Reduce[a x^2 == 1, x, Reals]",
             "Reduce[Equal[Times[a, Power[x, 2]], 1], x, Reals]");
    /* Invalid variable specification -> Reduce::ivar (stderr) + unevaluated. */
    run_test("Reduce[x == 1, 5]", "Reduce[Equal[x, 1], 5]");
    run_test("Reduce[x == 1, {}]", "Reduce[Equal[x, 1], List[]]");
}

/* ------------------------------------------------------------------ *
 *  White-box: drive the DNF layer on raw (un-evaluated) input         *
 * ------------------------------------------------------------------ */

/* Parse `input` WITHOUT evaluating, normalise to DNF over variable x, simplify,
 * emit, and compare the FullForm against `expected`. */
static void wb_test(const char* input, const char* expected) {
    Expr* e  = parse_expression(input);
    Expr* vx = parse_expression("x");   /* interned symbol x */
    ASSERT(e && vx);
    Expr* vars[1] = { vx };

    bool ok = true;
    RForm* f = reduce_form_from_expr(e, vars, 1, &ok);
    ASSERT(ok);
    rform_simplify(f, vars, 1);
    Expr* out = rform_to_expr(f, vars, 1);
    char* got = expr_to_string_fullform(out);
    if (strcmp(got, expected) != 0) {
        printf("WB FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        free(got); expr_free(out); rform_free(f); expr_free(e); expr_free(vx);
        ASSERT(0); return;
    }
    printf("WB PASS: %s -> %s\n", input, got);
    free(got); expr_free(out); rform_free(f); expr_free(e); expr_free(vx);
}

/* Constant atoms are decided; >/>= are canonicalised to </<=. */
static void test_wb_constant_atoms(void) {
    wb_test("3 < 2", "False");
    wb_test("5 > 2", "True");
    wb_test("2 <= 2", "True");
    wb_test("3 != 3", "False");
}

/* And/Or fold constants correctly around a symbolic atom. */
static void test_wb_logic_fold(void) {
    wb_test("3 < 2 || 5 > 2", "True");        /* False v True  = True  */
    wb_test("1 < 2 && 3 < 2", "False");       /* True  ^ False = False */
    /* A true constant conjunct drops, leaving the symbolic atom. */
    wb_test("5 > 2 && x^2 == 4", "Equal[Plus[-4, Power[x, 2]], 0]");
    /* A false constant conjunct kills the whole conjunction. */
    wb_test("3 < 2 && x^2 == 4", "False");
}

/* Assert that the DNF layer rejects an unsupported construct (ok == false),
 * the *ok=false path the front-end turns into an unevaluated Reduce[...]. */
static void wb_unsupported(const char* input) {
    Expr* e  = parse_expression(input);
    Expr* vx = parse_expression("x");
    ASSERT(e && vx);
    Expr* vars[1] = { vx };
    bool ok = true;
    RForm* f = reduce_form_from_expr(e, vars, 1, &ok);
    if (ok) {
        printf("WB FAIL: %s\n  expected: ok=false (unsupported)\n  got:      ok=true\n",
               input);
        if (f) rform_free(f); expr_free(e); expr_free(vx); ASSERT(0); return;
    }
    printf("WB PASS: %s -> ok=false (unsupported)\n", input);
    if (f) rform_free(f); expr_free(e); expr_free(vx);
}

/* A single symbolic atom canonicalises to `poly REL 0`.  >/>= are flipped to
 * </<= by swapping operands; equations clear a fully-cancelled denominator. */
static void test_wb_atom_emit(void) {
    wb_test("x^2 == 4", "Equal[Plus[-4, Power[x, 2]], 0]");
    wb_test("x != 1", "Unequal[Plus[-1, x], 0]");
    wb_test("x > 1", "Less[Plus[1, Times[-1, x]], 0]");
    wb_test("x >= 1", "LessEqual[Plus[1, Times[-1, x]], 0]");
    /* Swapped-sides equation. */
    wb_test("4 == x^2", "Equal[Plus[4, Times[-1, Power[x, 2]]], 0]");
    wb_test("2 - x == 0", "Equal[Plus[2, Times[-1, x]], 0]");
    wb_test("-x + 3 > 0", "Less[Plus[-3, x], 0]");
    /* Denominator clearing is sound for an equation: 1/x == 0 -> 1 == 0. */
    wb_test("1/x == 0", "False");
}

/* Higher-level logical connectives expand into the DNF of relational atoms:
 * Implies -> De Morgan, Xor -> distribution, chained Inequality splits, Not
 * pushes through And/Or/atoms. */
static void test_wb_logic_expand(void) {
    wb_test("Implies[x^2 == 1, x == 1]",
            "Or[Unequal[Plus[-1, Power[x, 2]], 0], Equal[Plus[-1, x], 0]]");
    wb_test("Xor[x == 1, x == 2]",
            "Or[And[Equal[Plus[-1, x], 0], Unequal[Plus[-2, x], 0]], "
            "And[Unequal[Plus[-1, x], 0], Equal[Plus[-2, x], 0]]]");
    /* A chained Inequality splits into a conjunction of two atoms. */
    wb_test("1 < x < 3",
            "And[Less[Plus[1, Times[-1, x]], 0], Less[Plus[-3, x], 0]]");
    /* Not of an equation atom. */
    wb_test("!(x^2 == 4)", "Unequal[Plus[-4, Power[x, 2]], 0]");
    /* De Morgan over a conjunction / disjunction. */
    wb_test("!(x == 1 && x == 2)",
            "Or[Unequal[Plus[-1, x], 0], Unequal[Plus[-2, x], 0]]");
    wb_test("!(x == 1 || x == 2)",
            "And[Unequal[Plus[-1, x], 0], Unequal[Plus[-2, x], 0]]");
}

/* Constructs the DNF layer does not support -> ok=false (Reduce stays
 * unevaluated).  Negating an Element predicate is one such case. */
static void test_wb_unsupported(void) {
    wb_unsupported("Not[Element[x, Integers]]");
}

/* ------------------------------------------------------------------ *
 *  Phase 6: multivariate nonlinear inequalities over the reals (CAD)  *
 * ------------------------------------------------------------------ */

/* The exact FullForm strings below are pinned from the built binary -- they are
 * a correct, complete description of each solution set (verified semantically by
 * the corpus sample-point oracle), NOT necessarily Mathematica's cosmetic form
 * (radicals stay in Solve's `(1/2) Sqrt[4 - 4 x^2]` shape).  The boundary-merge
 * pass fuses a closed region's limit sections into the adjacent sector, so a
 * non-strict region closes its x-range (`-1 <= x <= 1`) instead of listing the
 * boundary points separately. */
static void test_cad_real(void) {
    /* Unsatisfiable / universal statements collapse to False / True. */
    run_test("Reduce[x^2 + y^2 < 0, {x, y}, Reals]", "False");
    run_test("Reduce[x^2 + y^2 >= 0, {x, y}, Reals]", "True");

    /* Product of two lines: the two open quadrant-like sectors. */
    run_test("Reduce[x y > 0, {x, y}, Reals]",
             "Or[And[Less[x, 0], Less[y, 0]], And[Greater[x, 0], Greater[y, 0]]]");

    /* Closed unit disk: the boundary sections x == +/-1 merge into the sector,
     * closing the x-range to -1 <= x <= 1. */
    run_test("Reduce[x^2 + y^2 <= 1, {x, y}, Reals]",
             "And[Inequality[-1, LessEqual, x, LessEqual, 1], "
             "Inequality[Times[Rational[-1, 2], Power[Plus[4, Times[-4, Power[x, 2]]], Rational[1, 2]]], "
             "LessEqual, y, LessEqual, "
             "Times[Rational[1, 2], Power[Plus[4, Times[-4, Power[x, 2]]], Rational[1, 2]]]]]");

    /* Closed unit circle (the curve): the two symbolic arcs, x-range closed. */
    run_test("Reduce[x^2 + y^2 == 1, {x, y}, Reals]",
             "And[Inequality[-1, LessEqual, x, LessEqual, 1], "
             "Or[Equal[y, Times[Rational[-1, 2], Power[Plus[4, Times[-4, Power[x, 2]]], Rational[1, 2]]]], "
             "Equal[y, Times[Rational[1, 2], Power[Plus[4, Times[-4, Power[x, 2]]], Rational[1, 2]]]]]]");

    /* Open region (strict <): sections vanish, so the answer merges cleanly. */
    run_test("Reduce[x^2 - y^2 > 1, {x, y}, Reals]",
             "Or[And[Less[x, -1], "
             "Inequality[Times[-1, Power[Plus[-1, Power[x, 2]], Rational[1, 2]]], Less, y, Less, "
             "Power[Plus[-1, Power[x, 2]], Rational[1, 2]]]], "
             "And[Greater[x, 1], "
             "Inequality[Times[-1, Power[Plus[-1, Power[x, 2]], Rational[1, 2]]], Less, y, Less, "
             "Power[Plus[-1, Power[x, 2]], Rational[1, 2]]]]]");

    /* A pure-x conjunct prunes half the base cells (partial CAD); the closed
     * half-disk merges to 0 <= x <= 1. */
    run_test("Reduce[x^2 + y^2 <= 1 && x >= 0, {x, y}, Reals]",
             "And[Inequality[0, LessEqual, x, LessEqual, 1], "
             "Inequality[Times[Rational[-1, 2], Power[Plus[4, Times[-4, Power[x, 2]]], Rational[1, 2]]], "
             "LessEqual, y, LessEqual, "
             "Times[Rational[1, 2], Power[Plus[4, Times[-4, Power[x, 2]]], Rational[1, 2]]]]]");

    /* A parametric coefficient the linear (Fourier-Motzkin) engine cannot take,
     * but whose single base cell is sign-invariant, is solved by CAD. */
    run_test("Reduce[a x + y < 1, {x, y}, Reals]",
             "Less[y, Plus[1, Times[-1, Times[a, x]]]]");

    /* A nominally-2-var problem whose second variable is absent is really the
     * univariate sign diagram. */
    run_test("Reduce[x^2 > 1, {x, y}, Reals]",
             "Or[Less[x, -1], Greater[x, 1]]");

    /* Nonlinear multivariate equation over Complexes still needs the (deferred)
     * equational CAD route: declines. */
    run_test("Reduce[x^2 + y^2 == 1, {x, y}]",
             "Reduce[Equal[Plus[Power[x, 2], Power[y, 2]], 1], List[x, y]]");
}

/* Phase 6d: n-variable (nu>=3) CAD over the Reals (recursive McCallum
 * projection + lift).  Output for STRICT inequalities is the clean nested form
 * (no boundary sections); CLOSED regions currently emit a correct, verbose Or of
 * cells (the n-D boundary merge that would close an outer range is deferred).
 * The semantic corpus (reduce_corpus.m) certifies the closed/verbose cases by
 * sample-point equivalence; here we pin the stable clean forms and the
 * True/False/decline invariants. */
static void test_cad_nvar(void) {
    /* Unsatisfiable / universal collapse. */
    run_test("Reduce[x^2 + y^2 + z^2 < 0, {x, y, z}, Reals]", "False");
    run_test("Reduce[x^2 + y^2 + z^2 >= 0, {x, y, z}, Reals]", "True");

    /* Strict unit ball: one clean nested conjunction, bounds symbolic in the
     * outer variables. */
    run_test("Reduce[x^2 + y^2 + z^2 < 1, {x, y, z}, Reals]",
             "And[Inequality[-1, Less, x, Less, 1], "
             "Inequality[Times[Rational[-1, 2], Power[Plus[4, Times[-4, Power[x, 2]]], Rational[1, 2]]], "
             "Less, y, Less, "
             "Times[Rational[1, 2], Power[Plus[4, Times[-4, Power[x, 2]]], Rational[1, 2]]]], "
             "Inequality[Times[Rational[-1, 2], Power[Times[-4, Plus[-1, Power[x, 2], Power[y, 2]]], Rational[1, 2]]], "
             "Less, z, Less, "
             "Power[Plus[1, Times[-1, Power[x, 2]], Times[-1, Power[y, 2]]], Rational[1, 2]]]]");

    /* Closed unit ball: the n-D boundary merge (Stage B) closes each range to a
     * non-strict nested conjunction (-1 <= x <= 1 && ... && ...). */
    run_test("Reduce[x^2 + y^2 + z^2 <= 1, {x, y, z}, Reals]",
             "And[Inequality[-1, LessEqual, x, LessEqual, 1], "
             "Inequality[Times[Rational[-1, 2], Power[Plus[4, Times[-4, Power[x, 2]]], Rational[1, 2]]], "
             "LessEqual, y, LessEqual, "
             "Times[Rational[1, 2], Power[Plus[4, Times[-4, Power[x, 2]]], Rational[1, 2]]]], "
             "Inequality[Times[Rational[-1, 2], Power[Times[-4, Plus[-1, Power[x, 2], Power[y, 2]]], Rational[1, 2]]], "
             "LessEqual, z, LessEqual, "
             "Power[Plus[1, Times[-1, Power[x, 2]], Times[-1, Power[y, 2]]], Rational[1, 2]]]]");

    /* Open axis-aligned box: three independent one-variable ranges. */
    run_test("Reduce[x^2 < 1 && y^2 < 1 && z^2 < 1, {x, y, z}, Reals]",
             "And[Inequality[-1, Less, x, Less, 1], Inequality[-1, Less, y, Less, 1], "
             "Inequality[-1, Less, z, Less, 1]]");

    /* Sign product of three lines: the positive octants, factored by the sign
     * of x (the merge groups each x-sector's y/z sub-decomposition). */
    run_test("Reduce[x y z > 0, {x, y, z}, Reals]",
             "Or[And[Less[x, 0], Or[And[Less[y, 0], Greater[z, 0]], And[Greater[y, 0], Less[z, 0]]]], "
             "And[Greater[x, 0], Or[And[Less[y, 0], Less[z, 0]], And[Greater[y, 0], Greater[z, 0]]]]]");

    /* Product of two planes: two open wedges. */
    run_test("Reduce[(x - y) (y - z) > 0, {x, y, z}, Reals]",
             "Or[And[Less[y, x], Less[z, y]], And[Greater[y, x], Greater[z, y]]]");

    /* Irrational base breakpoint (+/-Sqrt[2]) is outside the v1 rational-fibre
     * regime: declines (Phase 6b), never a wrong answer. */
    run_test("Reduce[x^2 + y^2 + z^2 <= 2, {x, y, z}, Reals]",
             "Reduce[LessEqual[Plus[Power[x, 2], Power[y, 2], Power[z, 2]], 2], "
             "List[x, y, z], Reals]");
}

/* ------------------------------------------------------------------ *
 *  Phase 9: elementary real functions (radicals, Abs, Log,           *
 *  inverse-trig, Floor/Mod) over the Reals via the general sign      *
 *  diagram + preprocessing.  Sample-verified in reduce_corpus.m;     *
 *  a few representative outputs are pinned here.                     *
 * ------------------------------------------------------------------ */

static void test_real_functions(void) {
    /* Nested radicals: the equation is an identity on an interval. */
    run_test("Reduce[Sqrt[x + 3 - 4 Sqrt[x - 1]] + Sqrt[x + 8 - 6 Sqrt[x - 1]] == 1, x, Reals]",
             "Inequality[5, LessEqual, x, LessEqual, 10]");
    /* A squared polynomial factor contributes an isolated solution point. */
    run_test("Reduce[(2 x - 1)^2 (Sqrt[x + 4 - 4 Sqrt[x]] + Sqrt[x + 9 - 6 Sqrt[x]] - 1) == 0, x, Reals]",
             "Or[Equal[x, Rational[1, 2]], Inequality[4, LessEqual, x, LessEqual, 9]]");
    /* Nested Abs. */
    run_test("Reduce[Abs[Abs[x] - 2] + Abs[Abs[x] - 5] == 3, x, Reals]",
             "Or[Inequality[-5, LessEqual, x, LessEqual, -2], "
             "Inequality[2, LessEqual, x, LessEqual, 5]]");
    /* Abs numerator over a pole. */
    run_test("Reduce[(Abs[x + 3] + Abs[x - 3])/x == 2, x, Reals]", "GreaterEqual[x, 3]");
    /* Radical domain intersection. */
    run_test("Reduce[Sqrt[x^2 - 4] == Sqrt[x - 2] Sqrt[x + 2], x, Reals]", "GreaterEqual[x, 2]");
    /* Sqrt of a perfect square -> |.|. */
    run_test("Reduce[Sqrt[(x^2 - 4)^2] == 4 - x^2, x, Reals]",
             "Inequality[-2, LessEqual, x, LessEqual, 2]");
    /* Floor via defining inequalities. */
    run_test("Reduce[Floor[2 x - 1] == 3, x, Reals]",
             "Inequality[2, LessEqual, x, Less, Rational[5, 2]]");
    /* Mod -> Floor isolation (rational modulus). */
    run_test("Reduce[Mod[x, 4] == x, x, Reals]", "Inequality[0, LessEqual, x, Less, 4]");
    /* Mod with a transcendental modulus: transcendental breakpoints. */
    run_test("Reduce[Mod[x, 2 Pi] == x - 2 Pi, x, Reals]",
             "Inequality[Times[2, Pi], LessEqual, x, Less, Times[4, Pi]]");
    /* Inverse-trig identity, on its bounded domain. */
    run_test("Reduce[ArcSin[x] + ArcCos[x] == Pi/2, x, Reals]",
             "Inequality[-1, LessEqual, x, LessEqual, 1]");
    /* Log identity valid only on the real domain (x < 0). */
    run_test("Reduce[Log[x^2] == 2 Log[-x], x, Reals]", "Less[x, 0]");
    /* Abs inequality (previously declined). */
    run_test("Reduce[Abs[x] < 1, x, Reals]", "Inequality[-1, Less, x, Less, 1]");
    /* Simple radical equation / inequality (previously echoed unsolved). */
    run_test("Reduce[Sqrt[x - 1] == 2, x, Reals]", "Equal[x, 5]");
    run_test("Reduce[Sqrt[x - 1] < 5, x, Reals]", "Inequality[1, LessEqual, x, Less, 26]");

    /* Per-conjunct domain gate: a radical/Log domain from one Abs sign-branch must
     * not exclude the mutually-exclusive other branch.  Both were WRONG before the
     * fix (x==0, and False, respectively). */
    run_test("Reduce[Sqrt[Abs[x]] < 1, x, Reals]", "Inequality[-1, Less, x, Less, 1]");
    run_test("Reduce[Log[Abs[x]] < 0, x, Reals]",
             "Or[Inequality[-1, Less, x, Less, 0], Inequality[0, Less, x, Less, 1]]");

    /* Soundness: an out-of-domain point where the identity holds in C must be
     * EXCLUDED (ArcSin[2]+ArcCos[2]==Pi/2 is True in C but x=2 is not real). */
    run_test("Reduce[ArcSin[x] + ArcCos[x] == Pi/2 && x > 3/2, x, Reals]", "False");
    /* Soundness: a free parameter declines rather than guessing. */
    run_test("Reduce[Sqrt[x] == a, x, Reals]",
             "Reduce[Equal[Power[x, Rational[1, 2]], a], x, Reals]");
}

/* ------------------------------------------------------------------ *
 *  Piecewise functions (Sign/UnitStep/Ramp/Clip/Piecewise/Boole/     *
 *  HeavisideTheta/IntegerPart) case-split into polynomial branches;   *
 *  a polynomial-in-Floor with an unbounded integer set emits rays.    *
 * ------------------------------------------------------------------ */

static void test_piecewise_functions(void) {
    run_test("Reduce[Sign[x - 1] < 0, x, Reals]", "Less[x, 1]");
    run_test("Reduce[IntegerPart[x] == 2, x, Reals]",
             "Inequality[2, LessEqual, x, Less, 3]");
    run_test("Reduce[Clip[x, {-2, 2}] < 1, x, Reals]", "Less[x, 1]");
    run_test("Reduce[Boole[x > 0] + Boole[x > 1] == 2, x, Reals]", "Greater[x, 1]");
    /* HeavisideTheta is verified here (not in the sampling corpus): its numeric
     * value stays symbolic, so back-sampling has no judgeable grid points. */
    run_test("Reduce[HeavisideTheta[x - 2] == 1, x, Reals]", "Greater[x, 2]");
    /* Unbounded polynomial-in-Floor -> one-sided rays via the integer sub-solve. */
    run_test("Reduce[Floor[x]^2 > 5, x, Reals]",
             "Or[Less[x, -2], GreaterEqual[x, 3]]");
}

/* ------------------------------------------------------------------ *
 *  Phase 7: quantifier elimination (Exists / ForAll / Resolve)        *
 * ------------------------------------------------------------------ */

/* Case A -- fully-quantified sentences decide to True / False (a real-closed-
 * field decision procedure), reusing the whole engine over the bound vars. */
static void test_quantifiers_decision(void) {
    run_test("Resolve[Exists[x, x^2 == 4], Reals]", "True");
    run_test("Resolve[Exists[x, x^2 == -1], Reals]", "False");
    run_test("Resolve[ForAll[x, x^2 >= 0], Reals]", "True");
    run_test("Resolve[ForAll[x, x^2 > 0], Reals]", "False");
    /* Multivariate bound block, no free variables. */
    run_test("Resolve[Exists[{x, y}, x^2 + y^2 < 1], Reals]", "True");
    run_test("Resolve[ForAll[{x, y}, x^2 + y^2 >= 0], Reals]", "True");
    /* The body reduces to False -> Exists is False (short-circuit). */
    run_test("Resolve[Exists[x, x + 1 == x], Reals]", "False");
    /* A bound value that folds (HoldAll keeps x symbolic, not the constants). */
    run_test("Resolve[Exists[y, 1 + 1 < y], Reals]", "True");
}

/* Case B -- one free variable: parametric QE via CAD, emitted as a 1-D formula.
 * (The exact FullForm is pinned from the built binary and certified sound.) */
static void test_quantifiers_parametric(void) {
    run_test("Reduce[Exists[y, x^2 + y^2 < 1], {x}, Reals]",
             "Inequality[-1, Less, x, Less, 1]");
    run_test("Reduce[ForAll[y, x^2 + y^2 >= 1], {x}, Reals]",
             "Or[LessEqual[x, -1], GreaterEqual[x, 1]]");
    run_test("Resolve[Exists[x, x^2 == a], Reals]", "GreaterEqual[a, 0]");
    /* Cofinite ForAll: true for all x except the single excluded point. */
    run_test("Reduce[ForAll[y, x^2*(1 + y^2) > 0], {x}, Reals]", "Unequal[x, 0]");
    run_test("Reduce[ForAll[y, y^2 + x^2 > 0], {x}, Reals]", "Unequal[x, 0]");
    /* Free variable absent from the body -> a plain decision. */
    run_test("Reduce[Exists[y, y^2 == 2], {x}, Reals]", "True");
    /* Bound variable absent from the body -> the condition on the free one. */
    run_test("Reduce[Exists[y, x != 0], {x}, Reals]", "Unequal[x, 0]");
    /* nbound == 0: Exists[{}, g] == g -- the quantifier is stripped. */
    run_test("Reduce[Exists[{}, x > 0], {x}, Reals]", "Greater[x, 0]");
    /* A listed variable absent from the body is unconstrained (not a decline):
     * the complete description is the condition on the variable that appears. */
    run_test("Reduce[Exists[y, x^2 + y^2 < 1], {x, z}, Reals]",
             "Inequality[-1, Less, x, Less, 1]");
}

/* Soundness net: cases the v1 engine does not cover must stay UNEVALUATED. */
static void test_quantifiers_decline(void) {
    /* >=2 genuine free variables (algebraic boundary, needs Phase 6b). */
    run_contains("Resolve[Exists[x, x^2 + b*x + c == 0], Reals]",
                 "Resolve[Exists[x,");
    /* Alternating quantifier prefix. */
    run_contains("Resolve[ForAll[x, Exists[y, x + y == 0]], Reals]",
                 "Resolve[ForAll[x,");
    /* Explicit non-Reals domain with a quantifier. */
    run_contains("Reduce[Exists[y, x^2 + y^2 < 1], {x}, Complexes]",
                 "Reduce[Exists[y,");
    /* Resolve on a non-quantified statement is left alone in v1. */
    run_contains("Resolve[x^2 > 0, Reals]", "Resolve[");
    /* Bare Exists / ForAll are inert (stay symbolic). */
    run_contains("Exists[x, x > 0]", "Exists[x,");
    run_contains("ForAll[x, x^2 >= 0]", "ForAll[x,");
}

/* ------------------------------------------------------------------ *
 *  Phase 8: LogicalExpand (+ the NotElement head)                     *
 *                                                                     *
 * LogicalExpand distributes to disjunctive normal form over OPAQUE    *
 * Boolean atoms with idempotence/complementation/absorption           *
 * contractions, collapsing to True/False when the statement decides.  *
 * FullForm strings are pinned from the built binary; the True/False   *
 * collapses are certified by the fact that a DNF is unsatisfiable iff  *
 * it distributes to zero clauses.                                     *
 * ------------------------------------------------------------------ */
static void test_logical_expand(void) {
    /* De Morgan + distribution. */
    run_test("LogicalExpand[p && !(q || r)]", "And[p, Not[q], Not[r]]");
    run_test("LogicalExpand[!(a || b)]", "And[Not[a], Not[b]]");
    run_test("LogicalExpand[!(a && b)]", "Or[Not[a], Not[b]]");
    run_test("LogicalExpand[a && (b || c)]", "Or[And[a, b], And[a, c]]");
    run_test("LogicalExpand[(a || b) && !(c || d || e)]",
             "Or[And[a, Not[c], Not[d], Not[e]], And[b, Not[c], Not[d], Not[e]]]");

    /* Tautology / contradiction collapse. */
    run_test("LogicalExpand[a || !a]", "True");
    run_test("LogicalExpand[a && !a]", "False");
    run_test("LogicalExpand[(a || b) && !a && !b]", "False");
    run_test("LogicalExpand[Implies[a && b, a]]", "True");
    /* Multi-variable tautology proved by the Not[phi]-empty test. */
    run_test("LogicalExpand[Implies[(p || !r) && s, Implies[r && s, p && q || p]]]",
             "True");

    /* Xor expands to its odd-parity minterms. */
    run_test("LogicalExpand[Xor[p, q, r]]",
             "Or[And[p, q, r], And[p, Not[q], Not[r]], "
             "And[Not[p], q, Not[r]], And[Not[p], Not[q], r]]");

    /* Absorption over relational atoms (treated as opaque Booleans). */
    run_test("LogicalExpand[x == a && y == b || x == a || y == b]",
             "Or[Equal[x, a], Equal[y, b]]");

    /* Negation folds into the complementary relation head (no stray Not). */
    run_test("LogicalExpand[!(x == a)]", "Unequal[x, a]");
    run_test("LogicalExpand[!(x < 1)]", "GreaterEqual[x, 1]");

    /* Element over Alternatives distributes; its negation uses NotElement. */
    run_test("LogicalExpand[Element[x | y, Reals]]",
             "And[Element[x, Reals], Element[y, Reals]]");
    run_test("LogicalExpand[!Element[x | y | z, Reals]]",
             "Or[NotElement[x, Reals], NotElement[y, Reals], NotElement[z, Reals]]");

    /* A bare atom expands to itself (strips the wrapper, never declines). */
    run_test("LogicalExpand[x == a]", "Equal[x, a]");
    run_test("LogicalExpand[5]", "5");
    run_test("LogicalExpand[a && a]", "a");

    /* NotElement decides when membership decides, else stays symbolic. */
    run_test("NotElement[3, Reals]", "False");
    run_test("NotElement[I, Reals]", "True");
    run_contains("NotElement[x, Reals]", "NotElement[x, Reals]");
}

/* ------------------------------------------------------------------ *
 *  Options: Backsubstitution, Cubics, GeneratedParameters, Method,    *
 *  Modulus, Quartics, WorkingPrecision                                *
 * ------------------------------------------------------------------ */

/* Options[Reduce] reports the seven Mathematica-compatible options. */
static void test_options_registered(void) {
    run_test("Options[Reduce]",
        "List[Rule[Backsubstitution, False], Rule[Cubics, False], "
        "Rule[GeneratedParameters, C], Rule[Method, Automatic], "
        "Rule[Modulus, 0], Rule[Quartics, False], "
        "Rule[WorkingPrecision, Infinity]]");
}

/* Cubics: default emits Root[] for an irreducible cubic; -> True gives radicals. */
static void test_option_cubics(void) {
    run_contains("Reduce[x^3 + x + 1 == 0, x]", "Root");
    run_not_contains("Reduce[x^3 + x + 1 == 0, x, Cubics -> True]", "Root");
    run_contains("Reduce[x^3 + x + 1 == 0, x, Cubics -> True]", "Or[Equal[x");
}

/* Quartics: default emits Root[] for an irreducible quartic; -> True gives radicals. */
static void test_option_quartics(void) {
    run_contains("Reduce[x^4 - x - 1 == 0, x]", "Root");
    run_not_contains("Reduce[x^4 - x - 1 == 0, x, Quartics -> True]", "Root");
}

/* Modulus: residue enumeration over Z/pZ; 0 is characteristic 0; a
 * non-modular / symbolic modulus declines (unevaluated). */
static void test_option_modulus(void) {
    run_test("Reduce[x^2 == 2, x, Modulus -> 7]", "Or[Equal[x, 3], Equal[x, 4]]");
    run_test("Reduce[x^2 == 2, x, Modulus -> 5]", "False");
    run_test("Reduce[3 x == 1, x, Modulus -> 7]", "Equal[x, 5]");
    run_test("Reduce[x^2 + x + 1 == 0, x, Modulus -> 7]",
             "Or[Equal[x, 2], Equal[x, 4]]");
    run_test("Reduce[x^2 == 4, x, Modulus -> 0]", "Or[Equal[x, -2], Equal[x, 2]]");
    run_test("Reduce[Sin[x] == 0, x, Modulus -> 7]",
             "Reduce[Equal[Sin[x], 0], x, Rule[Modulus, 7]]");
    run_test("Reduce[x^2 == 2, x, Modulus -> p]",
             "Reduce[Equal[Power[x, 2], 2], x, Rule[Modulus, p]]");
}

/* GeneratedParameters: renames the generated free-parameter head (default C). */
static void test_option_generatedparameters(void) {
    run_test("Reduce[x + y == 3, {x, y}, Integers]",
             "And[Element[C[1], Integers], Equal[x, C[1]], "
             "Equal[y, Plus[3, Times[-1, C[1]]]]]");
    run_test("Reduce[x + y == 3, {x, y}, Integers, GeneratedParameters -> k]",
             "And[Element[k[1], Integers], Equal[x, k[1]], "
             "Equal[y, Plus[3, Times[-1, k[1]]]]]");
    run_test("Reduce[x + y == 3, {x, y}, Integers, GeneratedParameters -> m]",
             "And[Element[m[1], Integers], Equal[x, m[1]], "
             "Equal[y, Plus[3, Times[-1, m[1]]]]]");
}

/* Backsubstitution: accepted + honored.  The current linear engine emits the
 * fully-solved (grafted) form, which is the default (-> False) and also what
 * -> True requests, so both agree and neither echoes unevaluated. */
static void test_option_backsubstitution(void) {
    const char* solved =
        "And[Unequal[a, 1], Equal[x, Power[Plus[-1, a], -1]], "
        "Equal[y, Times[-1, Power[Plus[-1, a], -1]]]]";
    run_test("Reduce[a x + y == 1 && x + y == 0, {x, y}]", solved);
    run_test("Reduce[a x + y == 1 && x + y == 0, {x, y}, Backsubstitution -> True]", solved);
    run_test("Reduce[a x + y == 1 && x + y == 0, {x, y}, Backsubstitution -> False]", solved);
}

/* WorkingPrecision: does not change an exact symbolic answer -- Infinity
 * (default) and a finite value both give the exact result. */
static void test_option_workingprecision(void) {
    run_test("Reduce[x^2 > 1, x, WorkingPrecision -> 30]",
             "Or[Less[x, -1], Greater[x, 1]]");
    run_test("Reduce[x^2 > 1, x, WorkingPrecision -> Infinity]",
             "Or[Less[x, -1], Greater[x, 1]]");
}

/* Method: Automatic is accepted and inert (Automatic is the only method). */
static void test_option_method(void) {
    run_test("Reduce[x^2 == 4, x, Method -> Automatic]",
             "Or[Equal[x, -2], Equal[x, 2]]");
}

/* Peeling: options work with and without an explicit domain and in multiples;
 * an unknown trailing option leaves Reduce unevaluated (Reduce::optx). */
static void test_option_peeling(void) {
    run_test("Reduce[x^2 == 4, x, Cubics -> False]", "Or[Equal[x, -2], Equal[x, 2]]");
    run_test("Reduce[x^2 == 4, x, Reals, Quartics -> False]",
             "Or[Equal[x, -2], Equal[x, 2]]");
    run_test("Reduce[x^2 == 4, x, Method -> Automatic, Cubics -> False]",
             "Or[Equal[x, -2], Equal[x, 2]]");
    run_test("Reduce[x^2 == 4, x, Bogus -> 1]",
             "Reduce[Equal[Power[x, 2], 4], x, Rule[Bogus, 1]]");
}

/* ------------------------------------------------------------------ *
 *  FindInstance (Phase 8 companion)                                   *
 *                                                                     *
 *  Every returned instance is verified against the input by the       *
 *  builtin itself, so these pin the deterministic witness chosen for  *
 *  each shape.  {} is returned ONLY when the set is provably empty;    *
 *  a shape it can neither witness nor refute stays unevaluated.        *
 * ------------------------------------------------------------------ */
static void test_find_instance(void) {
    /* Complexes (default): univariate equation -> one root */
    run_test("FindInstance[x^2 == 2, x]",
             "List[List[Rule[x, Times[-1, Power[2, Rational[1, 2]]]]]]");
    run_test("FindInstance[x^3 - 2 x + 1 == 0, x]", "List[List[Rule[x, 1]]]");
    /* Complexes system Reduce declines on -> Solve fallback, free vars -> 0 */
    run_test("FindInstance[x^2 - y z == 1, {x, y, z}]",
             "List[List[Rule[x, -1], Rule[y, 0], Rule[z, 0]]]");

    /* Reals: inequality sign diagram, disk interior, linear system */
    run_test("FindInstance[x^5 - 2 x + 1 < 0, x, Reals]", "List[List[Rule[x, -3]]]");
    run_test("FindInstance[x^2 + y^2 <= 1, {x, y}, Reals]",
             "List[List[Rule[x, 0], Rule[y, 0]]]");
    run_test("FindInstance[2 x + 3 y - 5 z == 1 && 3 x - 4 y + 7 z == 3, {x, y, z}, Reals]",
             "List[List[Rule[x, 0], Rule[y, 22], Rule[z, 13]]]");
    run_test("FindInstance[x^2 >= 0, x, Reals]", "List[List[Rule[x, 0]]]");   /* tautology */

    /* Integers / Rationals */
    run_test("FindInstance[x^2 == 4, x, Integers]", "List[List[Rule[x, -2]]]");
    run_test("FindInstance[x^2 < 10 && x > 0, x, Integers, 3]",
             "List[List[Rule[x, 1]], List[Rule[x, 2]], List[Rule[x, 3]]]");
    run_test("FindInstance[x^2 - 3 y^2 == 1 && 10 < x < 100, {x, y}, Integers]",
             "List[List[Rule[x, 26], Rule[y, -15]]]");
    run_test("FindInstance[3 x == 2, x, Rationals]", "List[List[Rule[x, Rational[2, 3]]]]");

    /* Booleans: SAT via DNF; the assignment satisfies the formula */
    run_test("FindInstance[Xor[a, b, c, d] && (a || b) && ! (c || d), {a, b, c, d}, Booleans]",
             "List[List[Rule[a, True], Rule[b, False], Rule[c, False], Rule[d, False]]]");

    /* Modulus option routes over Z/pZ */
    run_test("FindInstance[x^3 - 2 x + 1 == 0, x, Modulus -> 5]", "List[List[Rule[x, 1]]]");

    /* {} ONLY when provably empty: unsatisfiable system, no rational sqrt, n=0 */
    run_test("FindInstance[x^2 + y^3 == 3 && x + 2 y >= 4 && x y == 5, {x, y}, Reals]", "List[]");
    run_test("FindInstance[x^2 == 2, x, Rationals]", "List[]");
    run_test("FindInstance[x^2 == 2, x, 0]", "List[]");

    /* Structured sampling: branch-cut disequations Reduce/Solve decline. */
    run_test("FindInstance[Sqrt[z^2] != z, z, Complexes]", "List[List[Rule[z, -1]]]");
    run_test("FindInstance[Log[x y] != Log[x] + Log[y], {x, y}, Complexes]",
             "List[List[Rule[x, -1], Rule[y, -1]]]");
    /* Rational equation with an excluded pole: Reduce yields the roots; the exact
     * zero-test verifier confirms frac==0 and root!=1 at the complex radical. */
    run_contains("FindInstance[(x^3 - 1)/(x - 1) == 0 && x != 1, x, Complexes]", "Complex[0, 1]");
    /* Sum-of-squares bounds the integer box; the bilinear equation is then checked. */
    run_contains("FindInstance[a b + b c + c d + d e == 0 && a^2 + b^2 + c^2 + d^2 + e^2 == 5, "
                 "{a, b, c, d, e}, Integers]", "Rule[a");
    /* 1-variable real transcendental root far out: Tan[x]==x above 10^6. */
    run_contains("FindInstance[Tan[x] == x && x > 10^6, x, Reals]", "1000001.9");
    /* Inexact-input feasibility (damped oscillation), two distinct points + a>0. */
    run_contains("FindInstance[E^(-a x) Cos[b x] == 0.1 && E^(-a y) Cos[b y] == 0.1 "
                 "&& x != y && a > 0, {a, b, x, y}, Reals]", "Rule[a");

    /* Structured sampling finds a witness Reduce/Solve decline (x != 0 over the
     * default Complexes: -1 is the first verified candidate). */
    run_test("FindInstance[x != 0, x]", "List[List[Rule[x, -1]]]");

    /* Sound declines (stay unevaluated -- never {} unless provably empty) */
    run_contains("FindInstance[x^2 + y z == 1 && x + 2 y <= 3 z + 1 && x y z > 7, {x, y, z}, Reals]",
                 "FindInstance");
    /* unknown option -> unevaluated */
    run_contains("FindInstance[x^2 == 2, x, Foo -> 1]", "FindInstance");
    /* bad variable spec -> unevaluated */
    run_contains("FindInstance[x^2 == 2, 3]", "FindInstance");

    /* Parametric Diophantine family: the generated parameter C[1] is instantiated
     * to the fundamental Pell solution (grid hits C[1] == 1). */
    run_test("FindInstance[x^2 - 61 y^2 == 1 && x > 0 && y > 0, {x, y}, Integers]",
             "List[List[Rule[x, 1766319049], Rule[y, 226153980]]]");

    /* Periodic transcendental instance over the Reals: Solve gives x == 1/(k Pi);
     * the parameter is solved against the window (k == 15916 -> 31831 = 2 k - 1). */
    run_test("FindInstance[Sin[1/x] == 0 && 0 < x < 10^-5, x, Reals]",
             "List[List[Rule[x, Times[Rational[1, 31831], Power[Pi, -1]]]]]");

    /* Bounded integer search: a reachable sum of cubes (3^3+4^3+5^3 == 6^3). */
    run_test("FindInstance[a^3 + b^3 + c^3 == d^3 && a > 0 && b > 0 && c > 0 && d > 0, {a, b, c, d}, Integers]",
             "List[List[Rule[a, 5], Rule[b, 4], Rule[c, 3], Rule[d, 6]]]");

    /* Indexed variables c[i] are accepted; finite box gives multiple witnesses. */
    run_test("FindInstance[c[1] + 2 c[2] == 3 && 0 <= c[1] <= 1 && 0 <= c[2] <= 1, {c[1], c[2]}, Integers]",
             "List[List[Rule[c[1], 1], Rule[c[2], 1]]]");
    run_test("FindInstance[c[1]^2 + c[2]^2 == 25 && c[1] > 0 && c[2] > 0, {c[1], c[2]}, Integers, 2]",
             "List[List[Rule[c[1], 3], Rule[c[2], 4]], List[Rule[c[1], 4], Rule[c[2], 3]]]");

    /* Finite-domain emptiness: 0/1 knapsack whose target (500) exceeds the sum of
     * all 15 primes (328) -> provably {} via the linear reach-range check. */
    run_test("FindInstance[Total[Array[c, 15]*Prime[Range[15]]] == 500 "
             "&& And @@ Thread[0 <= Array[c, 15] <= 1], Array[c, 15], Integers]", "List[]");

    /* Best-effort search declines out-of-reach systems (a^4+b^4+c^4==d^4, smallest
     * solution ~ 4*10^5) and transcendental emptiness (no real Exp==PolyGamma). */
    run_contains("FindInstance[a^4 + b^4 + c^4 == d^4 && a > 0 && b > 0 && c > 0 && d > 0, {a, b, c, d}, Integers]",
                 "FindInstance");
    run_contains("FindInstance[Exp[x] == PolyGamma[0, x] && x > 0, x, Reals]", "FindInstance");

    /* Booleans: Equivalent is now expanded by the DNF engine (was left opaque). */
    run_test("FindInstance[Xor[p, q] && Implies[q, r] && Not[Equivalent[p, r]], {p, q, r}, Booleans]",
             "List[List[Rule[p, True], Rule[q, False], Rule[r, False]]]");

    /* Numerical witness for a transcendental+inexact Real system Reduce declares
     * False unsoundly; assert the returned instance actually satisfies it. */
    run_test("(0 < x < 0.001 && Sin[1/x] > 0.999) /. "
             "First[FindInstance[0 < x < 0.001 && Sin[1/x] > 0.999, x, Reals]]", "True");

    /* Groebner emptiness: 2x2 nilpotent M^2==0 forces det==0, so det!=0 is empty. */
    run_test("FindInstance[a^2 + b c == 0 && a b + b d == 0 && a c + c d == 0 "
             "&& b c + d^2 == 0 && a d - b c != 0, {a, b, c, d}, Reals]", "List[]");

    /* Guards (must stay correct): Fermat n>=5 empty (FLT), and p*q with a prime
     * target empty -- both sound {} from Reduce, not the new mechanisms. */
    run_test("FindInstance[a^5 + b^5 == c^5 && a > 0 && b > 0 && c > 0, {a, b, c}, Integers]", "List[]");
    run_test("FindInstance[p q == 4172535013 && p > 1 && q > 1, {p, q}, Integers]", "List[]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_decides_true);
    TEST(test_decides_false);
    TEST(test_equations);
    TEST(test_equations_decline);
    TEST(test_real_inequalities);
    TEST(test_rational_inequalities);
    TEST(test_linear_systems);
    TEST(test_parametric_systems);
    TEST(test_parametric_systems_decline);
    TEST(test_integer_domain);
    TEST(test_rational_domain);
    TEST(test_inequality_defaults_reals);
    TEST(test_decline_soundness);
    TEST(test_real_functions);
    TEST(test_piecewise_functions);
    TEST(test_cad_real);
    TEST(test_cad_nvar);
    TEST(test_quantifiers_decision);
    TEST(test_quantifiers_parametric);
    TEST(test_quantifiers_decline);
    TEST(test_logical_expand);
    TEST(test_find_instance);
    TEST(test_wb_constant_atoms);
    TEST(test_wb_logic_fold);
    TEST(test_wb_atom_emit);
    TEST(test_wb_logic_expand);
    TEST(test_wb_unsupported);

    TEST(test_options_registered);
    TEST(test_option_cubics);
    TEST(test_option_quartics);
    TEST(test_option_modulus);
    TEST(test_option_generatedparameters);
    TEST(test_option_backsubstitution);
    TEST(test_option_workingprecision);
    TEST(test_option_method);
    TEST(test_option_peeling);

    printf("All reduce tests passed!\n");
    return 0;
}
