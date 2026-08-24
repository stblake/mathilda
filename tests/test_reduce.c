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
    /* !=  over Complexes has no ordering to reason about -> declines. */
    run_test("Reduce[x != 0, x]", "Reduce[Unequal[x, 0], x]");
    run_test("Reduce[x^2 != 4, x]", "Reduce[Unequal[Power[x, 2], 4], x]");
    /* Conjunction / disjunction of equations routes to the linear-only system
     * engine, which declines on the nonlinear polynomials. */
    run_test("Reduce[x^2 == 4 && x^3 == 8, x]",
             "Reduce[And[Equal[Power[x, 2], 4], Equal[Power[x, 3], 8]], x]");
    run_test("Reduce[x^2 == 4 || x == 5, x]",
             "Reduce[Or[Equal[Power[x, 2], 4], Equal[x, 5]], x]");
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

/* Phase-4 boundary: nonlinear systems are not handled (linear-only engine). */
static void test_parametric_systems_decline(void) {
    run_test("Reduce[x y == 1 && x + y == 3, {x, y}]",
             "Reduce[And[Equal[Times[x, y], 1], Equal[Plus[x, y], 3]], List[x, y]]");
    run_test("Reduce[x^2 + y^2 == 1 && x == 0, {x, y}]",
             "Reduce[And[Equal[Plus[Power[x, 2], Power[y, 2]], 1], Equal[x, 0]], "
             "List[x, y]]");
    /* A non-linear multivariate equation over Complexes needs CAD (phase 6). */
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
    /* Unbounded integer sets are not yet emitted as x>=1 style forms. */
    run_test("Reduce[x > 0, x, Integers]", "Reduce[Greater[x, 0], x, Integers]");
    run_test("Reduce[x < 5, x, Integers]", "Reduce[Less[x, 5], x, Integers]");
    run_test("Reduce[x^2 > 1, x, Integers]",
             "Reduce[Greater[Power[x, 2], 1], x, Integers]");
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

    /* Soundness: an out-of-domain point where the identity holds in C must be
     * EXCLUDED (ArcSin[2]+ArcCos[2]==Pi/2 is True in C but x=2 is not real). */
    run_test("Reduce[ArcSin[x] + ArcCos[x] == Pi/2 && x > 3/2, x, Reals]", "False");
    /* Soundness: a free parameter declines rather than guessing. */
    run_test("Reduce[Sqrt[x] == a, x, Reals]",
             "Reduce[Equal[Power[x, Rational[1, 2]], a], x, Reals]");
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
    TEST(test_cad_real);
    TEST(test_cad_nvar);
    TEST(test_wb_constant_atoms);
    TEST(test_wb_logic_fold);
    TEST(test_wb_atom_emit);
    TEST(test_wb_logic_expand);
    TEST(test_wb_unsupported);

    printf("All reduce tests passed!\n");
    return 0;
}
