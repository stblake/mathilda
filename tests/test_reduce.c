/*
 * test_reduce.c
 *
 * Unit tests for `Reduce` (src/solve/reduce.c) and its internal logical
 * normal-form layer (src/solve/reduce_form.{h,c}, reduce_atom.c).
 *
 * Phase 0 scope: the front-end (argument parsing, True/False short-circuit,
 * Reduce::ivar bad-variable handling) and the DNF layer (atom canonicalisation,
 * constant-atom decision, And/Or/Not construction, emission back to Expr).
 * The per-domain solving engines are added in later phases; here anything that
 * does not fully decide is expected to stay unevaluated.
 *
 * Black-box results are compared against FullForm strings.  White-box tests
 * drive the DNF layer directly on RAW (parse-only, un-evaluated) input, since
 * the evaluator pre-decides most constant logical combinations before Reduce
 * ever sees them.
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

/* Statements that fully decide collapse to True / False. */
static void test_decides_true(void) {
    run_test("Reduce[True, x]", "True");
    run_test("Reduce[1 < 2, x]", "True");
    run_test("Reduce[x == x, x]", "True");
    run_test("Reduce[2 <= 2, x]", "True");
    run_test("Reduce[Implies[True, 1 < 2], x]", "True");
    run_test("Reduce[! (3 < 2), x]", "True");
}

static void test_decides_false(void) {
    run_test("Reduce[False, x]", "False");
    run_test("Reduce[3 < 2, x]", "False");
    run_test("Reduce[2 > 5, x]", "False");
    run_test("Reduce[1 < 2 && 3 < 2, x]", "False");
}

/* Statements that need an engine not yet wired (inequalities, mixed systems,
 * multivariate) stay unevaluated. */
static void test_unevaluated(void) {
    run_test("Reduce[x > 0, x]", "Reduce[Greater[x, 0], x]");
    run_test("Reduce[x^2 == 4 && x > 0, x]",
             "Reduce[And[Equal[Power[x, 2], 4], Greater[x, 0]], x]");
    /* A non-linear multivariate system over Complexes is not yet handled (CAD). */
    run_test("Reduce[x^2 + y^2 == 1, {x, y}]",
             "Reduce[Equal[Plus[Power[x, 2], Power[y, 2]], 1], List[x, y]]");
}

/* Phase 1: complete univariate polynomial equation over Complexes. */
static void test_equations(void) {
    /* Parametric linear: the complete set keeps the degenerate a==0 branch. */
    run_test("Reduce[a x == b, x]",
             "Or[And[Unequal[a, 0], Equal[x, Times[Power[a, -1], b]]], "
             "And[Equal[a, 0], Equal[b, 0]]]");
    /* Numeric leading coefficient -> terminal, no case split. */
    run_test("Reduce[x^2 == 4, x]", "Or[Equal[x, -2], Equal[x, 2]]");
    run_test("Reduce[2 x == 6, x]", "Equal[x, 3]");
    run_test("Reduce[x - 5 == 0, x]", "Equal[x, 5]");
    run_test("Reduce[(x - 1)(x - 2) == 0, x]", "Or[Equal[x, 1], Equal[x, 2]]");
    /* Default domain is Complexes: complex roots are kept. */
    run_test("Reduce[x^2 == -1, x]",
             "Or[Equal[x, Complex[0, -1]], Equal[x, Complex[0, 1]]]");
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

/* Phase 2: univariate real sign diagram (equations + inequalities over Reals). */
static void test_real_inequalities(void) {
    run_test("Reduce[x^2 > 1, x, Reals]",  "Or[Less[x, -1], Greater[x, 1]]");
    run_test("Reduce[x^2 >= 1, x, Reals]", "Or[LessEqual[x, -1], GreaterEqual[x, 1]]");
    run_test("Reduce[x^2 < 1, x, Reals]",  "Inequality[-1, Less, x, Less, 1]");
    run_test("Reduce[(x - 1)(x - 2)(x - 3) > 0, x, Reals]",
             "Or[Inequality[1, Less, x, Less, 2], Greater[x, 3]]");
    /* Cofinite: complement of finitely many points. */
    run_test("Reduce[x^2 != 1, x, Reals]", "And[Unequal[x, -1], Unequal[x, 1]]");
    /* Equations over Reals fall out of the same sign diagram. */
    run_test("Reduce[x^2 == 4, x, Reals]", "Or[Equal[x, -2], Equal[x, 2]]");
    /* No real breakpoints -> the whole line decides. */
    run_test("Reduce[x^2 + 1 > 0, x, Reals]", "True");
    run_test("Reduce[x^2 + 1 < 0, x, Reals]", "False");
    /* A conjunction of bounds. */
    run_test("Reduce[x > 0 && x < 1, x, Reals]", "Inequality[0, Less, x, Less, 1]");
    /* Algebraic (radical) breakpoints, ordered and signed via the qqbar oracle. */
    run_test("Reduce[x^2 < 2, x, Reals]",
             "Inequality[Times[-1, Power[2, Rational[1, 2]]], Less, x, Less, "
             "Power[2, Rational[1, 2]]]");
}

/* Phase 3: multivariate linear systems over Reals (Fourier-Motzkin). */
static void test_linear_systems(void) {
    /* The plan's flagship: a triangular description of the feasible region. */
    run_test("Reduce[x + y < 1 && x > 0 && y > 0, {x, y}, Reals]",
             "And[Inequality[0, Less, x, Less, 1], "
             "Inequality[0, Less, y, Less, Plus[1, Times[-1, x]]]]");
    /* Infeasible system -> False. */
    run_test("Reduce[x > 1 && x < 0 && y > 0, {x, y}, Reals]", "False");
    /* An equation is re-detected as `==` in the triangular output. */
    run_test("Reduce[x + y == 1 && x > 0, {x, y}, Reals]",
             "And[Greater[x, 0], Equal[y, Plus[1, Times[-1, x]]]]");
    /* Rational bound coefficients. */
    run_test("Reduce[2 x + 3 y <= 6 && x >= 0 && y >= 0, {x, y}, Reals]",
             "And[Inequality[0, LessEqual, x, LessEqual, 3], "
             "Inequality[0, LessEqual, y, LessEqual, Plus[2, Times[Rational[-2, 3], x]]]]");
    /* A free variable is simply omitted. */
    run_test("Reduce[x > 0, {x, y}, Reals]", "Greater[x, 0]");
    /* Disjunction: each conjunct solved and OR-ed. */
    run_test("Reduce[x < 0 || x > 1, {x, y}, Reals]",
             "Or[Less[x, 0], Greater[x, 1]]");
}

/* Phase 5: Integers / Rationals domain. */
static void test_integer_domain(void) {
    /* Equation -> finite solution set. */
    run_test("Reduce[x^2 == 4, x, Integers]", "Or[Equal[x, -2], Equal[x, 2]]");
    run_test("Reduce[x^2 == 2, x, Integers]", "False");
    run_test("Reduce[x^2 == 4, x, Rationals]", "Or[Equal[x, -2], Equal[x, 2]]");
    /* Bounded inequality -> integer enumeration (Solve declines, we fall back). */
    run_test("Reduce[x^2 < 10 && x > 0, x, Integers]",
             "Or[Equal[x, 1], Equal[x, 2], Equal[x, 3]]");
    run_test("Reduce[1 <= x <= 3, x, Integers]",
             "Or[Equal[x, 1], Equal[x, 2], Equal[x, 3]]");
    /* Bounded system with an equation. */
    run_test("Reduce[x + y == 5 && x > 0 && y > 0, {x, y}, Integers]",
             "Or[And[Equal[x, 1], Equal[y, 4]], And[Equal[x, 2], Equal[y, 3]], "
             "And[Equal[x, 3], Equal[y, 2]], And[Equal[x, 4], Equal[y, 1]]]");
    /* Parametric linear Diophantine -> a C[k] family with an Element condition. */
    run_test("Reduce[2 x + 3 y == 1, {x, y}, Integers]",
             "And[Element[C[1], Integers], Equal[x, Plus[-1, Times[3, C[1]]]], "
             "Equal[y, Plus[1, Times[-2, C[1]]]]]");
}

/* Phase 4: parametric linear systems over Complexes (case analysis). */
static void test_parametric_systems(void) {
    /* Numeric determined system. */
    run_test("Reduce[x + y == 3 && x - y == 1, {x, y}]",
             "And[Equal[y, 1], Equal[x, 2]]");
    /* Underdetermined: a variable stays free. */
    run_test("Reduce[x + y == 1, {x, y}]", "Equal[x, Plus[1, Times[-1, y]]]");
    /* Parametric square system: genericity condition + Cramer solution. */
    run_test("Reduce[a x + y == 1 && x + y == 0, {x, y}]",
             "And[Unequal[Plus[1, Times[-1, a]], 0], Equal[x, Power[Plus[-1, a], -1]], "
             "Equal[y, Times[-1, Power[Plus[-1, a], -1]]]]");
    /* Overdetermined: a consistency condition on the parameter. */
    run_test("Reduce[a x == 1 && x == 2, {x}]",
             "And[Equal[Plus[-1, Times[2, a]], 0], Equal[x, 2]]");
    /* Three variables. */
    run_test("Reduce[x + y + z == 6 && x - y == 0 && z == 2, {x, y, z}]",
             "And[Equal[z, 2], Equal[y, 2], Equal[x, 2]]");
    /* A non-linear system is declined (stays unevaluated). */
    run_test("Reduce[x y == 1 && x + y == 3, {x, y}]",
             "Reduce[And[Equal[Times[x, y], 1], Equal[Plus[x, y], 3]], List[x, y]]");
}

/* Invalid variable spec -> Reduce::ivar (stderr) + unevaluated. */
static void test_bad_vars(void) {
    run_test("Reduce[x == 1, 5]", "Reduce[Equal[x, 1], 5]");
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
    wb_test("3 < 2 || 5 > 2", "True");        /* False ∨ True  = True  */
    wb_test("1 < 2 && 3 < 2", "False");       /* True  ∧ False = False */
    /* A true constant conjunct drops, leaving the symbolic atom. */
    wb_test("5 > 2 && x^2 == 4", "Equal[Plus[-4, Power[x, 2]], 0]");
    /* A false constant conjunct kills the whole conjunction. */
    wb_test("3 < 2 && x^2 == 4", "False");
}

/* A single symbolic atom canonicalises to `poly REL 0`. */
static void test_wb_atom_emit(void) {
    wb_test("x^2 == 4", "Equal[Plus[-4, Power[x, 2]], 0]");
    wb_test("x != 1", "Unequal[Plus[-1, x], 0]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_decides_true);
    TEST(test_decides_false);
    TEST(test_unevaluated);
    TEST(test_equations);
    TEST(test_real_inequalities);
    TEST(test_linear_systems);
    TEST(test_parametric_systems);
    TEST(test_integer_domain);
    TEST(test_bad_vars);
    TEST(test_wb_constant_atoms);
    TEST(test_wb_logic_fold);
    TEST(test_wb_atom_emit);

    printf("All reduce tests passed!\n");
    return 0;
}
