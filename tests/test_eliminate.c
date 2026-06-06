/* Unit tests for Eliminate (src/poly/eliminate.c).
 *
 * Eliminate is a thin driver over the lex-order Buchberger engine in
 * src/poly/groebner.c with an elimination block, plus a principal-branch
 * inverse-function pre-pass and a balanced `Equal[lhs, rhs]` reconstruction.
 *
 * The expected outputs are taken from a smoke run of the implementation
 * (cross-checked for mathematical equivalence with Mathematica's reference
 * answers for the documented examples).  Matching FullForm byte-for-byte
 * is the strongest correctness signal we have.
 *
 * Run binary directly: ./eliminate_tests
 * (per MEMORY.md note: ctest is not configured in tests/CMakeLists.txt). */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

/* FullForm string match. */
static void check_eq(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, expected) != 0) {
        fprintf(stdout, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* "expr evaluates to True" assertion (avoids hard-coded long FullForm). */
static void check_true(const char* input) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, "True") != 0) {
        fprintf(stdout, "FAIL: %s\n  expected: True\n  got:      %s\n",
                input, got);
        ASSERT_STR_EQ(got, "True");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Silence stderr once for diagnostic-emitting tests. */
static void mute_stderr_once(void) {
    static int done = 0;
    if (!done) { freopen("/dev/null", "w", stderr); done = 1; }
}

/* ------------------------------------------------------------------ */
/* 1. Linear pair: headline example from the prompt                    */
/* ------------------------------------------------------------------ */

static void test_linear_pair(void) {
    /* Eliminate[{x == 2+y, y == z}, y] -> 2 + z == x */
    check_eq("Eliminate[{x == 2+y, y == z}, y]",
             "Equal[Plus[2, z], x]");
}

/* ------------------------------------------------------------------ */
/* 2. Linear system, eliminate z                                       */
/* ------------------------------------------------------------------ */

static void test_linear_system_eliminate_one(void) {
    /* Eliminate[{2x+3y+4z == 1, 9x+8y+7z == 2}, z]
     * Mathematica: 1 - 11 y == 22 x.  Sign-equivalent FullForm below. */
    check_eq("Eliminate[{2x+3y+4z == 1, 9x+8y+7z == 2}, z]",
             "Equal[Plus[Times[22, x], Times[11, y]], 1]");
}

/* ------------------------------------------------------------------ */
/* 3. Symmetric-function elimination                                   */
/* ------------------------------------------------------------------ */

static void test_symmetric_functions(void) {
    /* Eliminate[{f == x^5+y^5, a == x+y, b == x*y}, {x, y}]
     * Mathematica: f == a^5 - 5 a^3 b + 5 a b^2 (sign-balanced equiv). */
    check_eq("Eliminate[{f == x^5+y^5, a == x+y, b == x*y}, {x, y}]",
             "Equal[Plus[Power[a, 5], Times[5, Times[a, Power[b, 2]]]], "
                   "Plus[Times[5, Times[Power[a, 3], b]], f]]");
}

/* ------------------------------------------------------------------ */
/* 4. Polynomial system, eliminate z -> two equations remain           */
/* ------------------------------------------------------------------ */

static void test_poly_system_eliminate_z(void) {
    /* Eliminate[{x^2+y^2+z^2 == 1, x-y+z == 2, x^3-y^2 == z+1}, z]
     * One purely-polynomial constraint on x plus a y-eq parametric in x. */
    check_eq("Eliminate[{x^2+y^2+z^2 == 1, x-y+z == 2, x^3-y^2 == z+1}, z]",
             "And[Equal[Plus[27, Times[4, Power[x, 2]], Times[8, Power[x, 4]], "
                            "Times[4, Power[x, 5]], Times[4, Power[x, 6]]], "
                       "Plus[Times[18, x], Times[28, Power[x, 3]]]], "
                  "Equal[Plus[12, Times[2, x], Times[5, Power[x, 2]], y], "
                       "Plus[Times[8, Power[x, 3]], Times[4, Power[x, 4]], "
                            "Times[2, Power[x, 5]]]]]");
}

/* ------------------------------------------------------------------ */
/* 5. Same system, two-variable elimination -> single equation         */
/* ------------------------------------------------------------------ */

static void test_poly_system_eliminate_y_z(void) {
    check_eq("Eliminate[{x^2+y^2+z^2 == 1, x-y+z == 2, x^3-y^2 == z+1}, {y, z}]",
             "Equal[Plus[27, Times[4, Power[x, 2]], Times[8, Power[x, 4]], "
                        "Times[4, Power[x, 5]], Times[4, Power[x, 6]]], "
                   "Plus[Times[18, x], Times[28, Power[x, 3]]]]");
}

/* ------------------------------------------------------------------ */
/* 6. && form input is equivalent to {} form                           */
/* ------------------------------------------------------------------ */

static void test_and_form_equivalent(void) {
    check_true(
        "Eliminate[(2x+3y+4z == 1) && (9x+8y+7z == 2), z] === "
        "Eliminate[{2x+3y+4z == 1, 9x+8y+7z == 2}, z]");
}

/* ------------------------------------------------------------------ */
/* 7. Single-symbol var shorthand: y vs {y}                            */
/* ------------------------------------------------------------------ */

static void test_single_symbol_var_shorthand(void) {
    check_true(
        "Eliminate[{x == 2+y, y == z}, y] === "
        "Eliminate[{x == 2+y, y == z}, {y}]");
}

/* ------------------------------------------------------------------ */
/* 8. Resultant cross-check: incompatible roots -> False               */
/* ------------------------------------------------------------------ */

static void test_resultant_incompatible(void) {
    /* x^2 - 1 = 0 and x^2 - 4 = 0 share no root. */
    check_eq("Eliminate[{x^2 - 1 == 0, x^2 - 4 == 0}, x]",
             "False");
}

/* ------------------------------------------------------------------ */
/* 9. Common-root condition: a == b                                    */
/* ------------------------------------------------------------------ */

static void test_common_root_condition(void) {
    /* x - a = 0 and x - b = 0 share a root iff a == b. */
    check_eq("Eliminate[{x - a == 0, x - b == 0}, x]",
             "Equal[b, a]");
}

/* ------------------------------------------------------------------ */
/* 10. Empty elimination ideal -> True                                 */
/* ------------------------------------------------------------------ */

static void test_empty_ideal_returns_true(void) {
    /* One equation, one unknown to eliminate -> y solvable for any x. */
    check_eq("Eliminate[x + y == 0, y]", "True");
}

/* ------------------------------------------------------------------ */
/* 11. Numeric contradiction -> False                                  */
/* ------------------------------------------------------------------ */

static void test_contradiction_returns_false(void) {
    /* `1 == 2` evaluates to `False` before Eliminate runs; we fold that
     * sentinel rather than emit Eliminate::eqf. */
    check_eq("Eliminate[1 == 2, x]", "False");
}

/* ------------------------------------------------------------------ */
/* 12. Tautology -> True                                               */
/* ------------------------------------------------------------------ */

static void test_tautology_returns_true(void) {
    /* `x == x` evaluates to `True`; same sentinel-folding path. */
    check_eq("Eliminate[x == x, x]", "True");
}

/* ------------------------------------------------------------------ */
/* 13. Transcendental pre-pass: principal-branch inverse rewrite       */
/* ------------------------------------------------------------------ */

static void test_inverse_function_prepass(void) {
    /* `Eliminate[Sin[y] == 1, y]` peels Sin via ArcSin (ifun fires),
     * leaves y == Pi/2; eliminating y from a single solvable equation
     * yields no residual constraint -> True. */
    mute_stderr_once();
    check_eq("Eliminate[Sin[y] == 1, y]", "True");
}

static void test_inverse_function_two_vars(void) {
    /* Sin[x + y] == 1 -> x + y == Pi/2, then eliminating y is trivial. */
    mute_stderr_once();
    check_eq("Eliminate[Sin[x + y] == 1, y]", "True");
}

/* ------------------------------------------------------------------ */
/* 14. Empty variable list -> equations pass through                   */
/* ------------------------------------------------------------------ */

static void test_empty_var_list(void) {
    /* No variables to eliminate -> just return the equation unchanged. */
    check_eq("Eliminate[{x == y + 1}, {}]", "Equal[x, Plus[1, y]]");
}

/* ------------------------------------------------------------------ */
/* 15. Attribute check: Eliminate is Protected                         */
/* ------------------------------------------------------------------ */

static void test_eliminate_protected(void) {
    /* Setting Eliminate[x_] := ... must fail on a Protected head. */
    mute_stderr_once();
    check_true("MemberQ[Attributes[Eliminate], Protected]");
}

/* ------------------------------------------------------------------ */
/* 16. Memory smoke -- run the symmetric-function spec example         */
/* ------------------------------------------------------------------ */

static void test_memory_smoke(void) {
    for (int i = 0; i < 25; i++) {
        Expr* e = parse_expression(
            "Eliminate[{f == x^5+y^5, a == x+y, b == x*y}, {x, y}]");
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(e);
        expr_free(r);
    }
}

/* ------------------------------------------------------------------ */
/* 17. Memory smoke -- larger polynomial system                        */
/* ------------------------------------------------------------------ */

static void test_memory_smoke_larger(void) {
    for (int i = 0; i < 25; i++) {
        Expr* e = parse_expression(
            "Eliminate[{x^2+y^2+z^2 == 1, x-y+z == 2, x^3-y^2 == z+1}, "
            "{y, z}]");
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(e);
        expr_free(r);
    }
}

/* ------------------------------------------------------------------ */
/* 18-23. Algebraisation pre-pass (radical / rational-power inputs)    */
/* ------------------------------------------------------------------ */

/* Each Sqrt-bearing or x^(p/q)-bearing equation is rewritten as a
 * polynomial in a fresh aux symbol, plus the constraint `aux^L == base`.
 * The expected output is what Mathematica returns (the
 * cross-multiplied generic consequence; sign / branch info is dropped
 * along with the `Eliminate::alg` diagnostic). */

static void test_radical_square(void) {
    /* Eliminate[{Sqrt[x] == y, x == z^2}, x]  ->  y^2 == z^2 */
    mute_stderr_once();
    check_eq("Eliminate[{Sqrt[x] == y, x == z^2}, x]",
             "Equal[Power[z, 2], Power[y, 2]]");
}

static void test_radical_inv_pair(void) {
    /* u == Sqrt[x^2+1] and v == 1/Sqrt[x^2+1]  ->  u*v == 1 (since
     * Sqrt[..] * 1/Sqrt[..] = 1 generically).  Both rationals share
     * the same base, so a single aux variable is introduced. */
    mute_stderr_once();
    check_eq("Eliminate[{u == Sqrt[x^2+1], v == 1/Sqrt[x^2+1]}, x]",
             "Equal[Times[u, v], 1]");
}

static void test_radical_cuberoot_pair(void) {
    /* u == x^(1/3), v == x^(2/3)  ->  v == u^2 (Power[x, 1/3] and
     * Power[x, 2/3] collapse onto Power[aux, 1] and Power[aux, 2]). */
    mute_stderr_once();
    check_eq("Eliminate[{u == x^(1/3), v == x^(2/3)}, x]",
             "Equal[v, Power[u, 2]]");
}

static void test_radical_dt_y_arc(void) {
    /* Headline example from the user prompt:
     *   Dt[y] == x^3 / Sqrt[x^2+1] Dt[x]
     *   u     == x^2 + 1
     *   Dt[u] == 2 x Dt[x]
     * Eliminating {Dt[x], x} yields a polynomial relation among
     * u, Dt[u], Dt[y]; the answer matches Mathematica's modulo a sign
     * rearrangement.  Mathematica reports:
     *   u^2 Dt[u]^2 + u (-2 Dt[u]^2 - 4 Dt[y]^2) == -Dt[u]^2
     * which equals the form below after moving negative terms across. */
    mute_stderr_once();
    check_eq("Eliminate[{Dt[y] == x^3/Sqrt[x^2+1] Dt[x], u == x^2 + 1,"
             " Dt[u] == 2 x Dt[x]}, {Dt[x], x}]",
             "Equal[Plus[Power[Dt[u], 2], Times[Power[u, 2], Power[Dt[u], 2]]], "
                   "Plus[Times[2, Times[u, Power[Dt[u], 2]]], "
                        "Times[4, Times[u, Power[Dt[y], 2]]]]]");
}

static void test_radical_nested_sqrt(void) {
    /* Nested radical: Sqrt[x + Sqrt[x]] in one equation, 1/Sqrt[x] in
     * another.  Two aux variables are needed (one for Sqrt[x], one for
     * Sqrt[x + aux_inner]).  Mathematica returns:
     *   4 u^4 Dt[u]^2 - 2 u Dt[u] Dt[y] - 8 u^3 Dt[u] Dt[y]
     *     + 4 u^2 Dt[y]^2 == -Dt[y]^2 */
    mute_stderr_once();
    check_eq("Eliminate[{Dt[y] == Sqrt[x + Sqrt[x]] Dt[x],"
             " u == Sqrt[x + Sqrt[x]],"
             " Dt[u] == (1 + 1/2/Sqrt[x]) Dt[x]}, {Dt[x], x}]",
             "Equal[Plus[Times[4, Times[Power[u, 4], Power[Dt[u], 2]]], "
                       "Power[Dt[y], 2], "
                       "Times[4, Times[Power[u, 2], Power[Dt[y], 2]]]], "
                   "Plus[Times[2, Times[u, Dt[u], Dt[y]]], "
                        "Times[8, Times[Power[u, 3], Dt[u], Dt[y]]]]]");
}

static void test_radical_quartic_with_d(void) {
    /* `D[x + 1/x, x]` is evaluated to `1 - 1/x^2` before Eliminate
     * runs; the Sqrt[x^4+1] is then the sole algebraic atom.  The
     * Numerator[Together[...]] normalisation also has to clear the
     * `1/x` denominator inside `u == x + 1/x`.  Mathematica returns:
     *   -2 u^2 Dt[y]^2 + u^4 Dt[y]^2 == Dt[u]^2 */
    mute_stderr_once();
    check_eq("Eliminate[{Dt[y] == (x^2 - 1)/((x^2 + 1) Sqrt[x^4 + 1]) Dt[x],"
             " u == x + 1/x,"
             " Dt[u] == (D[x + 1/x, x]) Dt[x]}, {Dt[x], x}]",
             "Equal[Plus[Power[Dt[u], 2], "
                       "Times[2, Times[Power[u, 2], Power[Dt[y], 2]]]], "
                   "Times[Power[u, 4], Power[Dt[y], 2]]]");
}

static void test_radical_sqrtxp1_with_d(void) {
    /* Single Sqrt[x+1] appearing in two equations, plus a 1/(x^2+1)
     * rational coefficient.  Confirms the monomial-factor strip
     * runs (without it the answer carries an extra factor of u, since
     * `u == Sqrt[x+1]` introduces a u-multiplier that survives
     * Buchberger).  Mathematica returns:
     *   u^4 (2 Dt[u] - Dt[y]) + u^2 (-4 Dt[u] + 2 Dt[y]) == 2 Dt[y] */
    mute_stderr_once();
    check_eq("Eliminate[{Dt[y] == (x^2 - 1)/((x^2 + 1) Sqrt[x + 1]) Dt[x],"
             " u == Sqrt[x + 1],"
             " Dt[u] == (D[Sqrt[x + 1], x]) Dt[x]}, {Dt[x], x}]",
             "Equal[Plus[Times[2, Times[Power[u, 4], Dt[u]]], "
                       "Times[2, Times[Power[u, 2], Dt[y]]]], "
                   "Plus[Times[4, Times[Power[u, 2], Dt[u]]], "
                        "Times[2, Dt[y]], "
                        "Times[Power[u, 4], Dt[y]]]]");
}

static void test_radical_memory_smoke(void) {
    /* Re-run the headline radical example repeatedly to surface any
     * AlgState / aux-symbol cleanup leaks. */
    mute_stderr_once();
    for (int i = 0; i < 25; i++) {
        Expr* e = parse_expression(
            "Eliminate[{Dt[y] == x^3/Sqrt[x^2+1] Dt[x],"
            " u == x^2 + 1, Dt[u] == 2 x Dt[x]}, {Dt[x], x}]");
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(e);
        expr_free(r);
    }
}

/* ------------------------------------------------------------------ */
/* 24. Truly non-polynomial: Sin[x*y] still triggers nlin             */
/* ------------------------------------------------------------------ */

static void test_nlin_for_transcendental_arg(void) {
    /* `Sin[x*y]` has an elim variable inside a non-polynomial head;
     * the algebraisation pass cannot help, so collect_main_vars'
     * post-validation must reject the input with `Eliminate::nlin`
     * and return Eliminate[...] unevaluated. */
    mute_stderr_once();
    check_eq("Eliminate[{Sin[x*y] + Cos[x] == 1, x + y == 2}, x]",
             "Eliminate[List[Equal[Plus[Cos[x], Sin[Times[x, y]]], 1], "
                            "Equal[Plus[x, y], 2]], x]");
}

/* ------------------------------------------------------------------ */
/* 25-32. Trigonometric algebraisation pre-pass                        */
/* ------------------------------------------------------------------ */

/* When an elim variable sits inside circular/hyperbolic trig functions,
 * each `Sin[theta]`/`Cos[theta]` (etc.) is replaced by fresh aux symbols
 * with the Pythagorean constraint, after a TrigExpand that reduces
 * multiple/sum angles to atomic Sin/Cos.  Expected outputs are the
 * cross-multiplied generic consequence (sign/branch info dropped along
 * with the `Eliminate::ifun` diagnostic), cross-checked for equivalence
 * with Mathematica. */

static void test_trig_integration_substitution(void) {
    /* Headline example: substitution method for
     *   Integrate[Cos[x] Sqrt[1 - Sin[x]], x]  with  u = Sin[x].
     * Mathematica returns  -Dt[y]^2 == (-1 + u) Dt[u]^2, i.e.
     *   Dt[y]^2 == (1 - u) Dt[u]^2, the balanced form below. */
    mute_stderr_once();
    check_eq("Eliminate[{Dt[y] == Cos[x] Sqrt[1 - Sin[x]] Dt[x], u == Sin[x],"
             " Dt[u == Sin[x]]}, {x, Dt[x]}]",
             "Equal[Plus[Times[u, Power[Dt[u], 2]], Power[Dt[y], 2]], "
                   "Power[Dt[u], 2]]");
}

static void test_trig_pythagorean(void) {
    /* u == Sin[x], v == Cos[x]  ->  u^2 + v^2 == 1. */
    mute_stderr_once();
    check_eq("Eliminate[{u == Sin[x], v == Cos[x]}, x]",
             "Equal[Plus[Power[u, 2], Power[v, 2]], 1]");
}

static void test_trig_sin_squared(void) {
    /* Two occurrences of Sin[x] (one squared) collapse onto one aux. */
    mute_stderr_once();
    check_eq("Eliminate[{u == Sin[x], w == Sin[x]^2}, x]",
             "Equal[w, Power[u, 2]]");
}

static void test_trig_double_angle(void) {
    /* Sin[2x] is TrigExpand'd to 2 Sin[x] Cos[x] onto the SAME atomic
     * angle x, so the relation v^2 == 4 u^2 (1 - u^2) survives. */
    mute_stderr_once();
    check_eq("Eliminate[{u == Sin[x], v == Sin[2 x]}, x]",
             "Equal[Plus[Times[4, Power[u, 4]], Power[v, 2]], "
                   "Times[4, Power[u, 2]]]");
}

static void test_trig_tangent(void) {
    /* Tan[x] is rewritten through Sin/Cos; the cleared denominator gives
     *   v^2 + u^2 v^2 == 1   (i.e. Tan^2 == Sec^2 - 1). */
    mute_stderr_once();
    check_eq("Eliminate[{u == Tan[x], v == Cos[x]}, x]",
             "Equal[Plus[Power[v, 2], Times[Power[u, 2], Power[v, 2]]], 1]");
}

static void test_trig_hyperbolic(void) {
    /* a == Cosh[x], b == Sinh[x]  ->  a^2 - b^2 == 1 (hyperbolic
     * Pythagorean: Cosh^2 - Sinh^2 == 1). */
    mute_stderr_once();
    check_eq("Eliminate[{a == Cosh[x], b == Sinh[x]}, x]",
             "Equal[Plus[1, Power[b, 2]], Power[a, 2]]");
}

static void test_trig_free_var_returns_true(void) {
    /* a == Cos[x] Dt[x], b == Sin[x], eliminate {x, Dt[x]}: Dt[x] is a
     * free differential, so a is unconstrained -> no relation -> True.
     * Exercises Gate B's poly-atom test (x inside the elim atom Dt[x]
     * must NOT count as a free occurrence). */
    mute_stderr_once();
    check_true("Eliminate[{a == Cos[x] Dt[x], b == Sin[x]}, {x, Dt[x]}]");
}

static void test_trig_mixed_arg_nlin(void) {
    /* Bare `x` alongside `Sin[x]` severs the x <-> Sin[x] link; Gate B
     * rejects it as `nlin` and returns the input unevaluated rather than
     * emit an unsound relation. */
    mute_stderr_once();
    check_eq("Eliminate[{u == Sin[x], x + y == 2}, x]",
             "Eliminate[List[Equal[u, Sin[x]], Equal[Plus[x, y], 2]], x]");
}

/* ------------------------------------------------------------------ */
/* 33-37. Exponential / logarithmic algebraisation                     */
/* ------------------------------------------------------------------ */

/* Exp (`b^x`) and Log kernels are handled by the same pass with a single
 * algebraically-free aux per kernel; multiple/sum exponents and product
 * logs are reduced onto a common atomic kernel before substitution. */

static void test_exp_double_angle(void) {
    /* Exp[2x] = Exp[x]^2 collapses onto one aux -> v == u^2. */
    mute_stderr_once();
    check_eq("Eliminate[{u == Exp[x], v == Exp[2 x]}, x]",
             "Equal[v, Power[u, 2]]");
}

static void test_exp_general_base(void) {
    /* Non-E base: 2^(3x) = (2^x)^3 -> v == u^3. */
    mute_stderr_once();
    check_eq("Eliminate[{u == 2^x, v == 2^(3 x)}, x]",
             "Equal[v, Power[u, 3]]");
}

static void test_exp_integration_substitution(void) {
    /* Substitution method for Integrate[Exp[x]/(1 + Exp[x]), x] with
     * u = Exp[x]:  Dt[y] == Dt[u]/(1 + u), the balanced form below. */
    mute_stderr_once();
    check_eq("Eliminate[{Dt[y] == Exp[x]/(1 + Exp[x]) Dt[x], u == Exp[x],"
             " Dt[u] == Exp[x] Dt[x]}, {x, Dt[x]}]",
             "Equal[Dt[u], Plus[Dt[y], Times[u, Dt[y]]]]");
}

static void test_log_power(void) {
    /* Log[x^3] = 3 Log[x] -> v == 3 u. */
    mute_stderr_once();
    check_eq("Eliminate[{u == Log[x], v == Log[x^3]}, x]",
             "Equal[v, Times[3, u]]");
}

static void test_log_squared(void) {
    /* Two occurrences of Log[x] (one squared) collapse onto one aux. */
    mute_stderr_once();
    check_eq("Eliminate[{u == Log[x], w == Log[x]^2}, x]",
             "Equal[w, Power[u, 2]]");
}

static void test_trig_memory_smoke(void) {
    /* Re-run the headline trig example to surface TrigState / aux-symbol
     * cleanup leaks across the multi-exit cleanup paths. */
    mute_stderr_once();
    for (int i = 0; i < 25; i++) {
        Expr* e = parse_expression(
            "Eliminate[{Dt[y] == Cos[x] Sqrt[1 - Sin[x]] Dt[x], u == Sin[x],"
            " Dt[u] == Cos[x] Dt[x]}, {x, Dt[x]}]");
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(e);
        expr_free(r);
    }
}

/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_linear_pair);
    TEST(test_linear_system_eliminate_one);
    TEST(test_symmetric_functions);
    TEST(test_poly_system_eliminate_z);
    TEST(test_poly_system_eliminate_y_z);
    TEST(test_and_form_equivalent);
    TEST(test_single_symbol_var_shorthand);
    TEST(test_resultant_incompatible);
    TEST(test_common_root_condition);
    TEST(test_empty_ideal_returns_true);
    TEST(test_contradiction_returns_false);
    TEST(test_tautology_returns_true);
    TEST(test_inverse_function_prepass);
    TEST(test_inverse_function_two_vars);
    TEST(test_empty_var_list);
    TEST(test_eliminate_protected);
    TEST(test_memory_smoke);
    TEST(test_memory_smoke_larger);
    TEST(test_radical_square);
    TEST(test_radical_inv_pair);
    TEST(test_radical_cuberoot_pair);
    TEST(test_radical_dt_y_arc);
    TEST(test_radical_nested_sqrt);
    TEST(test_radical_quartic_with_d);
    TEST(test_radical_sqrtxp1_with_d);
    TEST(test_radical_memory_smoke);
    TEST(test_nlin_for_transcendental_arg);
    TEST(test_trig_integration_substitution);
    TEST(test_trig_pythagorean);
    TEST(test_trig_sin_squared);
    TEST(test_trig_double_angle);
    TEST(test_trig_tangent);
    TEST(test_trig_hyperbolic);
    TEST(test_trig_free_var_returns_true);
    TEST(test_trig_mixed_arg_nlin);
    TEST(test_exp_double_angle);
    TEST(test_exp_general_base);
    TEST(test_exp_integration_substitution);
    TEST(test_log_power);
    TEST(test_log_squared);
    TEST(test_trig_memory_smoke);

    printf("All Eliminate tests passed!\n");
    return 0;
}
