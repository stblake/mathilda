/* test_risch_logderiv_radical.c — Risch`LogarithmicDerivativeOfRadical
 * (Bronstein §5.12 / eq. 7.44, the structure-theorem radical test Cor. 9.3.1/
 * 9.3.2 (ii)).
 *
 * Given f and a monomial tower, the builtin decides whether f is the logarithmic
 * derivative of a K-radical: it returns {n, u} with n a positive integer and u a
 * witness radical satisfying D[u]/u == n f exactly, or False.  Only the Exp/Log
 * (E u L) monomials participate (Cor. 9.3.1(ii)); the witness base is the log
 * ARGUMENT a_i for t_i = log(a_i) (supplied as the optional 4th monomial element)
 * and the monomial t_i itself for an exponential.
 *
 * The load-bearing check is the radical IDENTITY D[u]/u - n f == 0.  We verify it
 * end-to-end by substituting each tower symbol back to its concrete kernel
 * (t -> Log[x], t2 -> Exp[x], ...) so ordinary D[] differentiates it — an exact,
 * surface-form-independent certificate.
 *
 * Monomial derivatives (as in the real-structure suite):
 *   Log t=log(a):  Dt = Da/a  (t=log x -> 1/x)    Exp t=exp(a): Dt = Da·t (t=exp x -> t)
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0)
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Assert the radical identity D[u]/u - n f == 0 for the returned {n,u}, with the
 * tower symbols substituted back to concrete kernels `subs` so D[] closes it. */
static void assert_radical_identity(const char* f, const char* mons, const char* subs) {
    char buf[2048];
    snprintf(buf, sizeof buf,
        "Simplify[With[{radr = Risch`LogarithmicDerivativeOfRadical[%s, x, %s], "
        "rads = %s}, D[(Last[radr]) /. rads, x]/((Last[radr]) /. rads) "
        "- (First[radr]) ((%s) /. rads)]]",
        f, mons, subs, f);
    run_test(buf, "0");
}

/* Assert the builtin returns exactly `expected` (structural witness pin). */
static void assert_witness(const char* f, const char* mons, const char* expected) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Risch`LogarithmicDerivativeOfRadical[%s, x, %s]", f, mons);
    run_test(buf, expected);
}

/* ---- Witness values (structural) ----------------------------------------- */
static void test_witness_values(void) {
    /* D[a]/a with a = x^2: f = 2/x -> n=1, u=x^2. */
    assert_witness("2/x", "{{t, \"Log\", 1/x, x}}", "List[1, Power[x, 2]]");
    /* Fractional coefficient forces n: f = 1/(2x) -> n=2, u=x. */
    assert_witness("1/(2 x)", "{{t, \"Log\", 1/x, x}}", "List[2, x]");
    /* Pure exponential: f = 1 = Dt/t (t=exp x) -> n=1, u=t. */
    assert_witness("1", "{{t, \"Exp\", t}}", "List[1, t]");
    /* Zero f -> the trivial radical {1, 1}. */
    assert_witness("0", "{{t, \"Log\", 1/x, x}}", "List[1, 1]");
}

/* ---- The radical identity D[u]/u == n f (end-to-end diff-back) ------------ */
static void test_radical_identity(void) {
    assert_radical_identity("2/x", "{{t, \"Log\", 1/x, x}}", "{t -> Log[x]}");
    assert_radical_identity("1/(2 x)", "{{t, \"Log\", 1/x, x}}", "{t -> Log[x]}");
    assert_radical_identity("5/x", "{{t, \"Log\", 1/x, x}}", "{t -> Log[x]}");
    assert_radical_identity("1", "{{t, \"Exp\", t}}", "{t -> Exp[x]}");
    assert_radical_identity("2", "{{t, \"Exp\", t}}", "{t -> Exp[x]}");
    /* Mixed E u L tower: f = 3/x + 2 -> u = t2^2 x^3 (t1=log x, t2=exp x). */
    assert_radical_identity("3/x + 2", "{{t1, \"Log\", 1/x, x}, {t2, \"Exp\", t2}}",
                            "{t1 -> Log[x], t2 -> Exp[x]}");
    /* Two logs t1=log x, t2=log(x+1): f = 1/x + 1/(x+1) -> u = x (x+1). */
    assert_radical_identity("1/x + 1/(x + 1)",
        "{{t1, \"Log\", 1/x, x}, {t2, \"Log\", 1/(x + 1), x + 1}}",
        "{t1 -> Log[x], t2 -> Log[x + 1]}");
    /* Mixed signs / fractions: f = 2/x - 1/(3(x+1)). */
    assert_radical_identity("2/x - 1/(3 (x + 1))",
        "{{t1, \"Log\", 1/x, x}, {t2, \"Log\", 1/(x + 1), x + 1}}",
        "{t1 -> Log[x], t2 -> Log[x + 1]}");
}

/* ---- Non-radical declines (False) ---------------------------------------- */
static void test_declines(void) {
    /* 1/(x+1) is Dt of a NEW log log(x+1), not a radical of x. */
    assert_witness("1/(x + 1)", "{{t, \"Log\", 1/x, x}}", "False");
    /* f = x is not a Q-combination of any generator (not simple). */
    assert_witness("x", "{{t, \"Log\", 1/x, x}}", "False");
    /* An exp exponent that is not constant: f = 2x is not r·(Dt/t)=r. */
    assert_witness("2 x", "{{t, \"Exp\", t}}", "False");
}

/* ---- Scope boundary: the §5.12-only radical (log-monomial as a factor) ---- */
static void test_section512_only_scope(void) {
    /* Bronstein Example 7.3.1 finds 2(f-3w) = D(x^5 log x)/(x^5 log x), i.e.
     * f0 = 5/x + 1/(x log x) IS a logarithmic derivative of a radical whose
     * witness x^5·log(x) contains the LOG MONOMIAL itself as a factor.  The
     * structure-theorem test admits radicals built only from log ARGUMENTS and
     * exp monomials, so it (soundly) declines this §5.12-only case rather than
     * emitting a wrong witness.  Pinned so a future §5.12 recursion flips it.
     * (f0 has the term 1/(x t) = Dt/t for a LOGARITHM t, absent from eq. 7.44.) */
    assert_witness("5/x + 1/(x t)", "{{t, \"Log\", 1/x, x}}", "False");
}

/* ---- Robustness ---------------------------------------------------------- */
static void test_logderiv_robustness(void) {
    /* 3-tuple (no log argument): still decides; witness base falls back to t. */
    run_test("Head[Risch`LogarithmicDerivativeOfRadical[1/x]] "
             "=== Risch`LogarithmicDerivativeOfRadical", "True");    /* arity */
    /* Pure-exponential 3-tuple works fully (base = the exp monomial): the
     * coefficient 2 is already integral, so n=1 and u=t^2 (D[t^2]/t^2 = 2). */
    assert_witness("2", "{{t, \"Exp\", t}}", "List[1, Power[t, 2]]");
}

int main(void) {
    core_init();

    TEST(test_witness_values);
    TEST(test_radical_identity);
    TEST(test_declines);
    TEST(test_section512_only_scope);
    TEST(test_logderiv_robustness);

    printf("All risch_logderiv_radical tests passed.\n");
    return 0;
}
