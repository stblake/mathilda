/* test_int_rnb.c — RischNormanBlake: parallel Risch-Norman over a simple
 * radical extension L = K(y), y^m = q(x) (int_rnb.c).
 *
 * Correctness is checked by a high-precision NUMERICAL diff-back at an interior
 * rational point: a plain Simplify diff-back is unusable here because it is
 * pathologically slow on the nested radical fields (Sqrt[3/2] over Sqrt[x^2+1])
 * these integrals live in.  The engine itself gates on the same numerical check,
 * so a non-declined result is already verified; the assertions re-confirm it and
 * pin the four worked examples of Blake, "Parallel Integration over Simple
 * Radical Extensions".
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
#include <stdbool.h>

static bool eval_is(const char* s, const char* expected) {
    Expr* e = parse_expression(s);
    Expr* r = evaluate(e);
    char* got = expr_to_string_fullform(r);
    bool ok = got && strcmp(got, expected) == 0;
    if (!ok) printf("  [%s] -> %s (expected %s)\n", s, got ? got : "?", expected);
    free(got); expr_free(r); expr_free(e);
    return ok;
}

/* f integrates (non-declined) and D[result]-f vanishes numerically at x=53/100. */
static void assert_rnb(const char* f) {
    char buf[1600];
    snprintf(buf, sizeof(buf),
        "With[{r = Integrate`RischNormanBlake[%s, x]}, "
        "Head[r] =!= Integrate`RischNormanBlake && "
        "Abs[N[(D[r, x] - (%s)) /. x -> 53/100, 30]] < 1/10^12]", f, f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: RischNormanBlake diff-back nonzero/declined", f);
}

/* f is out of scope / non-elementary: the engine declines cleanly. */
static void assert_declines(const char* f) {
    char buf[800];
    snprintf(buf, sizeof(buf),
        "Head[Integrate`RischNormanBlake[%s, x]] === Integrate`RischNormanBlake", f);
    ASSERT_MSG(eval_is(buf, "True"), "%s: expected a clean decline", f);
}

/* Field setup: the Trager basis E_i and the derivation Lam_i (Info debug). */
static void test_rnb_field(void) {
    /* y^2 = x^2+1: E=(1,1), Lam=(0, x/(1+x^2)), Norm[w1] = -(1+x^2). */
    ASSERT_MSG(eval_is("Integrate`RNB`Info[1/Sqrt[x^2+1], x][[4]]", "List[1, 1]"),
        "y^2=x^2+1: E_i");
    ASSERT_MSG(eval_is("Simplify[Integrate`RNB`Info[1/Sqrt[x^2+1], x][[6]] + (1+x^2)]",
        "0"), "y^2=x^2+1: Norm[w1]");
    /* y^3 = x (from x^(2/3)): E=(1,1,1), Lam=(0, 1/(3x), 2/(3x)). */
    ASSERT_MSG(eval_is("Integrate`RNB`Info[1/x^(2/3), x][[4]]", "List[1, 1, 1]"),
        "y^3=x: E_i");
    ASSERT_MSG(eval_is("Simplify[Integrate`RNB`Info[1/x^(2/3), x][[5]] "
                       "- {0, 1/(3 x), 2/(3 x)}]", "List[0, 0, 0]"),
        "y^3=x: Lam_i");
}

/* The algebraic (rational-part-only) tier: no logand needed. */
static void test_rnb_algebraic(void) {
    assert_rnb("x/Sqrt[x^2+1]");        /* = Sqrt[x^2+1]                     */
    assert_rnb("2 x Sqrt[x^2+1]");      /* = (2/3)(x^2+1)^(3/2)              */
    assert_rnb("1/(x^2 Sqrt[x^2+1])");  /* = -Sqrt[x^2+1]/x                  */
    assert_rnb("x/Sqrt[1-x^2]");        /* = -Sqrt[1-x^2]  (lead c=-1 sign)  */
    assert_rnb("x^3/Sqrt[x^4+1]");
}

/* Units at infinity (m=2 continued fraction / Pell). */
static void test_rnb_units(void) {
    assert_rnb("1/Sqrt[x^2+1]");        /* = Log[x + Sqrt[x^2+1]]  (flagship) */
    assert_rnb("Sqrt[x^2+1]/x^2");      /* algebraic part + Log unit          */
    assert_rnb("1/Sqrt[2 x^2+3]");      /* non-unit leading coeff (c=2)       */
}

/* The four worked examples of the paper (residue-divisor / torsion logands). */
static void test_rnb_examples(void) {
    /* Example: genus 0, y^2=x^2+1, new constants Sqrt[2],Sqrt[3],Sqrt[6]. */
    assert_rnb("(2 x^2-5)/((x^2-2)(x^2-3) Sqrt[x^2+1]) + x/Sqrt[x^2+1] "
               "+ 1/(x^2 Sqrt[x^2+1])");
    /* Example: genus 1, y^2=x^4+1, a unit at infinity and a 2-torsion divisor. */
    assert_rnb("(x^2+1)/(x Sqrt[x^4+1]) + x^3/Sqrt[x^4+1] "
               "+ (x^4-1)/(x^2 Sqrt[x^4+1])");
    /* Example: m=3, y^3=x^2, nontrivial integral basis, three logs with omega. */
    assert_rnb("1/(x^(2/3)(x-1)) + x^(-4/3)");
    /* A single separated branch place (residue divisor over x^2-2). */
    assert_rnb("1/((x^2-2) Sqrt[x^2+1])");
}

/* Non-elementary / out-of-scope integrands decline cleanly (never wrong). */
static void test_rnb_declines(void) {
    assert_declines("1/Sqrt[x^4+1]");   /* lemniscatic elliptic (non-elementary) */
    assert_declines("x/Sqrt[x^3+1]");   /* elliptic (non-elementary)             */
    assert_declines("Sin[x]");          /* no radical of x -> not our business   */
}

int main(void) {
    core_init();
    TEST(test_rnb_field);
    TEST(test_rnb_algebraic);
    TEST(test_rnb_units);
    TEST(test_rnb_examples);
    TEST(test_rnb_declines);
    printf("All RischNormanBlake tests passed.\n");
    return 0;
}
