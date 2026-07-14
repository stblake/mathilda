/* test_risch_structure_real.c — the REAL Risch structure theorem (Bronstein
 * §9.3, Corollary 9.3.2, eqs. 9.14/9.15) and the four-way disjoint index-set
 * partitioning it requires.
 *
 * The complex structure theorem (Cor. 9.3.1) decides reducibility of a new
 * logarithm / exponential over a tower via a Q-linear span of the monomial
 * generators.  The REAL theorem (Cor. 9.3.2) adds two more decisions for a real
 * tower carrying tangents (T) and arc-tangents (A), and crucially splits the
 * tower monomials into FOUR DISJOINT index sets E/L/T/A: the log/exp decisions
 * (9.12/9.13) use only E∪L generators, while the tan/arctan decisions
 * (9.14/9.15) use only T∪A generators.
 *
 * Monomial derivatives used below (t = f(a), a the inner argument):
 *   Log t=log(a):     Dt = Da/a          e.g. t=log x  -> Dt = 1/x
 *   Exp t=exp(a):     Dt = Da·t          e.g. t=exp x  -> Dt = t
 *   Tan t=tan(a):     Dt = Da·(t^2+1)    e.g. t=tan x  -> Dt = 1+t^2
 *   ArcTan t=atan(a): Dt = Da/(a^2+1)    e.g. t=atan x -> Dt = 1/(1+x^2)
 *
 * Generators (Cor. 9.3.2):  Log,ArcTan -> Dt ; Exp -> Dt/t ; Tan -> Dt/(t^2+1).
 *
 * Verifies: Risch`TanReducible (eq. 9.15, target Db), Risch`ArcTanReducible
 * (eq. 9.14, target Db/(b^2+1)); the disjoint-index partitioning (a tan query
 * ignores log/exp generators and vice versa); back-compatibility of
 * Log/ExpReducible on pure log/exp towers; and 4-tuple decode tolerance.
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

/* ---- Tangent reducibility (Cor. 9.3.2 iv, eq. 9.15; target Db) ------------ */
static void test_tan_reducible(void) {
#define T1 "{{t1, \"Tan\", 1 + t1^2}}"                 /* t1 = tan x, generator 1 */
    /* tan(2x) = reducible with rational coefficient 2 (Db = 2 = 2·1). */
    run_test("Risch`TanReducible[2 x, x, " T1 "]", "List[2]");
    run_test("Risch`TanReducible[3 x, x, " T1 "]", "List[3]");
    run_test("Risch`TanReducible[x, x, " T1 "]", "List[1]");
    run_test("Risch`TanReducible[x/2, x, " T1 "]", "List[Rational[1, 2]]");
    /* tan(c) for constant c: Db = 0 -> the zero combination. */
    run_test("Risch`TanReducible[5, x, " T1 "]", "List[0]");
    /* Genuinely new: tan(x^2) has Db = 2x, not a constant multiple of 1. */
    run_test("Risch`TanReducible[x^2, x, " T1 "]", "False");
    run_test("Risch`TanReducible[x^3 + x, x, " T1 "]", "False");
#undef T1
}

/* ---- Arc-tangent reducibility (Cor. 9.3.2 iii, eq. 9.14; target Db/(b^2+1)) */
static void test_arctan_reducible(void) {
#define A1 "{{a1, \"ArcTan\", 1/(1 + x^2)}}"           /* a1 = arctan x, generator 1/(1+x^2) */
    /* arctan(x): Db/(b^2+1) = 1/(1+x^2) = 1·generator. */
    run_test("Risch`ArcTanReducible[x, x, " A1 "]", "List[1]");
    /* arctan(c) constant -> zero combination. */
    run_test("Risch`ArcTanReducible[7, x, " A1 "]", "List[0]");
    /* Genuinely new: arctan(2x) target = 2/(1+4x^2) is not r/(1+x^2). */
    run_test("Risch`ArcTanReducible[2 x, x, " A1 "]", "False");
    run_test("Risch`ArcTanReducible[x^2, x, " A1 "]", "False");
#undef A1
}

/* ---- Disjoint index sets: tan/arctan queries ignore log/exp generators ---- */
static void test_disjoint_index_tan(void) {
    /* Mixed real tower {t1=log x (L), t2=exp x (E), t3=tan x (T)}.  A tangent
     * query uses ONLY the T∪A generators: tan(x) reduces onto t3 alone. */
#define M3 "{{t1, \"Log\", 1/x}, {t2, \"Exp\", t2}, {t3, \"Tan\", 1 + t3^2}}"
    run_test("Risch`TanReducible[x, x, " M3 "]", "List[0, 0, 1]");
    run_test("Risch`TanReducible[2 x, x, " M3 "]", "List[0, 0, 2]");
    /* tan(x^2) still new even amid log/exp monomials. */
    run_test("Risch`TanReducible[x^2, x, " M3 "]", "False");
#undef M3
}

static void test_disjoint_index_arctan(void) {
    /* {t1=tan x (T), a1=arctan x (A)}: an arctan query uses the T∪A generators.
     * arctan(x): 1/(1+x^2) = 0·(tan gen 1) + 1·(arctan gen 1/(1+x^2)). */
#define TA "{{t1, \"Tan\", 1 + t1^2}, {a1, \"ArcTan\", 1/(1 + x^2)}}"
    run_test("Risch`ArcTanReducible[x, x, " TA "]", "List[0, 1]");
    /* A tangent query over the same tower: tan(x) -> t1 (slot 0). */
    run_test("Risch`TanReducible[x, x, " TA "]", "List[1, 0]");
#undef TA
}

static void test_disjoint_index_log(void) {
    /* A log/exp query must IGNORE tan/arctan generators (they slot to 0). */
#define LTA "{{t1, \"Log\", 1/x}, {t2, \"Tan\", 1 + t2^2}, {a1, \"ArcTan\", 1/(1 + x^2)}}"
    /* log(x^2) = 2 log x: only the log generator participates. */
    run_test("Risch`LogReducible[x^2, x, " LTA "]", "List[2, 0, 0]");
    /* exp(x) has no log/exp reduction here (Db = 1 not in span of {1/x}). */
    run_test("Risch`ExpReducible[x^2, x, " LTA "]", "False");
#undef LTA
}

/* ---- Back-compatibility: pure log/exp towers unchanged -------------------- */
static void test_backcompat_logexp(void) {
    run_test("Risch`LogReducible[x^2, x, {{t1, \"Log\", 1/x}}]", "List[2]");
    run_test("Risch`LogReducible[2 x, x, {{t1, \"Log\", 1/x}}]", "List[1]");
    run_test("Risch`LogReducible[x + 1, x, {{t1, \"Log\", 1/x}}]", "False");
    run_test("Risch`ExpReducible[2 x, x, {{t1, \"Exp\", t1}}]", "List[2]");
    run_test("Risch`ExpReducible[x^2, x, {{t1, \"Exp\", t1}}]", "False");
    /* Mixed log+exp (E∪L both participate). */
    run_test("Risch`ExpReducible[x + t1, x, {{t1, \"Log\", 1/x}, {t2, \"Exp\", t2}}]",
             "List[1, 1]");
}

/* ---- Deeper real towers (stress) ----------------------------------------- */
static void test_deep_real_tower(void) {
    /* Two tangents t1=tan x, t2=tan(3x) is degenerate (t2 reducible), so use two
     * INDEPENDENT tangents t1=tan x, t2=tan(exp x) with an exp monomial e=exp x.
     *   e = exp x: De = e ; t2 = tan(e): Dt2 = De·(t2^2+1) = e (t2^2+1). */
#define REAL4 "{{e, \"Exp\", e}, {t1, \"Tan\", 1 + t1^2}, {t2, \"Tan\", e (1 + t2^2)}, {a1, \"ArcTan\", 1/(1 + x^2)}}"
    /* tan(x) -> t1 only (slot 1 among {e,t1,t2,a1}). */
    run_test("Risch`TanReducible[x, x, " REAL4 "]", "List[0, 1, 0, 0]");
    /* tan(exp x) -> t2 (slot 2): b = e is the field generator for exp x, Db = e
     * matches t2's generator; the exp monomial itself is NOT a T∪A generator. */
    run_test("Risch`TanReducible[e, x, " REAL4 "]", "List[0, 0, 1, 0]");
    /* arctan(x) -> a1 (slot 3). */
    run_test("Risch`ArcTanReducible[x, x, " REAL4 "]", "List[0, 0, 0, 1]");
    /* A genuinely new tangent tan(x^2) declines. */
    run_test("Risch`TanReducible[x^2, x, " REAL4 "]", "False");
#undef REAL4
}

/* ---- Robustness: 4-tuple decode tolerance + arity ------------------------ */
static void test_real_robustness(void) {
    /* The optional 4th monomial element (the log argument, for radical witness)
     * is accepted and ignored by the reducibility decisions. */
    run_test("Risch`LogReducible[x^2, x, {{t1, \"Log\", 1/x, x}}]", "List[2]");
    run_test("Risch`TanReducible[2 x, x, {{t1, \"Tan\", 1 + t1^2, x}}]", "List[2]");
    /* Arity: wrong number of arguments stays unevaluated. */
    run_test("Head[Risch`TanReducible[x]] === Risch`TanReducible", "True");
    run_test("Head[Risch`ArcTanReducible[x, x]] === Risch`ArcTanReducible", "True");
}

int main(void) {
    core_init();

    TEST(test_tan_reducible);
    TEST(test_arctan_reducible);
    TEST(test_disjoint_index_tan);
    TEST(test_disjoint_index_arctan);
    TEST(test_disjoint_index_log);
    TEST(test_backcompat_logexp);
    TEST(test_deep_real_tower);
    TEST(test_real_robustness);

    printf("All risch_structure_real tests passed.\n");
    return 0;
}
