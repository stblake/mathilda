/* test_findintegernullvector.c -- unit tests for FindIntegerNullVector
 * (PSLQ / integer-relation detection).
 *
 * Relations are canonical only up to sign and (for complex input) a
 * Gaussian-unit multiple, so found relations are validated by the residual
 * a . x rather than by matching a fixed vector:
 *   - exact input   -> PossibleZeroQ[a . x] == True  (or Expand == 0),
 *   - inexact input -> Abs[a . x] below the input tolerance.
 * A few stable minimal-polynomial vectors are pinned exactly.
 *
 * Coverage: exact real / complex relations, inexact relations, minimal
 * polynomials, the norm-bound and no-relation diagnostics (returned
 * unevaluated), the WorkingPrecision option, $MaxExtraPrecision, and the
 * argument-validation edge cases.
 */
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

/* Evaluate `input`, assert its FullForm equals `expected`. */
static void expect(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); }
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
        ASSERT(0);
    } else {
        printf("PASS: %s -> %s\n", input, s);
    }
    free(s);
    expr_free(res);
    expr_free(e);
}

/* Exact real relations -- validated by PossibleZeroQ on the residual. */
static void test_exact_real(void) {
    expect("PossibleZeroQ[FindIntegerNullVector[{Log[2],Log[4]}] . {Log[2],Log[4]}]", "True");
    expect("PossibleZeroQ[FindIntegerNullVector[{Log[2],Log[3],Log[72]}] . {Log[2],Log[3],Log[72]}]", "True");
    expect("PossibleZeroQ[FindIntegerNullVector[{Pi,ArcTan[1/5],ArcTan[1/239]}] . {Pi,ArcTan[1/5],ArcTan[1/239]}]", "True");
    expect("PossibleZeroQ[FindIntegerNullVector[{Sin[1]^3,Sin[1],Sin[3]}] . {Sin[1]^3,Sin[1],Sin[3]}]", "True");
    expect("PossibleZeroQ[FindIntegerNullVector[{Cos[4],Cos[1]^4,Cos[1]^2,1}] . {Cos[4],Cos[1]^4,Cos[1]^2,1}]", "True");
}

/* Minimal polynomials of algebraic numbers (exact, stable vectors). */
static void test_minimal_polynomial(void) {
    /* coefficients of the minimal polynomial of Sqrt[2] + 3^(1/3) */
    expect("PossibleZeroQ[Module[{a=Sqrt[2]+3^(1/3)}, FindIntegerNullVector[a^Range[0,6]] . a^Range[0,6]]]", "True");
    /* 1, Sin[Pi/8], ... , Sin[Pi/8]^4 : a stable degree-4 relation */
    expect("FindIntegerNullVector[{1,Sin[Pi/8],Sin[Pi/8]^2,Sin[Pi/8]^3,Sin[Pi/8]^4}]",
           "List[-1, 0, 8, 0, -8]");
}

/* Complex input -> Gaussian-integer relation. */
static void test_complex(void) {
    expect("Expand[FindIntegerNullVector[{1,2I+Sqrt[3],(2I+Sqrt[3])^2}] . {1,2I+Sqrt[3],(2I+Sqrt[3])^2}]", "0");
    expect("Length[FindIntegerNullVector[{1,2I+Sqrt[3],(2I+Sqrt[3])^2}]]", "3");
}

/* Inexact input: the relation holds to the precision of the input. */
static void test_inexact(void) {
    expect("Abs[N[{Log[2],Log[4]}] . FindIntegerNullVector[N[{Log[2],Log[4]}]]] < 1/100000", "True");
    /* A 20-digit approximation of Sqrt[2]+Sqrt[3] recovers its quartic. */
    expect("Module[{b=N[Sqrt[2]+Sqrt[3],20]}, FindIntegerNullVector[b^Range[0,4]]]",
           "List[1, 0, -10, 0, 1]");
}

/* Norm-bound argument. */
static void test_norm_bound(void) {
    /* sqrt(5) ~ 2.236 <= 3: the relation is returned. */
    expect("PossibleZeroQ[FindIntegerNullVector[{Log[2],Log[4]},3] . {Log[2],Log[4]}]", "True");
    expect("Head[FindIntegerNullVector[{Log[2],Log[4]},3]]", "List");
    /* d = 1 is below the certified bound: no relation, left unevaluated. */
    expect("Head[FindIntegerNullVector[{Log[2],Log[4]},1]]", "FindIntegerNullVector");
    /* E and Pi: no relation up to a large bound. */
    expect("Head[FindIntegerNullVector[{E,Pi},1000000]]", "FindIntegerNullVector");
}

/* No integer relation -> returned unevaluated (rnfu). */
static void test_no_relation(void) {
    expect("Head[FindIntegerNullVector[{E,Pi}]]", "FindIntegerNullVector");
    /* Sqrt[2]+3^(1/3) is degree 6: {1, a, a^2, a^3, a^4} has no relation. */
    expect("Head[Module[{a=Sqrt[2]+3^(1/3)}, FindIntegerNullVector[a^Range[0,4]]]]",
           "FindIntegerNullVector");
}

/* WorkingPrecision option and $MaxExtraPrecision. */
static void test_precision(void) {
    /* default precision is not enough for this degree-30 relation (rnfu) ... */
    expect("Head[FindIntegerNullVector[Table[(2^(1/6)+3^(1/5))^i,{i,0,30}]]]",
           "FindIntegerNullVector");
    /* ... but WorkingPrecision -> 300 finds it (a length-31 vector). */
    expect("Length[FindIntegerNullVector[Table[(2^(1/6)+3^(1/5))^i,{i,0,30}],WorkingPrecision->300]]",
           "31");
    /* the new system variable exists and defaults to 50 */
    expect("$MaxExtraPrecision == 50", "True");
}

/* Argument validation: every malformed call is left unevaluated. */
static void test_errors(void) {
    expect("Head[FindIntegerNullVector[]]", "FindIntegerNullVector");
    expect("Head[FindIntegerNullVector[{}]]", "FindIntegerNullVector");
    expect("Head[FindIntegerNullVector[{Log[2]}]]", "FindIntegerNullVector");   /* length 1 */
    expect("Head[FindIntegerNullVector[{1,x}]]", "FindIntegerNullVector");      /* symbolic */
}

int main(void) {
    symtab_init();
    core_init();

    test_exact_real();
    test_minimal_polynomial();
    test_complex();
    test_inexact();
    test_norm_bound();
    test_no_relation();
    test_precision();
    test_errors();

    printf("\nAll FindIntegerNullVector tests passed.\n");
    return 0;
}
