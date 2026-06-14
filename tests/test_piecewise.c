#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void test_piecewise() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        // Integer inputs
        {"Round[5]", "5"},
        {"Floor[5]", "5"},
        {"Ceiling[5]", "5"},
        {"IntegerPart[5]", "5"},
        {"FractionalPart[5]", "0"},
        
        // Real inputs
        {"Floor[2.4]", "2"},
        {"Ceiling[2.4]", "3"},
        {"IntegerPart[2.4]", "2"},
        {"FractionalPart[2.4]", "0.4"},
        {"Round[2.4]", "2"},
        {"Round[2.5]", "2"},
        {"Round[3.5]", "4"},
        {"Floor[-2.4]", "-3"},
        {"Ceiling[-2.4]", "-2"},
        {"IntegerPart[-2.4]", "-2"},
        {"FractionalPart[-2.4]", "-0.4"},
        {"Round[-2.5]", "-2"},
        {"Round[-3.5]", "-4"},
        
        // Arbitrary-precision (MPFR) Real inputs — must reduce to exact integers
        // for Floor/Ceiling (the sibling MPFR coverage follows in later commits).
        {"Floor[N[3.7, 50]]", "3"},
        {"Floor[N[-3.7, 50]]", "-4"},
        {"Floor[N[5, 50]]", "5"},
        {"Ceiling[N[3.2, 50]]", "4"},
        {"Ceiling[N[-3.2, 50]]", "-3"},
        {"Ceiling[N[5, 50]]", "5"},
        {"Round[N[3.5, 50]]", "4"},
        {"Round[N[2.5, 50]]", "2"},
        {"Round[N[-3.5, 50]]", "-4"},
        {"Round[N[3.4, 50]]", "3"},
        {"IntegerPart[N[3.7, 50]]", "3"},
        {"IntegerPart[N[-3.7, 50]]", "-3"},
        {"IntegerPart[N[5, 50]]", "5"},
        // FractionalPart on MPFR preserves precision, so the exact bit pattern
        // is implementation-specific. Verify the value via a tolerance check
        // (and an exact integer-result case).
        {"Abs[FractionalPart[N[37/10, 50]] - 7/10] < 0.0001", "True"},
        {"Abs[FractionalPart[N[-37/10, 50]] + 7/10] < 0.0001", "True"},
        {"FractionalPart[N[5, 50]] == 0", "True"},

        // GMP BigInt inputs (overflow int64). All five ops must reduce.
        {"Floor[10^50]", "100000000000000000000000000000000000000000000000000"},
        {"Ceiling[10^50]", "100000000000000000000000000000000000000000000000000"},
        {"Round[10^50]", "100000000000000000000000000000000000000000000000000"},
        {"IntegerPart[10^50]", "100000000000000000000000000000000000000000000000000"},
        {"FractionalPart[10^50]", "0"},

        // Rational[BigInt, _] — exact GMP path. The old is_rational[int64]
        // gate refused to extract these, leaving the call symbolic.
        {"Floor[10^50/7]",          "14285714285714285714285714285714285714285714285714"},
        {"Ceiling[10^50/7]",        "14285714285714285714285714285714285714285714285715"},
        {"Round[10^50/7]",          "14285714285714285714285714285714285714285714285714"},
        {"IntegerPart[10^50/7]",    "14285714285714285714285714285714285714285714285714"},
        {"FractionalPart[10^50/7]", "Rational[2, 7]"},
        {"Floor[-(10^50)/7]",       "-14285714285714285714285714285714285714285714285715"},
        {"Ceiling[-(10^50)/7]",     "-14285714285714285714285714285714285714285714285714"},

        // Round banker's tie-breaking on BigInt rationals (half-to-even).
        {"Round[(10^50 + 1)/2]", "50000000000000000000000000000000000000000000000000"},
        {"Round[(10^50 - 1)/2]", "50000000000000000000000000000000000000000000000000"},
        {"Round[(10^50 + 3)/2]", "50000000000000000000000000000000000000000000000002"},

        // Rational inputs
        {"Floor[5/2]", "2"},
        {"Ceiling[5/2]", "3"},
        {"IntegerPart[5/2]", "2"},
        {"FractionalPart[5/2]", "Rational[1, 2]"},
        {"Round[5/2]", "2"},
        {"Round[7/2]", "4"},
        {"Floor[-5/2]", "-3"},
        {"Ceiling[-5/2]", "-2"},
        {"IntegerPart[-5/2]", "-2"},
        {"FractionalPart[-5/2]", "Rational[-1, 2]"},
        {"Round[-5/2]", "-2"},
        {"Round[-7/2]", "-4"},
        
        // Complex inputs
        {"Floor[Complex[2.4, -2.4]]", "Complex[2, -3]"},
        {"Round[Complex[2.5, -3.5]]", "Complex[2, -4]"},
        {"IntegerPart[Complex[5/2, -5/2]]", "Complex[2, -2]"},
        {"FractionalPart[Complex[5/2, -5/2]]", "Complex[Rational[1, 2], Rational[-1, 2]]"},
        {"FractionalPart[Complex[2, 3]]", "0"}, // Real 0, Im 0 -> integer 0
        
        // Sign extraction (symbolic)
        {"Floor[-x]", "Times[-1, Ceiling[x]]"},
        {"Ceiling[-x]", "Times[-1, Floor[x]]"},
        {"Round[-x]", "Times[-1, Round[x]]"},
        {"Floor[-2 x]", "Times[-1, Ceiling[Times[2, x]]]"},
        {"Ceiling[Times[-3, x, y]]", "Times[-1, Floor[Times[3, x, y]]]"},
        {"Floor[-x] + Ceiling[x]", "0"},

        // Idempotency / composition (symbolic)
        {"Floor[Floor[x]]", "Floor[x]"},
        {"Ceiling[Ceiling[x]]", "Ceiling[x]"},
        {"Round[Round[x]]", "Round[x]"},
        {"Floor[Ceiling[x]]", "Ceiling[x]"},
        {"Ceiling[Floor[x]]", "Floor[x]"},
        {"Floor[Round[x]]", "Round[x]"},
        {"Ceiling[Round[x]]", "Round[x]"},
        {"Round[Floor[x]]", "Floor[x]"},
        {"Round[Ceiling[x]]", "Ceiling[x]"},
        {"Ceiling[Floor[Ceiling[x]]]", "Ceiling[x]"},

        // IntegerPart / FractionalPart do NOT participate in these rules
        {"IntegerPart[-x]", "IntegerPart[Times[-1, x]]"},
        {"FractionalPart[Floor[x]]", "FractionalPart[Floor[x]]"},

        // Exact real numeric quantities (Pi, E, surds) — numericalized to
        // certified precision and reduced to the exact integer. Previously
        // these fell through unevaluated.
        {"Round[10000000*3^(2/3)]", "20800838"},
        {"Round[25000000000000000000 Pi]", "78539816339744830962"},
        {"Floor[Pi]", "3"},
        {"Ceiling[Pi]", "4"},
        {"Round[E]", "3"},
        {"Floor[100 E]", "271"},
        {"Ceiling[100 E]", "272"},
        {"IntegerPart[10000000*3^(2/3)]", "20800838"},
        {"IntegerPart[-(7 Pi)]", "-21"},
        {"Round[Sqrt[2]]", "1"},
        // Round of a huge exact multiple of Pi (51-digit answer).
        {"Round[10^50 * Pi]", "314159265358979323846264338327950288419716939937511"},
        // FractionalPart stays exact: x - IntegerPart[x].
        {"FractionalPart[10000000*3^(2/3)]", "Plus[-20800838, Times[10000000, Power[3, Rational[2, 3]]]]"},
        // Sign extraction composes with the numeric path.
        {"Round[-(3^(2/3))]", "-2"},
        {"Floor[-(3^(2/3))]", "-3"},
        // Genuinely symbolic / complex arguments stay unevaluated.
        {"Round[2.5 + y]", "Round[Plus[2.5, y]]"},
        {"Round[(-3)^(2/3)]", "Round[Power[-3, Rational[2, 3]]]"},

        // Two argument forms
        {"Floor[7, 3]", "6"},
        {"Ceiling[7, 3]", "9"},
        {"Round[7, 3]", "6"},
        {"Round[2.4, 0.5]", "2.5"},
        
        // Infinity inputs
        {"Floor[Infinity]", "Infinity"},
        {"Ceiling[Infinity]", "Infinity"},
        {"Round[Infinity]", "Infinity"},
        {"IntegerPart[Infinity]", "Infinity"},
        {"FractionalPart[Infinity]", "0"},
        {"Floor[-Infinity]", "Times[-1, Infinity]"},
        {"Ceiling[-Infinity]", "Times[-1, Infinity]"},
        {"Round[-Infinity]", "Times[-1, Infinity]"},
        {"IntegerPart[-Infinity]", "Times[-1, Infinity]"},
        {"FractionalPart[-Infinity]", "0"},
        
        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "Forward %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(e);
        expr_free(res);
    }
}

void test_unitstep() {
    struct {
        const char* input;
        const char* expected;
    } cases[] = {
        // --- Integer arguments. Value at 0 is 1. ---
        {"UnitStep[1]", "1"},
        {"UnitStep[0]", "1"},
        {"UnitStep[-1]", "0"},
        {"UnitStep[5]", "1"},
        {"UnitStep[-123]", "0"},

        // --- Real arguments: always an exact 0 or 1. ---
        {"UnitStep[.8]", "1"},
        {"UnitStep[0.0]", "1"},
        {"UnitStep[-0.001]", "0"},
        {"UnitStep[2.4]", "1"},
        {"UnitStep[-1.6]", "0"},
        {"UnitStep[3.200000000000]", "1"},

        // --- Rational arguments. ---
        {"UnitStep[1/2]", "1"},
        {"UnitStep[-3/7]", "0"},
        {"UnitStep[7/2]", "1"},

        // --- BigInt and BigInt-rational arguments (exact GMP path). ---
        {"UnitStep[10^50]", "1"},
        {"UnitStep[-(10^50)]", "0"},
        {"UnitStep[10^50/7]", "1"},
        {"UnitStep[-(10^50)/7]", "0"},

        // --- Exact symbolic real arguments resolved by certification. ---
        {"UnitStep[Pi]", "1"},
        {"UnitStep[-Pi]", "0"},
        {"UnitStep[Sqrt[2]]", "1"},
        {"UnitStep[E - 3]", "0"},
        {"UnitStep[3 - E]", "1"},
        {"UnitStep[Pi - 4]", "0"},
        // Tight boundary: Sqrt[2] - 99/70 ~ -6.4*10^-5 (certification must
        // separate it from zero and return 0, not guess).
        {"UnitStep[Sqrt[2] - 99/70]", "0"},
        {"UnitStep[99/70 - Sqrt[2]]", "1"},

        // --- Infinities (only the real positive point at infinity is >= 0). ---
        {"UnitStep[Infinity]", "1"},
        {"UnitStep[-Infinity]", "0"},

        // --- The empty form. ---
        {"UnitStep[]", "1"},

        // --- Multidimensional: 1 iff none are negative. ---
        {"UnitStep[1, 2, 3]", "1"},
        {"UnitStep[1, Pi, 5.3]", "1"},
        {"UnitStep[1, -2, 3]", "0"},
        {"UnitStep[-1, Pi]", "0"},
        {"UnitStep[Pi, Sqrt[2], 1/2]", "1"},
        // A negative arg wins even when others are unknown.
        {"UnitStep[x, -2, y]", "0"},

        // --- Mixed: proven-non-negative args are dropped (factor of 1). ---
        {"UnitStep[1, x]", "UnitStep[x]"},
        {"UnitStep[x, 2, y]", "UnitStep[x, y]"},

        // --- Purely symbolic / non-real arguments stay unevaluated. ---
        {"UnitStep[x]", "UnitStep[x]"},
        {"UnitStep[I]", "UnitStep[Complex[0, 1]]"},
        {"UnitStep[2 + 3 I]", "UnitStep[Complex[2, 3]]"},

        // --- Listable threading over a list argument. ---
        {"UnitStep[{-1.6, 3.200000000000}]", "List[0, 1]"},
        {"UnitStep[{-1, 0, 1}]", "List[0, 1, 1]"},

        // --- High-precision arguments resolve via their MPFR sign directly. ---
        {"UnitStep[1/7`100]", "1"},
        {"UnitStep[-1/7`100]", "0"},

        // --- Derivatives (product rule -> Piecewise). ---
        {"D[UnitStep[x], x]",
         "Piecewise[List[List[Indeterminate, Equal[x, 0]]], 0]"},
        {"D[UnitStep[x, y, z], z]",
         "Times[UnitStep[x, y], Piecewise[List[List[Indeterminate, Equal[z, 0]]], 0]]"},
        {"D[UnitStep[x^2], x]",
         "Times[2, x, Piecewise[List[List[Indeterminate, Equal[Power[x, 2], 0]]], 0]]"},
        // Derivative w.r.t. a variable the argument does not contain is 0.
        {"D[UnitStep[x], y]", "0"},

        {NULL, NULL}
    };

    for (int i = 0; cases[i].input != NULL; i++) {
        Expr* e = parse_expression(cases[i].input);
        Expr* res = evaluate(e);
        char* s = expr_to_string_fullform(res);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0, "UnitStep %s: expected %s, got %s", cases[i].input, cases[i].expected, s);
        free(s);
        expr_free(e);
        expr_free(res);
    }
}

int main() {
    symtab_init();
    core_init();

    TEST(test_piecewise);
    TEST(test_unitstep);

    return 0;
}
