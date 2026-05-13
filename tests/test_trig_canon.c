/* Unit tests for the trig / hyperbolic ratio canonicalization performed
 * inside builtin_times and builtin_power.
 *
 * Coverage:
 *   - solo reciprocal rewrites: 1/Cos -> Sec, 1/Tan -> Cot, etc.
 *   - solo positive forms left alone: Sin[x], Tan[x], Sec[x]^2, ...
 *   - pair collapses: Sin/Cos -> Tan, Cos/Sin -> Cot
 *   - cancellations: Sin*Csc -> 1, Tan*Cot -> 1
 *   - product simplifications: Cos*Tan -> Sin, Tan*Csc -> Sec
 *   - powers: Sin^2/Cos^2 -> Tan^2, Cos^-3 -> Sec^3
 *   - leftovers: Sin^3/Cos -> Sin^2 * Tan, Sin/Cos^3 -> Tan * Sec^2
 *   - hyperbolic family in isolation (mirror of trig)
 *   - mixed args don't merge: Sin[x] * Cos[y] stays
 *   - non-integer exponents left alone
 *   - symbolic exponents left alone
 *   - numeric coefficients pass through cleanly
 *   - compound arguments (Sin[x+y]/Cos[x+y] -> Tan[x+y])
 *   - trig and hyperbolic in the same product remain separate
 */

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

typedef struct { const char* input; const char* expected; } Case;

static void run_cases(const char* group, const Case* cases) {
    for (int i = 0; cases[i].input; i++) {
        Expr* e = parse_expression(cases[i].input);
        ASSERT_MSG(e != NULL, "[%s] parse failed: %s", group, cases[i].input);
        Expr* r = evaluate(e);
        char* s = expr_to_string(r);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0,
                   "[%s] %s -> %s (expected %s)",
                   group, cases[i].input, s, cases[i].expected);
        free(s);
        expr_free(r);
        expr_free(e);
    }
}

static void run_cases_full(const char* group, const Case* cases) {
    for (int i = 0; cases[i].input; i++) {
        Expr* e = parse_expression(cases[i].input);
        ASSERT_MSG(e != NULL, "[%s] parse failed: %s", group, cases[i].input);
        Expr* r = evaluate(e);
        char* s = expr_to_string_fullform(r);
        ASSERT_MSG(strcmp(s, cases[i].expected) == 0,
                   "[%s/Full] %s -> %s (expected %s)",
                   group, cases[i].input, s, cases[i].expected);
        free(s);
        expr_free(r);
        expr_free(e);
    }
}

/* --------------------------------------------------------------------- */
/* Solo reciprocal rewrites                                              */

static void test_solo_reciprocal_trig(void) {
    Case cs[] = {
        {"1/Sin[x]", "Csc[x]"},
        {"1/Cos[x]", "Sec[x]"},
        {"1/Tan[x]", "Cot[x]"},
        {"1/Cot[x]", "Tan[x]"},
        {"1/Sec[x]", "Cos[x]"},
        {"1/Csc[x]", "Sin[x]"},
        {NULL, NULL}
    };
    run_cases("solo_reciprocal_trig", cs);
}

static void test_solo_reciprocal_hyp(void) {
    Case cs[] = {
        {"1/Sinh[x]", "Csch[x]"},
        {"1/Cosh[x]", "Sech[x]"},
        {"1/Tanh[x]", "Coth[x]"},
        {"1/Coth[x]", "Tanh[x]"},
        {"1/Sech[x]", "Cosh[x]"},
        {"1/Csch[x]", "Sinh[x]"},
        {NULL, NULL}
    };
    run_cases("solo_reciprocal_hyp", cs);
}

/* --------------------------------------------------------------------- */
/* Solo reciprocal with higher power                                     */

static void test_solo_reciprocal_powers(void) {
    Case cs[] = {
        {"1/Cos[x]^2",  "Sec[x]^2"},
        {"1/Sin[x]^3",  "Csc[x]^3"},
        {"1/Tan[x]^4",  "Cot[x]^4"},
        {"1/Cot[x]^2",  "Tan[x]^2"},
        {"1/Sec[x]^5",  "Cos[x]^5"},
        {"1/Csc[x]^2",  "Sin[x]^2"},
        {"1/Cosh[x]^3", "Sech[x]^3"},
        {"1/Tanh[x]^2", "Coth[x]^2"},
        {NULL, NULL}
    };
    run_cases("solo_reciprocal_powers", cs);
}

/* --------------------------------------------------------------------- */
/* Solo positive forms must remain unchanged                             */

static void test_solo_positive_unchanged(void) {
    Case cs[] = {
        {"Sin[x]",     "Sin[x]"},
        {"Cos[x]",     "Cos[x]"},
        {"Tan[x]",     "Tan[x]"},
        {"Cot[x]",     "Cot[x]"},
        {"Sec[x]",     "Sec[x]"},
        {"Csc[x]",     "Csc[x]"},
        {"Sin[x]^2",   "Sin[x]^2"},
        {"Cos[x]^3",   "Cos[x]^3"},
        {"Tan[x]^4",   "Tan[x]^4"},
        {"Sec[x]^5",   "Sec[x]^5"},
        {"Sinh[x]",    "Sinh[x]"},
        {"Cosh[x]^2",  "Cosh[x]^2"},
        {NULL, NULL}
    };
    run_cases("solo_positive_unchanged", cs);
}

/* --------------------------------------------------------------------- */
/* Pair collapses (canonical Tan / Cot)                                  */

static void test_pair_collapse_tan(void) {
    Case cs[] = {
        {"Sin[x]/Cos[x]",     "Tan[x]"},
        {"Cos[x]/Sin[x]",     "Cot[x]"},
        {"Sin[x] Cos[x]^-1",  "Tan[x]"},
        {"Sin[x]^2/Cos[x]^2", "Tan[x]^2"},
        {"Cos[x]^3/Sin[x]^3", "Cot[x]^3"},
        {NULL, NULL}
    };
    run_cases("pair_collapse_tan", cs);
}

static void test_pair_collapse_hyp(void) {
    Case cs[] = {
        {"Sinh[x]/Cosh[x]",     "Tanh[x]"},
        {"Cosh[x]/Sinh[x]",     "Coth[x]"},
        {"Sinh[x]^2/Cosh[x]^2", "Tanh[x]^2"},
        {NULL, NULL}
    };
    run_cases("pair_collapse_hyp", cs);
}

/* --------------------------------------------------------------------- */
/* Full cancellations                                                    */

static void test_full_cancel(void) {
    Case cs[] = {
        {"Sin[x] Csc[x]",     "1"},
        {"Cos[x] Sec[x]",     "1"},
        {"Tan[x] Cot[x]",     "1"},
        {"Sin[x]^3 Csc[x]^3", "1"},
        {"Cos[x]^2 Sec[x]^2", "1"},
        {"Tan[x]^4 Cot[x]^4", "1"},
        /* Hyperbolic */
        {"Sinh[x] Csch[x]",   "1"},
        {"Cosh[x] Sech[x]",   "1"},
        {"Tanh[x] Coth[x]",   "1"},
        {NULL, NULL}
    };
    run_cases("full_cancel", cs);
}

/* --------------------------------------------------------------------- */
/* Productive simplifications                                            */

static void test_product_simplify(void) {
    Case cs[] = {
        /* Sin*Sec = Sin/Cos = Tan */
        {"Sin[x] Sec[x]",     "Tan[x]"},
        /* Cos*Csc = Cos/Sin = Cot */
        {"Cos[x] Csc[x]",     "Cot[x]"},
        /* Tan*Csc = (Sin/Cos)/Sin = Sec */
        {"Tan[x] Csc[x]",     "Sec[x]"},
        /* Cot*Sec = (Cos/Sin)/Cos = Csc */
        {"Cot[x] Sec[x]",     "Csc[x]"},
        /* Sin*Cot = Sin*Cos/Sin = Cos */
        {"Sin[x] Cot[x]",     "Cos[x]"},
        /* Cos*Tan = Cos*Sin/Cos = Sin */
        {"Cos[x] Tan[x]",     "Sin[x]"},
        /* Tan*Sec = Sin/Cos^2 -- canonical via Tan family with leftover Sec.
         * Times has ATTR_ORDERLESS so the printed order is canonical (Sec
         * sorts before Tan: shorter+lexically-earlier head). */
        {"Tan[x] Sec[x]",     "Sec[x] Tan[x]"},
        {NULL, NULL}
    };
    run_cases("product_simplify", cs);
}

/* --------------------------------------------------------------------- */
/* Leftover (non-zero remainder) cases                                   */

static void test_leftover_powers(void) {
    Case cs[] = {
        /* Sin^3/Cos: a=3, b=-1, s=2. Sin^2 * Tan */
        {"Sin[x]^3/Cos[x]",  "Sin[x]^2 Tan[x]"},
        /* Sin/Cos^3: a=1, b=-3, s=-2. Tan * Sec^2 (sorted: Sec^2 first). */
        {"Sin[x]/Cos[x]^3",  "Sec[x]^2 Tan[x]"},
        /* Cos^3/Sin: a=-1, b=3, s=2. Cot * Cos^2 (sort: Cos^2 before Cot). */
        {"Cos[x]^3/Sin[x]",  "Cos[x]^2 Cot[x]"},
        /* Cos/Sin^3: a=-3, b=1, s=-2. Cot * Csc^2 */
        {"Cos[x]/Sin[x]^3",  "Cot[x] Csc[x]^2"},
        /* Hyperbolic mirror of one case */
        {"Sinh[x]^3/Cosh[x]", "Sinh[x]^2 Tanh[x]"},
        {NULL, NULL}
    };
    run_cases("leftover_powers", cs);
}

/* --------------------------------------------------------------------- */
/* Different arguments must NOT merge                                    */

static void test_different_args_kept_separate(void) {
    Case cs[] = {
        /* Different symbol args. Cos with exponent -1 still rewrites solo
         * to Sec[y]; Sin[x] just stays. Times then sorts canonically:
         * under the polynomial-aware comparator Sin[x] (deg 0 in y) sorts
         * before Cos[y]/Sec[y] (deg infinity in y). */
        {"Sin[x] Cos[y]",     "Sin[x] Cos[y]"},
        {"Sin[x]/Cos[y]",     "Sin[x] Sec[y]"},
        {"Sin[x] Sin[y]",     "Sin[x] Sin[y]"},
        /* Mixing trig and hyperbolic of the same symbol must NOT merge,
         * even though both have Sin/Cos building blocks. They live in
         * separate buckets (kind=trig vs kind=hyp). */
        {"Sin[x] Cosh[x]",    "Cosh[x] Sin[x]"},
        {"Sin[x]/Cosh[x]",    "Sech[x] Sin[x]"},
        /* Different compound args don't merge either. The standalone
         * Cos[x]^-1 still rewrites to Sec[x]. */
        {"Sin[x + 1]/Cos[x]", "Sec[x] Sin[1 + x]"},
        {NULL, NULL}
    };
    run_cases("different_args_kept_separate", cs);
}

/* --------------------------------------------------------------------- */
/* Compound arguments (function call as the trig argument)               */

static void test_compound_argument(void) {
    Case cs[] = {
        {"Sin[x + y]/Cos[x + y]", "Tan[x + y]"},
        {"Cos[2 x]/Sin[2 x]",     "Cot[2 x]"},
        {"1/Cos[x^2]",            "Sec[x^2]"},
        {"Sin[a + b] Csc[a + b]", "1"},
        {NULL, NULL}
    };
    run_cases("compound_argument", cs);
}

/* --------------------------------------------------------------------- */
/* Numeric coefficients combine cleanly                                  */

static void test_with_numeric_coefficients(void) {
    Case cs[] = {
        {"2 Sin[x]/Cos[x]",   "2 Tan[x]"},
        {"Sin[x]/(2 Cos[x])", "1/2 Tan[x]"},
        {"3/Cos[x]",          "3 Sec[x]"},
        {"3/Cos[x]^2",        "3 Sec[x]^2"},
        {"-Sin[x]/Cos[x]",    "-Tan[x]"},
        {NULL, NULL}
    };
    run_cases("with_numeric_coefficients", cs);
}

/* --------------------------------------------------------------------- */
/* Non-integer exponents must be left alone                              */

static void test_non_integer_exponents_skipped(void) {
    Case cs[] = {
        /* Half-integer: parsed as Power[Cos[x], 1/2]; not eligible. */
        {"Power[Cos[x], 1/2]",  "Sqrt[Cos[x]]"},
        {"Power[Cos[x], -1/2]", "1/Sqrt[Cos[x]]"},
        /* Mixed: integer part collapses, rational part stays. Sin[x] /
         * Cos[x]^(3/2) -- exponent is rational so trig_canon doesn't fire
         * on Cos[x]; Sin stays as Sin. The printer shows the negative
         * exponent as a denominator. */
        {"Sin[x]/Cos[x]^(3/2)", "Sin[x]/Cos[x]^(3/2)"},
        {NULL, NULL}
    };
    run_cases("non_integer_exponents_skipped", cs);
}

/* --------------------------------------------------------------------- */
/* Symbolic / non-numeric exponents must be left alone                   */

static void test_symbolic_exponent_skipped(void) {
    Case cs[] = {
        {"Cos[x]^n",     "Cos[x]^n"},
        {"1/Cos[x]^n",   "Cos[x]^(-n)"},
        {"Sin[x]/Cos[x]^k", "Sin[x] Cos[x]^(-k)"},
        {NULL, NULL}
    };
    run_cases("symbolic_exponent_skipped", cs);
}

/* --------------------------------------------------------------------- */
/* FullForm sanity: the rewritten exprs are real Sec/Csc/Cot nodes,      */
/* not pretty-print trickery.                                            */

static void test_fullform_structure(void) {
    Case cs[] = {
        {"1/Cos[x]",       "Sec[x]"},
        {"1/Sin[x]^2",     "Power[Csc[x], 2]"},
        {"Sin[x]/Cos[x]",  "Tan[x]"},
        {"1/Tanh[x]",      "Coth[x]"},
        {"Sin[x]/Cos[x]^3","Times[Power[Sec[x], 2], Tan[x]]"},
        {NULL, NULL}
    };
    run_cases_full("fullform_structure", cs);
}

/* --------------------------------------------------------------------- */
/* Power[Sec[x] + 1, -1] etc. must not be touched (base is not a trig    */
/* function, just contains one).                                         */

static void test_non_trig_base_skipped(void) {
    Case cs[] = {
        {"1/(1 + Cos[x])", "1/(1 + Cos[x])"},
        {"1/(Sin[x] Cos[x])", "Csc[x] Sec[x]"},  /* product reciprocal IS rewritten via Times */
        {NULL, NULL}
    };
    run_cases("non_trig_base_skipped", cs);
}

/* --------------------------------------------------------------------- */
/* Numeric arguments still simplify (the rewrite is symbolic; Sin[1] is  */
/* held symbolic so Sin[1]/Cos[1] -> Tan[1]).                            */

static void test_numeric_argument(void) {
    Case cs[] = {
        {"Sin[1]/Cos[1]", "Tan[1]"},
        {"1/Cos[1]",      "Sec[1]"},
        {NULL, NULL}
    };
    run_cases("numeric_argument", cs);
}

/* --------------------------------------------------------------------- */
/* Multiple independent buckets in the same Times                        */

static void test_multiple_independent_buckets(void) {
    Case cs[] = {
        /* Two separate args, each independently collapses. Times sorts the
         * resulting Tan[x] Tan[y] alphabetically. */
        {"(Sin[x]/Cos[x]) (Sin[y]/Cos[y])", "Tan[x] Tan[y]"},
        /* Trig + Hyperbolic on same arg, both independently collapse. */
        {"(Sin[x]/Cos[x]) (Sinh[x]/Cosh[x])", "Tan[x] Tanh[x]"},
        {NULL, NULL}
    };
    run_cases("multiple_independent_buckets", cs);
}

/* --------------------------------------------------------------------- */
/* No spurious rewrites at top-level outside Times/Power                 */

static void test_no_rewrite_outside_times_power(void) {
    Case cs[] = {
        /* List elements are NOT inside a Times. */
        {"{Sin[x], Cos[x], Tan[x]}",  "{Sin[x], Cos[x], Tan[x]}"},
        /* Pythagorean identity is a Plus, not a Times. Stays as-is. */
        {"Sin[x]^2 + Cos[x]^2",       "Cos[x]^2 + Sin[x]^2"},
        {NULL, NULL}
    };
    run_cases("no_rewrite_outside_times_power", cs);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_solo_reciprocal_trig);
    TEST(test_solo_reciprocal_hyp);
    TEST(test_solo_reciprocal_powers);
    TEST(test_solo_positive_unchanged);
    TEST(test_pair_collapse_tan);
    TEST(test_pair_collapse_hyp);
    TEST(test_full_cancel);
    TEST(test_product_simplify);
    TEST(test_leftover_powers);
    TEST(test_different_args_kept_separate);
    TEST(test_compound_argument);
    TEST(test_with_numeric_coefficients);
    TEST(test_non_integer_exponents_skipped);
    TEST(test_symbolic_exponent_skipped);
    TEST(test_fullform_structure);
    TEST(test_non_trig_base_skipped);
    TEST(test_numeric_argument);
    TEST(test_multiple_independent_buckets);
    TEST(test_no_rewrite_outside_times_power);

    printf("All trig_canon tests passed.\n");
    return 0;
}
