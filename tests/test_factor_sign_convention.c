/* test_factor_sign_convention.c
 *
 * Factor writes every polynomial factor with a positive leading coefficient,
 * leading meaning lexicographic over the variables in order of first
 * appearance in the input (the order the reference's Variables[] reports),
 * and puts the sign into the numeric content; a lone irreducible factor of
 * unit content keeps its form.  Mathilda kept the constant term positive and
 * passed an already-factored product through unchanged, so FeynCalc's traces
 * came back as 4 (...) where the reference gives -4 (...) (DiracTrace,
 * PauliTrace, Tdec, FourSeries: 11 tests).  Expected values computed by the
 * reference kernel (compared in Mathilda's output form; the values agree).
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>

static void test_factors_normalised(void) {
    assert_eval_eq("Factor[(3 - D) s]", "-(-3 + D) s", 0);
    assert_eval_eq("Factor[-2 (1 - D) k]", "2 (-1 + D) k", 0);
    assert_eval_eq("Factor[(1 - 4 x) (1 + 3 y)]", "-(-1 + 4 x) (1 + 3 y)", 0);
    assert_eval_eq("Factor[x (1 - y) (2 - z)]", "x (-1 + y) (-2 + z)", 0);
    assert_eval_eq("Factor[4/3 (1 - D) (2 - D) (3 - D)]", "-4/3 (-3 + D) (-2 + D) (-1 + D)", 0);
    assert_eval_eq("Factor[1 - x^2]", "-(-1 + x) (1 + x)", 0);
    assert_eval_eq("Factor[12 - 16 eps]", "-4 (-3 + 4 eps)", 0);
    assert_eval_eq("Factor[(2 - D)/(3 - D)]", "(-2 + D)/(-3 + D)", 0);
    assert_eval_eq("Factor[(1 - x)/(1 + x)]", "-(-1 + x)/(1 + x)", 0);
}

static void test_lead_variable_is_first_appearance(void) {
    assert_eval_eq("Factor[-2 a c + 2 b]", "2 (b - a c)", 0);
    assert_eval_eq("Factor[-3 x^2 + 3 y - 3 w z]", "-3 (x^2 - y + w z)", 0);
    assert_eval_eq("Factor[-4 SP[t, z] LC[][c, d, k, p] + 4 SP[k, p] LC[][c, d, t, z]"
                   " + 4 SP[c, t] LC[][d, k, p, z]]",
                   "-4 (SP[t, z] LC[][c, d, k, p] - SP[k, p] LC[][c, d, t, z] - SP[c, t] LC[][d, k, p, z])", 0);
}

static void test_lone_irreducible_kept(void) {
    assert_eval_eq("Factor[1 - x]", "1 - x", 0);
    assert_eval_eq("Factor[-a - b]", "-a - b", 0);
    assert_eval_eq("Factor[-x + y^2]", "-x + y^2", 0);
    assert_eval_eq("Factor[Sqrt[2] - x]", "Sqrt[2] - x", 0);
    assert_eval_eq("Factor[-2 a + 2 b]", "-2 (a - b)", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_factors_normalised);
    TEST(test_lead_variable_is_first_appearance);
    TEST(test_lone_irreducible_kept);
    printf("All Factor sign-convention tests passed.\n");
    return 0;
}
