/* Symbolic-construction tests for Graphics[]/Show[]/Plot[]. These never
 * call graphics_show()'s real Raylib path -- MATHILDA_NO_GRAPHICS_WINDOW
 * forces the windowing call to no-op regardless of how USE_GRAPHICS
 * resolved for this build, so the suite stays headless everywhere. */
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"
#include <stdlib.h>
#include <stdio.h>

void test_graphics_literal_construction(void) {
    assert_eval_eq("FullForm[Graphics[{Point[{0,0}]}]]",
                    "Graphics[List[Point[List[0, 0]]]]", 0);
    assert_eval_eq("Head[Graphics[{Point[{0,0}]}]]", "Graphics", 0);
    assert_eval_eq("Graphics[{Point[{0,0}]}]", "-Graphics-", 0);
}

void test_plot_returns_graphics_head(void) {
    assert_eval_eq("Head[Plot[Sin[x], {x, 0, 2 Pi}]]", "Graphics", 0);
    assert_eval_eq("Head[Plot[x^2, {x, -1, 1}]]", "Graphics", 0);
}

void test_plot_invalid_args_stay_unevaluated(void) {
    assert_eval_eq("Plot[Sin[x], {x, a, b}]", "Plot[Sin[x], {x, a, b}]", 0);
    assert_eval_eq("Plot[Sin[x], {x, 0, 1}, PlotPoints -> 1]",
                    "Plot[Sin[x], {x, 0, 1}, PlotPoints -> 1]", 0);
}

void test_plot_honors_plot_points_option(void) {
    /* PlotPoints sets the *initial* sample count; MaxRecursion -> 0
     * disables all adaptive refinement, so the rendered Line[...] should
     * contain exactly that many points. */
    assert_eval_eq(
        "g = Plot[x, {x, 0, 1}, PlotPoints -> 7, MaxRecursion -> 0]; "
        "Length[g[[1]][[1]][[1]]]",
        "7", 0);
}

void test_show_requires_graphics_argument(void) {
    assert_eval_eq("Show[5]", "Show[5]", 0);
}

void test_show_merges_options(void) {
    assert_eval_eq("Show[Graphics[{Point[{0,0}]}], Axes -> True][[2]]",
                    "Axes -> True", 0);
}

int main(void) {
    setenv("MATHILDA_NO_GRAPHICS_WINDOW", "1", 1);
    symtab_init();
    core_init();

    TEST(test_graphics_literal_construction);
    TEST(test_plot_returns_graphics_head);
    TEST(test_plot_invalid_args_stay_unevaluated);
    TEST(test_plot_honors_plot_points_option);
    TEST(test_show_requires_graphics_argument);
    TEST(test_show_merges_options);

    printf("All graphics tests passed!\n");
    symtab_clear();
    return 0;
}
