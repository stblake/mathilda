/* Structural tests for StreamPlot[]'s StreamAnimate option. These never
 * open a real window -- they only inspect the returned Graphics[...] tree,
 * so no MATHILDA_NO_GRAPHICS_WINDOW guard is needed for these specific
 * assertions (StreamPlot itself never opens a window). */
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"
#include <stdlib.h>
#include <stdio.h>

void test_streamplot_default_emits_line_not_animated(void) {
    assert_eval_eq(
        "Length[Cases[StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, "
        "StreamPoints -> 3], _Line, Infinity]] > 0",
        "True", 0);
    assert_eval_eq(
        "Cases[StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, "
        "StreamPoints -> 3], _AnimatedStreamline, Infinity]",
        "{}", 0);
}

void test_streamplot_animate_true_emits_animated_not_line(void) {
    assert_eval_eq(
        "Length[Cases[StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, "
        "StreamPoints -> 3, StreamAnimate -> True], "
        "_AnimatedStreamline, Infinity]] > 0",
        "True", 0);
    assert_eval_eq(
        "Cases[StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, "
        "StreamPoints -> 3, StreamAnimate -> True], _Line, Infinity]",
        "{}", 0);
}

void test_streamplot_animate_false_matches_default(void) {
    assert_eval_eq(
        "Cases[StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, "
        "StreamPoints -> 3, StreamAnimate -> False], "
        "_AnimatedStreamline, Infinity]",
        "{}", 0);
}

void test_animatedstreamline_is_registered(void) {
    assert_eval_eq("Head[AnimatedStreamline[{{0,0},{1,1}}]]", "AnimatedStreamline", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_streamplot_default_emits_line_not_animated);
    TEST(test_streamplot_animate_true_emits_animated_not_line);
    TEST(test_streamplot_animate_false_matches_default);
    TEST(test_animatedstreamline_is_registered);

    printf("All streamplot tests passed.\n");
    return 0;
}
