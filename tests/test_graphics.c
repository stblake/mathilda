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
#ifdef USE_GRAPHICS
#include "render.h"
#endif
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

/* Extract the AspectRatio value Plot embeds in the returned Graphics. */
#define ASPECT(plot) \
    "First[Cases[" plot ", (AspectRatio -> v_) -> v]]"
#define IMGSIZE(plot) \
    "First[Cases[" plot ", (ImageSize -> v_) -> v]]"

void test_plot_aspect_ratio_default(void) {
    /* Plot defaults to AspectRatio -> 1/GoldenRatio, injected as a real. */
    assert_eval_eq(
        "Abs[" ASPECT("Plot[Sin[x], {x, 0, 1}]") " - N[1/GoldenRatio]] < 10^-6",
        "True", 0);
}

void test_plot_aspect_ratio_explicit_number(void) {
    /* An explicit ratio is numericalized to a real (integers and rationals
     * alike), so the renderer -- which has no evaluator -- gets a number. */
    assert_eval_eq(
        "Abs[" ASPECT("Plot[Sin[x], {x, 0, 1}, AspectRatio -> 2]") " - 2] < 10^-6",
        "True", 0);
    assert_eval_eq(
        "Abs[" ASPECT("Plot[Sin[x], {x, 0, 1}, AspectRatio -> 3/4]") " - 0.75] < 10^-6",
        "True", 0);
    assert_eval_eq(
        "Abs[" ASPECT("Plot[Sin[x], {x, 0, 1}, AspectRatio -> 0.4]") " - 0.4] < 10^-6",
        "True", 0);
}

void test_plot_aspect_ratio_symbolic_constant(void) {
    /* Symbolic-numeric ratios (1/GoldenRatio, GoldenRatio) resolve too. */
    assert_eval_eq(
        "Abs[" ASPECT("Plot[Sin[x], {x, 0, 1}, AspectRatio -> 1/GoldenRatio]")
        " - N[1/GoldenRatio]] < 10^-6", "True", 0);
    assert_eval_eq(
        "Abs[" ASPECT("Plot[Sin[x], {x, 0, 1}, AspectRatio -> GoldenRatio]")
        " - N[GoldenRatio]] < 10^-6", "True", 0);
}

void test_plot_aspect_ratio_symbol_settings(void) {
    /* Automatic and Full are interpreted by the renderer; they pass through
     * the option list verbatim rather than being numericalized. */
    assert_eval_eq(ASPECT("Plot[Sin[x], {x, 0, 1}, AspectRatio -> Automatic]"),
                    "Automatic", 0);
    assert_eval_eq(ASPECT("Plot[Sin[x], {x, 0, 1}, AspectRatio -> Full]"),
                    "Full", 0);
}

void test_plot_image_size_passthrough(void) {
    /* ImageSize is copied onto the Graphics unchanged for the renderer. */
    assert_eval_eq(IMGSIZE("Plot[Sin[x], {x, 0, 1}, ImageSize -> 600]"),
                    "600", 0);
    assert_eval_eq(IMGSIZE("Plot[Sin[x], {x, 0, 1}, ImageSize -> {400, 300}]"),
                    "{400, 300}", 0);
}

#ifdef USE_GRAPHICS
/* The window-shaping policy gfx_window_height() lives in render.c, which is
 * only compiled with a live Raylib; exercise it directly here. data_w/data_h
 * model a wide, short curve (Sin over a wide x-range). */
void test_window_height_policy(void) {
    const double dw = 12.0, dh = 2.0; /* data extent: ratio 1/6 */

    /* Explicit ratio shapes the window: height = round(width * a). */
    ASSERT(gfx_window_height(800, 600, 1.0, false, false, dw, dh) == 800);
    ASSERT(gfx_window_height(800, 600, 0.5, false, false, dw, dh) == 400);
    ASSERT(gfx_window_height(1000, 600, 0.618034, false, false, dw, dh) == 618);

    /* Automatic (aspect <= 0) follows the data ratio data_h/data_w. */
    ASSERT(gfx_window_height(1200, 600, -1.0, false, false, dw, dh) == 200);

    /* Full keeps the ImageSize box untouched (data stretches to fill it). */
    ASSERT(gfx_window_height(800, 600, 0.75, true, false, dw, dh) == 600);

    /* A pinned height (ImageSize -> {w,h}) wins over AspectRatio. */
    ASSERT(gfx_window_height(400, 400, 2.0, false, true, dw, dh) == 400);

    /* Extreme ratios clamp to the screen-friendly [100, 2000] band. */
    ASSERT(gfx_window_height(800, 600, 0.01, false, false, dw, dh) == 100);
    ASSERT(gfx_window_height(800, 600, 50.0, false, false, dw, dh) == 2000);
}
#endif

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
    TEST(test_plot_aspect_ratio_default);
    TEST(test_plot_aspect_ratio_explicit_number);
    TEST(test_plot_aspect_ratio_symbolic_constant);
    TEST(test_plot_aspect_ratio_symbol_settings);
    TEST(test_plot_image_size_passthrough);
#ifdef USE_GRAPHICS
    TEST(test_window_height_policy);
#endif

    printf("All graphics tests passed!\n");
    symtab_clear();
    return 0;
}
