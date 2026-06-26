/* Symbolic-construction tests for Plot3D[]/Graphics3D[]. These never call
 * graphics3d_show()'s real Raylib path -- MATHILDA_NO_GRAPHICS_WINDOW forces
 * the windowing call to no-op regardless of how USE_GRAPHICS resolved for
 * this build, so the suite stays headless everywhere (mirrors test_graphics.c). */
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"
#include <stdlib.h>
#include <stdio.h>

void test_plot3d_returns_graphics3d_head(void) {
    assert_eval_eq("Head[Plot3D[x + y, {x, -1, 1}, {y, -1, 1}]]", "Graphics3D", 0);
    assert_eval_eq("Head[Plot3D[x^2 + y^2, {x, -1, 1}, {y, -1, 1}]]", "Graphics3D", 0);
}

void test_plot3d_invalid_args_stay_unevaluated(void) {
    /* Symbolic (non-numeric) bounds: declines, exactly like Plot. */
    assert_eval_eq("Plot3D[Sin[x + y], {x, a, b}, {y, -1, 1}]",
                    "Plot3D[Sin[x + y], {x, a, b}, {y, -1, 1}]", 0);
    /* PlotPoints < 2 is rejected the same way split_options rejects it for Plot. */
    assert_eval_eq("Plot3D[Sin[x + y], {x, 0, 1}, {y, 0, 1}, PlotPoints -> 1]",
                    "Plot3D[Sin[x + y], {x, 0, 1}, {y, 0, 1}, PlotPoints -> 1]", 0);
    /* Plot3D[{}, ...] -- nothing to draw. */
    assert_eval_eq("Plot3D[{}, {x, -1, 1}, {y, -1, 1}]",
                    "Plot3D[{}, {x, -1, 1}, {y, -1, 1}]", 0);
}

void test_plot3d_honors_plot_points_and_max_recursion(void) {
    /* PlotPoints sets the *initial* per-axis grid resolution; MaxRecursion
     * -> 0 disables all whole-grid refinement, so a PlotPoints -> 4 grid
     * has exactly (4-1)^2 = 9 quad cells -> 9 Polygon[] primitives (no
     * Mesh, no ColorFunction, so nothing else is emitted). */
    assert_eval_eq(
        "Length[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 4,"
        " MaxRecursion -> 0, Mesh -> None][[1]]]",
        "9", 0);
    assert_eval_eq(
        "Length[Cases[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 4,"
        " MaxRecursion -> 0, Mesh -> None][[1]], Polygon[___]]]",
        "9", 0);
}

void test_plot3d_mesh_default_and_toggle(void) {
    /* Mesh defaults to True (unlike Plot's None): each of the 9 cells above
     * also gets 4 wireframe Line[] edges, prefixed by one GrayLevel[]
     * directive -- 9 polygons + 1 color directive + 9*4 lines = 46. */
    assert_eval_eq(
        "Length[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 4,"
        " MaxRecursion -> 0][[1]]]",
        "46", 0);
    /* The mesh color directive sits right after all 9 polygon fills. */
    assert_eval_eq(
        "Head[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 4,"
        " MaxRecursion -> 0][[1]][[10]]]",
        "GrayLevel", 0);
    /* Mesh -> None removes both the directive and the wireframe lines. */
    assert_eval_eq(
        "Cases[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 4,"
        " MaxRecursion -> 0, Mesh -> None][[1]], Line[___]]",
        "{}", 0);
}

void test_plot3d_color_function_rainbow(void) {
    /* Each cell is prefixed by ColorFunction's resolved color, exactly as
     * Plot's per-segment ColorFunction does. "Rainbow" resolves to Hue[...]. */
    assert_eval_eq(
        "Head[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, ColorFunction -> \"Rainbow\","
        " Mesh -> None, PlotPoints -> 2, MaxRecursion -> 0][[1]][[1]]]",
        "Hue", 0);
    /* One color directive per cell -> 2x as many primitives as cells. */
    assert_eval_eq(
        "Length[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, ColorFunction -> \"Rainbow\","
        " Mesh -> None, PlotPoints -> 2, MaxRecursion -> 0][[1]]]",
        "2", 0);
}

void test_plot3d_plot_style_and_axes_defaults(void) {
    assert_eval_eq(
        "First[Cases[Plot3D[x, {x, -1, 1}, {y, -1, 1}], (PlotStyle -> v_) -> v]]",
        "RGBColor[0.2, 0.4, 0.8]", 0);
    assert_eval_eq(
        "First[Cases[Plot3D[x, {x, -1, 1}, {y, -1, 1}, PlotStyle -> Red],"
        " (PlotStyle -> v_) -> v]]",
        "RGBColor[1, 0, 0]", 0);
    assert_eval_eq(
        "First[Cases[Plot3D[x, {x, -1, 1}, {y, -1, 1}], (Axes -> v_) -> v]]",
        "True", 0);
}

void test_plot3d_region_function_filters_cells(void) {
    /* RegionFunction[x,y,z] (3-arg, Plot3D's own convention): a region
     * containing the whole sampled cell keeps it (no decline). */
    assert_eval_eq(
        "Head[Plot3D[x + y, {x, -0.5, 0.5}, {y, -0.5, 0.5},"
        " RegionFunction -> Function[{x, y, z}, x^2 + y^2 < 1],"
        " PlotPoints -> 2, MaxRecursion -> 0]]",
        "Graphics3D", 0);
    /* A region excluding every grid corner leaves no valid cells, so the
     * whole call declines -- exactly like Plot's RegionFunction-starves-
     * every-point case. */
    assert_eval_eq(
        "Plot3D[x + y, {x, -2, 2}, {y, -2, 2},"
        " RegionFunction -> Function[{x, y, z}, x^2 + y^2 < 1],"
        " PlotPoints -> 2, MaxRecursion -> 0]",
        "Plot3D[x + y, {x, -2, 2}, {y, -2, 2},"
        " RegionFunction -> Function[{x, y, z}, x^2 + y^2 < 1],"
        " PlotPoints -> 2, MaxRecursion -> 0]", 0);
    /* RegionFunction also accepts Plot's 2-arg f[x,y] form as a fallback. */
    assert_eval_eq(
        "Head[Plot3D[x + y, {x, 0.1, 1}, {y, -1, 1},"
        " RegionFunction -> Function[{x, y}, x > 0],"
        " PlotPoints -> 2, MaxRecursion -> 0]]",
        "Graphics3D", 0);
}

void test_plot3d_multi_surface_palette(void) {
    /* Two surfaces, one cell each (PlotPoints -> 2, MaxRecursion -> 0,
     * Mesh -> None): each surface is one palette color directive + one
     * Polygon[], in that order. */
    assert_eval_eq(
        "Length[Plot3D[{x + y, x - y}, {x, -1, 1}, {y, -1, 1}, Mesh -> None,"
        " PlotPoints -> 2, MaxRecursion -> 0][[1]]]",
        "4", 0);
    assert_eval_eq(
        "Head[Plot3D[{x + y, x - y}, {x, -1, 1}, {y, -1, 1}, Mesh -> None,"
        " PlotPoints -> 2, MaxRecursion -> 0][[1]][[1]]]",
        "RGBColor", 0);
    assert_eval_eq(
        "Head[Plot3D[{x + y, x - y}, {x, -1, 1}, {y, -1, 1}, Mesh -> None,"
        " PlotPoints -> 2, MaxRecursion -> 0][[1]][[3]]]",
        "RGBColor", 0);
}

void test_plot3d_polygon_vertices_are_3d(void) {
    /* Every Polygon[] vertex is a {x,y,z} triple, not {x,y} -- the one
     * structural difference from Plot's 2D Polygon usage (Filling). */
    assert_eval_eq(
        "Length[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 2,"
        " MaxRecursion -> 0, Mesh -> None][[1]][[1]][[1]][[1]]]",
        "3", 0);
}

int main(void) {
    setenv("MATHILDA_NO_GRAPHICS_WINDOW", "1", 1);
    symtab_init();
    core_init();

    TEST(test_plot3d_returns_graphics3d_head);
    TEST(test_plot3d_invalid_args_stay_unevaluated);
    TEST(test_plot3d_honors_plot_points_and_max_recursion);
    TEST(test_plot3d_mesh_default_and_toggle);
    TEST(test_plot3d_color_function_rainbow);
    TEST(test_plot3d_plot_style_and_axes_defaults);
    TEST(test_plot3d_region_function_filters_cells);
    TEST(test_plot3d_multi_surface_palette);
    TEST(test_plot3d_polygon_vertices_are_3d);

    printf("All Plot3D tests passed!\n");
    symtab_clear();
    return 0;
}
