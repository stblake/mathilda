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
    /* Mesh defaults to True (unlike Plot's None).  The interior-only mesh
     * draws each shared edge exactly once: for a 4x4 grid (9 cells) there
     * are 2x3 = 6 horizontal interior edges and 3x2 = 6 vertical interior
     * edges (domain perimeter edges are suppressed), so:
     *   9 polygons + 1 GrayLevel directive + 12 Line edges = 22. */
    assert_eval_eq(
        "Length[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 4,"
        " MaxRecursion -> 0][[1]]]",
        "22", 0);
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

void test_plot3d_mesh_interior_only(void) {
    /* Verify that no mesh lines are emitted at the domain perimeter.
     * For a 3x3 grid (2x2 = 4 cells): interior lines =
     *   bottom edges (j=1): 2 (i=0,1)
     *   left  edges (i=1): 2 (j=0,1)
     * Total = 4.  The 8 perimeter edges are never emitted. */
    assert_eval_eq(
        "Length[Cases[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 3,"
        " MaxRecursion -> 0][[1]], Line[___]]]",
        "4", 0);
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

void test_plot3d_color_function_custom_3arg(void) {
    /* A 3-arg ColorFunction f[x,y,z] is tried first (Plot3D's own convention).
     * Here we pass Function[{x,y,z}, Hue[z]] which returns Hue[...]. */
    assert_eval_eq(
        "Head[First[Cases[Plot3D[x + y, {x, 0, 1}, {y, 0, 1},"
        " ColorFunction -> Function[{x, y, z}, Hue[z]],"
        " Mesh -> None, PlotPoints -> 2, MaxRecursion -> 0][[1]],"
        " Hue[___]]]]",
        "Hue", 0);
    /* ColorFunction -> GrayLevel colours by height (z). Produces GrayLevel[...]. */
    assert_eval_eq(
        "Head[First[Cases[Plot3D[x^2 + y^2, {x, 0, 1}, {y, 0, 1},"
        " ColorFunction -> GrayLevel,"
        " Mesh -> None, PlotPoints -> 2, MaxRecursion -> 0][[1]],"
        " GrayLevel[___]]]]",
        "GrayLevel", 0);
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

void test_plot3d_plot_style_list(void) {
    /* PlotStyle -> {color1, color2} assigns one palette entry per surface.
     * For two surfaces with Mesh -> None, PlotPoints -> 2, MaxRecursion -> 0,
     * the prim list is: [color1, poly1, color2, poly2] = 4 elements. */
    assert_eval_eq(
        "Length[Plot3D[{x + y, x - y}, {x, -1, 1}, {y, -1, 1},"
        " PlotStyle -> {RGBColor[1, 0, 0], RGBColor[0, 0, 1]},"
        " Mesh -> None, PlotPoints -> 2, MaxRecursion -> 0][[1]]]",
        "4", 0);
    /* First prim is the first explicit PlotStyle color. */
    assert_eval_eq(
        "Plot3D[{x + y, x - y}, {x, -1, 1}, {y, -1, 1},"
        " PlotStyle -> {RGBColor[1, 0, 0], RGBColor[0, 0, 1]},"
        " Mesh -> None, PlotPoints -> 2, MaxRecursion -> 0][[1]][[1]]",
        "RGBColor[1, 0, 0]", 0);
    /* Third prim is the second explicit PlotStyle color. */
    assert_eval_eq(
        "Plot3D[{x + y, x - y}, {x, -1, 1}, {y, -1, 1},"
        " PlotStyle -> {RGBColor[1, 0, 0], RGBColor[0, 0, 1]},"
        " Mesh -> None, PlotPoints -> 2, MaxRecursion -> 0][[1]][[3]]",
        "RGBColor[0, 0, 1]", 0);
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

void test_plot3d_boundary_clipping_polygon_count(void) {
    /* Sutherland-Hodgman clipping: with a 5x5 grid on [-2,2]x[-2,2] and
     * RegionFunction x^2+y^2 < 4, the 3x3 inner block of grid vertices is
     * valid, giving 2x2 = 4 fully-valid cells AND 12 boundary cells (the
     * remaining 16 - 4 = 12 cells each have >= 1 valid corner).  All 12
     * boundary cells are clipped into partial polygons (nc >= 3), so the
     * total polygon count is 4 + 12 = 16. */
    assert_eval_eq(
        "Length[Cases[Plot3D[x + y, {x, -2, 2}, {y, -2, 2},"
        " RegionFunction -> Function[{x,y,z}, x^2 + y^2 < 4],"
        " PlotPoints -> 5, MaxRecursion -> 0, Mesh -> None][[1]],"
        " Polygon[___]]]",
        "16", 0);
    /* Without RegionFunction, all 16 cells are fully valid -- same count. */
    assert_eval_eq(
        "Length[Cases[Plot3D[x + y, {x, -2, 2}, {y, -2, 2},"
        " PlotPoints -> 5, MaxRecursion -> 0, Mesh -> None][[1]],"
        " Polygon[___]]]",
        "16", 0);
}

void test_plot3d_boundary_clipping_triangle_cells(void) {
    /* The 4 corner boundary cells (e.g. (0,0) with only corner (-1,-1)
     * inside the circle) produce triangular (3-vertex) clipped polygons.
     * Min over all polygon vertex counts must be 3. */
    assert_eval_eq(
        "Min[Map[Length[#[[1]]] &,"
        " Cases[Plot3D[x + y, {x, -2, 2}, {y, -2, 2},"
        " RegionFunction -> Function[{x,y,z}, x^2 + y^2 < 4],"
        " PlotPoints -> 5, MaxRecursion -> 0, Mesh -> None][[1]],"
        " Polygon[___]]]]",
        "3", 0);
    /* Fully-valid interior cells produce 4-vertex quads; the max
     * vertex count over all polygons is 4. */
    assert_eval_eq(
        "Max[Map[Length[#[[1]]] &,"
        " Cases[Plot3D[x + y, {x, -2, 2}, {y, -2, 2},"
        " RegionFunction -> Function[{x,y,z}, x^2 + y^2 < 4],"
        " PlotPoints -> 5, MaxRecursion -> 0, Mesh -> None][[1]],"
        " Polygon[___]]]]",
        "4", 0);
}

void test_plot3d_exclusion_style_smooth_lines(void) {
    /* Each boundary cell contributes exactly one Line[] segment between the
     * two crossing points of its clipped polygon — these chord segments
     * trace the RegionFunction boundary smoothly.  For 12 boundary cells
     * this gives 12 ExclusionStyle Line[] primitives (and no staircase lines
     * from the fully-valid pass, because fn_ok=true everywhere for x+y). */
    assert_eval_eq(
        "Length[Cases[Plot3D[x + y, {x, -2, 2}, {y, -2, 2},"
        " RegionFunction -> Function[{x,y,z}, x^2 + y^2 < 4],"
        " PlotPoints -> 5, MaxRecursion -> 0, Mesh -> None][[1]],"
        " Line[___]]]",
        "12", 0);
    /* Without RegionFunction there is no ExclusionStyle pass at all. */
    assert_eval_eq(
        "Cases[Plot3D[x + y, {x, -2, 2}, {y, -2, 2},"
        " PlotPoints -> 5, MaxRecursion -> 0, Mesh -> None][[1]],"
        " Line[___]]",
        "{}", 0);
}

void test_plot3d_exclusion_style_custom_color(void) {
    /* When ExclusionStyle is given, its color directive appears in the prim
     * list before the ExclusionStyle Line[] segments. */
    assert_eval_eq(
        "Length[Cases[Plot3D[x + y, {x, -2, 2}, {y, -2, 2},"
        " RegionFunction -> Function[{x,y,z}, x^2 + y^2 < 4],"
        " ExclusionStyle -> RGBColor[1, 0.5, 0],"
        " PlotPoints -> 5, MaxRecursion -> 0, Mesh -> None][[1]],"
        " RGBColor[1, 0.5, 0]]]",
        "1", 0);
    /* Without ExclusionStyle the boundary color defaults to GrayLevel[0.35]. */
    assert_eval_eq(
        "Length[Cases[Plot3D[x + y, {x, -2, 2}, {y, -2, 2},"
        " RegionFunction -> Function[{x,y,z}, x^2 + y^2 < 4],"
        " PlotPoints -> 5, MaxRecursion -> 0, Mesh -> None][[1]],"
        " GrayLevel[0.35`]]]",
        "1", 0);
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

void test_plot3d_clipped_polygon_vertices_are_3d(void) {
    /* Clipped boundary cell polygons also use {x,y,z} triples (not {x,y}).
     * Check the first polygon (which, for this test, is a boundary cell). */
    assert_eval_eq(
        "Length[First[Cases[Plot3D[x + y, {x, -2, 2}, {y, -2, 2},"
        " RegionFunction -> Function[{x,y,z}, x^2 + y^2 < 4],"
        " PlotPoints -> 5, MaxRecursion -> 0, Mesh -> None][[1]],"
        " Polygon[___]]][[1]][[1]]]",
        "3", 0);
}

void test_plot3d_region_function_2arg_form(void) {
    /* Plot's 2-arg RegionFunction f[x,y] is accepted by Plot3D via the
     * eval_region3 fallback path. */
    assert_eval_eq(
        "Head[Plot3D[x + y, {x, 0, 2}, {y, 0, 2},"
        " RegionFunction -> Function[{x, y}, x + y < 3],"
        " PlotPoints -> 3, MaxRecursion -> 0, Mesh -> None]]",
        "Graphics3D", 0);
}

void test_plot3d_plot_range_stored_in_options(void) {
    /* PlotRange -> {zmin, zmax} is accepted and stored in the options list. */
    assert_eval_eq(
        "First[Cases[Plot3D[x + y, {x, -1, 1}, {y, -1, 1},"
        " PlotRange -> {-1, 1}], (PlotRange -> v_) -> v]]",
        "{-1, 1}", 0);
}

void test_plot3d_plot_label_stored_in_options(void) {
    assert_eval_eq(
        "First[Cases[Plot3D[x + y, {x, -1, 1}, {y, -1, 1},"
        " PlotLabel -> \"My Surface\"], (PlotLabel -> v_) -> v]]",
        "\"My Surface\"", 0);
}

int main(void) {
    setenv("MATHILDA_NO_GRAPHICS_WINDOW", "1", 1);
    symtab_init();
    core_init();

    TEST(test_plot3d_returns_graphics3d_head);
    TEST(test_plot3d_invalid_args_stay_unevaluated);
    TEST(test_plot3d_honors_plot_points_and_max_recursion);
    TEST(test_plot3d_mesh_default_and_toggle);
    TEST(test_plot3d_mesh_interior_only);
    TEST(test_plot3d_color_function_rainbow);
    TEST(test_plot3d_color_function_custom_3arg);
    TEST(test_plot3d_plot_style_and_axes_defaults);
    TEST(test_plot3d_plot_style_list);
    TEST(test_plot3d_region_function_filters_cells);
    TEST(test_plot3d_boundary_clipping_polygon_count);
    TEST(test_plot3d_boundary_clipping_triangle_cells);
    TEST(test_plot3d_exclusion_style_smooth_lines);
    TEST(test_plot3d_exclusion_style_custom_color);
    TEST(test_plot3d_multi_surface_palette);
    TEST(test_plot3d_polygon_vertices_are_3d);
    TEST(test_plot3d_clipped_polygon_vertices_are_3d);
    TEST(test_plot3d_region_function_2arg_form);
    TEST(test_plot3d_plot_range_stored_in_options);
    TEST(test_plot3d_plot_label_stored_in_options);

    printf("All Plot3D tests passed!\n");
    symtab_clear();
    return 0;
}
