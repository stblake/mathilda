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

void test_graphics_options_registered(void) {
    /* Options[Graphics] must list the options the renderer honours, not {}. */
    assert_eval_eq("Length[Options[Graphics]] > 0", "True", 0);
    assert_eval_eq("Options[Graphics, Axes]", "{Axes -> False}", 0);
    assert_eval_eq("Options[Graphics, PlotRange]",
                    "{PlotRange -> Automatic}", 0);
    assert_eval_eq("OptionValue[Graphics, Frame]", "False", 0);
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

void test_show_merges_frame_option(void) {
    /* Frame, like Axes, is copied verbatim onto the Show'd Graphics. */
    assert_eval_eq("Show[Graphics[{Point[{0,0}]}], Frame -> True][[2]]",
                    "Frame -> True", 0);
}

void test_plot_frame_passthrough(void) {
    /* Frame and its companion options pass through to the Graphics intact. */
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, Frame -> True], (Frame -> v_) -> v]]",
        "True", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, Frame -> True, FrameTicks -> None],"
        " (FrameTicks -> v_) -> v]]",
        "None", 0);
}

void test_plot_frame_suppresses_axes_default(void) {
    /* A frame replaces the interior axes cross, so Plot's Axes -> True default
     * is withheld when Frame -> True is supplied (matching Wolfram). */
    assert_eval_eq(
        "Cases[Plot[Sin[x], {x, 0, 1}, Frame -> True], (Axes -> _)]",
        "{}", 0);
    /* But Frame -> False leaves the axes default in place. */
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, Frame -> False], (Axes -> v_) -> v]]",
        "True", 0);
    /* And a plain plot keeps Axes -> True. */
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}], (Axes -> v_) -> v]]",
        "True", 0);
}

void test_axes_origin_passthrough(void) {
    /* Absent by default (the renderer auto-computes the origin then). */
    assert_eval_eq("Cases[Plot[Sin[x], {x, 0, 1}], (AxesOrigin -> _)]", "{}", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, AxesOrigin -> {1, 2}], (AxesOrigin -> v_) -> v]]",
        "{1, 2}", 0);
}

void test_axes_style_and_ticks_style_passthrough(void) {
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, AxesStyle -> RGBColor[1, 0, 0]],"
        " (AxesStyle -> v_) -> v]]",
        "RGBColor[1, 0, 0]", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, TicksStyle -> GrayLevel[0.5]],"
        " (TicksStyle -> v_) -> v]]",
        "GrayLevel[0.5]", 0);
}

void test_frame_label_and_rotate_label_passthrough(void) {
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, Frame -> True, FrameLabel -> {\"t\", \"y\"}],"
        " (FrameLabel -> v_) -> v]]",
        "{\"t\", \"y\"}", 0);
    /* RotateLabel defaults to absent (the renderer treats unset as True). */
    assert_eval_eq("Cases[Plot[Sin[x], {x, 0, 1}], (RotateLabel -> _)]", "{}", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, RotateLabel -> False], (RotateLabel -> v_) -> v]]",
        "False", 0);
}

void test_plot_range_padding_passthrough(void) {
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, PlotRangePadding -> None],"
        " (PlotRangePadding -> v_) -> v]]",
        "None", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, PlotRangePadding -> 0.2],"
        " (PlotRangePadding -> v_) -> v]]",
        "0.2", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, PlotRangePadding -> {0.1, 0.3}],"
        " (PlotRangePadding -> v_) -> v]]",
        "{0.1, 0.3}", 0);
}

void test_grid_lines_passthrough(void) {
    /* Absent (no grid) by default. */
    assert_eval_eq("Cases[Plot[Sin[x], {x, 0, 1}], (GridLines -> _)]", "{}", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, GridLines -> Automatic], (GridLines -> v_) -> v]]",
        "Automatic", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, GridLines -> {{0, 0.5, 1}, None}],"
        " (GridLines -> v_) -> v]]",
        "{{0, 0.5, 1}, None}", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, GridLinesStyle -> GrayLevel[0.8]],"
        " (GridLinesStyle -> v_) -> v]]",
        "GrayLevel[0.8]", 0);
}

void test_prolog_epilog_passthrough(void) {
    assert_eval_eq("Cases[Plot[Sin[x], {x, 0, 1}], (Prolog -> _) | (Epilog -> _)]", "{}", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, Prolog -> {GrayLevel[0.9], Disk[{0, 0}, 1]}],"
        " (Prolog -> v_) -> v]]",
        "{GrayLevel[0.9], Disk[{0, 0}, 1]}", 0);
    /* Plot evaluates each option's RHS once before storing it (it must,
     * since it's HoldAll), so a named color constant like Red resolves to
     * its RGBColor[...] value here -- same as it would for a non-Held
     * Graphics[]'s own arguments. */
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, Epilog -> {Red, Point[{0, 0}]}],"
        " (Epilog -> v_) -> v]]",
        "{RGBColor[1, 0, 0], Point[{0, 0}]}", 0);
}

void test_named_color_constants_resolve(void) {
    /* Bare Graphics[] isn't Held, so named colors resolve at construction
     * time regardless of Plot's evaluate-once-in-split_options fix. */
    assert_eval_eq("FullForm[Red]", "RGBColor[1, 0, 0]", 0);
    assert_eval_eq("FullForm[Graphics[{Blue, Point[{0,0}]}]]",
                    "Graphics[List[RGBColor[0, 0, 1], Point[List[0, 0]]]]", 0);
    /* And inside a Plot option, via the evaluate-once fix. */
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, PlotStyle -> Green], (PlotStyle -> v_) -> v]]",
        "RGBColor[0, 1, 0]", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, Background -> LightGray], (Background -> v_) -> v]]",
        "GrayLevel[0.85]", 0);
    /* ?Red inspects the symbol, not its value (Information is HoldFirst). */
    assert_eval_eq("StringTake[Information[Red], 3]", "\"Red\"", 0);
}

void test_hue_color_directive(void) {
    assert_eval_eq("Graphics[{Hue[0.5], Point[{0,0}]}]", "-Graphics-", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, PlotStyle -> Hue[0.5]], (PlotStyle -> v_) -> v]]",
        "Hue[0.5]", 0);
}

void test_color_function_builds_per_segment_colors(void) {
    /* Without ColorFunction, a smooth curve over a small range is a single
     * Line[] run. With it, each segment gets its own color directive +
     * 2-point Line[] pair, so the primitive count grows well past 1. */
    assert_eval_eq("Length[Plot[Sin[x], {x, 0, 1}][[1]]]", "1", 0);
    assert_eval_eq(
        "Length[Plot[Sin[x], {x, 0, 1}, ColorFunction -> \"Rainbow\"][[1]]] > 1",
        "True", 0);
    /* "Rainbow" resolves to a Hue[...] directive (the first primitive). */
    assert_eval_eq(
        "Head[Plot[Sin[x], {x, 0, 1}, ColorFunction -> \"Rainbow\"][[1]][[1]]]",
        "Hue", 0);
}

void test_filling_builds_polygon(void) {
    /* Default Filling -> Axis: Opacity[0.3], then the fill Polygon[], then
     * (since ColorFunction isn't also set) a colour restore, before the
     * curve's own Line[] outline. */
    assert_eval_eq("Head[Plot[Sin[x], {x, 0, 1}, Filling -> Axis][[1]][[1]]]", "Opacity", 0);
    assert_eval_eq("Head[Plot[Sin[x], {x, 0, 1}, Filling -> Axis][[1]][[2]]]", "Polygon", 0);
    /* No Filling -> no Polygon at all. */
    assert_eval_eq("Cases[Plot[Sin[x], {x, 0, 1}][[1]], Polygon[___]]", "{}", 0);
    /* An explicit FillingStyle is used directly (no Opacity bracketing). */
    assert_eval_eq(
        "Head[Plot[Sin[x], {x, 0, 1}, Filling -> Axis, FillingStyle -> Red][[1]][[1]]]",
        "RGBColor", 0);
}

/* Regression test for a real bug: a fill segment that crosses the baseline
 * (e.g. one sample point just above Filling -> Axis's y=0, the next just
 * below) used to become a single self-intersecting "bowtie" quad, which
 * render.c's triangle-fan Polygon fill turned into a stray sliver right at
 * the crossing (visible as a small misplaced triangle in an actual
 * screenshot). build_fill_quads now splits any baseline-crossing segment
 * into two plain triangles instead of one quad. */
void test_filling_splits_at_baseline_crossing(void) {
    /* x in [-0.1, 0.1] straddles Sin's zero at x=0: with just 2 plot
     * points (one negative y, one positive y) and recursion disabled, the
     * single sampled segment crosses the baseline, so it must become two
     * 3-vertex triangles, not one 4-vertex quad. */
    assert_eval_eq(
        "Cases[Plot[Sin[x], {x, -0.1, 0.1}, Filling -> Axis, MaxRecursion -> 0,"
        " PlotPoints -> 2][[1]], Polygon[pts_] :> Length[pts]]",
        "{3, 3}", 0);
    /* x in [0.1, 0.3] stays entirely positive: the single segment doesn't
     * cross the baseline, so it's the plain 4-vertex quad. */
    assert_eval_eq(
        "Cases[Plot[Sin[x], {x, 0.1, 0.3}, Filling -> Axis, MaxRecursion -> 0,"
        " PlotPoints -> 2][[1]], Polygon[pts_] :> Length[pts]]",
        "{4}", 0);
}

void test_plot_legends_metadata(void) {
    /* Absent by default. */
    assert_eval_eq("Cases[Plot[Sin[x], {x, 0, 1}], $PlotLegendData[___]]", "{}", 0);
    /* "Expressions": label derived from the function itself. */
    assert_eval_eq(
        "Cases[Plot[Sin[x], {x, 0, 1}, PlotLegends -> \"Expressions\"], $PlotLegendData[___]]",
        "{$PlotLegendData[{RGBColor[0.2, 0.4, 0.8], \"Sin[x]\"}]}", 0);
    /* An explicit label list is used as given, paired with the palette
     * colors for a multi-curve plot. */
    assert_eval_eq(
        "Cases[Plot[{Sin[x], Cos[x]}, {x, 0, 1}, PlotLegends -> {\"s\", \"c\"}], $PlotLegendData[___]]",
        "{$PlotLegendData[{RGBColor[0.368417, 0.506779, 0.709798], \"s\"},"
        " {RGBColor[0.880722, 0.611041, 0.142051], \"c\"}]}", 0);
}

void test_region_function_and_exclusions_split_domain(void) {
    /* RegionFunction excludes the middle band -- two disjoint runs. */
    assert_eval_eq(
        "Length[Plot[Sin[x], {x, -3, 3}, RegionFunction -> (Abs[#1] > 1 &)][[1]]]",
        "2", 0);
    /* Exclusions forces a break at x=0 even though Sin is perfectly smooth
     * there. */
    assert_eval_eq(
        "Length[Plot[Sin[x], {x, -3, 3}, Exclusions -> {0}][[1]]]",
        "2", 0);
    /* Exclusions -> {x == a} (an equation) is accepted too. */
    assert_eval_eq(
        "Length[Plot[Sin[x], {x, -3, 3}, Exclusions -> {x == 0}][[1]]]",
        "2", 0);
}

void test_label_style_passthrough(void) {
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, LabelStyle -> Red], (LabelStyle -> v_) -> v]]",
        "RGBColor[1, 0, 0]", 0);
    /* An explicit AxesStyle still passes through distinctly -- the
     * fallback-vs-override resolution itself happens in render.c at draw
     * time, not in this structural passthrough. */
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, LabelStyle -> Red, AxesStyle -> Blue],"
        " (AxesStyle -> v_) -> v]]",
        "RGBColor[0, 0, 1]", 0);
}

/* ---- ListPlot: symbolic Graphics[...] construction ---- */

void test_listplot_returns_graphics_head(void) {
    assert_eval_eq("Head[ListPlot[{1, 4, 9}]]", "Graphics", 0);
    /* A non-list / all-missing argument leaves ListPlot unevaluated. */
    assert_eval_eq("ListPlot[x]", "ListPlot[x]", 0);
    assert_eval_eq("ListPlot[{}]", "ListPlot[{}]", 0);
    assert_eval_eq("ListPlot[{a, b, c}]", "ListPlot[{a, b, c}]", 0);
}

void test_listplot_heights_form(void) {
    /* {y1,...,yn} -> points {i, yi}. */
    assert_eval_eq("FullForm[ListPlot[{1, 4, 9}][[1]]]",
        "List[Point[List[List[1.0, 1.0], List[2.0, 4.0], List[3.0, 9.0]]]]", 0);
    /* Data is evaluated (ListPlot is not HoldAll), so Table works. */
    assert_eval_eq("Length[ListPlot[Table[i^2, {i, 1, 5}]][[1, 1, 1]]]", "5", 0);
    /* A non-numeric height is missing; its index slot is skipped. */
    assert_eval_eq("FullForm[ListPlot[{1, x, 3}][[1]]]",
        "List[Point[List[List[1.0, 1.0], List[3.0, 3.0]]]]", 0);
}

void test_listplot_pairs_form(void) {
    /* {{x,y},...} -> the given points (one Point primitive). */
    assert_eval_eq("FullForm[ListPlot[{{0, 0}, {1, 1}, {2, 4}}][[1]]]",
        "List[Point[List[List[0.0, 0.0], List[1.0, 1.0], List[2.0, 4.0]]]]", 0);
}

void test_listplot_datarange_maps_heights(void) {
    /* DataRange -> {xmin, xmax} spreads heights uniformly across the range. */
    assert_eval_eq("FullForm[ListPlot[{1, 4, 9}, DataRange -> {0, 1}][[1]]]",
        "List[Point[List[List[0.0, 1.0], List[0.5, 4.0], List[1.0, 9.0]]]]", 0);
}

void test_listplot_joined_emits_line(void) {
    /* Joined -> True draws a Line polyline instead of a Point cloud. */
    assert_eval_eq("Head[ListPlot[{{1, 1}, {2, 4}}, Joined -> True][[1, 1]]]",
        "Line", 0);
    assert_eval_eq("Head[ListPlot[{{1, 1}, {2, 4}}][[1, 1]]]", "Point", 0);
}

void test_listplot_multiple_datasets_get_palette(void) {
    /* A list of sublists that aren't 2-pairs is several datasets, each
     * prefixed by a distinct palette colour directive. */
    assert_eval_eq(
        "Cases[ListPlot[{{1, 2, 3}, {4, 5, 6}}][[1]], RGBColor[___]]",
        "{RGBColor[0.368417, 0.506779, 0.709798],"
        " RGBColor[0.880722, 0.611041, 0.142051]}", 0);
    /* DataRange -> All forces a flat list of pairs to be read as datasets. */
    assert_eval_eq(
        "Length[Cases[ListPlot[{{1, 1}, {2, 4}}, DataRange -> All][[1]], Point[___]]]",
        "2", 0);
}

void test_listplot_filling_builds_stems(void) {
    /* Filling -> Axis draws one vertical stem Line per point to y = 0,
     * wrapped in Opacity, then restores the curve colour before the points. */
    assert_eval_eq(
        "Length[Cases[ListPlot[{1, 4, 9}, Filling -> Axis][[1]], Line[___]]]",
        "3", 0);
    assert_eval_eq(
        "MemberQ[ListPlot[{1, 4, 9}, Filling -> Axis][[1]], Opacity[0.3]]",
        "True", 0);
}

void test_listplot_default_options_injected(void) {
    /* ListPlot injects Axes -> True and AspectRatio -> 1/GoldenRatio. */
    assert_eval_eq("First[Cases[ListPlot[{1, 2, 3}], (Axes -> v_) -> v]]", "True", 0);
    assert_eval_eq(
        "Abs[First[Cases[ListPlot[{1, 2, 3}], (AspectRatio -> v_) -> v]]"
        " - N[1/GoldenRatio]] < 10^-6", "True", 0);
    /* An explicit PlotStyle is passed through and suppresses the default. */
    assert_eval_eq(
        "First[Cases[ListPlot[{1, 2}, PlotStyle -> Red], (PlotStyle -> v_) -> v]]",
        "RGBColor[1, 0, 0]", 0);
}

void test_listplot_options_registered(void) {
    assert_eval_eq("Options[ListPlot, Joined]", "{Joined -> False}", 0);
    assert_eval_eq("Options[ListPlot, DataRange]", "{DataRange -> Automatic}", 0);
}

void test_listplot_legends_metadata(void) {
    assert_eval_eq("Cases[ListPlot[{1, 2, 3}], $PlotLegendData[___]]", "{}", 0);
    assert_eval_eq(
        "Cases[ListPlot[{{1, 2, 3}, {4, 5, 6}}, PlotLegends -> {\"a\", \"b\"}],"
        " $PlotLegendData[___]]",
        "{$PlotLegendData[{RGBColor[0.368417, 0.506779, 0.709798], \"a\"},"
        " {RGBColor[0.880722, 0.611041, 0.142051], \"b\"}]}", 0);
}

#ifdef USE_GRAPHICS
/* The frame minor-tick subdivision policy frame_minor_divs() lives in
 * render.c; exercise it directly. A "nice" step (1/2/5 x 10^k) chooses round
 * minor spacings: 1 -> fifths, 2 -> quarters, 5 -> fifths, across magnitudes. */
void test_frame_minor_divs_policy(void) {
    ASSERT(frame_minor_divs(1.0)   == 5);
    ASSERT(frame_minor_divs(2.0)   == 4);
    ASSERT(frame_minor_divs(5.0)   == 5);
    ASSERT(frame_minor_divs(10.0)  == 5);
    ASSERT(frame_minor_divs(20.0)  == 4);
    ASSERT(frame_minor_divs(50.0)  == 5);
    ASSERT(frame_minor_divs(0.1)   == 5);
    ASSERT(frame_minor_divs(0.2)   == 4);
    ASSERT(frame_minor_divs(0.5)   == 5);
    ASSERT(frame_minor_divs(0.0)   == 5); /* degenerate guard */
}

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

/* Regression test for a real bug: Polygon[] silently rendered nothing for
 * a clockwise vertex list (e.g. {{0,0},{0,1},{1,1},{1,0}}, the natural
 * reading order for a square's corners) because raylib's DrawTriangleFan
 * requires counter-clockwise winding, while Mathematica's Polygon[]
 * imposes no winding convention on the caller. polygon_signed_area's sign
 * is what the renderer checks to decide whether to reverse the vertex
 * list before drawing -- confirmed by an actual screenshot during
 * development (a un-reversed clockwise square rendered as a blank
 * window), not just reasoned about. */
void test_polygon_signed_area_winding_detection(void) {
    /* Coordinates here are post-y-negation draw space (what render.c
     * actually feeds polygon_signed_area), matching the exact failing case:
     * Polygon[{{0,0},{0,1},{1,1},{1,0}}] -- the natural reading order for a
     * unit square's corners -- becomes (0,0),(0,-1),(1,-1),(1,0) in draw
     * space, which has positive signed area (the renderer reverses it). */
    double cw_x[4] = { 0, 0, 1, 1 };
    double cw_y[4] = { 0, -1, -1, 0 };
    ASSERT(polygon_signed_area(cw_x, cw_y, 4) > 0.0);

    /* The reverse vertex ordering of the same square: negative area --
     * already correctly wound, the renderer leaves it alone. */
    double ccw_x[4] = { 1, 1, 0, 0 };
    double ccw_y[4] = { 0, -1, -1, 0 };
    ASSERT(polygon_signed_area(ccw_x, ccw_y, 4) < 0.0);

    /* A degenerate (zero-area) "polygon" -- e.g. all points collinear --
     * is neither winding; must not crash or loop. */
    double line_x[3] = { 0, 1, 2 };
    double line_y[3] = { 0, 0, 0 };
    ASSERT(polygon_signed_area(line_x, line_y, 3) == 0.0);
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
    TEST(test_graphics_options_registered);
    TEST(test_show_merges_options);
    TEST(test_show_merges_frame_option);
    TEST(test_plot_frame_passthrough);
    TEST(test_plot_frame_suppresses_axes_default);
    TEST(test_plot_aspect_ratio_default);
    TEST(test_plot_aspect_ratio_explicit_number);
    TEST(test_plot_aspect_ratio_symbolic_constant);
    TEST(test_plot_aspect_ratio_symbol_settings);
    TEST(test_plot_image_size_passthrough);
    TEST(test_axes_origin_passthrough);
    TEST(test_axes_style_and_ticks_style_passthrough);
    TEST(test_frame_label_and_rotate_label_passthrough);
    TEST(test_plot_range_padding_passthrough);
    TEST(test_grid_lines_passthrough);
    TEST(test_prolog_epilog_passthrough);
    TEST(test_named_color_constants_resolve);
    TEST(test_hue_color_directive);
    TEST(test_color_function_builds_per_segment_colors);
    TEST(test_filling_builds_polygon);
    TEST(test_filling_splits_at_baseline_crossing);
    TEST(test_plot_legends_metadata);
    TEST(test_region_function_and_exclusions_split_domain);
    TEST(test_label_style_passthrough);
    TEST(test_listplot_returns_graphics_head);
    TEST(test_listplot_heights_form);
    TEST(test_listplot_pairs_form);
    TEST(test_listplot_datarange_maps_heights);
    TEST(test_listplot_joined_emits_line);
    TEST(test_listplot_multiple_datasets_get_palette);
    TEST(test_listplot_filling_builds_stems);
    TEST(test_listplot_default_options_injected);
    TEST(test_listplot_options_registered);
    TEST(test_listplot_legends_metadata);
#ifdef USE_GRAPHICS
    TEST(test_frame_minor_divs_policy);
    TEST(test_window_height_policy);
    TEST(test_polygon_signed_area_winding_detection);
#endif

    printf("All graphics tests passed!\n");
    symtab_clear();
    return 0;
}
