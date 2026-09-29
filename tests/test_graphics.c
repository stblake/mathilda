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
#include "plot_common.h"     /* gfx_*: the style sizes both back ends share */
#ifdef USE_GRAPHICS
#include "render.h"
#include "render_common.h"   /* frame_minor_divs: shared tick-spacing policy */
#include "graphics_export.h" /* graphics_render_rgba: offscreen raster */
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

void test_export_graphics_pdf(void) {
    /* PDF export is the dependency-free, headless vector path, so it runs in
     * the suite (PNG/JPEG need a real GL context and are not exercised here).
     * Export returns the path on success and $Failed on a claimed-but-failed
     * write, so the returned path is itself proof the file was written. */
    assert_eval_eq(
        "StringQ[Export[\"/tmp/mathilda_test_gx.pdf\", Plot[Sin[x], {x, 0, 6}]]]",
        "True", 0);
    assert_eval_eq("FileExistsQ[\"/tmp/mathilda_test_gx.pdf\"]", "True", 0);
    assert_eval_eq(
        "FileExistsQ[Export[\"/tmp/mathilda_test_gx2.pdf\", "
        "Graphics[{Line[{{0, 0}, {1, 1}, {2, 0}}]}]]]",
        "True", 0);
    /* Graphics3D export is not supported yet -> $Failed. */
    assert_eval_eq(
        "Export[\"/tmp/mathilda_test_gx3.pdf\", Graphics3D[{Point[{0, 0, 0}]}]]",
        "$Failed", 0);
}

/* ---- Show[g1, g2, ...] and the export option/directive fixes ---------- */

/* Read a whole file into a NUL-terminated buffer (caller frees). The PDF
 * writer emits an uncompressed ASCII content stream, so its drawing
 * operators can be asserted on directly. */
static char* slurp_file(const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) return NULL;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    char* buf = malloc((size_t)n + 1);
    size_t got = fread(buf, 1, (size_t)n, f);
    buf[got] = '\0';
    fclose(f);
    return buf;
}

static size_t count_substr(const char* hay, const char* needle) {
    size_t c = 0, n = strlen(needle);
    for (const char* p = strstr(hay, needle); p; p = strstr(p + n, needle)) c++;
    return c;
}

/* Export `expr_src` (a Graphics expression) to `path` and return the PDF
 * text (caller frees). */
static char* export_pdf_text(const char* expr_src, const char* path) {
    char cmd[2048];
    snprintf(cmd, sizeof(cmd), "Export[\"%s\", %s]", path, expr_src);
    struct Expr* parsed = parse_expression(cmd);
    ASSERT(parsed != NULL);
    Expr* r = evaluate(parsed);
    expr_free(parsed);
    ASSERT(r && r->type == EXPR_STRING);
    expr_free(r);
    char* text = slurp_file(path);
    ASSERT(text != NULL);
    return text;
}

/* Bug: Show[g1, g2] stayed unevaluated -- Show accepted only one graphic. */
void test_show_combines_graphics(void) {
    assert_eval_eq("Show[Graphics[Line[{{0,0},{1,1}}]], Graphics[Point[{2,2}]]]",
                   "-Graphics-", 0);
    /* One List scope per input, so directives cannot leak between them. */
    assert_eval_eq("Show[Graphics[{Red, Line[{{0,0},{1,1}}]}], Graphics[Point[{2,2}]]][[1]]",
                   "{{RGBColor[1, 0, 0], Line[{{0, 0}, {1, 1}}]}, {Point[{2, 2}]}}", 0);
    /* The list form, nested lists included, is the same call. */
    assert_eval_eq("Show[{Graphics[Point[{0,0}]], {Graphics[Point[{1,1}]]}}] === "
                   "Show[Graphics[Point[{0,0}]], Graphics[Point[{1,1}]]]", "True", 0);
    assert_eval_eq("Show[{Graphics[Point[{0,0}]]}] === Graphics[Point[{0,0}]]", "True", 0);
    /* Options come from the first graphic... */
    assert_eval_eq("Cases[Show[Graphics[{}, Axes -> True], Graphics[{}, Frame -> True]],"
                   " (h:Axes|Frame -> v_) :> {h, v}]", "{{Axes, True}}", 0);
    /* ...unless Show overrides them. */
    assert_eval_eq("Cases[Show[Graphics[{}, Axes -> True], Graphics[{}], Axes -> False],"
                   " (Axes -> v_) :> v]", "{False}", 0);
    /* Two Plot outputs combine too (the reported case). */
    assert_eval_eq("Head[Show[Plot[Sin[x], {x, 0, 1}], Plot[Cos[x], {x, 0, 1}]]]",
                   "Graphics", 0);
    /* 2D and 3D do not mix; non-graphics arguments are rejected. */
    assert_eval_eq("Head[Show[Graphics[Point[{0,0}]], Graphics3D[Point[{0,0,0}]]]]",
                   "Show", 0);
    assert_eval_eq("Head[Show[Graphics[Point[{0,0}]], 5]]", "Show", 0);
    assert_eval_eq("Head[Show[Axes -> True]]", "Show", 0);
}

void test_show_merges_plot_range(void) {
    /* Union of an explicit range with another input's extent. */
    assert_eval_eq("Cases[Show[Graphics[Line[{{0,0},{1,1}}], PlotRange -> {{0,1},{0,1}}],"
                   " Graphics[Point[{2,3}]]], (PlotRange -> v_) :> v]",
                   "{{{0.0, 2.0}, {0.0, 3.0}}}", 0);
    /* Circle[] (unit circle by default) contributes its extent. */
    assert_eval_eq("Cases[Show[Graphics[Circle[]], Graphics[Disk[{3, 0}],"
                   " PlotRange -> {{2, 4}, {-1, 1}}]], (PlotRange -> v_) :> v]",
                   "{{{-1.0, 4.0}, {-1.0, 1.0}}}", 0);
    /* No input fixes a range: stays automatic (no PlotRange option). */
    assert_eval_eq("Cases[Show[Graphics[Point[{0,0}]], Graphics[Point[{1,1}]]],"
                   " (PlotRange -> _)]", "{}", 0);
    /* An explicit PlotRange in Show wins. */
    assert_eval_eq("Cases[Show[Graphics[Point[{0,0}]], Graphics[Point[{1,1}]],"
                   " PlotRange -> {{0, 5}, {0, 5}}], (PlotRange -> v_) :> v]",
                   "{{{0, 5}, {0, 5}}}", 0);
}

void test_show_keeps_each_plot_style(void) {
    /* A single-curve Plot is drawn in its PlotStyle option; Show keeps only
     * g1's options, so each input's style is baked into its own scope. */
    assert_eval_eq("g = Show[Plot[x, {x, 0, 1}, PlotStyle -> Red, PlotPoints -> 2,"
                   " MaxRecursion -> 0], Plot[x^2, {x, 0, 1}, PlotStyle -> Blue,"
                   " PlotPoints -> 2, MaxRecursion -> 0]]; {g[[1, 1, 1]], g[[1, 2, 1]]}",
                   "{RGBColor[1, 0, 0], RGBColor[0, 0, 1]}", 0);
    assert_eval_eq("Cases[g, (PlotStyle -> _)]", "{}", 0);
    /* $PlotResample would re-sample only g1's curves; it is dropped. */
    assert_eval_eq("FreeQ[g, $PlotResample]", "True", 0);
    /* Legends are concatenated. */
    assert_eval_eq("List @@ Select[List @@ Show[Plot[x, {x, 0, 1}, PlotLegends -> {\"a\"}],"
                   " Plot[x^2, {x, 0, 1}, PlotLegends -> {\"b\"}]],"
                   " Head[#] === $PlotLegendData &][[1, All, 2]]", "{\"a\", \"b\"}", 0);
}

void test_style_directive_symbols(void) {
    assert_eval_eq("Dashed", "Dashing[{Small, Small}]", 0);
    assert_eval_eq("Dotted", "Dashing[{0, Small}]", 0);
    assert_eval_eq("DotDashed", "Dashing[{0, Small, Small, Small}]", 0);
    assert_eval_eq("Thick", "Thickness[Large]", 0);
    assert_eval_eq("Thin", "Thickness[Tiny]", 0);
    assert_eval_eq("Directive[Red, Dashed]", "Directive[RGBColor[1, 0, 0], Dashing[{Small, Small}]]", 0);

    /* The shared size resolution both back ends use (plot width 500). */
    double v, d[GFX_MAX_DASH];
    int nd;
    Expr* e = evaluate(parse_expression("PointSize[0.1]"));
    ASSERT(gfx_point_radius_pts(e, 500.0, &v) && fabs(v - 25.0) < 1e-9);  /* diameter 50 */
    expr_free(e);
    e = evaluate(parse_expression("Thickness[0.01]"));
    ASSERT(gfx_thickness_pts(e, 500.0, &v) && fabs(v - 5.0) < 1e-9);
    expr_free(e);
    e = evaluate(parse_expression("AbsoluteThickness[2]"));
    ASSERT(gfx_thickness_pts(e, 500.0, &v) && fabs(v - 2.0) < 1e-9);
    expr_free(e);
    e = evaluate(parse_expression("Dashing[{1/100, 1/50}]"));
    ASSERT(gfx_dash_pts(e, 500.0, d, GFX_MAX_DASH, &nd) && nd == 2
           && fabs(d[0] - 5.0) < 1e-9 && fabs(d[1] - 10.0) < 1e-9);
    expr_free(e);
    e = evaluate(parse_expression("Dashing[{}]"));
    ASSERT(gfx_dash_pts(e, 500.0, d, GFX_MAX_DASH, &nd) && nd == 0);
    expr_free(e);
}

/* PlotStyle -> {s1, s2, ...} was ignored for multi-curve plots (palette
 * colours only); each style is now baked into its curve, non-colour styles in
 * their own List scope, and legends show the styled colours. */
void test_plot_style_per_curve(void) {
    assert_eval_eq("g = Plot[{x, 2 x}, {x, 0, 1}, PlotPoints -> 2, MaxRecursion -> 0,"
                   " PlotStyle -> {Red, Directive[Blue, Dashed]}]; g[[1, 1]]",
                   "RGBColor[1, 0, 0]", 0);
    assert_eval_eq("g[[1, 3, 1]]", "Directive[RGBColor[0, 0, 1], Dashing[{Small, Small}]]", 0);
    assert_eval_eq("Head[g[[1, 3, 2]]]", "Line", 0);
    /* A style without a colour keeps the curve's palette colour. */
    assert_eval_eq("Head[Plot[{x, 2 x}, {x, 0, 1}, PlotStyle -> {Thick, Dashed}][[1, 1, 1, 1]]]",
                   "RGBColor", 0);
    assert_eval_eq("List @@ Plot[{x, 2 x}, {x, 0, 1}, PlotStyle -> {Red, Blue},"
                   " PlotLegends -> {\"a\", \"b\"}][[-1, All, 1]]",
                   "{RGBColor[1, 0, 0], RGBColor[0, 0, 1]}", 0);
    assert_eval_eq("ListPlot[{{1, 2, 3}, {3, 4, 5}}, PlotStyle -> {Red, Blue}][[1, {1, 3}]]",
                   "{RGBColor[1, 0, 0], RGBColor[0, 0, 1]}", 0);
    assert_eval_eq("ListPlot[{{1, 2, 3}, {3, 4, 5}}, PlotStyle -> {Red, PointSize[Large]}][[1, 3, 1]]",
                   "Directive[RGBColor[0.880722, 0.611041, 0.142051], PointSize[Large]]", 0);
    assert_eval_eq("ParametricPlot[{{t, t}, {t, 2 t}}, {t, 0, 1},"
                   " PlotStyle -> {Green, Red}][[1, 1]]", "RGBColor[0, 1, 0]", 0);
}

/* Bug: Line[{{0,0},{1/2,1}}] was dropped from the PDF (the writer read only
 * Integer/Real coordinates). Every primitive now reads exact values too. */
void test_export_pdf_exact_coordinates(void) {
    char* t = export_pdf_text("Graphics[Line[{{0, 0}, {1/2, 1}}]]", "/tmp/mathilda_test_rat.pdf");
    ASSERT(count_substr(t, " m\n") == 1 && count_substr(t, " l\n") >= 1);
    free(t);
    /* Point/Disk/Circle/Polygon/Arrow/Text with Rational and Pi/Sqrt coordinates. */
    t = export_pdf_text("Graphics[{Point[{1/3, 2/3}], Disk[{Pi/4, 1/4}, 1/10],"
                        " Circle[{3/4, 3/4}, 1/8], Polygon[{{0, 1/2}, {1/4, 1/2}, {1/8, 3/4}}],"
                        " Arrow[{{0, 0}, {Sqrt[2]/2, 1/2}}], Text[\"hi\", {Pi/4, 1/2}]}]",
                        "/tmp/mathilda_test_rat2.pdf");
    ASSERT(count_substr(t, " c\n") >= 12);          /* point + disk + circle: 3 x 4 Beziers */
    ASSERT(count_substr(t, "h f\n") >= 2);          /* polygon + arrowhead */
    ASSERT(strstr(t, "(hi) Tj") != NULL);
    free(t);
    /* Circle[] is the unit circle. */
    t = export_pdf_text("Graphics[Circle[]]", "/tmp/mathilda_test_circ.pdf");
    ASSERT(count_substr(t, " c\n") == 4);
    free(t);
}

/* Bug: the PDF ignored PlotStyle, Epilog/Prolog, AxesLabel, PlotLabel,
 * PlotLegends and Dashing. */
void test_export_pdf_options_and_directives(void) {
    char* t = export_pdf_text("Plot[x, {x, 0, 1}, PlotStyle -> Red]", "/tmp/mathilda_test_ps.pdf");
    ASSERT(strstr(t, "1.0000 0.0000 0.0000 RG") != NULL);
    free(t);
    t = export_pdf_text("Plot[x, {x, 0, 1}, PlotStyle -> Directive[Dashed, AbsoluteThickness[3]]]",
                        "/tmp/mathilda_test_dash.pdf");
    ASSERT(strstr(t, "[4.000 4.000] 0 d") != NULL);
    ASSERT(strstr(t, "3.000 w") != NULL);
    free(t);
    t = export_pdf_text("Graphics[{Dotted, Line[{{0, 0}, {1, 1}}]}]", "/tmp/mathilda_test_dot.pdf");
    ASSERT(strstr(t, "[0.000 4.000] 0 d 1 J") != NULL);   /* dots need round caps */
    free(t);
    t = export_pdf_text("Graphics[{Line[{{0, 0}, {1, 1}}]}, Epilog -> {Green, Line[{{0, 1}, {1, 0}}]},"
                        " Prolog -> {Blue, Rectangle[{0, 0}, {1, 1}]}]", "/tmp/mathilda_test_epi.pdf");
    ASSERT(strstr(t, "0.0000 1.0000 0.0000 RG") != NULL);
    ASSERT(strstr(t, "0.0000 0.0000 1.0000 rg") != NULL);
    /* Prolog draws beneath the primitives, Epilog over them. */
    ASSERT(strstr(t, "0.0000 0.0000 1.0000 rg") < strstr(t, "0.0000 1.0000 0.0000 RG"));
    free(t);
    t = export_pdf_text("Plot[{x, x^2}, {x, 0, 1}, AxesLabel -> {\"xlab\", \"ylab\"},"
                        " PlotLabel -> \"Title\", PlotLegends -> {\"alpha\", \"beta\"}]",
                        "/tmp/mathilda_test_lbl.pdf");
    ASSERT(strstr(t, "(xlab) Tj") && strstr(t, "(ylab) Tj") && strstr(t, "(Title) Tj"));
    ASSERT(strstr(t, "(alpha) Tj") && strstr(t, "(beta) Tj"));
    free(t);
    /* A directive set inside a List does not leak past it. */
    t = export_pdf_text("Graphics[{{Dashed, Line[{{0, 0}, {1, 1}}]}, Line[{{0, 1}, {1, 0}}]}]",
                        "/tmp/mathilda_test_scope.pdf");
    ASSERT(strstr(t, "[] 0 d") != NULL);
    free(t);
    /* Show[g1, g2] exports both inputs, each in its own style. */
    t = export_pdf_text("Show[Plot[x, {x, 0, 1}, PlotStyle -> Red], Plot[1 - x, {x, 0, 1}, PlotStyle -> Green]]",
                        "/tmp/mathilda_test_show.pdf");
    ASSERT(strstr(t, "1.0000 0.0000 0.0000 RG") && strstr(t, "0.0000 1.0000 0.0000 RG"));
    free(t);
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
    /* Black/White are grey levels, not RGB (matching WL's InputForm). */
    assert_eval_eq("FullForm[Black]", "GrayLevel[0]", 0);
    assert_eval_eq("FullForm[White]", "GrayLevel[1]", 0);
    /* The light variants resolve to their RGBColor[...] literals too. */
    assert_eval_eq("FullForm[LightRed]", "RGBColor[1, 0.85, 0.85]", 0);
    assert_eval_eq("FullForm[LightPurple]", "RGBColor[0.94, 0.88, 0.94]", 0);
    /* ?Red inspects the symbol, not its value (Information is HoldFirst). */
    assert_eval_eq("StringTake[Information[Red], 3]", "\"Red\"", 0);
    /* Docstrings are now sourced from info.c. */
    assert_eval_eq("StringTake[Information[LightOrange], 11]", "\"LightOrange\"", 0);
}

void test_hue_color_directive(void) {
    assert_eval_eq("Graphics[{Hue[0.5], Point[{0,0}]}]", "-Graphics-", 0);
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, PlotStyle -> Hue[0.5]], (PlotStyle -> v_) -> v]]",
        "Hue[0.5]", 0);
}

void test_cmykcolor_directive(void) {
    /* CMYKColor is an inert, protected style head: it stays unevaluated and
     * carries the Protected attribute, with its docstring set in info.c. */
    assert_eval_eq("Attributes[CMYKColor]", "{Protected}", 0);
    assert_eval_eq("FullForm[CMYKColor[0.1, 0.2, 0.3, 0.4]]",
                    "CMYKColor[0.1, 0.2, 0.3, 0.4]", 0);
    assert_eval_eq("StringTake[Information[CMYKColor], 9]", "\"CMYKColor\"", 0);
    /* All its surface forms are accepted inside a Graphics[] without error. */
    assert_eval_eq("Head[Graphics[{CMYKColor[0, 1, 1, 0], Disk[]}]]", "Graphics", 0);
    assert_eval_eq("Head[Graphics[{CMYKColor[0, 1, 1], Disk[]}]]", "Graphics", 0);
    assert_eval_eq("Head[Graphics[{CMYKColor[0, 1, 1, 0, 0.5], Disk[]}]]", "Graphics", 0);
    assert_eval_eq("Head[Graphics[{CMYKColor[{0, 1, 1, 0}], Disk[]}]]", "Graphics", 0);
    assert_eval_eq("Head[Graphics[{CMYKColor[{0, 1, 1, 0, 0.5}], Disk[]}]]", "Graphics", 0);
    /* As a Plot directive it passes through verbatim (no auto RGB conversion). */
    assert_eval_eq(
        "First[Cases[Plot[Sin[x], {x, 0, 1}, PlotStyle -> CMYKColor[1, 0, 0, 0]],"
        " (PlotStyle -> v_) -> v]]",
        "CMYKColor[1, 0, 0, 0]", 0);
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

/* ColorFunctionScaling -> False must hand the colour function the RAW value.
 *
 * The clamp that belongs to the scaling branch used to be applied to both, so
 * every cell below 0 (or above 1) got the colour of the nearest endpoint. A
 * signed field then renders as a uniform block -- which looks exactly like a
 * correct plot of a field that really is constant, so nothing about the picture
 * gives the mistake away. */
void test_densityplot_unscaled_color_function_gets_raw_value(void) {
    assert_eval_eq(
        "Union[Cases[DensityPlot[-0.5, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 2,"
        " ColorFunctionScaling -> False,"
        " ColorFunction -> Function[z, RGBColor[0., 0., -z]]],"
        " RGBColor[__], Infinity]]",
        "{RGBColor[0.0, 0.0, 0.5]}", 0);
    /* Above 1 is raw too, not pinned to the top of the ramp. */
    assert_eval_eq(
        "Union[Cases[DensityPlot[3.0, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 2,"
        " ColorFunctionScaling -> False,"
        " ColorFunction -> Function[z, GrayLevel[z/4]]],"
        " GrayLevel[_], Infinity]]",
        "{GrayLevel[0.75]}", 0);
    /* Default (scaling on) still normalises: a constant field maps to t = 0. */
    assert_eval_eq(
        "Union[Cases[DensityPlot[-0.5, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 2,"
        " ColorFunction -> Function[z, RGBColor[0., 0., z]]],"
        " RGBColor[__], Infinity]]",
        "{RGBColor[0.0, 0.0, 0.0]}", 0);
}

/* ArrayPlot: default ColorFunction is the shared "Greyscale" ramp (white at
 * the minimum, black at the maximum) -- Mathematica's own ArrayPlot default,
 * reusing named_color_ramp rather than a bespoke palette. */
void test_arrayplot_returns_graphics_head(void) {
    assert_eval_eq("Head[ArrayPlot[{{1,2},{3,4}}]]", "Graphics", 0);
}

void test_arrayplot_default_greyscale_color(void) {
    assert_eval_eq(
        "Union[Cases[ArrayPlot[{{1,2},{3,4}}], GrayLevel[__], Infinity]]",
        "{GrayLevel[0.0], GrayLevel[0.333333], GrayLevel[0.666667], GrayLevel[1.0]}", 0);
}

/* A matrix whose entries are already colour literals paints each cell that
 * colour directly, in row-major order, instead of deriving one from
 * ColorFunction -- lets ArrayPlot double as a raw pixel-grid renderer. */
void test_arrayplot_color_matrix_passthrough(void) {
    assert_eval_eq(
        "Cases[ArrayPlot[{{Red, Blue}, {Blue, Red}}], RGBColor[__], Infinity]",
        "{RGBColor[1, 0, 0], RGBColor[0, 0, 1], RGBColor[0, 0, 1], RGBColor[1, 0, 0]}", 0);
}

/* Numeric and colour cells freely mix within the same array: a cell that is
 * already a colour literal paints directly (bypassing ColorFunction) while
 * the rest of the array still follows the normal greyscale heatmap. */
void test_arrayplot_mixed_numeric_and_color_cells(void) {
    assert_eval_eq(
        "Union[Cases[ArrayPlot[{{1, 0, Pink}, {0, 1, Red}}], RGBColor[__], Infinity]]",
        "{RGBColor[1, 0, 0], RGBColor[1, 0.5, 0.5]}", 0);
    assert_eval_eq(
        "Union[Cases[ArrayPlot[{{1, 0, Pink}, {0, 1, Red}}], GrayLevel[__], Infinity]]",
        "{GrayLevel[0.0], GrayLevel[1.0]}", 0);
}

void test_arrayplot_invalid_args_stay_unevaluated(void) {
    assert_eval_eq("ArrayPlot[foo]", "ArrayPlot[foo]", 0);
    /* A cell that is neither numeric nor a colour literal leaves the whole
     * call unevaluated, matching the "can't evaluate this" convention. */
    assert_eval_eq("ArrayPlot[{{1, 2, foo}}]", "ArrayPlot[{{1, 2, foo}}]", 0);
}

/* Mesh -> All draws (rows+1) horizontal + (cols+1) vertical grid Lines. */
void test_arrayplot_mesh_option(void) {
    assert_eval_eq(
        "Length[Cases[ArrayPlot[{{1,2},{3,4}}, Mesh -> All], Line[__], Infinity]]",
        "6", 0);
}

/* AspectRatio defaults to rows/cols so cells render square, unless the
 * caller overrides it. */
void test_arrayplot_aspect_ratio_matches_dimensions(void) {
    assert_eval_eq(
        "Cases[ArrayPlot[{{1,2,3},{4,5,6}}], Rule[AspectRatio, r_] :> r]",
        "{0.666667}", 0);
}

/* ColorRules gives an explicit colour to cells matching a named value,
 * checked before ColorFunction; cells matching no rule still fall through
 * to the normal scaled ColorFunction colour (not left grey/uncoloured). */
void test_arrayplot_color_rules_override_matching_cells(void) {
    assert_eval_eq(
        "Union[Cases[ArrayPlot[{{1,0},{0,1}}, ColorRules -> {1 -> Pink, 0 -> Yellow}],"
        " RGBColor[__], Infinity]]",
        "{RGBColor[1, 0.5, 0.5], RGBColor[1, 1, 0]}", 0);
    /* A value with no matching rule keeps the default greyscale mapping. */
    assert_eval_eq(
        "Union[Cases[ArrayPlot[{{1,0,0.5},{0,1,0.5}}, ColorRules -> {1 -> Pink, 0 -> Yellow}],"
        " GrayLevel[__], Infinity]]",
        "{GrayLevel[0.5]}", 0);
    /* ColorRules -> None is equivalent to omitting the option. */
    assert_eval_eq(
        "ArrayPlot[{{1,2},{3,4}}, ColorRules -> None] === ArrayPlot[{{1,2},{3,4}}]",
        "True", 0);
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

void test_listplot_joined_filling_builds_polygons(void) {
    /* Joined -> True with Filling fills the continuous region under the
     * connecting curve with Polygon[] quads (one per segment), not isolated
     * per-point stem Lines -- so there are no fill Lines, only the single
     * joined curve Line, and at least one Polygon. */
    assert_eval_eq(
        "Length[Cases[ListPlot[{0, 1, 2, 3}, Joined -> True, Filling -> Axis][[1]],"
        " Polygon[___]]] >= 1", "True", 0);
    /* Three points -> one connecting Line, and no stem Lines were emitted. */
    assert_eval_eq(
        "Length[Cases[ListPlot[{0, 1, 2, 3}, Joined -> True, Filling -> Axis][[1]],"
        " Line[___]]]", "1", 0);
    /* A segment crossing the baseline (y: 1 -> -1) splits into two triangles. */
    assert_eval_eq(
        "Length[Cases[ListPlot[{1, -1}, Joined -> True, Filling -> Axis][[1]],"
        " Polygon[pts_] :> Length[pts]]]", "2", 0);
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

/* Regression test for a real bug (reported via screenshot): ArrayPlot's
 * default AspectRatio -> rows/cols, combined with Frame's fixed-pixel
 * margins (not proportional between width and height), made the coloured
 * grid letterbox inside its own frame instead of filling it edge-to-edge.
 * gfx_window_height alone sizes the *window* to match the aspect ratio;
 * gfx_window_height_fit_region additionally corrects for the margins so
 * the *region* the data actually renders into hits it. */
void test_window_height_fit_region_no_letterbox(void) {
    const double dw = 20.0, dh = 8.0;   /* an 8-row x 20-col ArrayPlot */
    double a = dh / dw;                 /* AspectRatio -> rows/cols = 0.4 */
    long width = 640;

    long naive_h = gfx_window_height(width, 480, a, false, false, dw, dh);
    long fit_h   = gfx_window_height_fit_region(width, naive_h, a, false, false,
                                                dw, dh, /*frame=*/true, /*axes=*/false,
                                                /*frame_label=*/false);

    float mL, mR, mT, mB;
    gfx_horizontal_margins((float)width, true, false, false, &mL, &mR);
    double reg_w = width - mL - mR;

    /* The corrected height's region matches the target ratio closely. */
    gfx_vertical_margins((float)fit_h, true, false, false, &mT, &mB);
    double reg_h = fit_h - mT - mB;
    ASSERT(fabs(reg_h / reg_w - a) < 0.01);

    /* The naive (pre-fix) height's region does not -- confirms this case
     * actually exercises the bug rather than one margins already got right. */
    gfx_vertical_margins((float)naive_h, true, false, false, &mT, &mB);
    double naive_reg_h = naive_h - mT - mB;
    ASSERT(fabs(naive_reg_h / reg_w - a) > 0.05);

    /* height_pinned / aspect_full: no-op, matching gfx_window_height. */
    ASSERT(gfx_window_height_fit_region(width, 480, a, true, false, dw, dh, true, false, false) == 480);
    ASSERT(gfx_window_height_fit_region(width, 480, a, false, true, dw, dh, true, false, false) == 480);

    /* Neither Frame nor Axes: nothing to correct for. */
    ASSERT(gfx_window_height_fit_region(width, naive_h, a, false, false, dw, dh, false, false, false) == naive_h);
}

/* Regression test for a real bug (reported via screenshot): raster export
 * (Export["f.png", ArrayPlot[...]]/DensityPlot) sized the offscreen canvas to
 * a fixed 4:3 shape, ignoring the plot's AspectRatio, so graphics_render_in_region
 * letterboxed every square/aspect-driven plot with wide left/right padding
 * while the top and bottom filled -- the "horizontal padding far too large"
 * report. graphics_raster_dims now sizes the canvas from AspectRatio exactly as
 * the on-screen window does, so the margin-reduced DATA REGION matches the data
 * aspect and the plot fills its frame edge-to-edge. */
void test_raster_dims_no_letterbox(void) {
    /* A square 3x3 ArrayPlot: Frame -> True, Axes -> False, AspectRatio -> 1
     * (rows/cols), PlotRange -> {{0,3},{0,3}}. */
    struct Expr* parsed = parse_expression("ArrayPlot[{{1,2,3},{4,5,6},{7,8,9}}]");
    ASSERT(parsed != NULL);
    Expr* g = evaluate(parsed);
    expr_free(parsed);
    ASSERT(g != NULL && g->type == EXPR_FUNCTION);

    int w = 0, h = 0;
    graphics_raster_dims(g, &w, &h);
    ASSERT(w > 0 && h > 0);

    /* The margin-reduced region (same Frame margins the renderer uses) must
     * match the data aspect of 1 -- i.e. no horizontal (or vertical) padding. */
    float mL, mR, mT, mB;
    gfx_horizontal_margins((float)w, /*frame=*/true, /*axes=*/false, /*label=*/false, &mL, &mR);
    gfx_vertical_margins((float)h, /*frame=*/true, /*axes=*/false, /*label=*/false, &mT, &mB);
    double reg_w = (double)w - mL - mR;
    double reg_h = (double)h - mT - mB;
    ASSERT(reg_w > 0 && reg_h > 0);
    ASSERT(fabs(reg_h / reg_w - 1.0) < 0.02);

    expr_free(g);

    /* A pinned ImageSize -> {w,h} must win over AspectRatio (the user asked for
     * that exact box), matching Mathematica. */
    parsed = parse_expression(
        "ArrayPlot[{{1,2,3},{4,5,6},{7,8,9}}, ImageSize -> {400, 300}]");
    ASSERT(parsed != NULL);
    g = evaluate(parsed);
    expr_free(parsed);
    graphics_raster_dims(g, &w, &h);
    ASSERT(w == 400 && h == 300);
    expr_free(g);
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

/* The CMYK->RGB conversion cmyk_to_rgb() lives in render.c (USE_GRAPHICS only);
 * exercise the subtractive-model math and the [0,1] input clipping directly. */
void test_cmyk_to_rgb_conversion(void) {
    double r, g, b;
    /* Pure black ink (k=1) is black regardless of c/m/y. */
    cmyk_to_rgb(0, 0, 0, 1, &r, &g, &b);
    ASSERT(r == 0.0 && g == 0.0 && b == 0.0);
    /* No ink at all is white. */
    cmyk_to_rgb(0, 0, 0, 0, &r, &g, &b);
    ASSERT(r == 1.0 && g == 1.0 && b == 1.0);
    /* Full cyan -> pure red's complement: (0,1,1) in RGB. */
    cmyk_to_rgb(1, 0, 0, 0, &r, &g, &b);
    ASSERT(r == 0.0 && g == 1.0 && b == 1.0);
    /* Full magenta -> (1,0,1); full yellow -> (1,1,0). */
    cmyk_to_rgb(0, 1, 0, 0, &r, &g, &b);
    ASSERT(r == 1.0 && g == 0.0 && b == 1.0);
    cmyk_to_rgb(0, 0, 1, 0, &r, &g, &b);
    ASSERT(r == 1.0 && g == 1.0 && b == 0.0);
    /* Half black darkens an otherwise white pixel by half. */
    cmyk_to_rgb(0, 0, 0, 0.5, &r, &g, &b);
    ASSERT(r == 0.5 && g == 0.5 && b == 0.5);
    /* Out-of-range inputs clip to [0,1] (negative c -> 0, k>1 -> 1 -> black). */
    cmyk_to_rgb(-1, 0, 0, 0, &r, &g, &b);
    ASSERT(r == 1.0 && g == 1.0 && b == 1.0);
    cmyk_to_rgb(0, 0, 0, 2, &r, &g, &b);
    ASSERT(r == 0.0 && g == 0.0 && b == 0.0);
}

/* Render `src` offscreen at w x h; NULL (test skipped) unless opted in.
 * Opt-in (MATHILDA_TEST_RASTER=1) rather than automatic: Raylib's InitWindow
 * segfaults inside GLFW when the session exists but no monitor is usable
 * (display asleep or locked), which gui_session_available() cannot detect,
 * and that would take the whole suite down on such a machine. */
static unsigned char* render_src(const char* src, int w, int h, int* ow, int* oh) {
    const char* opt = getenv("MATHILDA_TEST_RASTER");
    if (!opt || !opt[0] || !gui_session_available()) return NULL;
    struct Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* g = evaluate(parsed);
    expr_free(parsed);
    unsigned char* px = graphics_render_rgba(g, w, h, ow, oh);
    expr_free(g);
    return px;
}

/* Bug: the raster renderer read PointSize[d] as a RADIUS IN PLOT
 * COORDINATES, so over a large coordinate range (0..2000 here) PointSize[0.1]
 * was a 0.02-pixel dot and vanished from PNG exports, while the PDF (reading
 * it as a fraction of the plot width) drew it. Both now use Mathematica's
 * meaning: the diameter as a fraction of the plot width. */
void test_raster_point_size_fraction_of_width(void) {
    int w = 0, h = 0;
    unsigned char* px = render_src(
        "Graphics[{Red, PointSize[0.1], Point[{1000, 1000}]},"
        " PlotRange -> {{0, 2000}, {0, 2000}}]", 200, 200, &w, &h);
    if (!px) { printf("  (skipped: set MATHILDA_TEST_RASTER=1 with a live display)\n"); return; }
    /* Diameter 0.1 x 200 px = 20 px: the centre and a point 6 px off it are red. */
    const unsigned char* c = px + ((size_t)(h / 2) * w + w / 2) * 4;
    ASSERT(c[0] > 200 && c[1] < 80 && c[2] < 80);
    c = px + ((size_t)(h / 2) * w + w / 2 + 6) * 4;
    ASSERT(c[0] > 200 && c[1] < 80 && c[2] < 80);
    free(px);
}

/* Dashing reaches the raster renderer: a dashed horizontal line alternates
 * ink and paper along its row. */
void test_raster_dashing(void) {
    int w = 0, h = 0;
    unsigned char* px = render_src(
        "Graphics[{Black, AbsoluteThickness[3], AbsoluteDashing[{10, 10}],"
        " Line[{{0, 1}, {10, 1}}]}, PlotRange -> {{0, 10}, {0, 2}}]", 300, 100, &w, &h);
    if (!px) { printf("  (skipped: set MATHILDA_TEST_RASTER=1 with a live display)\n"); return; }
    int runs = 0, prev = -1;
    for (int x = 10; x < w - 10; x++) {
        const unsigned char* c = px + ((size_t)(h / 2) * w + x) * 4;
        int ink = c[0] < 128;
        if (ink != prev) { runs++; prev = ink; }
    }
    ASSERT(runs >= 10);            /* ~14 dashes + gaps across 280 px */
    free(px);
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
    TEST(test_export_graphics_pdf);
    TEST(test_show_combines_graphics);
    TEST(test_show_merges_plot_range);
    TEST(test_show_keeps_each_plot_style);
    TEST(test_style_directive_symbols);
    TEST(test_plot_style_per_curve);
    TEST(test_export_pdf_exact_coordinates);
    TEST(test_export_pdf_options_and_directives);
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
    TEST(test_cmykcolor_directive);
    TEST(test_color_function_builds_per_segment_colors);
    TEST(test_filling_builds_polygon);
    TEST(test_filling_splits_at_baseline_crossing);
    TEST(test_plot_legends_metadata);
    TEST(test_densityplot_unscaled_color_function_gets_raw_value);
    TEST(test_arrayplot_returns_graphics_head);
    TEST(test_arrayplot_default_greyscale_color);
    TEST(test_arrayplot_color_matrix_passthrough);
    TEST(test_arrayplot_mixed_numeric_and_color_cells);
    TEST(test_arrayplot_invalid_args_stay_unevaluated);
    TEST(test_arrayplot_mesh_option);
    TEST(test_arrayplot_aspect_ratio_matches_dimensions);
    TEST(test_arrayplot_color_rules_override_matching_cells);
    TEST(test_region_function_and_exclusions_split_domain);
    TEST(test_label_style_passthrough);
    TEST(test_listplot_returns_graphics_head);
    TEST(test_listplot_heights_form);
    TEST(test_listplot_pairs_form);
    TEST(test_listplot_datarange_maps_heights);
    TEST(test_listplot_joined_emits_line);
    TEST(test_listplot_multiple_datasets_get_palette);
    TEST(test_listplot_filling_builds_stems);
    TEST(test_listplot_joined_filling_builds_polygons);
    TEST(test_listplot_default_options_injected);
    TEST(test_listplot_options_registered);
    TEST(test_listplot_legends_metadata);
#ifdef USE_GRAPHICS
    TEST(test_frame_minor_divs_policy);
    TEST(test_window_height_policy);
    TEST(test_window_height_fit_region_no_letterbox);
    TEST(test_raster_dims_no_letterbox);
    TEST(test_polygon_signed_area_winding_detection);
    TEST(test_cmyk_to_rgb_conversion);
    TEST(test_raster_point_size_fraction_of_width);
    TEST(test_raster_dashing);
#endif

    printf("All graphics tests passed!\n");
    symtab_clear();
    return 0;
}
