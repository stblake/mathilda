/* Symbolic-construction tests for ParametricPlot[]/Graphics[]. Never opens a
 * Raylib window -- MATHILDA_NO_GRAPHICS_WINDOW (set as an env-var in main())
 * forces the windowing call to a no-op, keeping the suite headless everywhere. */
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* ---- Head / structural ---- */

void test_parametricplot_returns_graphics_head(void) {
    assert_eval_eq("Head[ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi}]]",
                   "Graphics", 0);
    assert_eval_eq("Head[ParametricPlot[{t, t^2}, {t, -2, 2}]]",
                   "Graphics", 0);
}

void test_parametricplot_invalid_args_unevaluated(void) {
    /* Symbolic bounds: cannot numericize → stays unevaluated. */
    assert_eval_eq("ParametricPlot[{t, t}, {t, a, b}]",
                   "ParametricPlot[{t, t}, {t, a, b}]", 0);
    /* Wrong curve spec — 3-element list is not {fx, fy}. */
    assert_eval_eq("ParametricPlot[{t, t, t}, {t, 0, 1}]",
                   "ParametricPlot[{t, t, t}, {t, 0, 1}]", 0);
    /* Missing iterator spec. */
    assert_eval_eq("ParametricPlot[{t, t}]",
                   "ParametricPlot[{t, t}]", 0);
    /* PlotPoints < 2 is invalid. */
    assert_eval_eq("ParametricPlot[{t, t}, {t, 0, 1}, PlotPoints -> 1]",
                   "ParametricPlot[{t, t}, {t, 0, 1}, PlotPoints -> 1]", 0);
}

/* ---- Default options in the returned Graphics ---- */

/* Default AspectRatio is 1 (square), unlike Plot's 1/GoldenRatio. */
void test_parametricplot_default_aspectratio_one(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{t, t}, {t, 0, 1}]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    bool found = false;
    for (size_t i = 1; i < g->data.function.arg_count; i++) {
        Expr* opt = g->data.function.args[i];
        if (opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2) continue;
        Expr* lhs = opt->data.function.args[0];
        Expr* rhs = opt->data.function.args[1];
        if (lhs->type != EXPR_SYMBOL || strcmp(lhs->data.symbol, "AspectRatio") != 0)
            continue;
        double v = 0.0;
        if      (rhs->type == EXPR_REAL)    v = rhs->data.real;
        else if (rhs->type == EXPR_INTEGER) v = (double)rhs->data.integer;
        found = (v > 0.99 && v < 1.01);
        break;
    }
    ASSERT(found); /* AspectRatio -> 1 must be present */
    expr_free(g);
}

/* Default Axes -> True so the coordinate grid shows. */
void test_parametricplot_default_axes_true(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{t, t}, {t, 0, 1}]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    bool found = false;
    for (size_t i = 1; i < g->data.function.arg_count; i++) {
        Expr* opt = g->data.function.args[i];
        if (opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2) continue;
        Expr* lhs = opt->data.function.args[0];
        Expr* rhs = opt->data.function.args[1];
        if (lhs->type == EXPR_SYMBOL && strcmp(lhs->data.symbol, "Axes") == 0
            && rhs->type == EXPR_SYMBOL && strcmp(rhs->data.symbol, "True") == 0) {
            found = true; break;
        }
    }
    ASSERT(found); /* Axes -> True must be present */
    expr_free(g);
}

/* ---- Primitive list structure: single curve ---- */

/* The primitive list (Graphics[[1]]) must contain at least one Line[]. */
void test_parametricplot_single_curve_has_line(void) {
    assert_eval_eq(
        "Length[Cases[ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi}][[1]],"
        " Line[___]]] >= 1",
        "True", 0);
}

/* MaxRecursion -> 0 means no adaptive subdivision: exactly PlotPoints samples
 * are produced.  The single Line[...] for {t, t} (a straight line) has its
 * coordinate list at [[1, 1, 1]]. */
void test_parametricplot_maxrecursion0_exact_count(void) {
    assert_eval_eq(
        "Length[ParametricPlot[{t, t}, {t, 0, 1},"
        " PlotPoints -> 10, MaxRecursion -> 0][[1, 1, 1]]]",
        "10", 0);
    assert_eval_eq(
        "Length[ParametricPlot[{t, t}, {t, 0, 1},"
        " PlotPoints -> 20, MaxRecursion -> 0][[1, 1, 1]]]",
        "20", 0);
}

/* The adaptive sampler refines curved paths: a circle gets more than the
 * initial PlotPoints=10 samples. */
void test_parametricplot_adaptive_refines_circle(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi},"
        " PlotPoints -> 10, MaxRecursion -> 6]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    Expr* prim_list = g->data.function.args[0];
    size_t npts = 0;
    for (size_t i = 0; i < prim_list->data.function.arg_count; i++) {
        Expr* e = prim_list->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "Line") == 0)
            npts += e->data.function.args[0]->data.function.arg_count;
    }
    ASSERT(npts > 10); /* adaptive refinement adds more points */
    expr_free(g);
}

/* Each {x,y} point in the Line must be a 2-element List. */
void test_parametricplot_points_are_2d(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi},"
        " PlotPoints -> 8, MaxRecursion -> 0]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    Expr* prim_list = g->data.function.args[0];
    for (size_t i = 0; i < prim_list->data.function.arg_count; i++) {
        Expr* e = prim_list->data.function.args[i];
        if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL
            || strcmp(e->data.function.head->data.symbol, "Line") != 0) continue;
        Expr* pts = e->data.function.args[0];
        for (size_t j = 0; j < pts->data.function.arg_count; j++) {
            Expr* pt = pts->data.function.args[j];
            ASSERT(pt->type == EXPR_FUNCTION && pt->data.function.arg_count == 2);
        }
    }
    expr_free(g);
}

/* ---- Multi-curve ---- */

/* Two palette RGBColor directives and two Line[] runs when given two curves. */
void test_parametricplot_multi_curve_two_palette_colors(void) {
    assert_eval_eq(
        "Length[Cases[ParametricPlot[{{Cos[t], Sin[t]}, {2 Cos[t], Sin[t]}},"
        " {t, 0, 2 Pi}][[1]], _RGBColor]]",
        "2", 0);
    assert_eval_eq(
        "Length[Cases[ParametricPlot[{{Cos[t], Sin[t]}, {2 Cos[t], Sin[t]}},"
        " {t, 0, 2 Pi}][[1]], Line[___]]] >= 2",
        "True", 0);
}

/* ---- PlotPoints controls density ---- */

/* Lower PlotPoints + MaxRecursion -> 0 yields fewer points. */
void test_parametricplot_plotpoints_controls_density(void) {
    /* MaxRecursion -> 0 means exactly PlotPoints points; verified by
     * comparing the point counts of two distinct PlotPoints settings. */
    Expr* g8 = evaluate(parse_expression(
        "ParametricPlot[{t, t}, {t, 0, 1},"
        " PlotPoints -> 8, MaxRecursion -> 0]"));
    Expr* g50 = evaluate(parse_expression(
        "ParametricPlot[{t, t}, {t, 0, 1},"
        " PlotPoints -> 50, MaxRecursion -> 0]"));
    ASSERT(g8  && g8->type  == EXPR_FUNCTION);
    ASSERT(g50 && g50->type == EXPR_FUNCTION);

    size_t n8 = 0, n50 = 0;
    Expr* pl8  = g8->data.function.args[0];
    Expr* pl50 = g50->data.function.args[0];
    for (size_t i = 0; i < pl8->data.function.arg_count; i++) {
        Expr* e = pl8->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "Line") == 0)
            n8 += e->data.function.args[0]->data.function.arg_count;
    }
    for (size_t i = 0; i < pl50->data.function.arg_count; i++) {
        Expr* e = pl50->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "Line") == 0)
            n50 += e->data.function.args[0]->data.function.arg_count;
    }
    ASSERT(n8 < n50);
    expr_free(g8); expr_free(g50);
}

/* ---- Mesh ---- */

/* Mesh -> All adds a PointSize[] and a Point[] primitive. */
void test_parametricplot_mesh_adds_pointsize_point(void) {
    assert_eval_eq(
        "Length[Cases[ParametricPlot[{t, t}, {t, 0, 1}, Mesh -> All][[1]],"
        " _PointSize]]",
        "1", 0);
    assert_eval_eq(
        "Length[Cases[ParametricPlot[{t, t}, {t, 0, 1}, Mesh -> All][[1]],"
        " _Point]]",
        "1", 0);
}

/* ---- ColorFunction ---- */

/* ColorFunction -> (Hue[#] &) produces alternating color + 2-pt Line pairs.
 * With MaxRecursion -> 0 and PlotPoints -> 5 there are 4 segments. */
void test_parametricplot_colorfn_segment_pairs(void) {
    assert_eval_eq(
        "Length[Cases[ParametricPlot[{t, t}, {t, 0, 1},"
        " ColorFunction -> (Hue[#] &), PlotPoints -> 5, MaxRecursion -> 0"
        "][[1]], Line[___]]]",
        "4", 0); /* 5 points → 4 segments → 4 Line[] primitives */
}

/* Every segment Line[] must have exactly 2 points. */
void test_parametricplot_colorfn_lines_are_2pt(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi},"
        " ColorFunction -> (Hue[#] &), PlotPoints -> 8, MaxRecursion -> 0]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    Expr* prim_list = g->data.function.args[0];
    int line_count = 0;
    for (size_t i = 0; i < prim_list->data.function.arg_count; i++) {
        Expr* e = prim_list->data.function.args[i];
        if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL
            || strcmp(e->data.function.head->data.symbol, "Line") != 0) continue;
        size_t npt = e->data.function.args[0]->data.function.arg_count;
        ASSERT(npt == 2);
        line_count++;
    }
    ASSERT(line_count > 0);
    expr_free(g);
}

/* ---- RegionFunction ---- */

/* RegionFunction -> Function[{x,y}, x > 0] restricts to the right half-plane.
 * The Graphics is still returned; all x-coordinates in Line[] points are >= 0. */
void test_parametricplot_regionfunction_right_half(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi},"
        " RegionFunction -> Function[{x, y}, x > 0]]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    Expr* prim_list = g->data.function.args[0];
    ASSERT(prim_list->data.function.arg_count > 0);

    for (size_t i = 0; i < prim_list->data.function.arg_count; i++) {
        Expr* e = prim_list->data.function.args[i];
        if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL
            || strcmp(e->data.function.head->data.symbol, "Line") != 0) continue;
        Expr* pts = e->data.function.args[0];
        for (size_t j = 0; j < pts->data.function.arg_count; j++) {
            Expr* pt = pts->data.function.args[j];
            double xv = 0.0;
            if (pt->data.function.args[0]->type == EXPR_REAL)
                xv = pt->data.function.args[0]->data.real;
            ASSERT(xv >= -0.1); /* small tolerance for boundary bisection */
        }
    }
    expr_free(g);
}

/* ---- PlotStyle pass-through ---- */

/* A supplied RGBColor PlotStyle should appear in the Graphics option list. */
void test_parametricplot_plotstyle_passthrough(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{t, t^2}, {t, -1, 1},"
        " PlotStyle -> RGBColor[1, 0, 0]]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    bool found_rgb_style = false;
    for (size_t i = 1; i < g->data.function.arg_count; i++) {
        Expr* opt = g->data.function.args[i];
        if (opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2) continue;
        Expr* lhs = opt->data.function.args[0];
        Expr* rhs = opt->data.function.args[1];
        if (lhs->type != EXPR_SYMBOL || strcmp(lhs->data.symbol, "PlotStyle") != 0)
            continue;
        found_rgb_style = (rhs->type == EXPR_FUNCTION
            && rhs->data.function.head->type == EXPR_SYMBOL
            && strcmp(rhs->data.function.head->data.symbol, "RGBColor") == 0);
        break;
    }
    ASSERT(found_rgb_style);
    expr_free(g);
}

/* ---- Custom AspectRatio ---- */

void test_parametricplot_custom_aspectratio(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{t, t}, {t, 0, 1}, AspectRatio -> 2]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    bool found = false;
    for (size_t i = 1; i < g->data.function.arg_count; i++) {
        Expr* opt = g->data.function.args[i];
        if (opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2) continue;
        Expr* lhs = opt->data.function.args[0];
        Expr* rhs = opt->data.function.args[1];
        if (lhs->type != EXPR_SYMBOL || strcmp(lhs->data.symbol, "AspectRatio") != 0)
            continue;
        double v = 0.0;
        if      (rhs->type == EXPR_REAL)    v = rhs->data.real;
        else if (rhs->type == EXPR_INTEGER) v = (double)rhs->data.integer;
        found = (v > 1.99 && v < 2.01);
        break;
    }
    ASSERT(found);
    expr_free(g);
}

/* ---- Lissajous figure: basic sanity ---- */

/* A Lissajous curve {Sin[2t], Sin[3t]} is smooth and closed.  Check that
 * ParametricPlot produces a non-trivial number of sample points (adaptive
 * refinement fires on the high-curvature regions near the self-crossings). */
void test_parametricplot_lissajous_sufficient_points(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{Sin[2 t], Sin[3 t]}, {t, 0, 2 Pi},"
        " PlotPoints -> 25]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    Expr* prim_list = g->data.function.args[0];
    size_t total_pts = 0;
    for (size_t i = 0; i < prim_list->data.function.arg_count; i++) {
        Expr* e = prim_list->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "Line") == 0)
            total_pts += e->data.function.args[0]->data.function.arg_count;
    }
    ASSERT(total_pts >= 25); /* at minimum the initial grid */
    expr_free(g);
}

/* ---- Computed (non-literal-list) body in 1-iterator form ---- */

/* Body is Times[2, List[...]] — head is Times, not List, but it evaluates to
 * a {x,y} list.  Old code would return NULL here; verify we get Graphics. */
void test_parametricplot_computed_body_1iter(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[2 {Cos[t], Sin[t]}, {t, 0, 2 Pi}]"));
    ASSERT(g && g->type == EXPR_FUNCTION
           && g->data.function.head->type == EXPR_SYMBOL
           && strcmp(g->data.function.head->data.symbol, "Graphics") == 0);
    Expr* prim_list = g->data.function.args[0];
    bool has_line = false;
    for (size_t i = 0; i < prim_list->data.function.arg_count; i++) {
        Expr* e = prim_list->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "Line") == 0) {
            has_line = true; break;
        }
    }
    ASSERT(has_line);
    expr_free(g);
}

/* ---- 2-iterator (filled region) form ---- */

/* Basic: body is a literal {r*Cos[t], r*Sin[t]} with two iterators.
 * Expect Graphics[] head and at least one Polygon primitive. */
void test_parametricplot_two_iter_returns_graphics(void) {
    assert_eval_eq(
        "Head[ParametricPlot[{r Cos[t], r Sin[t]}, {t, 0, 2 Pi}, {r, 1, 2}]]",
        "Graphics", 0);
}

void test_parametricplot_two_iter_has_polygon(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{r Cos[t], r Sin[t]}, {t, 0, 2 Pi}, {r, 1, 2},"
        " PlotPoints -> 5]"));
    ASSERT(g && g->type == EXPR_FUNCTION);
    Expr* prim_list = g->data.function.args[0];
    bool has_polygon = false;
    for (size_t i = 0; i < prim_list->data.function.arg_count; i++) {
        Expr* e = prim_list->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "Polygon") == 0) {
            has_polygon = true; break;
        }
    }
    ASSERT(has_polygon);
    expr_free(g);
}

/* Computed body with two iterators: r^2 * {Sqrt[t]*Cos[t], Sin[t]}.
 * This is the exact expression the user reported as broken. */
void test_parametricplot_two_iter_computed_body(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[r^2 {Sqrt[t] Cos[t], Sin[t]},"
        " {t, 0, 3 Pi/2}, {r, 1, 2}, PlotPoints -> 5]"));
    ASSERT(g && g->type == EXPR_FUNCTION
           && g->data.function.head->type == EXPR_SYMBOL
           && strcmp(g->data.function.head->data.symbol, "Graphics") == 0);
    /* Must produce at least one Polygon (some cells may be invalid for small t). */
    Expr* prim_list = g->data.function.args[0];
    bool has_polygon = false;
    for (size_t i = 0; i < prim_list->data.function.arg_count; i++) {
        Expr* e = prim_list->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "Polygon") == 0) {
            has_polygon = true; break;
        }
    }
    ASSERT(has_polygon);
    expr_free(g);
}

/* Each Polygon in the 2-iter form must have exactly 4 vertices (quad mesh). */
void test_parametricplot_two_iter_polygon_has_4_verts(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{r Cos[t], r Sin[t]}, {t, 0, 2 Pi}, {r, 1, 2},"
        " PlotPoints -> 4]"));
    ASSERT(g && g->type == EXPR_FUNCTION);
    Expr* prim_list = g->data.function.args[0];
    for (size_t i = 0; i < prim_list->data.function.arg_count; i++) {
        Expr* e = prim_list->data.function.args[i];
        if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL
            || strcmp(e->data.function.head->data.symbol, "Polygon") != 0)
            continue;
        /* Polygon[verts_list] — verts_list has 4 elements. */
        ASSERT(e->data.function.arg_count == 1);
        ASSERT(e->data.function.args[0]->data.function.arg_count == 4);
    }
    expr_free(g);
}

/* Mesh -> All in the 2-iter form adds Line primitives. */
void test_parametricplot_two_iter_mesh_adds_lines(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{r Cos[t], r Sin[t]}, {t, 0, 2 Pi}, {r, 1, 2},"
        " PlotPoints -> 4, Mesh -> All]"));
    ASSERT(g && g->type == EXPR_FUNCTION);
    Expr* prim_list = g->data.function.args[0];
    bool has_line = false;
    for (size_t i = 0; i < prim_list->data.function.arg_count; i++) {
        Expr* e = prim_list->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "Line") == 0) {
            has_line = true; break;
        }
    }
    ASSERT(has_line);
    expr_free(g);
}

/* ---- 2-iterator multi-surface form ---- */

/* {{body1},{body2}} with two iterators: each sub-body gets its own
 * color-keyed region.  This was the user-reported broken case. */
void test_parametricplot_two_iter_multi_surface_has_polygon(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{{2 r Cos[t], r Sin[t]}, {r Cos[t], 2 r Sin[t]}},"
        " {t, 0, 2 Pi}, {r, 0, 1}, PlotPoints -> 4]"));
    ASSERT(g && g->type == EXPR_FUNCTION
           && strcmp(g->data.function.head->data.symbol, "Graphics") == 0);
    Expr* prim_list = g->data.function.args[0];
    bool has_polygon = false;
    for (size_t i = 0; i < prim_list->data.function.arg_count; i++) {
        Expr* e = prim_list->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "Polygon") == 0) {
            has_polygon = true; break;
        }
    }
    ASSERT(has_polygon);
    expr_free(g);
}

/* Two surfaces → two distinct palette color directives in the prim list. */
void test_parametricplot_two_iter_multi_surface_two_colors(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{{2 r Cos[t], r Sin[t]}, {r Cos[t], 2 r Sin[t]}},"
        " {t, 0, 2 Pi}, {r, 0, 1}, PlotPoints -> 4]"));
    ASSERT(g && g->type == EXPR_FUNCTION);
    Expr* prim_list = g->data.function.args[0];
    size_t ncolors = 0;
    for (size_t i = 0; i < prim_list->data.function.arg_count; i++) {
        Expr* e = prim_list->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "RGBColor") == 0)
            ncolors++;
    }
    ASSERT(ncolors == 2);
    expr_free(g);
}

/* Mesh -> False is accepted (no PFAIL). */
void test_parametricplot_two_iter_multi_mesh_false(void) {
    assert_eval_eq(
        "Head[ParametricPlot[{{2 r Cos[t], r Sin[t]}, {r Cos[t], 2 r Sin[t]}},"
        " {t, 0, 2 Pi}, {r, 0, 1}, Mesh -> False]]",
        "Graphics", 0);
}

/* ---- PlotLegends ---- */

/* PlotLegends -> Automatic embeds $PlotLegendData in the Graphics object. */
void test_parametricplot_plotlegends_automatic(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi},"
        " PlotLegends -> Automatic]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    bool has_legend = false;
    for (size_t i = 1; i < g->data.function.arg_count; i++) {
        Expr* e = g->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "$PlotLegendData") == 0) {
            has_legend = true; break;
        }
    }
    ASSERT(has_legend);
    expr_free(g);
}

/* PlotLegends -> None produces no $PlotLegendData. */
void test_parametricplot_plotlegends_none(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi},"
        " PlotLegends -> None]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    bool has_legend = false;
    for (size_t i = 1; i < g->data.function.arg_count; i++) {
        Expr* e = g->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "$PlotLegendData") == 0) {
            has_legend = true; break;
        }
    }
    ASSERT(!has_legend);
    expr_free(g);
}

/* Multi-curve with PlotLegends: $PlotLegendData has one entry per curve. */
void test_parametricplot_plotlegends_multicurve(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{{Cos[t], Sin[t]}, {2 Cos[t], Sin[t]}},"
        " {t, 0, 2 Pi}, PlotLegends -> Automatic]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    bool found = false;
    for (size_t i = 1; i < g->data.function.arg_count; i++) {
        Expr* e = g->data.function.args[i];
        if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && strcmp(e->data.function.head->data.symbol, "$PlotLegendData") == 0) {
            ASSERT(e->data.function.arg_count == 2); /* two curves */
            found = true; break;
        }
    }
    ASSERT(found);
    expr_free(g);
}

/* Explicit list of labels: $PlotLegendData entries use those labels. */
void test_parametricplot_plotlegends_explicit_labels(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi},"
        " PlotLegends -> {\"my curve\"}]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    bool found_label = false;
    for (size_t i = 1; i < g->data.function.arg_count; i++) {
        Expr* e = g->data.function.args[i];
        if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL
            || strcmp(e->data.function.head->data.symbol, "$PlotLegendData") != 0) continue;
        /* $PlotLegendData[{color, label}]: check label is the string "my curve". */
        if (e->data.function.arg_count == 1) {
            Expr* entry = e->data.function.args[0];
            if (entry->type == EXPR_FUNCTION && entry->data.function.arg_count == 2) {
                Expr* lbl = entry->data.function.args[1];
                found_label = (lbl->type == EXPR_STRING
                               && strcmp(lbl->data.string, "my curve") == 0);
            }
        }
        break;
    }
    ASSERT(found_label);
    expr_free(g);
}

/* ---- Display options pass through ---- */

/* PlotLabel, GridLines, Epilog pass through to the Graphics result. */
void test_parametricplot_display_options_passthrough(void) {
    Expr* g = evaluate(parse_expression(
        "ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi},"
        " PlotLabel -> \"Circle\", GridLines -> Automatic]"));
    ASSERT(g && g->type == EXPR_FUNCTION);

    bool has_label = false, has_grid = false;
    for (size_t i = 1; i < g->data.function.arg_count; i++) {
        Expr* opt = g->data.function.args[i];
        if (opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2) continue;
        Expr* lhs = opt->data.function.args[0];
        if (lhs->type != EXPR_SYMBOL) continue;
        if (strcmp(lhs->data.symbol, "PlotLabel")  == 0) has_label = true;
        if (strcmp(lhs->data.symbol, "GridLines")  == 0) has_grid  = true;
    }
    ASSERT(has_label);
    ASSERT(has_grid);
    expr_free(g);
}

/* ---- runner ---- */

int main(void) {
    setenv("MATHILDA_NO_GRAPHICS_WINDOW", "1", 1);
    symtab_init();
    core_init();

    printf("=== ParametricPlot tests ===\n");

    TEST(test_parametricplot_returns_graphics_head);
    TEST(test_parametricplot_invalid_args_unevaluated);
    TEST(test_parametricplot_default_aspectratio_one);
    TEST(test_parametricplot_default_axes_true);
    TEST(test_parametricplot_single_curve_has_line);
    TEST(test_parametricplot_maxrecursion0_exact_count);
    TEST(test_parametricplot_adaptive_refines_circle);
    TEST(test_parametricplot_points_are_2d);
    TEST(test_parametricplot_multi_curve_two_palette_colors);
    TEST(test_parametricplot_plotpoints_controls_density);
    TEST(test_parametricplot_mesh_adds_pointsize_point);
    TEST(test_parametricplot_colorfn_segment_pairs);
    TEST(test_parametricplot_colorfn_lines_are_2pt);
    TEST(test_parametricplot_regionfunction_right_half);
    TEST(test_parametricplot_plotstyle_passthrough);
    TEST(test_parametricplot_custom_aspectratio);
    TEST(test_parametricplot_lissajous_sufficient_points);
    TEST(test_parametricplot_computed_body_1iter);
    TEST(test_parametricplot_two_iter_returns_graphics);
    TEST(test_parametricplot_two_iter_has_polygon);
    TEST(test_parametricplot_two_iter_computed_body);
    TEST(test_parametricplot_two_iter_polygon_has_4_verts);
    TEST(test_parametricplot_two_iter_mesh_adds_lines);
    TEST(test_parametricplot_two_iter_multi_surface_has_polygon);
    TEST(test_parametricplot_two_iter_multi_surface_two_colors);
    TEST(test_parametricplot_two_iter_multi_mesh_false);
    TEST(test_parametricplot_plotlegends_automatic);
    TEST(test_parametricplot_plotlegends_none);
    TEST(test_parametricplot_plotlegends_multicurve);
    TEST(test_parametricplot_plotlegends_explicit_labels);
    TEST(test_parametricplot_display_options_passthrough);

    printf("All ParametricPlot tests passed.\n");
    symtab_clear();
    return 0;
}
