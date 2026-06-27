/* graphics_init.c — registers the graphics engine's symbols.
 *
 * Primitives (Point, Line, Rectangle, Circle, Disk, Polygon, Text) and
 * style directives (RGBColor, GrayLevel, Opacity, Thickness, PointSize)
 * are, like in Mathematica, inert heads: they get no builtin C function,
 * so the evaluator leaves e.g. Point[{1,2}] unevaluated automatically.
 * Graphics[...] is the same -- nothing computes until Show[]/Plot[]
 * actually renders it. Only Show and Plot need real C logic. */

#include "graphics.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdlib.h>

/* Register a symbol as an inert, protected structural head with a
 * docstring but no builtin function -- mirrors how Rule/List/etc. are
 * "just data" in this codebase. */
static void register_inert(const char* name, const char* doc) {
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(name, doc);
}

/* Mathematica's named color constants are themselves OwnValues that
 * evaluate to an RGBColor[...]/GrayLevel[...] literal -- not inert heads.
 * That's what makes them resolve anywhere a real color literal would: a
 * bare (non-Held) Graphics[]'s arguments evaluate normally, and Plot
 * (HoldAll) now evaluates each of its own option values once before
 * storing them (see plot.c's split_options), so e.g. `Epilog -> {Red, ...}`
 * sees `RGBColor[1, 0, 0]`, exactly as it would for a literal Graphics[]. */
static void register_color(const char* name, Expr* value) {
    Expr* sym = expr_new_symbol(name);
    symtab_add_own_value(name, sym, value);
    expr_free(sym);
    expr_free(value);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    /* Docstrings for the named colours live centrally in info.c (info_init),
     * mirroring CMYKColor; info_init runs before graphics_init. */
}

/* `v` exactly 0 or 1 prints as a plain integer (matching Mathematica's own
 * `FullForm[Red]` -> `RGBColor[1, 0, 0]`); anything else is a real. */
static Expr* color_component(double v) {
    if (v == 0.0) return expr_new_integer(0);
    if (v == 1.0) return expr_new_integer(1);
    return expr_new_real(v);
}

static Expr* rgb_color(double r, double g, double b) {
    Expr* a[3] = { color_component(r), color_component(g), color_component(b) };
    return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
}

static Expr* gray_color(double g) {
    Expr* a[1] = { color_component(g) };
    return expr_new_function(expr_new_symbol(SYM_GrayLevel), a, 1);
}

void graphics_init(void) {
    register_inert("Point",
        "Point[{x, y}]\n\tA graphics primitive: a single point.\n"
        "Point[{{x1,y1}, {x2,y2}, ...}]\n\tA list of points.");
    register_inert("Line",
        "Line[{{x1,y1}, {x2,y2}, ...}]\n\tA graphics primitive: a polyline "
        "through the given points.");
    register_inert("Rectangle",
        "Rectangle[{xmin,ymin}, {xmax,ymax}]\n\tA graphics primitive: an "
        "axis-aligned filled rectangle.");
    register_inert("Circle",
        "Circle[{x,y}, r]\n\tA graphics primitive: a circle outline of "
        "radius r centered at {x,y}. Circle[{x,y}] takes radius 1 and "
        "Circle[] is the unit circle at the origin.");
    register_inert("Disk",
        "Disk[{x,y}, r]\n\tA graphics primitive: a filled disk of radius r "
        "centered at {x,y}. Disk[{x,y}] takes radius 1 and Disk[] is the "
        "unit disk at the origin.");
    register_inert("Polygon",
        "Polygon[{{x1,y1}, {x2,y2}, ...}]\n\tA graphics primitive: a filled "
        "closed polygon through the given vertices.");
    register_inert("Text",
        "Text[expr, {x,y}]\n\tA graphics primitive: expr rendered as text "
        "centered at {x,y}.");

    register_inert("RGBColor",
        "RGBColor[r, g, b]\n\tA style directive: sets the color (each "
        "component in [0,1]) of subsequent primitives.");
    register_inert("GrayLevel",
        "GrayLevel[g]\n\tA style directive: sets the color of subsequent "
        "primitives to a shade of gray (g in [0,1]).");
    register_inert("Hue",
        "Hue[h]\n\tA style directive: sets the color of subsequent "
        "primitives to the fully saturated hue h (in [0,1], wrapping).\n"
        "Hue[h, s, b]\n\tHue h, saturation s, brightness b (each in [0,1]).\n"
        "Hue[h, s, b, a]\n\tAs above, with opacity a.");
    register_inert("Opacity",
        "Opacity[a]\n\tA style directive: sets the opacity (a in [0,1]) of "
        "subsequent primitives.");
    /* CMYKColor is an inert, protected style directive like the colors above;
     * its docstring is set centrally in info.c (info_init). */
    symtab_get_def("CMYKColor")->attributes |= ATTR_PROTECTED;
    register_inert("Thickness",
        "Thickness[t]\n\tA style directive: sets the line thickness (in "
        "plot coordinates) of subsequent Line/Circle primitives.");

    /* Named colour constants -> RGBColor[...]/GrayLevel[...] literals.
     * Black/White/Gray and the grey LightGray use GrayLevel; the rest use
     * RGBColor, matching Mathematica's own InputForm values. */
    register_color("Black",        gray_color(0));
    register_color("White",        gray_color(1));
    register_color("Gray",         gray_color(0.5));
    register_color("Red",          rgb_color(1, 0, 0));
    register_color("Green",        rgb_color(0, 1, 0));
    register_color("Blue",         rgb_color(0, 0, 1));
    register_color("Cyan",         rgb_color(0, 1, 1));
    register_color("Magenta",      rgb_color(1, 0, 1));
    register_color("Yellow",       rgb_color(1, 1, 0));
    register_color("Brown",        rgb_color(0.6, 0.4, 0.2));
    register_color("Orange",       rgb_color(1, 0.5, 0));
    register_color("Pink",         rgb_color(1, 0.5, 0.5));
    register_color("Purple",       rgb_color(0.5, 0, 0.5));
    register_color("LightRed",     rgb_color(1, 0.85, 0.85));
    register_color("LightGreen",   rgb_color(0.88, 1, 0.88));
    register_color("LightBlue",    rgb_color(0.87, 0.94, 1));
    register_color("LightGray",    gray_color(0.85));
    register_color("LightCyan",    rgb_color(0.9, 1, 1));
    register_color("LightMagenta", rgb_color(1, 0.9, 1));
    register_color("LightYellow",  rgb_color(1, 1, 0.85));
    register_color("LightBrown",   rgb_color(0.94, 0.91, 0.88));
    register_color("LightOrange",  rgb_color(1, 0.9, 0.8));
    register_color("LightPink",    rgb_color(1, 0.925, 0.925));
    register_color("LightPurple",  rgb_color(0.94, 0.88, 0.94));
    register_inert("PointSize",
        "PointSize[s]\n\tA style directive: sets the radius (in plot "
        "coordinates) of subsequent Point primitives.");

    register_inert("Graphics",
        "Graphics[primitives, opts...]\n\tA symbolic 2D graphics object. "
        "primitives is a (possibly nested) list of graphics primitives "
        "(Point, Line, Rectangle, Circle, Disk, Polygon, Text) and style "
        "directives (RGBColor, GrayLevel, Hue, Opacity, Thickness, "
        "PointSize). Rendered on demand by Show[]. Prints as -Graphics-. "
        "Options: AspectRatio, Axes, AxesLabel, AxesOrigin, AxesStyle, "
        "Background, Epilog, Frame, FrameLabel, FrameStyle, FrameTicks, "
        "GridLines, GridLinesStyle, ImageSize, LabelStyle, PlotLabel, "
        "PlotRange, PlotRangePadding, PlotStyle, Prolog, RotateLabel, "
        "TicksStyle.");

    /* Internal: Plot embeds $PlotResample[var, {bodies}, {opts...}] inside
     * the Graphics it returns so the renderer can re-sample on zoom. HoldAll
     * keeps the bodies and held option values unevaluated through the
     * Graphics re-evaluation. */
    symtab_get_def("$PlotResample")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("$PlotResample",
        "$PlotResample[var, {f...}, {plotPoints, maxRecursion, maxPlotPoints, "
        "mesh, regionFunction, exclusions, colorFunction, "
        "colorFunctionScaling, filling, fillingStyle}]\n"
        "\tInternal Plot metadata used by the renderer to re-sample curves "
        "at the current zoom. Not intended for direct use.");

    /* Internal: Plot embeds $PlotLegendData[{color,label}, ...] when
     * PlotLegends is given, so the renderer's draw_legend() can show a
     * legend box matching each curve's actual drawn colour. */
    symtab_get_def("$PlotLegendData")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("$PlotLegendData",
        "$PlotLegendData[{color1, label1}, ...]\n"
        "\tInternal Plot metadata read by the renderer to draw a legend box. "
        "Not intended for direct use.");

    symtab_add_builtin("Show", builtin_show);
    symtab_get_def("Show")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Show",
        "Show[graphics, opts...]\n\tDisplays graphics (a Graphics[...] "
        "object) in an interactive window and returns it, merging any "
        "given options into its option list.");

    /* ListPlot's data must evaluate (ListPlot[Range[10]] / Table[...]), so —
     * unlike Plot — it is a plain protected builtin, not HoldAll. The option
     * keywords below are inert so they print and ?-document, but carry no
     * value of their own. */
    register_inert("Joined",
        "Joined\n\tA ListPlot option: True connects the points with a line "
        "instead of drawing markers (default False).");
    register_inert("DataRange",
        "DataRange\n\tA ListPlot option: the x-range {xmin, xmax} to assume "
        "for a list of heights (default Automatic = 1 to n). DataRange -> All "
        "reads a flat list of pairs as several height datasets.");
    register_inert("PlotMarkers",
        "PlotMarkers\n\tA ListPlot option naming the markers used for points "
        "(default None). Currently accepted but rendered with the default "
        "marker.");

    symtab_add_builtin("ListPlot", builtin_listplot);
    symtab_get_def("ListPlot")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ListPlot",
        "ListPlot[{y1, ..., yn}, opts...]\n\tPlots the values as points "
        "{i, yi} (a scatter/point plot). ListPlot[{{x1,y1}, ...}] plots the "
        "given coordinate pairs; ListPlot[{data1, data2, ...}] overlays each "
        "dataset in a distinct palette colour. Returns a Graphics[...] object. "
        "Options: Joined (connect points; default False), DataRange (x-range "
        "for heights), Filling (Axis/Bottom/Top/a number — draws stems), "
        "FillingStyle, PlotMarkers, PlotStyle, PlotLegends, and the Graphics "
        "options PlotRange, Axes (default True), AspectRatio (default "
        "1/GoldenRatio), Frame, AxesLabel, GridLines, ImageSize, Background, "
        "PlotLabel, Prolog, Epilog.");

    symtab_add_builtin("Plot", builtin_plot);
    symtab_get_def("Plot")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Plot",
        "Plot[f, {x, xmin, xmax}, opts...]\n\tAdaptively samples f over "
        "[xmin, xmax], displays the resulting curve in an interactive "
        "window, and returns it as a Graphics[...] object. A list of "
        "functions Plot[{f1, f2, ...}, {x, xmin, xmax}] draws each on the "
        "same axes in a distinct palette colour. Options: "
        "PlotPoints (initial sample count, default 50), MaxRecursion "
        "(adaptive refinement depth, default 6), MaxPlotPoints (overall "
        "point cap, default Infinity), Mesh (All overlays the evaluation "
        "points as dots; default None), PlotRange, PlotRangePadding, "
        "AspectRatio, PlotStyle, Axes, AxesLabel, AxesOrigin, AxesStyle, "
        "TicksStyle, LabelStyle, Frame, FrameLabel, FrameStyle, FrameTicks, "
        "RotateLabel, GridLines, GridLinesStyle, Prolog, Epilog, PlotLabel, "
        "Background, ImageSize, ColorFunction (a function, or \"Rainbow\"), "
        "ColorFunctionScaling (default True), Filling (Axis/Bottom/Top/a "
        "number), FillingStyle, PlotLegends (Automatic/\"Expressions\"/an "
        "explicit list), RegionFunction, Exclusions.");
}
