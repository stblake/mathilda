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
        "primitives is a (possibly nested) list of graphics primitives and "
        "style directives. Rendered on demand by Show[]. Prints as "
        "-Graphics-.");
    register_inert("Graphics3D",
        "Graphics3D[primitives, opts...]\n\tA symbolic 3D graphics object, "
        "as built by Plot3D[]. primitives is a (possibly nested) list of "
        "the same primitives Graphics[] uses (Polygon, Line, color "
        "directives, ...) but with 3-coordinate {x,y,z} points instead of "
        "2-coordinate {x,y}. Rendered on demand by Show[] in an orbit-"
        "camera window: drag to rotate, scroll to zoom, right-drag to pan, "
        "R to reset, S to save a screenshot. Prints as -Graphics3D-.");

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

    /* Internal: StreamPlot embeds $StreamColorBar[spd_min, spd_max] when
     * PlotLegends -> Automatic is set and speed-based coloring is active.
     * The renderer draws a vertical gradient color scale bar from spd_min
     * (dark blue-purple) to spd_max (bright yellow). */
    symtab_get_def("$StreamColorBar")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("$StreamColorBar",
        "$StreamColorBar[speed_min, speed_max]\n"
        "\tInternal StreamPlot metadata. Instructs the renderer to draw a "
        "vertical speed color scale bar (dark blue = slow, yellow = fast). "
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
        "Background, ImageSize, ColorFunction (a function, or named ramp: "
        "\"Rainbow\"/\"CoolTones\"/\"WarmTones\"/\"Greyscale\"/\"Temperature\"), "
        "ColorFunctionScaling (default True), Filling (Axis/Bottom/Top/a "
        "number), FillingStyle, PlotLegends (Automatic/\"Expressions\"/an "
        "explicit list), RegionFunction, Exclusions.");

    symtab_add_builtin("ParametricPlot", builtin_parametricplot);
    symtab_get_def("ParametricPlot")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("ParametricPlot",
        "ParametricPlot[{fx, fy}, {t, tmin, tmax}, opts...]\n"
        "\tAdaptively samples the parametric curve (fx(t), fy(t)) over "
        "[tmin, tmax] and returns a Graphics[...] object (auto-displayed). "
        "The body may be any expression that evaluates to a 2-element {x,y} "
        "list (not just a literal {fx,fy}). Multiple curves: "
        "ParametricPlot[{{fx1,fy1}, ...}, {t,...}]. Two-iterator (filled "
        "region) form: ParametricPlot[body, {t,...}, {r,...}] samples a "
        "PlotPoints x PlotPoints grid and emits Polygon[] quads. "
        "Default AspectRatio -> 1 (both axes equally important). "
        "Options: PlotPoints (default 25), MaxRecursion (default 6), "
        "MaxPlotPoints, Mesh (All: dots for curves, grid lines for regions), "
        "PlotLegends (Automatic/\"Expressions\"/{labels...}: draws a legend), "
        "ColorFunction (\"Rainbow\" or f[t] / f[t,r]), ColorFunctionScaling "
        "(default True), RegionFunction (f[x,y] mask), PlotStyle, "
        "AspectRatio, Axes, PlotRange, PlotRangePadding, AxesLabel, "
        "AxesOrigin, Frame, FrameLabel, GridLines, Prolog, Epilog, "
        "PlotLabel, Background, ImageSize (all passed through to Graphics).");

    symtab_add_builtin("ParametricPlot3D", builtin_parametricplot3d);
    symtab_get_def("ParametricPlot3D")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("ParametricPlot3D",
        "ParametricPlot3D[{fx, fy, fz}, {t, tmin, tmax}, opts...]\n"
        "\tAdaptively samples the parametric 3D space curve (fx(t), fy(t), fz(t)) "
        "over [tmin, tmax] and returns a Graphics3D[...] object (auto-displayed in "
        "an orbit-camera window). The body may be any expression that evaluates to "
        "a 3-element {x,y,z} list. Multiple curves: "
        "ParametricPlot3D[{{fx1,fy1,fz1}, ...}, {t,...}].\n"
        "ParametricPlot3D[{fx, fy, fz}, {t, tmin, tmax}, {u, umin, umax}, opts...]\n"
        "\tTwo-iterator form: samples a PlotPoints x PlotPoints grid of (t,u) "
        "pairs, maps each to {x,y,z}, and emits Polygon[] quads — a parametric "
        "3D surface patch. Options: "
        "PlotPoints (initial sample count/grid size, default 25), "
        "MaxRecursion (adaptive refinement depth for curves, default 6), "
        "MaxPlotPoints (overall point cap, default Infinity), "
        "Mesh (All/True: overlays sample dots for curves or grid lines for "
        "surfaces; default None), "
        "PlotLegends (Automatic/\"Expressions\"/{labels...}), "
        "ColorFunction (\"Rainbow\" or f[x,y,z] receiving scaled spatial coords, "
        "or f[x,z] / f[z] for height-based coloring), "
        "ColorFunctionScaling (default True), "
        "RegionFunction (f[x,y,z] mask; falls back to f[x,y] forms), "
        "PlotStyle, Axes, PlotRange, AxesLabel, PlotLabel, Background, ImageSize "
        "(all passed through to Graphics3D). "
        "Lighting -> None disables shading (flat colors); default is Automatic "
        "(Lambertian shading, same as Plot3D).");

    /* Lighting is an inert option keyword recognised by the renderer. */
    symtab_get_def("Lighting")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Lighting",
        "Lighting\n\tA Graphics3D / Plot3D / ParametricPlot3D option controlling "
        "surface shading. Lighting -> Automatic (default): per-face Lambertian "
        "(flat) shading with a fixed directional light; ambient 0.3, diffuse 0.7. "
        "Lighting -> None (or False): disables shading and draws surfaces in their "
        "raw PlotStyle/ColorFunction color.");

    symtab_add_builtin("ContourPlot", builtin_contourplot);
    symtab_get_def("ContourPlot")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("ContourPlot",
        "ContourPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]\n"
        "\tGenerates iso-contour lines of f(x,y) using the marching squares algorithm\n"
        "\tand returns a Graphics[...] object (auto-displayed). ContourPlot is HoldAll:\n"
        "\tf is held unevaluated until x and y are bound to numeric values.\n"
        "\n"
        "\tOptions:\n"
        "\t  Contours         - Integer n (n evenly spaced auto levels, default 10), or\n"
        "\t                     {c1, c2, ...} (explicit contour values).\n"
        "\t  ContourStyle     - Style directive(s) for the contour lines. A single\n"
        "\t                     directive is applied to all levels; a List cycles through\n"
        "\t                     the levels. Automatic (default) colours by height.\n"
        "\t                     None/False suppresses lines (leaves only shading).\n"
        "\t  ContourLabels    - True: draw the z value at the midpoint of each level's\n"
        "\t                     first visible segment. Default False.\n"
        "\t  ContourShading   - True: fill each grid cell by its z value (via\n"
        "\t                     ColorFunction or the built-in thermal gradient).\n"
        "\t                     False/None: lines only. Automatic (default): shade when\n"
        "\t                     ColorFunction is set, otherwise lines only.\n"
        "\t  ColorFunction    - A function f[t] → color (t in [0,1] after scaling), or\n"
        "\t                     a named ramp string: \"Rainbow\", \"Temperature\",\n"
        "\t                     \"CoolTones\", \"WarmTones\", \"Greyscale\".\n"
        "\t                     Applied to shading and auto line colors.\n"
        "\t  ColorFunctionScaling - True (default): normalise z to [0,1] before calling\n"
        "\t                         ColorFunction. False: pass raw z.\n"
        "\t  PlotPoints       - Grid resolution per axis (default 25; increase for\n"
        "\t                     smoother contours).\n"
        "\t  RegionFunction   - f[x,y] mask: cells where the function is False are\n"
        "\t                     skipped (neither shaded nor contoured).\n"
        "\t  Standard Graphics options (Axes, AspectRatio, Frame, PlotRange,\n"
        "\t  AxesLabel, GridLines, ImageSize, Background, PlotLabel, Prolog,\n"
        "\t  Epilog, ...) pass through to the Graphics[...] result.");

    /* ContourStyle, Contours, ContourLabels, ContourShading are inert option
     * keywords that print and ?-document but carry no value of their own. */
    symtab_get_def("Contours")->attributes     |= ATTR_PROTECTED;
    symtab_set_docstring("Contours",
        "Contours\n\tContourPlot option: integer count of auto-levels or explicit "
        "list of contour values.");
    symtab_get_def("ContourStyle")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ContourStyle",
        "ContourStyle\n\tContourPlot option: style directive(s) for contour lines. "
        "A list cycles through the levels; Automatic colours by height; None suppresses lines.");
    symtab_get_def("ContourLabels")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ContourLabels",
        "ContourLabels\n\tContourPlot option: True draws z-value text labels at the "
        "first visible point of each contour level. Default False.");
    symtab_get_def("ContourShading")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ContourShading",
        "ContourShading\n\tContourPlot option: True/False/Automatic — fills cells "
        "by z value using ColorFunction or the built-in blue-cyan-yellow-red thermal "
        "ramp. Automatic enables shading only when ColorFunction is set.");

    symtab_add_builtin("PolarPlot", builtin_polarplot);
    symtab_get_def("PolarPlot")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("PolarPlot",
        "PolarPlot[r, {theta, tmin, tmax}, opts...]\n"
        "\tPlots the polar curve r(theta) by converting to Cartesian\n"
        "\tcoordinates {r*Cos[theta], r*Sin[theta]} and sampling adaptively\n"
        "\tover [tmin, tmax]. Returns a Graphics[...] object (auto-displayed).\n"
        "PolarPlot[{r1, r2, ...}, {theta, tmin, tmax}, opts...]\n"
        "\tMultiple polar curves in distinct palette colours.\n"
        "\n"
        "\tNegative r values are plotted in the opposite direction (standard\n"
        "\tpolar convention). Default AspectRatio -> 1 (equal axes).\n"
        "\n"
        "\tOptions (same as ParametricPlot):\n"
        "\t  PlotPoints          - initial sample count per curve (default 75)\n"
        "\t  MaxRecursion        - adaptive refinement depth (default 6)\n"
        "\t  MaxPlotPoints       - total point cap (default Infinity)\n"
        "\t  Mesh                - All/True: overlay evaluation dots; None (default)\n"
        "\t  ColorFunction       - f[t] or \"Rainbow\" (sweeps scaled theta)\n"
        "\t  ColorFunctionScaling - True (default): normalise theta to [0,1]\n"
        "\t  RegionFunction      - f[x,y] mask\n"
        "\t  PlotStyle           - color/style directive(s)\n"
        "\t  PlotLegends         - Automatic / \"Expressions\" / label list\n"
        "\t  PolarAxes           - option keyword (accepted; polar grid overlay\n"
        "\t                        not yet rendered)\n"
        "\t  Standard Graphics options (AspectRatio, Axes, PlotRange,\n"
        "\t  AxesLabel, Frame, GridLines, PlotLabel, Background, ImageSize,\n"
        "\t  Prolog, Epilog) pass through to the Graphics[...] result.");

    /* PolarAxes is an inert option keyword -- recognised and documented but
     * actual polar grid overlay is not yet rendered. */
    symtab_get_def("PolarAxes")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PolarAxes",
        "PolarAxes\n\tPolarPlot option: True requests a polar grid overlay\n"
        "\t(radial circles + angle labels). Currently accepted but not yet\n"
        "\trendered; Cartesian axes are drawn instead.");

    symtab_add_builtin("DensityPlot", builtin_densityplot);
    symtab_get_def("DensityPlot")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("DensityPlot",
        "DensityPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]\n"
        "\tRenders f(x,y) as a heatmap: each grid cell is coloured by its\n"
        "\tfunction value via ColorFunction (default: thermal blue→yellow ramp).\n"
        "\tDensityPlot is HoldAll: f is held unevaluated until x and y are\n"
        "\tbound to numeric values. Returns a Graphics[...] object.\n"
        "\n"
        "\tOptions:\n"
        "\t  PlotPoints          grid resolution per axis (default 50)\n"
        "\t  ColorFunction       named ramp string or f[t]→color (t in [0,1]).\n"
        "\t                      Ramps: \"Rainbow\", \"CoolTones\", \"WarmTones\",\n"
        "\t                      \"Greyscale\", \"Temperature\" (all keyed to\n"
        "\t                      normalised z value, t∈[0,1])\n"
        "\t  ColorFunctionScaling True (default): normalise z to [0,1] before\n"
        "\t                       calling ColorFunction; False: pass raw z\n"
        "\t  RegionFunction      f[x,y] mask; excluded cells are not drawn\n"
        "\t  PlotLegends         Automatic: attach a vertical color scale bar\n"
        "\t  Standard Graphics options (Axes, AspectRatio→1, Frame, PlotRange,\n"
        "\t  AxesLabel, GridLines, ImageSize, Background, PlotLabel, …) pass\n"
        "\t  through to the Graphics[...] result.");

    symtab_add_builtin("BarChart", builtin_barchart);
    symtab_get_def("BarChart")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("BarChart",
        "BarChart[{v1, v2, ..., vn}, opts...]\n"
        "\tDraws a vertical bar chart: n bars at x = 1..n with heights v1..vn.\n"
        "BarChart[{{v1,...}, {w1,...}, ...}, opts...]\n"
        "\tMultiple grouped datasets, each in a distinct palette colour.\n"
        "\n"
        "\tOptions:\n"
        "\t  ChartStyle    color/style list cycling through bars (default: palette)\n"
        "\t  ChartLabels   list of x-axis tick labels\n"
        "\t  BarSpacing    gap fraction of bar width (default 0.2)\n"
        "\t  Standard Graphics options (Axes, AspectRatio, Frame, PlotRange,\n"
        "\t  PlotLabel, Background, ImageSize, …) pass through.");

    symtab_add_builtin("Histogram", builtin_histogram);
    symtab_get_def("Histogram")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Histogram",
        "Histogram[data, opts...]\n"
        "\tBins the numeric values in data and draws a frequency histogram.\n"
        "\tBin count defaults to Sturges' rule: ceil(Log2[n]) + 1.\n"
        "Histogram[data, k, opts...]\n"
        "\tk equal-width bins.\n"
        "Histogram[data, {step}, opts...]\n"
        "\tBins of width step.\n"
        "Histogram[data, {min, max, step}, opts...]\n"
        "\tExplicit range and width.\n"
        "\n"
        "\tOptions:\n"
        "\t  ChartStyle   color/style list cycling through bins\n"
        "\t  BarSpacing   gap fraction of bin width (default 0.2)\n"
        "\t  Standard Graphics options pass through.");

    /* Option keyword inerts for BarChart/Histogram. */
    symtab_get_def("BarSpacing")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("BarSpacing",
        "BarSpacing\n\tBarChart/Histogram option: gap between bars as a fraction\n"
        "\tof bar width (default 0.2). 0 = touching bars; 1 = all gap.");
    symtab_get_def("ChartStyle")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ChartStyle",
        "ChartStyle\n\tBarChart/Histogram option: color or list of colors cycling\n"
        "\tthrough bars. Defaults to the standard palette.");
    symtab_get_def("ChartLabels")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ChartLabels",
        "ChartLabels\n\tBarChart option: list of label expressions drawn below\n"
        "\teach bar on the x-axis.");

    symtab_add_builtin("VectorPlot", builtin_vectorplot);
    symtab_get_def("VectorPlot")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("VectorPlot",
        "VectorPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}, opts...]\n"
        "\tDraws a grid of arrows showing the direction (and optionally\n"
        "\tmagnitude) of the vector field {vx, vy} at each grid point.\n"
        "\tVectorPlot is HoldAll: vx, vy are held unevaluated until x and y\n"
        "\tare bound to numeric values. Returns a Graphics[...] object.\n"
        "\n"
        "\tOptions:\n"
        "\t  VectorPoints   integer n → n×n grid (default 15); Automatic = 15\n"
        "\t  VectorScale    Automatic: equal-length arrows (direction only)\n"
        "\t                 None: proportional to magnitude\n"
        "\t                 real f: arrow length = f × grid spacing\n"
        "\t  VectorStyle    style directive(s) applied to all arrows\n"
        "\t  ColorFunction  named ramp string (keyed to speed) or\n"
        "\t                 f[vx,vy,speed]/f[speed]→color.\n"
        "\t                 Ramps: \"Rainbow\", \"CoolTones\", \"WarmTones\",\n"
        "\t                 \"Greyscale\", \"Temperature\"\n"
        "\t  ColorFunctionScaling  True (default): normalise speed to [0,1]\n"
        "\t  RegionFunction f[x,y] mask: skip grid points outside the region\n"
        "\t  Standard Graphics options (Axes, AspectRatio→1, Frame, PlotRange,\n"
        "\t  AxesLabel, GridLines, ImageSize, Background, PlotLabel, …) pass\n"
        "\t  through to the Graphics[...] result.");

    /* Option keyword inerts for VectorPlot. */
    symtab_get_def("VectorPoints")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VectorPoints",
        "VectorPoints\n\tVectorPlot option: integer n specifies an n×n seed grid\n"
        "\t(default 15). Automatic also uses 15.");
    symtab_get_def("VectorScale")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VectorScale",
        "VectorScale\n\tVectorPlot option: controls arrow length.\n"
        "\t  Automatic (default): all arrows drawn at equal length (direction only)\n"
        "\t  None: length proportional to field magnitude\n"
        "\t  real f: arrow length = f × grid spacing");
    symtab_get_def("VectorStyle")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("VectorStyle",
        "VectorStyle\n\tVectorPlot option: style directive(s) (RGBColor, Thickness, …)\n"
        "\tapplied globally to all arrows. Overrides per-arrow ColorFunction.");

    symtab_get_def("ScalingFunctions")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ScalingFunctions",
        "ScalingFunctions\n"
        "\tOption for Plot, ListPlot, DensityPlot, ContourPlot, VectorPlot,\n"
        "\tStreamPlot: applies a coordinate transform to one or both axes.\n"
        "\n"
        "\tForms:\n"
        "\t  ScalingFunctions -> \"Log\"         both axes: natural log\n"
        "\t  ScalingFunctions -> \"Log10\"        both axes: log base 10\n"
        "\t  ScalingFunctions -> \"Log2\"         both axes: log base 2\n"
        "\t  ScalingFunctions -> \"Reverse\"      both axes: mirror (negate)\n"
        "\t  ScalingFunctions -> {\"Log\", None}  x-axis log, y-axis linear\n"
        "\t  ScalingFunctions -> {None, \"Log\"}  x-axis linear, y-axis log\n"
        "\t  ScalingFunctions -> None           identity (default)\n"
        "\t  ScalingFunctions -> Automatic      identity (default)\n"
        "\n"
        "\tWhen a log scale is active, tick labels show original data-space\n"
        "\tvalues (e.g. 1, 2, 5, 10, 20, …) at decade-based positions.\n"
        "\tNon-positive values on a log-scaled axis are suppressed (mapped to\n"
        "\t-1e30 in world space).");

    symtab_add_builtin("StreamPlot", builtin_streamplot);
    symtab_get_def("StreamPlot")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("StreamPlot",
        "StreamPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}, opts...]\n"
        "\tTraces streamlines of the 2-D vector field {vx, vy} by RK4 integration\n"
        "\tfrom a grid of seed points, and returns a Graphics[{Arrow[...], ...}, opts]\n"
        "\tobject (auto-displayed). StreamPlot is HoldAll: vx, vy, and the iterator\n"
        "\tspecs are held unevaluated until x and y are given numeric values.\n"
        "\n"
        "\tOptions:\n"
        "\t  StreamPoints  - Integer n (n x n seed grid) or Automatic (default 15 x 15).\n"
        "\t  StreamScale   – Automatic (8%% of domain diagonal, default), None (full run),\n"
        "\t                   or a real fraction of the domain diagonal.\n"
        "\t  StreamStyle   – Style directive(s) applied to all streams.\n"
        "\t  StreamColorFunction / ColorFunction\n"
        "\t                 – f[x,y,vx,vy,speed] (or fewer args) returning a color,\n"
        "\t                   or a named ramp: \"Rainbow\", \"CoolTones\",\n"
        "\t                   \"WarmTones\", \"Greyscale\", \"Temperature\"\n"
        "\t                   (all keyed to scaled speed).\n"
        "\t  RegionFunction – f[x,y] mask; seeds outside the region are skipped.\n"
        "\t  PlotLegends   – Automatic / \"Expressions\" / explicit label list.\n"
        "\t  Standard Graphics options (PlotRange, Axes, AspectRatio, Frame, …)\n"
        "\t                   pass through to the Graphics[...] result.");

    register_inert("Arrow",
        "Arrow[{{x1,y1}, {x2,y2}, ...}]\n"
        "\tA graphics primitive: a directed polyline with an arrowhead at its\n"
        "\tlast point. Used by StreamPlot to draw streamlines.");

    symtab_add_builtin("Plot3D", builtin_plot3d);
    symtab_get_def("Plot3D")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Plot3D",
        "Plot3D[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]\n\tSamples f "
        "over a uniform grid on [xmin,xmax] x [ymin,ymax], displays the "
        "resulting surface in an interactive orbit-camera window, and "
        "returns it as a Graphics3D[...] object. A list of functions "
        "Plot3D[{f1, f2, ...}, {x,...}, {y,...}] draws each surface in a "
        "distinct palette colour. Shares Plot's option semantics where "
        "they apply: PlotPoints (per-axis grid resolution, default 25), "
        "MaxRecursion (doubles the whole grid's resolution while a "
        "flatness check fails, default 2 -- a global, crack-free analogue "
        "of Plot's adaptive bisection), Mesh (overlay the grid wireframe; "
        "default True, unlike Plot's None), PlotStyle, ColorFunction (a "
        "function of scaled-x and z, or \"Rainbow\"), ColorFunctionScaling "
        "(default True), RegionFunction (f[x,y,z], or Plot's f[x,y]/f[x] "
        "forms), PlotRange (an explicit {zmin,zmax} z-band), Axes, "
        "PlotLabel, Background, ImageSize, "
        "Lighting (Automatic (default, Lambertian shading) or None to disable).");
}
