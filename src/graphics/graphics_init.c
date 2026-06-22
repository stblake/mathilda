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

/* Register a symbol as an inert, protected structural head with a
 * docstring but no builtin function -- mirrors how Rule/List/etc. are
 * "just data" in this codebase. */
static void register_inert(const char* name, const char* doc) {
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(name, doc);
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
        "radius r centered at {x,y}.");
    register_inert("Disk",
        "Disk[{x,y}, r]\n\tA graphics primitive: a filled disk of radius r "
        "centered at {x,y}.");
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
    register_inert("Opacity",
        "Opacity[a]\n\tA style directive: sets the opacity (a in [0,1]) of "
        "subsequent primitives.");
    register_inert("Thickness",
        "Thickness[t]\n\tA style directive: sets the line thickness (in "
        "plot coordinates) of subsequent Line/Circle primitives.");
    register_inert("PointSize",
        "PointSize[s]\n\tA style directive: sets the radius (in plot "
        "coordinates) of subsequent Point primitives.");

    register_inert("Graphics",
        "Graphics[primitives, opts...]\n\tA symbolic 2D graphics object. "
        "primitives is a (possibly nested) list of graphics primitives and "
        "style directives. Rendered on demand by Show[]. Prints as "
        "-Graphics-.");

    symtab_add_builtin("Show", builtin_show);
    symtab_get_def("Show")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Show",
        "Show[graphics, opts...]\n\tDisplays graphics (a Graphics[...] "
        "object) in an interactive window and returns it, merging any "
        "given options into its option list.");

    symtab_add_builtin("Plot", builtin_plot);
    symtab_get_def("Plot")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;
    symtab_set_docstring("Plot",
        "Plot[f, {x, xmin, xmax}, opts...]\n\tAdaptively samples f over "
        "[xmin, xmax], displays the resulting curve in an interactive "
        "window, and returns it as a Graphics[...] object. A list of "
        "functions Plot[{f1, f2, ...}, {x, xmin, xmax}] draws each on the "
        "same axes in a distinct palette colour. Options: "
        "PlotPoints (initial sample count, default 25), MaxRecursion "
        "(adaptive refinement depth, default 6), MaxPlotPoints (overall "
        "point cap, default Infinity), PlotRange, AspectRatio, PlotStyle, "
        "Axes, AxesLabel, PlotLabel, Background, ImageSize.");
}
