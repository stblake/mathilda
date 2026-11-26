/* Builtin-audit regressions: graphics (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* Plot3D#1 (high): a huge PlotPoints crashed the kernel -- build_grid's
 * n*n GridPt malloc wrapped (2^32 squared is 0 mod 2^64) and was written
 * through.  The grid is checked first: past memory aborts with General::nomem. */
static void test_plot3d_huge_plot_points(void) {
    audit_chk("Plot3D#1",
              "Quiet[CheckAbort[Head[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 2^32]], aborted]]",
              "aborted");
    audit_chk("Plot3D#1",
              "Quiet[Check[CheckAbort[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 2^31], 0], "
              "fired, General::nomem]]",
              "fired");
    audit_chk("Plot3D#1", "Head[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 10]]", "Graphics3D");
}

/* Histogram#1 (high): with an explicit {min, max, dw} bin spec, data outside the bins
 * was clamped into the edge bins (15 counted in [5, 10), 4 in [3, 4)).  The bins are
 * left-closed [min + k dw, min + (k+1) dw) and data in none of them is not counted. */
static void test_histogram_explicit_range_excludes(void) {
    audit_chk("Histogram#1",
              "Cases[Histogram[{1, 5, 9, 15}, {0, 10, 5}], Rectangle[{a_, _}, {b_, h_}] :> {a, b, h}, Infinity]",
              "{{0.0, 5.0, 1.0}, {5.0, 10.0, 2.0}}");
    audit_chk("Histogram#1",
              "Cases[Histogram[{1, 2, 3, 4}, {0, 4, 1}], Rectangle[{a_, _}, {b_, h_}] :> {a, b, h}, Infinity]",
              "{{1.0, 2.0, 1.0}, {2.0, 3.0, 1.0}, {3.0, 4.0, 1.0}}");
    audit_chk("Histogram#1",
              "Cases[Histogram[{-3, 0, 2, 2.5, 7}, {0, 5, 2.5}], Rectangle[{a_, _}, {b_, h_}] :> {a, b, h}, Infinity]",
              "{{0.0, 2.5, 2.0}, {2.5, 5.0, 1.0}}");
    /* the counted total is the number of data in [min, max) */
    audit_true("Histogram#1",
               "SeedRandom[9]; Module[{d = RandomReal[{-1, 2}, 500]},"
               " Total[Cases[Histogram[d, {0, 1, 0.1}], Rectangle[_, {_, h_}] :> h, Infinity]] =="
               " Count[d, v_ /; 0 <= v < 1]]");
    /* sanity: a derived range still counts every datum, the maximum included */
    audit_true("Histogram#1",
               "Total[Cases[Histogram[{1, 2, 3, 4}], Rectangle[_, {_, h_}] :> h, Infinity]] == 4");
    audit_true("Histogram#1",
               "Total[Cases[Histogram[{1, 2, 3, 4, 9}, 3], Rectangle[_, {_, h_}] :> h, Infinity]] == 5");
}

/* ================================================================== */
/* Wave 4 G3: graphics                                                */
/* ================================================================== */

/* A body that evaluates to a list plots each element, as the literal list
 * would (plot_redispatch_list_body, src/graphics/plot_common.c). */
static void test_wave4_list_valued_bodies(void) {
    audit_chk("ParametricPlot#1", "ParametricPlot[{2{Cos[t],Sin[t]},{3Cos[t],Sin[t]}},{t,0,2Pi}] // Head", "Graphics");
    audit_chk("ParametricPlot#2", "f={{Cos[t],Sin[t]},{2Cos[t],Sin[t]}}; Head[ParametricPlot[f,{t,0,1}]]", "Graphics");
    audit_chk("Plot#1", "Head[Plot[Table[x^n, {n, 3}], {x, 0, 1}]]", "Graphics");
    audit_true("Plot#1", "Length[Cases[Plot[Table[x^n, {n, 3}], {x, 0, 1}], _Line, Infinity]] == 3");
    audit_true("Plot#1", "Head[Plot[Table[x^n, {n, 3}], {x, 0, 1}][[1]]] === List");
    audit_chk("ParametricPlot3D#1", "L={{Cos[t],Sin[t],0},{Cos[t],0,Sin[t]}}; Head[ParametricPlot3D[L,{t,0,1}]]", "Graphics3D");
    audit_chk("PolarPlot#1", "f = {1, 2}; Head[PolarPlot[f, {t, 0, 2 Pi}]]", "Graphics");
    audit_chk("StreamPlot#1", "v={y,-x}; Head[TimeConstrained[StreamPlot[v,{x,-3,3},{y,-3,3}],10]]", "Graphics");
    audit_chk("StreamPlot#1", "ff[x_,y_]:={y,-x}; Head[TimeConstrained[StreamPlot[ff[x,y],{x,-3,3},{y,-3,3}],10]]", "Graphics");
    audit_chk("VectorPlot#2", "v = {-y, x}; Head[VectorPlot[v, {x, -1, 1}, {y, -1, 1}, VectorPoints -> 2]]", "Graphics");
    /* a global value of the plot variable is irrelevant */
    audit_chk("Plot#1", "x = 3; r = Head[Plot[x^2, {x, 0, 1}]]; Clear[x]; r", "Graphics");
    audit_chk("ParametricPlot#3",
              "Head[ParametricPlot[{Cos[t],Sin[t]},{t,0,2Pi},RegionFunction->Function[{x,y,u},u<1]]]", "Graphics");
    audit_true("ParametricPlot#3",
               "Max[Cases[ParametricPlot[{t, t}, {t, 0, 2}, RegionFunction -> Function[{x, y, u}, u < 1]],"
               " Line[p_] :> p[[All, 1]], Infinity]] <= 1");
}

/* An explicit Frame (False included) suppresses the injected Frame -> True. */
static void test_wave4_frame_option(void) {
    audit_true("ArrayPlot#1", "Count[ArrayPlot[{{1,0},{0,1}}, Frame -> False], HoldPattern[Frame -> _]] == 1");
    audit_true("ArrayPlot#1", "MemberQ[ArrayPlot[{{1,0},{0,1}}, Frame -> False], Frame -> False]");
    audit_true("ComplexPlot#1",
               "Count[ComplexPlot[z, {z, -1-I, 1+I}, PlotPoints -> 2, Frame -> False], HoldPattern[Frame -> _]] == 1");
}

/* PlotPoints -> {nx, ny} and Automatic; Plot3D Mesh forms. */
static void test_wave4_plot_points_pair(void) {
    audit_true("ComplexPlot#2", "Length[Cases[ComplexPlot[z, {z, -1-I, 1+I}, PlotPoints -> {2,3}], _Rectangle, Infinity]] == 6");
    audit_chk("ContourPlot#1", "ContourPlot[x+y,{x,0,1},{y,0,1},PlotPoints->{10,20}]//Head", "Graphics");
    audit_chk("DensityPlot#1", "g=DensityPlot[x+y,{x,0,1},{y,0,1},PlotPoints->{2,3}]; Head[g]", "Graphics");
    audit_true("DensityPlot#1", "Length[Cases[DensityPlot[x+y,{x,0,1},{y,0,1},PlotPoints->{2,3}], _Rectangle, Infinity]] == 6");
    audit_chk("Plot3D#2", "Head[Plot3D[x+y,{x,0,1},{y,0,1},PlotPoints->Automatic]]", "Graphics3D");
    audit_chk("Plot3D#2", "Head[Plot3D[x+y,{x,0,1},{y,0,1},MaxRecursion->Automatic]]", "Graphics3D");
    audit_chk("Plot3D#3", "Head[Plot3D[x+y,{x,0,1},{y,0,1},PlotPoints->{5,10}]]", "Graphics3D");
    audit_chk("Plot3D#4", "Head[Plot3D[x+y,{x,0,1},{y,0,1},Mesh->5]]", "Graphics3D");
    audit_chk("Plot3D#4", "Head[Plot3D[x+y,{x,0,1},{y,0,1},Mesh->{3,4}]]", "Graphics3D");
    audit_chk("Plot3D#4", "Head[Plot3D[x+y,{x,0,1},{y,0,1},Mesh->Full]]", "Graphics3D");
    audit_true("ComplexPlot3D#1",
               "g = ComplexPlot3D[z^2, {z, -2-2I, 2+2I}, PlotPoints -> 20, PlotRange -> All];"
               " Max[Cases[g, Polygon[p_] :> p[[All, 3]], Infinity]] == 8.");
}

static void test_wave4_bar_histogram(void) {
    audit_true("BarChart#1", "Cases[BarChart[{1, x, 3}], Rectangle[__], Infinity] === "
                             "{Rectangle[{0.6, 0.}, {1.4, 1.}], Rectangle[{2.6, 0.}, {3.4, 3.}]}");
    audit_true("BarChart#1", "Cases[BarChart[{1, Missing[], 3}], Rectangle[__], Infinity] === "
                             "{Rectangle[{0.6, 0.}, {1.4, 1.}], Rectangle[{2.6, 0.}, {3.4, 3.}]}");
    audit_true("BarChart#2", "Cases[BarChart[{1, 2, 3}, BarSpacing -> 1], Rectangle[{a_, _}, {b_, _}] :> b - a, Infinity] == {0.5, 0.5, 0.5}");
    audit_true("BarChart#2", "And @@ Cases[BarChart[{1, 2}, BarSpacing -> 2], Rectangle[{a_, _}, {b_, _}] :> b > a, Infinity]");
    audit_true("BarChart#3", "Cases[BarChart[{{1,3,2},{4,2,5}}], Rectangle[_, {_, h_}] :> h, Infinity] == {1., 3., 2., 4., 2., 5.}");
    audit_true("BarChart#3", "Cases[BarChart[{{1,3,2},{4,2,5}}], Rectangle[{a_, _}, _] :> a, Infinity][[4]] > 1.5");
    audit_true("BarChart#1", "Cases[BarChart[{1, 2, 3}], Rectangle[__], Infinity] === "
                             "{Rectangle[{0.6, 0.}, {1.4, 1.}], Rectangle[{1.6, 0.}, {2.4, 2.}], Rectangle[{2.6, 0.}, {3.4, 3.}]}");
    audit_true("Histogram#2", "Cases[Histogram[{0.5, 4, 9.5}, {0, 10, 3}], Rectangle[{a_, _}, {b_, h_}] :> {a, b, h}, Infinity] == "
                              "{{0., 3., 1.}, {3., 6., 1.}, {9., 12., 1.}}");
    audit_true("Histogram#2", "Cases[Histogram[{1.5, 2.5, 7.2}, {2}], Rectangle[{a_, _}, {b_, h_}] :> {a, b, h}, Infinity] == "
                              "{{0., 2., 1.}, {2., 4., 1.}, {6., 8., 1.}}");
    audit_true("Histogram#3", "Cases[Histogram[{1/2, 3/2, Pi}, {0, 4, 1}], Rectangle[{a_, _}, {b_, h_}] :> {a, b, h}, Infinity] == "
                              "{{0., 1., 1.}, {1., 2., 1.}, {3., 4., 1.}}");
    audit_true("Histogram#3", "Cases[Histogram[{1, 2, E}, {0, 4, 1}], Rectangle[{a_, _}, {b_, h_}] :> {a, b, h}, Infinity] == "
                              "{{1., 2., 1.}, {2., 3., 2.}}");
}

static void test_wave4_listplot_missing_pair(void) {
    audit_true("ListPlot#1", "Cases[ListPlot[{{1,2},{3,x}}], Point[p_] :> p, Infinity] == {{{1., 2.}}}");
    audit_true("ListPlot#1", "Cases[ListPlot[{{1,2},{3,4}}], Point[p_] :> p, Infinity] == {{{1., 2.}, {3., 4.}}}");
}

static void test_wave4_rgbcolor_list(void) {
    audit_true("RGBColor#1", "RGBColor[{0.2, 0.4, 0.6, 0.5}] === RGBColor[0.2, 0.4, 0.6, 0.5]");
    audit_true("RGBColor#1", "RGBColor[{1, 0, 0}] === RGBColor[1, 0, 0] === Red");
}

/* PDF export: Arrow[{line1, line2}], Opacity[a, color], Text of any expression. */
static void test_wave4_export_primitives(void) {
    audit_run("g3f = FileNameJoin[{$TemporaryDirectory, \"audit_g3_wave4.pdf\"}]");
    audit_true("Arrow#1", "Export[g3f, Graphics[{Arrow[{{{0,0},{1,1}},{{1,0},{0,1}}}]}]]; StringCount[ReadString[g3f], \" h f\"] == 2");
    audit_true("Opacity#1", "Export[g3f, Graphics[{Opacity[0.5, Red], Disk[]}]]; StringContainsQ[ReadString[g3f], \"1.0000 0.0000 0.0000 rg\"]");
    audit_true("Text#1", "Export[g3f, Graphics[{Text[1/2, {0,0}], Text[10^30, {1,1}], Text[7, {2,2}]}]];"
                         " SubsetQ[StringCases[ReadString[g3f], RegularExpression[\"\\\\(([^)]*)\\\\) Tj\"] -> \"$1\"],"
                         " {\"1/2\", \"1000000000000000000000000000000\", \"7\"}]");
    audit_run("DeleteFile[g3f]");
}

/* Show's union reads PlotRange -> r, partial ranges, GraphicsComplex, BezierCurve, Sphere. */
static void test_wave4_show_ranges(void) {
    audit_true("Show#1", "Cases[Show[Graphics[Point[{0,0}], PlotRange->3], Graphics[Point[{1,1}]]], HoldPattern[PlotRange -> r_] :> r, Infinity] == "
                         "{{{-3., 3.}, {-3., 3.}}}");
    audit_true("Show#1", "Length[Cases[Show[Graphics[Point[{0,0}], PlotRange->{{0,1},All}], Graphics[Point[{5,5}]]], HoldPattern[PlotRange -> _], Infinity]] == 1");
    audit_true("Show#2", "Cases[Show[Graphics[Point[{0,0}], PlotRange->{{0,1},{0,1}}], Graphics[GraphicsComplex[{{5,5},{6,6}}, Line[{1,2}]]]], "
                         "HoldPattern[PlotRange -> r_] :> r, Infinity] == {{{0., 6.}, {0., 6.}}}");
    audit_true("Show#2", "Cases[Show[Graphics3D[Sphere[{0,0,0}], PlotRange->{{0,1},{0,1},{0,1}}], Graphics3D[Sphere[{5,5,5}]]], "
                         "HoldPattern[PlotRange -> r_] :> r, Infinity] == {{{0., 6.}, {0., 6.}, {0., 6.}}}");
    audit_true("Show#2", "Cases[Show[Graphics[Point[{0,0}], PlotRange->{{0,1},{0,1}}], Graphics[BezierCurve[{{0,0},{2,3},{4,0}}]]], "
                         "HoldPattern[PlotRange -> r_] :> r, Infinity] == {{{0., 4.}, {0., 3.}}}");
}

/* VectorPlot's ColorFunction gets f[x, y, vx, vy, n], scaled to [0, 1]. */
static void test_wave4_vectorplot_colorfunction(void) {
    audit_true("VectorPlot#1", "Cases[VectorPlot[{1, 0}, {x, 0, 1}, {y, 0, 1}, VectorPoints -> 2,"
                               " ColorFunction -> Function[{x, y, vx, vy, n}, Hue[x]]][[1]], _Hue] == "
                               "{Hue[0.], Hue[1.], Hue[0.], Hue[1.]}");
    audit_true("VectorPlot#1", "Head[VectorPlot[{1, 0}, {x, 0, 1}, {y, 0, 1}, VectorPoints -> 2,"
                               " ColorFunction -> Function[{s}, Hue[s]]]] === Graphics");
}

static void test_wave4_complexplot3d_list(void) {
    audit_true("ComplexPlot3D#2", "g = ComplexPlot3D[{z, z^2}, {z, -1-I, 1+I}, PlotPoints -> 4]; Head[g] === Graphics3D && "
                                  "Length[Cases[g, _Polygon, Infinity]] == 2 Length[Cases[ComplexPlot3D[z, {z, -1-I, 1+I}, PlotPoints -> 4], _Polygon, Infinity]]");
}

/* ================================================================== */
/* Wave 5 G3: plot range messages, empty plots, option defaults, ...   */
/* ================================================================== */

/* Range specs: pllim / plln / plld / nonopt, each leaving the call alone;
 * src/graphics/plot_common.c plot_range_bounds, src/arity_table.c. */
static void test_wave5_plot_range_messages(void) {
    audit_chk("Plot#5", "Check[Plot[x, {x, 0, 1, 0.5}], \"fired\", Plot::pllim]", "\"fired\"");
    audit_chk("Plot#6", "Check[Plot[x, {x, 0, 0}], \"fired\", Plot::plld]", "\"fired\"");
    audit_chk("Plot#6", "Check[Plot[x, {x, a, b}], \"fired\", Plot::plln]", "\"fired\"");
    audit_chk("Plot#6", "Check[Plot[x, {x, 0, 1}, foo], \"fired\", Plot::nonopt]", "\"fired\"");
    audit_chk("Plot#6", "Quiet[Plot[x, {x, 0, 0}]]", "Plot[x, {x, 0, 0}]");
    audit_chk("ParametricPlot#5", "Check[ParametricPlot[{Cos[t],Sin[t]},{t,0,0}], \"fired\", ParametricPlot::plld]", "\"fired\"");
    audit_chk("ParametricPlot#6", "Check[ParametricPlot[{Cos[t],Sin[t]},{t,0,1},foo], \"fired\", ParametricPlot::nonopt]", "\"fired\"");
    audit_chk("Plot3D#6", "Check[Plot3D[x+y,{x,1,1},{y,0,1}], \"fired\", Plot3D::plld]", "\"fired\"");
    audit_chk("Plot3D#6", "Quiet[Plot3D[x+y,{x,1,1},{y,0,1}]]", "Plot3D[x + y, {x, 1, 1}, {y, 0, 1}]");
}

/* Nothing plottable is an empty Graphics / Graphics3D, not an unevaluated call. */
static void test_wave5_empty_plots(void) {
    audit_chk("Plot#2", "Plot[I x, {x, 0, 1}][[1]]", "{}");
    audit_chk("Plot#2", "Head[Plot[Indeterminate, {x, 0, 1}]]", "Graphics");
    audit_chk("Plot#3", "Plot[{}, {x,0,1}][[1]]", "{}");
    audit_chk("ParametricPlot#4", "Head[ParametricPlot[{x,y},{t,0,1}]]", "Graphics");
    audit_chk("ParametricPlot3D#2", "ParametricPlot3D[{x,y,z}, {t, 0, 1}][[1]]", "{}");
    audit_chk("Plot3D#5", "Head[Plot3D[f,{x,0,1},{y,0,1}]]", "Graphics3D");
    audit_chk("Plot3D#5", "Plot3D[Sqrt[-1-x^2],{x,0,1},{y,0,1}][[1]]", "{}");
    audit_chk("Plot3D#5", "Plot3D[{},{x,0,1},{y,0,1}][[1]]", "{}");
    audit_chk("StreamPlot#3", "Head[TimeConstrained[StreamPlot[{0,0},{x,-1,1},{y,-1,1}],10]]", "Graphics");
    audit_true("StreamPlot#3", "Cases[StreamPlot[{0,0},{x,-1,1},{y,-1,1}], HoldPattern[PlotRange -> r_] :> r, Infinity] == {{{-1., 1.}, {-1., 1.}}}");
}

/* Exclusions -> Automatic splits at Floor / UnitStep jumps. */
static void test_wave5_plot_exclusions(void) {
    audit_chk("Plot#4", "Length[Cases[Plot[Floor[x], {x, 0, 3}, Exclusions -> Automatic], _Line, Infinity]]", "3");
    audit_chk("Plot#4", "Length[Cases[Plot[UnitStep[x - 1/2], {x, 0, 1}], _Line, Infinity]]", "2");
    audit_chk("Plot#4", "Length[Cases[Plot[Floor[x], {x, 0, 3}, Exclusions -> None], _Line, Infinity]]", "1");
    audit_true("Plot#4", "Max[Cases[Plot[Floor[x], {x, 0, 3}], Line[p_] :> Max[Abs[Differences[p[[All, 2]]]]], Infinity]] == 0");
}

/* Options[...] report the defaults the code uses. */
static void test_wave5_option_defaults(void) {
    audit_chk("Plot#7", "Options[Plot, {AspectRatio, MaxPlotPoints, Mesh}]",
              "{AspectRatio -> 1/GoldenRatio, MaxPlotPoints -> Infinity, Mesh -> None}");
    audit_chk("Graphics#1", "Options[Graphics, {PlotRange, Background, AxesStyle, FrameStyle, TicksStyle, GridLinesStyle}]",
              "{PlotRange -> All, Background -> None, AxesStyle -> {}, FrameStyle -> {}, TicksStyle -> {}, GridLinesStyle -> {}}");
    audit_chk("Graphics3D#1", "Options[Graphics3D, Boxed]", "{Boxed -> True}");
    audit_chk("Graphics3D#1", "Options[Graphics3D, ViewPoint]", "{ViewPoint -> {1.3, -2.4, 2.0}}");
    audit_true("ParametricPlot3D#4", "Length[Options[ParametricPlot3D]] > 5 && Options[ParametricPlot3D, PlotPoints] === {PlotPoints -> 25}");
    audit_true("PolarPlot#3", "Length[Options[PolarPlot]] > 5 && Options[PolarPlot, PlotPoints] === {PlotPoints -> 75}");
}

/* Docstrings that state the real behaviour. */
static void test_wave5_docstrings(void) {
    audit_true("ParametricPlot#9", "!StringContainsQ[Information[ParametricPlot], \"PlotPoints (default 25)\"] && "
                                   "StringContainsQ[Information[ParametricPlot], \"default 50 for curves, 75\"]");
    audit_true("ParametricPlot#9", "Length[Cases[ParametricPlot[{Cos[t],Sin[t]},{t,0,2Pi},MaxRecursion->0],Line[p_]:>p,Infinity][[1]]] == 50");
    audit_true("ContourPlot#3", "!StringContainsQ[Information[ContourPlot], \"default 25\"] && StringContainsQ[Information[ContourPlot], \"default 75\"]");
    audit_true("GrayLevel#1", "StringContainsQ[Information[GrayLevel], \"GrayLevel[g, a]\"]");
    audit_true("Opacity#2", "StringContainsQ[Information[Opacity], \"Opacity[a, color]\"]");
    audit_true("StreamPlot#4", "With[{s = Information[StreamPlot]}, !StringContainsQ[s, \"Arrow[\"] && !StringContainsQ[s, \"%%\"] && "
                               "!StringContainsQ[s, \"15 x 15\"] && StringContainsQ[s, \"default 25\"]]");
    audit_true("Histogram#6", "!StringContainsQ[Information[Histogram], \"BarSpacing\"]");
    audit_chk("Histogram#6", "Cases[Histogram[{1, 2, 3}, 2, BarSpacing -> 0.5], Rectangle[{a_, _}, {b_, h_}] :> {a, b, h}, Infinity]",
              "{{1.0, 2.0, 1.0}, {2.0, 3.0, 2.0}}");
    audit_true("StreamPlot#5", "g = StreamPlot[{y,-x},{x,-3,3},{y,-3,3}]; {Length[Cases[g, _Arrow, Infinity]], Options[g, Frame], Options[g, Axes]} === {0, {Frame -> True}, {Axes -> False}}");
    audit_chk("ContourPlot#4", "Rest[List @@ ContourPlot[x+y,{x,0,1},{y,0,1}]]",
              "{Frame -> True, Axes -> False, AspectRatio -> 1, PlotRange -> {{0.0, 1.0}, {0.0, 1.0}}}");
    audit_true("VectorPlot#6", "First[Cases[VectorPlot[{x, 0}, {x, -1, 1}, {y, -1, 1}, VectorPoints -> 3, VectorScale -> None][[1]], _RGBColor]] === RGBColor[0.993248, 0.906157, 0.143936]");
    audit_chk("ComplexPlot#4", "Cases[ComplexPlot[z, {z, -1-I, 1+I}, PlotPoints -> 2], HoldPattern[(Frame | Axes) -> _], Infinity]",
              "{Frame -> True, Axes -> False}");
}

/* ParametricPlot: Mesh -> n, ColorFunction f[x, y, u]; muted singular samples. */
static void test_wave5_parametric(void) {
    audit_chk("ParametricPlot#7", "Cases[ParametricPlot[{Cos[t],Sin[t]},{t,0,2Pi},Mesh->5], Point[p_] :> Length[p], Infinity]", "{5}");
    audit_chk("ParametricPlot#8", "Cases[ParametricPlot[{2 t,0},{t,0,1},ColorFunction->(Hue[#]&), ColorFunctionScaling->False, "
                                  "MaxRecursion->0, PlotPoints->3],_Hue,Infinity]", "{Hue[0.5], Hue[1.5]}");
    audit_chk("ParametricPlot3D#3", "Check[Head[ParametricPlot3D[{t,1/t,t},{t,-1,1}]], \"fired\"]", "Graphics3D");
    audit_chk("PolarPlot#2", "Check[Cases[PolarPlot[1/t, {t, 0, 2 Pi}], Line[p_] :> Length[p], Infinity], \"fired\"]", "{159}");
}

/* ComplexPlot / ComplexPlot3D defaults and list options. */
static void test_wave5_complexplot(void) {
    audit_chk("ComplexPlot#3", "Cases[ComplexPlot[z, {z, -1-I, 1+I}, PlotPoints -> 2, Axes -> True], HoldPattern[(Frame | Axes) -> _], Infinity]",
              "{Axes -> True, Frame -> True}");
    audit_chk("ComplexPlot#5", "Head[ComplexPlot[z, {z, -1-I, 1+I}, {PlotPoints -> 2}]]", "Graphics");
    audit_chk("ComplexPlot3D#3", "Rest[List @@ ComplexPlot3D[z, {z, -1-I, 1+I}, PlotPoints -> 2, Lighting -> None]][[;; 2]]",
              "{Lighting -> None, Axes -> True}");
    audit_true("ComplexPlot3D#3", "FreeQ[ComplexPlot3D[z, {z, -1-I, 1+I}, PlotPoints -> 2], Frame | AspectRatio]");
}

/* ContourPlot, DensityPlot, VectorPlot option handling. */
static void test_wave5_field_plots(void) {
    audit_chk("ContourPlot#2", "Length[Cases[ContourPlot[x+y,{x,0,1},{y,0,1},Contours->{}],_Line,Infinity]]", "0");
    audit_chk("DensityPlot#2", "CheckAbort[DensityPlot[x,{x,0,1},{y,0,1},PlotPoints->2^32+2], \"aborted\"]", "\"aborted\"");
    audit_chk("DensityPlot#3", "Cases[DensityPlot[z,{x,0,1},{y,0,1},PlotPoints->2,PlotLegends->Automatic], _$StreamColorBar, Infinity]",
              "{$StreamColorBar[0.0, 1.0, Automatic]}");
    audit_true("DensityPlot#4", "Min[Cases[DensityPlot[x,{x,-1,1},{y,-1,1},PlotPoints->2,ScalingFunctions->\"Log\"], Rectangle[a_, b_] :> Min[a, b], Infinity]] > -10");
    audit_chk("VectorPlot#3", "VectorPlot[{1, 0}, {x, -1, 1}, {y, -1, 1}, VectorPoints -> 3, VectorScale -> 0.5][[1, 3]]",
              "Arrow[{{-1.25, -1.0}, {-0.75, -1.0}}]");
    audit_chk("VectorPlot#4", "Length[Cases[VectorPlot[{-y, x}, {x, -1, 1}, {y, -1, 1}, VectorPoints -> {2, 3}], _Arrow, Infinity]]", "6");
    audit_true("VectorPlot#5", "Min[Cases[VectorPlot[{1, 1}, {x, 0, 10}, {y, 0, 1}, VectorPoints -> 2, ScalingFunctions -> \"Log\"], Arrow[p_] :> Min[p], Infinity]] > -20");
}

/* StreamPlot of several fields. */
static void test_wave5_streamplot_fields(void) {
    audit_true("StreamPlot#2", "g = TimeConstrained[StreamPlot[{{y,-x},{1,1}},{x,-3,3},{y,-3,3}],30]; "
                               "Head[g] === Graphics && Length[Cases[g, _Line, Infinity]] > Length[Cases[StreamPlot[{y,-x},{x,-3,3},{y,-3,3}], _Line, Infinity]]");
}

/* BarChart / Histogram. */
static void test_wave5_charts(void) {
    audit_chk("BarChart#4", "{Head[BarChart[{1, 2}, BarSpacing -> None]], Head[BarChart[{1, 2}, BarSpacing -> Automatic]]}", "{Graphics, Graphics}");
    audit_chk("BarChart#4", "Cases[BarChart[{1, 2}, BarSpacing -> None], Rectangle[__], Infinity]",
              "{Rectangle[{0.5, 0.0}, {1.5, 1.0}], Rectangle[{1.5, 0.0}, {2.5, 2.0}]}");
    audit_chk("BarChart#5", "Cases[BarChart[{1, {2, 3}}], Rectangle[_, {_, h_}] :> h, Infinity]", "{1.0, 2.0, 3.0}");
    audit_chk("BarChart#5", "Cases[BarChart[{{1, 2}, 3}], Rectangle[_, {_, h_}] :> h, Infinity]", "{1.0, 2.0, 3.0}");
    audit_true("BarChart#7", "Cases[BarChart[{1, 2}, ChartStyle -> \"Pastel\"], _String, Infinity] === {} && "
                             "Length[Cases[BarChart[{1, 2}, ChartStyle -> \"Rainbow\"], _RGBColor, Infinity]] == 2");
    audit_chk("BarChart#8", "Head[BarChart[<|\"a\" -> 1, \"b\" -> 2|>]]", "Graphics");
    audit_chk("BarChart#8", "Cases[BarChart[<|\"a\" -> 1, \"b\" -> 2|>], Rectangle[_, {_, h_}] :> h, Infinity]", "{1.0, 2.0}");
    audit_true("BarChart#8", "(List @@ First[Cases[BarChart[<|\"a\" -> 1, \"b\" -> 2|>], _$BarChartLabels, Infinity]])[[All, 2]] === {\"a\", \"b\"}");
    audit_chk("Histogram#4", "CheckAbort[Histogram[{1,2,3,4}, 2^32+2], \"aborted\"]", "\"aborted\"");
    audit_chk("Histogram#5", "Length[Union[Cases[Histogram[{1,2,3,4,5,6,7,8}], _RGBColor, Infinity]]]", "1");
}

/* ListPlot: messages and the Joined gap. */
static void test_wave5_listplot(void) {
    audit_chk("ListPlot#2", "Cases[ListPlot[{1,x,3,4},Joined->True], _Line, Infinity]", "{Line[{{3.0, 3.0}, {4.0, 4.0}}]}");
    audit_chk("ListPlot#2", "Length[Cases[ListPlot[{1,2,x,4,5},Joined->True], _Line, Infinity]]", "2");
    audit_chk("ListPlot#3", "Check[ListPlot[{1,2},5], \"fired\", ListPlot::nonopt]", "\"fired\"");
    audit_chk("ListPlot#3", "Check[ListPlot[{1,2,3},DataRange->{1,2,3}], \"fired\", ListPlot::drng]", "\"fired\"");
    audit_chk("ListPlot#3", "Check[ListPlot[{1,2},Foo->3], \"fired\", ListPlot::optx]", "\"fired\"");
    audit_chk("ListPlot#3", "Quiet[ListPlot[{1,2},Foo->3]]", "ListPlot[{1, 2}, Foo -> 3]");
    audit_chk("ListPlot#3", "Head[ListPlot[{1,2}, Joined -> True, ImageSize -> 300, AxesLabel -> {\"a\", \"b\"}]]", "Graphics");
}

/* Show: elliptical radii, Rectangle[{x, y}], gtype. */
static void test_wave5_show(void) {
    audit_true("Show#3", "Cases[Show[Graphics[Point[{0,0}], PlotRange->{{0,1},{0,1}}], Graphics[Disk[{0,0},{5,2}]]], "
                         "HoldPattern[PlotRange -> r_] :> r, Infinity] == {{{-5., 5.}, {-2., 2.}}}");
    audit_true("Show#4", "Cases[Show[Graphics[Point[{0,0}], PlotRange->{{0,1},{0,1}}], Graphics[Rectangle[{2,2}]]], "
                         "HoldPattern[PlotRange -> r_] :> r, Infinity] == {{{0., 3.}, {0., 3.}}}");
    audit_chk("Show#5", "Check[Show[Graphics[Point[{0,0}]], 1], \"fired\", Show::gtype]", "\"fired\"");
}

/* Styles and primitives: RGBColor hex, ArrayPlot cells, the PDF writer. */
static void test_wave5_styles_export(void) {
    audit_chk("RGBColor#2", "RGBColor[\"#FF0000\"]", "RGBColor[1.0, 0.0, 0.0]");
    audit_chk("RGBColor#2", "RGBColor[\"#00ff0080\"] === RGBColor[0., 1., 0., 128/255.]", "True");
    audit_chk("ArrayPlot#2", "Cases[ArrayPlot[{{RGBColor[1,0,0,0.5], RGBColor[0,0,1,0.5]}}], _Rectangle, Infinity]",
              "{Rectangle[{0.0, 0.0}, {1.0, 1.0}], Rectangle[{1.0, 0.0}, {2.0, 1.0}]}");
    audit_run("g3w5 = FileNameJoin[{$TemporaryDirectory, \"audit_g3_wave5.pdf\"}]");
    audit_run("g3w5b = FileNameJoin[{$TemporaryDirectory, \"audit_g3_wave5b.pdf\"}]");
    audit_true("AbsoluteDashing#1", "Export[g3w5, Graphics[{AbsoluteDashing[{1,2,3,4,5,6,7,8,9,10}], Line[{{0,0},{1,1}}]}]]; "
                                    "StringContainsQ[ReadString[g3w5], \"[1.000 2.000 3.000 4.000 5.000 6.000 7.000 8.000 9.000 10.000] 0 d\"]");
    audit_true("Arrow#2", "Export[g3w5, Graphics[{Arrow[{{0,0},{1,1}}]}]]; Export[g3w5b, Graphics[{Arrow[{{0,0},{1,1}}, 0.3]}]]; "
                          "ReadString[g3w5] =!= ReadString[g3w5b]");
    audit_true("Arrow#3", "Export[g3w5, Graphics[{Arrow[BezierCurve[{{0,0},{1,1},{2,0}}]]}]]; StringCount[ReadString[g3w5], \" h f\"] == 1");
    audit_true("Directive#1", "Export[g3w5, Plot[{x, 2x}, {x,0,1}, PlotStyle -> {Directive[Red, Thickness[0.02]], Directive[Blue]}]]; "
                              "StringContainsQ[ReadString[g3w5], \"0.0000 0.0000 1.0000 RG\\n1.000 w\"]");
    audit_run("DeleteFile[{g3w5, g3w5b}]");
}

/* Manipulate / Animate controls (built without a window:
 * MATHILDA_NO_GRAPHICS_WINDOW, as tests/test_manipulate.c does). */
static void test_wave5_controls(void) {
    setenv("MATHILDA_NO_GRAPHICS_WINDOW", "1", 1);
    audit_chk("Manipulate#2", "Manipulate[x, {x, Range[33]}]", "Null");
    audit_chk("Manipulate#3", "Manipulate[x^2, {{x, 1, \"lbl\"}, 0, 1}]", "Null");
    audit_chk("Animate#3", "Animate[t, {t, {1, 2, 3}}]", "Null");
    audit_chk("Manipulate#2", "Manipulate[x, {x, {}}]", "Manipulate[x, {x, {}}]");
}

/* Wave 6 G4: Animate's stepped iterator {t, tmin, tmax, dt} is the discrete
 * value list Range[tmin, tmax, dt] (built without a window). */
static void test_wave6_g4_graphics(void) {
    audit_chk("Animate#2", "Animate[t, {t, 0, 1, 0.1}]", "Null");
    audit_chk("Animate#2", "Animate[t, {t, 0, 1, 1/4}]", "Null");
    audit_chk("Animate#2", "Animate[t, {t, 0, 1, -0.1}]", "Animate[t, {t, 0, 1, -0.1}]");
}

int main(void) {
    audit_init();
    TEST(test_plot3d_huge_plot_points);
    TEST(test_histogram_explicit_range_excludes);
    TEST(test_wave4_list_valued_bodies);
    TEST(test_wave4_frame_option);
    TEST(test_wave4_plot_points_pair);
    TEST(test_wave4_bar_histogram);
    TEST(test_wave4_listplot_missing_pair);
    TEST(test_wave4_rgbcolor_list);
    TEST(test_wave4_export_primitives);
    TEST(test_wave4_show_ranges);
    TEST(test_wave4_vectorplot_colorfunction);
    TEST(test_wave4_complexplot3d_list);
    TEST(test_wave5_plot_range_messages);
    TEST(test_wave5_empty_plots);
    TEST(test_wave5_plot_exclusions);
    TEST(test_wave5_option_defaults);
    TEST(test_wave5_docstrings);
    TEST(test_wave5_parametric);
    TEST(test_wave5_complexplot);
    TEST(test_wave5_field_plots);
    TEST(test_wave5_streamplot_fields);
    TEST(test_wave5_charts);
    TEST(test_wave5_listplot);
    TEST(test_wave5_show);
    TEST(test_wave5_styles_export);
    TEST(test_wave5_controls);
    TEST(test_wave6_g4_graphics);
    return audit_finish("audit_graphics");
}
