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
    return audit_finish("audit_graphics");
}
