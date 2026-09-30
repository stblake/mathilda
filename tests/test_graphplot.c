/* test_graphplot.c - GraphPlot / HypergraphPlot and the layout engine
 * (src/graph/graphplot.c, hyp_plot.c, glayout.c).
 *
 * Determinism (the same call gives an identical Graphics expression, and the
 * layout engine gives bit-identical coordinates), every option (arrows for
 * directed edges, red + thicker highlighting, VertexCoordinates in both forms,
 * labels, styles, edge weights, each GraphLayout), the frame (no axes, equal
 * aspect ratio), the degenerate inputs (empty graph, one vertex, disconnected
 * and edgeless graphs), and 500-vertex graphs in bounded time.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"
#include "graph.h"
#include "glayout.h"
#include "test_utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#define PETERSEN "PetersenGraph[]"
#define H1 "Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]"

static void test_basic_shape(void) {
    assert_eval_eq("Head[GraphPlot[" PETERSEN "]]", "Graphics", 0);
    /* one Line per undirected edge, one Disk per vertex, no labels by default */
    assert_eval_eq("Count[GraphPlot[" PETERSEN "], _Line, Infinity]", "15", 0);
    assert_eval_eq("Count[GraphPlot[" PETERSEN "], _Disk, Infinity]", "10", 0);
    assert_eval_eq("Count[GraphPlot[" PETERSEN "], _Text, Infinity]", "0", 0);
    /* frame: no axes, equal aspect ratio, explicit PlotRange and ImageSize */
    assert_eval_eq("Axes /. Rest[List @@ GraphPlot[" PETERSEN "]]", "False", 0);
    assert_eval_eq("AspectRatio /. Rest[List @@ GraphPlot[" PETERSEN "]]", "Automatic", 0);
    assert_eval_eq("MatchQ[PlotRange /. Rest[List @@ GraphPlot[" PETERSEN "]], "
                   "{{_Real, _Real}, {_Real, _Real}}]", "True", 0);
    assert_eval_eq("MatchQ[ImageSize /. Rest[List @@ GraphPlot[" PETERSEN "]], "
                   "{_Integer, _Integer}]", "True", 0);
    /* default colours: Mathematica's vertex blue */
    assert_eval_eq("MemberQ[GraphPlot[" PETERSEN "], RGBColor[0.368417, 0.506779, 0.709798], Infinity]",
                   "True", 0);
    /* a user option passes through to Graphics */
    assert_eval_eq("PlotLabel /. Rest[List @@ GraphPlot[CycleGraph[3], PlotLabel -> \"C3\"]]",
                   "\"C3\"", 0);
    /* non-graphs and non-option arguments stay unevaluated */
    assert_eval_eq("Head[GraphPlot[5]]", "GraphPlot", 0);
    assert_eval_eq("Head[GraphPlot[CycleGraph[3], 7]]", "GraphPlot", 0);
    /* a list of rules is taken as the graph it denotes */
    assert_eval_eq("Count[GraphPlot[{1 -> 2, 2 -> 3}], _Arrow, Infinity]", "2", 0);
}

static void test_determinism(void) {
    assert_eval_eq("GraphPlot[" PETERSEN "] === GraphPlot[" PETERSEN "]", "True", 0);
    assert_eval_eq("GraphPlot[GridGraph[{5, 4}]] === GraphPlot[GridGraph[{5, 4}]]", "True", 0);
    assert_eval_eq("GraphPlot[CompleteKaryTree[4, 3], VertexLabels -> \"Name\"] === "
                   "GraphPlot[CompleteKaryTree[4, 3], VertexLabels -> \"Name\"]", "True", 0);
    assert_eval_eq("(SeedRandom[7]; g = RandomGraph[{60, 120}]; "
                   "GraphPlot[g, GraphLayout -> \"SpringElectricalEmbedding\"] === "
                   "GraphPlot[g, GraphLayout -> \"SpringElectricalEmbedding\"])", "True", 0);
    assert_eval_eq("HypergraphPlot[" H1 "] === HypergraphPlot[" H1 "]", "True", 0);

    /* The engine itself: bit-identical coordinates on repeated calls, for
     * every method. */
    const int n = 12, m = 18;
    int eu[18], ev[18];
    unsigned char dir[18];
    for (int k = 0; k < m; k++) {
        eu[k] = k % n; ev[k] = (k * 5 + 1) % n;
        if (ev[k] == eu[k]) ev[k] = (ev[k] + 1) % n;
        dir[k] = 0;
    }
    GLMethod ms[] = { GL_AUTOMATIC, GL_CIRCULAR, GL_SPRING, GL_STRESS, GL_LAYERED,
                      GL_BIPARTITE, GL_GRID };
    for (size_t t = 0; t < sizeof(ms) / sizeof(ms[0]); t++) {
        double a[24], b[24];
        int ra = glayout_compute(n, m, eu, ev, dir, ms[t], a);
        int rb = glayout_compute(n, m, eu, ev, dir, ms[t], b);
        ASSERT_MSG(ra >= 0 && ra == rb, "method %d: %d vs %d", (int)ms[t], ra, rb);
        ASSERT_MSG(memcmp(a, b, sizeof(a)) == 0, "method %d not deterministic", (int)ms[t]);
        for (int i = 0; i < 2 * n; i++) ASSERT_MSG(isfinite(a[i]), "method %d: non-finite", (int)ms[t]);
    }
    /* n == 0 is a no-op success */
    ASSERT_MSG(glayout_compute(0, 0, NULL, NULL, NULL, GL_AUTOMATIC, NULL) >= 0, "empty layout");
}

static void test_directed(void) {
    /* directed edges are Arrows (with an Arrowheads directive), undirected Lines */
    assert_eval_eq("Count[GraphPlot[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]], _Arrow, Infinity]", "3", 0);
    assert_eval_eq("Count[GraphPlot[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]], _Line, Infinity]", "0", 0);
    assert_eval_eq("Count[GraphPlot[Graph[{1 -> 2, 2 <-> 3}]], _Arrow, Infinity]", "1", 0);
    assert_eval_eq("Count[GraphPlot[Graph[{1 -> 2, 2 <-> 3}]], _Line, Infinity]", "1", 0);
    assert_eval_eq("Count[GraphPlot[Graph[{1 -> 2}]], _Arrowheads, Infinity] >= 1", "True", 0);
    /* the arrow stops outside the target disk and starts outside the source */
    assert_eval_eq("Module[{g = GraphPlot[Graph[{1 -> 2}], VertexCoordinates -> {{0, 0}, {1, 0}}], r, s, t},"
                   " r = Cases[g, Disk[_, x_] :> x, Infinity][[1]];"
                   " {s, t} = Cases[g, Arrow[{a_, b_}] :> {a, b}, Infinity][[1]];"
                   " {s[[1]] >= r, 1 - t[[1]] > r}]", "{True, True}", 0);
    /* a mutual pair u->v, v->u is drawn as two separated arrows */
    assert_eval_eq("Module[{g = GraphPlot[Graph[{1 -> 2, 2 -> 1}], VertexCoordinates -> {{0, 0}, {1, 0}}]},"
                   " Length[Union[Cases[g, Arrow[{{_, y_}, _}] :> y, Infinity]]]]", "2", 0);
    /* a DAG defaults to a layered drawing: every edge points downwards */
    assert_eval_eq("Module[{g = GraphPlot[Graph[{1 -> 2, 1 -> 3, 2 -> 4, 3 -> 4, 4 -> 5, 2 -> 5}]]},"
                   " And @@ Cases[g, Arrow[{{_, y0_}, {_, y1_}}] :> y1 < y0, Infinity]]", "True", 0);
}

static void test_highlight_and_styles(void) {
    /* highlighted vertex: red, drawn larger */
    assert_eval_eq("Count[GraphPlot[CycleGraph[3]], RGBColor[1., 0., 0.], Infinity]", "0", 0);
    assert_eval_eq("MatchQ[First[GraphPlot[CycleGraph[3], GraphHighlight -> {1}]],"
                   " {___, RGBColor[1., 0., 0.], _Disk, ___}]", "True", 0);
    assert_eval_eq("Length[Union[Cases[GraphPlot[CycleGraph[3], GraphHighlight -> {1}],"
                   " Disk[_, r_] :> r, Infinity]]]", "2", 0);
    /* highlighted edge: red and thicker than the others */
    assert_eval_eq("MatchQ[First[GraphPlot[CycleGraph[3], GraphHighlight -> {1 <-> 2}]],"
                   " {___, RGBColor[1., 0., 0.], _Line, ___}]", "True", 0);
    assert_eval_eq("Module[{t = Cases[GraphPlot[CycleGraph[3], GraphHighlight -> {1 <-> 2}],"
                   " Thickness[x_] :> x, Infinity]}, Length[t] >= 2 && t[[2]] > 2 t[[1]]]", "True", 0);
    /* a highlighted directed edge (Rule syntax) */
    assert_eval_eq("MatchQ[First[GraphPlot[Graph[{1 -> 2, 2 -> 3}], GraphHighlight -> {2 -> 3}]],"
                   " {___, RGBColor[1., 0., 0.], _Arrow, ___}]", "True", 0);
    /* VertexStyle: per-vertex rules and a global colour */
    assert_eval_eq("MatchQ[First[GraphPlot[CycleGraph[3], VertexStyle -> {2 -> RGBColor[0, 1, 0]}]],"
                   " {___, RGBColor[0, 1, 0], _Disk, ___}]", "True", 0);
    assert_eval_eq("Count[First[GraphPlot[CycleGraph[4], VertexStyle -> RGBColor[0, 0, 1]]],"
                   " RGBColor[0, 0, 1]]", "4", 0);
    /* EdgeStyle rules */
    assert_eval_eq("MatchQ[First[GraphPlot[CycleGraph[4], EdgeStyle -> {(2 <-> 3) -> RGBColor[1, 0, 1]}]],"
                   " {___, RGBColor[1, 0, 1], _Line, ___}]", "True", 0);
}

static void test_labels(void) {
    assert_eval_eq("Count[GraphPlot[" PETERSEN ", VertexLabels -> \"Name\"], _Text, Infinity]", "10", 0);
    assert_eval_eq("Count[GraphPlot[" PETERSEN ", VertexLabels -> Automatic], _Text, Infinity]", "10", 0);
    assert_eval_eq("Count[GraphPlot[" PETERSEN ", VertexLabels -> None], _Text, Infinity]", "0", 0);
    assert_eval_eq("Cases[GraphPlot[CycleGraph[3], VertexLabels -> {2 -> \"two\"}],"
                   " Text[s_, __] :> s, Infinity]", "{\"two\"}", 0);
    /* labels are Strings with an offset, so the PDF writer can align them */
    assert_eval_eq("MatchQ[Cases[GraphPlot[CycleGraph[3], VertexLabels -> \"Name\"], _Text, Infinity],"
                   " {Text[_String, {_, _}, {_, _}] ..}]", "True", 0);
    /* labels fit inside the PlotRange */
    assert_eval_eq("Module[{g = GraphPlot[CompleteGraph[5], VertexLabels -> \"Name\"], pr, pts},"
                   " pr = PlotRange /. Rest[List @@ g];"
                   " pts = Cases[g, Text[_, p_, _] :> p, Infinity];"
                   " And @@ (pr[[1, 1]] < #[[1]] < pr[[1, 2]] && pr[[2, 1]] < #[[2]] < pr[[2, 2]] & /@ pts)]",
                   "True", 0);
    /* edge weights */
    assert_eval_eq("Cases[GraphPlot[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}, EdgeWeight -> {5, 9}],"
                   " EdgeLabels -> \"EdgeWeight\"], Text[s_, __] :> s, Infinity]", "{\"5\", \"9\"}", 0);
}

static void test_coordinates_and_layouts(void) {
    /* VertexCoordinates as a list, in VertexList order */
    assert_eval_eq("Cases[GraphPlot[PathGraph[{1, 2, 3}], VertexCoordinates -> {{0, 0}, {1, 0}, {2, 1}}],"
                   " Disk[p_, _] :> p, Infinity]", "{{0.0, 0.0}, {1.0, 0.0}, {2.0, 1.0}}", 0);
    /* ... and as rules for a subset (the rest from the layout) */
    assert_eval_eq("Cases[GraphPlot[CycleGraph[4], VertexCoordinates -> {3 -> {5, 5}}],"
                   " Disk[p_, _] :> p, Infinity][[3]]", "{5.0, 5.0}", 0);
    /* every GraphLayout gives a Graphics with all the vertices */
    const char* layouts[] = { "CircularEmbedding", "SpringElectricalEmbedding", "StressEmbedding",
                              "LayeredEmbedding", "BipartiteEmbedding", "GridEmbedding", NULL };
    for (int i = 0; layouts[i]; i++) {
        char buf[256];
        snprintf(buf, sizeof(buf), "Count[GraphPlot[CycleGraph[6], GraphLayout -> \"%s\"], _Disk, Infinity]",
                 layouts[i]);
        assert_eval_eq(buf, "6", 0);
    }
    /* CircularEmbedding: every vertex at the same distance from the centre */
    assert_eval_eq("Module[{p = Cases[GraphPlot[" PETERSEN ", GraphLayout -> \"CircularEmbedding\"],"
                   " Disk[q_, _] :> q, Infinity], c}, c = Mean[p];"
                   " Max[Norm[# - c] & /@ p] - Min[Norm[# - c] & /@ p] < 10^-4]", "True", 0);
    /* BipartiteEmbedding: two columns */
    assert_eval_eq("Length[Union[Cases[GraphPlot[CompleteGraph[{3, 4}], GraphLayout -> \"BipartiteEmbedding\"],"
                   " Disk[{x_, _}, _] :> x, Infinity]]]", "2", 0);
    /* a grid is drawn as a grid: 4 distinct x and 4 distinct y coordinates */
    assert_eval_eq("Module[{p = Cases[GraphPlot[GridGraph[{4, 4}]], Disk[q_, _] :> q, Infinity]},"
                   " {Length[Union[Round[p[[All, 1]], 0.001]]], Length[Union[Round[p[[All, 2]], 0.001]]]}]",
                   "{4, 4}", 0);
    /* trees default to a layered drawing with the root on top */
    assert_eval_eq("Module[{p = Cases[GraphPlot[CompleteKaryTree[3, 2]], Disk[q_, _] :> q, Infinity]},"
                   " p[[1, 2]] > Max[Rest[p][[All, 2]]]]", "True", 0);
    /* a cycle is a regular polygon: all edges the same length */
    assert_eval_eq("Module[{g = GraphPlot[CycleGraph[7]], l}, l = Cases[g, Line[{a_, b_}] :> Norm[a - b], Infinity];"
                   " Max[l] - Min[l] < 10^-3 Max[l]]", "True", 0);
}

static void test_degenerate(void) {
    assert_eval_eq("Head[GraphPlot[Graph[{}, {}]]]", "Graphics", 0);
    assert_eval_eq("Count[GraphPlot[Graph[{}, {}]], _Disk, Infinity]", "0", 0);
    assert_eval_eq("Count[GraphPlot[Graph[{a}, {}], VertexLabels -> \"Name\"], _Disk, Infinity]", "1", 0);
    /* an edgeless graph: distinct positions for every vertex */
    assert_eval_eq("Length[Union[Cases[GraphPlot[Graph[Range[10], {}]], Disk[p_, _] :> p, Infinity]]]", "10", 0);
    /* disconnected: components do not overlap (distinct positions) */
    assert_eval_eq("Module[{g = GraphDisjointUnion[CycleGraph[5], CompleteGraph[4]], p},"
                   " p = Cases[GraphPlot[g], Disk[q_, _] :> q, Infinity]; {Length[p], Length[Union[p]]}]",
                   "{9, 9}", 0);
    /* a two-vertex graph and a path */
    assert_eval_eq("Count[GraphPlot[Graph[{1 <-> 2}]], _Line, Infinity]", "1", 0);
    assert_eval_eq("Count[GraphPlot[PathGraph[Range[5]]], _Disk, Infinity]", "5", 0);
}

static double now_s(void) { return (double)clock() / CLOCKS_PER_SEC; }

static void test_large(void) {
    double t0 = now_s();
    assert_eval_eq("(SeedRandom[3]; Count[GraphPlot[RandomGraph[{500, 1500}]], _Disk, Infinity])", "500", 0);
    assert_eval_eq("Count[GraphPlot[GridGraph[{25, 20}]], _Disk, Infinity]", "500", 0);
    assert_eval_eq("Count[GraphPlot[CompleteKaryTree[9, 2]], _Disk, Infinity]", "511", 0);
    assert_eval_eq("Count[GraphPlot[Graph[Range[500], {}]], _Disk, Infinity]", "500", 0);
    assert_eval_eq("(SeedRandom[4]; Count[GraphPlot[RandomGraph[{500, 1500}],"
                   " GraphLayout -> \"SpringElectricalEmbedding\"], _Disk, Infinity])", "500", 0);
    assert_eval_eq("(SeedRandom[5]; Count[GraphPlot[RandomGraph[{1500, 3000}]], _Disk, Infinity])", "1500", 0);
    double dt = now_s() - t0;
    ASSERT_MSG(dt < 20.0, "large graphs took %.2fs", dt);
    printf("(%.2fs) ", dt);
}

static void test_hypergraph_plot(void) {
    assert_eval_eq("Head[HypergraphPlot[" H1 "]]", "Graphics", 0);
    /* one translucent shape per hyperedge (sizes 3, 2, 3 and 1), one disk per vertex */
    assert_eval_eq("Count[HypergraphPlot[" H1 "], _Polygon, Infinity]", "4", 0);
    assert_eval_eq("Count[HypergraphPlot[" H1 "], _Disk, Infinity]", "7", 0);
    assert_eval_eq("MemberQ[HypergraphPlot[" H1 "], Opacity[x_ /; x < 1], Infinity]", "True", 0);
    /* hyperedges get distinct colours */
    assert_eval_eq("Module[{p = First[HypergraphPlot[" H1 "]]},"
                   " Length[Union[Cases[p, c_RGBColor /; MemberQ[p, c], 1]]] >= 4]", "True", 0);
    assert_eval_eq("Axes /. Rest[List @@ HypergraphPlot[" H1 "]]", "False", 0);
    assert_eval_eq("Count[HypergraphPlot[" H1 ", VertexLabels -> \"Name\"], _Text, Infinity]", "7", 0);
    /* a plain list of hyperedges */
    assert_eval_eq("Count[HypergraphPlot[{{1, 2, 3}, {3, 4}}], _Polygon, Infinity]", "2", 0);
    /* VertexCoordinates: every rounded hull encloses its members */
    assert_eval_eq("Cases[HypergraphPlot[Hypergraph[{{1, 2}}], VertexCoordinates -> {{0, 0}, {3, 0}}],"
                   " Disk[p_, _] :> p, Infinity]", "{{0.0, 0.0}, {3.0, 0.0}}", 0);
    assert_eval_eq("Module[{g = HypergraphPlot[Hypergraph[{{1, 2, 3}}],"
                   " VertexCoordinates -> {{0, 0}, {2, 0}, {1, 2}}], poly},"
                   " poly = Cases[g, Polygon[p_] :> p, Infinity][[1]];"
                   " {Min[poly[[All, 1]]] < 0, Max[poly[[All, 1]]] > 2, Max[poly[[All, 2]]] > 2}]",
                   "{True, True, True}", 0);
    /* degenerate: empty and edgeless */
    assert_eval_eq("Head[HypergraphPlot[Hypergraph[{}, {}]]]", "Graphics", 0);
    assert_eval_eq("Count[HypergraphPlot[Hypergraph[{1, 2, 3}, {}]], _Disk, Infinity]", "3", 0);
    assert_eval_eq("Head[HypergraphPlot[5]]", "HypergraphPlot", 0);
    /* larger hypergraph in bounded time */
    double t0 = now_s();
    assert_eval_eq("(SeedRandom[2]; Count[HypergraphPlot[RandomHypergraph[{300, 200}, 3]], _Disk, Infinity])",
                   "300", 0);
    ASSERT_MSG(now_s() - t0 < 10.0, "hypergraph plot too slow");
}

static void test_pdf_export(void) {
    assert_eval_eq("Export[\"/tmp/mathilda_test_graphplot.pdf\", GraphPlot[" PETERSEN
                   ", VertexLabels -> \"Name\"]]", "\"/tmp/mathilda_test_graphplot.pdf\"", 0);
    assert_eval_eq("Export[\"/tmp/mathilda_test_hypergraphplot.pdf\", HypergraphPlot[" H1 "]]",
                   "\"/tmp/mathilda_test_hypergraphplot.pdf\"", 0);
    remove("/tmp/mathilda_test_graphplot.pdf");
    remove("/tmp/mathilda_test_hypergraphplot.pdf");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_basic_shape);
    TEST(test_determinism);
    TEST(test_directed);
    TEST(test_highlight_and_styles);
    TEST(test_labels);
    TEST(test_coordinates_and_layouts);
    TEST(test_degenerate);
    TEST(test_large);
    TEST(test_hypergraph_plot);
    TEST(test_pdf_export);
    printf("All graphplot tests passed!\n");
    return 0;
}
