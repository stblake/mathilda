#ifndef GLAYOUT_H
#define GLAYOUT_H

/* glayout.h - graph layout engine and shared drawing helpers for GraphPlot and
 * HypergraphPlot (src/graph/glayout.c, graphplot.c, hyp_plot.c).
 *
 * The layout half is pure numerics over integer edge arrays: it knows nothing
 * about Expr, so it can be driven from a Graph (GraphPlot), from the star
 * expansion of a Hypergraph (HypergraphPlot), or from a unit test. Every
 * algorithm is deterministic -- no random numbers, fixed iteration orders --
 * so the same graph always yields bit-identical coordinates.
 *
 * Coordinates come out in "edge units": a typical edge has length about 1.
 */

#include "expr.h"

typedef enum {
    GL_AUTOMATIC = 0,  /* trees/DAGs -> layered, paths/cycles/rest -> stress   */
    GL_CIRCULAR,       /* "CircularEmbedding": VertexList order on a circle     */
    GL_SPRING,         /* "SpringElectricalEmbedding": Hu's spring-electrical   */
    GL_STRESS,         /* "StressEmbedding": SMACOF on BFS distances            */
    GL_LAYERED,        /* "LayeredEmbedding" / "LayeredDigraphEmbedding"        */
    GL_BIPARTITE,      /* "BipartiteEmbedding": the two parts in two columns    */
    GL_GRID            /* "GridEmbedding": VertexList order on a square grid    */
} GLMethod;

/* Lays out n vertices joined by the m edges eu[k] -- ev[k] (0-based vertex
 * indices; dir[k] != 0 marks a directed edge eu -> ev, dir may be NULL for all
 * undirected). Writes vertex i's position to xy[2i], xy[2i+1]. Every connected
 * component is laid out separately and the components are then packed side by
 * side. Returns the method actually used (GL_AUTOMATIC resolved; an inapplicable
 * choice, e.g. BipartiteEmbedding of an odd cycle, falls back to GL_STRESS), or
 * -1 on allocation failure. n == 0 is a no-op success. */
int glayout_compute(int n, int m, const int* eu, const int* ev,
                    const unsigned char* dir, GLMethod method, double* xy);

/* Parses a GraphLayout option value ("StressEmbedding", Automatic, ...).
 * Returns 1 and sets *out on a recognised value, 0 otherwise. */
int glayout_parse_method(const Expr* v, GLMethod* out);

/* ---- Shared drawing helpers (graphplot.c) --------------------------------
 * Used by both GraphPlot and HypergraphPlot so the two draw vertices, labels,
 * margins and options identically. */

/* The value of option `name` (Rule or RuleDelayed with a Symbol lhs) among
 * res's arguments from position `first` on, or NULL. Later options do NOT
 * override earlier ones: the first occurrence wins, as in Mathematica. */
const Expr* gd_option(const Expr* res, size_t first, const char* name);

/* Applies a VertexCoordinates spec -- {{x,y}, ...} in VertexList order, or
 * {v -> {x,y}, ...} for any subset of vertices -- onto xy. vindex(v) maps a
 * vertex to its position (or -1). Returns the number of vertices set; a
 * malformed spec sets none. */
int gd_apply_vertex_coordinates(const Expr* spec, int n, double* xy,
                                int (*vindex)(const void* ctx, const Expr* v),
                                const void* ctx);

/* A growable array of owned primitive Exprs. */
typedef struct { Expr** p; size_t n, cap; int oom; } GDPrims;
void gd_push(GDPrims* P, Expr* e);        /* takes ownership (NULL -> oom)  */
void gd_prims_free(GDPrims* P);

/* Small constructors. */
Expr* gd_pt(double x, double y);                       /* {x, y}            */
Expr* gd_rgb(double r, double g, double b);            /* RGBColor[r, g, b] */
Expr* gd_head1(const char* head, Expr* a);             /* head[a]           */
Expr* gd_head2(const char* head, Expr* a, Expr* b);    /* head[a, b]        */

/* The style given for vertex/edge by a VertexStyle/EdgeStyle-like spec: a
 * List of rules (or one Rule) looked up with `match`, else the spec itself as a
 * global style. NULL when no style applies. Borrowed. */
const Expr* gd_style_for(const Expr* spec, const void* item,
                         int (*match)(const Expr* key, const void* item));

/* Label text for an arbitrary expression: a String's contents, else its
 * printed form. Caller frees. */
char* gd_label_text(const Expr* e);

/* Label policy decoded from a VertexLabels value. */
typedef enum { GD_LBL_NONE, GD_LBL_NAME, GD_LBL_RULES } GDLabelMode;
GDLabelMode gd_label_mode(const Expr* v);

/* The rhs of the first rule `key -> rhs` in the List `rules` whose lhs is
 * SameQ key, or NULL. */
const Expr* gd_rule_lookup(const Expr* rules, const Expr* key);

/* Frames the picture. On entry bb (xmin, xmax, ymin, ymax, world units) boxes
 * everything drawn except labels; on exit it also boxes the labels plus a small
 * margin, and *pw x *ph is the page size in points at which the world maps with
 * equal x/y scale. texts[i] (NULL = unlabelled; texts may be NULL) labels the
 * vertex at xy[2i], xy[2i+1] with disk radius r; the label goes on the side with
 * the widest angular gap between the directions ang[aoff[i] .. aoff[i+1]-1]
 * (incident edges), preferring the upper right. Label Text primitives are
 * appended to L. The page size honours an ImageSize option in res. */
void gd_frame(double* bb, const Expr* res, size_t first, int n,
              const double* xy, double r, char* const* texts,
              const int* aoff, const double* ang, GDPrims* L,
              double* pw, double* ph);

/* Assembles Graphics[prims, PlotRange -> bb, AspectRatio -> Automatic,
 * Axes -> False, ImageSize -> {pw, ph}, passthrough...]. Any option of res
 * (from `first` on) whose name is not in the NULL-terminated `consumed` list
 * is passed through to Graphics (PlotLabel, Background, ...). Consumes P. */
Expr* gd_finish(GDPrims* P, const double* bb, double pw, double ph,
                const Expr* res, size_t first, const char* const* consumed);

/* The Mathematica ColorData[97] palette, cycled by index. */
void gd_palette(int k, double* r, double* g, double* b);

/* Vertex disks (default colour, or vstyle[i] when non-NULL; hl[i] marks a
 * highlighted vertex, drawn red and 15% larger), each with a thin darker rim.
 * plot_w_pt is the plot width in points (Thickness is relative to it). */
void gd_emit_vertices(GDPrims* P, int n, const double* xy, double r,
                      const Expr* const* vstyle, const unsigned char* hl,
                      double plot_w_pt);

/* Label texts per a VertexLabels value for the vertices of the List verts
 * (NULL when labels are off; entries NULL for unlabelled vertices). */
char** gd_vertex_texts(const Expr* spec, int n, const Expr* verts);
void   gd_free_texts(char** t, int n);

/* Typical spacing of a drawing: the median edge length, else the median
 * nearest-neighbour distance, else 1. */
double gd_drawing_unit(int n, int m, const int* eu, const int* ev, const double* xy);

/* Vertex disk radius for a drawing of spacing `unit`; also writes the
 * vertices' bounding box (xmin, xmax, ymin, ymax) to bb. With nncap the
 * radius is also held under 0.3 x the median nearest-neighbour distance
 * (layered drawings, whose leaves sit closer than their edges are long). */
double gd_vertex_radius(int n, const double* xy, double unit, int nncap, double* bb);

/* A lone vertex (or a cluster smaller than one unit in both directions)
 * would fill the page: widen bb to a unit square about its centre. */
void gd_min_extent(double* bb, double unit);

/* expr_eq(key, (const Expr*)item): the vertex matcher for gd_style_for. */
int gd_vertex_match(const Expr* key, const void* item);

/* Is item in a GraphHighlight-like spec (a List of items, or one item)? */
int gd_in_highlight(const Expr* hl, const void* item,
                    int (*match)(const Expr*, const void*));

/* Default vertex colour (Mathematica's first ColorData[97] entry). */
#define GD_VERTEX_R 0.368417
#define GD_VERTEX_G 0.506779
#define GD_VERTEX_B 0.709798

#endif /* GLAYOUT_H */
