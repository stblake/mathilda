/* Mathilda -- planar triangular meshes for the finite element method.
 *
 * A region (Rectangle, Polygon, Triangle, Disk and its sectors, Annulus) is
 * described by its boundary: closed loops of vertices, each edge a straight
 * segment or an arc of a circle, oriented so that the region lies on the
 * LEFT of every edge (outer loops counter-clockwise, holes clockwise).
 * fem_mesh_build triangulates it by Delaunay refinement (Ruppert's
 * algorithm with Shewchuk's concentric-shell segment splitting):
 *
 *   1. the loop vertices are inserted into a Delaunay triangulation of a
 *      large enclosing triangle (Bowyer-Watson);
 *   2. every boundary segment that is missing from the triangulation, or
 *      that has a vertex inside its diametral circle ("encroached"), is
 *      split -- at its midpoint, at a power-of-two distance from an input
 *      corner (the concentric shells that keep small corner angles from
 *      cascading), or at its mid-angle when it is an arc, so that curved
 *      boundaries are resolved by the refinement itself;
 *   3. each triangle inside the region that is larger than the requested
 *      cell measure or whose circumradius-to-shortest-edge ratio exceeds
 *      sqrt(2) (every angle >= 20.7 degrees) gets a vertex at its
 *      circumcentre, unless that vertex would encroach a segment, which is
 *      then split instead.
 *
 * The result is a conforming Delaunay mesh of the region's interior: every
 * boundary edge is a mesh edge, no triangle is obtuse against the boundary,
 * and every angle is bounded away from zero (except near input corners
 * sharper than 60 degrees, where the guarantee does not hold and the
 * refinement is cut off by a minimum edge length).
 *
 * fem_mesh_p2 then adds the quadratic-element nodes: one node per edge, at
 * the edge midpoint, or on the arc (at its mid-angle) for a curved boundary
 * edge, which makes those elements isoparametric.
 */
#ifndef MATHILDA_FEM_MESH_H
#define MATHILDA_FEM_MESH_H

#include <stdbool.h>
#include <stdint.h>

/* A circle carrying curved boundary edges. */
typedef struct { double cx, cy, r; } FemArc;

#define FEM_MAX_ARCS 16

/* Boundary description. Loop k has loop_n[k] vertices starting at
 * vertex index loop_start[k] of xy (2 doubles per vertex); edge i of the
 * loop runs from its vertex i to vertex i+1 (cyclically) and lies on arc
 * edge_arc[global vertex index] (-1: straight). */
typedef struct {
    int nloop;
    int loop_start[8], loop_n[8];
    int64_t nv;
    double* xy;
    int* edge_arc;
    /* optional (NULL: unknown, treated as one boundary-condition type
     * everywhere): 1 when edge i carries a Dirichlet condition, 0 a natural
     * one; a corner where the type changes is singular (exponent
     * pi/(2 omega)) and is graded like a re-entrant corner */
    unsigned char* edge_dir;
    int narc;
    FemArc arc[FEM_MAX_ARCS];
} FemGeom;

void fem_geom_init(FemGeom* g);
void fem_geom_free(FemGeom* g);
/* Append a loop of n vertices; arcs[i] is the arc of edge i (NULL: all
 * straight). The loop is reoriented if needed: `hole` false wants it
 * counter-clockwise, true clockwise. False when there is no room. */
bool fem_geom_add_loop(FemGeom* g, const double* xy, const int* arcs, int n, bool hole);
int fem_geom_add_arc(FemGeom* g, double cx, double cy, double r);

/* The quadratic (P2) mesh. Node k is at (x[2k], x[2k+1]). Element e has
 * nodes tri[6e .. 6e+5]: the three vertices counter-clockwise, then the
 * midside nodes of edges (0,1), (1,2), (2,0). curved[e] is true when one of
 * its midside nodes was placed on an arc. A boundary edge b has vertex
 * nodes bedge[3b], bedge[3b+1] (the region on the left) and midside node
 * bedge[3b+2]. bnode[k] marks the nodes on the boundary. */
typedef struct {
    int64_t nn, ne, nb, nvert;
    double* x;
    int64_t* tri;
    unsigned char* curved;
    int64_t* bedge;
    unsigned char* bnode;
    double bbox[4];          /* xmin, xmax, ymin, ymax */
} FemMesh;

typedef struct {
    double max_area;         /* largest triangle area (MaxCellMeasure) */
    double max_bedge;        /* largest boundary edge length; <= 0: none */
    int64_t max_vertices;    /* safety cap */
    bool grade_corners;      /* refine towards re-entrant corners */
} FemMeshOpts;

/* Triangulate g into a P2 mesh. False (with *why set to a short reason)
 * when the region cannot be meshed. */
bool fem_mesh_build(const FemGeom* g, const FemMeshOpts* o, FemMesh* m, const char** why);
void fem_mesh_free(FemMesh* m);

/* Uniform refinement: every element split into four through its six P2
 * nodes (which become the refined mesh's vertices, same indices), the new
 * midside nodes by the parent's quadratic map -- on a curved boundary edge,
 * on the circle through that edge's three nodes. The refined mesh is nested
 * in the parent, so h is exactly halved element by element (what Richardson
 * extrapolation in h needs). False on allocation failure. */
bool fem_mesh_refine(const FemMesh* in, FemMesh* out);

/* Signed area of the region described by g (holes subtract). */
double fem_geom_area(const FemGeom* g);

#endif /* MATHILDA_FEM_MESH_H */
