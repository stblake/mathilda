/* Mathilda -- planar Delaunay-refinement meshing for the finite element method.
 *
 * See fem_mesh.h for the algorithm. Implementation notes:
 *
 *   - Triangles are kept in one array with neighbour links (nb[i] is the
 *     triangle across the edge opposite vertex v[i]); dead triangles are
 *     recycled through a free list. Every live triangle is counter-clockwise.
 *   - Insertion is Bowyer-Watson: the cavity of triangles whose circumcircle
 *     contains the new point, grown from the triangle that contains it and
 *     never across a boundary segment (so the triangulation stays
 *     constrained-Delaunay even when rounding makes an incircle test lie),
 *     re-triangulated as a fan from the new point. A cavity that is not
 *     star-shaped from the point (rounding again) is shrunk until it is.
 *   - Boundary segments are directed with the region on their left. Each
 *     triangle carries an inside flag: a new triangle inherits the flag of
 *     the cavity triangle whose edge it is built on, except that a triangle
 *     with a segment edge takes the flag from the side of that segment it
 *     lies on -- which keeps the flags exact while arcs are split (an arc
 *     midpoint is off the old chord, so a cavity can straddle it).
 *   - A segment is found encroached either directly (an apex of one of its
 *     two triangles sees it at an angle >= 90 degrees; for a Delaunay edge
 *     that is equivalent to a vertex in its diametral circle) or, for a
 *     candidate circumcentre, when the point location walk would cross it
 *     or it bounds the candidate's cavity and the candidate sees it at
 *     >= 90 degrees. The candidate is then not inserted; the segment is
 *     split and the triangle queued again.
 */
#include "fem_mesh.h"
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  Geometry description                                               */
/* ------------------------------------------------------------------ */

void fem_geom_init(FemGeom* g) { memset(g, 0, sizeof *g); }

void fem_geom_free(FemGeom* g) {
    free(g->xy); free(g->edge_arc); free(g->edge_dir);
    memset(g, 0, sizeof *g);
}

int fem_geom_add_arc(FemGeom* g, double cx, double cy, double r) {
    if (g->narc >= FEM_MAX_ARCS) return -1;
    g->arc[g->narc].cx = cx; g->arc[g->narc].cy = cy; g->arc[g->narc].r = r;
    return g->narc++;
}

bool fem_geom_add_loop(FemGeom* g, const double* xy, const int* arcs, int n, bool hole) {
    if (n < 3 || g->nloop >= 8) return false;
    double a2 = 0.0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        a2 += xy[2 * i] * xy[2 * j + 1] - xy[2 * j] * xy[2 * i + 1];
    }
    bool ccw = a2 > 0.0;
    bool flip = hole ? ccw : !ccw;
    g->xy = realloc(g->xy, sizeof(double) * 2 * (size_t)(g->nv + n));
    g->edge_arc = realloc(g->edge_arc, sizeof(int) * (size_t)(g->nv + n));
    int64_t s = g->nv;
    for (int i = 0; i < n; i++) {
        /* reversed: vertex order n-1 .. 0; edge i (v'_i -> v'_{i+1}) is the old
         * edge from vertex n-1-i-1 to n-1-i, i.e. arcs[n-2-i] (cyclically) */
        int src = flip ? (n - 1 - i) : i;
        g->xy[2 * (s + i)] = xy[2 * src];
        g->xy[2 * (s + i) + 1] = xy[2 * src + 1];
        int ea = -1;
        if (arcs) ea = flip ? arcs[((n - 2 - i) % n + n) % n] : arcs[i];
        g->edge_arc[s + i] = ea;
    }
    g->loop_start[g->nloop] = (int)s;
    g->loop_n[g->nloop] = n;
    g->nloop++;
    g->nv += n;
    return true;
}

double fem_geom_area(const FemGeom* g) {
    double area = 0.0;
    for (int k = 0; k < g->nloop; k++) {
        int s = g->loop_start[k], n = g->loop_n[k];
        for (int i = 0; i < n; i++) {
            const double* p = g->xy + 2 * (s + i);
            const double* q = g->xy + 2 * (s + (i + 1) % n);
            area += 0.5 * (p[0] * q[1] - q[0] * p[1]);
            int ai = g->edge_arc[s + i];
            if (ai >= 0) {
                const FemArc* A = &g->arc[ai];
                double ta = atan2(p[1] - A->cy, p[0] - A->cx);
                double tb = atan2(q[1] - A->cy, q[0] - A->cx);
                double d = tb - ta;
                while (d > M_PI) d -= 2 * M_PI;
                while (d <= -M_PI) d += 2 * M_PI;
                double th = fabs(d);
                double seg = 0.5 * A->r * A->r * (th - sin(th));
                /* the circular segment lies on the region side when the
                 * centre is on the left of the edge (convex boundary) */
                double side = (q[0] - p[0]) * (A->cy - p[1]) - (q[1] - p[1]) * (A->cx - p[0]);
                area += side > 0 ? seg : -seg;
            }
        }
    }
    return area;
}

/* ------------------------------------------------------------------ */
/*  Triangulation state                                                */
/* ------------------------------------------------------------------ */

typedef struct { int v[3]; int nb[3]; unsigned char alive, inside, noref; } Tri;
typedef struct { int a, b, arc; unsigned char alive, queued, force; } Seg;

typedef struct {
    double* px; double* py; unsigned char* corner; int* vtri;
    int np, pcap;
    Tri* t; int nt, tcap; int* freel; int nfree, fcap;
    int* stamp; int curstamp;
    Seg* s; int ns, scap;
    /* edge -> segment hash (open addressing, keys min*2^32+max) */
    int64_t* hkey; int* hval; int64_t hcap; int64_t hused;
    const FemGeom* g;
    unsigned int rng;
    /* scratch */
    int* cav; int ncav, ccap;
    int (*bnd)[4]; int nbnd, bcap;      /* a, b, owner, outer */
    int* newt; int nnew, ncap;
    int* segq; int nsegq, sqcap;
    int* triq; int ntriq, tqhead, tqcap;
} Mesh;

static double orient(const Mesh* M, int a, int b, int c) {
    return (M->px[b] - M->px[a]) * (M->py[c] - M->py[a])
         - (M->py[b] - M->py[a]) * (M->px[c] - M->px[a]);
}
static double orient_p(const Mesh* M, int a, int b, double x, double y) {
    return (M->px[b] - M->px[a]) * (y - M->py[a]) - (M->py[b] - M->py[a]) * (x - M->px[a]);
}
/* > 0 when (x, y) is inside the circumcircle of the CCW triangle t */
static double incircle(const Mesh* M, int t, double x, double y) {
    const int* v = M->t[t].v;
    double adx = M->px[v[0]] - x, ady = M->py[v[0]] - y;
    double bdx = M->px[v[1]] - x, bdy = M->py[v[1]] - y;
    double cdx = M->px[v[2]] - x, cdy = M->py[v[2]] - y;
    double ad = adx * adx + ady * ady, bd = bdx * bdx + bdy * bdy, cd = cdx * cdx + cdy * cdy;
    return adx * (bdy * cd - bd * cdy) - ady * (bdx * cd - bd * cdx) + ad * (bdx * cdy - bdy * cdx);
}

static int64_t ekey(int a, int b) {
    int lo = a < b ? a : b, hi = a < b ? b : a;
    return ((int64_t)lo << 32) | (int64_t)(unsigned int)hi;
}
static uint64_t hmix(int64_t k) {
    uint64_t x = (uint64_t)k * UINT64_C(0x9E3779B97F4A7C15);
    return x ^ (x >> 29);
}
static void hash_grow(Mesh* M);
static void hash_put(Mesh* M, int a, int b, int seg) {
    if (2 * (M->hused + 1) > M->hcap) hash_grow(M);
    int64_t k = ekey(a, b);
    uint64_t i = hmix(k) & (uint64_t)(M->hcap - 1);
    while (M->hkey[i] != -1 && M->hkey[i] != k && M->hkey[i] != -2) i = (i + 1) & (uint64_t)(M->hcap - 1);
    if (M->hkey[i] != k) M->hused++;
    M->hkey[i] = k; M->hval[i] = seg;
}
static int hash_get(const Mesh* M, int a, int b) {
    if (!M->hcap) return -1;
    int64_t k = ekey(a, b);
    uint64_t i = hmix(k) & (uint64_t)(M->hcap - 1);
    while (M->hkey[i] != -1) {
        if (M->hkey[i] == k) return M->hval[i];
        i = (i + 1) & (uint64_t)(M->hcap - 1);
    }
    return -1;
}
static void hash_del(Mesh* M, int a, int b) {
    int64_t k = ekey(a, b);
    uint64_t i = hmix(k) & (uint64_t)(M->hcap - 1);
    while (M->hkey[i] != -1) {
        if (M->hkey[i] == k) { M->hkey[i] = -2; return; }
        i = (i + 1) & (uint64_t)(M->hcap - 1);
    }
}
static void hash_grow(Mesh* M) {
    int64_t oc = M->hcap; int64_t* ok = M->hkey; int* ov = M->hval;
    M->hcap = oc ? 2 * oc : 1024;
    M->hkey = malloc(sizeof(int64_t) * (size_t)M->hcap);
    M->hval = malloc(sizeof(int) * (size_t)M->hcap);
    for (int64_t i = 0; i < M->hcap; i++) M->hkey[i] = -1;
    M->hused = 0;
    for (int64_t i = 0; i < oc; i++)
        if (ok[i] >= 0) {
            int a = (int)(ok[i] >> 32), b = (int)(ok[i] & 0xffffffff);
            hash_put(M, a, b, ov[i]);
        }
    free(ok); free(ov);
}

/* the live segment on edge (a, b), or -1 */
static int seg_of(const Mesh* M, int a, int b) {
    int s = hash_get(M, a, b);
    return (s >= 0 && M->s[s].alive) ? s : -1;
}

static int add_point(Mesh* M, double x, double y, bool corner) {
    if (M->np == M->pcap) {
        M->pcap = M->pcap ? 2 * M->pcap : 256;
        M->px = realloc(M->px, sizeof(double) * (size_t)M->pcap);
        M->py = realloc(M->py, sizeof(double) * (size_t)M->pcap);
        M->corner = realloc(M->corner, (size_t)M->pcap);
        M->vtri = realloc(M->vtri, sizeof(int) * (size_t)M->pcap);
    }
    M->px[M->np] = x; M->py[M->np] = y; M->corner[M->np] = corner; M->vtri[M->np] = -1;
    return M->np++;
}

static int new_tri(Mesh* M) {
    int t;
    if (M->nfree) t = M->freel[--M->nfree];
    else {
        if (M->nt == M->tcap) {
            M->tcap = M->tcap ? 2 * M->tcap : 512;
            M->t = realloc(M->t, sizeof(Tri) * (size_t)M->tcap);
            M->stamp = realloc(M->stamp, sizeof(int) * (size_t)M->tcap);
        }
        t = M->nt++;
    }
    M->stamp[t] = 0;
    M->t[t].alive = 1; M->t[t].inside = 0; M->t[t].noref = 0;
    return t;
}
static void kill_tri(Mesh* M, int t) {
    M->t[t].alive = 0;
    if (M->nfree == M->fcap) {
        M->fcap = M->fcap ? 2 * M->fcap : 256;
        M->freel = realloc(M->freel, sizeof(int) * (size_t)M->fcap);
    }
    M->freel[M->nfree++] = t;
}

static int add_seg(Mesh* M, int a, int b, int arc) {
    if (M->ns == M->scap) {
        M->scap = M->scap ? 2 * M->scap : 256;
        M->s = realloc(M->s, sizeof(Seg) * (size_t)M->scap);
    }
    int s = M->ns++;
    M->s[s].a = a; M->s[s].b = b; M->s[s].arc = arc; M->s[s].alive = 1; M->s[s].queued = 0; M->s[s].force = 0;
    hash_put(M, a, b, s);
    return s;
}

static void push_segq(Mesh* M, int s) {
    if (!M->s[s].alive || M->s[s].queued) return;
    if (M->nsegq == M->sqcap) {
        M->sqcap = M->sqcap ? 2 * M->sqcap : 256;
        M->segq = realloc(M->segq, sizeof(int) * (size_t)M->sqcap);
    }
    M->s[s].queued = 1;
    M->segq[M->nsegq++] = s;
}
static void push_triq(Mesh* M, int t) {
    if (M->ntriq == M->tqcap) {
        if (M->tqhead > 0) {           /* compact the consumed prefix */
            memmove(M->triq, M->triq + M->tqhead, sizeof(int) * (size_t)(M->ntriq - M->tqhead));
            M->ntriq -= M->tqhead; M->tqhead = 0;
        }
        if (M->ntriq == M->tqcap) {
            M->tqcap = M->tqcap ? 2 * M->tqcap : 1024;
            M->triq = realloc(M->triq, sizeof(int) * (size_t)M->tqcap);
        }
    }
    M->triq[M->ntriq++] = t;
}

/* Walk from triangle t towards (x, y). Returns the containing triangle, or
 * -1 on failure; when stop_seg and the walk must cross a segment, returns -2
 * with *hit set to it. */
static int locate(Mesh* M, int t, double x, double y, bool stop_seg, int* hit) {
    int steps = 0, maxsteps = 4 * M->nt + 100;
    while (t >= 0 && steps++ < maxsteps) {
        const Tri* T = &M->t[t];
        int start = (int)(M->rng = M->rng * 1103515245u + 12345u) >> 16;
        int moved = 0;
        for (int k = 0; k < 3; k++) {
            int i = (start + k) % 3;
            if (i < 0) i += 3;
            int a = T->v[(i + 1) % 3], b = T->v[(i + 2) % 3];
            if (orient_p(M, a, b, x, y) < 0) {
                if (stop_seg) {
                    int s = seg_of(M, a, b);
                    if (s >= 0) { *hit = s; return -2; }
                }
                t = T->nb[i];
                moved = 1;
                break;
            }
        }
        if (!moved) return t;
    }
    return -1;
}

/* Compute the cavity of (x, y) from triangle t0 into M->cav / M->bnd.
 * False when the point coincides with a vertex. */
static bool cavity(Mesh* M, int t0, double x, double y) {
    for (int k = 0; k < 3; k++) {
        int v = M->t[t0].v[k];
        if (M->px[v] == x && M->py[v] == y) return false;
    }
    /* triangles excluded from this cavity (rounding broke star-shapedness)
     * carry the stamp -excl for every round of this call */
    int excl = ++M->curstamp;
    for (;;) {
        M->curstamp++;
        M->ncav = 0; M->nbnd = 0;
        if (M->ccap < 16) { M->ccap = 64; M->cav = realloc(M->cav, sizeof(int) * (size_t)M->ccap); }
        M->cav[M->ncav++] = t0; M->stamp[t0] = M->curstamp;
        for (int q = 0; q < M->ncav; q++) {
            int t = M->cav[q];
            for (int i = 0; i < 3; i++) {
                int n = M->t[t].nb[i];
                int a = M->t[t].v[(i + 1) % 3], b = M->t[t].v[(i + 2) % 3];
                bool inc = false;
                if (n >= 0 && M->stamp[n] == M->curstamp) continue;      /* interior edge */
                if (n >= 0 && M->stamp[n] != -excl && seg_of(M, a, b) < 0
                    && incircle(M, n, x, y) > 0) inc = true;
                if (inc) {
                    if (M->ncav == M->ccap) {
                        M->ccap *= 2; M->cav = realloc(M->cav, sizeof(int) * (size_t)M->ccap);
                    }
                    M->stamp[n] = M->curstamp;
                    M->cav[M->ncav++] = n;
                } else {
                    if (M->nbnd == M->bcap) {
                        M->bcap = M->bcap ? 2 * M->bcap : 64;
                        M->bnd = realloc(M->bnd, sizeof(int[4]) * (size_t)M->bcap);
                    }
                    M->bnd[M->nbnd][0] = a; M->bnd[M->nbnd][1] = b;
                    M->bnd[M->nbnd][2] = t; M->bnd[M->nbnd][3] = n;
                    M->nbnd++;
                }
            }
        }
        /* edges recorded before their neighbour joined the cavity are not
         * boundary edges: drop them */
        int w = 0;
        for (int k = 0; k < M->nbnd; k++) {
            int n = M->bnd[k][3];
            if (n >= 0 && M->stamp[n] == M->curstamp) continue;
            if (w != k) memcpy(M->bnd[w], M->bnd[k], sizeof(int[4]));
            w++;
        }
        M->nbnd = w;
        /* star-shapedness: every boundary edge must see the point on its left */
        int bad = -1;
        for (int k = 0; k < M->nbnd && bad < 0; k++) {
            double o = orient_p(M, M->bnd[k][0], M->bnd[k][1], x, y);
            if (!(o > 0) && M->bnd[k][2] != t0) bad = M->bnd[k][2];
        }
        if (bad < 0) {
            for (int k = 0; k < M->nbnd; k++)
                if (!(orient_p(M, M->bnd[k][0], M->bnd[k][1], x, y) > 0)) return false;
            return true;
        }
        /* exclude the offending triangle and recompute */
        M->stamp[bad] = -excl;
    }
}

/* Replace the current cavity with the fan from point p. */
static void commit(Mesh* M, int p) {
    M->nnew = 0;
    if (M->ncap < M->nbnd) {
        M->ncap = 2 * M->nbnd + 16;
        M->newt = realloc(M->newt, sizeof(int) * (size_t)M->ncap);
    }
    unsigned char* oldin = malloc((size_t)M->nbnd);
    for (int k = 0; k < M->nbnd; k++) oldin[k] = M->t[M->bnd[k][2]].inside;
    for (int k = 0; k < M->ncav; k++) kill_tri(M, M->cav[k]);
    for (int k = 0; k < M->nbnd; k++) {
        int a = M->bnd[k][0], b = M->bnd[k][1], n = M->bnd[k][3], owner = M->bnd[k][2];
        int t = new_tri(M);
        Tri* T = &M->t[t];
        T->v[0] = p; T->v[1] = a; T->v[2] = b;
        T->nb[0] = n; T->nb[1] = -1; T->nb[2] = -1;
        T->inside = oldin[k];
        if (n >= 0) {
            Tri* N = &M->t[n];
            for (int j = 0; j < 3; j++) if (N->nb[j] == owner) {
                int na = N->v[(j + 1) % 3], nb2 = N->v[(j + 2) % 3];
                if ((na == b && nb2 == a) || (na == a && nb2 == b)) { N->nb[j] = t; break; }
            }
        }
        M->vtri[a] = t; M->vtri[b] = t;
        M->newt[M->nnew++] = t;
    }
    free(oldin);
    M->vtri[p] = M->newt[0];
    /* link the fan: T = (p, a, b); across edge (b, p) (opposite a, nb[1]) is
     * the fan triangle starting at b; across (p, a) (nb[2]) the one ending at a */
    for (int k = 0; k < M->nnew; k++) {
        Tri* T = &M->t[M->newt[k]];
        for (int j = 0; j < M->nnew; j++) {
            if (j == k) continue;
            Tri* U = &M->t[M->newt[j]];
            if (U->v[1] == T->v[2]) T->nb[1] = M->newt[j];
            if (U->v[2] == T->v[1]) T->nb[2] = M->newt[j];
        }
    }
    /* inside flags from segment sides */
    for (int k = 0; k < M->nnew; k++) {
        int t = M->newt[k];
        Tri* T = &M->t[t];
        for (int i = 0; i < 3; i++) {
            int a = T->v[(i + 1) % 3], b = T->v[(i + 2) % 3];
            int s = seg_of(M, a, b);
            if (s >= 0) {
                int w = T->v[i];
                T->inside = orient(M, M->s[s].a, M->s[s].b, w) > 0;
                break;
            }
        }
    }
}

/* Locate and insert (x, y) as a new point. In refine mode the point is
 * rejected (returns -1) when it encroaches segments, which are queued. */
static int insert_point(Mesh* M, int start, double x, double y, bool corner, bool refine) {
    int hit = -1;
    int t = locate(M, start, x, y, refine, &hit);
    if (t == -2) { M->s[hit].force = 1; push_segq(M, hit); return -1; }
    if (t < 0) return -1;
    if (refine && !M->t[t].inside) return -1;
    if (!cavity(M, t, x, y)) return -1;
    if (refine) {
        bool enc = false;
        for (int k = 0; k < M->nbnd; k++) {
            int a = M->bnd[k][0], b = M->bnd[k][1];
            int s = seg_of(M, a, b);
            if (s < 0) continue;
            double dot = (M->px[a] - x) * (M->px[b] - x) + (M->py[a] - y) * (M->py[b] - y);
            if (dot <= 0) { M->s[s].force = 1; push_segq(M, s); enc = true; }
        }
        if (enc) return -1;
    }
    int p = add_point(M, x, y, corner);
    commit(M, p);
    return p;
}

/* The triangle containing edge (a, b) and the index of that edge, by
 * rotation around a. */
static int find_edge(const Mesh* M, int a, int b, int* idx) {
    int t0 = M->vtri[a];
    if (t0 < 0 || !M->t[t0].alive) return -1;
    for (int dir = 0; dir < 2; dir++) {
        int t = t0, guard = 0;
        while (t >= 0 && guard++ < 10000) {
            const Tri* T = &M->t[t];
            int k = T->v[0] == a ? 0 : T->v[1] == a ? 1 : T->v[2] == a ? 2 : -1;
            if (k < 0) break;
            for (int i = 0; i < 3; i++) {
                int ea = T->v[(i + 1) % 3], eb = T->v[(i + 2) % 3];
                if ((ea == a && eb == b) || (ea == b && eb == a)) { *idx = i; return t; }
            }
            t = dir == 0 ? T->nb[(k + 2) % 3] : T->nb[(k + 1) % 3];
            if (t == t0) break;
        }
    }
    return -1;
}

/* Is segment s missing from the triangulation, or encroached by an apex? */
static bool seg_needs_split(const Mesh* M, int s, double max_len) {
    int a = M->s[s].a, b = M->s[s].b, i;
    double dx = M->px[b] - M->px[a], dy = M->py[b] - M->py[a];
    if (max_len > 0 && dx * dx + dy * dy > max_len * max_len) return true;
    int t = find_edge(M, a, b, &i);
    if (t < 0) return true;
    int tri2[2] = { t, M->t[t].nb[i] };
    for (int k = 0; k < 2; k++) {
        int u = tri2[k];
        if (u < 0) continue;
        int w = -1;
        for (int j = 0; j < 3; j++) {
            int v = M->t[u].v[j];
            if (v != a && v != b) w = v;
        }
        if (w < 3) continue;                     /* the enclosing triangle */
        double dot = (M->px[a] - M->px[w]) * (M->px[b] - M->px[w])
                   + (M->py[a] - M->py[w]) * (M->py[b] - M->py[w]);
        if (dot <= 0) return true;
    }
    return false;
}

static void split_seg(Mesh* M, int s) {
    int a = M->s[s].a, b = M->s[s].b, arc = M->s[s].arc;
    double x, y;
    if (arc >= 0) {
        const FemArc* A = &M->g->arc[arc];
        double ta = atan2(M->py[a] - A->cy, M->px[a] - A->cx);
        double tb = atan2(M->py[b] - A->cy, M->px[b] - A->cx);
        double d = tb - ta;
        while (d > M_PI) d -= 2 * M_PI;
        while (d <= -M_PI) d += 2 * M_PI;
        x = A->cx + A->r * cos(ta + 0.5 * d);
        y = A->cy + A->r * sin(ta + 0.5 * d);
    } else {
        double dx = M->px[b] - M->px[a], dy = M->py[b] - M->py[a];
        double L = sqrt(dx * dx + dy * dy), f = 0.5;
        bool ca = M->corner[a], cb = M->corner[b];
        if (ca != cb && L > 0) {
            /* concentric shells: a power-of-two distance from the corner */
            double d = ldexp(1.0, (int)floor(log2(L / 2) + 0.5));
            f = d / L;
            if (cb) f = 1.0 - f;
        }
        x = M->px[a] + f * dx; y = M->py[a] + f * dy;
    }
    M->s[s].alive = 0;
    hash_del(M, a, b);
    int s1 = M->ns, s2 = M->ns + 1;
    int start = M->vtri[a];
    /* the new segments are registered after the point exists */
    int p = add_point(M, x, y, false);
    int hit = -1;
    int t = locate(M, start, x, y, false, &hit);
    if (t < 0 || !cavity(M, t, x, y)) {
        /* cannot place the point: restore the segment */
        M->np--;
        M->s[s].alive = 1; hash_put(M, a, b, s);
        M->s[s].queued = 0;
        return;
    }
    add_seg(M, a, p, arc);
    add_seg(M, p, b, arc);
    commit(M, p);
    (void)s1; (void)s2;
    push_segq(M, M->ns - 2);
    push_segq(M, M->ns - 1);
}

/* After an insertion: queue the segments the new triangles may encroach,
 * and (refinement phase) the new triangles themselves. */
static void after_insert(Mesh* M, bool queue_tris) {
    for (int k = 0; k < M->nnew; k++) {
        int t = M->newt[k];
        const Tri* T = &M->t[t];
        for (int i = 0; i < 3; i++) {
            int s = seg_of(M, T->v[(i + 1) % 3], T->v[(i + 2) % 3]);
            if (s >= 0 && !M->s[s].queued) push_segq(M, s);
        }
        if (queue_tris) push_triq(M, t);
    }
}

static void process_segments(Mesh* M, double max_len, int64_t maxv) {
    while (M->nsegq > 0) {
        int s = M->segq[--M->nsegq];
        M->s[s].queued = 0;
        if (!M->s[s].alive) continue;
        if (M->np >= maxv) continue;
        /* forced: a rejected circumcentre encroaches it */
        bool forced = M->s[s].force;
        M->s[s].force = 0;
        if (forced || seg_needs_split(M, s, max_len)) {
            split_seg(M, s);
            if (!M->s[s].alive) after_insert(M, true);
        }
    }
}

static void tri_quality(const Mesh* M, int t, double* area, double* ratio, double* lmin,
                        double* cx, double* cy) {
    const int* v = M->t[t].v;
    double ax = M->px[v[0]], ay = M->py[v[0]];
    double bx = M->px[v[1]] - ax, by = M->py[v[1]] - ay;
    double qx = M->px[v[2]] - ax, qy = M->py[v[2]] - ay;
    double d = 2.0 * (bx * qy - by * qx);
    double b2 = bx * bx + by * by, q2 = qx * qx + qy * qy;
    double ux = (qy * b2 - by * q2) / d, uy = (bx * q2 - qx * b2) / d;
    *cx = ax + ux; *cy = ay + uy;
    *area = 0.25 * d;
    double ex = M->px[v[2]] - M->px[v[1]], ey = M->py[v[2]] - M->py[v[1]];
    double l1 = b2, l2 = q2, l3 = ex * ex + ey * ey;
    double m = l1 < l2 ? l1 : l2; if (l3 < m) m = l3;
    *lmin = sqrt(m);
    *ratio = sqrt(ux * ux + uy * uy) / (*lmin > 0 ? *lmin : 1e-300);
}

/* ------------------------------------------------------------------ */
/*  Re-entrant corner grading                                          */
/* ------------------------------------------------------------------ */
/* At a corner of interior angle omega > pi the solution of an elliptic
 * problem behaves like r^lambda, lambda = pi/omega < 1 (pi/(2 omega) where a
 * Dirichlet edge meets a natural one, < 1 already for omega > pi/2), and a
 * quasi-uniform
 * mesh converges only as h^(2 lambda) in the eigenvalue (h^(4/3) at an
 * L-shape's corner) instead of h^4. Grading the mesh as
 *     h(r) ~ h0 (r/R)^(1 - lambda/2)
 * within R (the shorter adjacent edge) restores the P2 rate (Babuska): the
 * area limit is min(maxA, R^2/4) (r/R)^(2 - lambda), floored at 1e-8 of it.
 * At a mixed-condition corner R is half that: on SIAM Problem 14's pentagon
 * the full radius cost 1.75x the nodes (2.4 s against 1.1 s a solve) for a
 * 2.5x smaller error (1.1e-9 against 2.7e-9), not worth it there; at a
 * re-entrant corner of a Koch snowflake the full radius halves the error. */
typedef struct { double x, y, R, e; } Corner;
typedef struct {
    Corner* c; int nc;
    double cell, x0, y0; int nx, ny;
    int* start; int* list;
} Grading;

static void grading_init(Grading* G, const FemGeom* g, bool on) {
    memset(G, 0, sizeof *G);
    if (!on) return;
    int cap = 16;
    G->c = malloc(sizeof(Corner) * (size_t)cap);
    for (int k = 0; k < g->nloop; k++) {
        int s = g->loop_start[k], n = g->loop_n[k];
        for (int i = 0; i < n; i++) {
            int gi = s + i, gp = s + (i + n - 1) % n, gq = s + (i + 1) % n;
            if (g->edge_arc[gp] >= 0 || g->edge_arc[gi] >= 0) continue;   /* straight corners */
            double ax = g->xy[2 * gi] - g->xy[2 * gp], ay = g->xy[2 * gi + 1] - g->xy[2 * gp + 1];
            double bx = g->xy[2 * gq] - g->xy[2 * gi], by = g->xy[2 * gq + 1] - g->xy[2 * gi + 1];
            double cr = ax * by - ay * bx, dt = ax * bx + ay * by;
            double omega = M_PI - atan2(cr, dt);          /* region on the left */
            bool mixed = g->edge_dir && g->edge_dir[gp] != g->edge_dir[gi];
            double lam = mixed ? M_PI / (2 * omega) : M_PI / omega;
            if (!(lam < 1.0 - 1e-9)) continue;           /* u in H^2 there: no grading */
            double la = sqrt(ax * ax + ay * ay), lb = sqrt(bx * bx + by * by);
            if (G->nc == cap) { cap *= 2; G->c = realloc(G->c, sizeof(Corner) * (size_t)cap); }
            Corner* C = &G->c[G->nc++];
            C->x = g->xy[2 * gi]; C->y = g->xy[2 * gi + 1];
            C->R = (mixed ? 0.5 : 1.0) * (la < lb ? la : lb);
            C->e = 2.0 - lam;
        }
    }
    if (G->nc == 0) return;
    double xmin = INFINITY, xmax = -INFINITY, ymin = INFINITY, ymax = -INFINITY, rmax = 0;
    for (int i = 0; i < G->nc; i++) {
        const Corner* C = &G->c[i];
        if (C->x < xmin) xmin = C->x;
        if (C->x > xmax) xmax = C->x;
        if (C->y < ymin) ymin = C->y;
        if (C->y > ymax) ymax = C->y;
        if (C->R > rmax) rmax = C->R;
    }
    G->cell = rmax > 0 ? rmax : 1.0;
    G->x0 = xmin - G->cell; G->y0 = ymin - G->cell;
    G->nx = (int)((xmax - G->x0) / G->cell) + 2;
    G->ny = (int)((ymax - G->y0) / G->cell) + 2;
    if ((double)G->nx * G->ny > 4e6) {        /* corners too spread for a grid: scan */
        G->nx = G->ny = 0;
        return;
    }
    int ncell = G->nx * G->ny;
    G->start = calloc((size_t)ncell + 1, sizeof(int));
    G->list = malloc(sizeof(int) * (size_t)G->nc);
    for (int i = 0; i < G->nc; i++) {
        int cx = (int)((G->c[i].x - G->x0) / G->cell), cy = (int)((G->c[i].y - G->y0) / G->cell);
        G->start[cy * G->nx + cx + 1]++;
    }
    for (int c = 0; c < ncell; c++) G->start[c + 1] += G->start[c];
    int* fill = calloc((size_t)ncell, sizeof(int));
    for (int i = 0; i < G->nc; i++) {
        int cx = (int)((G->c[i].x - G->x0) / G->cell), cy = (int)((G->c[i].y - G->y0) / G->cell);
        int c = cy * G->nx + cx;
        G->list[G->start[c] + fill[c]++] = i;
    }
    free(fill);
}

static void grading_free(Grading* G) { free(G->c); free(G->start); free(G->list); }

/* The corner's area limit at (x, y): the local size at r = R is the
 * smaller of h0 and R/2 (a corner between short edges -- a fractal-like
 * boundary -- is graded relative to its own scale, not the global one). */
static double corner_limit(const Corner* C, double maxA, double x, double y) {
    double dx = x - C->x, dy = y - C->y, r = sqrt(dx * dx + dy * dy);
    if (r >= C->R) return INFINITY;
    double base = 0.25 * C->R * C->R;
    if (maxA < base) base = maxA;
    double f = pow(r / C->R, C->e);
    return base * (f < 1e-8 ? 1e-8 : f);
}

/* The area limit (<= maxA) at (x, y). */
static double grading_limit(const Grading* G, double maxA, double x, double y) {
    double f = maxA;
    if (G->nx == 0) {
        for (int i = 0; i < G->nc; i++) { double v = corner_limit(&G->c[i], maxA, x, y); if (v < f) f = v; }
        return f;
    }
    int cx = (int)floor((x - G->x0) / G->cell), cy = (int)floor((y - G->y0) / G->cell);
    for (int j = cy - 1; j <= cy + 1; j++) {
        if (j < 0 || j >= G->ny) continue;
        for (int i = cx - 1; i <= cx + 1; i++) {
            if (i < 0 || i >= G->nx) continue;
            int c = j * G->nx + i;
            for (int p = G->start[c]; p < G->start[c + 1]; p++) {
                double v = corner_limit(&G->c[G->list[p]], maxA, x, y);
                if (v < f) f = v;
            }
        }
    }
    return f;
}

static double cx_of(const Mesh* M, int t) {
    const int* v = M->t[t].v;
    return (M->px[v[0]] + M->px[v[1]] + M->px[v[2]]) / 3.0;
}
static double cy_of(const Mesh* M, int t) {
    const int* v = M->t[t].v;
    return (M->py[v[0]] + M->py[v[1]] + M->py[v[2]]) / 3.0;
}

/* ------------------------------------------------------------------ */
/*  Driver                                                             */
/* ------------------------------------------------------------------ */

static void mesh_free_state(Mesh* M) {
    free(M->px); free(M->py); free(M->corner); free(M->vtri);
    free(M->t); free(M->freel); free(M->stamp); free(M->s);
    free(M->hkey); free(M->hval); free(M->cav); free(M->bnd); free(M->newt);
    free(M->segq); free(M->triq);
}

void fem_mesh_free(FemMesh* m) {
    free(m->x); free(m->tri); free(m->curved); free(m->bedge); free(m->bnode);
    memset(m, 0, sizeof *m);
}

/* ------------------------------------------------------------------ */
/*  Uniform refinement                                                 */
/* ------------------------------------------------------------------ */

/* The parent element's quadratic map at reference point (s, r). */
static void p2_map(const double* X, const double* Y, double s, double r, double* x, double* y) {
    double l0 = 1.0 - s - r, l1 = s, l2 = r;
    double N[6] = { l0 * (2 * l0 - 1), l1 * (2 * l1 - 1), l2 * (2 * l2 - 1),
                    4 * l0 * l1, 4 * l1 * l2, 4 * l2 * l0 };
    *x = 0; *y = 0;
    for (int k = 0; k < 6; k++) { *x += N[k] * X[k]; *y += N[k] * Y[k]; }
}

/* The point of the circle through a, m, b (a boundary edge of the mesh with
 * its midside node m) halfway in angle between a and m; false when the three
 * are collinear (a straight edge). */
static bool arc_mid(double ax, double ay, double mx, double my, double bx, double by,
                    double* x, double* y) {
    double d = 2.0 * (ax * (my - by) + mx * (by - ay) + bx * (ay - my));
    double L2 = (bx - ax) * (bx - ax) + (by - ay) * (by - ay);
    double dev = fabs((mx - ax) * (by - ay) - (my - ay) * (bx - ax));   /* |cross| */
    if (!(L2 > 0) || dev <= 1e-12 * L2) return false;
    double a2 = ax * ax + ay * ay, m2 = mx * mx + my * my, b2 = bx * bx + by * by;
    double cx = (a2 * (my - by) + m2 * (by - ay) + b2 * (ay - my)) / d;
    double cy = (a2 * (bx - mx) + m2 * (ax - bx) + b2 * (mx - ax)) / d;
    double r = sqrt((ax - cx) * (ax - cx) + (ay - cy) * (ay - cy));
    double ta = atan2(ay - cy, ax - cx), tm = atan2(my - cy, mx - cx);
    double dt = tm - ta;
    while (dt > M_PI) dt -= 2 * M_PI;
    while (dt <= -M_PI) dt += 2 * M_PI;
    *x = cx + r * cos(ta + 0.5 * dt);
    *y = cy + r * sin(ta + 0.5 * dt);
    return true;
}

bool fem_mesh_refine(const FemMesh* in, FemMesh* out) {
    memset(out, 0, sizeof *out);
    int64_t ne = 4 * in->ne, nvert = in->nn;
    int64_t cap = nvert + 3 * ne;
    double* X = malloc(sizeof(double) * 2 * (size_t)cap);
    int64_t* tri = malloc(sizeof(int64_t) * 6 * (size_t)ne);
    unsigned char* curved = calloc((size_t)ne, 1);
    int64_t* bedge = malloc(sizeof(int64_t) * 3 * (size_t)(2 * in->nb + 1));
    if (!X || !tri || !curved || !bedge) { free(X); free(tri); free(curved); free(bedge); return false; }
    memcpy(X, in->x, sizeof(double) * 2 * (size_t)in->nn);
    /* boundary edges of the parent: their child edges' midside nodes go on
     * the circle through the parent edge's three nodes */
    int64_t bcap = 1; while (bcap < 4 * in->nb + 4) bcap <<= 1;
    int64_t* bk = malloc(sizeof(int64_t) * (size_t)bcap);
    for (int64_t i = 0; i < bcap; i++) bk[i] = -1;
    for (int64_t b = 0; b < in->nb; b++) {
        int64_t key = ekey((int)in->bedge[3 * b], (int)in->bedge[3 * b + 1]);
        uint64_t h = hmix(key) & (uint64_t)(bcap - 1);
        while (bk[h] != -1 && bk[h] != key) h = (h + 1) & (uint64_t)(bcap - 1);
        bk[h] = key;
    }
    int64_t ecap = 1; while (ecap < 4 * (3 * ne)) ecap <<= 1;
    int64_t* ek = malloc(sizeof(int64_t) * (size_t)ecap);
    int64_t* ev = malloc(sizeof(int64_t) * (size_t)ecap);
    for (int64_t i = 0; i < ecap; i++) ek[i] = -1;
    /* reference coordinates of the parent's six nodes */
    static const double RS[6][2] = { {0, 0}, {1, 0}, {0, 1}, {0.5, 0}, {0.5, 0.5}, {0, 0.5} };
    /* children as parent-node indices, counter-clockwise */
    static const int CH[4][3] = { {0, 3, 5}, {3, 1, 4}, {5, 4, 2}, {3, 4, 5} };
    /* the parent edge (vertex pair) a half-edge lies on, -1 interior */
    int64_t nn = nvert, e2 = 0;
    for (int64_t e = 0; e < in->ne; e++) {
        const int64_t* t = in->tri + 6 * e;
        double PX[6], PY[6];
        for (int k = 0; k < 6; k++) { PX[k] = in->x[2 * t[k]]; PY[k] = in->x[2 * t[k] + 1]; }
        for (int c = 0; c < 4; c++, e2++) {
            int64_t* ct = tri + 6 * e2;
            for (int k = 0; k < 3; k++) ct[k] = t[CH[c][k]];
            bool cv = false;
            for (int k = 0; k < 3; k++) {
                int pa = CH[c][k], pb = CH[c][(k + 1) % 3];
                int64_t a = t[pa], b = t[pb];
                int64_t key = ekey((int)a, (int)b);
                uint64_t h = hmix(key) & (uint64_t)(ecap - 1);
                while (ek[h] != -1 && ek[h] != key) h = (h + 1) & (uint64_t)(ecap - 1);
                int64_t node;
                if (ek[h] == key) node = ev[h];
                else {
                    node = nn++;
                    ek[h] = key; ev[h] = node;
                    double x, y;
                    p2_map(PX, PY, 0.5 * (RS[pa][0] + RS[pb][0]), 0.5 * (RS[pa][1] + RS[pb][1]), &x, &y);
                    /* a half of a parent boundary edge: (vertex, midside) */
                    int pv = pa < 3 ? pa : pb, pm = pa < 3 ? pb : pa;
                    if (pv < 3 && pm >= 3) {
                        int q = pm - 3;                    /* parent edge (q, q+1) */
                        int64_t va = t[q], vb = t[(q + 1) % 3];
                        int64_t pk = ekey((int)va, (int)vb);
                        uint64_t hb = hmix(pk) & (uint64_t)(bcap - 1);
                        while (bk[hb] != -1 && bk[hb] != pk) hb = (hb + 1) & (uint64_t)(bcap - 1);
                        if (bk[hb] == pk) {
                            /* the circle through the edge's three nodes, from vertex pv */
                            int64_t vo = t[pv] == va ? vb : va;
                            double ax = in->x[2 * t[pv]], ay = in->x[2 * t[pv] + 1];
                            double mx = in->x[2 * t[pm]], my = in->x[2 * t[pm] + 1];
                            double bx = in->x[2 * vo], by = in->x[2 * vo + 1];
                            double qx, qy;
                            if (arc_mid(ax, ay, mx, my, bx, by, &qx, &qy)) { x = qx; y = qy; }
                        }
                    }
                    X[2 * node] = x; X[2 * node + 1] = y;
                }
                ct[3 + k] = node;
                double ax = X[2 * a], ay = X[2 * a + 1], bx = X[2 * b], by = X[2 * b + 1];
                double dx = X[2 * node] - 0.5 * (ax + bx), dy = X[2 * node + 1] - 0.5 * (ay + by);
                double L2 = (bx - ax) * (bx - ax) + (by - ay) * (by - ay);
                if (dx * dx + dy * dy > 1e-24 * L2) cv = true;
            }
            curved[e2] = cv;
        }
    }
    /* boundary edges: each parent edge (a, b, m) -> (a, m) and (m, b) */
    int64_t nb = 0;
    for (int64_t b = 0; b < in->nb; b++) {
        int64_t va = in->bedge[3 * b], vb = in->bedge[3 * b + 1], vm = in->bedge[3 * b + 2];
        int64_t pr[2][2] = { { va, vm }, { vm, vb } };
        for (int j = 0; j < 2; j++) {
            int64_t key = ekey((int)pr[j][0], (int)pr[j][1]);
            uint64_t h = hmix(key) & (uint64_t)(ecap - 1);
            while (ek[h] != -1 && ek[h] != key) h = (h + 1) & (uint64_t)(ecap - 1);
            if (ek[h] != key) continue;
            bedge[3 * nb] = pr[j][0]; bedge[3 * nb + 1] = pr[j][1]; bedge[3 * nb + 2] = ev[h];
            nb++;
        }
    }
    free(ek); free(ev); free(bk);
    out->nn = nn; out->ne = ne; out->nb = nb; out->nvert = nvert;
    out->x = realloc(X, sizeof(double) * 2 * (size_t)nn);
    out->tri = tri; out->curved = curved; out->bedge = bedge;
    out->bnode = calloc((size_t)nn, 1);
    for (int64_t b = 0; b < nb; b++)
        for (int k = 0; k < 3; k++) out->bnode[bedge[3 * b + k]] = 1;
    memcpy(out->bbox, in->bbox, sizeof out->bbox);
    return true;
}

bool fem_mesh_build(const FemGeom* g, const FemMeshOpts* o, FemMesh* out, const char** why) {
    memset(out, 0, sizeof *out);
    *why = NULL;
    if (g->nv < 3) { *why = "the region has no boundary"; return false; }
    Mesh M; memset(&M, 0, sizeof M);
    M.g = g; M.rng = 12345u;

    double xmin = INFINITY, xmax = -INFINITY, ymin = INFINITY, ymax = -INFINITY;
    for (int64_t i = 0; i < g->nv; i++) {
        double x = g->xy[2 * i], y = g->xy[2 * i + 1];
        if (!isfinite(x) || !isfinite(y)) { *why = "the region has a non-numeric vertex"; return false; }
        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;
        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
    }
    for (int k = 0; k < g->narc; k++) {        /* arcs may bulge past their vertices */
        const FemArc* A = &g->arc[k];
        if (A->cx - A->r < xmin) xmin = A->cx - A->r;
        if (A->cx + A->r > xmax) xmax = A->cx + A->r;
        if (A->cy - A->r < ymin) ymin = A->cy - A->r;
        if (A->cy + A->r > ymax) ymax = A->cy + A->r;
    }
    double w = xmax - xmin, h = ymax - ymin, sz = w > h ? w : h;
    if (!(sz > 0)) { *why = "the region is degenerate"; return false; }
    double mx = 0.5 * (xmin + xmax), my = 0.5 * (ymin + ymax);
    /* enclosing triangle */
    add_point(&M, mx - 20 * sz, my - 10 * sz, false);
    add_point(&M, mx + 20 * sz, my - 10 * sz, false);
    add_point(&M, mx, my + 20 * sz, false);
    int t0 = new_tri(&M);
    M.t[t0].v[0] = 0; M.t[t0].v[1] = 1; M.t[t0].v[2] = 2;
    M.t[t0].nb[0] = M.t[t0].nb[1] = M.t[t0].nb[2] = -1;
    M.vtri[0] = M.vtri[1] = M.vtri[2] = t0;

    /* 1. input vertices */
    int* vid = malloc(sizeof(int) * (size_t)g->nv);
    for (int k = 0; k < g->nloop; k++) {
        int s = g->loop_start[k], n = g->loop_n[k];
        for (int i = 0; i < n; i++) {
            int gi = s + i, prev = s + (i + n - 1) % n;
            bool corner = !(g->edge_arc[gi] >= 0 && g->edge_arc[gi] == g->edge_arc[prev]);
            int start = M.np > 3 ? M.vtri[M.np - 1] : t0;
            int p = insert_point(&M, start, g->xy[2 * gi], g->xy[2 * gi + 1], corner, false);
            if (p < 0) {
                /* a repeated vertex: reuse the existing one */
                int found = -1;
                for (int q = 3; q < M.np && found < 0; q++)
                    if (M.px[q] == g->xy[2 * gi] && M.py[q] == g->xy[2 * gi + 1]) found = q;
                if (found < 0) { free(vid); mesh_free_state(&M); *why = "a boundary vertex could not be inserted"; return false; }
                p = found;
            }
            vid[gi] = p;
        }
    }
    /* 2. segments */
    for (int k = 0; k < g->nloop; k++) {
        int s = g->loop_start[k], n = g->loop_n[k];
        for (int i = 0; i < n; i++) {
            int a = vid[s + i], b = vid[s + (i + 1) % n];
            if (a == b) continue;
            add_seg(&M, a, b, g->edge_arc[s + i]);
        }
    }
    free(vid);
    for (int s = 0; s < M.ns; s++) push_segq(&M, s);
    int64_t maxv = o->max_vertices > 0 ? o->max_vertices : 4000000;
    process_segments(&M, o->max_bedge, maxv);

    /* inside flags: from the segment sides, flooded across other edges */
    {
        for (int t = 0; t < M.nt; t++) M.t[t].inside = 2;      /* unknown */
        int* st = malloc(sizeof(int) * (size_t)(M.nt + 1)); int nst = 0;
        for (int t = 0; t < M.nt; t++) {
            if (!M.t[t].alive) continue;
            for (int i = 0; i < 3; i++) {
                int s = seg_of(&M, M.t[t].v[(i + 1) % 3], M.t[t].v[(i + 2) % 3]);
                if (s >= 0) {
                    M.t[t].inside = orient(&M, M.s[s].a, M.s[s].b, M.t[t].v[i]) > 0;
                    st[nst++] = t;
                    break;
                }
            }
        }
        while (nst > 0) {
            int t = st[--nst];
            for (int i = 0; i < 3; i++) {
                int n = M.t[t].nb[i];
                if (n < 0 || M.t[n].inside != 2) continue;
                if (seg_of(&M, M.t[t].v[(i + 1) % 3], M.t[t].v[(i + 2) % 3]) >= 0) continue;
                M.t[n].inside = M.t[t].inside;
                st[nst++] = n;
            }
        }
        free(st);
        for (int t = 0; t < M.nt; t++) if (M.t[t].inside == 2) M.t[t].inside = 0;
    }

    /* 3. refinement */
    double maxA = o->max_area > 0 ? o->max_area : INFINITY;
    double hmin = isfinite(maxA) ? 1e-3 * sqrt(maxA) : 1e-6 * sz;
    Grading G; grading_init(&G, g, isfinite(maxA) && o->grade_corners);
    const double B = 1.4142135623730951;
    for (int t = 0; t < M.nt; t++) if (M.t[t].alive && M.t[t].inside) push_triq(&M, t);
    while (M.tqhead < M.ntriq) {
        int t = M.triq[M.tqhead++];
        if (!M.t[t].alive || !M.t[t].inside || M.t[t].noref) continue;
        if (M.np >= maxv) break;
        double area, ratio, lmin, cx, cy;
        tri_quality(&M, t, &area, &ratio, &lmin, &cx, &cy);
        bool big = area > maxA;
        /* the minimum edge for quality refinement scales with the graded
         * size, so graded triangles keep the angle bound */
        double hloc = hmin;
        if (G.nc > 0) {
            double lim = grading_limit(&G, maxA, cx_of(&M, t), cy_of(&M, t));
            if (!big) big = area > lim;
            if (lim < maxA) hloc = 1e-3 * sqrt(lim);
        }
        bool skinny = ratio > B && lmin > hloc;
        if (!big && !skinny) continue;
        int p = insert_point(&M, t, cx, cy, false, true);
        if (p >= 0) {
            after_insert(&M, true);
        } else if (M.nsegq > 0) {
            process_segments(&M, o->max_bedge, maxv);
            if (M.t[t].alive) push_triq(&M, t);
        } else {
            M.t[t].noref = 1;       /* cannot be improved (rounding) */
        }
        if (M.nsegq > 0) process_segments(&M, o->max_bedge, maxv);
    }

    grading_free(&G);

    /* 4. extract the inside triangles and build the P2 mesh */
    int* newid = malloc(sizeof(int) * (size_t)M.np);
    for (int i = 0; i < M.np; i++) newid[i] = -1;
    int64_t ne = 0, nvert = 0;
    for (int t = 0; t < M.nt; t++) {
        if (!M.t[t].alive || !M.t[t].inside) continue;
        if (M.t[t].v[0] < 3 || M.t[t].v[1] < 3 || M.t[t].v[2] < 3) continue;
        ne++;
        for (int k = 0; k < 3; k++) if (newid[M.t[t].v[k]] < 0) newid[M.t[t].v[k]] = (int)nvert++;
    }
    if (ne == 0) { free(newid); mesh_free_state(&M); *why = "the region has no interior"; return false; }
    /* edges: hash (lo, hi) -> midside node */
    int64_t ecap = 1; while (ecap < 4 * (3 * ne)) ecap <<= 1;
    int64_t* ek = malloc(sizeof(int64_t) * (size_t)ecap);
    int64_t* ev = malloc(sizeof(int64_t) * (size_t)ecap);
    for (int64_t i = 0; i < ecap; i++) ek[i] = -1;
    int64_t nn = nvert, cap = nvert + 3 * ne;
    double* X = malloc(sizeof(double) * 2 * (size_t)cap);
    for (int i = 0; i < M.np; i++)
        if (newid[i] >= 0) { X[2 * newid[i]] = M.px[i]; X[2 * newid[i] + 1] = M.py[i]; }
    int64_t* tri = malloc(sizeof(int64_t) * 6 * (size_t)ne);
    unsigned char* curved = calloc((size_t)ne, 1);
    int64_t* bedge = malloc(sizeof(int64_t) * 3 * (size_t)(3 * ne));
    int64_t nb = 0, e = 0;
    for (int t = 0; t < M.nt; t++) {
        if (!M.t[t].alive || !M.t[t].inside) continue;
        if (M.t[t].v[0] < 3 || M.t[t].v[1] < 3 || M.t[t].v[2] < 3) continue;
        const int* v = M.t[t].v;
        for (int k = 0; k < 3; k++) tri[6 * e + k] = newid[v[k]];
        /* midside node of edge (v[k], v[k+1]) is entry 3 + k */
        for (int k = 0; k < 3; k++) {
            int a = v[k], b = v[(k + 1) % 3];
            int64_t key = ekey(a, b);
            uint64_t hi = hmix(key) & (uint64_t)(ecap - 1);
            while (ek[hi] != -1 && ek[hi] != key) hi = (hi + 1) & (uint64_t)(ecap - 1);
            int64_t node;
            if (ek[hi] == key) node = ev[hi];
            else {
                node = nn++;
                ek[hi] = key; ev[hi] = node;
                double x = 0.5 * (M.px[a] + M.px[b]), y = 0.5 * (M.py[a] + M.py[b]);
                int s = seg_of(&M, a, b);
                if (s >= 0 && M.s[s].arc >= 0) {
                    const FemArc* A = &g->arc[M.s[s].arc];
                    double ta = atan2(M.py[a] - A->cy, M.px[a] - A->cx);
                    double tb = atan2(M.py[b] - A->cy, M.px[b] - A->cx);
                    double d = tb - ta;
                    while (d > M_PI) d -= 2 * M_PI;
                    while (d <= -M_PI) d += 2 * M_PI;
                    x = A->cx + A->r * cos(ta + 0.5 * d);
                    y = A->cy + A->r * sin(ta + 0.5 * d);
                }
                X[2 * node] = x; X[2 * node + 1] = y;
            }
            tri[6 * e + 3 + k] = node;
            int s = seg_of(&M, a, b);
            if (s >= 0) {
                if (M.s[s].arc >= 0) curved[e] = 1;
                bedge[3 * nb] = newid[M.s[s].a]; bedge[3 * nb + 1] = newid[M.s[s].b];
                bedge[3 * nb + 2] = node;
                nb++;
            }
        }
        e++;
    }
    free(ek); free(ev); free(newid);
    out->nn = nn; out->ne = ne; out->nb = nb; out->nvert = nvert;
    out->x = realloc(X, sizeof(double) * 2 * (size_t)nn);
    out->tri = tri; out->curved = curved;
    out->bedge = bedge;
    out->bnode = calloc((size_t)nn, 1);
    for (int64_t b = 0; b < nb; b++)
        for (int k = 0; k < 3; k++) out->bnode[bedge[3 * b + k]] = 1;
    out->bbox[0] = INFINITY; out->bbox[1] = -INFINITY; out->bbox[2] = INFINITY; out->bbox[3] = -INFINITY;
    for (int64_t i = 0; i < nn; i++) {
        double x = out->x[2 * i], y = out->x[2 * i + 1];
        if (x < out->bbox[0]) out->bbox[0] = x;
        if (x > out->bbox[1]) out->bbox[1] = x;
        if (y < out->bbox[2]) out->bbox[2] = y;
        if (y > out->bbox[3]) out->bbox[3] = y;
    }
    mesh_free_state(&M);
    return true;
}
