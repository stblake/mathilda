/* FindClusters[list] / [list, n] / [list, UpTo[n]] -- partition a 1D numeric
 * list into clusters of nearby elements.
 *
 *   FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}]  -> {{1, 2, 3, 1}, {10, 12, 13}, {25}}
 *   FindClusters[{1, 2, 3, 5, 8, 9, 10}, 2]     -> {{1, 2, 3, 5}, {8, 9, 10}}
 *
 * WE DO NOT REPRODUCE MATHEMATICA'S OUTPUT, and this is not a shortfall to be
 * fixed later -- it is not reachable. Mathematica auto-selects a distance
 * function and preprocesses the data by rules it does not publish; its own
 * messages name the choice ("... method Agglomerate and distance function
 * HammingDistance ..."). The visible consequence is that even single-linkage
 * with an explicit count disagrees with the textbook algorithm:
 *
 *   FindClusters[{1,4,9,16,25,36}, 3, Method->"Agglomerate"] -> {{36},{1,4,9},{16,25}}
 *
 * on gaps 3,5,7,9,11 it cut 11 and 7 rather than the two largest, 11 and 9.
 * Its Automatic cluster count is likewise not a gap rule -- {1,4,9,16,25,36}
 * (linear gap growth) splits three ways while {1,2,4,8,16,32,64} (geometric,
 * far more extreme gaps) gives one. So: this file implements the textbook
 * algorithm for each named method with the semantics stated below, and the
 * acceptance table in tests/test_list.c is the specification.
 *
 * SEMANTICS
 *
 *   Distance          |a - b| on the line. Never via Abs -- see EXACTNESS.
 *   Cluster order     by first occurrence of any member in the input.
 *                     Mathematica's is unstable (Automatic and n=3 disagree on
 *                     the same input); first occurrence matches Gather.
 *   Element order     input order, within every cluster.
 *   Count             three modes, see FcCountMode. `n` is capped at the
 *                     distinct-value count, never exceeded.
 *
 * EXACTNESS
 *
 * The sorted order and the fixed-count gap selection are computed on the
 * ELEMENTS THEMSELVES via list_numeric_cmp, so exact input is ordered exactly:
 * FindClusters[{1/10^25, 2/10^25, 1}, 2] works. Nearest cannot do this,
 * because it routes its distance through Abs, and builtin_abs declines on a
 * rational with a bigint component (complex.c:418-421). A sorted 1D pass has no
 * need of Abs at all -- the gap between sorted neighbours is non-negative by
 * construction -- so the landmine is avoided rather than inherited. Rewriting a
 * gap as Abs[b - a] would reintroduce it, and an acceptance row fails if anyone
 * does.
 *
 * The Automatic threshold and the inherently-inexact methods (KMeans, the
 * density family, GaussianMixture, Spectral) work on a machine-double
 * projection instead, which is correct: a mean, a kernel and an eigenvector are
 * inexact by definition, and the Automatic threshold carries a fitted constant.
 *
 * COST
 *
 * Every method is O(n log n), dominated by the initial sort, except Spectral,
 * which is O(n^2) in memory and declines above FC_SPECTRAL_MAX_N. The density
 * family stays linear only because the array is sorted -- the same
 * neighbourhood queries in general dimension are the expensive part of those
 * algorithms.
 *
 * SCOPE
 *
 * 1D numeric lists. Mathematica also clusters symbolic elements as nominal
 * features (FindClusters[{1, a, 3}, 2] gives {{a}, {1, 3}}); we decline. The
 * rule forms, Association input, CriterionFunction, PerformanceGoal, Weights
 * and the FeatureX options are not implemented. */

#include "list_common.h"
#include "internal.h"
#include "find_clusters.h"
#include "distance.h"
#include "../ndarray.h"   /* is_ndarray, ndarray_to_nested_list */

#include <math.h>       /* isnan, fabs, exp, sqrt, log -- all C99 */

#ifndef M_PI            /* POSIX, not C99: glibc hides it under -std=c99 */
#define M_PI 3.14159265358979323846
#endif

/* The one tunable behind the Automatic cluster count: cut a sorted-adjacent gap
 * when it exceeds this multiple of the median gap. Fitted to reproduce
 * Mathematica on the cases probed while writing this (see the plan), NOT a
 * recovered internal index. Named here so the acceptance rows that move when it
 * is tuned are obvious. */
/* Integer, so the Automatic threshold (FC_GAP_FACTOR * median gap) stays exact
 * when the gaps are exact. The density defaults cast it to double. */
#define FC_GAP_FACTOR        3

#define FC_MAX_ITER          100    /* Lloyd / mean-shift / EM iteration cap */
#define FC_SPECTRAL_MAX_N    2000   /* Spectral builds an n x n matrix */

/* MeanShift and NeighborhoodContraction rescan the whole sample for every
 * point on every iteration, so they are Theta(FC_MAX_ITER * n^2) -- measured
 * 0.48 s at n = 2000 and 1.73 s at 4000, which extrapolates to roughly a
 * quarter of an hour at 10^5 inside a single uninterruptible builtin call.
 * Until they are reworked onto fc_eps_window (the flat kernel is exactly an
 * eps-window mean) they decline above this, the way Spectral already declines
 * above its own limit. Refusing is strictly better than appearing to hang. */
#define FC_SHIFT_MAX_N       4000
/* Vector input is quadratic in n (Prim, and the neighbourhood queries a sort
 * makes linear on a line), so it carries the same order of cap as the other
 * quadratic methods. */
#define FC_NDIM_MAX_N        2000
/* The machine-precision point builder is over two orders of magnitude faster
 * than the exact one -- 2000 2-D points went from 1.49 s to 6.8 ms -- so it
 * earns a correspondingly larger ceiling. Still quadratic: 20000 points is
 * roughly 0.7 s. */
#define FC_NDIM_MACHINE_MAX_N 20000
/* Lloyd in n dimensions is linear in n but bilinear in (n, k), so its ceiling is
 * on n * k * dim rather than on n -- see fc_lloyd_ndim. Admits 20000 points in ten
 * dimensions at any sensible k, and refuses k on the order of n at that size,
 * which would be 100 quadratic passes over the whole sample. */
#define FC_LLOYD_MAX_WORK    20000000u

/* Merge tolerance for the shift methods, in units of the data scale (one
 * median gap). Strictly 1.0 is too tight: on evenly spaced data the interior
 * points are stationary but the ones near each end drift inward slightly, so a
 * few adjacent spacings land just above one scale and the run fragments --
 * Range[40] came back as 9 clusters at a tolerance of exactly 1.0. Genuine
 * mode separation is many scales wide, so the slack costs nothing. */
#define FC_MERGE_SLACK       1.5

/* ------------------------------------------------------------------------- */
/* Methods and count modes                                                   */
/* ------------------------------------------------------------------------- */

/* Order must match the rows of FC_ALLOWED. */
typedef enum {
    FC_AGGLOMERATE = 0,
    FC_SPANNINGTREE,
    FC_KMEANS,
    FC_KMEDOIDS,
    FC_SPECTRAL,
    FC_DBSCAN,
    FC_GAUSSIANMIXTURE,
    FC_JARVISPATRICK,
    FC_MEANSHIFT,
    FC_NEIGHBORHOODCONTRACTION,
    FC_METHOD_COUNT
} FcMethod;

/* Three modes, not two. Mathematica gives each its own error class
 * (wrgmthundef / wrgmthbound / wrgmthdef) and the capability matrix differs
 * between them -- Spectral accepts UpTo[n] but rejects a bare n. */
typedef enum {
    FC_COUNT_AUTOMATIC = 0,   /* count omitted, or Automatic */
    FC_COUNT_BOUNDED   = 1,   /* UpTo[n]: at most n, MAY BE FEWER */
    FC_COUNT_FIXED     = 2    /* n: exactly n, capped at the distinct count */
} FcCountMode;

typedef struct { FcCountMode mode; size_t n; } FcCount;

/* [method][mode]. Transcribed from the allowed-lists in Mathematica's three
 * error messages, which are authoritative -- the documentation bullets place
 * Spectral in neither constraint list, implying it takes a fixed count, and the
 * runtime rejects it. */
static const bool FC_ALLOWED[FC_METHOD_COUNT][3] = {
    /*                            Automatic  Bounded  Fixed */
    /* Agglomerate             */ { true,    true,    true  },
    /* SpanningTree            */ { true,    true,    true  },
    /* KMeans                  */ { false,   true,    true  },
    /* KMedoids                */ { false,   true,    true  },
    /* Spectral                */ { true,    true,    false },
    /* DBSCAN                  */ { true,    false,   false },
    /* GaussianMixture         */ { true,    false,   false },
    /* JarvisPatrick           */ { true,    false,   false },
    /* MeanShift               */ { true,    false,   false },
    /* NeighborhoodContraction */ { true,    false,   false },
};

static const struct { const char* name; FcMethod m; } FC_METHOD_NAMES[] = {
    { "Agglomerate",             FC_AGGLOMERATE             },
    { "SpanningTree",            FC_SPANNINGTREE            },
    { "KMeans",                  FC_KMEANS                  },
    { "KMedoids",                FC_KMEDOIDS                },
    { "Spectral",                FC_SPECTRAL                },
    { "DBSCAN",                  FC_DBSCAN                  },
    { "GaussianMixture",         FC_GAUSSIANMIXTURE         },
    { "JarvisPatrick",           FC_JARVISPATRICK           },
    { "MeanShift",               FC_MEANSHIFT               },
    { "NeighborhoodContraction", FC_NEIGHBORHOODCONTRACTION },
};

/* Which metric the tree edges are weighted by.
 *
 * Euclidean and SquaredEuclidean are deliberately ONE case, not two. Squaring is
 * monotone on non-negatives, so it preserves the ranking of edges, and the
 * Automatic threshold compares against a multiple of the median -- where
 * d > 3 * median(d) if and only if d^2 > 9 * median(d^2). Both the ranking and
 * the threshold test are therefore identical under the two, so they always
 * produce the same partition, and collapsing them avoids taking a square root of
 * an exact value: Sqrt[2] is not rational, and an exact ordering that has to
 * compare irrationals is a different and much harder problem than this needs.
 *
 * Manhattan genuinely differs and is the one that had to be wired.
 *
 * FC_DIST_AUTOMATIC is kept distinct from FC_DIST_SQUAREDEUCLIDEAN even though
 * they behave identically, so that a future Automatic that chooses a metric from
 * the data does not have to be told apart from an explicit request. */
typedef enum {
    FC_DIST_AUTOMATIC,
    FC_DIST_EUCLIDEAN,
    FC_DIST_SQUAREDEUCLIDEAN,
    FC_DIST_MANHATTAN
} FcDistance;

/* True when this metric's weights are squared distances, which the Automatic
 * threshold factor has to match -- see FcData.thresh_factor. */
static bool fc_dist_is_squared(FcDistance dist) {
    return dist == FC_DIST_AUTOMATIC || dist == FC_DIST_EUCLIDEAN
        || dist == FC_DIST_SQUAREDEUCLIDEAN;
}

typedef struct {
    FcMethod method;
    double   radius;          bool radius_given;
    long     min_points;      bool min_points_given;
    long     neighbor_count;  bool neighbor_count_given;
    FcDistance dist;
} FcOpts;

/* ------------------------------------------------------------------------- */
/* Decoded input                                                             */
/* ------------------------------------------------------------------------- */

/* What the elements ARE, which decides both the distance and the code path.
 *
 *   SCALAR   real numbers. Sorted, so the spanning tree is the adjacency chain
 *            and edge weights are exact differences. The original path.
 *   POINT    a List or a colour, i.e. a compound expression whose arguments are
 *            numeric coordinates. Weights are exact SQUARED distances.
 *   SEQUENCE strings. Weights are exact integer edit distances. There are no
 *            coordinates at all, which is fine: the gap methods only ever ask
 *            for pairwise distances.
 *
 * Rational and Complex are also compound expressions, which is why POINT is
 * restricted to an explicit head list -- treating Rational[1, 2] as a 2-D point
 * would silently cluster 1/2 by its numerator and denominator. */
typedef enum { FC_KIND_SCALAR, FC_KIND_POINT, FC_KIND_SEQUENCE } FcKind;

/* Everything the methods share, computed once.
 *
 * `order` is the EXACT sorted permutation and `gap` the EXACT adjacent
 * differences; `val` is the machine projection, for the methods whose output is
 * inexact by definition. Keeping both is what lets the default gap method stay
 * exact while KMeans and the density family use doubles. */
typedef struct {
    Expr**  elem;      /* borrowed: the input elements, input order */
    size_t  n;
    size_t  dim;       /* component count for FC_KIND_POINT; 1 otherwise */
    FcKind  kind;      /* WHICH code path. Never infer this from `dim`:
                        * {{1}, {2}, {100}} is a list of 1-component points, so
                        * dim == 1 while the kind is POINT. Conflating the two
                        * sent 1-vectors down the scalar path, where the elements
                        * are still Lists and every arithmetic step threaded over
                        * them. `dim` only ever counts components. */
    size_t* order;     /* owned: sorted permutation of 0..n-1. dim == 1 ONLY --
                        * there is no total order on vectors, so this is NULL in
                        * higher dimensions and every consumer of it is a
                        * 1D-only path. */
    double* val;       /* owned: val[i] is the double projection of elem[i].
                        * dim == 1 only; see `coord` otherwise. */
    double* coord;     /* owned: row-major n x dim machine projection, dim > 1 */

    /* THE SPANNING TREE. `gap[j]` is edge j's exact weight and eu/ev its
     * endpoints. In 1D this is the sorted adjacency chain -- edge j is
     * (order[j], order[j+1]) weighted by their exact difference, i.e. literally
     * the adjacent-gap array this field used to be -- because the MST of points
     * on a line IS that chain. In higher dimensions it is a real MST over exact
     * SQUARED distances, built by Prim.
     *
     * Storing it as edges rather than gaps is what lets one implementation serve
     * both: single-linkage clustering is "cut the heaviest tree edges", and the
     * ranking and threshold code below never needs to know which case it is in.
     *
     * eu[j] is always the endpoint that was already in the tree when ev[j] joined
     * (in 1D, the earlier sorted position). Parent-before-child ordering is
     * relied upon by the equal-elements fold, which would not be transitive
     * without it. */
    Expr**  gap;       /* owned: exact weight of edge j */
    size_t* eu;        /* owned: parent endpoint of edge j */
    size_t* ev;        /* owned: child endpoint of edge j */
    size_t  n_gap;     /* n - 1, or 0 when n == 0 */

    /* Multiplier for the Automatic threshold, applied to the weights above.
     * FC_GAP_FACTOR when the weights are plain distances (1D differences), and
     * its SQUARE when they are squared distances (n-D), which is the same test:
     * squaring is monotone on non-negatives so the median commutes with it, and
     * d^2 > 9 * median(d^2) if and only if d > 3 * median(d). Getting this wrong
     * would silently apply a factor of sqrt(3) in higher dimensions. */
    size_t  thresh_factor;

    /* The metric the edge weights above were built with. Carried here rather than
     * passed alongside FcData because every distance consumer already takes a
     * FcData and nothing else needs the option struct. Meaningful for
     * FC_KIND_POINT only: on a line all four accepted metrics agree up to the
     * monotone transform the threshold factor cancels, and for strings the only
     * meaningful metric is edit distance. */
    FcDistance dist;

    size_t  n_distinct;
    /* owned: bnd[j] is true when edge j's endpoints hold EXACTLY different
     * values. Derived by comparing the elements themselves, never from the
     * double projection -- `2^60` and `2^60 + 1` are distinct but project to the
     * same double, and a boundary set computed in double space silently loses
     * them. Every consumer of "is this a real boundary" must read this, not
     * compare val[]. */
    bool*   bnd;
} FcData;

static void fc_data_free(FcData* d) {
    if (d->gap) {
        for (size_t j = 0; j < d->n_gap; j++) expr_free(d->gap[j]);
        free(d->gap);
    }
    free(d->order);
    free(d->val);
    free(d->coord);
    free(d->eu);
    free(d->ev);
    free(d->bnd);
    d->gap = NULL; d->order = NULL; d->val = NULL; d->bnd = NULL;
    d->coord = NULL; d->eu = NULL; d->ev = NULL;
}

/* Machine projection of a real number. Only for the inexact methods and the
 * Automatic threshold; ordering never goes through this. */
static double fc_to_double(Expr* e) {
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    if (e->type == EXPR_REAL)    return e->data.real;
    if (e->type == EXPR_BIGINT)  return mpz_get_d(e->data.bigint);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
#endif
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2) {
        double num = fc_to_double(e->data.function.args[0]);
        double den = fc_to_double(e->data.function.args[1]);
        if (den != 0.0) return num / den;
    }
    return 0.0;
}

/* ------------------------------------------------------------------------- */
/* Input shape                                                                */
/* ------------------------------------------------------------------------- */

/* Component c of element i, treating a scalar as a 1-vector. Borrowed. */
static Expr* fc_comp(const FcData* d, size_t i, size_t c) {
    Expr* e = d->elem[i];
    return (d->kind == FC_KIND_POINT) ? e->data.function.args[c] : e;
}

/* Heads whose arguments are coordinates. Colours are points in their own space,
 * which is all clustering needs; RGBColor[r, g, b] is the same shape as a
 * 3-vector. Deliberately a closed list -- see the FcKind comment for why any
 * compound head would be wrong. */
static bool fc_is_point_head(Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    const char* nm = h->data.symbol.name;
    return strcmp(nm, "List") == 0 || strcmp(nm, "RGBColor") == 0 ||
           strcmp(nm, "GrayLevel") == 0 || strcmp(nm, "Hue") == 0 ||
           strcmp(nm, "CMYKColor") == 0;
}

/* Exact equality of two whole points. Component-wise through the exact
 * comparator, so this is the vector generalisation of the 1D distinctness test
 * and inherits its refusal to decide through doubles. */
static bool fc_elem_equal(const FcData* d, size_t i, size_t j, bool* ok) {
    /* Strings have no numeric components; structural equality IS value equality
     * for them, and it cannot fail to decide. */
    if (d->kind == FC_KIND_SEQUENCE) {
        *ok = true;
        return expr_eq(d->elem[i], d->elem[j]);
    }
    for (size_t c = 0; c < d->dim; c++) {
        if (list_numeric_cmp(fc_comp(d, i, c), fc_comp(d, j, c), ok) != 0) return false;
        if (!*ok) return false;
    }
    return true;
}

/* Exact distance between two elements, by kind. Caller owns the result; NULL
 * means the pair is not comparable, which the caller turns into an unevaluated
 * call rather than a guess. */
static Expr* fc_pair_distance(const FcData* d, size_t i, size_t j) {
    if (d->kind == FC_KIND_SEQUENCE) return distance_edit(d->elem[i], d->elem[j]);
    if (d->dist == FC_DIST_MANHATTAN) return distance_manhattan(d->elem[i], d->elem[j]);
    return distance_squared_euclidean(d->elem[i], d->elem[j]);
}

/* Decide (n, dim) for the input, or fail.
 *
 * Two accepted shapes: every element a real scalar (dim 1, the original path),
 * or every element a List of the SAME length with every component a real number
 * (dim k). Everything else declines -- ragged rows, a mix of scalars and
 * vectors, depth over 2, a non-real component, a visible NDArray (not a List, so
 * it never reaches here as one). Declining rather than guessing follows the rest
 * of the file: a symbolic element is not silently dropped or reinterpreted as a
 * nominal feature. */
static bool fc_probe_shape(Expr** elem, size_t n, size_t* dim, FcKind* kind) {
    *dim = 1;

    /* Scalars first, so a Rational -- itself a compound expression -- is read as
     * the number it is rather than as a pair of coordinates. */
    if (list_real_number_q(elem[0])) {
        for (size_t i = 0; i < n; i++)
            if (!list_real_number_q(elem[i])) return false;
        *kind = FC_KIND_SCALAR;
        return true;
    }

    if (elem[0]->type == EXPR_STRING) {
        for (size_t i = 0; i < n; i++)
            if (elem[i]->type != EXPR_STRING) return false;
        *kind = FC_KIND_SEQUENCE;
        return true;
    }

    if (!fc_is_point_head(elem[0])) return false;
    Expr* head0 = elem[0]->data.function.head;
    size_t k = elem[0]->data.function.arg_count;
    if (k == 0) return false;                  /* points of no dimension */
    for (size_t i = 0; i < n; i++) {
        if (!fc_is_point_head(elem[i])) return false;                /* mixed */
        /* One head for the whole list: a red RGBColor and a 3-vector are not
         * points in a common space, whatever their arity suggests. */
        if (!expr_eq(elem[i]->data.function.head, head0)) return false;
        if (elem[i]->data.function.arg_count != k) return false;      /* ragged */
        for (size_t c = 0; c < k; c++) {
            Expr* comp = elem[i]->data.function.args[c];
            if (is_listq(comp)) return false;                         /* depth > 2 */
            if (!list_real_number_q(comp)) return false;
        }
    }
    *dim = k;
    *kind = FC_KIND_POINT;
    return true;
}

/* ------------------------------------------------------------------------- */
/* Exact stable merge sort of an index array                                 */
/* ------------------------------------------------------------------------- */

/* qsort cannot carry a failure flag out of its comparator, and list_numeric_cmp
 * can decline. A bottom-up merge sort keeps O(n log n), is stable (so equal
 * values retain input order, which the cluster-order contract depends on), and
 * threads the flag. */
typedef struct {
    Expr** key;      /* key[i] is the value ranked for index i */
    bool   ok;
} FcSortCtx;

static void fc_merge_sort(size_t* idx, size_t n, FcSortCtx* c) {
    if (n < 2) return;
    size_t* tmp = malloc(sizeof(size_t) * n);
    if (!tmp) { c->ok = false; return; }

    for (size_t width = 1; width < n; width *= 2) {
        for (size_t lo = 0; lo < n; lo += 2 * width) {
            size_t mid = lo + width;      if (mid > n) mid = n;
            size_t hi  = lo + 2 * width;  if (hi  > n) hi  = n;
            size_t i = lo, j = mid, k = lo;
            while (i < mid && j < hi) {
                int r = list_numeric_cmp(c->key[idx[i]], c->key[idx[j]], &c->ok);
                if (!c->ok) { free(tmp); return; }
                /* <= keeps the left run first on a tie: stability, which is how
                 * equal elements keep input order. */
                tmp[k++] = (r <= 0) ? idx[i++] : idx[j++];
            }
            while (i < mid) tmp[k++] = idx[i++];
            while (j < hi)  tmp[k++] = idx[j++];
        }
        memcpy(idx, tmp, sizeof(size_t) * n);
    }
    free(tmp);
}

/* ------------------------------------------------------------------------- */
/* Exact minimum spanning tree (dim > 1)                                      */
/* ------------------------------------------------------------------------- */

/* Prim's algorithm over exact SQUARED distances, O(n^2) time and O(n) memory.
 *
 * Squared rather than true distance for one reason: it is rational for rational
 * input, so the whole tree -- and therefore the whole partition -- is decided by
 * exact comparisons. Squaring is monotone on non-negatives, so the MST of the
 * squared metric is an MST of the true metric; only the recorded weights differ,
 * which `thresh_factor` accounts for.
 *
 * O(n^2) is the honest cost of general-dimension neighbourhood work; the 1D path
 * gets an O(n log n) sort instead and never comes here. Callers cap n.
 *
 * Edges are emitted in Prim insertion order with eu = the endpoint already in the
 * tree, which gives the parent-before-child property the equal-elements fold
 * needs. Ties between equal-weight edges therefore break by insertion order --
 * deterministic for a given input, and with no prior behaviour to preserve since
 * higher dimensions were previously rejected outright. */
/* Is every coordinate already a machine number?
 *
 * If so, the input carries no precision beyond a double, and computing its
 * distances through exact Expr arithmetic preserves nothing -- it just allocates
 * n^2 expressions to reach the same answer a double would. An exact Rational, a
 * bigint or an MPFR value is different: there the exact path is the only one that
 * can order the points correctly, and it is kept.
 *
 * Integers are required to be within 2^53 so that squaring differences stays
 * exact in a double; beyond that, doubles start losing integers and the fast path
 * would silently disagree with the exact one. */
static bool fc_all_machine(const FcData* d) {
    for (size_t i = 0; i < d->n; i++) {
        for (size_t c = 0; c < d->dim; c++) {
            Expr* e = fc_comp(d, i, c);
            if (e->type == EXPR_REAL) continue;
            if (e->type == EXPR_INTEGER) {
                int64_t v = e->data.integer;
                if (v > -9007199254740992LL && v < 9007199254740992LL) continue;
            }
            return false;
        }
    }
    return true;
}

/* Distance over the machine projection, in whichever metric the edge weights are
 * being built with.
 *
 * Must agree with fc_pair_distance on the SAME input, or the two spanning-tree
 * builders would disagree and the partition would depend on whether the input
 * happened to be machine-precision. */
static double fc_sqdist(const FcData* d, size_t i, size_t j) {
    const double* a = d->coord + i * d->dim;
    const double* b = d->coord + j * d->dim;
    double s = 0.0;
    if (d->dist == FC_DIST_MANHATTAN) {
        for (size_t c = 0; c < d->dim; c++) s += fabs(a[c] - b[c]);
        return s;
    }
    for (size_t c = 0; c < d->dim; c++) {
        double t = a[c] - b[c];
        s += t * t;
    }
    return s;
}

/* Prim over doubles, for input that is machine-precision to begin with.
 *
 * Same algorithm and the same O(n^2) comparison count as the exact builder; the
 * difference is that a distance costs a few flops instead of allocating and
 * evaluating a chain of expressions. Only the n-1 chosen edge weights become
 * Exprs, so allocation drops from O(n^2) to O(n).
 *
 * Everything downstream is unaffected: bnd[] and n_distinct are derived by
 * comparing the ELEMENTS with the exact comparator, never these weights, so
 * distinctness stays exact on this path too. */
static bool fc_build_mst_machine(FcData* d) {
    size_t n = d->n;
    double* best = malloc(sizeof(double) * n);
    size_t* from = malloc(sizeof(size_t) * n);
    bool*   in   = calloc(n, sizeof(bool));
    if (!best || !from || !in) { free(best); free(from); free(in); return false; }

    in[0] = true;
    for (size_t i = 1; i < n; i++) { best[i] = fc_sqdist(d, 0, i); from[i] = 0; }

    bool ok = true;
    for (size_t added = 1; added < n; added++) {
        size_t pick = SIZE_MAX;
        double bw = 0.0;
        for (size_t i = 0; i < n; i++) {
            if (in[i]) continue;
            if (pick == SIZE_MAX || best[i] < bw) { pick = i; bw = best[i]; }
        }
        if (pick == SIZE_MAX || isnan(bw)) { ok = false; break; }

        in[pick] = true;
        size_t j = added - 1;
        d->eu[j] = from[pick];
        d->ev[j] = pick;
        d->gap[j] = expr_new_real(bw);
        if (!d->gap[j]) { ok = false; break; }

        for (size_t i = 0; i < n; i++) {
            if (in[i]) continue;
            double cand = fc_sqdist(d, pick, i);
            if (cand < best[i]) { best[i] = cand; from[i] = pick; }
        }
    }

    free(best);
    free(from);
    free(in);
    return ok;
}

static bool fc_build_mst(FcData* d) {
    size_t n = d->n;
    Expr** best = calloc(n, sizeof(Expr*));    /* exact weight to the tree */
    size_t* from = malloc(sizeof(size_t) * n); /* which tree vertex realises it */
    bool*   in   = calloc(n, sizeof(bool));
    bool ok = best && from && in;

    if (ok) {
        in[0] = true;
        for (size_t i = 1; i < n && ok; i++) {
            best[i] = fc_pair_distance(d, 0, i);
            from[i] = 0;
            if (!best[i] || !list_real_number_q(best[i])) ok = false;
        }
    }

    for (size_t added = 1; added < n && ok; added++) {
        /* Cheapest vertex not yet in the tree, by exact comparison. */
        size_t pick = SIZE_MAX;
        for (size_t i = 0; i < n && ok; i++) {
            if (in[i]) continue;
            if (pick == SIZE_MAX) { pick = i; continue; }
            if (list_numeric_cmp(best[i], best[pick], &ok) < 0) pick = i;
        }
        if (!ok || pick == SIZE_MAX) { ok = false; break; }

        in[pick] = true;
        size_t j = added - 1;                  /* this edge's index */
        d->eu[j] = from[pick];
        d->ev[j] = pick;
        d->gap[j] = best[pick];
        best[pick] = NULL;                     /* ownership moves to d->gap */

        /* Relax: does joining `pick` bring anyone closer to the tree? */
        for (size_t i = 0; i < n && ok; i++) {
            if (in[i]) continue;
            Expr* cand = fc_pair_distance(d, pick, i);
            if (!cand || !list_real_number_q(cand)) { expr_free(cand); ok = false; break; }
            if (list_numeric_cmp(cand, best[i], &ok) < 0 && ok) {
                expr_free(best[i]);
                best[i] = cand;
                from[i] = pick;
            } else {
                expr_free(cand);
            }
        }
    }

    for (size_t i = 0; i < n; i++) expr_free(best[i]);
    free(best);
    free(from);
    free(in);
    return ok;
}

/* ------------------------------------------------------------------------- */
/* Result construction                                                        */
/* ------------------------------------------------------------------------- */

/* Build {{...}, {...}} from a per-element cluster assignment.
 *
 * assign[i] is the cluster id of input element i, in 0..k-1. Clusters come out
 * ordered by FIRST OCCURRENCE in the input and elements within a cluster in
 * input order -- both fall out of two ascending passes over the input, with no
 * comparator and no sort.
 *
 * Every method produces an assignment array and calls this, so result
 * construction and its two OOM unwinds are written exactly once. Sizes are
 * exact after the counting pass, so the vectors are pre-sized and never grow
 * (the Split idiom, split.c:15, not the SubsetBuf one). */
static Expr* fc_emit_clusters(Expr** elem, size_t n, const size_t* assign, size_t k) {
    if (k == 0) return NULL;

    size_t* count = calloc(k, sizeof(size_t));
    size_t* slot  = malloc(sizeof(size_t) * k);   /* raw id -> emission slot */
    Expr**  out   = malloc(sizeof(Expr*) * k);
    if (!count || !slot || !out) { free(count); free(slot); free(out); return NULL; }

    for (size_t j = 0; j < k; j++) slot[j] = SIZE_MAX;

    size_t nslots = 0;
    for (size_t i = 0; i < n; i++) {
        size_t a = assign[i];
        count[a]++;
        if (slot[a] == SIZE_MAX) slot[a] = nslots++;   /* first occurrence wins */
    }

    /* Per-slot fill cursors, and the element buffer for each cluster. */
    Expr*** buf  = calloc(nslots, sizeof(Expr**));
    size_t* fill = calloc(nslots, sizeof(size_t));
    size_t* cap  = calloc(nslots, sizeof(size_t));
    if (!buf || !fill || !cap) {
        free(count); free(slot); free(out); free(buf); free(fill); free(cap);
        return NULL;
    }
    for (size_t a = 0; a < k; a++)
        if (slot[a] != SIZE_MAX) cap[slot[a]] = count[a];

    bool ok = true;
    for (size_t s = 0; s < nslots && ok; s++) {
        buf[s] = malloc(sizeof(Expr*) * (cap[s] ? cap[s] : 1));
        if (!buf[s]) ok = false;
    }

    if (ok) {
        for (size_t i = 0; i < n; i++) {
            size_t s = slot[assign[i]];
            buf[s][fill[s]++] = expr_copy(elem[i]);
        }
    }

    /* Wrap each cluster. The head is a named local so it can be freed if
     * expr_new_function fails -- it adopts nothing on failure. */
    size_t built = 0;
    for (size_t s = 0; s < nslots && ok; s++) {
        Expr* head = expr_new_symbol(SYM_List);
        Expr* cl   = expr_new_function(head, buf[s], fill[s]);
        if (!cl) {
            for (size_t i = 0; i < fill[s]; i++) expr_free(buf[s][i]);
            expr_free(head);
            ok = false;
            break;
        }
        out[built++] = cl;
    }

    if (!ok) {
        /* Free the clusters already wrapped, then any element buffers whose
         * contents were filled but never wrapped.
         *
         * The second loop starts at built + 1, not built: when the wrap fails
         * at slot s we have built == s and the `if (!cl)` arm above has ALREADY
         * freed that slot's elements. Starting at built would free them a
         * second time -- and expr_copy is a refcount bump (expr.c:623), not a
         * deep copy, so buf[s][i] aliases an element of the caller's input.
         * The double release would recycle a node `res` still owns and the
         * evaluator would then walk freed memory. */
        for (size_t s = 0; s < built; s++) expr_free(out[s]);
        for (size_t s = built + 1; s < nslots; s++) {
            if (!buf[s]) continue;
            for (size_t i = 0; i < fill[s]; i++) expr_free(buf[s][i]);
        }
        for (size_t s = 0; s < nslots; s++) free(buf[s]);
        free(buf); free(fill); free(cap); free(count); free(slot); free(out);
        return NULL;
    }

    for (size_t s = 0; s < nslots; s++) free(buf[s]);
    free(buf); free(fill); free(cap); free(count); free(slot);

    Expr* head   = expr_new_symbol(SYM_List);
    Expr* result = expr_new_function(head, out, nslots);
    if (!result) {
        for (size_t s = 0; s < nslots; s++) expr_free(out[s]);
        expr_free(head);
    }
    free(out);
    return result;
}

/* ------------------------------------------------------------------------- */
/* Count selection                                                            */
/* ------------------------------------------------------------------------- */

/* The Automatic cluster count for the gap family: one more than the number of
 * sorted-adjacent gaps exceeding FC_GAP_FACTOR times the median gap. A run of
 * equal or near-equal gaps yields no cut and therefore one cluster, which is
 * what keeps {1,...,8} and {7,7,7,7} whole.
 *
 * Computed EXACTLY, on the gap expressions. Two reasons, both found by review
 * rather than by reasoning:
 *
 *   - The median is the LOWER of the two middle gaps for an even-length list,
 *     not their average. Averaging makes a split arithmetically impossible for
 *     any 2- or 3-element input: with one or two gaps the largest gap is itself
 *     one of the values averaged, so median >= max/2 and the threshold
 *     3 * median > max can never be met. FindClusters[{0, 1, 10^12}] came back
 *     as a single cluster however far apart the values were.
 *
 *   - Comparing through the double projection loses whole scales. A gap above
 *     DBL_MAX saturates to +inf, so the threshold becomes unsatisfiable and
 *     FindClusters[{0, 1, 10^400, 2*10^400}] collapsed to one cluster, while
 *     the same shape at 10^25 worked -- a scale cliff, not a design choice.
 *
 * Sets *ok false only if it cannot decide a comparison, which the caller turns
 * into an unevaluated result rather than a fabricated count. */
static size_t fc_automatic_gap_count(const FcData* d, bool* ok) {
    *ok = true;
    if (d->n_gap == 0) return 1;

    size_t* gi = malloc(sizeof(size_t) * d->n_gap);
    if (!gi) { *ok = false; return 1; }
    for (size_t j = 0; j < d->n_gap; j++) gi[j] = j;

    FcSortCtx c = { d->gap, true };
    fc_merge_sort(gi, d->n_gap, &c);
    if (!c.ok) { free(gi); *ok = false; return 1; }

    /* Lower median: index (m - 1) / 2 is the middle for odd m and the lower of
     * the two middles for even m. */
    Expr* median = d->gap[gi[(d->n_gap - 1) / 2]];
    free(gi);

    Expr* mul[2] = { expr_new_integer((long)d->thresh_factor), expr_copy(median) };
    Expr* thresh = eval_and_free(internal_times(mul, 2));
    if (!list_real_number_q(thresh)) { expr_free(thresh); *ok = false; return 1; }

    size_t cuts = 0;
    for (size_t j = 0; j < d->n_gap; j++) {
        bool cok = true;
        if (list_numeric_cmp(d->gap[j], thresh, &cok) > 0) cuts++;
        if (!cok) { expr_free(thresh); *ok = false; return 1; }
    }
    expr_free(thresh);
    return cuts + 1;
}

/* Reduce a method's natural count to what the count mode requires.
 *   Automatic : the natural count, unchanged
 *   Bounded   : min(natural, n)   -- MAY BE FEWER than n, that is the point
 *   Fixed     : exactly n
 * Always capped at the distinct-value count: no method can separate two equal
 * elements, so asking for more clusters than distinct values yields fewer. */
static size_t fc_target_count(size_t natural, FcCount spec, size_t n_distinct) {
    size_t want;
    switch (spec.mode) {
        case FC_COUNT_BOUNDED: want = (natural < spec.n) ? natural : spec.n; break;
        case FC_COUNT_FIXED:   want = spec.n;                                break;
        default:               want = natural;                               break;
    }
    if (want > n_distinct) want = n_distinct;
    if (want < 1) want = 1;
    return want;
}

/* ------------------------------------------------------------------------- */
/* Method: Agglomerate == SpanningTree                                        */
/* ------------------------------------------------------------------------- */

/* Single-linkage clustering equals cutting the largest edges of the minimum
 * spanning tree, and in 1D the MST of points on a line IS the sorted adjacency
 * chain. So both method names reduce to the same operation -- sort, cut the
 * widest gaps -- and share this one implementation rather than drifting into
 * two. Confirmed empirically: Mathematica's Agglomerate and SpanningTree agree
 * on every probe, with and without a count.
 *
 * Gap selection is EXACT (list_numeric_cmp over the gap expressions), so the
 * fixed-count path orders 1/10^25-scale rationals correctly. Only the Automatic
 * threshold uses the double projection, because it carries a fitted constant. */
static bool fc_method_gap(const FcData* d, FcCount spec, const FcOpts* o,
                          size_t* assign, size_t* k) {
    (void)o;
    size_t n = d->n;

    bool nok = true;
    size_t natural = fc_automatic_gap_count(d, &nok);
    if (!nok) return false;
    size_t target  = fc_target_count(natural, spec, d->n_distinct);

    /* Rank gap positions by width, widest first, earlier position winning a
     * tie. Sorting ascending by an exact comparator and then taking from the
     * end keeps the comparison exact and the tie-break stable. */
    size_t ncut = (target > 0) ? target - 1 : 0;
    if (ncut > d->n_gap) ncut = d->n_gap;

    bool* cut = calloc(d->n_gap ? d->n_gap : 1, sizeof(bool));
    size_t* gidx = malloc(sizeof(size_t) * (d->n_gap ? d->n_gap : 1));
    if (!cut || !gidx) { free(cut); free(gidx); return false; }
    for (size_t j = 0; j < d->n_gap; j++) gidx[j] = j;

    if (ncut > 0) {
        FcSortCtx c = { d->gap, true };
        fc_merge_sort(gidx, d->n_gap, &c);
        if (!c.ok) { free(cut); free(gidx); return false; }

        /* gidx is now ascending by width and, being a stable sort, ascending by
         * position within a run of equal widths. We want the widest first and,
         * among equal widths, the EARLIEST position -- Mathematica breaks the
         * tie that way (on {1,2,3,5,8,9,10} with n=4 it cuts the first of the
         * four unit gaps, not the last). So walk the runs from the last
         * backwards, but each run left to right. */
        size_t taken = 0;
        size_t hi = d->n_gap;
        while (hi > 0 && taken < ncut) {
            /* [lo, hi) is the trailing run of equal width. */
            size_t lo = hi - 1;
            while (lo > 0) {
                bool eq_ok = true;
                if (list_numeric_cmp(d->gap[gidx[lo - 1]], d->gap[gidx[hi - 1]], &eq_ok) != 0)
                    break;
                if (!eq_ok) { free(cut); free(gidx); return false; }
                lo--;
            }
            for (size_t t = lo; t < hi && taken < ncut; t++, taken++)
                cut[gidx[t]] = true;
            hi = lo;
        }
    }
    free(gidx);

    /* Connected components of the tree minus the cut edges.
     *
     * This replaces a walk along the sorted chain, which only existed because in
     * 1D the tree IS a chain and its components are contiguous runs. Union-find
     * over the surviving edges says the same thing there and is the only version
     * that also works on a general tree. Component ids come out in an arbitrary
     * order, which costs nothing: fc_emit_clusters orders clusters by first
     * occurrence in input order regardless. */
    size_t* parent = malloc(sizeof(size_t) * n);
    size_t* rank   = calloc(n ? n : 1, sizeof(size_t));
    if (!parent || !rank) { free(cut); free(parent); free(rank); return false; }
    for (size_t i = 0; i < n; i++) parent[i] = i;

    /* Path halving plus union by rank, which is not a micro-optimisation here.
     * A naive find that walks to the root is O(n) per query when the surviving
     * edges form one long chain -- exactly the shape of a small cluster count on
     * sorted 1-D data -- and made the whole method quadratic: an explicit count
     * of 10 over 100,000 points took 2.24 s against 0.087 s for Automatic on the
     * same input, and timings grew 4x per doubling. Automatic looked fine only
     * because many cuts leave many short components. */
    for (size_t j = 0; j < d->n_gap; j++) {
        if (cut[j]) continue;
        size_t a = d->eu[j], b = d->ev[j];
        while (parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; }
        while (parent[b] != b) { parent[b] = parent[parent[b]]; b = parent[b]; }
        if (a == b) continue;
        if (rank[a] < rank[b]) { size_t t = a; a = b; b = t; }
        parent[b] = a;
        if (rank[a] == rank[b]) rank[a]++;
    }
    free(cut);
    free(rank);

    /* Number the roots by first appearance in input order. */
    size_t* label = malloc(sizeof(size_t) * n);
    if (!label) { free(parent); return false; }
    for (size_t i = 0; i < n; i++) label[i] = SIZE_MAX;

    size_t id = 0;
    for (size_t i = 0; i < n; i++) {
        size_t r = i;
        while (parent[r] != r) { parent[r] = parent[parent[r]]; r = parent[r]; }
        if (label[r] == SIZE_MAX) label[r] = id++;
        assign[i] = label[r];
    }
    free(parent);
    free(label);

    *k = id;
    return true;
}

/* ------------------------------------------------------------------------- */
/* Sorted-neighbourhood kernel                                                */
/* ------------------------------------------------------------------------- */

/* The primitives the density family shares. Every one is O(1), O(k) or O(n)
 * ONLY because the array is sorted -- the same neighbourhood queries in general
 * dimension are the expensive part of DBSCAN, mean shift and Jarvis-Patrick.
 * `sv` is always the values in sorted order, i.e. sv[j] = val[order[j]]. */

/* The eps-window that used to live here went when DBSCAN was unified onto the
 * general kernel: it was that method's only caller, and the sorted window is
 * precisely the part with no meaning off a line. */

/* Half-open [*lo, *hi) covering sv[i] and its k nearest neighbours. In 1D they
 * are contiguous around i, so expanding whichever side is closer is exact. */
static void fc_knn_window(const double* sv, size_t n, size_t i, size_t k,
                          size_t* lo, size_t* hi) {
    size_t a = i, b = i + 1;          /* [a, b) currently holds just i */
    size_t got = 0;
    while (got < k && (a > 0 || b < n)) {
        if (a == 0)              { b++; }
        else if (b == n)         { a--; }
        else {
            double dl = sv[i] - sv[a - 1];
            double dr = sv[b] - sv[i];
            if (dl <= dr) a--; else b++;
        }
        got++;
    }
    *lo = a; *hi = b;
}

/* There is deliberately no kernel-density helper here. The plan called for one,
 * but mean shift's update is the kernel-weighted MEAN of the sample, which
 * never needs the density value itself -- fc_shift_cluster computes the
 * weighted mean directly in one pass. Adding a fc_kde() that nothing calls
 * would be dead code that -Werror=unused-function rightly rejects. */

/* Median of the sorted-adjacent gaps: the default length scale for eps and for
 * kernel bandwidths, so every density method shares one notion of "close". */
static double fc_median_gap(const double* sv, size_t n) {
    if (n < 2) return 0.0;
    size_t m = n - 1;
    double* g = malloc(sizeof(double) * m);
    if (!g) return 0.0;
    for (size_t j = 0; j < m; j++) g[j] = sv[j + 1] - sv[j];
    for (size_t i = 1; i < m; i++) {
        double v = g[i]; size_t j = i;
        while (j > 0 && g[j - 1] > v) { g[j] = g[j - 1]; j--; }
        g[j] = v;
    }
    double med = (m % 2) ? g[m / 2] : 0.5 * (g[m / 2 - 1] + g[m / 2]);
    free(g);
    return med;
}

/* Sorted values, allocated by the caller of a density method. */
static double* fc_sorted_values(const FcData* d) {
    double* sv = malloc(sizeof(double) * d->n);
    if (!sv) return NULL;
    for (size_t j = 0; j < d->n; j++) sv[j] = d->val[d->order[j]];
    return sv;
}

/* Turn "cluster id per SORTED position" into the input-indexed assignment every
 * method must return, renumbering ids to be contiguous from 0. */
static void fc_scatter(const FcData* d, const size_t* sorted_id, size_t* assign,
                       size_t* k) {
    size_t* remap = malloc(sizeof(size_t) * d->n);
    if (!remap) { *k = 0; return; }
    for (size_t j = 0; j < d->n; j++) remap[j] = SIZE_MAX;

    size_t next = 0;
    for (size_t j = 0; j < d->n; j++) {
        size_t id = sorted_id[j];
        if (remap[id] == SIZE_MAX) remap[id] = next++;
        assign[d->order[j]] = remap[id];
    }
    free(remap);
    *k = next;
}

/* ------------------------------------------------------------------------- */
/* Dimension-general helpers                                                  */
/*                                                                            */
/* Everything here works for FC_KIND_SCALAR and FC_KIND_POINT alike, so the    */
/* methods built on them need one implementation rather than two. The scalar   */
/* case is not a special case: d->val is a row-major n x 1 buffer indexed by   */
/* INPUT position, which is exactly the layout d->coord has for dim == 1.      */
/* ------------------------------------------------------------------------- */

/* The coordinate buffer, whichever field holds it. Both are n x dim row-major
 * and input-indexed; only their names differ. NULL for FC_KIND_SEQUENCE, which
 * has no coordinates at all. */
static const double* fc_points(const FcData* d) {
    if (d->kind == FC_KIND_POINT)  return d->coord;
    if (d->kind == FC_KIND_SCALAR) return d->val;
    return NULL;
}

/* LINEAR distance from an arbitrary position to data point t.
 *
 * Linear, never squared, unlike fc_sqdist -- a kernel bandwidth and a merge
 * tolerance are lengths, and mixing a squared distance into either silently
 * changes the scale by an exponent rather than a factor. */
static double fc_dist_to_point(const FcData* d, const double* pts,
                               const double* p, size_t t) {
    const double* b = pts + t * d->dim;
    double s = 0.0;
    if (d->dist == FC_DIST_MANHATTAN) {
        for (size_t c = 0; c < d->dim; c++) s += fabs(p[c] - b[c]);
        return s;
    }
    for (size_t c = 0; c < d->dim; c++) { double t2 = p[c] - b[c]; s += t2 * t2; }
    return sqrt(s);
}

/* LINEAR distance between two arbitrary positions. */
static double fc_dist_pos(const FcData* d, const double* a, const double* b) {
    double s = 0.0;
    if (d->dist == FC_DIST_MANHATTAN) {
        for (size_t c = 0; c < d->dim; c++) s += fabs(a[c] - b[c]);
        return s;
    }
    for (size_t c = 0; c < d->dim; c++) { double t = a[c] - b[c]; s += t * t; }
    return sqrt(s);
}

/* The length scale of the data, in any dimension: the median spanning-tree edge
 * weight, as a LINEAR distance.
 *
 * This is not merely analogous to fc_median_gap, it GENERALISES it exactly. On a
 * line the minimum spanning tree IS the sorted adjacency chain, so the tree's
 * edge weights are precisely the adjacent gaps that fc_median_gap takes the
 * median of. Two consequences, both relied on below:
 *
 *   - the median edge weight reproduces fc_median_gap on scalar input, so one
 *     scale serves both dimensionalities and the 1-D answers do not move;
 *   - the MEAN edge weight reproduces (max - min) / (n - 1) there, because the
 *     adjacent gaps of a sorted line sum to its range. That is exactly the
 *     fallback the 1-D code used when the median gap came out zero, so the
 *     fallback generalises for free too.
 *
 * Weights are converted to lengths BEFORE the median is taken. Taking the median
 * of squared weights and rooting afterwards happens to give the same answer,
 * since a root is monotone, but the mean does not commute that way and the
 * fallback would be wrong. */
static double fc_scale_ndim(const FcData* d) {
    size_t m = d->n_gap;
    if (m == 0) return 0.0;
    double* w = malloc(sizeof(double) * m);
    if (!w) return 0.0;

    bool squared = fc_dist_is_squared(d->dist) && d->kind == FC_KIND_POINT;
    double sum = 0.0;
    for (size_t j = 0; j < m; j++) {
        double v = fc_to_double(d->gap[j]);
        if (v < 0.0) v = 0.0;                 /* a weight is a distance */
        w[j] = squared ? sqrt(v) : v;
        sum += w[j];
    }

    for (size_t i = 1; i < m; i++) {          /* insertion sort, as fc_median_gap */
        double v = w[i]; size_t j = i;
        while (j > 0 && w[j - 1] > v) { w[j] = w[j - 1]; j--; }
        w[j] = v;
    }
    double med = (m % 2) ? w[m / 2] : 0.5 * (w[m / 2 - 1] + w[m / 2]);
    free(w);

    /* A zero median does NOT mean the points coincide -- it means at least half
     * the tree edges have zero length, which any tie-heavy input has. Collapsing
     * to one cluster there would discard real structure, so fall back to the mean
     * edge length, which on a line is the range over n-1. */
    if (med <= 0.0) med = sum / (double)m;
    return med;
}

/* Union-find over converged positions, with path halving. Replaces the 1-D
 * adjacent-difference merge, which cannot generalise: it relies on a sorted order
 * that vectors do not have.
 *
 * It also reproduces that merge exactly on a line. Merging every pair within tol
 * is transitive chaining, and in one dimension a point lying between two others
 * is closer to each than they are to one another -- so pairs within tol of each
 * other are always reachable through adjacent steps, and the two agree. */
static size_t fc_uf_find(size_t* parent, size_t x) {
    while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
}

static void fc_uf_union(size_t* parent, size_t a, size_t b) {
    a = fc_uf_find(parent, a);
    b = fc_uf_find(parent, b);
    if (a != b) parent[b] = a;
}

/* Turn union-find roots into the input-indexed assignment every method returns,
 * numbering clusters by first appearance in INPUT order -- the same order
 * fc_emit_clusters and the equal-elements fold use, so the three agree. */
static void fc_assign_from_uf(size_t* parent, size_t n, size_t* assign, size_t* k) {
    size_t* label = malloc(sizeof(size_t) * n);
    if (!label) { *k = 0; return; }
    for (size_t i = 0; i < n; i++) label[i] = SIZE_MAX;
    size_t next = 0;
    for (size_t i = 0; i < n; i++) {
        size_t r = fc_uf_find(parent, i);
        if (label[r] == SIZE_MAX) label[r] = next++;
        assign[i] = label[r];
    }
    free(label);
    *k = next;
}

/* fc_merge_modes lived here: an adjacent-difference pass over sorted, monotone
 * converged positions. It was correct and it was one-dimensional to its core --
 * it needed a total order on positions, which vectors do not have. The union-find
 * merge above replaces it and agrees with it on a line, so it was deleted rather
 * than kept alongside; -Werror=unused-function makes leaving it a build failure,
 * which is the right pressure. */

/* ------------------------------------------------------------------------- */
/* Method: KMeans and KMedoids                                                */
/* ------------------------------------------------------------------------- */

/* 1D Lloyd, shared by both: `use_medoid` swaps the centre from the mean of a
 * run to its median element, which is the cost-minimising medoid of a
 * contiguous run in 1D and so is an O(1) index lookup after the sort.
 *
 * Note on KMedoids: this is Voronoi iteration (alternating k-medoids), not PAM.
 * Both are "partitioning around medoids" in the loose sense, but PAM's
 * swap-based search escapes local optima that the alternating form does not, so
 * this can settle on a worse partition than PAM would -- measured on
 * {1,2,10,12,3,1,13,25} with n=3, where it keeps 25 with {12,13} rather than
 * isolating it. Deterministic and a legitimate named algorithm; simply not the
 * stronger one. Upgrading to PAM is a follow-up, not a silent expectation.
 *
 * Initialisation is deterministic (quantile), not random, so results are
 * reproducible run to run and RandomSeeding is not needed. Because points and
 * centres are both sorted, a cluster is always a contiguous run of the sorted
 * array, and the assignment step is a merge walk: O(n + k), not O(nk). */
static bool fc_lloyd(const FcData* d, FcCount spec, size_t* assign, size_t* k,
                     bool use_medoid) {
    size_t n = d->n;
    double* sv = fc_sorted_values(d);
    if (!sv) return false;

    /* R5: `natural` must be a DATA-DRIVEN count, not n_distinct. Passing the
     * maximum possible value made min(natural, spec.n) always spec.n, so
     * UpTo[n] was indistinguishable from a bare n -- and since FC_ALLOWED
     * denies these two methods FC_COUNT_AUTOMATIC, UpTo[n] is the only
     * data-driven form a caller can even write for them. Use the same gap rule
     * the Automatic methods use. */
    size_t natural = d->n_distinct;
    if (spec.mode == FC_COUNT_BOUNDED) {
        bool nok = true;
        natural = fc_automatic_gap_count(d, &nok);
        if (!nok) { free(sv); return false; }
    }
    size_t target = fc_target_count(natural, spec, d->n_distinct);
    if (target < 1) target = 1;

    double* c   = malloc(sizeof(double) * target);
    size_t* cut = malloc(sizeof(size_t) * (target + 1));
    size_t* id  = malloc(sizeof(size_t) * n);
    if (!c || !cut || !id) { free(sv); free(c); free(cut); free(id); return false; }

    /* Seed from the DISTINCT values, not from raw sorted positions. Quantile
     * init over raw positions puts two centroids on the same value whenever the
     * data is tie-heavy, and a centroid that never wins a point leaves an empty
     * cluster -- so fewer clusters come back than were asked for. Measured
     * before the fix: FindClusters[{1,1,1,1,2,3}, 3, Method -> "KMeans"]
     * returned two clusters despite three distinct values. */
    double* uniq = malloc(sizeof(double) * n);
    if (!uniq) { free(sv); free(c); free(cut); free(id); return false; }
    size_t nu = 0;
    for (size_t j = 0; j < n; j++)
        if (j == 0 || d->bnd[j - 1]) uniq[nu++] = sv[j];   /* exact distinctness */
    for (size_t j = 0; j < target; j++)
        c[j] = uniq[(size_t)(((double)j + 0.5) * (double)nu / (double)target)];
    free(uniq);

    for (int it = 0; it < FC_MAX_ITER; it++) {
        /* Assignment: cluster j owns the sorted points nearer to c[j] than to
         * c[j+1], so the boundaries are the midpoints, found by one walk. */
        cut[0] = 0; cut[target] = n;
        size_t p = 0;
        for (size_t j = 0; j + 1 < target; j++) {
            double mid = 0.5 * (c[j] + c[j + 1]);
            while (p < n && sv[p] < mid) p++;
            cut[j + 1] = p;
        }
        for (size_t j = 1; j < target; j++)
            if (cut[j] < cut[j - 1]) cut[j] = cut[j - 1];

        bool moved = false;
        for (size_t j = 0; j < target; j++) {
            size_t lo = cut[j], hi = cut[j + 1];
            if (lo >= hi) continue;                 /* empty: leave the centre */
            double nc;
            if (use_medoid) {
                nc = sv[lo + (hi - lo) / 2];
            } else {
                double s = 0.0;
                for (size_t t = lo; t < hi; t++) s += sv[t];
                nc = s / (double)(hi - lo);
            }
            if (nc != c[j]) { c[j] = nc; moved = true; }
        }
        if (!moved) break;
    }

    /* Guarantee exactly `target` non-empty clusters.
     *
     * Better seeding makes an empty cluster rare, but Lloyd can still strand one
     * mid-iteration, and "exactly n" is what the docstring promises for a fixed
     * count. A boundary is only meaningful at a value boundary -- splitting a
     * run of equal values would put identical elements in different clusters --
     * so keep the boundaries Lloyd chose that sit on one, then top up from the
     * widest unused value boundaries, widest first with the earlier position
     * winning a tie, exactly as fc_method_gap does. n_distinct >= target holds
     * because fc_target_count caps there, so enough boundaries always exist. */
    size_t* bnd = malloc(sizeof(size_t) * (target > 0 ? target : 1));
    if (!bnd) { free(sv); free(c); free(cut); free(id); return false; }
    size_t nb = 0;
    for (size_t j = 1; j < target && nb + 1 < target; j++) {
        size_t p = cut[j];
        if (p == 0 || p >= n) continue;
        if (!d->bnd[p - 1]) continue;                   /* inside a tie run */
        if (nb > 0 && bnd[nb - 1] == p) continue;       /* already taken */
        bnd[nb++] = p;
    }

    while (nb + 1 < target) {
        size_t best_p = 0; double best_g = -1.0;
        for (size_t p = 1; p < n; p++) {
            if (!d->bnd[p - 1]) continue;
            bool used = false;
            for (size_t t = 0; t < nb; t++) if (bnd[t] == p) { used = true; break; }
            if (used) continue;
            double g = sv[p] - sv[p - 1];
            if (g > best_g) { best_g = g; best_p = p; }  /* > keeps earliest on a tie */
        }
        if (best_g < 0.0) break;                         /* no boundary left */
        bnd[nb++] = best_p;
        for (size_t t = nb - 1; t > 0 && bnd[t - 1] > bnd[t]; t--) {
            size_t tmp = bnd[t - 1]; bnd[t - 1] = bnd[t]; bnd[t] = tmp;
        }
    }

    size_t cur = 0, next_b = 0;
    for (size_t t = 0; t < n; t++) {
        while (next_b < nb && bnd[next_b] == t) { cur++; next_b++; }
        id[t] = cur;
    }
    free(bnd);

    fc_scatter(d, id, assign, k);
    free(sv); free(c); free(cut); free(id);
    return *k > 0;
}

/* Lloyd's algorithm in any dimension.
 *
 * A SECOND implementation rather than a generalisation of fc_lloyd above, and
 * deliberately so. Iteration 3 was able to unify MeanShift because its n-D form
 * provably reproduced the 1-D answers -- the median spanning-tree edge weight IS
 * the median adjacent gap on a line, and union-find merging IS the
 * adjacent-difference pass there. Nothing like that holds here: the 1-D kernel
 * seeds at quantiles of the sorted distinct values and this one seeds
 * farthest-first, and two initialisations are two algorithms that settle in
 * different local optima. Routing one dimension through the other would change
 * answers that are pinned, so unifying is a deliberate behaviour change and not a
 * refactor -- and it is not this iteration's.
 *
 * What survives the port is the shape: assign, move the centres, repeat. What does
 * not is every mechanism the 1-D kernel uses to do it quickly -- a merge walk over
 * midpoints of sorted centres, a scalar mean over a contiguous run, and an
 * empty-cluster repair that picks value boundaries in sorted order. All three need
 * a total order that vectors do not have.
 *
 * Initialisation is farthest-first (Gonzalez): the point nearest the centroid,
 * then repeatedly the point farthest from everything chosen so far. Three reasons:
 *
 *   - It is DETERMINISTIC, which the quantile seeding it replaces also was. The
 *     docs promise reproducibility without RandomSeeding and there is no
 *     RandomVariate in the tree yet, so keeping that property costs nothing here.
 *   - Starting from the centroid-nearest point rather than from index 0 makes the
 *     result independent of input order, which a k-means has no business depending
 *     on. Ties break to the lowest index so it stays deterministic.
 *   - Every centre is a distinct data point, so each cluster starts owning at
 *     least itself, which is what makes empty clusters rare rather than routine.
 *     Distinctness is guaranteed: target is capped at n_distinct, so while fewer
 *     than target centres are chosen some point is still at nonzero distance from
 *     all of them.
 */
static bool fc_lloyd_ndim(const FcData* d, FcCount spec, size_t* assign, size_t* k,
                          bool medoid) {
    size_t n = d->n, dim = d->dim;
    const double* pts = fc_points(d);
    if (!pts) return false;

    /* Same count rule as the 1-D kernel. Both helpers read only d->gap, d->bnd and
     * d->n_distinct, all of which are kind-agnostic -- which is why the gap methods
     * already worked in n-D and why these need no porting. */
    size_t natural = d->n_distinct;
    if (spec.mode == FC_COUNT_BOUNDED) {
        bool nok = true;
        natural = fc_automatic_gap_count(d, &nok);
        if (!nok) return false;
    }
    size_t target = fc_target_count(natural, spec, d->n_distinct);
    if (target < 1) target = 1;
    if (target > n) target = n;

    /* Θ(FC_MAX_ITER * n * target * dim), so the cap is on that PRODUCT and not on
     * n, which is the mistake the shift methods' cap invites by analogy. Lloyd is
     * LINEAR in n -- cheaper than the spanning tree already built for this input --
     * so capping n would decline work the builder had just finished paying for:
     * 20000 machine points in ten dimensions is admitted by FC_NDIM_MACHINE_MAX_N
     * and there is no reason for k-means to be the step that refuses it.
     *
     * What can still run away is target approaching n, where the assignment scan
     * becomes quadratic. At the ceiling this budget is a few seconds of assignment
     * work, the same order as the tree build that preceded it. */
    if (target != 0 && dim != 0 &&
        n > FC_LLOYD_MAX_WORK / target / dim) return false;
    /* KMedoids carries a SECOND, tighter ceiling, because its update step is a
     * different complexity class: the medoid search compares every member against
     * every other member of its own cluster, which sums to O(n^2 * dim) per
     * iteration however small k is. KMeans' mean is O(n * dim). So the bilinear
     * budget above, which admits 20000 points at a small k, is the wrong bound
     * here and would admit a run of minutes. */
    if (medoid && dim != 0 && n != 0 && n > FC_LLOYD_MAX_WORK / n / dim) return false;

    double* c     = malloc(sizeof(double) * target * dim);
    size_t* owner = malloc(sizeof(size_t) * n);
    /* One buffer, two jobs across the two phases: the distance to the nearest
     * chosen centre while seeding, then the distance to its own centre once the
     * assignment exists. The second is what the empty-cluster repair reads. */
    double* nd    = malloc(sizeof(double) * n);
    size_t* count = malloc(sizeof(size_t) * target);
    double* acc   = malloc(sizeof(double) * target * dim);
    if (!c || !owner || !nd || !count || !acc) {
        free(c); free(owner); free(nd); free(count); free(acc); return false;
    }

    /* ---- Seed 0: the point nearest the centroid ---- */
    double* mean = calloc(dim, sizeof(double));
    if (!mean) { free(c); free(owner); free(nd); free(count); return false; }
    for (size_t i = 0; i < n; i++)
        for (size_t comp = 0; comp < dim; comp++) mean[comp] += pts[i * dim + comp];
    for (size_t comp = 0; comp < dim; comp++) mean[comp] /= (double)n;

    size_t seed = 0; double bestd = -1.0;
    for (size_t i = 0; i < n; i++) {
        double dd = fc_dist_to_point(d, pts, mean, i);
        if (bestd < 0.0 || dd < bestd) { bestd = dd; seed = i; }
    }
    free(mean);
    memcpy(c, pts + seed * dim, sizeof(double) * dim);

    /* ---- Seeds 1..target-1: farthest from everything chosen ---- */
    for (size_t i = 0; i < n; i++) nd[i] = fc_dist_to_point(d, pts, c, i);
    for (size_t j = 1; j < target; j++) {
        size_t pick = 0; double far = -1.0;
        for (size_t i = 0; i < n; i++)
            if (nd[i] > far) { far = nd[i]; pick = i; }   /* > keeps lowest index */
        memcpy(c + j * dim, pts + pick * dim, sizeof(double) * dim);
        for (size_t i = 0; i < n; i++) {
            double dd = fc_dist_to_point(d, pts, c + j * dim, i);
            if (dd < nd[i]) nd[i] = dd;
        }
    }

    /* ---- Lloyd ---- */
    for (int it = 0; it < FC_MAX_ITER; it++) {
        for (size_t i = 0; i < n; i++) {
            size_t best = 0; double bd = -1.0;
            for (size_t j = 0; j < target; j++) {
                double dd = fc_dist_to_point(d, pts, c + j * dim, i);
                if (bd < 0.0 || dd < bd) { bd = dd; best = j; }  /* lowest j on tie */
            }
            owner[i] = best;
            nd[i] = bd;
        }

        for (size_t j = 0; j < target; j++) count[j] = 0;
        for (size_t i = 0; i < n; i++) count[owner[i]]++;

        /* An empty cluster is reseeded at the point currently worst served by its
         * own centre. The 1-D kernel instead split at the widest unused value
         * boundary, which needs an order over positions; "farthest from its
         * centre" is the order-free statement of the same idea, and it cannot
         * loop, since the reseeded centre then owns at least that point. */
        for (size_t j = 0; j < target; j++) {
            if (count[j] > 0) continue;
            size_t pick = SIZE_MAX; double far = -1.0;
            for (size_t i = 0; i < n; i++) {
                if (count[owner[i]] < 2) continue;      /* do not empty another */
                if (nd[i] > far) { far = nd[i]; pick = i; }
            }
            if (pick == SIZE_MAX) break;                /* nothing to give */
            memcpy(c + j * dim, pts + pick * dim, sizeof(double) * dim);
            count[owner[pick]]--;
            owner[pick] = j;
            count[j] = 1;
            nd[pick] = 0.0;
        }

        bool moved = false;
        if (medoid) {
            /* KMedoids: the new centre must BE one of the members, the one whose
             * total distance to the rest of its cluster is least. The 1-D kernel
             * could take the middle element of a contiguous sorted run in O(1);
             * off a line there is no such shortcut, so this is the real search,
             * and it is what makes the method quadratic rather than linear in n.
             *
             * Restricting centres to data points is the whole point of the method:
             * it is what lets a medoid resist an outlier that would drag a mean,
             * and what lets the method run on a metric with no notion of average. */
            for (size_t j = 0; j < target; j++) {
                if (count[j] == 0) continue;
                size_t best = SIZE_MAX; double bestsum = 0.0;
                for (size_t a = 0; a < n; a++) {
                    if (owner[a] != j) continue;
                    double s = 0.0;
                    for (size_t b = 0; b < n; b++) {
                        if (owner[b] != j) continue;
                        s += fc_dist_pos(d, pts + a * dim, pts + b * dim);
                    }
                    if (best == SIZE_MAX || s < bestsum) { bestsum = s; best = a; }
                }
                if (best == SIZE_MAX) continue;
                for (size_t comp = 0; comp < dim; comp++) {
                    double v = pts[best * dim + comp];
                    if (v != c[j * dim + comp]) { c[j * dim + comp] = v; moved = true; }
                }
            }
        } else {
            /* KMeans: the component-wise mean of the members. That the mean
             * generalises component-wise is the whole reason this step is
             * mechanical rather than a new derivation. */
            for (size_t j = 0; j < target * dim; j++) acc[j] = 0.0;
            for (size_t i = 0; i < n; i++)
                for (size_t comp = 0; comp < dim; comp++)
                    acc[owner[i] * dim + comp] += pts[i * dim + comp];
            for (size_t j = 0; j < target; j++) {
                if (count[j] == 0) continue;
                for (size_t comp = 0; comp < dim; comp++) {
                    double v = acc[j * dim + comp] / (double)count[j];
                    if (v != c[j * dim + comp]) { c[j * dim + comp] = v; moved = true; }
                }
            }
        }
        if (!moved) break;
    }

    /* Number clusters by first appearance in INPUT order, matching what
     * fc_emit_clusters and the equal-elements fold do, and dropping any cluster
     * that ended up empty so the count reported is the count returned. */
    size_t* label = malloc(sizeof(size_t) * target);
    if (!label) {
        free(c); free(owner); free(nd); free(count); free(acc); return false;
    }
    for (size_t j = 0; j < target; j++) label[j] = SIZE_MAX;
    size_t next = 0;
    for (size_t i = 0; i < n; i++) {
        if (label[owner[i]] == SIZE_MAX) label[owner[i]] = next++;
        assign[i] = label[owner[i]];
    }
    *k = next;

    free(label); free(c); free(owner); free(nd); free(count); free(acc);
    return *k > 0;
}

static bool fc_method_kmeans(const FcData* d, FcCount spec, const FcOpts* o,
                             size_t* assign, size_t* k) {
    (void)o;
    if (d->kind == FC_KIND_POINT) return fc_lloyd_ndim(d, spec, assign, k, false);
    return fc_lloyd(d, spec, assign, k, false);
}

static bool fc_method_kmedoids(const FcData* d, FcCount spec, const FcOpts* o,
                               size_t* assign, size_t* k) {
    (void)o;
    if (d->kind == FC_KIND_POINT) return fc_lloyd_ndim(d, spec, assign, k, true);
    return fc_lloyd(d, spec, assign, k, true);
}

/* ------------------------------------------------------------------------- */
/* Method: DBSCAN                                                             */
/* ------------------------------------------------------------------------- */

/* Density-based, Automatic count only. A point is a core point when its
 * eps-window holds at least MinPoints members; core points whose windows
 * overlap join. A point in no dense region is NOISE, and rather than dropping
 * it -- which would lose an input element -- it becomes its own singleton
 * cluster, so the result is always a partition. */
/* DBSCAN in any dimension, and the ONLY implementation: the textbook
 * formulation, with real eps-neighbourhoods instead of a sorted window and
 * union-find instead of a left-to-right sweep.
 *
 * Unlike KMeans, which needs a second implementation because two initialisations
 * settle in different local optima, this one replaced the 1-D kernel outright.
 * The old kernel linked only ADJACENT sorted pairs -- "consecutive points, one of
 * them core" -- where DBSCAN specifies any pair within eps joined through a core
 * point. The two are not obviously the same rule, and at the default MinPoints of
 * 2 they provably collapse together, since a point with any neighbour inside eps
 * is then automatically core and both reduce to single linkage at eps. Above that
 * default the argument runs out, so the question was put to the pin suite instead
 * of settled by reasoning: all 22 one-dimensional pins pass unchanged through this
 * kernel, so the general rule is answer-preserving on a line and there is no case
 * for keeping a second copy.
 *
 * fc_eps_window went with it -- DBSCAN was its only caller, and
 * -Werror=unused-function turns leaving a dead 1-D helper behind into a build
 * failure, which is the right pressure.
 */
static bool fc_dbscan_ndim(const FcData* d, const FcOpts* o,
                           size_t* assign, size_t* k) {
    size_t n = d->n, dim = d->dim;
    const double* pts = fc_points(d);
    if (!pts) return false;

    /* Quadratic in n and linear in dim, like every neighbourhood method here: a
     * sort is what made this linear on a line, and vectors have no sort. The
     * point builder's own ceiling already bounds n for this input. */
    if (dim != 0 && n > FC_LLOYD_MAX_WORK / dim) return false;

    double eps = o->radius_given ? o->radius
                                 : (double)FC_GAP_FACTOR * fc_scale_ndim(d);
    size_t minpts = o->min_points_given ? (size_t)o->min_points : 2;

    bool*   core = calloc(n ? n : 1, sizeof(bool));
    size_t* uf   = malloc(sizeof(size_t) * (n ? n : 1));
    if (!core || !uf) { free(core); free(uf); return false; }
    for (size_t i = 0; i < n; i++) uf[i] = i;

    /* Core test counts the point itself, matching the 1-D window which spans
     * [lo, hi) around j inclusive. */
    for (size_t i = 0; i < n; i++) {
        size_t cnt = 0;
        for (size_t j = 0; j < n; j++)
            if (fc_dist_to_point(d, pts, pts + j * dim, i) <= eps) cnt++;
        core[i] = (cnt >= minpts);
    }

    /* Link every eps-close pair in which at least one side is core -- the same
     * predicate the 1-D sweep applies to adjacent pairs, lifted off the line.
     * A point in no dense region is joined to nothing and so falls out as its own
     * singleton, which is how noise stays in the partition instead of being
     * dropped and losing an input element. */
    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++)
            if ((core[i] || core[j]) &&
                fc_dist_pos(d, pts + i * dim, pts + j * dim) <= eps)
                fc_uf_union(uf, i, j);

    fc_assign_from_uf(uf, n, assign, k);
    free(core); free(uf);
    return *k > 0;
}

static bool fc_method_dbscan(const FcData* d, FcCount spec, const FcOpts* o,
                             size_t* assign, size_t* k) {
    (void)spec;
    /* One path for both dimensionalities. SEQUENCE has no coordinates and so still
     * declines; everything else goes through the general kernel. */
    return fc_points(d) ? fc_dbscan_ndim(d, o, assign, k) : false;
}

/* ------------------------------------------------------------------------- */
/* Method: MeanShift and NeighborhoodContraction                              */
/* ------------------------------------------------------------------------- */

/* Both shift every point toward higher density and then merge points that
 * arrived at the same place; only the update rule differs, so the merge is
 * shared (fc_merge_modes). MeanShift takes a kernel-weighted mean over the
 * whole sample; NeighborhoodContraction takes a flat mean over a fixed radius.
 * Both preserve order -- a monotone update cannot make two sorted points cross
 * -- so the merge can be a single adjacent-difference pass. */
static bool fc_shift_cluster(const FcData* d, const FcOpts* o, size_t* assign,
                             size_t* k, bool flat_kernel) {
    size_t n = d->n, dim = d->dim;

    /* Strings have no coordinates to shift toward anything, so there is no
     * meaningful update rule -- declined rather than approximated. */
    const double* pts = fc_points(d);
    if (!pts) return false;

    /* Quadratic in n and linear in dim. The 1-D cap was chosen for
     * Theta(FC_MAX_ITER * n^2); dim multiplies that, so the cap is divided by the
     * dimension rather than left to grow with it. */
    size_t cap = (dim > 1) ? (FC_SHIFT_MAX_N / dim) : FC_SHIFT_MAX_N;
    if (cap < 2) cap = 2;
    if (n > cap) return false;

    /* Two scales, deliberately separate -- conflating them is what produced
     * both of the bugs found here.
     *
     *   scale : the resolution of the data, one median gap. Two converged
     *           points closer than this are the same mode.
     *   h     : the kernel bandwidth. The flat kernel sees nothing outside its
     *           radius so it needs FC_GAP_FACTOR * scale (measured: a flat
     *           radius of one scale left {1,2,3} split three ways); the
     *           Gaussian has infinite support and one scale is right.
     *
     * Tying the merge tolerance to h instead of to scale broke both ends. At
     * h/2 with a one-scale Gaussian, evenly spaced data never merged --
     * FindClusters[Range[40], Method -> "MeanShift"] gave 34 clusters, because
     * every interior point of a flat density is already stationary and h/2
     * cannot span the one scale between neighbours. Widening the Gaussian to
     * FC_GAP_FACTOR * scale fixed that and immediately over-merged the other
     * direction: {1,1,1,1,100} came back whole. */
    double scale = fc_scale_ndim(d);

    double h = o->radius_given ? o->radius
                               : (flat_kernel ? (double)FC_GAP_FACTOR : 1.0) * scale;
    if (o->radius_given) scale = o->radius;

    if (h <= 0.0) {                     /* genuinely no spread: one cluster */
        for (size_t i = 0; i < n; i++) assign[i] = 0;
        *k = 1;
        return true;
    }

    /* Positions are dim-vectors now, so every buffer is n * dim and every
     * "distance" goes through the metric rather than fabs of a difference. */
    double* pos = malloc(sizeof(double) * n * dim);
    double* nxt = malloc(sizeof(double) * n * dim);
    size_t* parent = malloc(sizeof(size_t) * n);
    if (!pos || !nxt || !parent) { free(pos); free(nxt); free(parent); return false; }
    memcpy(pos, pts, sizeof(double) * n * dim);

    for (int it = 0; it < FC_MAX_ITER; it++) {
        double delta = 0.0;
        for (size_t j = 0; j < n; j++) {
            double* pj = pos + j * dim;
            double* nj = nxt + j * dim;
            /* The update is a weighted MEAN of the sample, which generalises
             * component-wise: one shared weight per point, applied to every
             * coordinate. That is the whole of what makes this port mechanical. */
            double den = 0.0;
            for (size_t c = 0; c < dim; c++) nj[c] = 0.0;
            for (size_t t = 0; t < n; t++) {
                double dist = fc_dist_to_point(d, pts, pj, t);
                double w;
                if (flat_kernel) {
                    if (dist > h) continue;              /* sees nothing outside */
                    w = 1.0;
                } else {
                    double z = dist / h;
                    w = exp(-0.5 * z * z);
                }
                const double* bt = pts + t * dim;
                for (size_t c = 0; c < dim; c++) nj[c] += w * bt[c];
                den += w;
            }
            if (den > 0.0) { for (size_t c = 0; c < dim; c++) nj[c] /= den; }
            else           { for (size_t c = 0; c < dim; c++) nj[c] = pj[c]; }

            double mv = fc_dist_pos(d, nj, pj);
            if (mv > delta) delta = mv;
        }
        memcpy(pos, nxt, sizeof(double) * n * dim);
        if (delta < h * 1e-6) break;
    }

    /* Merge every pair of converged positions within tolerance. The 1-D code
     * walked adjacent sorted positions instead, which vectors have no analogue
     * of; on a line the two agree, because a point between two others is closer
     * to each of them than they are to one another. */
    double tol = scale * FC_MERGE_SLACK;
    for (size_t i = 0; i < n; i++) parent[i] = i;
    for (size_t i = 0; i < n; i++)
        for (size_t j = i + 1; j < n; j++)
            if (fc_dist_pos(d, pos + i * dim, pos + j * dim) <= tol)
                fc_uf_union(parent, i, j);

    fc_assign_from_uf(parent, n, assign, k);
    free(pos); free(nxt); free(parent);
    return *k > 0;
}

static bool fc_method_meanshift(const FcData* d, FcCount spec, const FcOpts* o,
                                size_t* assign, size_t* k) {
    (void)spec; return fc_shift_cluster(d, o, assign, k, false);
}

static bool fc_method_neighborhood(const FcData* d, FcCount spec, const FcOpts* o,
                                   size_t* assign, size_t* k) {
    (void)spec; return fc_shift_cluster(d, o, assign, k, true);
}

/* ------------------------------------------------------------------------- */
/* Method: JarvisPatrick                                                      */
/* ------------------------------------------------------------------------- */

/* Shared-nearest-neighbour: two points join when each is in the other's k-NN
 * list and they share at least `t` neighbours. In 1D a k-NN list is the
 * contiguous window fc_knn_window returns, so both tests are index arithmetic
 * and the whole pass is O(nk). Only adjacent sorted points can be mutual
 * neighbours here, which keeps clusters contiguous. */
/* Jarvis-Patrick in any dimension: real k-NN lists and a real set intersection,
 * where the 1-D kernel had index arithmetic on contiguous windows.
 *
 * One detail carried over deliberately rather than reinvented: the 1-D window
 * spans sv[i] AND its k nearest, so a point is a member of its own neighbour
 * list and counts toward the shared total. The lists built here include self for
 * the same reason, which is also the classical formulation -- two points that are
 * each other's neighbours already share those two members.
 *
 * Lists are kept sorted by INDEX rather than by distance, so the intersection is
 * a linear merge instead of a nested scan. Selection breaks ties toward the lower
 * index, so the result does not depend on scan order.
 *
 * This one does NOT replace the 1-D kernel, and finding that out corrected the
 * procedure DBSCAN established. Routing scalars through here leaves all 22
 * one-dimensional pins passing -- which is what licensed the DBSCAN unification --
 * and yet moves an answer that `list_tests` covers and the pin file does not:
 *
 *   FindClusters[{1, 2, 10, 12, 3, 1, 13, 25},
 *                Method -> {"JarvisPatrick", "NeighborCount" -> 2}]
 *     1-D kernel: {{1, 2, 1}, {10, 12, 13}, {3}, {25}}
 *     this kernel: {{1, 2, 3, 1}, {10, 12, 13}, {25}}
 *
 * So the unify test is "does any covered 1-D answer move", not "do the pins
 * pass". The two rules differ because the 1-D kernel counts shared neighbours as
 * the OVERLAP OF TWO CONTIGUOUS WINDOWS and restricts linkage to adjacent sorted
 * pairs, where this one takes a true set intersection over all pairs -- the
 * textbook formulation. The general rule is very likely the better one, but
 * adopting it on a line moves a checked answer, so it is a deliberate behaviour
 * change to argue separately rather than a side effect of a port.
 */
static bool fc_jarvis_ndim(const FcData* d, const FcOpts* o,
                           size_t* assign, size_t* k) {
    size_t n = d->n, dim = d->dim;
    const double* pts = fc_points(d);
    if (!pts) return false;
    /* Quadratic in n, as the neighbourhood methods all are off a line; the point
     * builder's ceiling is what actually bounds n for this input. */
    if (dim != 0 && n > FC_LLOYD_MAX_WORK / dim) return false;

    size_t kn = o->neighbor_count_given ? (size_t)o->neighbor_count : 5;
    if (kn > n - 1) kn = (n > 1) ? n - 1 : 0;
    size_t thresh = (kn + 1) / 2;
    if (thresh < 1) thresh = 1;

    size_t w = kn + 1;                    /* list width, self included */
    size_t* nbr = malloc(sizeof(size_t) * n * w);
    double* nds = malloc(sizeof(double) * w);
    size_t* uf  = malloc(sizeof(size_t) * (n ? n : 1));
    if (!nbr || !nds || !uf) { free(nbr); free(nds); free(uf); return false; }
    for (size_t i = 0; i < n; i++) uf[i] = i;

    for (size_t i = 0; i < n; i++) {
        size_t* row = nbr + i * w;
        size_t cnt = 0;
        /* Insertion selection of the kn nearest, held sorted by distance. Ties
         * keep the incumbent, which is the lower index because j ascends. */
        for (size_t j = 0; j < n; j++) {
            if (j == i) continue;
            double dd = fc_dist_to_point(d, pts, pts + i * dim, j);
            if (cnt < kn) {
                size_t p = cnt++;
                while (p > 0 && nds[p - 1] > dd) {
                    nds[p] = nds[p - 1]; row[p] = row[p - 1]; p--;
                }
                nds[p] = dd; row[p] = j;
            } else if (kn > 0 && dd < nds[kn - 1]) {
                size_t p = kn - 1;
                while (p > 0 && nds[p - 1] > dd) {
                    nds[p] = nds[p - 1]; row[p] = row[p - 1]; p--;
                }
                nds[p] = dd; row[p] = j;
            }
        }
        row[cnt++] = i;                   /* self */
        /* Re-sort the row by index so intersections can merge. cnt <= kn + 1 and
         * kn is small, so an insertion sort is the right tool. */
        for (size_t a = 1; a < cnt; a++) {
            size_t v = row[a], b = a;
            while (b > 0 && row[b - 1] > v) { row[b] = row[b - 1]; b--; }
            row[b] = v;
        }
        for (size_t a = cnt; a < w; a++) row[a] = (size_t)-1;   /* pad */
    }

    for (size_t i = 0; i < n; i++) {
        const size_t* ri = nbr + i * w;
        for (size_t j = i + 1; j < n; j++) {
            const size_t* rj = nbr + j * w;
            /* Mutual membership first: it is O(w) and rejects most pairs. */
            bool in_i = false, in_j = false;
            for (size_t a = 0; a < w; a++) {
                if (ri[a] == j) in_i = true;
                if (rj[a] == i) in_j = true;
            }
            if (!(in_i && in_j)) continue;
            size_t a = 0, b = 0, shared = 0;
            while (a < w && b < w && ri[a] != (size_t)-1 && rj[b] != (size_t)-1) {
                if (ri[a] == rj[b]) { shared++; a++; b++; }
                else if (ri[a] < rj[b]) a++;
                else b++;
            }
            if (shared >= thresh) fc_uf_union(uf, i, j);
        }
    }

    fc_assign_from_uf(uf, n, assign, k);
    free(nbr); free(nds); free(uf);
    return *k > 0;
}

static bool fc_method_jarvispatrick(const FcData* d, FcCount spec, const FcOpts* o,
                                    size_t* assign, size_t* k) {
    (void)spec;
    if (d->kind == FC_KIND_POINT) return fc_jarvis_ndim(d, o, assign, k);
    size_t n = d->n;
    double* sv = fc_sorted_values(d);
    if (!sv) return false;

    size_t kn = o->neighbor_count_given ? (size_t)o->neighbor_count : 5;
    if (kn > n - 1) kn = (n > 1) ? n - 1 : 0;
    size_t thresh = (kn + 1) / 2;
    if (thresh < 1) thresh = 1;

    size_t* lo = malloc(sizeof(size_t) * n);
    size_t* hi = malloc(sizeof(size_t) * n);
    size_t* id = malloc(sizeof(size_t) * n);
    if (!lo || !hi || !id) { free(sv); free(lo); free(hi); free(id); return false; }
    for (size_t j = 0; j < n; j++) fc_knn_window(sv, n, j, kn, &lo[j], &hi[j]);

    size_t cur = 0;
    id[0] = 0;
    for (size_t j = 1; j < n; j++) {
        bool mutual = (j >= lo[j - 1] && j < hi[j - 1]) &&
                      (j - 1 >= lo[j] && j - 1 < hi[j]);
        /* Shared neighbours = the overlap of the two contiguous windows. */
        size_t a = (lo[j] > lo[j - 1]) ? lo[j] : lo[j - 1];
        size_t b = (hi[j] < hi[j - 1]) ? hi[j] : hi[j - 1];
        size_t shared = (b > a) ? b - a : 0;
        if (!(mutual && shared >= thresh)) cur++;
        id[j] = cur;
    }

    fc_scatter(d, id, assign, k);
    free(sv); free(lo); free(hi); free(id);
    return *k > 0;
}

/* ------------------------------------------------------------------------- */
/* Method: GaussianMixture                                                    */
/* ------------------------------------------------------------------------- */

/* 1D EM over k Gaussians, with BIC choosing k from 1..k_max. Deterministic
 * quantile initialisation, as in fc_lloyd. Each point is assigned to its
 * highest-responsibility component.
 *
 * The variance floor is not defensive decoration: identical points give zero
 * variance and a singular Gaussian, so {7,7,7,7} would divide by zero without
 * it. It also prunes collapsed components, which is what makes a mixture
 * collapse to a single cluster on data that does not separate. */
static double fc_gmm_fit(const double* sv, size_t n, size_t k, double vfloor,
                         double* mu, double* var, double* w, size_t* id) {
    for (size_t j = 0; j < k; j++) {
        mu[j]  = sv[(size_t)(((double)j + 0.5) * (double)n / (double)k)];
        var[j] = vfloor;
        w[j]   = 1.0 / (double)k;
    }

    double* r = malloc(sizeof(double) * n * k);
    if (!r) return -INFINITY;

    double loglik = -INFINITY;
    for (int it = 0; it < FC_MAX_ITER; it++) {
        loglik = 0.0;
        for (size_t i = 0; i < n; i++) {
            double tot = 0.0;
            for (size_t j = 0; j < k; j++) {
                double z = sv[i] - mu[j];
                double p = w[j] * exp(-0.5 * z * z / var[j]) / sqrt(2.0 * M_PI * var[j]);
                r[i * k + j] = p;
                tot += p;
            }
            if (tot <= 0.0) { for (size_t j = 0; j < k; j++) r[i * k + j] = 1.0 / (double)k; tot = 1.0; }
            for (size_t j = 0; j < k; j++) r[i * k + j] /= tot;
            loglik += log(tot);
        }
        for (size_t j = 0; j < k; j++) {
            double sw = 0.0, sm = 0.0, sv2 = 0.0;
            for (size_t i = 0; i < n; i++) sw += r[i * k + j];
            if (sw <= 1e-12) { w[j] = 0.0; var[j] = vfloor; continue; }
            for (size_t i = 0; i < n; i++) sm += r[i * k + j] * sv[i];
            double m = sm / sw;
            for (size_t i = 0; i < n; i++) { double z = sv[i] - m; sv2 += r[i * k + j] * z * z; }
            mu[j]  = m;
            var[j] = sv2 / sw;
            if (var[j] < vfloor) var[j] = vfloor;
            w[j]   = sw / (double)n;
        }
    }

    for (size_t i = 0; i < n; i++) {
        size_t best = 0; double bp = -1.0;
        for (size_t j = 0; j < k; j++)
            if (r[i * k + j] > bp) { bp = r[i * k + j]; best = j; }
        id[i] = best;
    }
    free(r);
    return loglik;
}

static bool fc_method_gaussianmixture(const FcData* d, FcCount spec, const FcOpts* o,
                                      size_t* assign, size_t* k) {
    (void)spec; (void)o;
    size_t n = d->n;
    double* sv = fc_sorted_values(d);
    if (!sv) return false;

    /* Zero spread means every point is identical and one component is the only
     * sensible answer -- and the variance floor below would be zero. An empty
     * input (n == 0) has no sv[n-1] to read, so it folds into the same branch. */
    double spread = (n > 0) ? sv[n - 1] - sv[0] : 0.0;
    if (spread <= 0.0) {
        size_t* id0 = calloc(n, sizeof(size_t));
        if (!id0) { free(sv); return false; }
        fc_scatter(d, id0, assign, k);
        free(sv); free(id0);
        return *k > 0;
    }

    /* The variance floor is the squared median gap: a component may not be
     * narrower than the typical spacing between samples.
     *
     * This is load-bearing, not a guard against division by zero (though it is
     * that too, for identical points). A Gaussian mixture is unbounded above --
     * a component collapsing onto a single point drives its density, and hence
     * the likelihood, to infinity -- so with a floor set merely "small" the BIC
     * search buys arbitrarily many near-singular spikes and reports one cluster
     * per point. Measured: a floor of (spread/n)^2 * 1e-4 selected SIX
     * components for eight points. Flooring at the sampling scale says the
     * honest thing instead: structure finer than the point spacing is not
     * resolvable from this data. */
    double mg = fc_median_gap(sv, n);
    double vfloor = (mg > 0.0) ? mg * mg : (spread / (double)n) * (spread / (double)n);
    if (vfloor <= 0.0) vfloor = 1e-300;

    size_t kmax = d->n_distinct < 10 ? d->n_distinct : 10;
    double best_bic = INFINITY;
    size_t* best_id = malloc(sizeof(size_t) * n);
    size_t* id      = malloc(sizeof(size_t) * n);
    double* mu      = malloc(sizeof(double) * kmax);
    double* var     = malloc(sizeof(double) * kmax);
    double* w       = malloc(sizeof(double) * kmax);
    if (!best_id || !id || !mu || !var || !w) {
        free(sv); free(best_id); free(id); free(mu); free(var); free(w); return false;
    }
    for (size_t i = 0; i < n; i++) best_id[i] = 0;

    for (size_t kk = 1; kk <= kmax; kk++) {
        double ll = fc_gmm_fit(sv, n, kk, vfloor, mu, var, w, id);
        if (ll == -INFINITY) continue;
        double params = 3.0 * (double)kk - 1.0;          /* mu, var, weight */
        double bic = params * log((double)n) - 2.0 * ll;
        if (bic < best_bic) { best_bic = bic; memcpy(best_id, id, sizeof(size_t) * n); }
    }

    fc_scatter(d, best_id, assign, k);
    free(sv); free(best_id); free(id); free(mu); free(var); free(w);
    return *k > 0;
}

/* ------------------------------------------------------------------------- */
/* Method: Spectral                                                           */
/* ------------------------------------------------------------------------- */

/* Similarity graph, normalised Laplacian, cluster on the leading non-trivial
 * eigenvector.
 *
 * Honest note: in 1D a Gaussian similarity graph is nearly a path graph, so its
 * Fiedler vector largely recovers the sorted order and spectral clustering
 * degenerates toward the gap cuts of fc_method_gap. It is implemented as
 * specified rather than aliased, and an acceptance row pins whatever it does
 * rather than asserting it differs.
 *
 * The similarity matrix is n x n, so the method DECLINES above
 * FC_SPECTRAL_MAX_N rather than allocating gigabytes. Inverse iteration against
 * the known-constant first eigenvector is enough for the Fiedler vector here
 * and avoids pulling a full eigensolver onto this path. */
static bool fc_method_spectral(const FcData* d, FcCount spec, const FcOpts* o,
                               size_t* assign, size_t* k) {
    size_t n = d->n;
    if (n > FC_SPECTRAL_MAX_N) return false;
    if (n < 3) {                        /* nothing for a spectrum to say */
        size_t* id0 = calloc(n, sizeof(size_t));
        if (!id0) return false;
        fc_scatter(d, id0, assign, k);
        free(id0);
        return *k > 0;
    }

    double* sv = fc_sorted_values(d);
    if (!sv) return false;
    double h = o->radius_given ? o->radius : fc_median_gap(sv, n);
    if (h <= 0.0) h = 1.0;

    /* Fiedler vector of the path-like similarity graph by power iteration on
     * (I - L_sym), deflated against the constant vector. */
    double* deg = calloc(n, sizeof(double));
    double* v   = malloc(sizeof(double) * n);
    double* t   = malloc(sizeof(double) * n);
    size_t* id  = malloc(sizeof(size_t) * n);
    if (!deg || !v || !t || !id) { free(sv); free(deg); free(v); free(t); free(id); return false; }

    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++) {
            double z = (sv[i] - sv[j]) / h;
            deg[i] += exp(-0.5 * z * z);
        }

    for (size_t i = 0; i < n; i++) v[i] = ((i % 2) ? 1.0 : -1.0) + 0.001 * (double)i;

    for (int it = 0; it < FC_MAX_ITER; it++) {
        for (size_t i = 0; i < n; i++) {
            double acc = 0.0;
            for (size_t j = 0; j < n; j++) {
                double z = (sv[i] - sv[j]) / h;
                double wij = exp(-0.5 * z * z);
                acc += wij * v[j] / sqrt(deg[i] * deg[j]);
            }
            t[i] = acc;
        }
        /* Deflate the constant (trivial) eigenvector, then normalise. */
        double dot = 0.0, nrm = 0.0, cn = 0.0;
        for (size_t i = 0; i < n; i++) cn += deg[i];
        cn = sqrt(cn);
        for (size_t i = 0; i < n; i++) dot += t[i] * sqrt(deg[i]) / cn;
        for (size_t i = 0; i < n; i++) t[i] -= dot * sqrt(deg[i]) / cn;
        for (size_t i = 0; i < n; i++) nrm += t[i] * t[i];
        nrm = sqrt(nrm);
        if (nrm < 1e-300) break;
        for (size_t i = 0; i < n; i++) v[i] = t[i] / nrm;
    }

    /* The embedding is monotone along the line for a path-like graph, so the
     * natural cut is its widest jump -- the eigengap heuristic, applied to the
     * coordinates rather than the spectrum. */
    size_t natural = 1;
    double med = 0.0;
    {
        double* g = malloc(sizeof(double) * (n - 1));
        if (!g) { free(sv); free(deg); free(v); free(t); free(id); return false; }
        for (size_t j = 0; j + 1 < n; j++) g[j] = fabs(v[j + 1] - v[j]);
        double* s = malloc(sizeof(double) * (n - 1));
        if (!s) { free(g); free(sv); free(deg); free(v); free(t); free(id); return false; }
        memcpy(s, g, sizeof(double) * (n - 1));
        for (size_t i = 1; i + 1 < n; i++) {
            double x = s[i]; size_t j = i;
            while (j > 0 && s[j - 1] > x) { s[j] = s[j - 1]; j--; }
            s[j] = x;
        }
        size_t m = n - 1;
        med = (m % 2) ? s[m / 2] : 0.5 * (s[m / 2 - 1] + s[m / 2]);
        for (size_t j = 0; j + 1 < n; j++) if (g[j] > (double)FC_GAP_FACTOR * med) natural++;
        free(g); free(s);
    }

    size_t target = fc_target_count(natural, spec, d->n_distinct);

    /* Cut the target-1 widest embedding jumps. */
    {
        size_t m = n - 1;
        double* g = malloc(sizeof(double) * m);
        size_t* gi = malloc(sizeof(size_t) * m);
        bool* cut = calloc(m ? m : 1, sizeof(bool));
        if (!g || !gi || !cut) { free(g); free(gi); free(cut); free(sv); free(deg); free(v); free(t); free(id); return false; }
        for (size_t j = 0; j < m; j++) { g[j] = fabs(v[j + 1] - v[j]); gi[j] = j; }
        for (size_t i = 1; i < m; i++) {
            size_t x = gi[i]; size_t j = i;
            while (j > 0 && g[gi[j - 1]] > g[x]) { gi[j] = gi[j - 1]; j--; }
            gi[j] = x;
        }
        size_t ncut = target > 0 ? target - 1 : 0;
        if (ncut > m) ncut = m;
        for (size_t s2 = 0; s2 < ncut; s2++) cut[gi[m - 1 - s2]] = true;
        size_t cur = 0;
        id[0] = 0;
        for (size_t j = 1; j < n; j++) { if (cut[j - 1]) cur++; id[j] = cur; }
        free(g); free(gi); free(cut);
    }

    fc_scatter(d, id, assign, k);
    free(sv); free(deg); free(v); free(t); free(id);
    return *k > 0;
}

/* ------------------------------------------------------------------------- */
/* Option and count parsing                                                   */
/* ------------------------------------------------------------------------- */

static bool fc_is_known_option_name(const char* s) {
    return s == SYM_Method
        || s == SYM_DistanceFunction
        || s == SYM_CriterionFunction
        || s == SYM_PerformanceGoal;
}

static bool fc_is_option_arg(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* lhs = e->data.function.args[0];
    if (lhs->type != EXPR_SYMBOL) return false;
    return fc_is_known_option_name(lhs->data.symbol.name);
}

static bool fc_lookup_method(const char* s, FcMethod* out) {
    for (size_t i = 0; i < sizeof(FC_METHOD_NAMES) / sizeof(FC_METHOD_NAMES[0]); i++)
        if (strcmp(s, FC_METHOD_NAMES[i].name) == 0) { *out = FC_METHOD_NAMES[i].m; return true; }
    return false;
}

/* Read a real-valued suboption into a double. */
static bool fc_read_double(Expr* v, double* out) {
    if (!list_real_number_q(v)) return false;
    *out = fc_to_double(v);
    return true;
}

/* Method -> "Name" or Method -> {"Name", subopt -> v, ...}. */
static bool fc_parse_method(Expr* rhs, FcOpts* o) {
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol.name == SYM_Automatic) {
        o->method = FC_AGGLOMERATE;    /* see fc_resolve_automatic_method */
        return true;
    }
    if (rhs->type == EXPR_STRING) return fc_lookup_method(rhs->data.string, &o->method);

    if (is_listq(rhs) && rhs->data.function.arg_count >= 1) {
        Expr* nm = rhs->data.function.args[0];
        if (nm->type != EXPR_STRING) return false;
        if (!fc_lookup_method(nm->data.string, &o->method)) return false;
        for (size_t i = 1; i < rhs->data.function.arg_count; i++) {
            Expr* r = rhs->data.function.args[i];
            if (r->type != EXPR_FUNCTION || r->data.function.arg_count != 2) return false;
            const char* h = r->data.function.head->data.symbol.name;
            if (h != SYM_Rule && h != SYM_RuleDelayed) return false;
            Expr* key = r->data.function.args[0];
            Expr* val = r->data.function.args[1];
            if (key->type != EXPR_SYMBOL && key->type != EXPR_STRING) return false;
            const char* kn = (key->type == EXPR_SYMBOL) ? key->data.symbol.name
                                                        : key->data.string;
            if (strcmp(kn, "NeighborhoodRadius") == 0) {
                if (!fc_read_double(val, &o->radius) || o->radius <= 0.0) return false;
                o->radius_given = true;
            } else if (strcmp(kn, "MinPoints") == 0) {
                if (val->type != EXPR_INTEGER || val->data.integer < 1) return false;
                o->min_points = (long)val->data.integer;
                o->min_points_given = true;
            } else if (strcmp(kn, "NeighborCount") == 0) {
                if (val->type != EXPR_INTEGER || val->data.integer < 1) return false;
                o->neighbor_count = (long)val->data.integer;
                o->neighbor_count_given = true;
            } else {
                return false;               /* unknown suboption: decline */
            }
        }
        return true;
    }
    return false;
}

/* On a LINE every accepted metric is a monotone transform of |a - b|, so all four
 * induce the same ordering of gaps and the same partition -- which is why this
 * used to validate the name and store nothing at all.
 *
 * That stops being true the moment the points have more than one component:
 * Manhattan and Euclidean rank pairs differently in the plane, and the option
 * silently doing nothing there would be a wrong answer rather than a missing
 * feature. The choice is now recorded and applied for FC_KIND_POINT. */
static bool fc_parse_distance_function(Expr* rhs, FcOpts* o) {
    if (rhs->type == EXPR_SYMBOL) {
        const char* s = rhs->data.symbol.name;
        if (s == SYM_Automatic)                              { o->dist = FC_DIST_AUTOMATIC;        return true; }
        if (strcmp(s, "EuclideanDistance") == 0)             { o->dist = FC_DIST_EUCLIDEAN;        return true; }
        if (strcmp(s, "SquaredEuclideanDistance") == 0)      { o->dist = FC_DIST_SQUAREDEUCLIDEAN; return true; }
        if (strcmp(s, "ManhattanDistance") == 0)             { o->dist = FC_DIST_MANHATTAN;        return true; }
    }
    return false;
}

static bool fc_apply_option(Expr* rule, FcOpts* o) {
    const char* name = rule->data.function.args[0]->data.symbol.name;
    Expr* rhs = rule->data.function.args[1];

    if (name == SYM_Method) return fc_parse_method(rhs, o);
    if (name == SYM_DistanceFunction) return fc_parse_distance_function(rhs, o);
    /* CriterionFunction and PerformanceGoal are accepted and currently have no
     * effect; the docs state this rather than leaving it implicit. */
    if (name == SYM_CriterionFunction || name == SYM_PerformanceGoal) return true;
    return false;
}

/* Decode the optional count argument. Returns false for a spec this builtin
 * does not understand, which the caller turns into an unevaluated result. */
static bool fc_parse_count(Expr* e, FcCount* out) {
    if (e == NULL) { out->mode = FC_COUNT_AUTOMATIC; out->n = 0; return true; }
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Automatic) {
        out->mode = FC_COUNT_AUTOMATIC; out->n = 0; return true;
    }
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer < 1) return false;
        out->mode = FC_COUNT_FIXED;
        out->n = (size_t)e->data.integer;
        return true;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_UpTo &&
        e->data.function.arg_count == 1) {
        Expr* a = e->data.function.args[0];
        if (a->type != EXPR_INTEGER || a->data.integer < 1) return false;
        out->mode = FC_COUNT_BOUNDED;
        out->n = (size_t)a->data.integer;
        return true;
    }
    return false;
}

/* Method -> Automatic resolves to Agglomerate in every count mode: it is the
 * only family valid in all three, it is deterministic, and its single tunable
 * is FC_GAP_FACTOR. Deliberately NOT a criterion-driven search across methods
 * -- that is Mathematica's unpublished internal index, and this file does not
 * pretend to reproduce it. */
static FcMethod fc_resolve_automatic_method(FcCount spec) {
    (void)spec;
    return FC_AGGLOMERATE;
}

/* ------------------------------------------------------------------------- */
/* Builtin                                                                    */
/* ------------------------------------------------------------------------- */

typedef bool (*FcMethodFn)(const FcData*, FcCount, const FcOpts*, size_t*, size_t*);

static FcMethodFn fc_method_fn(FcMethod m) {
    switch (m) {
        case FC_AGGLOMERATE:
        case FC_SPANNINGTREE:              return fc_method_gap;
        case FC_KMEANS:                    return fc_method_kmeans;
        case FC_KMEDOIDS:                  return fc_method_kmedoids;
        case FC_SPECTRAL:                  return fc_method_spectral;
        case FC_DBSCAN:                    return fc_method_dbscan;
        case FC_GAUSSIANMIXTURE:           return fc_method_gaussianmixture;
        case FC_JARVISPATRICK:             return fc_method_jarvispatrick;
        case FC_MEANSHIFT:                 return fc_method_meanshift;
        case FC_NEIGHBORHOODCONTRACTION:   return fc_method_neighborhood;
        default:                           return NULL;
    }
}

/* The whole of FindClusters, with the data list passed in rather than read off
 * `res`.
 *
 * Split out purely for ownership. There are eighteen `return NULL` paths below,
 * and the visible-NDArray surface has to materialise its argument into a nested
 * List which must outlive fc_emit_clusters (that reads the elements to build the
 * result) and be freed on every one of those paths. Threading a free through
 * eighteen exits is how a leak or a double free gets in; one caller that owns the
 * temporary and one callee that only borrows it cannot go wrong. */
static Expr* fc_find_clusters(Expr* res, Expr* list) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    /* Split trailing options from the positional arguments. */
    FcOpts opts = { FC_AGGLOMERATE, 0.0, false, 0, false, 0, false, FC_DIST_AUTOMATIC };
    bool method_given = false;
    size_t n_pos = argc;
    for (size_t i = 1; i < argc; i++) {
        if (fc_is_option_arg(res->data.function.args[i])) {
            if (i < n_pos) n_pos = i;
        }
    }
    for (size_t i = n_pos; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        if (!fc_is_option_arg(a)) return NULL;    /* junk after the options */
        if (a->data.function.args[0]->data.symbol.name == SYM_Method) method_given = true;
        if (!fc_apply_option(a, &opts)) return NULL;
    }
    if (n_pos < 1 || n_pos > 2) return NULL;

    if (!is_listq(list)) return NULL;

    size_t n = list->data.function.arg_count;
    /* Empty in, unevaluated out -- matching Mathematica's ::mlmpty, and
     * defensible on its own terms: no clusters can be formed from nothing. */
    if (n == 0) return NULL;

    FcCount spec;
    if (!fc_parse_count(n_pos == 2 ? res->data.function.args[1] : NULL, &spec)) return NULL;

    if (!method_given) opts.method = fc_resolve_automatic_method(spec);
    if (!FC_ALLOWED[opts.method][spec.mode]) return NULL;

    FcMethodFn fn = fc_method_fn(opts.method);
    if (!fn) return NULL;

    Expr** elem = list->data.function.args;

    FcData d;
    memset(&d, 0, sizeof d);
    d.elem = elem;
    d.n = n;
    d.n_gap = n - 1;

    /* Shape gate. Every element must be a real scalar, or every element a vector
     * of the same length over real components. One symbolic element declines the
     * whole call rather than being dropped or treated as a nominal feature. */
    if (!fc_probe_shape(elem, n, &d.dim, &d.kind)) return NULL;

    /* Which methods are dimension-general, as opposed to still reading the sorted
     * projection. The rest are declined above one dimension rather than being
     * silently run on meaningless data -- d->val and d->order are both NULL there,
     * so "silently" would in fact be a crash.
     *
     * This list grows one method at a time as each is ported. Kept as an explicit
     * check on the function pointer rather than a flag on FcMethod so that a new
     * method cannot default into being considered ported. */
    bool fn_is_ndim = (fn == fc_method_gap)
                   || (fn == fc_method_meanshift)
                   || (fn == fc_method_neighborhood)
                   || (fn == fc_method_kmeans)
                   || (fn == fc_method_dbscan)
                   || (fn == fc_method_jarvispatrick)
                   || (fn == fc_method_kmedoids);
    /* The shift methods need coordinates; strings have none, so they stay
     * declined for those even though POINT input is now fine. */
    if (d.kind == FC_KIND_SEQUENCE && fn != fc_method_gap) return NULL;
    if (d.kind != FC_KIND_SCALAR && !fn_is_ndim) return NULL;

    /* Vector work is quadratic (Prim, and the neighbourhood queries that a sort
     * makes linear in 1D), so it is capped the way the other quadratic methods
     * already are. */
    /* Which point builder will run, decided here because it sets the ceiling:
     * both are quadratic, but their constants differ by more than two orders of
     * magnitude. Machine-precision input carries no precision a double would
     * lose, so nothing is given up by taking the fast path. */
    bool machine_points = (d.kind == FC_KIND_POINT) && fc_all_machine(&d);
    if (d.kind != FC_KIND_SCALAR) {
        size_t cap = machine_points ? FC_NDIM_MACHINE_MAX_N : FC_NDIM_MAX_N;
        if (n > cap) return NULL;
    }

    /* The metric the tree will be weighted by. Set before either builder runs,
     * since both read it, and only meaningful for POINT input -- the scalar path
     * below builds plain exact differences and the sequence path edit distances,
     * neither of which the option can change. */
    d.dist = (d.kind == FC_KIND_POINT) ? opts.dist : FC_DIST_AUTOMATIC;

    if (d.n_gap) {
        d.gap = calloc(d.n_gap, sizeof(Expr*));
        d.eu  = malloc(sizeof(size_t) * d.n_gap);
        d.ev  = malloc(sizeof(size_t) * d.n_gap);
        d.bnd = calloc(d.n_gap, sizeof(bool));
        if (!d.gap || !d.eu || !d.ev || !d.bnd) { fc_data_free(&d); return NULL; }
    }

    if (d.kind == FC_KIND_SCALAR) {
        /* The MST of points on a line is the sorted adjacency chain, so sorting
         * builds it in O(n log n) and the weights are plain exact differences --
         * no Abs needed, sorted neighbours differ non-negatively by
         * construction. This is the original code path, unchanged, and the reason
         * one-dimensional results are identical to before. */
        d.thresh_factor = FC_GAP_FACTOR;
        d.order = malloc(sizeof(size_t) * n);
        d.val   = malloc(sizeof(double) * n);
        if (!d.order || !d.val) { fc_data_free(&d); return NULL; }
        for (size_t i = 0; i < n; i++) { d.order[i] = i; d.val[i] = fc_to_double(elem[i]); }

        FcSortCtx sc = { elem, true };
        fc_merge_sort(d.order, n, &sc);
        if (!sc.ok) { fc_data_free(&d); return NULL; }

        for (size_t j = 0; j < d.n_gap; j++) {
            d.eu[j] = d.order[j];
            d.ev[j] = d.order[j + 1];
            Expr* a[2] = { expr_copy(elem[d.ev[j]]), expr_copy(elem[d.eu[j]]) };
            d.gap[j] = eval_and_free(internal_subtract(a, 2));
            if (!list_real_number_q(d.gap[j])) { fc_data_free(&d); return NULL; }
        }
    } else {
        /* POINT weights are SQUARED distances under the Euclidean metrics, so the
         * threshold factor is squared to match -- see the FcData comment. Under
         * Manhattan the weights are plain sums of absolute differences, already
         * linear, so the factor stays plain: squaring it there would silently
         * apply a threshold nine times too large. SEQUENCE weights are edit
         * distances, also linear. */
        d.thresh_factor = (d.kind == FC_KIND_POINT && fc_dist_is_squared(d.dist))
                        ? FC_GAP_FACTOR * FC_GAP_FACTOR
                        : FC_GAP_FACTOR;

        /* Coordinates exist only where there are coordinates. Strings have none,
         * and the gap methods never ask for any -- they work purely off the exact
         * pairwise distances in the tree. */
        if (d.kind == FC_KIND_POINT) {
            d.coord = malloc(sizeof(double) * n * d.dim);
            if (!d.coord) { fc_data_free(&d); return NULL; }
            for (size_t i = 0; i < n; i++)
                for (size_t c = 0; c < d.dim; c++)
                    d.coord[i * d.dim + c] = fc_to_double(fc_comp(&d, i, c));
        }

        /* Machine-precision points take the double builder; exact ones (a
         * Rational, a bigint, an MPFR value) keep the exact builder, which is the
         * only one that can order them correctly. */
        if (d.n_gap && !(machine_points ? fc_build_mst_machine(&d)
                                        : fc_build_mst(&d))) {
            fc_data_free(&d);
            return NULL;
        }
    }

    /* Exact boundary flags and the distinct count, one pass over the tree edges.
     *
     * Distinctness comes from comparing the ELEMENTS with the exact comparator,
     * never from the sign or magnitude of an edge weight. The 1D weights are
     * built with internal_subtract, which widens a mixed exact/inexact pair to a
     * double: (2^60 + 1) - 2.0^60 is 0.0, so a weight test would call two
     * different values equal, collapse n_distinct, and make the whole result
     * depend on input order.
     *
     * Counting distinct points this way is exact in any dimension. Zero-weight
     * edges are precisely those joining identical points, and a spanning tree
     * restricted to a class of m identical points uses m-1 of them, so
     * 1 + (number of nonzero edges) is the number of distinct points. */
    d.n_distinct = 1;
    for (size_t j = 0; j < d.n_gap; j++) {
        bool ok = true;
        if (!fc_elem_equal(&d, d.eu[j], d.ev[j], &ok)) {
            d.bnd[j] = true;
            d.n_distinct++;
        }
        if (!ok) { fc_data_free(&d); return NULL; }
    }

    size_t* assign = calloc(n, sizeof(size_t));
    if (!assign) { fc_data_free(&d); return NULL; }

    size_t k = 0;
    bool ok = fn(&d, spec, &opts, assign, &k);

    /* THE EQUAL-ELEMENTS INVARIANT, enforced once here rather than trusted to
     * ten independent methods.
     *
     * No method may put two equal elements in different clusters -- the
     * docstring states it unconditionally, and the distinct-value cap depends
     * on it. Three methods broke it in review: JarvisPatrick returned four
     * clusters for nine copies of 7 (its k-NN window walks left on a zero
     * distance, so a point is not its own neighbour's neighbour), and DBSCAN
     * split a duplicated pair whenever MinPoints exceeded its multiplicity.
     * Both are real bugs in those kernels, but patching each one leaves the
     * next method free to reintroduce it, and the property is cheap to enforce
     * globally: walk the sorted order and, wherever the exact gap is zero,
     * fold the later element into the earlier one's cluster.
     *
     * This can only ever REDUCE the cluster count, and only to the distinct
     * count, which is the documented cap. fc_scatter then renumbers, so the
     * ids stay contiguous. */
    if (ok && d.n_gap) {
        /* Ascending edge order with eu = parent means a point's representative is
         * already folded when its own edge is processed, so one pass is
         * transitive across a whole run of identical points. In 1D these are the
         * sorted-chain edges, so this is the same walk as before. */
        for (size_t j = 0; j < d.n_gap; j++)
            if (!d.bnd[j]) assign[d.ev[j]] = assign[d.eu[j]];
        size_t* renum = malloc(sizeof(size_t) * n);
        if (!renum) { ok = false; }
        else {
            /* Renumber contiguously by first occurrence in INPUT order, the
             * same ordering fc_emit_clusters applies. */
            for (size_t i = 0; i < n; i++) renum[i] = SIZE_MAX;
            size_t kk = 0;
            for (size_t i = 0; i < n; i++)
                if (renum[assign[i]] == SIZE_MAX) renum[assign[i]] = kk++;
            for (size_t i = 0; i < n; i++) assign[i] = renum[assign[i]];
            k = kk;
            free(renum);
        }
    }

    fc_data_free(&d);
    if (!ok || k == 0) { free(assign); return NULL; }

    Expr* result = fc_emit_clusters(elem, n, assign, k);
    free(assign);
    return result;
}

Expr* builtin_find_clusters(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 1) return NULL;
    Expr* a0 = res->data.function.args[0];

    /* A VISIBLE NDArray. Not a fast path -- a materialise guard, the same shape
     * Cases and FlattenAt use (see tools/check_packed_aware.py's EXEMPT list).
     *
     * Why materialise rather than read the buffer. Everything downstream is
     * Expr-centric for reasons that are not incidental: the exact MST orders
     * Rationals and bigints through internal_subtract, the boundary set compares
     * ELEMENTS because 2^60 and 2^60+1 are distinct yet project to the same
     * double, and the result is built from the input elements themselves. A
     * buffer path would have to give all of that up, and the exactness is the
     * reason one-dimensional answers are exact. Machine speed is already had:
     * every value in an NDArray is machine by construction, so fc_all_machine is
     * true and fc_build_mst_machine runs -- the same fast builder a machine List
     * takes.
     *
     * The packed-List surface needs nothing here. The transparency gate
     * (eval.c step 2.7) materialises it for any head not on pack.c's AWARE list,
     * and FindClusters is deliberately not on it. The gate tests only
     * is_packed_list, though, which is why the visible surface arrived here
     * untouched and returned unevaluated before this. */
    if (is_ndarray(a0)) {
        Expr* nested = ndarray_to_nested_list(a0);
        if (!nested) return NULL;
        Expr* out = fc_find_clusters(res, nested);
        expr_free(nested);
        return out;
    }
    return fc_find_clusters(res, a0);
}
