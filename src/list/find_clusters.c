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

typedef struct {
    FcMethod method;
    double   radius;          bool radius_given;
    long     min_points;      bool min_points_given;
    long     neighbor_count;  bool neighbor_count_given;
} FcOpts;

/* ------------------------------------------------------------------------- */
/* Decoded input                                                             */
/* ------------------------------------------------------------------------- */

/* Everything the methods share, computed once.
 *
 * `order` is the EXACT sorted permutation and `gap` the EXACT adjacent
 * differences; `val` is the machine projection, for the methods whose output is
 * inexact by definition. Keeping both is what lets the default gap method stay
 * exact while KMeans and the density family use doubles. */
typedef struct {
    Expr**  elem;      /* borrowed: the input elements, input order */
    size_t  n;
    size_t* order;     /* owned: sorted permutation of 0..n-1 */
    double* val;       /* owned: val[i] is the double projection of elem[i] */
    Expr**  gap;       /* owned: gap[j] = elem[order[j+1]] - elem[order[j]] */
    size_t  n_gap;     /* n - 1, or 0 when n == 0 */
    size_t  n_distinct;
    /* owned: bnd[j] is true when gap j is EXACTLY nonzero, i.e. sorted
     * positions j and j+1 hold different values. Derived from the gap Exprs,
     * never from the double projection -- `2^60` and `2^60 + 1` are distinct
     * but project to the same double, and a boundary set computed in double
     * space silently loses them. Every consumer of "is this a real boundary"
     * must read this, not compare val[]. */
    bool*   bnd;
} FcData;

static void fc_data_free(FcData* d) {
    if (d->gap) {
        for (size_t j = 0; j < d->n_gap; j++) expr_free(d->gap[j]);
        free(d->gap);
    }
    free(d->order);
    free(d->val);
    free(d->bnd);
    d->gap = NULL; d->order = NULL; d->val = NULL; d->bnd = NULL;
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

    Expr* mul[2] = { expr_new_integer(FC_GAP_FACTOR), expr_copy(median) };
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

    /* Walk the sorted order, opening a new cluster after every cut gap. */
    size_t id = 0;
    for (size_t j = 0; j < n; j++) {
        assign[d->order[j]] = id;
        if (j + 1 < n && cut[j]) id++;
    }
    free(cut);

    *k = id + 1;
    return true;
}

/* ------------------------------------------------------------------------- */
/* Sorted-neighbourhood kernel                                                */
/* ------------------------------------------------------------------------- */

/* The primitives the density family shares. Every one is O(1), O(k) or O(n)
 * ONLY because the array is sorted -- the same neighbourhood queries in general
 * dimension are the expensive part of DBSCAN, mean shift and Jarvis-Patrick.
 * `sv` is always the values in sorted order, i.e. sv[j] = val[order[j]]. */

/* Half-open [*lo, *hi) : the indices within `eps` of sv[i]. */
static void fc_eps_window(const double* sv, size_t n, size_t i, double eps,
                          size_t* lo, size_t* hi) {
    size_t a = i;
    while (a > 0 && sv[i] - sv[a - 1] <= eps) a--;
    size_t b = i + 1;
    while (b < n && sv[b] - sv[i] <= eps) b++;
    *lo = a; *hi = b;
}

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

/* Merge contiguous runs of a sorted, monotone "converged position" array into
 * cluster ids: two points join when their converged positions differ by less
 * than `tol`. Shared by MeanShift and NeighborhoodContraction, whose update
 * rules differ but whose merge step is identical. */
static void fc_merge_modes(const double* pos, size_t n, double tol, size_t* id) {
    size_t cur = 0;
    id[0] = 0;
    for (size_t j = 1; j < n; j++) {
        if (fabs(pos[j] - pos[j - 1]) > tol) cur++;
        id[j] = cur;
    }
}

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

static bool fc_method_kmeans(const FcData* d, FcCount spec, const FcOpts* o,
                             size_t* assign, size_t* k) {
    (void)o; return fc_lloyd(d, spec, assign, k, false);
}

static bool fc_method_kmedoids(const FcData* d, FcCount spec, const FcOpts* o,
                               size_t* assign, size_t* k) {
    (void)o; return fc_lloyd(d, spec, assign, k, true);
}

/* ------------------------------------------------------------------------- */
/* Method: DBSCAN                                                             */
/* ------------------------------------------------------------------------- */

/* Density-based, Automatic count only. A point is a core point when its
 * eps-window holds at least MinPoints members; core points whose windows
 * overlap join. A point in no dense region is NOISE, and rather than dropping
 * it -- which would lose an input element -- it becomes its own singleton
 * cluster, so the result is always a partition. */
static bool fc_method_dbscan(const FcData* d, FcCount spec, const FcOpts* o,
                             size_t* assign, size_t* k) {
    (void)spec;
    size_t n = d->n;
    double* sv = fc_sorted_values(d);
    if (!sv) return false;

    double eps = o->radius_given ? o->radius : (double)FC_GAP_FACTOR * fc_median_gap(sv, n);
    size_t minpts = o->min_points_given ? (size_t)o->min_points : 2;

    bool*   core = calloc(n, sizeof(bool));
    size_t* id   = malloc(sizeof(size_t) * n);
    size_t* lo   = malloc(sizeof(size_t) * n);
    size_t* hi   = malloc(sizeof(size_t) * n);
    if (!core || !id || !lo || !hi) { free(sv); free(core); free(id); free(lo); free(hi); return false; }

    for (size_t j = 0; j < n; j++) {
        fc_eps_window(sv, n, j, eps, &lo[j], &hi[j]);
        core[j] = (hi[j] - lo[j]) >= minpts;
    }

    /* One left-to-right sweep: stay in the current cluster while consecutive
     * points are eps-reachable through a core point. */
    size_t cur = 0;
    id[0] = 0;
    for (size_t j = 1; j < n; j++) {
        bool linked = (sv[j] - sv[j - 1] <= eps) && (core[j] || core[j - 1]);
        if (!linked) cur++;
        id[j] = cur;
    }

    fc_scatter(d, id, assign, k);
    free(sv); free(core); free(id); free(lo); free(hi);
    return *k > 0;
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
    size_t n = d->n;
    if (n > FC_SHIFT_MAX_N) return false;      /* quadratic; see FC_SHIFT_MAX_N */
    double* sv = fc_sorted_values(d);
    if (!sv) return false;

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
    double scale = fc_median_gap(sv, n);

    /* A zero median gap does NOT mean the points are identical -- it means at
     * least half the adjacent gaps are zero, which any tie-heavy input has.
     * Collapsing to one cluster there discarded real structure. */
    if (scale <= 0.0 && n > 1 && sv[n - 1] > sv[0])
        scale = (sv[n - 1] - sv[0]) / (double)(n - 1);

    double h = o->radius_given ? o->radius
                               : (flat_kernel ? (double)FC_GAP_FACTOR : 1.0) * scale;
    if (o->radius_given) scale = o->radius;

    if (h <= 0.0) {                     /* genuinely no spread: one cluster */
        size_t* id0 = calloc(n, sizeof(size_t));
        if (!id0) { free(sv); return false; }
        fc_scatter(d, id0, assign, k);
        free(sv); free(id0);
        return *k > 0;
    }

    double* pos = malloc(sizeof(double) * n);
    double* nxt = malloc(sizeof(double) * n);
    size_t* id  = malloc(sizeof(size_t) * n);
    if (!pos || !nxt || !id) { free(sv); free(pos); free(nxt); free(id); return false; }
    memcpy(pos, sv, sizeof(double) * n);

    for (int it = 0; it < FC_MAX_ITER; it++) {
        double delta = 0.0;
        for (size_t j = 0; j < n; j++) {
            double num = 0.0, den = 0.0;
            if (flat_kernel) {
                for (size_t t = 0; t < n; t++)
                    if (fabs(sv[t] - pos[j]) <= h) { num += sv[t]; den += 1.0; }
            } else {
                for (size_t t = 0; t < n; t++) {
                    double z = (pos[j] - sv[t]) / h;
                    double w = exp(-0.5 * z * z);
                    num += w * sv[t]; den += w;
                }
            }
            nxt[j] = (den > 0.0) ? num / den : pos[j];
            double mv = fabs(nxt[j] - pos[j]);
            if (mv > delta) delta = mv;
        }
        memcpy(pos, nxt, sizeof(double) * n);
        if (delta < h * 1e-6) break;
    }

    fc_merge_modes(pos, n, scale * FC_MERGE_SLACK, id);
    fc_scatter(d, id, assign, k);
    free(sv); free(pos); free(nxt); free(id);
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
static bool fc_method_jarvispatrick(const FcData* d, FcCount spec, const FcOpts* o,
                                    size_t* assign, size_t* k) {
    (void)spec;
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
     * sensible answer -- and the variance floor below would be zero. */
    double spread = sv[n - 1] - sv[0];
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

/* In 1D every accepted distance is a monotone transform of |a - b|, so all of
 * them induce the same ordering of gaps and therefore the same partition for
 * every distance-ranking method. Accepting four names is not four
 * implementations, and the docs say so. */
static bool fc_parse_distance_function(Expr* rhs) {
    if (rhs->type == EXPR_SYMBOL) {
        const char* s = rhs->data.symbol.name;
        return s == SYM_Automatic
            || strcmp(s, "EuclideanDistance") == 0
            || strcmp(s, "ManhattanDistance") == 0
            || strcmp(s, "SquaredEuclideanDistance") == 0;
    }
    return false;
}

static bool fc_apply_option(Expr* rule, FcOpts* o) {
    const char* name = rule->data.function.args[0]->data.symbol.name;
    Expr* rhs = rule->data.function.args[1];

    if (name == SYM_Method) return fc_parse_method(rhs, o);
    if (name == SYM_DistanceFunction) return fc_parse_distance_function(rhs);
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

Expr* builtin_find_clusters(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    /* Split trailing options from the positional arguments. */
    FcOpts opts = { FC_AGGLOMERATE, 0.0, false, 0, false, 0, false };
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

    Expr* list = res->data.function.args[0];
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

    /* Numeric gate: every element must be a real number we can order exactly.
     * One symbolic element declines the whole call rather than being dropped or
     * treated as a nominal feature. */
    for (size_t i = 0; i < n; i++)
        if (!list_real_number_q(elem[i])) return NULL;

    FcData d;
    memset(&d, 0, sizeof d);
    d.elem = elem;
    d.n = n;
    d.n_gap = n - 1;
    d.order = malloc(sizeof(size_t) * n);
    d.val   = malloc(sizeof(double) * n);
    if (!d.order || !d.val) { fc_data_free(&d); return NULL; }
    for (size_t i = 0; i < n; i++) { d.order[i] = i; d.val[i] = fc_to_double(elem[i]); }

    FcSortCtx sc = { elem, true };
    fc_merge_sort(d.order, n, &sc);
    if (!sc.ok) { fc_data_free(&d); return NULL; }

    /* Exact adjacent gaps over the sorted order. No Abs: sorted neighbours give
     * a non-negative difference by construction. */
    if (d.n_gap) {
        d.gap = calloc(d.n_gap, sizeof(Expr*));
        if (!d.gap) { fc_data_free(&d); return NULL; }
        for (size_t j = 0; j < d.n_gap; j++) {
            Expr* a[2] = { expr_copy(elem[d.order[j + 1]]), expr_copy(elem[d.order[j]]) };
            d.gap[j] = eval_and_free(internal_subtract(a, 2));
            if (!list_real_number_q(d.gap[j])) { fc_data_free(&d); return NULL; }
        }
    }

    /* Exact boundary flags and the distinct count, from one pass over the gap
     * signs. Both are exact by construction; nothing downstream needs to
     * rediscover "are these two values different" from the double projection. */
    d.n_distinct = 1;
    if (d.n_gap) {
        d.bnd = calloc(d.n_gap, sizeof(bool));
        if (!d.bnd) { fc_data_free(&d); return NULL; }
    }
    for (size_t j = 0; j < d.n_gap; j++) {
        /* Distinctness comes from comparing the ELEMENTS with the exact
         * comparator, never from the sign of their difference. The gap Exprs
         * are built with internal_subtract, which widens a mixed exact/inexact
         * pair to a double: (2^60 + 1) - 2.0^60 is 0.0, so a gap-sign test
         * would call two different values equal, collapse n_distinct, and make
         * the whole result depend on input order. */
        bool ok = true;
        if (list_numeric_cmp(elem[d.order[j + 1]], elem[d.order[j]], &ok) != 0) {
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
        for (size_t j = 0; j < d.n_gap; j++)
            if (!d.bnd[j]) assign[d.order[j + 1]] = assign[d.order[j]];
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
