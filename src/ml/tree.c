/* tree.c -- CART classification trees.
 *
 * GINI RATHER THAN ENTROPY, and it is not a coin toss. Both rank splits almost identically in
 * practice, but Gini is a sum of squares where entropy is a sum of x log x -- no logarithm per
 * candidate threshold, and no special case at p = 0 where log p is undefined. A splitter
 * evaluates the criterion O(dim * n) times per node, so the cheaper one with no domain edge is
 * the better default. The choice is visible in the docs, not hidden.
 *
 * THE SWEEP. For a node and a feature, sorting the node's points by that feature turns finding
 * the best threshold into a single pass: walk the sorted order moving one point at a time from
 * the right child to the left, updating both class histograms incrementally, and evaluate the
 * weighted impurity only where the feature value actually changes. That is O(n log n) per
 * feature per node from the sort, and O(n) for the sweep itself -- as against O(n) per
 * candidate threshold if the histograms were recomputed each time, which is the obvious
 * implementation and is quadratic.
 *
 * Thresholds are MIDPOINTS between consecutive distinct values, not the values themselves. A
 * threshold sitting exactly on a training value makes the `<=` boundary depend on floating
 * point equality with that value, so an unseen point equal to it lands by luck; the midpoint
 * puts the boundary in the gap where nothing lies.
 *
 * DETERMINISM IS A REQUIREMENT, NOT A COURTESY. Ties on impurity decrease are broken by lower
 * feature index, then lower threshold. Without that a tie could resolve differently between
 * runs (or between platforms, since qsort is not required to be stable), and then no test
 * could pin a tree at all -- which is the same reasoning that ruled out a softmax for
 * multi-class logistic regression, where the parameters are not unique either. To keep the
 * sort itself from deciding anything, its comparator falls back to the point INDEX when two
 * feature values are equal, making the order total and reproducible.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tree.h"

/* ---- the sort key: (value, index), so the order is total ---- */

typedef struct { double v; size_t i; } TreeKey;

static int treekey_cmp(const void* a, const void* b) {
    const TreeKey* p = (const TreeKey*)a;
    const TreeKey* q = (const TreeKey*)b;
    if (p->v < q->v) return -1;
    if (p->v > q->v) return 1;
    /* Equal values fall back to the index. qsort is not required to be stable, so without
     * this two runs could order a tie differently and the fitted tree would not be
     * reproducible. */
    if (p->i < q->i) return -1;
    if (p->i > q->i) return 1;
    return 0;
}

/* Gini impurity of a histogram, times its total -- the form the weighted sum needs.
 *
 * gini = 1 - sum (c_i/m)^2, so m * gini = m - (sum c_i^2)/m. Returning the weighted value
 * avoids a division per child per candidate and keeps the comparison in one scale. */
static double gini_weighted(const double* c, size_t k, double m) {
    if (m <= 0.0) return 0.0;
    double sq = 0.0;
    for (size_t i = 0; i < k; i++) sq += c[i] * c[i];
    return m - sq / m;
}

/* ---- the growing tree ---- */

typedef struct {
    MlTree* t;
    size_t  cap;
    /* Scratch reused across every node, allocated once. A node-local malloc per recursion
     * would dominate the run time on a deep tree. */
    TreeKey* key;
    double*  cl;        /* left histogram during the sweep  */
    double*  cr;        /* right histogram during the sweep */
    size_t*  idx;       /* the point indices this subtree owns, partitioned in place */
} TreeBuild;

static bool tree_reserve(TreeBuild* b, size_t need) {
    if (need <= b->cap) return true;
    size_t cap = b->cap ? b->cap * 2 : 32;
    while (cap < need) cap *= 2;
    int64_t* nf = realloc(b->t->feature, sizeof(int64_t) * cap);
    if (!nf) return false;
    b->t->feature = nf;
    double* nt = realloc(b->t->thresh, sizeof(double) * cap);
    if (!nt) return false;
    b->t->thresh = nt;
    size_t* nl = realloc(b->t->left, sizeof(size_t) * cap);
    if (!nl) return false;
    b->t->left = nl;
    size_t* nr = realloc(b->t->right, sizeof(size_t) * cap);
    if (!nr) return false;
    b->t->right = nr;
    double* nd = realloc(b->t->dist, sizeof(double) * cap * b->t->k);
    if (!nd) return false;
    b->t->dist = nd;
    b->cap = cap;
    return true;
}

/* Build the subtree owning idx[lo..hi) and return its node index, or SIZE_MAX on failure. */
static size_t tree_build(TreeBuild* b, const double* x, const size_t* y,
                         size_t dim, size_t lo, size_t hi, size_t depth,
                         size_t max_depth, size_t min_split) {
    MlTree* t = b->t;
    size_t k = t->k;
    if (!tree_reserve(b, t->count + 1)) return SIZE_MAX;
    size_t self = t->count++;

    double* hist = t->dist + self * k;
    for (size_t c = 0; c < k; c++) hist[c] = 0.0;
    for (size_t p = lo; p < hi; p++) hist[y[b->idx[p]]] += 1.0;

    size_t m = hi - lo;
    t->feature[self] = -1;          /* a leaf until a split is found */
    t->thresh[self]  = 0.0;
    t->left[self]    = 0;
    t->right[self]   = 0;

    /* Pure node, too small to split, or at the depth limit: stop. A pure node is the common
     * exit, and checking it first avoids the whole sweep for it. */
    size_t nonzero = 0;
    for (size_t c = 0; c < k; c++) if (hist[c] > 0.0) nonzero++;
    if (nonzero <= 1 || m < min_split || depth >= max_depth) return self;

    double parent = gini_weighted(hist, k, (double)m);

    size_t best_f = 0;
    double best_thr = 0.0, best_gain = 0.0;
    size_t best_nleft = 0;
    bool found = false;

    for (size_t f = 0; f < dim; f++) {
        for (size_t p = lo; p < hi; p++) {
            b->key[p - lo].v = x[b->idx[p] * dim + f];
            b->key[p - lo].i = b->idx[p];
        }
        qsort(b->key, m, sizeof(TreeKey), treekey_cmp);

        for (size_t c = 0; c < k; c++) { b->cl[c] = 0.0; b->cr[c] = hist[c]; }

        for (size_t s = 0; s + 1 <= m - 1; s++) {
            size_t pi = b->key[s].i;
            b->cl[y[pi]] += 1.0;
            b->cr[y[pi]] -= 1.0;
            /* Only a boundary between DISTINCT values is a real split; between equal values
             * the two children cannot be separated by any threshold. */
            if (!(b->key[s].v < b->key[s + 1].v)) continue;

            double nl = (double)(s + 1), nr = (double)(m - s - 1);
            double gain = parent - gini_weighted(b->cl, k, nl)
                                 - gini_weighted(b->cr, k, nr);
            /* Strictly greater, so the tie-break is "first found wins" -- and because f
             * ascends outermost and s ascends inside it, that is lower feature index then
             * lower threshold, deterministically. */
            if (gain > best_gain + 1e-12) {
                best_gain  = gain;
                best_f     = f;
                best_thr   = 0.5 * (b->key[s].v + b->key[s + 1].v);
                best_nleft = s + 1;
                found      = true;
            }
        }
    }

    /* No split improves impurity: identical feature rows carrying different classes reach
     * here, and a leaf holding the majority is the only honest answer. Not an error -- the
     * data genuinely does not determine the class. */
    if (!found || best_nleft == 0 || best_nleft == m) return self;

    /* Partition idx[lo..hi) in place around the chosen split. Rebuilt from the feature rather
     * than reusing the sorted key array, because the loop above left `key` sorted by the LAST
     * feature examined, not the best one. */
    size_t w = lo;
    for (size_t p = lo; p < hi; p++)
        if (x[b->idx[p] * dim + best_f] <= best_thr) {
            size_t tmp = b->idx[w]; b->idx[w] = b->idx[p]; b->idx[p] = tmp;
            w++;
        }
    if (w == lo || w == hi) return self;    /* defensive: nothing moved, stay a leaf */

    size_t l = tree_build(b, x, y, dim, lo, w, depth + 1, max_depth, min_split);
    if (l == SIZE_MAX) return SIZE_MAX;
    size_t r = tree_build(b, x, y, dim, w, hi, depth + 1, max_depth, min_split);
    if (r == SIZE_MAX) return SIZE_MAX;

    /* Written AFTER the children, because tree_reserve may have reallocated the arrays while
     * they were being built -- holding a pointer to this node's row across those calls would
     * be a use-after-free. This is why the fields are assigned through t->... here rather
     * than through a cached row pointer. */
    t->feature[self] = (int64_t)best_f;
    t->thresh[self]  = best_thr;
    t->left[self]    = l;
    t->right[self]   = r;
    return self;
}

MlTree* ml_tree_fit(const double* x, const size_t* y, size_t n, size_t dim, size_t k,
                    size_t max_depth, size_t min_split) {
    if (!x || !y || n == 0 || dim == 0 || k == 0) return NULL;
    if (min_split < 2) min_split = 2;

    MlTree* t = calloc(1, sizeof(MlTree));
    if (!t) return NULL;
    t->k = k;
    t->dim = dim;

    TreeBuild b;
    b.t = t; b.cap = 0;
    b.key = malloc(sizeof(TreeKey) * n);
    b.cl  = malloc(sizeof(double) * k);
    b.cr  = malloc(sizeof(double) * k);
    b.idx = malloc(sizeof(size_t) * n);
    if (!b.key || !b.cl || !b.cr || !b.idx) {
        free(b.key); free(b.cl); free(b.cr); free(b.idx);
        ml_tree_free(t);
        return NULL;
    }
    for (size_t i = 0; i < n; i++) b.idx[i] = i;

    size_t root = tree_build(&b, x, y, dim, 0, n, 0, max_depth, min_split);
    free(b.key); free(b.cl); free(b.cr); free(b.idx);
    if (root == SIZE_MAX || t->count == 0) { ml_tree_free(t); return NULL; }
    return t;
}

size_t ml_tree_route(const MlTree* t, const double* x) {
    if (!t || t->count == 0) return 0;
    size_t at = 0;
    /* Bounded by the node count as a backstop: a hand-built tree with a cycle in it would
     * otherwise spin here forever, and this walk runs on every prediction. */
    for (size_t guard = 0; guard <= t->count; guard++) {
        int64_t f = t->feature[at];
        if (f < 0) return at;
        size_t nxt = (x[(size_t)f] <= t->thresh[at]) ? t->left[at] : t->right[at];
        if (nxt >= t->count || nxt == at) return at;
        at = nxt;
    }
    return at;
}

size_t ml_tree_node_class(const MlTree* t, size_t node) {
    if (!t || node >= t->count) return 0;
    const double* h = t->dist + node * t->k;
    size_t best = 0;
    for (size_t c = 1; c < t->k; c++) if (h[c] > h[best]) best = c;
    return best;
}

void ml_tree_free(MlTree* t) {
    if (!t) return;
    free(t->feature); free(t->thresh);
    free(t->left); free(t->right); free(t->dist);
    free(t);
}
