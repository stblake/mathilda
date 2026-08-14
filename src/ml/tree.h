/* tree.h -- CART classification trees, at buffer level.
 *
 * The first thing in src/ml that reuses NOTHING from the other five families: no distance,
 * no density, no linear algebra, no label-index arithmetic beyond counting. A tree needs an
 * impurity criterion, a recursive splitter, a stopping rule and a node representation, none
 * of which any earlier family had a use for. That is exactly why it was recorded as deferred
 * for as long as it was, rather than being wedged into an iteration.
 *
 * Buffer level and in its own module for the same reason gmm.c was: there will be three
 * consumers -- Classify's DecisionTree method, LearnDistribution's, and RandomForest, which
 * is this same fit run many times over bootstrap resamples. Writing it inside classify.c
 * would guarantee the forest either duplicated it or reached across into it.
 *
 * THE NODE REPRESENTATION, and why the class histogram is stored at EVERY node rather than
 * only at leaves. Two flat arrays of the same length: the split (feature, threshold, left,
 * right) and the class counts. A leaf is `feature < 0`. Storing counts everywhere costs
 * nodes*k doubles and buys three things -- the predicted class is the arg-max of the node's
 * own histogram with no separate leaf table to keep in step, "Probabilities" comes out of
 * the same array for free, and an internal node's histogram is meaningful in its own right
 * (it is the subtree's distribution), so nothing about the layout is dead weight.
 */
#ifndef ML_TREE_H
#define ML_TREE_H

#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>

typedef struct {
    size_t   count;     /* nodes in use */
    size_t   k;         /* classes */
    size_t   dim;       /* features */
    int64_t* feature;   /* split feature, or -1 for a leaf */
    double*  thresh;    /* go left when x[feature] <= thresh */
    size_t*  left;
    size_t*  right;
    double*  dist;      /* count * k class COUNTS, one row per node */
} MlTree;

/* Grow a tree. `y` holds class indices in [0, k).
 *
 * `max_depth` bounds the tree; `min_split` is the fewest points a node may hold and still be
 * split. With min_split = 2 and a generous depth the tree grows until every leaf is pure or
 * unsplittable, which is what makes "reproduces every training label" an exact property to
 * assert rather than an accuracy figure to hope for.
 *
 * Returns NULL on allocation failure or degenerate input (n == 0, k == 0). */
MlTree* ml_tree_fit(const double* x, const size_t* y, size_t n, size_t dim, size_t k,
                    size_t max_depth, size_t min_split);

/* Index of the leaf a point routes to. */
size_t ml_tree_route(const MlTree* t, const double* x);

/* Arg-max of a node's histogram, ties to the lowest class index. */
size_t ml_tree_node_class(const MlTree* t, size_t node);

void ml_tree_free(MlTree* t);

#endif /* ML_TREE_H */
