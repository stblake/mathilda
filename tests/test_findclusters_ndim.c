/* test_findclusters_ndim.c -- FindClusters methods that work on vectors.
 *
 * Eight of the ten were one-dimensional by algorithm, reaching their data only
 * through the sorted projection (d->val indexed by d->order, both NULL for vector
 * input). They are being ported one at a time; this file grows with each, and the
 * guard in the builtin names the ported set explicitly so a new method cannot
 * default into being considered done.
 *
 * Ported so far: MeanShift, NeighborhoodContraction.
 *
 * The property that matters most here is not "n-D works" but "n-D and 1-D are the
 * same code". A dim-1 POINT is not a SCALAR -- its elements are Lists -- yet the
 * two must agree exactly, because they are the same points written differently.
 * Every ported method is asserted both ways.
 */
#include <stdio.h>
#include <string.h>
#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

/* Three well-separated square blobs in the plane, and three in five dimensions.
 * "Well separated" is doing real work: the gaps between blobs are an order of
 * magnitude larger than the spacing within one, so any method that has a notion
 * of scale at all must recover exactly these three groups. A method that returns
 * anything else on this input is broken rather than merely opinionated. */
#define BLOBS_2D "{{0,0},{1,0},{0,1},{1,1}, {20,20},{21,20},{20,21},{21,21}, " \
                 "{0,40},{1,40},{0,41},{1,41}}"
#define BLOBS_2D_OUT "{{{0, 0}, {1, 0}, {0, 1}, {1, 1}}, " \
                     "{{20, 20}, {21, 20}, {20, 21}, {21, 21}}, " \
                     "{{0, 40}, {1, 40}, {0, 41}, {1, 41}}}"

#define BLOBS_5D "{{0,0,0,0,0},{1,0,0,0,1},{0,1,0,1,0}, " \
                 "{30,30,30,30,30},{31,30,30,30,31},{30,31,30,31,30}, " \
                 "{0,0,60,0,0},{1,0,60,0,1},{0,1,60,1,0}}"
#define BLOBS_5D_OUT "{{{0, 0, 0, 0, 0}, {1, 0, 0, 0, 1}, {0, 1, 0, 1, 0}}, " \
                     "{{30, 30, 30, 30, 30}, {31, 30, 30, 30, 31}, {30, 31, 30, 31, 30}}, " \
                     "{{0, 0, 60, 0, 0}, {1, 0, 60, 0, 1}, {0, 1, 60, 1, 0}}}"

static void check(const char* method, const char* data, const char* expect) {
    char in[1024];
    snprintf(in, sizeof in, "FindClusters[%s, Method -> \"%s\"]", data, method);
    assert_eval_eq(in, expect, 0);
}

static void test_shift_methods_recover_blobs_2d(void) {
    check("MeanShift", BLOBS_2D, BLOBS_2D_OUT);
    check("NeighborhoodContraction", BLOBS_2D, BLOBS_2D_OUT);
}

static void test_shift_methods_recover_blobs_5d(void) {
    /* Five dimensions rather than three, so that nothing can pass by accidentally
     * treating the input as a plane. The third coordinate is what separates the
     * last blob, which also catches an implementation that only reads the first
     * two components. */
    check("MeanShift", BLOBS_5D, BLOBS_5D_OUT);
    check("NeighborhoodContraction", BLOBS_5D, BLOBS_5D_OUT);
}

static void test_scalar_and_dim1_point_agree(void) {
    /* The same numbers written as scalars and as one-component vectors. These take
     * different branches -- FC_KIND_SCALAR reads d->val, FC_KIND_POINT reads
     * d->coord -- and must agree, since they are the same points.
     *
     * The first pair is degenerate on purpose: with three points the median of two
     * spanning-tree edge weights averages 1 and 98 into a scale of 49.5, which
     * cannot resolve them, so both forms return one cluster. That is the method's
     * behaviour and not a porting artefact, which is exactly why it is asserted
     * for both surfaces rather than only for the one that looks right. */
    assert_eval_eq("FindClusters[{1, 2, 100}, Method -> \"MeanShift\"]",
                   "{{1, 2, 100}}", 0);
    assert_eval_eq("FindClusters[{{1}, {2}, {100}}, Method -> \"MeanShift\"]",
                   "{{{1}, {2}, {100}}}", 0);
    /* A fourth point gives the median something to work with, and both forms
     * resolve the split. */
    assert_eval_eq("FindClusters[{1, 2, 3, 100}, Method -> \"MeanShift\"]",
                   "{{1, 2, 3}, {100}}", 0);
    assert_eval_eq("FindClusters[{{1}, {2}, {3}, {100}}, Method -> \"MeanShift\"]",
                   "{{{1}, {2}, {3}}, {{100}}}", 0);
    assert_eval_eq("FindClusters[{1, 2, 3, 100}, Method -> \"NeighborhoodContraction\"]",
                   "{{1, 2, 3}, {100}}", 0);
    assert_eval_eq("FindClusters[{{1}, {2}, {3}, {100}}, "
                   "Method -> \"NeighborhoodContraction\"]",
                   "{{{1}, {2}, {3}}, {{100}}}", 0);
}

static void test_equal_points_are_never_split(void) {
    /* Enforced globally by the fold over zero-weight tree edges rather than by each
     * kernel. The fold is kind-agnostic, so it should cover the ports for free --
     * this asserts that it actually does, in the plane where three of the 1-D
     * kernels historically broke the same invariant.
     *
     * A duplicated point inside a blob, with enough points that the median edge
     * weight is nonzero. An earlier version of this test used four points with
     * three coincident, which is degenerate for a different reason worth recording:
     * a zero median makes the length scale fall back to the mean edge weight, i.e.
     * range / (n - 1), and the flat kernel's default radius is three times the
     * scale -- so at n == 4 the radius is EXACTLY the range, the outlier lands
     * precisely on the inclusion boundary, and NeighborhoodContraction merges
     * everything however far away it is put. That is pre-existing 1-D behaviour,
     * identical in both dimensionalities (checked), so it is recorded as a quirk
     * rather than asserted as though it were intended. */
    const char* dup2 = "{{0,0},{0,0},{1,0},{0,1},{1,1}, "
                       "{20,20},{21,20},{20,21},{21,21}}";
    const char* dup2_out = "{{{0, 0}, {0, 0}, {1, 0}, {0, 1}, {1, 1}}, "
                           "{{20, 20}, {21, 20}, {20, 21}, {21, 21}}}";
    check("MeanShift", dup2, dup2_out);
    check("NeighborhoodContraction", dup2, dup2_out);

    /* And the 1-D form of the same shape, to show the fold behaves the same way on
     * both surfaces. */
    assert_eval_eq("FindClusters[{0, 0, 1, 2, 3, 20, 21, 22, 23}, "
                   "Method -> \"MeanShift\"]",
                   "{{0, 0, 1, 2, 3}, {20, 21, 22, 23}}", 0);
    assert_eval_eq("FindClusters[{0, 0, 1, 2, 3, 20, 21, 22, 23}, "
                   "Method -> \"NeighborhoodContraction\"]",
                   "{{0, 0, 1, 2, 3}, {20, 21, 22, 23}}", 0);
}

static void test_strings_still_decline(void) {
    /* SEQUENCE input has no coordinates, so there is nothing for a shift method to
     * move toward higher density. Declined rather than approximated -- the gap
     * methods still handle strings, via edit distances in the tree. */
    assert_eval_eq("FindClusters[{\"aa\", \"ab\", \"zz\"}, Method -> \"MeanShift\"]",
                   "FindClusters[{\"aa\", \"ab\", \"zz\"}, Method -> \"MeanShift\"]", 0);
    assert_eval_eq("FindClusters[{\"aa\", \"ab\", \"zz\"}, "
                   "Method -> \"NeighborhoodContraction\"]",
                   "FindClusters[{\"aa\", \"ab\", \"zz\"}, "
                   "Method -> \"NeighborhoodContraction\"]", 0);
    /* Still works for the gap methods, which read only the tree. */
    assert_eval_eq("FindClusters[{\"aa\", \"ab\", \"zz\"}, 2]",
                   "{{\"aa\", \"ab\"}, {\"zz\"}}", 0);
}

static void test_unported_methods_still_decline_in_ndim(void) {
    /* The guard names the ported set explicitly, so the rest must still decline
     * above one dimension -- reading d->val there would dereference NULL, so a
     * premature relaxation is a crash rather than a wrong answer. These rows come
     * off the list as each method is ported, which makes the progress visible. */
    const char* unported[] = { "KMeans", "KMedoids", "DBSCAN",
                               "GaussianMixture", "JarvisPatrick" };
    for (size_t i = 0; i < sizeof(unported) / sizeof(unported[0]); i++) {
        char in[256], out[256];
        snprintf(in,  sizeof in,  "FindClusters[{{1, 1}, {9, 9}}, Method -> \"%s\"]",
                 unported[i]);
        snprintf(out, sizeof out, "FindClusters[{{1, 1}, {9, 9}}, Method -> \"%s\"]",
                 unported[i]);
        assert_eval_eq(in, out, 0);
    }
    /* Spectral takes no fixed count, so its n-D decline is checked in the form it
     * does accept. */
    assert_eval_eq("FindClusters[{{1, 1}, {9, 9}}, Method -> \"Spectral\"]",
                   "FindClusters[{{1, 1}, {9, 9}}, Method -> \"Spectral\"]", 0);
}

static void test_metric_applies_to_ported_methods(void) {
    /* The ported methods go through the same metric selection as the tree, so a
     * Manhattan bandwidth is a Manhattan bandwidth. Asserted on the blobs, where
     * both metrics agree on the answer -- the point is that neither declines nor
     * changes the grouping, not that they differ. */
    assert_eval_eq("FindClusters[" BLOBS_2D ", Method -> \"MeanShift\", "
                   "DistanceFunction -> ManhattanDistance]", BLOBS_2D_OUT, 0);
    assert_eval_eq("FindClusters[" BLOBS_2D ", Method -> \"MeanShift\", "
                   "DistanceFunction -> EuclideanDistance]", BLOBS_2D_OUT, 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_shift_methods_recover_blobs_2d);
    TEST(test_shift_methods_recover_blobs_5d);
    TEST(test_scalar_and_dim1_point_agree);
    TEST(test_equal_points_are_never_split);
    TEST(test_strings_still_decline);
    TEST(test_unported_methods_still_decline_in_ndim);
    TEST(test_metric_applies_to_ported_methods);

    printf("All FindClusters n-dimensional tests passed.\n");
    return 0;
}
