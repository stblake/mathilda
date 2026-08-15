/* test_findclusters_ndim.c -- FindClusters methods that work on vectors.
 *
 * Eight of the ten were one-dimensional by algorithm, reaching their data only
 * through the sorted projection (d->val indexed by d->order, both NULL for vector
 * input), so each was an algorithm to write rather than a check to relax.
 *
 * ALL TEN are ported. DBSCAN and GaussianMixture replaced their 1-D kernels
 * outright; KMeans, KMedoids, JarvisPatrick and Spectral kept two paths because
 * the general rule moved a covered 1-D answer in each case.
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

/* Separate from check() because the count-taking methods are a different shape,
 * not a variation: FC_ALLOWED denies KMeans and KMedoids the Automatic count, so
 * for them a count is mandatory rather than optional. */
static void check_k(const char* method, const char* data, const char* count,
                    const char* expect) {
    char in[1024];
    snprintf(in, sizeof in, "FindClusters[%s, %s, Method -> \"%s\"]",
             data, count, method);
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

static void test_kmeans_recovers_blobs(void) {
    check_k("KMeans", BLOBS_2D, "3", BLOBS_2D_OUT);
    check_k("KMeans", BLOBS_5D, "3", BLOBS_5D_OUT);
    /* UpTo is the only data-driven count KMeans accepts, so it is worth its own
     * row: asked for more clusters than the data has, it must return the data's
     * count and not the asked-for one. */
    check_k("KMeans", BLOBS_2D, "UpTo[3]", BLOBS_2D_OUT);
    check_k("KMeans", BLOBS_2D, "UpTo[9]", BLOBS_2D_OUT);
    /* ...and a FIXED count must be obeyed exactly, which is the opposite
     * behaviour on the same data. Twelve distinct points, nine asked for. */
    assert_eval_eq("Length[FindClusters[" BLOBS_2D ", 9, Method -> \"KMeans\"]]",
                   "9", 0);
}

static void test_dbscan_recovers_blobs(void) {
    check("DBSCAN", BLOBS_2D, BLOBS_2D_OUT);
    check("DBSCAN", BLOBS_5D, BLOBS_5D_OUT);
    /* An explicit radius must reach the n-D kernel too, not just the default
     * derived from the length scale. */
    assert_eval_eq("FindClusters[" BLOBS_2D ", "
                   "Method -> {\"DBSCAN\", \"NeighborhoodRadius\" -> 2}]",
                   BLOBS_2D_OUT, 0);
}

static void test_dbscan_keeps_noise_as_singletons(void) {
    /* DBSCAN's one structural difference from the other methods: a point in no
     * dense region is noise. Dropping it would lose an input element and stop the
     * result being a partition, so it becomes its own cluster. */
    assert_eval_eq("FindClusters[{{0, 0}, {1, 0}, {0, 1}, {1, 1}, {99, 99}}, "
                   "Method -> \"DBSCAN\"]",
                   "{{{0, 0}, {1, 0}, {0, 1}, {1, 1}}, {{99, 99}}}", 0);
    /* Raising MinPoints makes more points noise, which is the option doing its
     * job: at 3 neither of the first two is core, so neither joins anything and
     * both fall out separately, while the triple survives. */
    assert_eval_eq("FindClusters[{{0, 0}, {1, 0}, {50, 50}, {51, 50}, {52, 50}}, "
                   "Method -> {\"DBSCAN\", \"MinPoints\" -> 3}]",
                   "{{{0, 0}}, {{1, 0}}, {{50, 50}, {51, 50}, {52, 50}}}", 0);
}

/* Three blobs of EIGHT points, not four. Jarvis-Patrick's default NeighborCount is
 * 5, and a 5-NN list cannot fit inside a 4-point blob -- it must reach into a
 * neighbouring one, which then links them. That is the algorithm behaving as
 * specified on data too small for its default, not a defect, and the 1-D kernel
 * does the same thing (a k clamped to n-1 puts every point in one window). So the
 * blob test for this method needs blobs the default k fits inside; asserting the
 * merge instead would be pinning a data-sizing artefact. */
#define JP_BLOBS_2D "Join[Table[{Mod[t, 4], Quotient[t, 4]}, {t, 0, 7}], " \
                    "Table[{40 + Mod[t, 4], Quotient[t, 4]}, {t, 0, 7}], " \
                    "Table[{Mod[t, 4], 40 + Quotient[t, 4]}, {t, 0, 7}]]"
#define JP_BLOBS_5D "Join[" \
    "Table[{Mod[t, 4], Quotient[t, 4], 0, 0, 0}, {t, 0, 7}], " \
    "Table[{40 + Mod[t, 4], Quotient[t, 4], 0, 0, 0}, {t, 0, 7}], " \
    "Table[{Mod[t, 4], Quotient[t, 4], 0, 0, 60}, {t, 0, 7}]]"

static void test_jarvispatrick_recovers_blobs(void) {
    assert_eval_eq("Length[FindClusters[" JP_BLOBS_2D ", "
                   "Method -> \"JarvisPatrick\"]]", "3", 0);
    assert_eval_eq("Length[FindClusters[" JP_BLOBS_5D ", "
                   "Method -> \"JarvisPatrick\"]]", "3", 0);
    /* Membership, not just the count: the first blob must come back whole. */
    assert_eval_eq("First[FindClusters[" JP_BLOBS_2D ", Method -> \"JarvisPatrick\"]]",
                   "{{0, 0}, {1, 0}, {2, 0}, {3, 0}, "
                   "{0, 1}, {1, 1}, {2, 1}, {3, 1}}", 0);
    /* NeighborCount reaches the n-D kernel: on the small 4-point blobs the
     * default merges everything, and lowering k to fit recovers all three. This
     * pair is the option doing visible work. */
    assert_eval_eq("Length[FindClusters[" BLOBS_2D ", Method -> \"JarvisPatrick\"]]",
                   "1", 0);
    assert_eval_eq("FindClusters[" BLOBS_2D ", "
                   "Method -> {\"JarvisPatrick\", \"NeighborCount\" -> 3}]",
                   BLOBS_2D_OUT, 0);
}

static void test_kmedoids_recovers_blobs(void) {
    check_k("KMedoids", BLOBS_2D, "3", BLOBS_2D_OUT);
    check_k("KMedoids", BLOBS_5D, "3", BLOBS_5D_OUT);
    check_k("KMedoids", BLOBS_2D, "UpTo[3]", BLOBS_2D_OUT);
}

static void test_kmedoids_carries_a_tighter_ceiling_than_kmeans(void) {
    /* The two methods share fc_lloyd_ndim and differ in their update step, and that
     * difference is a complexity class: a mean is O(n * dim), a medoid search
     * compares every member against every other member of its cluster and is
     * O(n^2 * dim) however small k is. So KMedoids carries a second, tighter
     * ceiling -- and the cleanest proof that it is real is that the SAME input
     * separates the two methods. */
    const char* sample = "Table[Table[1.0 Mod[7 i + 13 j, 997], {j, 10}], {i, 2000}]";
    char in[512];
    snprintf(in, sizeof in, "Length[FindClusters[%s, 4, Method -> \"KMeans\"]]", sample);
    assert_eval_eq(in, "4", 0);
    snprintf(in, sizeof in, "Head[FindClusters[%s, 4, Method -> \"KMedoids\"]]", sample);
    assert_eval_eq(in, "FindClusters", 0);
}

static void test_kmedoids_ndim_finds_a_better_optimum_than_the_1d_kernel(void) {
    /* KMedoids is split like KMeans -- quantile seeding on a line, farthest-first
     * off it -- and here the two DISAGREE on the same seven numbers. That is not a
     * defect to hide behind an agreement assertion: the n-D answer is strictly
     * better by the method's own objective (total distance from each member to its
     * cluster's medoid), 4 against 16.
     *
     *   1-D kernel:  {{1, 2, 3}, {10}, {11, 12, 25}}   cost 16
     *   n-D kernel:  {{1, 2, 3}, {10, 11, 12}, {25}}   cost 4
     *
     * Both are pinned as they actually behave. Adopting farthest-first on a line
     * would improve the 1-D answer and move a pinned one, so it is a deliberate
     * change of its own -- recorded rather than smuggled in with a port. */
    assert_eval_eq("FindClusters[{1, 2, 3, 10, 11, 12, 25}, 3, Method -> \"KMedoids\"]",
                   "{{1, 2, 3}, {10}, {11, 12, 25}}", 0);
    assert_eval_eq("FindClusters[{{1}, {2}, {3}, {10}, {11}, {12}, {25}}, 3, "
                   "Method -> \"KMedoids\"]",
                   "{{{1}, {2}, {3}}, {{10}, {11}, {12}}, {{25}}}", 0);
}

static void test_spectral_recovers_blobs(void) {
    check("Spectral", BLOBS_2D, BLOBS_2D_OUT);
    check("Spectral", BLOBS_5D, BLOBS_5D_OUT);
}

static void test_spectral_respects_upto_by_merging_nearest_components(void) {
    /* Spectral takes Automatic and UpTo only, never a fixed count. Asked for fewer
     * clusters than the affinity graph has components, it must MERGE -- and which
     * pair is not something the spectrum decides, since on a disconnected graph
     * every partition respecting components has zero cut and all are optimal. The
     * nearest pair is the tie-break, taken from the spanning tree.
     *
     * The three blobs sit at separations 20 and 40, so UpTo[2] must join the two
     * twenty apart and leave the far one alone. A merge of an arbitrary pair would
     * still give two clusters, so the count alone would not catch it -- the
     * membership is the assertion that matters. */
    assert_eval_eq("Length[FindClusters[" BLOBS_2D ", UpTo[1], Method -> \"Spectral\"]]",
                   "1", 0);
    assert_eval_eq("FindClusters[" BLOBS_2D ", UpTo[2], Method -> \"Spectral\"]",
                   "{{{0, 0}, {1, 0}, {0, 1}, {1, 1}, "
                   "{20, 20}, {21, 20}, {20, 21}, {21, 21}}, "
                   "{{0, 40}, {1, 40}, {0, 41}, {1, 41}}}", 0);
    /* And it saturates: asked for more than the data supports, the natural count
     * wins rather than the request. */
    assert_eval_eq("Length[FindClusters[" BLOBS_2D ", UpTo[5], Method -> \"Spectral\"]]",
                   "3", 0);
}

/* GaussianMixture with a FULL covariance needs more points per component than
 * dimensions -- dim*(dim+1)/2 covariance entries plus dim means -- so the nine-point
 * five-dimensional probe the other methods use is under-determined for it: k_max
 * falls to n/(dim+1) == 1 and it correctly returns one cluster. That is BIC and the
 * parameter count doing their job, the same shape of finding as JarvisPatrick's
 * NeighborCount not fitting inside a 4-point blob. Its blob test therefore uses
 * samples a full covariance can afford. */
#define GMM_BLOBS_2D "Join[Table[{Mod[t, 4], Quotient[t, 4]}, {t, 0, 11}], " \
                     "Table[{40 + Mod[t, 4], Quotient[t, 4]}, {t, 0, 11}], " \
                     "Table[{Mod[t, 4], 40 + Quotient[t, 4]}, {t, 0, 11}]]"
#define GMM_BLOBS_5D "Join[" \
    "Table[{Mod[t, 4], Quotient[t, 4], 0, 0, 0}, {t, 0, 23}], " \
    "Table[{40 + Mod[t, 4], Quotient[t, 4], 0, 0, 0}, {t, 0, 23}], " \
    "Table[{Mod[t, 4], Quotient[t, 4], 0, 0, 60}, {t, 0, 23}]]"

static void test_gaussianmixture_recovers_blobs(void) {
    assert_eval_eq("Length[FindClusters[" GMM_BLOBS_2D ", "
                   "Method -> \"GaussianMixture\"]]", "3", 0);
    assert_eval_eq("Length[FindClusters[" GMM_BLOBS_5D ", "
                   "Method -> \"GaussianMixture\"]]", "3", 0);
    /* And the under-determined case is asserted as what it is, so that a later
     * change loosening k_max shows up here rather than silently. */
    assert_eval_eq("Length[FindClusters[" BLOBS_5D ", "
                   "Method -> \"GaussianMixture\"]]", "1", 0);
    /* Two dimensions with twelve points IS affordable, and is recovered exactly. */
    check("GaussianMixture", BLOBS_2D, BLOBS_2D_OUT);
}

static void test_gaussianmixture_survives_a_singular_component(void) {
    /* Identical points give a zero scatter matrix, which is singular: without the
     * variance ridge the Cholesky fails and without the fallback the whole fit dies.
     * One cluster is the right answer and it must arrive without a crash. */
    assert_eval_eq("FindClusters[{{7, 7}, {7, 7}, {7, 7}, {7, 7}}, "
                   "Method -> \"GaussianMixture\"]",
                   "{{{7, 7}, {7, 7}, {7, 7}, {7, 7}}}", 0);
    /* Collinear points in the plane are rank-deficient too -- a real covariance with
     * a zero eigenvalue -- and must be modelled, not refused. */
    assert_eval_eq("Head[FindClusters[{{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}}, "
                   "Method -> \"GaussianMixture\"]]", "List", 0);
}

static void test_kmeans_is_independent_of_input_order(void) {
    /* A k-means whose answer depends on the order its points were typed in is
     * seeding from index 0. This is why the n-D initialisation starts from the
     * point nearest the centroid rather than from the first point: the centroid is
     * a property of the SET. Compared as sorted-sorted so that cluster numbering,
     * which legitimately follows input order, is not what is being asserted. */
    const char* srt = "Sort[Map[Sort, FindClusters[%s, 3, Method -> \"KMeans\"]]]";
    char a[1024], b[1024], c[1024];
    snprintf(a, sizeof a, srt, BLOBS_2D);
    snprintf(b, sizeof b, srt,
             "{{20,21},{0,1},{21,21},{1,1},{0,40},{20,20},{1,40},{0,0},"
             "{1,41},{1,0},{21,20},{0,41}}");
    snprintf(c, sizeof c, srt, "Reverse[" BLOBS_2D "]");
    char eq[3072];
    snprintf(eq, sizeof eq, "(%s === %s) && (%s === %s)", a, b, a, c);
    assert_eval_eq(eq, "True", 0);
}

static void test_kmeans_automatic_is_refused_on_both_surfaces(void) {
    /* Not a porting gap. FC_ALLOWED denies KMeans FC_COUNT_AUTOMATIC, so a bare
     * call declines in one dimension too -- asserted here for both so that the
     * n-D decline is never mistaken for an unported path and "fixed". */
    assert_eval_eq("FindClusters[{1, 2, 3, 10, 11, 12, 25}, Method -> \"KMeans\"]",
                   "FindClusters[{1, 2, 3, 10, 11, 12, 25}, Method -> \"KMeans\"]", 0);
    assert_eval_eq("FindClusters[{{1, 1}, {2, 2}, {20, 20}}, Method -> \"KMeans\"]",
                   "FindClusters[{{1, 1}, {2, 2}, {20, 20}}, Method -> \"KMeans\"]", 0);
}

static void test_kmeans_declines_only_on_the_work_product(void) {
    /* The cap is on n * k * dim, not on n -- Lloyd is LINEAR in n and cheaper than
     * the spanning tree already built for the same input, so capping n would refuse
     * work just paid for. Both rows use the same 5000 x 10 machine sample and
     * differ only in k, which isolates the claim: a small k is admitted at a size
     * an n-cap would have refused, and a k in the thousands -- 100 near-quadratic
     * passes over the whole sample -- declines.
     *
     * The bound is deliberately conservative in one direction: k == n converges in
     * a single iteration, since farthest-first hands every point its own centre, so
     * it is refused by a budget that assumes the iteration count. Refusing a
     * degenerate request is cheaper than modelling it. */
    const char* sample = "Table[Table[1.0 Mod[7 i + 13 j, 997], {j, 10}], {i, 5000}]";
    char in[512];
    snprintf(in, sizeof in, "Length[FindClusters[%s, 8, Method -> \"KMeans\"]]", sample);
    assert_eval_eq(in, "8", 0);
    snprintf(in, sizeof in, "Head[FindClusters[%s, 1000, Method -> \"KMeans\"]]", sample);
    assert_eval_eq(in, "FindClusters", 0);
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

    /* KMeans keeps two implementations on purpose -- 1-D seeds at quantiles of the
     * sorted distinct values, n-D seeds farthest-first -- so this is the one place
     * the two are made to answer the same question, and the only guard that they
     * have not diverged. Unifying them would move the pinned 1-D answers, so it is
     * a deliberate behaviour change rather than a refactor; this row is what would
     * catch it being done by accident. */
    assert_eval_eq("FindClusters[{1, 2, 3, 10, 11, 12, 25}, 3, Method -> \"KMeans\"]",
                   "{{1, 2, 3}, {10, 11, 12}, {25}}", 0);
    assert_eval_eq("FindClusters[{{1}, {2}, {3}, {10}, {11}, {12}, {25}}, 3, "
                   "Method -> \"KMeans\"]",
                   "{{{1}, {2}, {3}}, {{10}, {11}, {12}}, {{25}}}", 0);

    /* DBSCAN went the other way from KMeans: one kernel now serves both
     * dimensionalities, because all 22 one-dimensional pins pass through the
     * general rule unchanged. So this pair is not two implementations agreeing but
     * one implementation reached by two representations -- still worth asserting,
     * since FC_KIND_SCALAR and a dim-1 FC_KIND_POINT read different fields. */
    assert_eval_eq("FindClusters[{1, 2, 3, 100}, Method -> \"DBSCAN\"]",
                   "{{1, 2, 3}, {100}}", 0);
    assert_eval_eq("FindClusters[{{1}, {2}, {3}, {100}}, Method -> \"DBSCAN\"]",
                   "{{{1}, {2}, {3}}, {{100}}}", 0);

    /* JarvisPatrick is split like KMeans, for a reason found rather than assumed:
     * its general rule leaves all 22 pins passing yet moves a list_tests answer at
     * NeighborCount -> 2, because the 1-D kernel counts shared neighbours as the
     * overlap of two contiguous windows and links only adjacent sorted pairs. Here
     * the two forms are asked the same question on data where they agree. */
    assert_eval_eq("FindClusters[{1, 2, 3, 100}, Method -> \"JarvisPatrick\"]",
                   "{{1, 2, 3, 100}}", 0);
    assert_eval_eq("FindClusters[{{1}, {2}, {3}, {100}}, Method -> \"JarvisPatrick\"]",
                   "{{{1}, {2}, {3}, {100}}}", 0);

    /* Spectral is split too, and the reason is the exact mirror of the KMedoids
     * story. Its n-D kernel thresholds embedding jumps against the MEAN because a
     * good embedding collapses within-cluster distance and drives the MEDIAN to
     * zero; on a line the 1-D kernel thresholds DATA gaps against the median, where
     * the mean is the one that misleads. Each statistic is right in its own domain,
     * so routing scalars through the n-D kernel under-counts -- it turned
     * {{1,2,3},{10,11,12},{25}} into {{1,2,3,10,11,12},{25}}. Here the two forms
     * agree, which is what makes the row worth keeping. */
    assert_eval_eq("FindClusters[{1, 2, 3, 100}, Method -> \"Spectral\"]",
                   "{{1, 2, 3}, {100}}", 0);
    assert_eval_eq("FindClusters[{{1}, {2}, {3}, {100}}, Method -> \"Spectral\"]",
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
    /* KMeans is the method most able to break this on its own: asked for more
     * clusters than there are distinct points it would have to split a duplicate
     * pair to comply, so the count must yield instead. */
    check_k("KMeans", dup2, "2", dup2_out);
    assert_eval_eq("FindClusters[{{5, 5}, {5, 5}, {5, 5}, {40, 40}, {41, 41}}, 4, "
                   "Method -> \"KMeans\"]",
                   "{{{5, 5}, {5, 5}, {5, 5}}, {{40, 40}}, {{41, 41}}}", 0);
    /* KMedoids has the same exposure as KMeans: asked for more clusters than there
     * are distinct points it would have to split a duplicate pair to comply. */
    assert_eval_eq("FindClusters[{{5, 5}, {5, 5}, {5, 5}, {40, 40}, {41, 41}}, 4, "
                   "Method -> \"KMedoids\"]",
                   "{{{5, 5}, {5, 5}, {5, 5}}, {{40, 40}}, {{41, 41}}}", 0);
    check("DBSCAN", dup2, dup2_out);
    /* JarvisPatrick needs blobs its default k fits inside (see JP_BLOBS_2D), so the
     * duplicate goes there: a repeated {0, 0} prepended to the first blob must come
     * back inside that blob's cluster, adjacent to its twin, and must not become a
     * fourth cluster of its own. */
    assert_eval_eq("Length[FindClusters[Join[{{0, 0}}, " JP_BLOBS_2D "], "
                   "Method -> \"JarvisPatrick\"]]", "3", 0);
    assert_eval_eq("First[FindClusters[Join[{{0, 0}}, " JP_BLOBS_2D "], "
                   "Method -> \"JarvisPatrick\"]]",
                   "{{0, 0}, {0, 0}, {1, 0}, {2, 0}, {3, 0}, "
                   "{0, 1}, {1, 1}, {2, 1}, {3, 1}}", 0);
    assert_eval_eq("FindClusters[{{5, 5}, {5, 5}, {5, 5}, {40, 40}, {41, 41}}, "
                   "Method -> \"DBSCAN\"]",
                   "{{{5, 5}, {5, 5}, {5, 5}}, {{40, 40}, {41, 41}}}", 0);

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
    /* DBSCAN too, and for the same reason: an eps-neighbourhood needs coordinates.
     * This is now the ONLY thing its single kernel refuses. */
    assert_eval_eq("FindClusters[{\"aa\", \"ab\", \"zz\"}, Method -> \"DBSCAN\"]",
                   "FindClusters[{\"aa\", \"ab\", \"zz\"}, Method -> \"DBSCAN\"]", 0);
    assert_eval_eq("FindClusters[{\"aa\", \"ab\", \"zz\"}, Method -> \"JarvisPatrick\"]",
                   "FindClusters[{\"aa\", \"ab\", \"zz\"}, "
                   "Method -> \"JarvisPatrick\"]", 0);
    /* Still works for the gap methods, which read only the tree. */
    assert_eval_eq("FindClusters[{\"aa\", \"ab\", \"zz\"}, 2]",
                   "{{\"aa\", \"ab\"}, {\"zz\"}}", 0);
}

static void test_every_method_now_answers_in_ndim(void) {
    /* The guard names the ported set explicitly, so the rest must still decline
     * above one dimension -- reading d->val there would dereference NULL, so a
     * premature relaxation is a crash rather than a wrong answer. These rows come
     * off the list as each method is ported, which makes the progress visible. */
    /* EMPTY at last: all ten methods cluster vectors. What replaces this list is the
     * opposite assertion -- that every method now answers, and that the only thing
     * still declined above one dimension is SEQUENCE input, which has no
     * coordinates. A method silently regressing to a decline would show up here. */
    const char* all[] = { "Agglomerate", "SpanningTree", "MeanShift",
                          "NeighborhoodContraction", "DBSCAN", "JarvisPatrick",
                          "Spectral", "GaussianMixture" };
    for (size_t i = 0; i < sizeof(all) / sizeof(all[0]); i++) {
        char in[256];
        snprintf(in, sizeof in,
                 "Head[FindClusters[{{1, 1}, {9, 9}}, Method -> \"%s\"]]", all[i]);
        assert_eval_eq(in, "List", 0);
    }
    /* KMeans and KMedoids need a count, so they are asked in the form they take. */
    assert_eval_eq("Head[FindClusters[{{1, 1}, {9, 9}}, 2, Method -> \"KMeans\"]]",
                   "List", 0);
    assert_eval_eq("Head[FindClusters[{{1, 1}, {9, 9}}, 2, Method -> \"KMedoids\"]]",
                   "List", 0);
    /* Spectral is ported; two points are below the n < 3 floor and come back as one
     * cluster. That floor must be handled on the point path itself -- falling
     * through to the 1-D branch reaches fc_scatter, which indexes d->order, NULL
     * here, and segfaulted before this row existed. */
    assert_eval_eq("FindClusters[{{1, 1}, {9, 9}}, Method -> \"Spectral\"]",
                   "{{{1, 1}, {9, 9}}}", 0);
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

    /* For KMeans the metric is asserted where the two answers genuinely DIFFER,
     * which is the only form of this test that proves anything. Three points: the
     * apex is nearer the origin under Euclidean (10 vs 11) and nearer the far point
     * under Manhattan (14 vs 11), so the pair that clusters flips. */
    assert_eval_eq("FindClusters[{{0, 0}, {0, 11}, {8, 6}}, 2, Method -> \"KMeans\", "
                   "DistanceFunction -> EuclideanDistance]",
                   "{{{0, 0}}, {{0, 11}, {8, 6}}}", 0);
    assert_eval_eq("FindClusters[{{0, 0}, {0, 11}, {8, 6}}, 2, Method -> \"KMeans\", "
                   "DistanceFunction -> ManhattanDistance]",
                   "{{{0, 0}, {0, 11}}, {{8, 6}}}", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_shift_methods_recover_blobs_2d);
    TEST(test_shift_methods_recover_blobs_5d);
    TEST(test_kmeans_recovers_blobs);
    TEST(test_dbscan_recovers_blobs);
    TEST(test_jarvispatrick_recovers_blobs);
    TEST(test_kmedoids_recovers_blobs);
    TEST(test_spectral_recovers_blobs);
    TEST(test_gaussianmixture_recovers_blobs);
    TEST(test_gaussianmixture_survives_a_singular_component);
    TEST(test_spectral_respects_upto_by_merging_nearest_components);
    TEST(test_kmedoids_carries_a_tighter_ceiling_than_kmeans);
    TEST(test_kmedoids_ndim_finds_a_better_optimum_than_the_1d_kernel);
    TEST(test_dbscan_keeps_noise_as_singletons);
    TEST(test_kmeans_is_independent_of_input_order);
    TEST(test_kmeans_automatic_is_refused_on_both_surfaces);
    TEST(test_kmeans_declines_only_on_the_work_product);
    TEST(test_scalar_and_dim1_point_agree);
    TEST(test_equal_points_are_never_split);
    TEST(test_strings_still_decline);
    TEST(test_every_method_now_answers_in_ndim);
    TEST(test_metric_applies_to_ported_methods);

    printf("All FindClusters n-dimensional tests passed.\n");
    return 0;
}
