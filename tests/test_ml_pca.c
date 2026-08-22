/* test_ml_pca.c -- Standardize and PrincipalComponents.
 *
 * These are the first NEW builtins of the machine-learning work rather than ports, so
 * unlike the FindClusters suites there is no prior behaviour to pin against. The
 * assertions are therefore chosen to be things that must hold for any correct
 * implementation, not things this one happens to produce:
 *
 *   - Standardize must agree with (x - Mean[x]) / StandardDeviation[x] written out by
 *     hand. That is the definition, and it also forces the n-1 divisor to match
 *     StandardDeviation's -- a divisor mismatch is invisible on a mean but not here.
 *   - A principal-component rotation is orthogonal, so it must preserve TOTAL
 *     variance while concentrating it in the leading component. Both halves matter: a
 *     projection that lost variance would still pass a "PC1 is largest" check.
 *   - Collinear data has rank 1, so the second component must be exactly zero.
 *   - Correlation and Covariance must differ where the columns have different scales,
 *     which is the only case in which the option means anything.
 */
#include <stdio.h>
#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

static void test_standardize_matches_its_definition(void) {
    /* Not an independently computed constant: the same quantity written two ways.
     * If the divisor were n instead of n - 1 this row would fail while a
     * mean-is-zero check would still pass. */
    assert_eval_eq("Standardize[{1., 2., 3., 4.}] == "
                   "((# - Mean[#])/StandardDeviation[#] &)[{1., 2., 3., 4.}]",
                   "True", 0);
    assert_eval_eq("Chop[Mean[Standardize[{1., 2., 3., 4., 9.}]]]", "0", 0);
    assert_eval_eq("StandardDeviation[Standardize[{1., 2., 3., 4., 9.}]]", "1.0", 0);
}

static void test_standardize_is_column_wise(void) {
    /* Rows are observations, columns are variables, so the two columns here -- which
     * differ by a factor of ten -- must standardise to the SAME values. A row-wise
     * implementation would give something else entirely. */
    assert_eval_eq("Standardize[{{1., 10.}, {2., 20.}, {3., 30.}}]",
                   "{{-1.0, -1.0}, {0.0, 0.0}, {1.0, 1.0}}", 0);
}

static void test_standardize_leaves_a_constant_column_at_zero(void) {
    /* Zero variance carries no information, so the honest standardised value is "no
     * deviation from the mean". Dividing by the zero standard deviation would give
     * Indeterminate and poison every reduction over the row. */
    assert_eval_eq("Standardize[{{1., 5.}, {2., 5.}, {3., 5.}}]",
                   "{{-1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}}", 0);
}

static void test_result_is_a_plain_list_not_a_visible_ndarray(void) {
    /* The machine bridge na_build_matrix returns a VISIBLE NDArray, whose head is
     * NDArray, so a result built that way compares False against the literal list a
     * user would write -- while Inverse, Dot and LinearSolve all compare True. This
     * row is what keeps the new builtins on the same surface as the old ones. */
    assert_eval_eq("Head[Standardize[{{1., 10.}, {2., 20.}, {3., 30.}}]]", "List", 0);
    assert_eval_eq("Standardize[{{1., 10.}, {2., 20.}, {3., 30.}}] === "
                   "{{-1.0, -1.0}, {0.0, 0.0}, {1.0, 1.0}}", "True", 0);
    assert_eval_eq("Head[PrincipalComponents[{{1., 2.}, {3., 5.}, {4., 4.}}]]",
                   "List", 0);
}

static void test_pca_rotation_preserves_total_variance(void) {
    /* An orthogonal rotation moves variance between components without creating or
     * destroying any, so the column variances must sum to the same total before and
     * after. This is the assertion that would catch a projection using the wrong
     * matrix, or eigenvectors that were not orthonormal. */
    assert_eval_eq("Module[{d = {{1., 2.}, {3., 5.}, {4., 4.}, {6., 9.}, {7., 8.}}, p},"
                   " p = PrincipalComponents[d];"
                   " Chop[(Variance[Map[First, p]] + Variance[Map[Last, p]]) - "
                   "      (Variance[Map[First, d]] + Variance[Map[Last, d]])]]",
                   "0", 0);
}

static void test_pca_orders_components_by_decreasing_variance(void) {
    assert_eval_eq("Module[{p = PrincipalComponents["
                   "{{1., 2.}, {3., 5.}, {4., 4.}, {6., 9.}, {7., 8.}}]},"
                   " Variance[Map[First, p]] > Variance[Map[Last, p]]]", "True", 0);
}

static void test_pca_puts_rank_one_data_entirely_in_one_component(void) {
    /* Points on a line have one non-zero eigenvalue, so every second coordinate must
     * be zero -- not merely small. Chop at 10^-10 rather than comparing to 0.0
     * because the projection is a floating-point contraction, but the tolerance is
     * far below anything a wrong answer would produce. */
    assert_eval_eq("Chop[Map[Last, PrincipalComponents["
                   "{{0., 0.}, {1., 1.}, {2., 2.}, {3., 3.}, {4., 4.}}]], 10.^-10]",
                   "{0, 0, 0, 0, 0}", 0);
}

static void test_pca_method_correlation_differs_from_covariance(void) {
    /* The columns here differ by two orders of magnitude, so covariance PCA is
     * dominated by the second and correlation PCA is not. Asserted as a difference
     * rather than against fixed numbers: the point of the option is that it changes
     * the answer, and a test that passed for both would prove nothing. */
    assert_eval_eq("PrincipalComponents[{{1., 100.}, {2., 300.}, {3., 200.}}] === "
                   "PrincipalComponents[{{1., 100.}, {2., 300.}, {3., 200.}}, "
                   "Method -> \"Correlation\"]", "False", 0);
    /* Correlation standardises to unit variance, so the total variance of the
     * rotated data is the number of variables. */
    assert_eval_eq("Module[{p = PrincipalComponents["
                   "{{1., 100.}, {2., 300.}, {3., 200.}, {4., 500.}, {5., 400.}}, "
                   "Method -> \"Correlation\"]},"
                   " Chop[(Variance[Map[First, p]] + Variance[Map[Last, p]]) - 2.0]]",
                   "0", 0);
}

static void test_bad_arguments_decline(void) {
    /* An unknown Method declines rather than silently picking one, so a typo is
     * visible instead of quietly changing the statistics. */
    assert_eval_eq("Head[PrincipalComponents[{{1., 2.}, {3., 4.}}, "
                   "Method -> \"Nope\"]]", "PrincipalComponents", 0);
    /* A flat list is n observations of ONE variable, which has no components to
     * rotate. */
    assert_eval_eq("Head[PrincipalComponents[{1., 2., 3.}]]",
                   "PrincipalComponents", 0);
    assert_eval_eq("Head[Standardize[x]]", "Standardize", 0);
    assert_eval_eq("Head[PrincipalComponents[{{1., a}, {2., 3.}}]]",
                   "PrincipalComponents", 0);
}

#define DR_DATA "{{1., 2., 3.}, {3., 5., 4.}, {4., 4., 8.}, {6., 9., 2.}, {7., 8., 9.}}"

static void test_dimension_reduce_is_a_prefix_of_the_full_rotation(void) {
    /* Reducing to k dimensions must give exactly the first k principal components, not
     * a separately-fitted k-component model. This is the assertion that would catch a
     * reducer that re-derived its axes from a truncated matrix. */
    assert_eval_eq("Chop[Map[Take[#, 2] &, PrincipalComponents[" DR_DATA "]] - "
                   "DimensionReduce[" DR_DATA ", 2], 10.^-10]",
                   "{{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}", 0);
    assert_eval_eq("Length[First[DimensionReduce[" DR_DATA ", 1]]]", "1", 0);
    assert_eval_eq("Head[DimensionReduce[" DR_DATA ", 2]]", "List", 0);
}

static void test_classical_mds_on_euclidean_distances_reproduces_pca(void) {
    /* The strongest check available here, and the reason MDS was worth implementing in
     * the same iteration: classical (Torgerson) scaling of EUCLIDEAN distances is
     * mathematically the same embedding as PCA, reached by a completely different
     * route -- an n x n double-centred distance matrix rather than a dim x dim
     * covariance. Agreement is therefore evidence about BOTH, in a way that no
     * self-consistency check could be.
     *
     * Compared on absolute values because an eigenvector's sign is arbitrary and the
     * two matrices have different shapes, so the canonicalisation in
     * ml_sym_eigen_desc does not force the same choice for both. */
    assert_eval_eq("Chop[Map[Abs, DimensionReduce[" DR_DATA ", 2, "
                   "Method -> \"MultidimensionalScaling\"]] - "
                   "Map[Abs, Map[Take[#, 2] &, PrincipalComponents[" DR_DATA "]]], "
                   "10.^-8]",
                   "{{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}", 0);
}

static void test_lsa_differs_from_pca_because_it_does_not_centre(void) {
    /* The single `if` that separates them is the whole mathematical difference between
     * principal components and a truncated SVD, so the two must not agree on data with
     * a non-zero mean -- if they did, the centring step would be dead code. */
    assert_eval_eq("DimensionReduce[" DR_DATA ", 2, "
                   "Method -> \"LatentSemanticAnalysis\"] === "
                   "DimensionReduce[" DR_DATA ", 2]", "False", 0);
}

static void test_dimension_reduce_declines_rather_than_padding(void) {
    /* Asking for more dimensions than the data supports has no answer, and padding
     * with zeros would look like a successful reduction to a caller checking only the
     * shape. */
    assert_eval_eq("Head[DimensionReduce[" DR_DATA ", 5]]", "DimensionReduce", 0);
    assert_eval_eq("Head[DimensionReduce[" DR_DATA ", 0]]", "DimensionReduce", 0);
    assert_eval_eq("Head[DimensionReduce[" DR_DATA ", 2, Method -> \"tSNE\"]]",
                   "DimensionReduce", 0);
    /* No target dimension: Wolfram picks one, this declines rather than guessing. */
    assert_eval_eq("Head[DimensionReduce[" DR_DATA "]]", "DimensionReduce", 0);
    assert_eval_eq("Head[DimensionReduce[{1., 2., 3.}, 1]]", "DimensionReduce", 0);
}

/* Training data with a deliberately OFF-CENTRE mean. This matters: the reducer must
 * centre new rows on the TRAINING means, and centring on the incoming batch's own mean
 * instead gives all zeros for a single point -- which is indistinguishable from correct
 * on data that happens to sit near the origin. Column means here are {104, 210}. */
#define RED_TRAIN "{{100., 200.}, {102., 205.}, {104., 210.}, {106., 215.}, {108., 220.}}"

static void test_reducer_reproduces_dimension_reduce_on_its_training_data(void) {
    /* The two builtins must agree where their domains overlap: DimensionReduce returns
     * the reduced training data, and the reducer applied to that same training data must
     * give the same thing. Disagreement would mean one of them centres or projects
     * differently. */
    assert_eval_eq("Chop[DimensionReduce[" RED_TRAIN ", 1] - "
                   "DimensionReduction[" RED_TRAIN ", 1][" RED_TRAIN "], 10.^-9]",
                   "{{0}, {0}, {0}, {0}, {0}}", 0);
}

static void test_reducer_centres_new_data_on_the_training_mean(void) {
    /* The discriminating pair, and the reason the training data is off-centre.
     *
     * A point AT the training mean must project to zero -- that identifies WHICH mean is
     * being subtracted. And a point away from it must project far from zero: if the
     * implementation centred each incoming batch on its own mean, a single point would
     * always land at zero, so this second row is what fails on the plausible bug. */
    assert_eval_eq("Chop[DimensionReduction[" RED_TRAIN ", 1][{104., 210.}], 10.^-9]",
                   "{0}", 0);
    assert_eval_eq("Chop[Abs[First[DimensionReduction[" RED_TRAIN ", 1][{110., 225.}]]] "
                   "> 10.]", "True", 0);
    /* Symmetry as a further constraint: points equally far either side of the training
     * mean must project to values of equal magnitude and opposite sign. */
    assert_eval_eq("Chop[First[DimensionReduction[" RED_TRAIN ", 1][{110., 225.}]] + "
                   "First[DimensionReduction[" RED_TRAIN ", 1][{98., 195.}]]]", "0", 0);
}

static void test_reducer_accepts_a_point_or_a_batch(void) {
    assert_eval_eq("DimensionReduction[" RED_TRAIN ", 1][{110., 225.}] === "
                   "First[DimensionReduction[" RED_TRAIN ", 1][{{110., 225.}}]]",
                   "True", 0);
    assert_eval_eq("Head[DimensionReduction[" RED_TRAIN ", 1][{{110., 225.}}]]",
                   "List", 0);
    assert_eval_eq("Length[DimensionReduction[" RED_TRAIN ", 2][{110., 225.}]]", "2", 0);
}

static void test_reducer_properties_and_declines(void) {
    assert_eval_eq("DimensionReduction[" RED_TRAIN ", 1][\"ReducedDimension\"]", "1", 0);
    assert_eval_eq("DimensionReduction[" RED_TRAIN ", 1][\"FeatureCount\"]", "2", 0);
    assert_eval_eq("Head[DimensionReduction[" RED_TRAIN ", 3]]", "DimensionReduction", 0);
    assert_eval_eq("Head[DimensionReduction[" RED_TRAIN ", 0]]", "DimensionReduction", 0);
    assert_eval_eq("Head[DimensionReduction[{1., 2., 3.}, 1]]", "DimensionReduction", 0);
    /* A wrongly-shaped input must not be projected. Asserted as "no list of the reduced
     * width came back" rather than by Head, since Head of an unevaluated composite
     * application is the whole reducer, not a symbol. */
    assert_eval_eq("MatrixQ[DimensionReduction[" RED_TRAIN ", 1][{1., 2., 3.}]] || "
                   "VectorQ[DimensionReduction[" RED_TRAIN ", 1][{1., 2., 3.}], NumberQ]",
                   "False", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_standardize_matches_its_definition);
    TEST(test_standardize_is_column_wise);
    TEST(test_standardize_leaves_a_constant_column_at_zero);
    TEST(test_result_is_a_plain_list_not_a_visible_ndarray);
    TEST(test_pca_rotation_preserves_total_variance);
    TEST(test_pca_orders_components_by_decreasing_variance);
    TEST(test_pca_puts_rank_one_data_entirely_in_one_component);
    TEST(test_pca_method_correlation_differs_from_covariance);
    TEST(test_bad_arguments_decline);
    TEST(test_dimension_reduce_is_a_prefix_of_the_full_rotation);
    TEST(test_classical_mds_on_euclidean_distances_reproduces_pca);
    TEST(test_lsa_differs_from_pca_because_it_does_not_centre);
    TEST(test_dimension_reduce_declines_rather_than_padding);
    TEST(test_reducer_reproduces_dimension_reduce_on_its_training_data);
    TEST(test_reducer_centres_new_data_on_the_training_mean);
    TEST(test_reducer_accepts_a_point_or_a_batch);
    TEST(test_reducer_properties_and_declines);

    printf("All ml PCA tests passed.\n");
    return 0;
}
