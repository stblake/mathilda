/* test_ml_predict.c -- Predict, LinearModelFit, and the fitted-model object.
 *
 * The strongest assertion available here is AGREEMENT WITH AN INDEPENDENT
 * IMPLEMENTATION: `Fit` already solves least squares by a different route, so a
 * matching answer is evidence about both. That is the same kind of check that
 * validated MDS against PCA, and it is worth more than any number chosen by hand.
 *
 * Beyond that, the properties asserted are ones any correct implementation must have:
 * exact recovery of an exactly-linear relationship, survival of storage and
 * re-application (the whole point of a trained-model object), and refusal on a
 * singular system rather than returning one of infinitely many answers.
 */
#include <stdio.h>
#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

/* y = 2x + 1 exactly, so the coefficients are not a fitting compromise and any
 * correct least squares must return them to machine precision. */
#define LINE "{1. -> 3., 2. -> 5., 3. -> 7., 4. -> 9.}"
/* y = 1 + 2a + 3b exactly, matrix form with the response last. */
#define PLANE "{{1., 1., 6.}, {2., 1., 8.}, {1., 2., 9.}, {3., 2., 13.}, {2., 3., 14.}}"

static void test_recovers_an_exact_relationship(void) {
    assert_eval_eq("Predict[" LINE "][\"Coefficients\"]", "{1.0, 2.0}", 0);
    assert_eval_eq("Predict[" PLANE "][\"Coefficients\"]", "{1.0, 2.0, 3.0}", 0);
}

static void test_agrees_with_the_existing_fit(void) {
    /* Fit is an independent least-squares implementation already in the tree, so this
     * is a cross-check rather than a self-consistency check: the two share no code
     * path. Compared through the symbolic form Fit returns, evaluated at a point, so
     * the comparison does not depend on how either presents its coefficients. */
    assert_eval_eq("Chop[(Fit[{{1., 3.}, {2., 5.}, {3., 7.}, {4., 9.}}, {1, x}, x] "
                   "/. x -> 7.) - Predict[" LINE "][7.]]", "0", 0);
    assert_eval_eq("Chop[(Fit[{{1., 3.}, {2., 5.}, {3., 7.}, {4., 9.}}, {1, x}, x] "
                   "/. x -> -2.5) - Predict[" LINE "][-2.5]]", "0", 0);
}

static void test_the_model_survives_storage_and_reapplication(void) {
    /* The reason a fitted model is an object at all. If it were re-fitted on each use,
     * or lost its parameters when assigned, this would fail. */
    assert_eval_eq("Module[{p = Predict[" LINE "], q}, q = p; q[5.]]", "11.0", 0);
    assert_eval_eq("Module[{p = Predict[" LINE "]}, "
                   "{p[10.], p[10.], p[{10.}]}]", "{21.0, 21.0, 21.0}", 0);
}

static void test_a_one_feature_model_accepts_a_bare_scalar(void) {
    /* p[3.] for a single-variable regression is the natural thing to write, so
     * requiring p[{3.}] would be pedantry. Both must give the same answer. */
    assert_eval_eq("Predict[" LINE "][3.] == Predict[" LINE "][{3.}]", "True", 0);
}

static void test_properties_are_readable(void) {
    assert_eval_eq("Predict[" LINE "][\"Method\"]", "\"LinearRegression\"", 0);
    assert_eval_eq("Predict[" PLANE "][\"FeatureCount\"]", "2", 0);
}

static void test_linear_model_fit_agrees_with_predict(void) {
    assert_eval_eq("LinearModelFit[" PLANE "][\"Coefficients\"] === "
                   "Predict[" PLANE "][\"Coefficients\"]", "True", 0);
}

static void test_rules_and_matrix_forms_agree(void) {
    /* The same data written both ways must fit the same model -- which form was used is
     * deliberately not remembered. */
    assert_eval_eq("Predict[{{1., 1.}, {2., 2.}, {3., 3.}, {4., 4.}}][\"Coefficients\"] "
                   "=== Predict[{1. -> 1., 2. -> 2., 3. -> 3., 4. -> 4.}][\"Coefficients\"]",
                   "True", 0);
}

static void test_singular_and_malformed_input_declines(void) {
    /* Perfectly collinear features: the normal equations are singular and there is no
     * unique fit. A pseudo-inverse would return one of infinitely many answers and look
     * like a successful fit, so this declines instead. */
    assert_eval_eq("Head[Predict[{{1., 2., 3.}, {2., 4., 6.}, "
                   "{3., 6., 9.}, {4., 8., 12.}}]]", "Predict", 0);
    /* Fewer observations than parameters is the other way to be underdetermined. */
    assert_eval_eq("Head[Predict[{{1., 1., 2.}}]]", "Predict", 0);
    /* An unimplemented method declines rather than silently linear-regressing. */
    assert_eval_eq("Head[Predict[" LINE ", Method -> \"RandomForest\"]]", "Predict", 0);
    assert_eval_eq("Head[Predict[{}]]", "Predict", 0);
    assert_eval_eq("Head[Predict[x]]", "Predict", 0);
}

static void test_wrong_shaped_application_stays_unevaluated(void) {
    /* A two-feature model given one feature must not guess. The application stays
     * unevaluated, and the right way to assert that is that no NUMBER came back --
     * Head of an unevaluated application is the whole composite head
     * PredictorFunction[...], not the symbol, so a Head test would be asserting the
     * wrong thing (and did, on the first attempt). */
    assert_eval_eq("NumberQ[Predict[" PLANE "][{1.}]]", "False", 0);
    assert_eval_eq("NumberQ[Predict[" PLANE "][{1., 2.}]]", "True", 0);
    assert_eval_eq("NumberQ[Predict[" LINE "][\"NoSuchProperty\"]]", "False", 0);
}

#define KNN_DATA "{1. -> 10., 2. -> 20., 3. -> 30., 10. -> 100., 11. -> 110.}"
#define KNN1 "Method -> {\"NearestNeighbors\", \"NeighborsNumber\" -> 1}"

static void test_knn_with_k_one_agrees_with_the_existing_nearest(void) {
    /* The independent cross-check for this method. `Nearest` is a separate
     * implementation already in the tree, so a 1-nearest predictor must return the
     * label of whatever Nearest picks. Queried on both sides of a midpoint (2.4 and
     * 2.6) so the comparison actually exercises the tie-breaking region rather than
     * points that are obviously nearest to one row.
     *
     * The limit of the check is worth stating: Nearest here supports only the
     * one-neighbour scalar form -- its k form and its point form decline -- so this
     * validates neighbour SELECTION at k = 1 and nothing about the averaging. */
    const char* qs[] = { "0.5", "2.4", "2.6", "7.", "9.9", "12." };
    for (size_t i = 0; i < sizeof(qs) / sizeof(qs[0]); i++) {
        char in[512];
        snprintf(in, sizeof in,
                 "Chop[Predict[" KNN_DATA ", " KNN1 "][%s] - "
                 "10. First[Nearest[{1., 2., 3., 10., 11.}, %s]]]", qs[i], qs[i]);
        assert_eval_eq(in, "0", 0);
    }
}

static void test_knn_averages_its_neighbours(void) {
    /* k = 3 around 2. takes 10, 20, 30 -- mean 20. k = 2 takes 10 and 20 -- mean 15.
     * Different k must give different answers, or the option is decoration. */
    assert_eval_eq("Predict[" KNN_DATA ", " KNN1 "][2.]", "20.0", 0);
    assert_eval_eq("Predict[" KNN_DATA ", Method -> {\"NearestNeighbors\", "
                   "\"NeighborsNumber\" -> 2}][2.]", "15.0", 0);
    assert_eval_eq("Predict[" KNN_DATA ", Method -> \"NearestNeighbors\"][2.]", "20.0", 0);
    /* Asking for every point makes the prediction the global mean, everywhere. */
    assert_eval_eq("Predict[" KNN_DATA ", Method -> {\"NearestNeighbors\", "
                   "\"NeighborsNumber\" -> 5}][2.] == Mean[{10., 20., 30., 100., 110.}]",
                   "True", 0);
}

static void test_knn_is_not_linear(void) {
    /* The two methods must genuinely differ -- and this test needs CURVED data to show
     * it. KNN_DATA above is exactly y = 10x, so a neighbour average at an interior
     * point equals the linear prediction there and the first version of this test
     * passed for the wrong reason, proving nothing. y = x^2 has curvature a line cannot
     * follow, so the two methods must disagree. */
    const char* quad = "{1. -> 1., 2. -> 4., 3. -> 9., 4. -> 16., 5. -> 25.}";
    char in[512];
    snprintf(in, sizeof in,
             "Predict[%s, Method -> \"NearestNeighbors\"][3.] == Predict[%s][3.]",
             quad, quad);
    assert_eval_eq(in, "False", 0);
    /* And the k-NN answer is the neighbour mean, which on x = 3 with k = 3 is the mean
     * of the labels at 2, 3 and 4: (4 + 9 + 16)/3. */
    snprintf(in, sizeof in,
             "Chop[Predict[%s, Method -> \"NearestNeighbors\"][3.] - (4. + 9. + 16.)/3.]",
             quad);
    assert_eval_eq(in, "0", 0);
}

static void test_knn_properties_and_multifeature(void) {
    assert_eval_eq("Predict[" KNN_DATA ", Method -> \"NearestNeighbors\"]"
                   "[\"NeighborCount\"]", "3", 0);
    assert_eval_eq("Predict[" KNN_DATA ", Method -> \"NearestNeighbors\"][\"Method\"]",
                   "\"NearestNeighbors\"", 0);
    /* Two features, matrix form with the response last: the single nearest row to
     * {10., 10.} is its own, so k = 1 returns that row's label exactly. */
    assert_eval_eq("Predict[{{0., 0., 1.}, {1., 0., 2.}, {0., 1., 3.}, {10., 10., 99.}}, "
                   KNN1 "][{10., 10.}]", "99.0", 0);
}

static void test_knn_option_errors_decline(void) {
    /* A neighbour count is meaningless for a regression, so it is refused rather than
     * accepted and ignored -- silently ignoring it would hide a real mistake. */
    assert_eval_eq("Head[Predict[" KNN_DATA ", Method -> {\"LinearRegression\", "
                   "\"NeighborsNumber\" -> 2}]]", "Predict", 0);
    assert_eval_eq("Head[Predict[" KNN_DATA ", Method -> {\"NearestNeighbors\", "
                   "\"Nope\" -> 2}]]", "Predict", 0);
    assert_eval_eq("Head[Predict[" KNN_DATA ", Method -> {\"NearestNeighbors\", "
                   "\"NeighborsNumber\" -> 0}]]", "Predict", 0);
}

static void test_fitted_models_print_abbreviated(void) {
    /* Note the quoting: this harness renders the method string WITH quotes, where the
     * REPL's Print does not -- the abbreviation itself is identical either way.
     *
     * A NearestNeighbors predictor carries its ENTIRE training set as its parameter
     * block, so the unabridged form is unreadable at any real size and scrolls the
     * useful answer off screen. Models print their method and elide the rest, following
     * InterpolatingFunction, which already does exactly this in print.c. */
    assert_eval_eq("Predict[" LINE "]", "PredictorFunction[\"LinearRegression\", <>]", 0);
    assert_eval_eq("Predict[" PLANE ", Method -> \"NearestNeighbors\"]",
                   "PredictorFunction[\"NearestNeighbors\", <>]", 0);
    assert_eval_eq("DimensionReduction[{{100., 200.}, {102., 205.}, {104., 210.}}, 1]",
                   "DimensionReducerFunction[\"PrincipalComponentsAnalysis\", <>]", 0);
}

static void test_abbreviation_hides_nothing_and_breaks_nothing(void) {
    /* Eliding is only safe because the information is one keystroke away rather than
     * gone -- FullForm must still reveal the parameters. */
    assert_eval_eq("FullForm[Predict[" LINE "]]",
                   "PredictorFunction[\"LinearRegression\", List[1.0, 2.0], 1, 0]", 0);
    /* And the abbreviation is a PRINTING change only: application and property access
     * must be untouched. A change that reached into the object rather than its rendering
     * would show up here. */
    assert_eval_eq("Predict[" LINE "][10.]", "21.0", 0);
    assert_eval_eq("Predict[" LINE "][\"Coefficients\"]", "{1.0, 2.0}", 0);
    assert_eval_eq("Length[Predict[" PLANE ", Method -> \"NearestNeighbors\"]"
                   "[\"TrainingData\"]]", "5", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_recovers_an_exact_relationship);
    TEST(test_agrees_with_the_existing_fit);
    TEST(test_the_model_survives_storage_and_reapplication);
    TEST(test_a_one_feature_model_accepts_a_bare_scalar);
    TEST(test_properties_are_readable);
    TEST(test_linear_model_fit_agrees_with_predict);
    TEST(test_rules_and_matrix_forms_agree);
    TEST(test_singular_and_malformed_input_declines);
    TEST(test_wrong_shaped_application_stays_unevaluated);
    TEST(test_knn_with_k_one_agrees_with_the_existing_nearest);
    TEST(test_knn_averages_its_neighbours);
    TEST(test_knn_is_not_linear);
    TEST(test_knn_properties_and_multifeature);
    TEST(test_knn_option_errors_decline);
    TEST(test_fitted_models_print_abbreviated);
    TEST(test_abbreviation_hides_nothing_and_breaks_nothing);

    printf("All ml Predict tests passed.\n");
    return 0;
}
