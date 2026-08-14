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

    printf("All ml Predict tests passed.\n");
    return 0;
}
