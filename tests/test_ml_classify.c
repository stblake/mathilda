/* test_ml_classify.c -- Classify, ClassifierFunction, and the label vocabulary.
 *
 * Two of these assertions are ABSOLUTE rather than comparative, which is what makes them
 * worth more than an accuracy figure:
 *
 *   - a k = 1 nearest-neighbour classifier must reproduce EVERY training label exactly,
 *     because each training point is its own nearest neighbour. Not "high accuracy" --
 *     exactly right, every row;
 *   - class probabilities must sum to 1. A classifier whose probabilities were each
 *     individually plausible but did not sum to 1 would pass any per-class check.
 *
 * The rest establish that a class really is an arbitrary expression -- strings, symbols
 * and integers all tested -- which is the point of the vocabulary this family added.
 */
#include <stdio.h>
#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

#define TR "{{0., 0.} -> \"red\", {1., 0.} -> \"red\", {0., 1.} -> \"red\", " \
           "{10., 10.} -> \"blue\", {11., 10.} -> \"blue\", {10., 11.} -> \"blue\"}"

static void test_k1_reproduces_every_training_label(void) {
    /* The exact, absolute check. Each training point is its own nearest neighbour, so a
     * k = 1 classifier that gets any row wrong has a bug in its distance, its vocabulary
     * or its vote -- there is no "close enough" here. */
    assert_eval_eq("Module[{tr = " TR ", c},"
                   " c = Classify[tr];"
                   " And @@ Table[c[First[tr[[i]]]] === Last[tr[[i]]], {i, Length[tr]}]]",
                   "True", 0);
    /* And it generalises the obvious way on either side of the gap. */
    assert_eval_eq("Classify[" TR "][{0.5, 0.5}]", "\"red\"", 0);
    assert_eval_eq("Classify[" TR "][{10.5, 10.5}]", "\"blue\"", 0);
}

static void test_class_probabilities_sum_to_one(void) {
    /* The second absolute check. Asserted at k = 5 as well, where the vote is genuinely
     * split (0.6/0.4) rather than degenerate -- a 1.0/0.0 split would sum to 1 even if the
     * normalisation were wrong. */
    assert_eval_eq("Chop[Total[Map[Last, Classify[" TR "][{0.5, 0.5}, \"Probabilities\"]]] - 1.]",
                   "0", 0);
    assert_eval_eq("Chop[Total[Map[Last, Classify[" TR ", "
                   "Method -> {\"NearestNeighbors\", \"NeighborsNumber\" -> 5}]"
                   "[{5., 5.}, \"Probabilities\"]]] - 1.]", "0", 0);
    /* The split at k = 5 must be non-degenerate, or the row above proves little. */
    assert_eval_eq("Module[{p = Classify[" TR ", Method -> {\"NearestNeighbors\", "
                   "\"NeighborsNumber\" -> 5}][{5., 5.}, \"Probabilities\"]},"
                   " 0. < Last[First[p]] < 1.]", "True", 0);
}

static void test_a_class_can_be_any_expression(void) {
    /* The whole point of the label vocabulary: earlier families took numeric responses,
     * and a class is not a number. Strings, symbols and integers must all work, and a
     * symbol must not be confused with the string of the same name. */
    assert_eval_eq("Classify[{1. -> a, 2. -> a, 9. -> b, 10. -> b}][\"Classes\"]",
                   "{a, b}", 0);
    assert_eval_eq("Classify[{1. -> a, 2. -> a, 9. -> b, 10. -> b}][1.5]", "a", 0);
    assert_eval_eq("Classify[{1. -> 0, 2. -> 0, 9. -> 1, 10. -> 1}][\"Classes\"]",
                   "{0, 1}", 0);
    assert_eval_eq("Classify[{1. -> 0, 2. -> 0, 9. -> 1, 10. -> 1}][9.5]", "1", 0);
    /* Comparison is structural, so "a" and a are two classes rather than one. */
    assert_eval_eq("Length[Classify[{1. -> \"a\", 2. -> a, 9. -> \"b\"}][\"Classes\"]]",
                   "3", 0);
}

static void test_class_order_is_first_appearance_and_deterministic(void) {
    /* Not sorted -- first appearance in the training data. What matters is that it is
     * deterministic and readable next to the data; reversing the data reverses the
     * vocabulary, which demonstrates the rule rather than an accident. */
    assert_eval_eq("Classify[" TR "][\"Classes\"]", "{\"red\", \"blue\"}", 0);
    assert_eval_eq("Classify[Reverse[" TR "]][\"Classes\"]", "{\"blue\", \"red\"}", 0);
    /* The prediction must not depend on that ordering. */
    assert_eval_eq("Classify[" TR "][{0.5, 0.5}] === "
                   "Classify[Reverse[" TR "]][{0.5, 0.5}]", "True", 0);
}

static void test_three_classes_and_properties(void) {
    assert_eval_eq("Classify[{{0.,0.} -> \"a\", {5.,5.} -> \"b\", "
                   "{10.,10.} -> \"c\"}][{4.9, 4.9}]", "\"b\"", 0);
    assert_eval_eq("Classify[" TR "][\"FeatureCount\"]", "2", 0);
    assert_eval_eq("Classify[" TR "][\"NeighborCount\"]", "1", 0);
    assert_eval_eq("Classify[" TR "][\"Method\"]", "\"NearestNeighbors\"", 0);
    /* A fitted classifier elides, like every other fitted model. */
    assert_eval_eq("Classify[" TR "]", "ClassifierFunction[\"NearestNeighbors\", <>]", 0);
}

static void test_malformed_input_declines(void) {
    /* Only the rule form is accepted, and that is not arbitrary: a matrix with the class
     * in its last column cannot work, because a class need not be a number. */
    assert_eval_eq("Head[Classify[{{1., 2.}, {3., 4.}}]]", "Classify", 0);
    assert_eval_eq("Head[Classify[{}]]", "Classify", 0);
    assert_eval_eq("Head[Classify[x]]", "Classify", 0);
    assert_eval_eq("Head[Classify[" TR ", Method -> \"RandomForest\"]]", "Classify", 0);
    assert_eval_eq("Head[Classify[" TR ", Method -> {\"NearestNeighbors\", "
                   "\"NeighborsNumber\" -> 0}]]", "Classify", 0);
    /* A wrongly-shaped query, and an unknown property, leave the application
     * unevaluated -- so no class comes back. */
    assert_eval_eq("StringQ[Classify[" TR "][{1.}]]", "False", 0);
    assert_eval_eq("StringQ[Classify[" TR "][\"NoSuchProperty\"]]", "False", 0);
}

#define NB "{{0.,0.} -> \"red\", {1.,0.} -> \"red\", {0.,1.} -> \"red\", " \
           "{1.,1.} -> \"red\", {10.,10.} -> \"blue\", {11.,10.} -> \"blue\", " \
           "{10.,11.} -> \"blue\", {11.,11.} -> \"blue\"}"
#define NB1 "{1. -> \"a\", 2. -> \"a\", 3. -> \"a\", " \
            "7. -> \"b\", 8. -> \"b\", 9. -> \"b\"}"

static void test_naivebayes_separates_well_separated_blobs(void) {
    /* Absolute, not a percentage: every training row, correct. */
    assert_eval_eq("Module[{tr = " NB ", b},"
                   " b = Classify[tr, Method -> \"NaiveBayes\"];"
                   " And @@ Table[b[First[tr[[i]]]] === Last[tr[[i]]], {i, Length[tr]}]]",
                   "True", 0);
    assert_eval_eq("Classify[" NB ", Method -> \"NaiveBayes\"][{0.5, 0.5}]", "\"red\"", 0);
    assert_eval_eq("Classify[" NB ", Method -> \"NaiveBayes\"][{10.5, 10.5}]",
                   "\"blue\"", 0);
}

static void test_naivebayes_posteriors_sum_to_one(void) {
    /* Asserted at the MIDPOINT as well, where the split is an exact 0.5/0.5 -- a
     * degenerate 1.0/0.0 split would sum to 1 even with a broken softmax. */
    assert_eval_eq("Chop[Total[Map[Last, Classify[" NB ", Method -> \"NaiveBayes\"]"
                   "[{5.5, 5.5}, \"Probabilities\"]]] - 1.]", "0", 0);
    assert_eval_eq("Module[{p = Classify[" NB ", Method -> \"NaiveBayes\"]"
                   "[{5.5, 5.5}, \"Probabilities\"]}, 0.4 < Last[First[p]] < 0.6]",
                   "True", 0);
    assert_eval_eq("Chop[Total[Map[Last, Classify[" NB ", Method -> \"NaiveBayes\"]"
                   "[{0.5, 0.5}, \"Probabilities\"]]] - 1.]", "0", 0);
}

static void test_naivebayes_matches_the_closed_form_gaussian_comparison(void) {
    /* The independent-implementation check, and the strongest available for this method.
     * With one feature and equal priors the decision is just "which prior-weighted
     * Gaussian density is larger", which PDF[NormalDistribution[...]] computes by a
     * completely separate path. Note the fit uses the ML (n) divisor per class, so the
     * comparison variance is Variance[...] * (n-1)/n -- getting that wrong would move the
     * boundary slightly and the straddling points below would catch it.
     *
     * Queried at 4.9 and 5.1, either side of the boundary at 5.0, so the test actually
     * exercises the decision rather than two obvious regions. */
    const char* xs[] = { "0.", "1.5", "3.5", "4.9", "5.1", "6.5", "8.5", "12." };
    for (size_t i = 0; i < sizeof(xs) / sizeof(xs[0]); i++) {
        char in[900];
        snprintf(in, sizeof in,
                 "Module[{nb, ma, mb, va, vb, la, lb},"
                 " nb = Classify[" NB1 ", Method -> \"NaiveBayes\"];"
                 " ma = Mean[{1.,2.,3.}]; mb = Mean[{7.,8.,9.}];"
                 " va = Variance[{1.,2.,3.}] 2/3; vb = Variance[{7.,8.,9.}] 2/3;"
                 " la = 0.5 PDF[NormalDistribution[ma, Sqrt[va]], %s];"
                 " lb = 0.5 PDF[NormalDistribution[mb, Sqrt[vb]], %s];"
                 " If[la > lb, \"a\", \"b\"] === nb[%s]]", xs[i], xs[i], xs[i]);
        assert_eval_eq(in, "True", 0);
    }
}

static void test_naivebayes_variance_floor_keeps_a_constant_feature_finite(void) {
    /* A class whose feature is constant has zero variance there and therefore INFINITE
     * density at that value -- it would win every comparison involving that feature. The
     * floor is a fraction of the feature's overall variance, so it is scale-invariant. The
     * observable consequence: probabilities stay finite and the informative feature still
     * decides. */
    const char* cz = "{{5., 1.} -> \"x\", {5., 2.} -> \"x\", {5., 3.} -> \"x\", "
                     "{9., 1.} -> \"y\", {9., 2.} -> \"y\", {9., 3.} -> \"y\"}";
    char in[700];
    snprintf(in, sizeof in, "Classify[%s, Method -> \"NaiveBayes\"][{5., 2.}]", cz);
    assert_eval_eq(in, "\"x\"", 0);
    snprintf(in, sizeof in, "Classify[%s, Method -> \"NaiveBayes\"][{9., 2.}]", cz);
    assert_eval_eq(in, "\"y\"", 0);
    snprintf(in, sizeof in,
             "Chop[Total[Map[Last, Classify[%s, Method -> \"NaiveBayes\"]"
             "[{7., 2.}, \"Probabilities\"]]] - 1.]", cz);
    assert_eval_eq(in, "0", 0);
}

static void test_naivebayes_option_errors_and_coexistence(void) {
    /* A neighbour count means nothing to a Bayes classifier, so it is refused rather than
     * ignored; and NeighborCount is not a property of one. */
    assert_eval_eq("Head[Classify[" NB ", Method -> {\"NaiveBayes\", "
                   "\"NeighborsNumber\" -> 3}]]", "Classify", 0);
    assert_eval_eq("StringQ[Classify[" NB ", Method -> \"NaiveBayes\"][\"NeighborCount\"]]",
                   "False", 0);
    /* Both methods on one head must coexist -- adding a method must not disturb the other. */
    assert_eval_eq("Classify[" NB "][{0.5, 0.5}]", "\"red\"", 0);
    assert_eval_eq("Classify[" NB ", Method -> \"NaiveBayes\"][\"Method\"]",
                   "\"NaiveBayes\"", 0);
    assert_eval_eq("Classify[" NB ", Method -> \"NaiveBayes\"]",
                   "ClassifierFunction[\"NaiveBayes\", <>]", 0);
}

#define LR1 "{1. -> \"a\", 2. -> \"a\", 3. -> \"a\", " \
            "7. -> \"b\", 8. -> \"b\", 9. -> \"b\"}"

static void test_logistic_separates_separable_data(void) {
    assert_eval_eq("Module[{tr = " LR1 ", L},"
                   " L = Classify[tr, Method -> \"LogisticRegression\"];"
                   " And @@ Table[L[First[tr[[i]]]] === Last[tr[[i]]], {i, Length[tr]}]]",
                   "True", 0);
    assert_eval_eq("Module[{t2 = {{0.,0.} -> \"n\", {1.,0.} -> \"n\", {0.,1.} -> \"n\","
                   " {9.,9.} -> \"p\", {10.,9.} -> \"p\", {9.,10.} -> \"p\"}, L},"
                   " L = Classify[t2, Method -> \"LogisticRegression\"];"
                   " And @@ Table[L[First[t2[[i]]]] === Last[t2[[i]]], {i, Length[t2]}]]",
                   "True", 0);
}

static void test_logistic_probability_is_exactly_half_on_the_boundary(void) {
    /* The EXACT identity, and much stronger than an accuracy figure. The fitted boundary is
     * where intercept + coef x = 0, i.e. x = -intercept/coef; the probability there must be
     * exactly 0.5 because that is what the logistic of zero is. Any error in how the
     * coefficients are stored, read back, or combined at application time shows up here,
     * because the identity ties the FIT and the APPLICATION together. */
    assert_eval_eq("Module[{L, b, x0},"
                   " L = Classify[" LR1 ", Method -> \"LogisticRegression\"];"
                   " b = Part[L, 2, 2];"                   /* the coefficient row */
                   " x0 = -First[b]/Last[b];"
                   " Abs[Last[Last[L[x0, \"Probabilities\"]]] - 0.5] < 1.*^-12]",
                   "True", 0);
    /* And that boundary sits in the gap between the classes, not outside the data. */
    assert_eval_eq("Module[{L, b, x0},"
                   " L = Classify[" LR1 ", Method -> \"LogisticRegression\"];"
                   " b = Part[L, 2, 2]; x0 = -First[b]/Last[b];"
                   " 3. < x0 < 7.]", "True", 0);
}

static void test_logistic_probability_is_monotone_and_normalised(void) {
    /* Monotonicity along the coefficient direction is what makes it a logistic model rather
     * than an arbitrary function fitted to the labels. */
    assert_eval_eq("Module[{L, ps},"
                   " L = Classify[" LR1 ", Method -> \"LogisticRegression\"];"
                   " ps = Table[Last[Last[L[x, \"Probabilities\"]]],"
                   " {x, {0., 2., 4., 5., 6., 8., 10.}}];"
                   " And @@ Table[ps[[i]] <= ps[[i + 1]], {i, Length[ps] - 1}]]", "True", 0);
    assert_eval_eq("Chop[Total[Map[Last, Classify[" LR1 ", "
                   "Method -> \"LogisticRegression\"][5., \"Probabilities\"]]] - 1.]",
                   "0", 0);
}

static void test_logistic_coefficients_stay_finite_on_separable_data(void) {
    /* THE reason the ridge exists. On linearly separable data the unpenalised likelihood is
     * UNBOUNDED -- driving the coefficients to infinity drives every fitted probability to 0
     * or 1 -- so plain Newton never converges. The ridge makes the penalised objective
     * strictly concave, so the fit is finite and unique. This row is what fails if the ridge
     * is ever removed as "unnecessary". */
    assert_eval_eq("Module[{b = Part[Classify[" LR1 ", "
                   "Method -> \"LogisticRegression\"], 2, 2]},"
                   " And @@ Map[(NumberQ[#] && Abs[#] < 1000.) &, b]]", "True", 0);
}

static void test_logistic_declines_and_coexists(void) {
    /* More than two classes would need one-vs-rest or a softmax; rather than guess, it
     * declines -- a silently binarised three-class problem would just be wrong. */
    assert_eval_eq("Head[Classify[{1. -> \"a\", 5. -> \"b\", 9. -> \"c\"}, "
                   "Method -> \"LogisticRegression\"]]", "Classify", 0);
    assert_eval_eq("Head[Classify[" LR1 ", Method -> {\"LogisticRegression\", "
                   "\"NeighborsNumber\" -> 3}]]", "Classify", 0);
    /* All three methods coexist on the one head. */
    assert_eval_eq("Classify[" LR1 "][2.5]", "\"a\"", 0);
    assert_eval_eq("Classify[" LR1 ", Method -> \"NaiveBayes\"][2.5]", "\"a\"", 0);
    assert_eval_eq("Classify[" LR1 ", Method -> \"LogisticRegression\"]",
                   "ClassifierFunction[\"LogisticRegression\", <>]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_k1_reproduces_every_training_label);
    TEST(test_class_probabilities_sum_to_one);
    TEST(test_a_class_can_be_any_expression);
    TEST(test_class_order_is_first_appearance_and_deterministic);
    TEST(test_three_classes_and_properties);
    TEST(test_malformed_input_declines);
    TEST(test_naivebayes_separates_well_separated_blobs);
    TEST(test_naivebayes_posteriors_sum_to_one);
    TEST(test_naivebayes_matches_the_closed_form_gaussian_comparison);
    TEST(test_naivebayes_variance_floor_keeps_a_constant_feature_finite);
    TEST(test_naivebayes_option_errors_and_coexistence);
    TEST(test_logistic_separates_separable_data);
    TEST(test_logistic_probability_is_exactly_half_on_the_boundary);
    TEST(test_logistic_probability_is_monotone_and_normalised);
    TEST(test_logistic_coefficients_stay_finite_on_separable_data);
    TEST(test_logistic_declines_and_coexists);

    printf("All ml Classify tests passed.\n");
    return 0;
}
