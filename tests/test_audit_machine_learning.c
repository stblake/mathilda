/* Builtin-audit regressions: machine learning (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"

/* Classify#1 (high, crash): "Classify[{Rule[]}] segfaults: the first rule's
 * lhs is read without checking its argument count".  Malformed training data
 * leaves the call unevaluated; Predict's reader had the same read. */
static void test_classify_malformed_rules(void) {
    audit_chk("Classify#1", "Classify[{Rule[]}]", "Classify[{Rule[]}]");
    audit_chk("Classify#1", "Classify[{Rule[1]}]", "Classify[{Rule[1]}]");
    audit_chk("Classify#1", "Classify[{Rule[1, 2, 3]}] === Unevaluated[Classify[{Rule[1, 2, 3]}]]", "True");
    audit_chk("Classify#1", "Classify[{RuleDelayed[]}]", "Classify[{RuleDelayed[]}]");
    audit_chk("Classify#1", "Classify[{}]", "Classify[{}]");
    audit_chk("Classify#1", "Predict[{Rule[]}]", "Predict[{Rule[]}]");
    audit_chk("Classify#1", "Predict[{Rule[1]}]", "Predict[{Rule[1]}]");
    /* sanity: well-formed data still trains */
    audit_chk("Classify#1",
              "Classify[{1. -> \"a\", 2. -> \"a\", 10. -> \"b\", 11. -> \"b\"}][10.5]", "\"b\"");
}

/* LearnDistribution#1 (high): two points in the plane have a rank-1 sample covariance,
 * but the exact Cholesky test (pivot > 0) accepted the roundoff residue of 0.5 - 0.5
 * and PDF gave a spike of 2e7.  n <= dim is singular by rank, and a pivot must carry
 * more than 64 eps of its own variance. */
static void test_learn_distribution_singular(void) {
    audit_chk("LearnDistribution#1", "LearnDistribution[{{1., 2.}, {2., 3.}}]",
              "LearnDistribution[{{1.0, 2.0}, {2.0, 3.0}}]");
    audit_chk("LearnDistribution#1", "LearnDistribution[{{1., 2., 3.}, {2., 3., 5.}, {0., 1., 7.}}]",
              "LearnDistribution[{{1.0, 2.0, 3.0}, {2.0, 3.0, 5.0}, {0.0, 1.0, 7.0}}]");
    /* sanity: three non-collinear points fit, and the density at the mean is the
     * closed form 1/(2 Pi Sqrt[Det[cov]]) */
    audit_true("LearnDistribution#1",
               "Module[{d = {{1., 2.}, {2., 3.}, {3., 5.}}, ld},"
               " ld = LearnDistribution[d];"
               " Abs[PDF[ld, Mean[d]] - 1/(2 Pi Sqrt[Det[Covariance[d]]])] < 1*^-12]");
}

/* RandomVariate on a learned distribution returned a constant 0.0 -- ml_draw had no
 * case for the Multinormal / GaussianMixture / SmoothKernel kinds (and gave a scalar
 * for 2-D data).  Each is now sampled from its definition; checked by moments of
 * 20000 seeded draws against the fitted parameters (standard errors are ~1% here,
 * the tolerances several of them). */
static void test_random_variate_learned(void) {
    const char* ids[] = { "LearnDistribution#2", "LearnedDistribution#1", "RandomVariate#1" };
    for (int k = 0; k < 3; k++) {
        /* N(2, 1) fitted to {1, 2, 3} */
        audit_true(ids[k],
                   "SeedRandom[7]; Module[{s = RandomVariate[LearnDistribution[{1., 2., 3.}], 20000]},"
                   " Length[s] == 20000 && Abs[Mean[s] - 2] < 0.05 && Abs[Variance[s] - 1] < 0.05]");
        /* 2-D Multinormal: vectors, with the fitted mean and covariance */
        audit_true(ids[k],
                   "SeedRandom[7]; Module[{d = {{1., 2.}, {2., 3.}, {3., 5.}, {4., 4.}, {5., 7.}, {6., 8.}}, s},"
                   " s = RandomVariate[LearnDistribution[d], 20000];"
                   " Dimensions[s] == {20000, 2} && Max[Abs[Mean[s] - Mean[d]]] < 0.1"
                   " && Max[Abs[Flatten[Covariance[s] - Covariance[d]]]] < 0.25]");
        audit_true(ids[k],
                   "VectorQ[RandomVariate[LearnDistribution[{{1., 2.}, {2., 3.}, {3., 5.}}]], NumericQ]"
                   " && Length[RandomVariate[LearnDistribution[{{1., 2.}, {2., 3.}, {3., 5.}}]]] == 2");
        /* a two-cluster mixture: draws land near the clusters, almost none between */
        audit_true(ids[k],
                   "SeedRandom[7]; Module[{g, s},"
                   " g = LearnDistribution[{1., 2., 3., 4., 10., 11., 12., 13.}, Method -> \"GaussianMixture\"];"
                   " s = RandomVariate[g, 20000];"
                   " Abs[Mean[s] - 7] < 0.2 && Count[s, v_ /; 6.5 < v < 7.5] < 400"
                   " && Count[s, v_ /; v < 7] > 9000 && Count[s, v_ /; v > 7] > 9000]");
    }
    /* SmoothKernelDistribution#1: the KDE is the uniform mixture of N(x_i, h), so its
     * mean is the sample mean and its variance the (biased) sample variance plus h^2. */
    audit_true("SmoothKernelDistribution#1",
               "SeedRandom[7]; Module[{x = {1, 2, 3, 4, 10}, k, h, s},"
               " k = SmoothKernelDistribution[x]; h = k[[2, 1, 1]];"
               " s = RandomVariate[k, 40000];"
               " Abs[Mean[s] - 4] < 0.1 && Abs[Variance[s] - (Variance[x] 4/5 + h^2)] < 0.6]");
    audit_true("SmoothKernelDistribution#1",
               "Dimensions[RandomVariate[SmoothKernelDistribution[{{1., 2.}, {2., 1.}, {3., 3.}}], 4]] == {4, 2}");
    audit_true("SmoothKernelDistribution#1",
               "SeedRandom[3]; Length[Union[RandomVariate[SmoothKernelDistribution[{1., 2., 3.}], 50]]] == 50");
    /* sanity: specified distributions and the empty request are unchanged */
    audit_true("RandomVariate#1",
               "SeedRandom[1]; Module[{s = RandomVariate[NormalDistribution[5., 2.], 20000]},"
               " Abs[Mean[s] - 5] < 0.05 && Abs[StandardDeviation[s] - 2] < 0.05]");
    audit_chk("RandomVariate#1", "RandomVariate[LearnDistribution[{1., 2., 3.}], 0]", "{}");
}

/* LinearModelFit#1 (high): the normal equations on raw columns squared in the offset;
 * data exactly on y = 2 x - 199999997 at x ~ 1e8 gave slope 1.68.  The problem is
 * centred first.  The residuals of the fit are orthogonal to the design columns, here
 * zero since the data are collinear. */
static void test_linear_model_fit_offset(void) {
    audit_true("LinearModelFit#1",
               "Module[{c = LinearModelFit[{{100000000., 3.}, {100000001., 5.}, {100000002., 7.},"
               " {100000003., 9.}}][\"Coefficients\"]},"
               " Abs[c[[2]] - 2] < 1*^-9 && Abs[c[[1]] + 199999997] < 1*^-3]");
    audit_true("LinearModelFit#1",
               "Module[{x = {1., 2., 3., 4., 5.}, y = {1.1, 1.9, 3.2, 3.9, 5.1}, c, r},"
               " c = LinearModelFit[Transpose[{x, y}]][\"Coefficients\"];"
               " r = y - (c[[1]] + c[[2]] x);"
               " Abs[Total[r]] < 1*^-12 && Abs[r . x] < 1*^-12]");
    /* sanity: the documented examples */
    audit_chk("LinearModelFit#1", "LinearModelFit[{1. -> 3., 2. -> 5., 3. -> 7., 4. -> 9.}][\"Coefficients\"]",
              "{1.0, 2.0}");
    audit_true("LinearModelFit#1",
               "Max[Abs[LinearModelFit[{{1., 1., 6.}, {2., 1., 8.}, {1., 2., 9.}, {3., 2., 13.},"
               " {2., 3., 14.}}][\"Coefficients\"] - {1, 2, 3}]] < 1*^-12");
}

/* LinearModelFit#2 (high): LinearModelFit[data, funs, vars] was unevaluated.  It fits
 * c0 + c . funs by least squares and returns a FittedModel. */
static void test_linear_model_fit_basis(void) {
    audit_true("LinearModelFit#2",
               "Max[Abs[LinearModelFit[{{1, 3}, {2, 5}, {3, 7}}, x, x][\"BestFitParameters\"] - {1, 2}]] < 1*^-12");
    audit_true("LinearModelFit#2",
               "Max[Abs[LinearModelFit[{3, 5, 7, 9}, x, x][\"BestFitParameters\"] - {1, 2}]] < 1*^-12");
    audit_true("LinearModelFit#2",
               "Max[Abs[LinearModelFit[{{1, 3}, {2, 5}, {3, 7}}, {x}, x][\"BestFitParameters\"] - {1, 2}]] < 1*^-12");
    audit_true("LinearModelFit#2",
               "Abs[LinearModelFit[{{1, 3}, {2, 5}, {3, 7}}, x, x][10] - 21] < 1*^-12");
    audit_true("LinearModelFit#2",
               "Head[LinearModelFit[{{1, 3}, {2, 5}, {3, 7}}, x, x]] === FittedModel");
    audit_chk("LinearModelFit#2", "LinearModelFit[{{1, 3}, {2, 5}, {3, 7}}, x, x][\"BasisFunctions\"]",
              "{1, x}");
    /* the best fit is a polynomial in x with the fitted coefficients */
    audit_true("LinearModelFit#2",
               "Module[{b = LinearModelFit[{{1, 3}, {2, 5}, {3, 7}}, x, x][\"BestFit\"]},"
               " Abs[Coefficient[b, x] - 2] < 1*^-12 && Abs[(b /. x -> 0) - 1] < 1*^-12]");
    /* an exact quadratic is recovered from the basis {x, x^2} */
    audit_true("LinearModelFit#2",
               "Max[Abs[LinearModelFit[Table[{t, 1 + 2 t + 3 t^2}, {t, 0, 5}], {x, x^2}, x]"
               "[\"BestFitParameters\"] - {1, 2, 3}]] < 1*^-9");
    /* two variables, rows {x, y, z} */
    audit_true("LinearModelFit#2",
               "Max[Abs[LinearModelFit[{{1, 1, 6}, {2, 1, 8}, {1, 2, 9}, {3, 2, 13}, {2, 3, 14}},"
               " {x, y}, {x, y}][\"BestFitParameters\"] - {1, 2, 3}]] < 1*^-12");
    /* least squares, not interpolation: residuals orthogonal to the basis columns 1 and x */
    audit_true("LinearModelFit#2",
               "Module[{d = {{0, 1}, {1, 0}, {3, 2}, {5, 4}}, m, r},"
               " m = LinearModelFit[d, x, x]; r = d[[All, 2]] - (m /@ d[[All, 1]]);"
               " Abs[Total[r]] < 1*^-12 && Abs[r . d[[All, 1]]] < 1*^-12]");
    /* the offset case through the basis form too */
    audit_true("LinearModelFit#2",
               "Abs[LinearModelFit[{{100000000., 3.}, {100000001., 5.}, {100000002., 7.}, {100000003., 9.}},"
               " x, x][\"BestFitParameters\"][[2]] - 2] < 1*^-9");
}

int main(void) {
    audit_init();
    TEST(test_classify_malformed_rules);
    TEST(test_learn_distribution_singular);
    TEST(test_random_variate_learned);
    TEST(test_linear_model_fit_offset);
    TEST(test_linear_model_fit_basis);
    return audit_finish("audit_machine_learning");
}
