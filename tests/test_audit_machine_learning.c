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

/* Wave 4: distribution closed forms, batch application, feature reading. */
static void test_wave4_machine_learning(void) {
    /* exact and symbolic PDF */
    audit_true("PDF#1", "PDF[NormalDistribution[], 0] === 1/Sqrt[2 Pi]");
    audit_chk("PDF#1", "{PDF[UniformDistribution[{0, 2}], 1], PDF[UniformDistribution[{0, 2}], 3]}", "{1/2, 0}");
    audit_true("NormalDistribution#1", "PossibleZeroQ[PDF[NormalDistribution[0, 1], 1] - E^(-1/2)/Sqrt[2 Pi]] && "
               "PossibleZeroQ[PDF[NormalDistribution[1/2, 3/2], 1/2] - Sqrt[2/Pi]/3]");
    audit_chk("UniformDistribution#1", "PDF[UniformDistribution[{0, 3}], 1]", "1/3");
    audit_true("NormalDistribution#2", "PossibleZeroQ[PDF[NormalDistribution[m, s], x] - "
               "E^(-(-m + x)^2/(2 s^2))/(Sqrt[2 Pi] s)]");
    audit_true("PDF#2", "PossibleZeroQ[PDF[NormalDistribution[], x] - E^(-x^2/2)/Sqrt[2 Pi]]");
    audit_chk("UniformDistribution#2", "PDF[UniformDistribution[{0, 2}], x]",
              "Piecewise[{{1/2, 0 <= x <= 2}}, 0]");
    audit_chk("PDF#2", "PDF[UniformDistribution[{0, 2}], x]", "Piecewise[{{1/2, 0 <= x <= 2}}, 0]");
    audit_true("PDF#2", "With[{p = PDF[NormalDistribution[], {0., x}]}, Abs[p[[1]] - 0.3989422804014327] < 1*^-15 && "
               "PossibleZeroQ[p[[2]] - E^(-x^2/2)/Sqrt[2 Pi]]]");
    /* moments, CDF, Quantile, Median */
    audit_chk("NormalDistribution#3", "{Mean[NormalDistribution[m, s]], Variance[NormalDistribution[m, s]], "
              "StandardDeviation[NormalDistribution[2, 3]]}", "{m, s^2, 3}");
    audit_chk("UniformDistribution#4", "{Mean[UniformDistribution[{0, 2}]], Variance[UniformDistribution[{0, 2}]], "
              "CDF[UniformDistribution[{0, 2}], 1], Quantile[UniformDistribution[{0, 4}], 1/4], "
              "Median[UniformDistribution[{1, 3}]]}", "{1, 1/3, 1/2, 1, 2}");
    /* numeric-constant parameters */
    audit_true("RandomVariate#2", "Head[RandomVariate[NormalDistribution[Sqrt[2], 1]]] === Real && "
               "Head[RandomVariate[NormalDistribution[Pi, E]]] === Real");
    audit_true("UniformDistribution#3", "With[{r = RandomVariate[UniformDistribution[{0, Pi}]]}, "
               "Head[r] === Real && 0 <= r <= Pi]");
    audit_chk("UniformDistribution#3", "{PDF[UniformDistribution[{0, Pi}], 1], PDF[UniformDistribution[{0, 2}], Pi]}",
              "{1/Pi, 0}");
    /* array shapes */
    audit_chk("RandomVariate#3", "{Dimensions[RandomVariate[NormalDistribution[], {2, 3}]], "
              "Length[RandomVariate[UniformDistribution[{0, 1}], {3}]]}", "{{2, 3}, 3}");
    /* univariate learned distributions take a scalar point */
    audit_run("auditM = LearnDistribution[{1., 2., 3., 4., 5., 6.}]; "
              "auditK = SmoothKernelDistribution[Table[1. i/10., {i, 0, 40}]]; "
              "auditK2 = SmoothKernelDistribution[{1, 2, 3, 4, 10}]");
    audit_true("LearnDistribution#3", "Abs[PDF[auditM, 3.5] - 0.213244] < 1*^-6");
    audit_true("LearnedDistribution#2", "PDF[auditM, 3.5] == PDF[auditM, {3.5}]");
    audit_true("SmoothKernelDistribution#2", "Abs[PDF[auditK, 2.0] - PDF[auditK, {2.0}]] < 1*^-15");
    audit_true("SmoothKernelDistribution#3", "Abs[Mean[auditK2] - 4.] < 1*^-12 && Variance[auditK2] > 10. && "
               "0 < CDF[auditK2, 3.0] < 1 && Abs[CDF[auditK2, Quantile[auditK2, 0.5]] - 0.5] < 1*^-9");
    audit_run("ClearAll[auditM, auditK, auditK2]");
    /* batch application */
    audit_run("auditC = Classify[{{0., 0.} -> \"red\", {1., 0.} -> \"red\", {0., 1.} -> \"red\", "
              "{10., 10.} -> \"blue\", {11., 10.} -> \"blue\", {10., 11.} -> \"blue\"}]; "
              "auditD = Classify[{1. -> \"a\", 2. -> \"a\", 8. -> \"b\", 9. -> \"b\"}]; "
              "auditP = Predict[{1. -> 3., 2. -> 5., 3. -> 7., 4. -> 9.}]");
    audit_chk("ClassifierFunction#1", "auditC[{{0.5, 0.5}, {10.5, 10.5}}]", "{\"red\", \"blue\"}");
    audit_chk("ClassifierFunction#1", "auditD[{1.5, 8.5}]", "{\"a\", \"b\"}");
    audit_chk("Classify#4", "auditC[{{0.5, 0.5}, {10.5, 10.5}}]", "{\"red\", \"blue\"}");
    audit_true("Predict#2", "Max[Abs[auditP[{1., 2.}] - {3., 5.}]] < 1*^-12");
    audit_true("PredictorFunction#1", "Max[Abs[auditP[{1., 2., 3.}] - {3., 5., 7.}]] < 1*^-12");
    audit_run("ClearAll[auditC, auditD, auditP]");
    /* features: exact numbers yes, arbitrary heads no */
    audit_chk("Classify#2", "Head[Classify[{1/2 -> \"a\", 3 -> \"b\"}]]", "ClassifierFunction");
    audit_chk("Classify#3", "Head[Classify[{g[9] -> \"a\", 1. -> \"b\"}]]", "Classify");
    audit_chk("Classify#3", "Head[Classify[{{1., 2.} -> \"a\", f[30, 40] -> \"b\"}]]", "Classify");
    /* a relative pivot tolerance */
    audit_true("Predict#1", "Max[Abs[Predict[{1.*^-7 -> 1., 2.*^-7 -> 2., 3.*^-7 -> 3., 4.*^-7 -> 4.}]"
               "[\"Coefficients\"] - {0., 1.*^7}]] < 1*^-3");
    audit_true("LinearModelFit#3", "Head[LinearModelFit[{{1.*^-7, 1.}, {2.*^-7, 2.}, {3.*^-7, 3.}}]] =!= LinearModelFit");
    /* Standardize: exact, and with location / scale functions */
    audit_chk("Standardize#1", "{Standardize[{1, 2, 3}], Standardize[{1, 2, 4}], Standardize[{10^30, 2*10^30, 3*10^30}]}",
              "{{-1, 0, 1}, {-4/Sqrt[21], -1/Sqrt[21], 5/Sqrt[21]}, {-1, 0, 1}}");
    audit_chk("Standardize#2", "{Standardize[{1, 2, 3}, Mean, StandardDeviation], Standardize[{1, 2, 3}, Median]}",
              "{{-1, 0, 1}, {-1, 0, 1}}");
    /* InputForm round-trips a trained model */
    audit_true("DimensionReducerFunction#1", "With[{r = DimensionReduction[{{100., 200.}, {102., 205.}, "
               "{104., 210.}, {106., 215.}, {108., 220.}}, 1]}, ToExpression[ToString[r, InputForm]] === r]");
}

/* Wave 5 (G6, low): classifier / predictor properties, learned and specified
 * distributions. */
static void test_wave5_ml(void) {
    audit_run("w5c = Classify[{1. -> \"a\", 2. -> \"a\", 8. -> \"b\", 9. -> \"b\"}]");
    audit_chk("ClassifierFunction#2", "w5c[1.5, \"Decision\"]", "\"a\"");
    audit_chk("ClassifierFunction#2", "w5c[1.5, \"TopProbabilities\"]", "{\"a\" -> 1.0}");
    audit_chk("ClassifierFunction#2", "w5c[8.5, {\"TopProbabilities\", 1}]", "{\"b\" -> 1.0}");
    audit_chk("Classify#5", "StringContainsQ[Information[Classify], \"RandomForest\"] && "
              "!StringContainsQ[Information[Classify], \"only method\"]", "True");
    audit_chk("PredictorFunction#2", "StringContainsQ[Information[PredictorFunction], \"featureCount, extra\"]", "True");
    audit_run("w5p = Predict[{1. -> 3., 2. -> 5., 3. -> 7., 4. -> 9.}]");
    audit_chk("PredictorFunction#3", "w5p[\"Properties\"]", "{\"Coefficients\", \"FeatureCount\", \"Method\", \"Properties\"}");
    audit_chk("Predict#4", "Predict[{1., 2., 3.} -> {1., 2., 3.}] === Predict[{1. -> 1., 2. -> 2., 3. -> 3.}]", "True");
    audit_chk("Predict#4", "Predict[<|1. -> 1., 2. -> 2., 3. -> 3.|>] === Predict[{1. -> 1., 2. -> 2., 3. -> 3.}]", "True");
    /* DimensionReduce / DimensionReduction */
    audit_chk("DimensionReduce#1", "DimensionReduce[{{1, 2, 3, 4}, {5, 6, 7, 9}}, 3]", "DimensionReduce[{{1, 2, 3, 4}, {5, 6, 7, 9}}, 3]");
    audit_chk("DimensionReduce#1", "Dimensions[DimensionReduce[{{1, 2, 3, 4}, {5, 6, 7, 9}}, 1]]", "{2, 1}");
    audit_chk("DimensionReduce#2",
              "DimensionReduce[{{1, 2}, {3, 4}, {5, 7}}, 1, Method -> Automatic] === DimensionReduce[{{1, 2}, {3, 4}, {5, 7}}, 1]",
              "True");
    audit_chk("DimensionReduction#1",
              "DimensionReduction[{{1, 2}, {3, 4}, {5, 7}}, 1, Method -> \"PrincipalComponentsAnalysis\"] === "
              "DimensionReduction[{{1, 2}, {3, 4}, {5, 7}}, 1]", "True");
    /* LearnDistribution */
    audit_chk("LearnDistribution#4",
              "LearnDistribution[{1., 2., 3.}, Method -> Automatic] === LearnDistribution[{1., 2., 3.}]", "True");
    audit_chk("LearnDistribution#5", "PDF[LearnDistribution[{\"a\", \"b\", \"a\"}], \"a\"] == 2/3", "True");
    audit_run("w5d = LearnDistribution[{\"r\", \"r\", \"r\", \"b\"}, Method -> \"ContingencyTable\"]");
    audit_chk("LearnedDistribution#3", "MemberQ[{\"r\", \"b\"}, RandomVariate[w5d]]", "True");
    audit_chk("LearnedDistribution#3",
              "SeedRandom[7]; Count[RandomVariate[w5d, 4000], \"r\"] > 2700", "True");
    /* specified distributions */
    audit_chk("NormalDistribution#4", "InverseCDF[NormalDistribution[], 0.975] == Quantile[NormalDistribution[], 0.975]", "True");
    audit_chk("NormalDistribution#4", "InverseCDF[NormalDistribution[], 1/2]", "0");
    audit_chk("NormalDistribution#6", "Quiet[Check[PDF[NormalDistribution[0, -1], 0.], fired, NormalDistribution::posprm]]", "fired");
    audit_chk("PDF#6", "Quiet[Check[PDF[NormalDistribution[0, -1], 0.], fired, NormalDistribution::posprm]]", "fired");
    audit_chk("PDF#6", "Quiet[Check[PDF[UniformDistribution[{2, 0}], 1.], fired, UniformDistribution::ordr]]", "fired");
    audit_chk("PDF#6", "Quiet[PDF[UniformDistribution[{2, 0}], 1.]]", "PDF[UniformDistribution[{2, 0}], 1.0]");
    audit_chk("PDF#3", "Precision[PDF[NormalDistribution[], 0.5`30]] >= 29", "True");
    audit_chk("PDF#3", "Abs[PDF[NormalDistribution[], 0.5`30] - N[Exp[-1/8]/Sqrt[2 Pi], 30]] < 10^-28", "True");
    audit_chk("PDF#5", "PDF[NormalDistribution[]][0.] == PDF[NormalDistribution[], 0.]", "True");
    audit_chk("PDF#5", "CDF[NormalDistribution[]][0]", "1/2");
    audit_chk("UniformDistribution#7", "PDF[UniformDistribution[{{0, 1}, {0, 1}}], {0.5, 0.5}]", "1");
    audit_chk("UniformDistribution#7", "PDF[UniformDistribution[{{0, 1}, {0, 2}}], {0.5, 3}]", "0");
    audit_chk("UniformDistribution#7",
              "SeedRandom[2]; p = RandomVariate[UniformDistribution[{{0, 1}, {0, 1}}]]; "
              "Length[p] == 2 && And @@ Thread[0 <= p <= 1]", "True");
}

int main(void) {
    audit_init();
    TEST(test_wave4_machine_learning);
    TEST(test_classify_malformed_rules);
    TEST(test_learn_distribution_singular);
    TEST(test_random_variate_learned);
    TEST(test_linear_model_fit_offset);
    TEST(test_linear_model_fit_basis);
    TEST(test_wave5_ml);
    return audit_finish("audit_machine_learning");
}
