/* test_ml_dist.c -- RandomVariate, PDF, and distribution objects.
 *
 * A sampler is awkward to test because its output is supposed to be unpredictable. The
 * assertions here are the three things that ARE deterministic about it:
 *
 *   - reproducibility under SeedRandom, which is the property a user relies on and the
 *     one a sampler with its own generator would silently break;
 *   - the MOMENTS of a large sample, which pin the distribution's shape without pinning
 *     any individual draw;
 *   - the PDF, which is a closed-form function and can be checked against the same
 *     formula written out symbolically -- an independent implementation, since the
 *     symbolic side shares no code with the sampler.
 */
#include <stdio.h>
#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

static void test_pdf_matches_the_closed_form(void) {
    /* The independent cross-check for this iteration: the density at zero is
     * 1/Sqrt[2 Pi], computed by the existing symbolic machinery, which shares no code
     * with the C density function. */
    assert_eval_eq("Chop[PDF[NormalDistribution[], 0.] - 1./Sqrt[2. Pi]]", "0", 0);
    /* And at a general point, against the formula written out longhand. */
    assert_eval_eq("Chop[PDF[NormalDistribution[3., 2.], 1.5] - "
                   "Exp[-0.5 ((1.5 - 3.)/2.)^2]/(2. Sqrt[2. Pi])]", "0", 0);
    /* Symmetry is a further constraint the formula must satisfy. */
    assert_eval_eq("Chop[PDF[NormalDistribution[], -1.] - PDF[NormalDistribution[], 1.]]",
                   "0", 0);
}

static void test_pdf_threads_and_handles_uniform_support(void) {
    assert_eval_eq("Length[PDF[NormalDistribution[], {-1., 0., 1.}]]", "3", 0);
    /* A uniform density is flat inside its support and exactly zero outside -- the
     * outside case is what distinguishes a real density from a bare formula. */
    assert_eval_eq("PDF[UniformDistribution[{0., 2.}], {1., 5.}]", "{0.5, 0.0}", 0);
    assert_eval_eq("PDF[UniformDistribution[{0., 2.}], -0.001]", "0.0", 0);
}

static void test_draws_are_reproducible_under_seedrandom(void) {
    /* The property that makes a sampler usable at all. It also verifies the sampler
     * draws from the SAME stream as RandomReal: a generator of its own would ignore
     * SeedRandom entirely. */
    assert_eval_eq("Module[{a, b}, SeedRandom[42]; "
                   "a = RandomVariate[NormalDistribution[], 5]; SeedRandom[42]; "
                   "b = RandomVariate[NormalDistribution[], 5]; a === b]", "True", 0);
    assert_eval_eq("Module[{a, c}, SeedRandom[42]; "
                   "a = RandomVariate[NormalDistribution[], 5]; SeedRandom[7]; "
                   "c = RandomVariate[NormalDistribution[], 5]; a =!= c]", "True", 0);
}

static void test_reseeding_clears_the_gaussian_cache(void) {
    /* THE non-obvious one. Box-Muller produces deviates in pairs and caches the spare,
     * so without clearing that cache on SeedRandom the first draw after reseeding comes
     * from the PREVIOUS stream. The bug only shows after an ODD number of draws -- an
     * even number leaves the cache empty and hides it -- so both are asserted.
     *
     * This failure mode is worth a test rather than a comment because it half-works:
     * reproducibility would hold for RandomReal and fail only for RandomVariate. */
    assert_eval_eq("Module[{x, y}, SeedRandom[42]; x = RandomVariate[NormalDistribution[]];"
                   " SeedRandom[42]; y = RandomVariate[NormalDistribution[]]; x === y]",
                   "True", 0);
    assert_eval_eq("Module[{x, z}, SeedRandom[42]; x = RandomVariate[NormalDistribution[]];"
                   " SeedRandom[42]; RandomVariate[NormalDistribution[]];"      /* one draw */
                   " SeedRandom[42]; z = RandomVariate[NormalDistribution[]]; x === z]",
                   "True", 0);
    /* Three draws, then reseed: an odd count again, and the pair boundary lands
     * differently, so this is not the same case as above. */
    assert_eval_eq("Module[{x, z}, SeedRandom[42]; x = RandomVariate[NormalDistribution[]];"
                   " SeedRandom[42]; RandomVariate[NormalDistribution[], 3];"
                   " SeedRandom[42]; z = RandomVariate[NormalDistribution[]]; x === z]",
                   "True", 0);
}

static void test_sample_moments_match_the_distribution(void) {
    /* Pins the SHAPE without pinning any draw. 20000 samples puts the standard error of
     * the mean at sigma/141, so a tolerance of 0.1 on a sigma of 2 is roughly seven
     * standard errors -- loose enough not to be flaky, tight enough that a wrong
     * variance or a missing mean shift fails it. */
    assert_eval_eq("Module[{s}, SeedRandom[1]; "
                   "s = RandomVariate[NormalDistribution[5., 2.], 20000]; "
                   "Abs[Mean[s] - 5.] < 0.1 && Abs[StandardDeviation[s] - 2.] < 0.1]",
                   "True", 0);
    /* A uniform sample must also respect its SUPPORT, which a normal sampler scaled by
     * mistake would not. */
    assert_eval_eq("Module[{s}, SeedRandom[1]; "
                   "s = RandomVariate[UniformDistribution[{10., 20.}], 20000]; "
                   "Min[s] >= 10. && Max[s] < 20. && Abs[Mean[s] - 15.] < 0.1]",
                   "True", 0);
}

static void test_distribution_objects_print_in_full(void) {
    /* Deliberately the OPPOSITE convention from a fitted model, which prints elided. A
     * distribution is SPECIFIED by its parameters, so they are the information; a
     * fitted model's are an implementation detail. Same mechanism, opposite visibility,
     * and this row is what stops a later change from unifying them by accident. */
    assert_eval_eq("NormalDistribution[0., 1.]", "NormalDistribution[0.0, 1.0]", 0);
    assert_eval_eq("UniformDistribution[{2., 5.}]", "UniformDistribution[{2.0, 5.0}]", 0);
}

static void test_invalid_parameters_decline(void) {
    /* A non-positive standard deviation is not a distribution, and an inverted range is
     * not an interval. Declining beats returning NaNs that propagate silently through a
     * whole sample and only surface as a strange plot much later. */
    assert_eval_eq("Head[RandomVariate[NormalDistribution[0., -1.]]]", "RandomVariate", 0);
    assert_eval_eq("Head[RandomVariate[NormalDistribution[0., 0.]]]", "RandomVariate", 0);
    assert_eval_eq("Head[RandomVariate[UniformDistribution[{5., 1.}]]]",
                   "RandomVariate", 0);
    assert_eval_eq("Head[RandomVariate[CauchyDistribution[0., 1.]]]", "RandomVariate", 0);
    assert_eval_eq("Head[PDF[NormalDistribution[0., -1.], 0.]]", "PDF", 0);
    assert_eval_eq("Head[RandomVariate[NormalDistribution[], -3]]", "RandomVariate", 0);
    /* Zero draws is a valid request with an empty answer, not an error. */
    assert_eval_eq("RandomVariate[NormalDistribution[], 0]", "{}", 0);
}

#define D1 "{1., 2., 3., 4., 5., 6.}"
#define D2 "{{1., 2.}, {2., 3.}, {3., 5.}, {4., 4.}, {5., 7.}, {6., 8.}}"

static void test_one_dimensional_multinormal_equals_the_normal_pdf(void) {
    /* THE cross-check for LearnDistribution, and a strong one: the fitted multinormal
     * reaches its density through a Cholesky factor and a Mahalanobis distance, while
     * PDF[NormalDistribution[...]] evaluates the scalar closed form. The two share no
     * code, so agreement is evidence about both.
     *
     * Checked in the far tail (x = 10, about 3.5 sigma out) as well as at the mean,
     * because the tail is where a wrong log-determinant or a missing 2 pi would show up
     * as a small relative error rather than an obvious one. */
    const char* xs[] = { "1.", "3.5", "6.", "10." };
    for (size_t i = 0; i < sizeof(xs) / sizeof(xs[0]); i++) {
        char in[512];
        snprintf(in, sizeof in,
                 "Chop[PDF[LearnDistribution[" D1 "], {%s}] - "
                 "PDF[NormalDistribution[Mean[" D1 "], StandardDeviation[" D1 "]], %s]]",
                 xs[i], xs[i]);
        assert_eval_eq(in, "0", 0);
    }
}

static void test_the_fitted_density_integrates_to_one(void) {
    /* A density that is merely proportional to the right shape would pass every
     * agreement test against another density with the same normalisation error. This
     * one pins the normalisation absolutely: a trapezoid sum over +/- 6 sigma must be
     * 1 to three decimals. */
    assert_eval_eq("Module[{m = LearnDistribution[" D1 "], mu, sd, h, g},"
                   " mu = Mean[" D1 "]; sd = StandardDeviation[" D1 "]; h = 0.01;"
                   " g = Table[mu - 6. sd + h i, {i, 0, Round[12. sd/h]}];"
                   " Abs[h Total[Map[PDF[m, {#}] &, g]] - 1.] < 0.001]", "True", 0);
}

static void test_multinormal_in_two_dimensions(void) {
    /* The density must be maximal at the fitted mean and negligible far away -- the
     * minimum any correct multivariate density satisfies. */
    assert_eval_eq("Module[{m = LearnDistribution[" D2 "]},"
                   " PDF[m, {3.5, 4.8333333}] > PDF[m, {6., 2.}] > PDF[m, {50., 50.}]]",
                   "True", 0);
    /* A MATRIX of points threads, giving one density each. Note the asymmetry with the
     * scalar case: for a multinormal a flat list is ONE point, because its argument is
     * itself a list -- reading it as many points would silently treat each coordinate
     * as a separate observation. */
    assert_eval_eq("Length[PDF[LearnDistribution[" D2 "], {{3.5, 4.8}, {50., 50.}}]]",
                   "2", 0);
    assert_eval_eq("NumberQ[PDF[LearnDistribution[" D2 "], {3.5, 4.8}]]", "True", 0);
}

static void test_learned_distributions_elide_but_specified_ones_do_not(void) {
    /* Both directions of the deliberate contrast. A fitted distribution's parameters are
     * derived, so it elides; a specified distribution's are what the user wrote, so it
     * prints them. These two rows together are what stops a later change from unifying
     * the conventions. */
    assert_eval_eq("LearnDistribution[" D1 "]",
                   "LearnedDistribution[\"Multinormal\", <>]", 0);
    assert_eval_eq("NormalDistribution[3.5, 2.]", "NormalDistribution[3.5, 2.0]", 0);
    /* And FullForm still reveals the fitted parameters. */
    assert_eval_eq("FullForm[LearnDistribution[" D1 "]]",
                   "LearnedDistribution[\"Multinormal\", "
                   "List[List[3.5], List[3.5]], 1, 0]", 0);
}

static void test_singular_fits_decline(void) {
    /* Perfectly collinear columns give a singular covariance, so no density exists --
     * declining is the honest answer, where a pseudo-inverse would invent one. */
    assert_eval_eq("Head[LearnDistribution[{{1., 2.}, {2., 4.}, {3., 6.}}]]",
                   "LearnDistribution", 0);
    /* One observation has no dispersion to fit. */
    assert_eval_eq("Head[LearnDistribution[{{1., 2.}}]]", "LearnDistribution", 0);
    assert_eval_eq("Head[LearnDistribution[" D1 ", Method -> \"Poisson\"]]",
                   "LearnDistribution", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_pdf_matches_the_closed_form);
    TEST(test_pdf_threads_and_handles_uniform_support);
    TEST(test_draws_are_reproducible_under_seedrandom);
    TEST(test_reseeding_clears_the_gaussian_cache);
    TEST(test_sample_moments_match_the_distribution);
    TEST(test_distribution_objects_print_in_full);
    TEST(test_invalid_parameters_decline);
    TEST(test_one_dimensional_multinormal_equals_the_normal_pdf);
    TEST(test_the_fitted_density_integrates_to_one);
    TEST(test_multinormal_in_two_dimensions);
    TEST(test_learned_distributions_elide_but_specified_ones_do_not);
    TEST(test_singular_fits_decline);

    printf("All ml distribution tests passed.\n");
    return 0;
}
