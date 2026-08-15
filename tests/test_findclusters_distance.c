/* test_findclusters_distance.c -- DistanceFunction actually selects a metric.
 *
 * It used to be validated and thrown away: fc_parse_distance_function checked the
 * name against four accepted spellings, returned true, and stored nothing. On a
 * LINE that was harmless and arguably correct -- every accepted metric is a
 * monotone transform of |a - b| there, so all four induce the same ordering of
 * gaps and the same partition. Above one dimension it is a wrong answer:
 * Manhattan and Euclidean rank pairs differently in the plane, and a user asking
 * for one and silently getting the other has no way to tell.
 *
 * Two properties are load-bearing and tested here rather than assumed:
 *
 *   Euclidean and SquaredEuclidean must ALWAYS agree. Squaring is monotone on
 *   non-negatives, so it preserves edge ranking, and the Automatic threshold
 *   compares against a multiple of the median where d > 3*median(d) exactly when
 *   d^2 > 9*median(d^2). They are one implementation for that reason -- which also
 *   avoids taking a square root of an exact value, since Sqrt[2] is irrational and
 *   an exact ordering that must compare irrationals is a much harder problem.
 *
 *   The exact and machine spanning-tree builders must agree under every metric.
 *   They are separate implementations (Expr arithmetic vs doubles) selected by
 *   whether the input is machine-precision, so a metric wired into one and not the
 *   other would make the partition depend on how the numbers were written.
 */
#include <stdio.h>
#include <string.h>
#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

static void test_metric_selection_changes_the_partition(void) {
    /* Three points chosen so the two metrics pick DIFFERENT heaviest MST edges,
     * with no ties to make the outcome arbitrary:
     *
     *   A=(0,0)  B=(0,11)  C=(8,6)
     *   Euclidean^2:  A-B=121  A-C=100  B-C=89   MST {89,100}, heaviest A-C
     *   Manhattan:    A-B=11   A-C=14   B-C=13   MST {11,13},  heaviest B-C
     *
     * Cutting the heaviest edge therefore isolates A under Euclidean and C under
     * Manhattan. An earlier attempt at this test used A=(0,0) B=(0,10) C=(7,7),
     * which ties the two Manhattan MST edges at 10 and lets both metrics return
     * the same partition -- it passed while proving nothing. */
    const char* p = "{{0, 0}, {0, 11}, {8, 6}}";
    (void)p;

    assert_eval_eq("FindClusters[{{0, 0}, {0, 11}, {8, 6}}, 2, "
                   "DistanceFunction -> EuclideanDistance]",
                   "{{{0, 0}}, {{0, 11}, {8, 6}}}", 0);
    assert_eval_eq("FindClusters[{{0, 0}, {0, 11}, {8, 6}}, 2, "
                   "DistanceFunction -> SquaredEuclideanDistance]",
                   "{{{0, 0}}, {{0, 11}, {8, 6}}}", 0);
    assert_eval_eq("FindClusters[{{0, 0}, {0, 11}, {8, 6}}, 2, "
                   "DistanceFunction -> ManhattanDistance]",
                   "{{{0, 0}, {0, 11}}, {{8, 6}}}", 0);
    /* Automatic is the Euclidean family, so it agrees with those two and not with
     * Manhattan. This row is what would have caught the option being ignored: it
     * used to equal the Manhattan row too. */
    assert_eval_eq("FindClusters[{{0, 0}, {0, 11}, {8, 6}}, 2]",
                   "{{{0, 0}}, {{0, 11}, {8, 6}}}", 0);
}

static void test_euclidean_and_squared_always_agree(void) {
    /* Several shapes, each asserted twice. If the two ever diverge, the threshold
     * factor has stopped tracking whether the weights are squared -- which is the
     * one way to break this pair, since d^2 > 9*median(d^2) is only equivalent to
     * d > 3*median(d) when the factor is squared to match. */
    const char* sets[] = {
        "{{0, 0}, {1, 4}, {9, 1}, {10, 6}, {2, 1}}",
        "{{0, 0}, {3, 3}, {8, 0}, {8, 5}}",
        "{{1, 1, 1}, {2, 1, 1}, {30, 30, 30}, {31, 30, 30}}",
    };
    for (size_t i = 0; i < sizeof(sets) / sizeof(sets[0]); i++) {
        char a[512], b[512];
        snprintf(a, sizeof a, "FindClusters[%s, DistanceFunction -> EuclideanDistance]", sets[i]);
        snprintf(b, sizeof b, "FindClusters[%s, DistanceFunction -> SquaredEuclideanDistance]", sets[i]);
        struct Expr* pa = parse_expression(a);
        struct Expr* pb = parse_expression(b);
        struct Expr* ea = evaluate(pa);
        struct Expr* eb = evaluate(pb);
        char* sa = expr_to_string(ea);
        char* sb = expr_to_string(eb);
        if (strcmp(sa, sb) != 0)
            fprintf(stderr, "FAIL: Euclidean and Squared disagree on %s\n  E: %s\n  S: %s\n",
                    sets[i], sa, sb);
        ASSERT(strcmp(sa, sb) == 0);
        free(sa); free(sb);
        expr_free(pa); expr_free(pb); expr_free(ea); expr_free(eb);
    }
}

static void test_one_dimensional_is_metric_invariant(void) {
    /* On a line all four agree, which is why storing the choice was unnecessary
     * before n-D mattered. Asserted so that wiring the option cannot quietly
     * change a 1-D answer -- the metric is deliberately not applied to scalar
     * input, whose weights are plain exact differences. */
    const char* expect = "{{1, 2, 3}, {10, 11, 12}, {25}}";
    assert_eval_eq("FindClusters[{1, 2, 3, 10, 11, 12, 25}]", expect, 0);
    assert_eval_eq("FindClusters[{1, 2, 3, 10, 11, 12, 25}, "
                   "DistanceFunction -> EuclideanDistance]", expect, 0);
    assert_eval_eq("FindClusters[{1, 2, 3, 10, 11, 12, 25}, "
                   "DistanceFunction -> SquaredEuclideanDistance]", expect, 0);
    assert_eval_eq("FindClusters[{1, 2, 3, 10, 11, 12, 25}, "
                   "DistanceFunction -> ManhattanDistance]", expect, 0);
}

static void test_exact_and_machine_builders_agree_under_manhattan(void) {
    /* The same points twice: as integers, which are machine-precision and take the
     * double builder, and divided by 3, which makes them exact rationals and takes
     * the Expr builder. The partition must be identical, because two builders for
     * one definition is a standing risk and Manhattan is newly wired into both.
     *
     * Scaling by a positive constant cannot change a partition under a metric that
     * is homogeneous of degree 1 -- which Manhattan is, and which is why the
     * comparison is meaningful rather than merely similar. */
    assert_eval_eq("FindClusters[{{0, 0}, {0, 11}, {8, 6}}, 2, "
                   "DistanceFunction -> ManhattanDistance]",
                   "{{{0, 0}, {0, 11}}, {{8, 6}}}", 0);
    assert_eval_eq("FindClusters[{{0, 0}, {0, 11/3}, {8/3, 2}}, 2, "
                   "DistanceFunction -> ManhattanDistance]",
                   "{{{0, 0}, {0, 11/3}}, {{8/3, 2}}}", 0);
    /* And the Euclidean family, same pair of builders. */
    assert_eval_eq("FindClusters[{{0, 0}, {0, 11}, {8, 6}}, 2, "
                   "DistanceFunction -> EuclideanDistance]",
                   "{{{0, 0}}, {{0, 11}, {8, 6}}}", 0);
    assert_eval_eq("FindClusters[{{0, 0}, {0, 11/3}, {8/3, 2}}, 2, "
                   "DistanceFunction -> EuclideanDistance]",
                   "{{{0, 0}}, {{0, 11/3}, {8/3, 2}}}", 0);
}

static void test_unknown_metric_declines(void) {
    /* An unrecognised name leaves the call unevaluated rather than falling back to
     * a default, so a typo cannot silently cluster by the wrong metric. */
    assert_eval_eq("FindClusters[{1, 2, 10}, DistanceFunction -> NoSuchDistance]",
                   "FindClusters[{1, 2, 10}, DistanceFunction -> NoSuchDistance]", 0);
    /* A quoted string is not a function: DistanceFunction takes a metric, the way
     * Wolfram's does, where Method takes a string. Getting these the wrong way
     * round is easy, so the distinction is pinned. */
    assert_eval_eq("FindClusters[{1, 2, 10}, DistanceFunction -> \"EuclideanDistance\"]",
                   "FindClusters[{1, 2, 10}, DistanceFunction -> \"EuclideanDistance\"]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_metric_selection_changes_the_partition);
    TEST(test_euclidean_and_squared_always_agree);
    TEST(test_one_dimensional_is_metric_invariant);
    TEST(test_exact_and_machine_builders_agree_under_manhattan);
    TEST(test_unknown_metric_declines);

    printf("All FindClusters DistanceFunction tests passed.\n");
    return 0;
}
