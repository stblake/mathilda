/* test_findclusters_scalar_pin.c -- what every FindClusters method returns on
 * ONE-DIMENSIONAL input, pinned before the n-dimensional ports touch anything
 * those methods share.
 *
 * WHY THIS FILE EXISTS, and why it is separate from test_list.c's FindClusters
 * coverage. Eight of the ten methods are one-dimensional by algorithm: they reach
 * their data only through fc_sorted_values, which dereferences
 * d->val[d->order[j]] -- the sorted permutation and the double projection, both
 * NULL for vector input. Making them work in n dimensions means giving them a
 * distance function, a length scale that is not the median adjacent gap, and a
 * cluster-assignment contract that is not "id per sorted position". All three are
 * shared with the 1-D paths that work today.
 *
 * A refactor of that kind fails quietly. The methods still return a partition, it
 * is still plausible, and nothing crashes -- so the only way to know the 1-D
 * answers did not move is to have written them down first. Every value below was
 * produced by running the built binary, not derived by hand.
 *
 * These are REGRESSION pins, not specifications. They record present behaviour so
 * a change to it has to be deliberate. Where a partition below is poor (KMedoids
 * splits {10,11,12,25} badly, JarvisPatrick merges two obvious groups), that is
 * the current kernel's behaviour and improving it is a separate, intentional
 * change -- which this file will then correctly flag.
 */
#include <stdio.h>
#include <string.h>
#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

static void test_find_clusters_scalar_pin(void) {
    /* One dataset for every method, chosen so the answer is not in doubt: three
     * groups at 1-3, 10-12 and 25, with gaps far wider than any within-group
     * spacing. A method that returns anything but those three groups is either
     * broken or has an opinion worth reading. */
    const char* data = "{1, 2, 3, 10, 11, 12, 25}";
    const char* three = "{{1, 2, 3}, {10, 11, 12}, {25}}";

    struct { const char* in; const char* out; } cases[] = {
        /* ---- The two dimension-general methods. Single linkage IS cut-the-MST,
         * so these share one implementation and must agree. ---- */
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, 3, Method -> \"Agglomerate\"]",  three},
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, 3, Method -> \"SpanningTree\"]", three},

        /* ---- Fixed and bounded count. ---- */
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, 3, Method -> \"KMeans\"]",       three},
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, UpTo[3], Method -> \"KMeans\"]", three},
        /* KMedoids' 1-D medoid is the median of a contiguous sorted run, an O(1)
         * index lookup, and it lands here rather than on the obvious grouping.
         * Pinned as-is: the n-D port introduces a real medoid search, and when
         * that changes this answer the change should be visible and argued. */
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, 3, Method -> \"KMedoids\"]",
         "{{1, 2, 3}, {10}, {11, 12, 25}}"},
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, UpTo[3], Method -> \"KMedoids\"]",
         "{{1, 2, 3}, {10}, {11, 12, 25}}"},

        /* ---- Spectral takes Automatic or UpTo[k] but NOT a fixed count, which is
         * transcribed from Mathematica's own error messages. The fixed-count form
         * declining is itself the behaviour under test. ---- */
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, Method -> \"Spectral\"]",        three},
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, UpTo[3], Method -> \"Spectral\"]", three},
        /* The echoed input keeps its quotes here, where Print would have stripped
         * them -- the expected string is the printed EXPRESSION, not Print output. */
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, 3, Method -> \"Spectral\"]",
         "FindClusters[{1, 2, 3, 10, 11, 12, 25}, 3, Method -> \"Spectral\"]"},

        /* ---- The density family: Automatic count only. ---- */
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, Method -> \"DBSCAN\"]",          three},
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, Method -> \"GaussianMixture\"]", three},
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, Method -> \"MeanShift\"]",       three},
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, Method -> \"NeighborhoodContraction\"]", three},
        /* JarvisPatrick merges the first two groups at its default k. Its shared-
         * neighbour test walks contiguous windows of the sorted order, so "near"
         * means "close in rank", not "close in value" -- exactly the assumption
         * the n-D port has to replace. */
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, Method -> \"JarvisPatrick\"]",
         "{{1, 2, 3, 10, 11, 12}, {25}}"},

        /* ---- Method sub-options reach the kernels and change the answer. ---- */
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, "
         "Method -> {\"DBSCAN\", \"NeighborhoodRadius\" -> 2, \"MinPoints\" -> 2}]", three},
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, "
         "Method -> {\"JarvisPatrick\", \"NeighborCount\" -> 3}]",
         "{{1, 2, 3, 10, 11, 12}, {25}}"},
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}, "
         "Method -> {\"MeanShift\", \"NeighborhoodRadius\" -> 3}]", three},

        /* ---- Automatic method and count. ---- */
        {"FindClusters[{1, 2, 3, 10, 11, 12, 25}]", three},

        /* ---- The three properties a shared refactor is most likely to break. ----
         *
         * Exact ordering: the sort and the gap selection compare the ELEMENTS, so
         * rationals with bigint components order correctly. A refactor that routes
         * ordering through the double projection passes every test above and fails
         * this one. */
        {"FindClusters[{1/10^25, 2/10^25, 1}, 2]",
         "{{1/10000000000000000000000000, 1/5000000000000000000000000}, {1}}"},
        /* A one-component POINT is not a SCALAR: the elements are still Lists, and
         * conflating the two sent 1-vectors down the scalar path where every
         * arithmetic step threaded over them. dim == 1 with kind == POINT is a real
         * combination and the ports must preserve it. */
        {"FindClusters[{{1}, {2}, {100}}]", "{{{1}, {2}}, {{100}}}"},
        /* Equal elements are never split, whatever the method. Enforced globally by
         * the fold over zero-weight tree edges rather than trusted to ten kernels
         * -- three of them broke it in review. The fold is kind-agnostic, so it
         * should keep holding through the ports for free; this asserts that it does. */
        {"FindClusters[{7, 7, 7, 7, 1, 100}]", "{{7, 7, 7, 7}, {1}, {100}}"},
    };

    (void)data;
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++)
        assert_eval_eq(cases[i].in, cases[i].out, 0);
}

int main(void) {
    /* The evaluator needs its symbol table and builtins before any of this can
     * parse, let alone evaluate. Omitting these two segfaults before the first
     * printf reaches the terminal, which is a confusing way to learn it. */
    symtab_init();
    core_init();

    printf("Running test: test_find_clusters_scalar_pin\n");
    test_find_clusters_scalar_pin();
    printf("All FindClusters scalar pins passed.\n");
    return 0;
}
