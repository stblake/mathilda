/*
 * test_geometry.c -- GEO-1: core computational-geometry builtins.
 *
 * Area / Perimeter / RegionCentroid / RegionMember over 2D Polygon, plus
 * ConvexHullRegion. Every expected string below is locked to a live Wolfram
 * kernel reference output (thoughts/shared/tickets/GEO-1/research.md, 2026-08-27),
 * rendered through mathilda's own printer (machine reals print 4.0 where WL
 * prints 4.; values identical).
 *
 * NOTE: geo_check deliberately does NOT use assert_eval_eq -- that helper
 * fails via libc assert(), which -DNDEBUG (any Release configuration)
 * compiles to a no-op, silently passing a failing test. ASSERT_STR_EQ
 * exit(1)s unconditionally, so CTest always sees a real verdict.
 */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

static void geo_check(const char* input, const char* expected) {
    struct Expr* parsed = parse_expression(input);
    ASSERT_MSG(parsed != NULL, "parse failed: %s", input);
    struct Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    evaluated = test_delist(evaluated);
    char* str = expr_to_string(evaluated);
    if (strcmp(str, expected) != 0)
        fprintf(stderr, "  input: %s\n", input);
    ASSERT_STR_EQ(str, expected);
    free(str);
    expr_free(evaluated);
}

/* Prefix form, for results containing a fully-expanded 400-digit integer whose
 * literal spelling would dwarf the assertion. Same loud-failure discipline. */
static void geo_check_prefix(const char* input, const char* prefix) {
    struct Expr* parsed = parse_expression(input);
    ASSERT_MSG(parsed != NULL, "parse failed: %s", input);
    struct Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    evaluated = test_delist(evaluated);
    char* str = expr_to_string(evaluated);
    ASSERT_MSG(strncmp(str, prefix, strlen(prefix)) == 0,
               "input: %s\n  expected prefix: %s\n  actual: %.120s...", input, prefix, str);
    free(str);
    expr_free(evaluated);
}

/* ---- Area (AC-1..AC-5, AC-25) ---- */
static void test_area(void) {
    geo_check("Area[Polygon[{{0,0},{1,0},{1/2,1/2}}]]", "1/4");                    /* AC-1 */
    geo_check("Area[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}]]", "10");             /* AC-2 concave */
    geo_check("Area[Polygon[{{0,0},{1.5,0},{1.5,1},{0,1}}]]", "1.5");              /* AC-3 machine */
    geo_check("Area[Polygon[{{0,0},{1,0},{1,1},{0,1},{0,0}}]]", "1");              /* AC-4 closed form */
    geo_check("Area[Polygon[{{0,0},{1,0}}]]", "Undefined");                        /* AC-5 degenerate */
    geo_check("Area[Polygon[NDArray[{{0.,0.},{2.,0.},{2.,2.},{0.,2.}}]]]", "4.0"); /* AC-25 */
}

/* ---- Perimeter (AC-6..AC-8, AC-26) ---- */
static void test_perimeter(void) {
    geo_check("Perimeter[Polygon[{{0,0},{1,0},{0,1}}]]", "2 + Sqrt[2]");           /* AC-6 exact */
    geo_check("Perimeter[Polygon[{{0,0},{3.,0},{3.,4.}}]]", "12.0");               /* AC-7 machine */
    geo_check("Perimeter[Polygon[{{0,0},{1,0}}]]", "Undefined");                   /* AC-8 degenerate */
    geo_check("Perimeter[Polygon[{{0,0},{1/2,0},{0,1}}]]", "3/2 + 1/2 Sqrt[5]");   /* AC-26 rational */
}

/* ---- RegionCentroid (AC-9, AC-10) ---- */
static void test_region_centroid(void) {
    geo_check("RegionCentroid[Polygon[{{0,0},{1,0},{0,1}}]]", "{1/3, 1/3}");       /* AC-9 */
    geo_check("RegionCentroid[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}]]", "{2, 7/5}"); /* AC-10 */
    /* zero-area polygon declines (documented deviation from WL) */
    geo_check("RegionCentroid[Polygon[{{0,0},{1,1},{2,2}}]]",
              "RegionCentroid[Polygon[{{0, 0}, {1, 1}, {2, 2}}]]");
}

/* ---- RegionMember (AC-11..AC-16) ---- */
static void test_region_member(void) {
    const char* sq = "Polygon[{{0,0},{2,0},{2,2},{0,2}}]";
    char buf[256];
    sprintf(buf, "RegionMember[%s, {1,1}]", sq);   geo_check(buf, "True");   /* AC-11 interior */
    sprintf(buf, "RegionMember[%s, {2,1}]", sq);   geo_check(buf, "True");   /* AC-12 edge */
    sprintf(buf, "RegionMember[%s, {0,0}]", sq);   geo_check(buf, "True");   /* AC-13 vertex */
    sprintf(buf, "RegionMember[%s, {3,1}]", sq);   geo_check(buf, "False");  /* AC-14 exterior */
    sprintf(buf, "RegionMember[%s, {1/2,1/2}]", sq); geo_check(buf, "True"); /* AC-16 rational */
    geo_check("RegionMember[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}], {2,3}]",
              "False");                                                      /* AC-15 notch */
    /* machine-path membership + NDArray query point */
    sprintf(buf, "RegionMember[%s, {0.5,0.5}]", sq); geo_check(buf, "True");
    sprintf(buf, "RegionMember[%s, NDArray[{1.,1.}]]", sq); geo_check(buf, "True");
}

/* ---- ConvexHullRegion (AC-17..AC-20, AC-24) ---- */
static void test_convex_hull(void) {
    geo_check("ConvexHullRegion[{{0,0},{2,0},{1,0},{2,2},{0,2},{1,1}}]",
              "Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}]");                  /* AC-17 dedup */
    geo_check("ConvexHullRegion[{{0,0},{1,0},{0,1},{2,2},{1/2,1/2}}]",
              "Polygon[{{0, 0}, {1, 0}, {2, 2}, {0, 1}}]");                  /* AC-18 CCW order */
    geo_check("ConvexHullRegion[{{0,0},{1,1},{2,2},{3,3}}]",
              "Line[{{0, 0}, {3, 3}}]");                                     /* AC-19 collinear */
    geo_check("ConvexHullRegion[{{1,2}}]", "Point[{1, 2}]");                 /* AC-20 single */
    geo_check("ConvexHullRegion[NDArray[{{0.,0.},{2.,0.},{2.,2.},{0.,2.},{1.,1.}}]]",
              "Polygon[{{0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}}]");  /* AC-24 */
}

/* ---- Composition + attributes (AC-22, AC-23) ---- */
static void test_composition_and_attributes(void) {
    geo_check("Area[ConvexHullRegion[{{0,0},{2,0},{1,0},{2,2},{0,2},{1,1}}]]", "4"); /* AC-22 */
    geo_check("Attributes[Area]", "{Protected}");                     /* AC-23 */
    geo_check("Attributes[ConvexHullRegion]", "{Protected}");
}

/* ---- Decline paths: out-of-scope input stays unevaluated (AC-21 shape) ---- */
static void test_declines(void) {
    geo_check("Area[Polygon[{{0,0},{a,0},{0,1}}]]",
              "Area[Polygon[{{0, 0}, {a, 0}, {0, 1}}]]");                    /* AC-21 symbolic */
    geo_check("Area[5]", "Area[5]");                                         /* wrong type */
    geo_check("Area[Polygon[{{0,0},{1,0},{0,1}}, {{0,0},{1,0}}]]",
              "Area[Polygon[{{0, 0}, {1, 0}, {0, 1}}, {{0, 0}, {1, 0}}]]");  /* holes: out of scope */
    geo_check("Perimeter[Polygon[{{0,0},{1,0},{0,x}}]]",
              "Perimeter[Polygon[{{0, 0}, {1, 0}, {0, x}}]]");
    geo_check("RegionMember[Polygon[{{0,0},{2,0},{2,2},{0,2}}], {x,1}]",
              "RegionMember[Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}], {x, 1}]");
    geo_check("ConvexHullRegion[{{0,0},{1,b}}]",
              "ConvexHullRegion[{{0, 0}, {1, b}}]");
    geo_check("ConvexHullRegion[{{0,0,0},{1,1,1}}]",
              "ConvexHullRegion[{{0, 0, 0}, {1, 1, 1}}]");                   /* 3D: out of scope */
    geo_check("RegionMember[Polygon[{{0,0},{1,0}}], {0,0}]",
              "RegionMember[Polygon[{{0, 0}, {1, 0}}], {0, 0}]");            /* degenerate poly */
}

/* ---- Findings from the independent adversarial review (2026-08-27) ----
 * Each of these was a wrong answer or a silent divergence before the fix; they
 * are pinned here so the fix cannot regress. */
static void test_review_findings(void) {
    /* Non-finite mirror: an exact coordinate too large for a double must never
     * be answered from the machine path (it produced True, the wrong answer,
     * via NaN cross products reading as "collinear"). */
    geo_check_prefix("RegionMember[Polygon[{{0,0},{10^400,0},{0,1}}], {1.,1.}]",
                     "RegionMember[Polygon[{{0, 0}, {1000000");
    /* ...while the fully exact spelling still answers, and correctly. */
    geo_check("RegionMember[Polygon[{{0,0},{10^400,0},{0,1}}], {1,1}]", "False");
    /* Mixed exact/machine with an unrepresentable coordinate: declined, rather
     * than a Line[{{0.0,0.0},{inf.0,0.0}}] whose "inf.0" is not re-parseable. */
    geo_check_prefix("ConvexHullRegion[{{0,0},{10^400,0},{0,1},{1.,1.}}]",
                     "ConvexHullRegion[{{0, 0}, {1000000");
    /* DISTINCT-vertex counting, as the docstring promises: a repeated vertex
     * does not make a polygon. Was a confident 0 / a there-and-back perimeter. */
    geo_check("Area[Polygon[{{0,0},{0,0},{1,1}}]]", "Undefined");
    geo_check("Perimeter[Polygon[{{0,0},{1,0},{1,0}}]]", "Undefined");
    /* The exact big-integer path is untouched by the finite check. */
    geo_check("Area[Polygon[{{0,0},{10^20,0},{0,10^20}}]]", "5000000000000000000000000000000000000000");
    /* Packed integer point lists stay EXACT (the packed-aware EXEMPT rationale). */
    geo_check("Area[Polygon[Table[{i, i^2}, {i, 0, 5}]]]", "20");
    /* A visible NDArray is float64 by construction, so a machine answer is
     * correct here -- matching Total[NDArray[{1,2,3}]] == 6.0. */
    geo_check("Area[Polygon[NDArray[{{0,0},{2,0},{2,2},{0,2}}]]]", "4.0");
}

/* ---- Repeated-evaluation loop ----
 *
 * This loop CANNOT detect a leak by itself -- an earlier version of it ran 700
 * evaluations asserting only `e != NULL` and was structurally incapable of
 * failing on the 320 bytes/call leak an independent review later found. It is
 * kept for two things it CAN do: catch a wrong answer under repetition (it now
 * asserts the values, not just non-NULL), and give an external leak checker a
 * workload to drive.
 *
 * The real leak gate is tests/scripts/geometry_leakcheck.sh (leaks --atExit on
 * macOS, valgrind on Linux, clean skip when neither exists). Run it after any
 * change to this module:
 *     bash tests/scripts/geometry_leakcheck.sh
 */
static void test_repeated_evaluation(void) {
    for (int i = 0; i < 100; i++) {
        /* Exact and machine paths, every head, plus a decline path. Values are
         * asserted so a repetition-dependent wrong answer fails here. */
        geo_check("Area[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}]]", "10");
        geo_check("Area[Polygon[{{0,0},{1.5,0},{1.5,1},{0,1}}]]", "1.5");
        geo_check("Perimeter[Polygon[{{0,0},{1/2,0},{0,1}}]]", "3/2 + 1/2 Sqrt[5]");
        geo_check("RegionCentroid[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}]]", "{2, 7/5}");
        geo_check("RegionMember[Polygon[{{0,0},{4,0},{4,4},{2,1},{0,4}}], {2,3}]", "False");
        geo_check("ConvexHullRegion[{{0,0},{2,0},{1,0},{2,2},{0,2},{1,1}}]",
                  "Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}]");
        geo_check("Area[Polygon[{{0,0},{a,0},{0,1}}]]",
                  "Area[Polygon[{{0, 0}, {a, 0}, {0, 1}}]]");
    }
}

/* ---- Scale and scratch reuse (invariants, not oracle rows) ----
 *
 * The exact cross-product predicate keeps its mpq_t temporaries for the whole
 * builtin call and reuses them for every evaluation inside it, so GMP grows
 * each limb buffer to that call's high-water mark and later, narrower
 * evaluations run on a buffer an earlier, wider one sized. Nothing above can
 * see a fault in that: every AC row is at most five vertices with small
 * coordinates, and a stale or aliased temporary needs one call that mixes
 * widths and runs the predicate many times.
 *
 * These cases build such a call and assert properties that hold for any
 * correct implementation -- containment, idempotence, hull-of-hull stability,
 * and exact/machine agreement -- so no Wolfram oracle is needed and the
 * assertions stay true if the vertex order convention is ever revisited. */

/* Deterministic LCG (Numerical Recipes constants); no rand(), so a failure
 * reproduces byte-for-byte on any machine. */
static unsigned long geo_lcg(unsigned long* st) {
    *st = *st * 1664525UL + 1013904223UL;
    return (*st >> 8) & 0xFFFFFFUL;
}

/* {{x,y},...} with `n` points whose coordinates span six orders of magnitude,
 * so one hull call runs the predicate over operands of very different widths.
 * Caller frees. */
static char* geo_wide_point_set(size_t n) {
    size_t cap = n * 64 + 8;
    char* buf = (char*)malloc(cap);
    ASSERT_MSG(buf != NULL, "out of memory building point set");
    size_t at = 0;
    unsigned long st = 20260827UL;
    at += (size_t)snprintf(buf + at, cap - at, "{");
    for (size_t i = 0; i < n; i++) {
        unsigned long a = geo_lcg(&st), b = geo_lcg(&st);
        /* Every fourth point is pushed out to ~10^12 so the set is not
         * uniformly narrow; the rest stay small. */
        long x = (long)(a % 1000), y = (long)(b % 1000);
        if (i % 4 == 0) { x *= 1000000000L; y *= 1000000L; }
        at += (size_t)snprintf(buf + at, cap - at, "%s{%ld,%ld}",
                               i ? "," : "", x, y);
    }
    snprintf(buf + at, cap - at, "}");
    return buf;
}

static void test_scale_and_scratch_reuse(void) {
    char* pts = geo_wide_point_set(240);
    size_t need = strlen(pts) * 2 + 256;
    char* expr = (char*)malloc(need);
    ASSERT_MSG(expr != NULL, "out of memory building expression");

    /* 1. Containment: every input point lies in or on its own hull. A predicate
     *    that returns a wrong sign for some operand width drops a point off the
     *    chain, and that point then reads as outside. */
    snprintf(expr, need,
             "Apply[And, Map[RegionMember[ConvexHullRegion[%s], #]&, %s]]",
             pts, pts);
    geo_check(expr, "True");

    /* 2. Idempotence: hulling the hull's own vertices reproduces it exactly. */
    snprintf(expr, need,
             "ConvexHullRegion[First[ConvexHullRegion[%s]]] === ConvexHullRegion[%s]",
             pts, pts);
    geo_check(expr, "True");

    /* 3. Repeat in the same process: the scratch block is per call, so a
     *    second call must not inherit anything from the first. */
    snprintf(expr, need,
             "ConvexHullRegion[%s] === ConvexHullRegion[%s]", pts, pts);
    geo_check(expr, "True");

    free(expr);
    free(pts);

    /* 4. Exact/machine agreement on one polygon written both ways -- the two
     *    paths share no code below geom_read_points, so this pins them to each
     *    other rather than to a printed constant. */
    geo_check("Area[Polygon[{{0,0},{1000000,0},{1000000,1000000},{0,1000000}}]]",
              "1000000000000");
    geo_check("Area[Polygon[{{0.,0.},{1000000.,0.},{1000000.,1000000.},{0.,1000000.}}]]",
              "1e+12");

    /* 5. Narrow evaluation after a wide one inside a SINGLE call: the small
     *    triangle's vertices are visited after the 10^18 vertex, on temporaries
     *    the wide vertex already grew. */
    geo_check("RegionMember[Polygon[{{0,0},{1000000000000000000,0},"
              "{1000000000000000000,1},{0,1}}], {1,1}]", "True");
    geo_check("RegionMember[Polygon[{{0,0},{1000000000000000000,0},"
              "{1000000000000000000,1},{0,1}}], {1,2}]", "False");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_area);
    TEST(test_perimeter);
    TEST(test_region_centroid);
    TEST(test_region_member);
    TEST(test_convex_hull);
    TEST(test_composition_and_attributes);
    TEST(test_declines);
    TEST(test_review_findings);
    TEST(test_repeated_evaluation);
    TEST(test_scale_and_scratch_reuse);

    printf("All geometry tests passed!\n");
    return 0;
}
