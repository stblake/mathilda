/* Builtin-audit regressions: geometry (review/issues/<Head>.md; tests/test_audit.h). */
#include "test_audit.h"

/* RegionMember#1 (high): "RegionMember[Annulus[...], {p1, p2}] returns a
 * nonsense inequality instead of a list of Booleans". */
static void test_region_member_point_list(void) {
    audit_chk("RegionMember#1", "RegionMember[Annulus[{0, 0}, {1, 2}], {{3/2, 0}, {0, 0}}]", "{True, False}");
    audit_chk("RegionMember#1", "RegionMember[Annulus[{0, 0}, {1, 2}], {3/2, 0}]", "True");
    audit_chk("RegionMember#1", "RegionMember[Disk[{0, 0}, 2], {{1, 1}, {3, 0}, {0, -2}}]", "{True, False, True}");
    audit_chk("RegionMember#1", "RegionMember[Polygon[{{0, 0}, {1, 0}, {0, 1}}], {{1/10, 1/10}, {2, 2}}]", "{True, False}");
    audit_chk("RegionMember#1", "RegionMember[Polygon[{{0, 0}, {1, 0}, {0, 1}}], {1/10, 1/10}]", "True");
}

/* RegionMember#2 (high): "RegionMember on a Disk is left unevaluated". */
static void test_region_member_disk(void) {
    audit_chk("RegionMember#2", "RegionMember[Disk[], {0, 0}]", "True");
    audit_chk("RegionMember#2", "RegionMember[Disk[{0, 0}, 2], {1, 1}]", "True");
    audit_chk("RegionMember#2", "RegionMember[Disk[{0, 0}, 2], {2, 0}]", "True");
    audit_chk("RegionMember#2", "RegionMember[Disk[{0, 0}, 2], {3/2, 3/2}]", "False");
    audit_chk("RegionMember#2", "RegionMember[Disk[{1, 1}], {1.5, 1.5}]", "True");
    audit_chk("RegionMember#2", "RegionMember[Disk[{1, 1}, {2, 1}], {{2.9, 1}, {1, 2.1}}]", "{True, False}");
    audit_chk("RegionMember#2", "RegionMember[Disk[{0, 0}, 1, {0, Pi/2}], {{1/2, 1/2}, {-1/2, 1/2}}]", "{True, False}");
    audit_chk("RegionMember#2", "RegionMember[Disk[], {qx, qy}]", "Element[qx | qy, Reals] && qx^2 + qy^2 <= 1");
}

/* ================================================================== */
/* Wave 4 G3: measures never negative, sectors, the further regions   */
/* (src/geometry_curves.c gc_span / gc_absnum; src/geometry_regions.c) */
/* ================================================================== */

/* A sector's span is |t2 - t1| capped at a full turn; a numeric radius
 * measures as |r|. */
static void test_wave4_sector_span(void) {
    audit_true("Annulus#1", "Area[Annulus[{0,0},{1,2},{Pi/2,0}]] === 3 Pi/4");
    audit_true("Annulus#2", "Area[Annulus[{0,0},{1,2},{0,4Pi}]] === Area[Annulus[{0,0},{1,2}]] === 3 Pi");
    audit_true("Area#4", "Area[Annulus[{0,0},{1,2},{Pi,0}]] === 3 Pi/2");
    audit_true("Area#4", "Area[Annulus[{0,0},{1,2},{0,3 Pi}]] === 3 Pi");
    audit_true("RegionMeasure#2", "RegionMeasure[Annulus[{0,0},{1,2},{0,4 Pi}]] === 3 Pi");
    audit_true("Circle#1", "ArcLength[Circle[{0,0},1,{Pi/2,0}]] === Pi/2");
    audit_true("Circle#1", "RegionMeasure[Circle[{0,0},1,{Pi,0}]] === Pi");
    audit_true("ArcLength#2", "ArcLength[Circle[{0,0},-2]] === 4 Pi");
    audit_true("ArcLength#2", "ArcLength[Circle[{0,0},1,{Pi/2,0}]] === Pi/2");
    audit_true("RegionMeasure#1", "RegionMeasure[Circle[{0,0},-1]] === 2 Pi");
    audit_true("RegionMeasure#1", "RegionMeasure[Disk[{0,0},-1]] === Pi");
    audit_true("Disk#3", "Perimeter[Disk[{0,0},-2]] === 4 Pi");
    audit_true("Disk#3", "Area[Disk[{0,0},{-2,3}]] === 6 Pi");
    audit_true("Disk#3", "Area[Disk[{0,0},-2]] === 4 Pi");
    /* a full-turn annular sector has no radial edges */
    audit_true("Annulus#2", "Perimeter[Annulus[{0,0},{1,2},{0,4 Pi}]] === 6 Pi");
    /* the reversed sector contains the same points */
    audit_true("Annulus#1", "RegionMember[Annulus[{0,0},{1,2},{Pi/2,0}], {1, 1}]");
    /* Perimeter[Disk[]] no longer rewrites the caller's Disk into a Circle */
    audit_true("Disk#3", "d = Disk[]; Perimeter[d]; d === Disk[]");
}

static void test_wave4_disk_sector(void) {
    audit_true("Area#1", "Area[Disk[{0,0},2,{0,Pi}]] === 2 Pi");
    audit_true("Area#1", "Area[Disk[{0,0},{1,2},{0,Pi}]] === Pi");
    audit_true("Disk#1", "Area[Disk[{0,0},1,{0,Pi/2}]] === Pi/4");
    audit_true("Disk#1", "Perimeter[Disk[{0,0},1,{0,Pi/2}]] === 2 + Pi/2");
    audit_true("Disk#1", "RegionCentroid[Disk[{0,0},1,{0,Pi}]] === {0, 4/(3 Pi)}");
    audit_true("Disk#1", "RegionCentroid[Annulus[{0,0},{1,2},{0,Pi}]] === {0, 28/(9 Pi)}");
}

/* ArcLength over an infinite range: Integrate's divergence verdict, never a
 * quadrature's last partial sum. */
static void test_wave4_arclength_infinite(void) {
    audit_chk("ArcLength#1", "Quiet[ArcLength[{t,t}, {t, 0, Infinity}]]", "Integrate[Sqrt[2], {t, 0, Infinity}]");
    audit_chk("ArcLength#1", "Check[ArcLength[{t,t}, {t, 0, Infinity}], fired, Integrate::idiv]", "fired");
}

static void test_wave4_rectangle_triangle(void) {
    audit_true("Rectangle#1", "Area[Rectangle[{0,0},{2,3}]] === 6");
    audit_true("Rectangle#1", "RegionMeasure[Rectangle[{0,0},{2,3}]] === 6");
    audit_true("Rectangle#1", "Area[Rectangle[]] === 1");
    audit_true("RegionCentroid#1", "RegionCentroid[Rectangle[{0,0},{2,4}]] === {1, 2}");
    audit_true("RegionCentroid#2", "RegionCentroid[Triangle[{{0,0},{3,0},{0,3}}]] === {1, 1}");
    audit_true("RegionCentroid#2", "Perimeter[Triangle[]] === 2 + Sqrt[2]");
}

static void test_wave4_compound_polygons(void) {
    audit_true("Area#2", "Area[Polygon[{{{0,0},{1,0},{0,1}},{{2,2},{3,2},{2,3}}}]] === 1");
    audit_true("Perimeter#1", "Perimeter[Polygon[{{{0,0},{1,0},{0,1}},{{2,2},{3,2},{2,3}}}]] === 4 + 2 Sqrt[2]");
    audit_true("Area#2", "RegionCentroid[Polygon[{{{0,0},{1,0},{0,1}},{{2,2},{3,2},{2,3}}}]] === {4/3, 4/3}");
    audit_true("Area#2", "Area[Polygon[{{0,0},{4,0},{4,4},{0,4}} -> {{{1,1},{2,1},{2,2},{1,2}}}]] === 15");
}

static void test_wave4_exact_irrational(void) {
    audit_true("Area#3", "Area[Polygon[{{0,0},{Sqrt[2],0},{0,Sqrt[2]}}]] === 1");
    audit_true("Area#3", "Area[Polygon[{{0,0},{Pi,0},{0,1}}]] === Pi/2");
    audit_true("Area#3", "RegionCentroid[Polygon[{{0,0},{Pi,0},{0,1}}]] === {Pi/3, 1/3}");
    audit_true("ConvexHullRegion#1",
               "ConvexHullRegion[{{0,0},{1,0},{0,1},{Sqrt[2], 0}}] === Polygon[{{0, 0}, {Sqrt[2], 0}, {0, 1}}]");
    audit_true("RegionMember#6", "RegionMember[Polygon[{{0,0},{2,0},{2,2},{0,2}}], {Sqrt[2],1}]");
    audit_true("RegionMember#6", "!RegionMember[Polygon[{{0,0},{2,0},{2,2},{0,2}}], {Pi,1}]");
}

static void test_wave4_region_member_more(void) {
    audit_true("RegionMember#3", "RegionMember[Rectangle[], {1/2,1/2}]");
    audit_true("RegionMember#3", "RegionMember[Triangle[], {1/4,1/4}]");
    audit_true("RegionMember#3", "RegionMember[Ball[{0,0,0},1], {0,0,0}]");
    audit_true("RegionMember#3", "RegionMember[Circle[], {1,0}]");
    audit_true("RegionMember#3", "RegionMember[Point[{0,0}], {0,0}]");
    audit_true("RegionMember#3", "RegionMember[Line[{{0,0},{1,1}}], {1/2,1/2}]");
    audit_true("RegionMember#3", "!RegionMember[Line[{{0,0},{1,1}}], {1/2,1/3}]");
    audit_true("RegionMember#3", "!RegionMember[Ball[{0,0,0},1], {1,1,0}]");
    audit_true("RegionMember#3", "RegionMember[Circle[], {{1,0},{0,0}}] === {True, False}");
    audit_true("RegionMember#4",
               "!RegionMember[Polygon[{{0,0},{4,0},{4,4},{0,4}} -> {{{1,1},{2,1},{2,2},{1,2}}}], {3/2,3/2}]");
    audit_true("RegionMember#4",
               "RegionMember[Polygon[{{0,0},{4,0},{4,4},{0,4}} -> {{{1,1},{2,1},{2,2},{1,2}}}], {1,3/2}]");
    audit_true("RegionMember#4", "RegionMember[Polygon[{{0,0,0},{1,0,0},{0,1,0}}], {0,0,0}]");
    audit_true("RegionMember#4", "!RegionMember[Polygon[{{0,0,0},{1,0,0},{0,1,0}}], {0,0,1}]");
    audit_true("RegionMember#4", "RegionMember[Polygon[{{{0,0},{1,0},{0,1}},{{2,2},{3,2},{2,3}}}], {5/2,5/2}]");
}

/* ================================================================== */
/* Wave 5 G3: point sets, lower/higher-dimensional centroids, Disk     */
/* validation, the elliptical arc, machine overflow                    */
/* ================================================================== */

/* Disk / Circle are 2-D primitives with a real radius (gc_round_ok). */
static void test_wave5_disk_validation(void) {
    audit_chk("Disk#4", "Area[Disk[{0,0,0},1]]", "Area[Disk[{0, 0, 0}, 1]]");
    audit_chk("Disk#4", "Area[Disk[1,2]]", "Area[Disk[1, 2]]");
    audit_chk("Disk#4", "RegionCentroid[Disk[1,2]]", "RegionCentroid[Disk[1, 2]]");
    audit_chk("Disk#4", "Head[Area[Disk[{0,0},\"a\"]]]", "Area");
    audit_chk("Disk#4", "Area[Disk[{0,0},I]]", "Area[Disk[{0, 0}, I]]");
    audit_chk("Area#6", "Area[Disk[{0,0},I]]", "Area[Disk[{0, 0}, I]]");
    audit_chk("RegionMeasure#3", "RegionMeasure[Disk[{0,0,0},1]]", "RegionMeasure[Disk[{0, 0, 0}, 1]]");
    audit_chk("Disk#4", "{Area[Disk[{1,2},3]], RegionCentroid[Disk[{1,2},3]], Area[Disk[]]}", "{9 Pi, {1, 2}, Pi}");
}

/* Point sets and higher-dimensional regions in ArcLength / Area / RegionMeasure. */
static void test_wave5_measures(void) {
    audit_chk("ArcLength#4", "{ArcLength[Polygon[{{0,0},{1,0},{0,1}}]], ArcLength[Disk[]]}", "{Undefined, Undefined}");
    audit_chk("Point#1", "{Area[Point[{1,2}]], ArcLength[Point[{1,2}]]}", "{Undefined, Undefined}");
    audit_chk("Point#2", "RegionMeasure[Point[{{0,0},{1,1}}]]", "2");
    audit_chk("Line#3", "RegionMeasure[Line[{{0,0},{1,0}}], 1]", "1");
    audit_chk("Line#3", "RegionMeasure[Polygon[{{0,0},{1,0},{0,1}}], 2]", "1/2");
    audit_chk("ArcLength#3", "ArcLength[Circle[{0,0},{1,2},{0,Pi}]]", "4 EllipticE[3/4]");
    audit_true("ArcLength#3", "Abs[N[ArcLength[Circle[{0,0},{1,2},{0,Pi/3}]]] - NIntegrate[Sqrt[Sin[t]^2 + 4 Cos[t]^2], {t, 0, Pi/3}]] < 10^-8");
    audit_true("Area#5", "Abs[Area[Polygon[{{0,0},{1,0},{0,1}}*1.*10^200]] / (5*10^399) - 1] < 10^-12");
    audit_true("Area#5", "Abs[Area[Polygon[{{0,0},{1,0},{0,1}}*1.*10^-200]] / (5*10^-401) - 1] < 10^-12");
}

/* RegionCentroid of points, curves and solids (src/geometry_regions.c). */
static void test_wave5_centroids(void) {
    audit_chk("Point#2", "{RegionCentroid[Point[{{0,0},{2,0},{1,3}}]], RegionMember[Point[{1,2}],{1,2}], RegionMember[Point[{1,2}],{1,3}]}",
              "{{1, 1}, True, False}");
    audit_chk("Line#2", "RegionCentroid[Line[{{0,0},{2,0}}]]", "{1, 0}");
    audit_chk("RegionCentroid#4", "RegionCentroid[Line[{{0,0},{2,0}}]]", "{1, 0}");
    audit_chk("RegionCentroid#4", "RegionCentroid[Line[{{0,0},{2,0},{2,2}}]]", "{3/2, 1/2}");
    audit_chk("Circle#3", "RegionCentroid[Circle[{1,2},3]]", "{1, 2}");
    audit_chk("Circle#3", "RegionCentroid[Circle[{0,0},1,{0,Pi}]]", "{0, 2/Pi}");
    audit_chk("RegionCentroid#5", "RegionCentroid[Ball[{1,2,3},1]]", "{1, 2, 3}");
    audit_chk("RegionCentroid#5", "RegionCentroid[Cuboid[{0,0,0},{2,4,6}]]", "{1, 2, 3}");
    audit_chk("RegionCentroid#5", "RegionCentroid[Polygon[{{0,0,0},{2,0,0},{2,2,0},{0,2,0}}]]", "{1, 1, 0}");
}

static void test_wave5_hull_member(void) {
    audit_chk("ConvexHullRegion#2", "ConvexHullRegion[Polygon[{{0,0},{2,0},{1,1},{2,2},{0,2}}]]", "Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}]");
    audit_chk("RegionMember#7", "RegionMember[Polygon[{{0,0},{2,0},{1,0}}], {1,0}]", "True");
}

int main(void) {
    audit_init();
    TEST(test_region_member_point_list);
    TEST(test_region_member_disk);
    TEST(test_wave4_sector_span);
    TEST(test_wave4_disk_sector);
    TEST(test_wave4_arclength_infinite);
    TEST(test_wave4_rectangle_triangle);
    TEST(test_wave4_compound_polygons);
    TEST(test_wave4_exact_irrational);
    TEST(test_wave4_region_member_more);
    TEST(test_wave5_disk_validation);
    TEST(test_wave5_measures);
    TEST(test_wave5_centroids);
    TEST(test_wave5_hull_member);
    return audit_finish("audit_geometry");
}
