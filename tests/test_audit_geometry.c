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

int main(void) {
    audit_init();
    TEST(test_region_member_point_list);
    TEST(test_region_member_disk);
    return audit_finish("audit_geometry");
}
