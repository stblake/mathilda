/* Builtin-audit regressions: image processing (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* RandomImage#1 (critical): w*h*ch overflowed (2^40 x 2^40 is 0 mod 2^64), so
 * the buffer was under-allocated (SIGBUS), or the request ran for minutes.
 * Checked, an impossible size aborts at once with General::nomem. */
static void test_random_image_huge_size(void) {
    audit_chk("RandomImage#1", "Quiet[CheckAbort[RandomImage[1, {2^40, 2^40}], aborted]]", "aborted");
    audit_chk("RandomImage#1",
              "Quiet[Check[CheckAbort[RandomImage[1, {4611686018427387904, 4}], 0], fired, General::nomem]]",
              "fired");
    audit_chk("RandomImage#1",
              "Quiet[CheckAbort[RandomImage[1, {2^31, 2^31}, ColorSpace -> \"RGB\"], aborted]]", "aborted");
    audit_chk("RandomImage#1", "ImageDimensions[RandomImage[1, {4, 3}]]", "{4, 3}");
    audit_chk("RandomImage#1", "ImageDimensions[RandomImage[2, 5]]", "{5, 5}");
}

int main(void) {
    audit_init();
    TEST(test_random_image_huge_size);
    return audit_finish("audit_image_processing");
}
