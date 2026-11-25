/* Builtin-audit regressions: packed arrays (review/issues/<Head>.md; tests/test_audit.h). */
#include "test_audit.h"

/* ================================================================== */
/* Wave 4, G2-structural (packed arrays): the typed PackedArrayQ and  */
/* ToPackedArray forms, and a single-precision buffer that cannot     */
/* hold a value.                                                      */
/* ================================================================== */

static void test_wave4_packed_arrays(void) {
    audit_chk("PackedArrayQ#1", "PackedArrayQ[ToPackedArray[{1,2,3}], Integer]", "True");
    audit_chk("PackedArrayQ#1", "PackedArrayQ[ToPackedArray[{1,2,3}], Real]", "False");
    audit_chk("PackedArrayQ#1", "PackedArrayQ[ToPackedArray[{{1,2},{3,4}}], Integer, 2]", "True");
    audit_chk("PackedArrayQ#1", "PackedArrayQ[ToPackedArray[{{1,2},{3,4}}], Integer, 1]", "False");
    audit_chk("PackedArrayQ#1", "PackedArrayQ[{1, 2, 3}, Integer]", "False");
    audit_chk("ToPackedArray#1", "ToPackedArray[{1, 2, 3}, Real]", "{1.0, 2.0, 3.0}");
    audit_chk("ToPackedArray#1", "PackedArrayQ[ToPackedArray[{1, 2, 3}, Real], Real]", "True");
    audit_chk("ToPackedArray#1", "PackedArrayQ[ToPackedArray[{1, 2}, Complex], Complex]", "True");
    audit_chk("ToPackedArray#1", "PackedArrayQ[ToPackedArray[{1.5, 2}, Integer]]", "False");
    audit_chk("ToNDArray#1", "ToNDArray[{1., 1.*^300}, DataType -> \"float32\"]", "{1.0, 1.0*^300}");
    audit_chk("ToNDArray#1", "NDArrayQ[ToNDArray[{1., 1.*^300}, DataType -> \"float32\"]]", "False");
    audit_chk("ToNDArray#1", "DataType[ToNDArray[{1., 2.}, DataType -> \"float32\"]]", "\"float32\"");
}

int main(void) {
    audit_init();
    TEST(test_wave4_packed_arrays);
    return audit_finish("audit_packed_arrays");
}
