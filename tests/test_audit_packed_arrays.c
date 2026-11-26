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

/* Wave 5 G2: DataType options, deep FromPackedArray, docstrings. */
static void test_wave5_packed_batch(void) {
    audit_true("DataType#1", "StringContainsQ[Information[DataType], \"int64\"] && "
                             "StringContainsQ[Information[DataType], \"bool\"]");
    audit_chk("FromPackedArray#1", "FromPackedArray[{Range[5], NDArray[{1, 2}]}]",
              "{{1, 2, 3, 4, 5}, {1.0, 2.0}}");
    audit_chk("FromPackedArray#1", "FromPackedArray[f[NDArray[{1, 2}]]]", "f[{1.0, 2.0}]");
    audit_true("NDArrayQ#2", "StringContainsQ[Information[NDArrayQ], \"packed List\"] && "
                             "NDArrayQ[ToNDArray[{1., 2., 3.}]]");
    audit_chk("ToNDArray#2", "Quiet[Check[ToNDArray[{1, 2}, DataType -> \"bogus\"], fired, ToNDArray::dtype]]",
              "fired");
    audit_chk("ToNDArray#2", "Quiet[Check[ToNDArray[{1, 2}, DataType -> \"int32\"], fired, ToNDArray::dtype]]",
              "fired");
    audit_chk("ToNDArray#2", "Quiet[ToNDArray[{300, 2}, DataType -> \"int8\"]]",
              "ToNDArray[{300, 2}, DataType -> \"int8\"]");
    audit_chk("ToNDArray#3", "DataType[ToNDArray[{1, 2}, DataType -> \"int64\", DataType -> \"float64\"]]",
              "\"int64\"");
    audit_chk("ToNDArray#6", "NDArrayQ[ToNDArray[{1, 2.5}]]", "True");
    audit_chk("ToPackedArray#2", "NDArrayQ[ToNDArray[{1, 2.5}]]", "True");
}

/* Wave 6 G4: ToNDArray / ToPackedArray widen a mixed Integer/Real list to
 * float64 and pack a machine complex list (the spec now says so). */
static void test_wave6_g4_packed(void) {
    audit_chk("NDArrayQ#1", "NDArrayQ[ToNDArray[{1, 2.5}]]", "True");
    audit_chk("NDArrayQ#1", "Head[ToNDArray[{1, 2.5}][[1]]]", "Real");
    audit_chk("NDArrayQ#1", "NDArrayQ[{1, 2.5, 3, 4}]", "False");
    audit_chk("ToPackedArray#2", "{NDArrayQ[ToNDArray[{1. + I, 2.}]], DataType[ToNDArray[{1. + I, 2.}]]}",
              "{True, \"complex64\"}");
    audit_chk("ToPackedArray#2", "NDArrayQ[ToNDArray[{1 + I, 2}]]", "False");
}

int main(void) {
    audit_init();
    TEST(test_wave4_packed_arrays);
    TEST(test_wave5_packed_batch);
    TEST(test_wave6_g4_packed);
    return audit_finish("audit_packed_arrays");
}
