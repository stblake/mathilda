# Packed arrays

4 built-in function(s) in this category.

- [`FromNDArray`](FromNDArray.md) — FromNDArray[expr] returns expr with any dense buffer storage undone: a packed List becomes an ordinary List of separate elements, and an NDArray[...] becomes the nested List of its entries. Anything else is returned unchanged. Inverse of ToNDArray.  _(Stable)_
- [`FromPackedArray`](FromPackedArray.md) — FromPackedArray[expr] is FromNDArray[expr]: it returns expr with any dense buffer storage undone, so a packed List becomes an ordinary List of separate elements and an NDArray[...] becomes the nested List of its entries. Provided under Mathematica's name for the same operation; see FromNDArray.  _(Stable)_
- [`ToNDArray`](ToNDArray.md) — ToNDArray[list] returns list stored as a dense machine-precision buffer. The result is still a List -- same Head, same printed form, same elements -- but NDArrayQ gives True for it. ToNDArray[list, DataType -> "float64"] forces the element type. Returns list unchanged when it is not rectangular, is empty, or holds anything other than uniformly Integer or uniformly Real machine values. Unlike automatic packing it ignores the size threshold.  _(Stable)_
- [`ToPackedArray`](ToPackedArray.md) — ToPackedArray[list] is ToNDArray[list]: it returns list stored as a dense machine-precision buffer. The result is still a List -- same Head, same printed form, same elements -- but NDArrayQ gives True for it. Provided under Mathematica's name for the same operation; see ToNDArray.  _(Stable)_
