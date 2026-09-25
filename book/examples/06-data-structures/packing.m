# 6.4 Packed arrays: invisible representation, exact contract, automatic packing
pv = ToNDArray[{1., 2., 3.}]
{Head[pv], ListQ[pv], AtomQ[pv]}
NDArrayQ[pv]
pv === {1., 2., 3.}
DataType[pv]
iv = ToNDArray[{1, 2, 3}];
DataType[iv]
Head[iv[[2]]]
Total[Range[1000]]
Range[10]/2
NDArrayQ[ToNDArray[{1, 1/2, 3}]]
NDArrayQ[Range[1., 300.]]
NDArrayQ[{1., 2., 3., 4.}]
$AutoArrayPacking = False; NDArrayQ[Range[1., 300.]]
$AutoArrayPacking = True; NDArrayQ[Range[1., 300.]]
FromPackedArray[pv]
NDArrayQ[FromPackedArray[pv]]
