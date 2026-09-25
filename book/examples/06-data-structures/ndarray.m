# 6.4 The visible NDArray, PackedArrayQ, and what stays packed
pv = ToNDArray[{1., 2., 3.}];
c = NDArray[{1., 2., 3.}]
Head[c]
c === {1., 2., 3.}
{NDArrayQ[c], NDArrayQ[pv]}
{PackedArrayQ[c], PackedArrayQ[pv]}
Normal[c]
DataType[ToPackedArray[{1, 2, 3.}]]
ToPackedArray[{1, 2, 3.}]
big = Range[1., 1000.];
NDArrayQ[big + 1]
NDArrayQ[Sin[big]]
NDArrayQ[N[Range[1000]]]
NDArray[{{1, 2}, {3, 4}}]
Dimensions[NDArray[{{1, 2}, {3, 4}}]]
