# 6.4 The transparency gate: packing never changes an answer
p = Range[1., 6.];
NDArrayQ[p]
Cases[p, x_ /; x > 3.]
Reverse[p]
NDArrayQ[Reverse[p]]
Total[p] == Total[FromPackedArray[p]]
