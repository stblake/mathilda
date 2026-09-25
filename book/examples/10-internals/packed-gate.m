# A packed array and the list it packs are indistinguishable -- bar NDArrayQ.
v = Range[5];
ListQ[v]
NDArrayQ[v]
Head[v]
# A visible NDArray is a different, atomic object:
w = NDArray[{1., 2., 3.}];
Head[w]
AtomQ[w]
ListQ[w]
NDArrayQ[w]
# The transparency gate keeps an integer sum exact even over a packed buffer:
Total[Range[5]]
Head[Total[Range[5]]]
