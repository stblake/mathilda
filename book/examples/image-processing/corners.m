# Corner detection: flat regions and straight edges score exactly zero; corners do not.
sq = Image[Table[If[5 <= y <= 12 && 5 <= x <= 12, 1., 0.], {y, 16}, {x, 16}]];
ImageCorners[sq]
edge = Image[Table[If[x > 8, 1., 0.], {y, 16}, {x, 16}]];
Max[ImageData[CornerFilter[edge]]]
Min[ImageData[CornerFilter[edge, 2, "Harris"]]] < 0
chk = Image[Table[If[Mod[Quotient[y - 1, 6] + Quotient[x - 1, 6], 2] == 0, 0., 1.], {y, 24}, {x, 24}]];
Length[ImageCorners[chk]]
Max[Abs[ImageData[CornerFilter[ImageRotate[chk]]] - ImageData[ImageRotate[CornerFilter[chk]]]]]
