# Flat morphology: a dot dilates into the element; the laws hold exactly.
dot = Image[ReplacePart[Table[0., {5}, {5}], {3, 3} -> 1.]];
ImageData[Dilation[dot, {{0, 1, 0}, {1, 1, 1}, {0, 1, 0}}]]
SameQ[Dilation[dot, 1], Dilation[dot, BoxMatrix[1]], Dilation[dot, 5 BoxMatrix[1]]]
SeedRandom[2]; f = RandomImage[1, {12, 9}];
Max[Abs[ImageData[Erosion[f, 1]] - (1 - ImageData[Dilation[Image[1 - ImageData[f]], 1]])]]
{e, o, c, d} = ImageData /@ {Erosion[f, 2], Opening[f, 2], Closing[f, 2], Dilation[f, 2]};
And @@ Flatten[MapThread[LessEqual, {e, o, ImageData[f], c, d}, 2]]
Opening[Opening[f, 2], 2] === Opening[f, 2]
