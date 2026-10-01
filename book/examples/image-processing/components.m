# Connected components: 8- versus 4-connectivity, and the two-pass labelling of a U.
MorphologicalComponents[Image[{{1., 0.}, {0., 1.}}]]
MorphologicalComponents[Image[{{1., 0.}, {0., 1.}}], CornerNeighbors -> False]
u = Image[{{1, 0, 1}, {1, 0, 1}, {1, 1, 1}}];
MorphologicalComponents[u, CornerNeighbors -> False]
four = Image[{{1, 0, 0, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {1, 0, 0, 1}}];
{Max[MorphologicalComponents[four]], Max[MorphologicalComponents[Dilation[four, 1]]]}
MorphologicalComponents[Image[{{0.2, 0.9}, {0.6, 0.1}}], 0.5]
