# format: png
# Featured: find the inner corners of a calibration chessboard; marks are dilated for visibility.
board = Image[Table[0.15 + 0.7 Mod[Quotient[x, 64] + Quotient[y, 64], 2], {y, 0, 383}, {x, 0, 511}]];
pts = ImageCorners[board, 2, 0.05, 32];
Length[pts]
marks = Image[ReplacePart[Table[0., {384}, {512}], pts -> 1.]];
fig = ImageAssemble[{board, Dilation[marks, 8]}]
