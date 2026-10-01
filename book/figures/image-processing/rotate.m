# format: png
# Figure: a free-angle rotation resamples bilinearly and fills the uncovered corners with black.
sq = Image[Table[0.3 + 0.6 Boole[Abs[x - 200] < 100 && Abs[y - 200] < 100], {y, 400}, {x, 400}]];
fig = ImageAssemble[{sq, ImageRotate[sq, Pi/4], ImageRotate[ImageRotate[sq, Pi/4], -Pi/4]}]
