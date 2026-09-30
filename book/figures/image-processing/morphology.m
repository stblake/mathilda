# format: png
# Figure: a speckled shape; top row erosion and dilation, bottom row opening and closing.
# Enlarged with "Nearest" so single pixels stay visible.
SeedRandom[9]; blob = Image[Table[Boole[(x - 48)^2 + (y - 48)^2 < 900 && RandomReal[] > 0.08 || RandomReal[] < 0.03], {y, 96}, {x, 96}]];
fig = ImageResize[ImageAssemble[{{blob, Erosion[blob, 1], Dilation[blob, 1]}, {blob, Opening[blob, 1], Closing[blob, 1]}}], 1296, Resampling -> "Nearest"]
