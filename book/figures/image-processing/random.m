# format: png
# Figure: uniform noise from RandomImage, grey and colour, reproducible via SeedRandom.
SeedRandom[7];
fig = ImageAssemble[{RandomImage[1, {400, 400}], RandomImage[1, {400, 400}, ColorSpace -> "RGB"]}]
