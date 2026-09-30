# Composition: centring, placement from the bottom left, opacity, and tiling.
base = Image[Table[0., {4}, {6}]];
patch = Image[{{1., 1.}, {1., 1.}}];
ImageData[ImageCompose[base, patch]]
ImageData[ImageCompose[base, {patch, 0.5}, {1, 1}]]
ImageChannels[ImageCompose[base, Image[{{{1., 0., 0.}}}]]]
ImageDimensions[ImageAssemble[{{patch, patch}, {patch, Image[{{1.}}]}}]]
