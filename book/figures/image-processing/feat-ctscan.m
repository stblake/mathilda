# format: png
# Featured: look inside a scanned volume. A hollow shell hides a small ball; one slice shows
# both, and the mean along the depth axis gives a simulated projection radiograph.
ct = Image3D[Table[Boole[196 < (x - 24)^2 + (y - 24)^2 + (z - 24)^2 < 400 || (x - 28)^2 + (y - 20)^2 + (z - 24)^2 < 16] + 0., {z, 48}, {y, 48}, {x, 48}]];
ImageDimensions[ct]
fig = ImageResize[ImageAssemble[{Image[ImageData[ct][[24]]], ImageAdjust[Image[Mean[ImageData[ct]]]]}], 864, Resampling -> "Nearest"]
