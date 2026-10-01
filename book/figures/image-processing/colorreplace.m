# format: png
# Figure: ColorReplace swaps one exact colour and leaves the rest alone.
dots = Image[Table[Which[(x - 128)^2 + (y - 128)^2 < 8000, {0.9, 0.1, 0.1}, (x - 384)^2 + (y - 128)^2 < 8000, {0.1, 0.2, 0.9}, True, {0.95, 0.95, 0.9}], {y, 256}, {x, 512}]];
fig = ImageAssemble[{{dots}, {ColorReplace[dots, RGBColor[0.9, 0.1, 0.1] -> RGBColor[1, 0.8, 0]]}}]
