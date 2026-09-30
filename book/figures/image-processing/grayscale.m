# format: png
# Figure: three pure primaries and their Rec. 601 luminances.
bars = Image[Table[Which[x <= 160, {1., 0, 0}, x <= 320, {0, 1., 0}, True, {0, 0, 1.}], {y, 240}, {x, 480}]];
fig = ImageAssemble[{bars, ColorConvert[bars, "Grayscale"]}]
