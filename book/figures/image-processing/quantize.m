# format: png
# Figure: median-cut colour quantisation to 2, 4 and 16 colours.
col = Image[Table[{x, y, 256 + (x - y)/2}/512., {y, 512}, {x, 512}]];
fig = ImageAssemble[{col, ColorQuantize[col, 2], ColorQuantize[col, 4], ColorQuantize[col, 16]}]
