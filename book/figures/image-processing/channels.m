# format: png
# Figure: an RGB image and its red and green channels shown as grey images.
col = Image[Table[{x, y, 256}/512., {y, 512}, {x, 512}]];
{red, green} = Image[ImageData[col][[All, All, #]]] & /@ {1, 2};
fig = ImageAssemble[{col, red, green}]
