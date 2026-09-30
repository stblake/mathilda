# format: png
# Figure: a zone plate, a standard test image built from one formula. Non-square on purpose.
zone = Image[Table[(1 + Cos[((x - 320)^2 + (y - 256)^2)/1600.])/2, {y, 512}, {x, 640}]];
ImageDimensions[zone]
fig = zone
