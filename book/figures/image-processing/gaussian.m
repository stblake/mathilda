# format: png
# Figure: GaussianFilter at radius 0, 8, 16 and 32.
zone = Image[Table[(1 + Cos[((x - 200)^2 + (y - 200)^2)/1040.])/2, {y, 400}, {x, 400}]];
fig = ImageAssemble[GaussianFilter[zone, #] & /@ {0, 8, 16, 32}]
