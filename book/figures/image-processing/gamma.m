# format: png
# Figure: a grey ramp under ImageAdjust gamma 1/2, 1 and 2 (top to bottom), then contrast +1.
ramp = Image[Table[x/511., {y, 48}, {x, 512}]];
fig = ImageAssemble[{{ImageAdjust[ramp, {0, 0, 0.5}]}, {ramp}, {ImageAdjust[ramp, {0, 0, 2}]}, {ImageAdjust[ramp, {1, 0}]}}]
