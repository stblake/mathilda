# format: png
# Figure: an asymmetric ramp triangle, its two reflections, and its half turn.
t = Image[Table[If[x < y, x/400., 0.], {y, 400}, {x, 400}]];
fig = ImageAssemble[ImagePad[#, 8, 1.] & /@ {t, ImageReflect[t], ImageReflect[t, Left], ImageRotate[t, Pi]}]
