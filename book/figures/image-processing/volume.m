# format: png
# Figure: three slices through a ball (top), its dilation (middle) and its Gaussian blur (bottom).
ball = Image3D[Table[Boole[(x - 16)^2 + (y - 16)^2 + (z - 16)^2 < 100] + 0., {z, 32}, {y, 32}, {x, 32}]];
slices[v_] := Image[ImageData[v][[#]]] & /@ {8, 12, 16};
fig = ImageResize[ImageAssemble[{slices[ball], slices[Dilation[ball, 2]], slices[GaussianFilter[ball, 3]]}], 1296, Resampling -> "Nearest"]
