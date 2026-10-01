# format: png
# Figure: the x derivative, the y derivative (both rescaled for display) and the gradient magnitude.
shapes = Image[Table[0.2 + 0.6 Boole[(x - 160)^2 + (y - 192)^2 < 11200] + 0.3 Boole[Abs[x - 400] < 80 && Abs[y - 192] < 120], {y, 384}, {x, 576}]];
fig = ImageAssemble[{{ImageAdjust[DerivativeFilter[shapes, {0, 1}]]}, {ImageAdjust[DerivativeFilter[shapes, {1, 0}]]}, {ImageAdjust[GradientFilter[shapes]]}}]
