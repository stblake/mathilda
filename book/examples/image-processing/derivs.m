# Derivative stencils on a ramp and a step, and Canny's one-pixel edge.
ramp = Image[Table[0.1 x, {y, 3}, {x, 6}]];
ImageData[DerivativeFilter[ramp, {0, 1}]][[2]]
ImageData[DerivativeFilter[ramp, {1, 0}]][[2]]
step = Image[Table[If[x > 4, 1., 0.], {y, 3}, {x, 8}]];
ImageData[DerivativeFilter[step, {0, 1}]][[2]]
ImageData[GradientFilter[step]][[2]]
ImageData[EdgeDetect[step, 0], "Bit"]
diag = Image[Table[If[x > y, 1., 0.], {y, 8}, {x, 8}]];
{dx, dy} = ImageData[DerivativeFilter[diag, #]][[4, 4]] & /@ {{0, 1}, {1, 0}}
{ImageData[GradientFilter[diag]][[4, 4]], Sqrt[dx^2 + dy^2], Abs[dx] + Abs[dy]}
