# format: png
# Featured: make a six-colour poster from a smooth colour field by median cut.
sunset = Image[Table[{0.9 - 0.4 y/384, 0.3 + 0.4 Sin[x/120.]^2, 0.2 + 0.7 y/384}, {y, 384}, {x, 576}]];
fig = ImageAssemble[{sunset, ColorQuantize[sunset, 6]}]
