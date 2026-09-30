# format: png
# Figure: the three padding rules of ImagePad on an off-centre disk.
d = Image[Table[0.2 + 0.7 Boole[(x - 80)^2 + (y - 120)^2 < 4800], {y, 256}, {x, 256}]];
fig = ImageAssemble[{ImagePad[d, 96], ImagePad[d, 96, "Fixed"], ImagePad[d, 96, "Reflected"]}]
