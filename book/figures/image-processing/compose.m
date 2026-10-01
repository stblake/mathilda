# format: png
# Figure: a half-transparent disk composited three times onto a colour gradient.
bg = Image[Table[{x/768., 0.3, 1 - x/768.}, {y, 384}, {x, 768}]];
disk = SetAlphaChannel[Image[Table[{1., 1., 0.3}, {192}, {192}]], Image[Table[Boole[(x - 96)^2 + (y - 96)^2 < 8000] 0.7, {y, 192}, {x, 192}]]];
fig = ImageCompose[ImageCompose[ImageCompose[bg, disk, {160, 192}], disk, {384, 240}], disk, {608, 144}]
