# format: png
# Figure: shrinking a zone plate 4x by point sampling and by area averaging, then enlarging both
# back with "Nearest" so the pixels stay visible.
zone = Image[Table[(1 + Cos[((x - 256)^2 + (y - 256)^2)/1600.])/2, {y, 512}, {x, 512}]];
up[i_] := ImageResize[i, 512, Resampling -> "Nearest"];
fig = ImageAssemble[{zone, up[ImageResize[zone, 128, Resampling -> "Nearest"]], up[ImageResize[zone, 128]]}]
