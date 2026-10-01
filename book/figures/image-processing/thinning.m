# format: png
# Figure: a thick ring with a spur, its Zhang-Suen skeleton, and the skeleton after Pruning.
ring = Image[Table[Boole[400 < (x - 48)^2 + (y - 48)^2 < 1300 || Abs[y - 48] < 4 && 60 < x < 90], {y, 96}, {x, 96}]];
skel = Thinning[ring];
fig = ImageResize[ImageAssemble[{ring, skel, Pruning[skel, 10]}], 1296, Resampling -> "Nearest"]
