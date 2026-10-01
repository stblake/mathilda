# format: png
# Featured: skeletonize a pen stroke, then prune the spur. The dot of the "i" thins to a point.
stroke = Image[Table[Boole[Abs[y - 48 - 25 Sin[x/18.]] < 7 && 12 < x < 180 || (x - 150)^2 + (y - 30)^2 < 60], {y, 96}, {x, 192}]];
skel = Thinning[stroke];
fig = ImageResize[ImageAssemble[{{stroke}, {skel}, {Pruning[skel, 8]}}], 768, Resampling -> "Nearest"]
