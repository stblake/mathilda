# format: png
# Featured: find a patch in a texture, and mark it. NCC locates the 11x11 crop; ImageCompose
# draws a red frame there, converting {row, column} to {x, y from the bottom}.
SeedRandom[4]; scene = GaussianFilter[RandomImage[1, {160, 96}], 1];
t = ImageData[scene][[40 ;; 50, 111 ;; 121]];
s = ImageData[ImageCorrelate[scene, t, "NormalizedCrossCorrelation"]];
{row, col} = First[Position[s, Max[s]]]
box = SetAlphaChannel[Image[Table[{1., 0.1, 0.1}, {15}, {15}]], Image[Table[Boole[Max[Abs[i - 8], Abs[j - 8]] >= 6] + 0., {i, 15}, {j, 15}]]];
fig = ImageResize[ImageCompose[scene, box, {col, 96 - row + 1}], 640, Resampling -> "Nearest"]
