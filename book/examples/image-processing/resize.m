# Resampling: area averaging, point sampling and bilinear enlargement.
chk = Image[Table[Mod[x + y, 2] + 0., {y, 4}, {x, 4}]];
ImageData[ImageResize[chk, {2, 2}]]
ImageData[ImageResize[chk, {2, 2}, Resampling -> "Nearest"]]
a = Image[{{1, 2, 3}, {4, 5, 6}}/10.];
ImageData[ImageResize[a, {6, 4}, Resampling -> "Nearest"]]
ImageData[ImageResize[a, {6, 4}]]
ImageDimensions[ImageResize[a, 12]]
SeedRandom[1]; big = RandomImage[1, {60, 40}];
{Mean[Flatten[ImageData[big]]], Mean[Flatten[ImageData[ImageResize[big, {40, 25}]]]]}
