# Import and Export: a round trip through a PNG file, and reproducible noise.
f = Export["/tmp/mathilda-tiny.png", Image[{{0, 128, 255}}]]
back = Import[f];
{ImageType[back], ImageData[back]}
SeedRandom[42]; a = RandomImage[1, {4, 3}];
ImageData[a]
SeedRandom[42]; ImageData[RandomImage[1, {4, 3}]] === ImageData[a]
ImageDimensions[RandomImage[]]
