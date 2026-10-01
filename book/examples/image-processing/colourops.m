# Colour operations: quantisation by median cut, replacement by distance, and the alpha channel.
ImageData[ColorQuantize[Image[{{0.1, 0.12, 0.8, 0.82}}], 2]]
rgb = Image[{{{1., 0., 0.}, {0.97, 0.03, 0.}}, {{0., 0., 1.}, {1., 1., 1.}}}];
ImageData[ColorReplace[rgb, RGBColor[1, 0, 0] -> RGBColor[0, 1, 0]]][[1]]
ImageData[ColorReplace[rgb, RGBColor[1, 0, 0] -> RGBColor[0, 1, 0], 0.05]][[1]]
ImageData[AlphaChannel[Image[{{0.2, 0.7}}]]]
half = SetAlphaChannel[Image[{{1., 1.}}], 0.5];
{ImageChannels[half], ImageData[half]}
ImageData[RemoveAlphaChannel[half, 0.]]
ImageData[RemoveAlphaChannel[half]]
