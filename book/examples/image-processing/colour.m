# Colour images: an interleaved height x width x channels array.
rgb = Image[{{{1., 0., 0.}, {0., 1., 0.}}, {{0., 0., 1.}, {1., 1., 1.}}}];
{ImageChannels[rgb], ImageDimensions[rgb], Dimensions[ImageData[rgb]]}
ImageData[rgb][[1, 2]]
ImageData[ColorConvert[rgb, "Grayscale"]]
0.299 + 0.587 + 0.114
ImageChannels[ColorConvert[rgb, "Grayscale"]]
