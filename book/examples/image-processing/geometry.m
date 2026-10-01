# Geometry at the pixel level, on a non-square 2x3 image.
a = Image[{{1, 2, 3}, {4, 5, 6}}/10.];
ImageData[ImageReflect[a]]
ImageData[ImageReflect[a, Left]]
ImageData[ImageRotate[a, Pi]]
{ImageDimensions[a], ImageDimensions[ImageRotate[a]]}
Nest[ImageRotate, a, 4] === a
ImageData[ImagePad[a, 1]]
ImageData[ImagePad[a, 1, "Fixed"]]
ImageCrop[ImagePad[a, 2], ImageDimensions[a]] === a
ImageData[ImageCrop[Image[{{1, 1, 1, 1}, {1, 0, 0.5, 1}, {1, 1, 1, 1}}]]]
