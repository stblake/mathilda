# Images as expressions: a 2x3 byte image, its canonical form and its accessors.
img = Image[{{0, 128, 255}, {64, 192, 32}}];
ImageQ[img]
{ImageDimensions[img], Dimensions[ImageData[img]]}
ImageType[img]
ImageData[img, "Byte"]
ImageData[img]
FullForm[img]
