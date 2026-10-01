# Volumes: Image3D stores slices outermost and reports its dimensions fully reversed.
vol = Image3D[Table[x + 10 y + 100 z, {z, 2}, {y, 3}, {x, 4}]/1000.];
{ImageDimensions[vol], Dimensions[ImageData[vol]]}
{ImageQ[vol], Image3DQ[vol]}
ImageData[vol][[2, 3, 4]]
ball = Image3D[Table[Boole[(x - 8)^2 + (y - 8)^2 + (z - 8)^2 < 25] + 0., {z, 15}, {y, 15}, {x, 15}]];
Total[ImageData[#], 3] & /@ {ball, GaussianFilter[ball, 1], Erosion[ball, 1], Dilation[ball, 1]}
Max[ImageData[DistanceTransform[ball]]]
ImageDimensions[ImagePad[ball, 2]]
ImageData[ImageReflect[vol, Front]][[1, 1]]
