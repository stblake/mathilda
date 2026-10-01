# The exact Euclidean distance transform: a 3-4-5 triangle comes out as exactly 5.
field = Image[ReplacePart[Table[1, {5}, {6}], {1, 1} -> 0]];
dt = ImageData[DistanceTransform[field]];
dt[[5, 4]]
dt[[2, 2]] == Sqrt[2.]
dt[[1 ;; 3, 1 ;; 4]]
disk = Image[Table[Boole[(x - 10)^2 + (y - 10)^2 < 50], {y, 19}, {x, 19}]];
Total[Flatten[ImageData[Erosion[disk, 1]]]] == Count[Flatten[ImageData[DistanceTransform[disk]]], v_ /; v >= 2]
