# format: png
# Featured: sharpen a soft picture with an unsharp mask, written as one convolution kernel
# (twice the identity minus a Gaussian).
zone = Image[Table[(1 + Cos[((x - 256)^2 + (y - 192)^2)/1440.])/2, {y, 384}, {x, 512}]];
soft = GaussianFilter[zone, 8];
k = 2 ReplacePart[Table[0., {25}, {25}], {13, 13} -> 1.] - GaussianMatrix[12];
fig = ImageAssemble[{soft, ImageConvolve[soft, k]}]
