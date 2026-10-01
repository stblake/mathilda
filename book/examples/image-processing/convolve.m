# Convolution reflects the kernel; correlation does not. Padding replicates the edge.
delta = Image[{{0., 1., 0.}}];
ImageData[ImageConvolve[delta, {{1, 2, 3}}]]
ImageData[ImageCorrelate[delta, {{1, 2, 3}}]]
ImageData[ImageConvolve[Image[{{0.3, 0.3, 0.3}, {0.3, 0.3, 0.3}}], {{0.25, 0.5, 0.25}}]]
ImageData[ImageConvolve[Image[{{0.3, 0.3, 0.3}}], {{1, 1, 1}}]]
ImageType[ImageConvolve[Image[{{0, 128, 255}}], {{1}}]]
