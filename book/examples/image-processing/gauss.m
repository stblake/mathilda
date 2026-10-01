# The Gaussian kernel, the box kernel, and the identities between the filters.
GaussianMatrix[1]
Total[GaussianMatrix[5], 2]
BoxMatrix[1]
d = Image[ReplacePart[Table[0., {5}, {5}], {3, 3} -> 1.]];
Max[Abs[ImageData[GaussianFilter[d, 1]][[2 ;; 4, 2 ;; 4]] - GaussianMatrix[1]]]
GaussianFilter[d, 1] === ImageConvolve[d, GaussianMatrix[1]]
ImageData[MeanFilter[d, 1]][[2 ;; 4, 2 ;; 4]]
Max[Abs[ImageData[MeanFilter[d, 1]] - ImageData[ImageConvolve[d, BoxMatrix[1]/9.]]]]
