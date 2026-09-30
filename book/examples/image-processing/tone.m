# Point operations: stretching, the contrast/brightness/gamma curve, and the histogram.
low = Image[{{0.40, 0.45}, {0.50, 0.55}}];
ImageData[ImageAdjust[low]]
ImageAdjust[ImageAdjust[low]] === ImageAdjust[low]
ImageData[ImageAdjust[Image[{{0.25, 0.5, 0.75}}], {1, 0}]]
ImageData[ImageAdjust[Image[{{0.25, 0.5, 0.75}}], {0, 0.1}]]
ImageData[ImageAdjust[Image[{{0.25, 0.5, 0.75}}], {0, 0, 2}]]
ImageLevels[Image[{{0, 1, 1}, {1, 0, 1}}]]
ImageLevels[Image[{{0.1, 0.2, 0.9, 0.95}}], 4]
ImageData[HistogramTransform[Image[{{0.40, 0.42}, {0.44, 0.60}}]]]
