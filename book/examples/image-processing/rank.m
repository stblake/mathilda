# The median removes an outlier exactly; the mean smears it. The median is not separable.
spike = Image[{{0.2, 0.2, 0.2}, {0.2, 1., 0.2}, {0.2, 0.2, 0.2}}];
ImageData[MedianFilter[spike, 1]]
ImageData[MeanFilter[spike, 1]]
w = {{1, 2, 9}, {3, 4, 5}, {6, 7, 8}};
{Median[Flatten[w]], Median[Median /@ w]}
ImageData[MedianFilter[Image[w/10.], 1]][[2, 2]]
