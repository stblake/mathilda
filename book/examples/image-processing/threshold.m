# Otsu's threshold, the strict comparison in Binarize, and a local threshold.
two = Image[{{0.1, 0.2, 0.15}, {0.8, 0.9, 0.85}}];
FindThreshold[two]
ImageData[Binarize[two], "Bit"]
ImageData[Binarize[Image[{{0.3, 0.5, 0.7}}], 0.5], "Bit"]
ImageType[Binarize[two]]
lit = Image[Table[(0.2 + 0.6 x/32) (0.6 + 0.4 Mod[Quotient[x - 1, 4] + Quotient[y - 1, 4], 2]), {y, 16}, {x, 32}]];
Union[Flatten[ImageData[Binarize[lit]][[All, 1 ;; 10]]]]
Union[Flatten[ImageData[LocalAdaptiveBinarize[lit, 4]][[All, 1 ;; 10]]]]
