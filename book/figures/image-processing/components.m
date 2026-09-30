# format: png
# Figure: connected components of a random blob pattern. Each label is drawn as its own grey
# level; background stays black.
SeedRandom[11]; spots = Binarize[GaussianFilter[RandomImage[1, {128, 96}], 5], 0.53];
labels = MorphologicalComponents[spots];
Max[labels]
fig = ImageResize[ImageAssemble[{spots, Image[Sign[labels] (0.3 + 0.7 labels/Max[labels])]}], 1024, Resampling -> "Nearest"]
