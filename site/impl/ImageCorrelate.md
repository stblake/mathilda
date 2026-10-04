---
source: src/imagefilter.c
---
**Algorithm.** `builtin_imagecorrelate` correlates an image with a kernel. Plain
correlation (`correlate_planes`) is identical to `ImageConvolve` except the
kernel is **not** reflected — correlation equals convolution with the kernel
reversed on both axes, so the two differ only on an asymmetric kernel, where a
delta image with `{{1,2,3}}` gives `{3,2,1}` here and `{1,2,3}` convolved. The
kernel may be a numeric matrix (`ker_load`) or a **single-channel `Image`**,
whose unit-scale pixels are used as the matrix (a multichannel template is
declined, since channel-against-channel vs. luminance is a choice the matrix form
never had to make). The third-argument form
`ImageCorrelate[image, template, "NormalizedCrossCorrelation"]` is template
matching (`ncc_planes`): at each position it subtracts the local mean and divides
by the local standard deviation, so it measures **shape** and is invariant to
brightness offset and contrast scale. Where the template is a crop of the image
the NCC score is exactly `1` and is the global maximum; a flat window has no
shape and scores `0` rather than dividing by zero. For the NCC form colour is
reduced to Rec. 601 luminance first (`img_grey_plane`), the same choice
`GradientFilter` makes. Padding and the interior/border split are inherited from
the shared convolution core.

**Data structures.** Flat `double` buffers: plain correlation keeps all channels
(`image_load` → height · width · channels), NCC works on a single luminance plane.
The kernel/template is a `kh × kw` `double` matrix. The result is a `"Real"`
image (`image_build_real`).

**Complexity / limits.** `O(width · height · channels · kw · kh)` for plain
correlation; the NCC pass adds per-window mean and variance (a constant number of
sums per tap). Plain correlation is maximised by brightness rather than
similarity — a white patch beats a correct but darker match — which is why the
NCC form exists for matching.
