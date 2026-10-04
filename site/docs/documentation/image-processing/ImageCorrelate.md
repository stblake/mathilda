# ImageCorrelate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageCorrelate[image, kernel] correlates image with kernel: the kernel is NOT reflected, which is the only difference from ImageConvolve. The two are related exactly -- correlation equals convolution with the kernel reversed on both axes -- and they agree on any symmetric kernel, so the distinction only shows on an asymmetric one, where a delta with {{1,2,3}} gives {3,2,1} here and {1,2,3} convolved. The kernel may also be a single-channel Image, whose unit-scale pixels are used as the matrix. ImageCorrelate[image, template, "NormalizedCrossCorrelation"] is template matching: it subtracts the local mean and divides by the local standard deviation, so it measures SHAPE and is invariant to brightness offset and contrast scale. Plain correlation is maximised by brightness rather than similarity -- a white patch beats a correct but darker match -- which is why raw correlation is a poor matcher. Where the template is a crop of the image the score is exactly 1 and is the global maximum. A flat window has no shape to compare and scores 0 rather than dividing by zero; scoring 1 would make every flat region match everything. Colour is reduced to luminance first for the NCC form.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (1)

```mathematica
In[1]:= ImageData[ImageCorrelate[Image[{{0., 1., 0.}, {1., 0., 1.}}], {{1.}}]]
Out[1]= {{0.0, 1.0, 0.0}, {1.0, 0.0, 1.0}}
```

### Scope (1)

```mathematica
In[2]:= ImageData[ImageCorrelate[Image[{{0., 1., 0.}}], Image[{{0.5, 1.}}]]] === ImageData[ImageCorrelate[Image[{{0., 1., 0.}}], {{0.5, 1.}}]]
Out[2]= True
```

### Applications (3)

```mathematica
In[3]:= ImageData[ImageCorrelate[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], {{1, 2, 3}}]]
Out[3]= {{0.0, 0.0, 0.0}, {3.0, 2.0, 1.0}, {0.0, 0.0, 0.0}}

In[4]:= ImageData[ImageConvolve[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], {{1, 2, 3}}]]
Out[4]= {{0.0, 0.0, 0.0}, {1.0, 2.0, 3.0}, {0.0, 0.0, 0.0}}

In[5]:= ImageData[ImageCorrelate[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], {{1, 1, 1}}]]
Out[5]= {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}}
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [ImageConvolve](../../image-processing/ImageConvolve/), [Image](../../image-processing/Image/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`ImageCorrelate[image, kernel]` correlates image with kernel — the kernel is
**not** reflected, which is the only difference from `ImageConvolve`. The two are
related exactly (correlation equals convolution with the kernel reversed on both
axes) and agree on any symmetric kernel, so the distinction shows only on an
asymmetric one. The kernel may also be a single-channel `Image`.

`ImageCorrelate[image, template, "NormalizedCrossCorrelation"]` is template
matching: it subtracts the local mean and divides by the local standard
deviation, so it measures shape and is invariant to brightness offset and
contrast scale. Where the template is a crop of the image the score is exactly
`1` and is the global maximum; a flat window scores `0` rather than dividing by
zero. Plain correlation is maximised by brightness, not similarity, which is why
raw correlation is a poor matcher. Colour is reduced to luminance first for the
NCC form.
