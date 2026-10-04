### Worked examples

```mathematica
(* correlation does NOT reflect: a delta with {{1,2,3}} comes back {3,2,1} *)
In[1]:= ImageData[ImageCorrelate[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], {{1, 2, 3}}]]
```

```mathematica
(* convolution, which does reflect, gives {1,2,3} on the same input *)
In[1]:= ImageData[ImageConvolve[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], {{1, 2, 3}}]]
```

```mathematica
(* a symmetric kernel makes correlation and convolution agree *)
In[1]:= ImageData[ImageCorrelate[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], {{1, 1, 1}}]]
```

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
