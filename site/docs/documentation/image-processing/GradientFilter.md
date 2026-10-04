# GradientFilter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GradientFilter[image] gives the gradient magnitude Sqrt[dx^2 + dy^2], using the normalised Sobel derivatives of DerivativeFilter. The magnitude rather than |dx| + |dy| because it is ROTATION INVARIANT: an edge at 45 degrees reports the same strength as one at 0, where the absolute sum would report it sqrt(2) times stronger and so bias every downstream threshold by orientation. A colour image is reduced to luminance first and differentiated once, rather than differentiated per channel and combined by some arbitrary rule.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (1)

```mathematica
In[1]:= ImageData[GradientFilter[Image[{{0., 0., 1., 1.}, {0., 0., 1., 1.}}]]]
Out[1]= {{0.0, 0.5, 0.5, 0.0}, {0.0, 0.5, 0.5, 0.0}}
```

### Applications (3)

A vertical step edge

```mathematica
In[2]:= edge = Image[{{0., 0., 1., 1.}, {0., 0., 1., 1.}, {0., 0., 1., 1.}}];
```

The response peaks across the step

```mathematica
In[3]:= Part[ImageData[GradientFilter[edge]], 2, 2]
Out[3]= 0.5
```

Same size as the input, single channel

```mathematica
In[4]:= ImageDimensions[GradientFilter[edge]]
Out[4]= {4, 3}
```

## Implementation notes

**Algorithm.** `builtin_gradientfilter` returns the gradient **magnitude**
`Sqrt[dx² + dy²]`. Colour is reduced to luminance **first** (`img_grey_plane`), then
differentiated once — not differentiated per channel and combined, which would need an
arbitrary rule (max? sum? norm?). The two first derivatives come from the same normalised Sobel
stencils `DerivativeFilter` uses (`deriv_kernel(0,1)` and `deriv_kernel(1,0)`, each through
`convolve_dispatch`), and are combined pixelwise as the Euclidean norm.

The magnitude rather than `|dx| + |dy|` is chosen for **rotation invariance**: the length of a
vector does not depend on which way the axes point, so an edge at 45° reports the same strength
as one at 0°, where the absolute-sum alternative would report it `√2` times stronger and bias
every downstream threshold by orientation. A volume (`gradient3_run`) takes three separable
derivatives and combines them with one square root at the end — `sqrt` is not additive, so
taking it per component would give the sum of absolute derivatives, a different and larger
quantity.

**Data structures.** A grey buffer plus two derivative buffers `dx`, `dy`; the result through
`image_build_real` as a packed single-channel `"Real"` image (every image head returns a packed
buffer — `make check-image-packing`).

**Complexity / limits.** `O(pixels)` — two separable 3-tap convolutions plus a per-pixel
`hypot`. The border is replicate, as for all convolutions; the output is single-channel.

**Attributes:** `Protected`.

## References

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`GradientFilter[image]` gives the gradient magnitude `Sqrt[dx^2 + dy^2]`, using the normalised
Sobel derivatives of `DerivativeFilter`. The magnitude rather than `|dx| + |dy|` because it is
**rotation invariant**: an edge at 45° reports the same strength as one at 0°, where the
absolute sum would report it `√2` times stronger and so bias every downstream threshold by
orientation.

A colour image is reduced to luminance first and differentiated once, rather than
differentiated per channel and combined by some arbitrary rule. The result is a single-channel
`"Real"` image.
