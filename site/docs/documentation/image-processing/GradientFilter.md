# GradientFilter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GradientFilter[image] gives the gradient magnitude Sqrt[dx^2 + dy^2], using the normalised Sobel derivatives of DerivativeFilter. The magnitude rather than |dx| + |dy| because it is ROTATION INVARIANT: an edge at 45 degrees reports the same strength as one at 0, where the absolute sum would report it sqrt(2) times stronger and so bias every downstream threshold by orientation. A colour image is reduced to luminance first and differentiated once, rather than differentiated per channel and combined by some arbitrary rule.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (1)

```mathematica
In[1]:= ImageData[GradientFilter[Image[{{0., 0., 1., 1.}, {0., 0., 1., 1.}}]]]
Out[1]= {{0.0, 0.5, 0.5, 0.0}, {0.0, 0.5, 0.5, 0.0}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
