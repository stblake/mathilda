# ImageConvolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageConvolve[image, kernel] convolves image with the rank-2 numeric kernel. This is true convolution: the kernel is REFLECTED before summing, so it differs from correlation on an asymmetric kernel (the two agree exactly on a symmetric one such as a Gaussian or a box). Out-of-range reads clamp to the nearest edge pixel, replicating the border, so a constant image convolved with a kernel summing to 1 comes back unchanged everywhere including the edges -- zero padding would darken them. The result is always a "Real" image of the same dimensions, since a filtered byte is not generally a byte. Each colour channel is convolved independently.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= ImageData[ImageConvolve[Image[{{0., 1., 0.}}], {{1, 2, 3}}]]
Out[1]= {{1.0, 2.0, 3.0}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ImageCorrelate](../../image-processing/ImageCorrelate/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
