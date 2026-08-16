# ImageCorrelate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageCorrelate[image, kernel] correlates image with kernel: the kernel is NOT reflected, which is the only difference from ImageConvolve. The two are related exactly -- correlation equals convolution with the kernel reversed on both axes -- and they agree on any symmetric kernel, so the distinction only shows on an asymmetric one, where a delta with {{1,2,3}} gives {3,2,1} here and {1,2,3} convolved. ImageCorrelate[image, template, "NormalizedCrossCorrelation"] is template matching: it subtracts the local mean and divides by the local standard deviation, so it measures SHAPE and is invariant to brightness offset and contrast scale. Plain correlation is maximised by brightness rather than similarity -- a white patch beats a correct but darker match -- which is why raw correlation is a poor matcher. Where the template is a crop of the image the score is exactly 1 and is the global maximum. A flat window has no shape to compare and scores 0 rather than dividing by zero; scoring 1 would make every flat region match everything. Colour is reduced to luminance first for the NCC form.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (1)

```mathematica
In[1]:= ImageData[ImageCorrelate[Image[{{0., 1., 0.}, {1., 0., 1.}}], {{1.}}]]
Out[1]= {{0.0, 1.0, 0.0}, {1.0, 0.0, 1.0}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ImageConvolve](../../image-processing/ImageConvolve/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
