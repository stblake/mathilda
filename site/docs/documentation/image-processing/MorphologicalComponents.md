# MorphologicalComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MorphologicalComponents[image] labels the connected components of the foreground, giving an INTEGER MATRIX with background 0 and components numbered 1..k in raster order of first appearance. MorphologicalComponents[image, t] takes pixels above t as foreground (default 0, so nonzero is foreground). CornerNeighbors -> False uses 4-connectivity instead of the default 8. Two pixels touching only at a corner are ONE component under 8 and TWO under 4, which is the property that distinguishes the two rules -- every other property holds under either. A matrix rather than an Image, deliberately: Image type inference would call a label array of 1..12 a "Byte" image and ImageData would then divide every label by 255. Labels are indices, not brightnesses. Contiguous labels in scan order mean Max of the result is the component count.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= MorphologicalComponents[Image[{{1., 0.}, {0., 1.}}]]
Out[1]= {{1, 0}, {0, 1}}
```

### Options (1)

```mathematica
In[2]:= MorphologicalComponents[Image[{{1., 0.}, {0., 1.}}], CornerNeighbors -> False]
Out[2]= {{1, 0}, {0, 2}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Image](../../image-processing/Image/), [ImageData](../../image-processing/ImageData/), [Max](../../data-structures/Max/), [List](../../other-advanced/List/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
