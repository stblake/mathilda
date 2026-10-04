# MorphologicalComponents

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MorphologicalComponents[image] labels the connected components of the foreground, giving an INTEGER MATRIX with background 0 and components numbered 1..k in raster order of first appearance. MorphologicalComponents[image, t] takes pixels above t as foreground (default 0, so nonzero is foreground). CornerNeighbors -> False uses 4-connectivity instead of the default 8. Two pixels touching only at a corner are ONE component under 8 and TWO under 4, which is the property that distinguishes the two rules -- every other property holds under either. A matrix rather than an Image, deliberately: Image type inference would call a label array of 1..12 a "Byte" image and ImageData would then divide every label by 255. Labels are indices, not brightnesses. Contiguous labels in scan order mean Max of the result is the component count.`**

## Examples (6)

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

### Applications (4)

```mathematica
In[3]:= MorphologicalComponents[Image[{{1, 0, 1}, {0, 0, 0}, {1, 0, 1}}]]
Out[3]= {{1, 0, 2}, {0, 0, 0}, {3, 0, 4}}

In[4]:= MorphologicalComponents[Image[{{1, 0}, {0, 1}}]]
Out[4]= {{1, 0}, {0, 1}}

In[5]:= MorphologicalComponents[Image[{{1, 0}, {0, 1}}], 0, CornerNeighbors -> False]
Out[5]= {{1, 0}, {0, 2}}

In[6]:= Max[MorphologicalComponents[Image[{{1, 0, 1}, {1, 0, 1}}]]]
Out[6]= 2
```

## Implementation notes

**Algorithm.** `builtin_morphologicalcomponents` labels the connected components
of the foreground by **two-pass union-find**, returning an integer matrix with
background `0` and components numbered `1..k` in raster order of first
appearance. `MorphologicalComponents[image, t]` takes pixels strictly above `t`
as foreground (default `0`); `CornerNeighbors -> False` switches from the default
8-connectivity to 4. The **first pass** walks in raster order and can only see
already-visited neighbours (W, NW, N, NE for 8-connectivity), assigning
provisional labels and recording equivalences with `cc_union`/`cc_find` (path
compression; union by lower index, kept deterministic rather than by rank since
the second pass relabels anyway). A U-shaped region is why one pass is
insufficient — its two arms get different labels that the base reveals equal. The
**second pass** resolves each pixel to its representative and **relabels to
`1..k`** in raster order of first appearance, so there are no gaps and
`Max[result]` is exactly the component count. Connectivity is the one
discriminating property: two pixels touching only at a corner are **one**
component under 8 and **two** under 4.

**Data structures.** A `size_t` parent array for union-find over the pixel grid
and a label matrix; the result is a plain nested-`List` **integer matrix**, *not*
an `Image` — deliberately, because `Image` type inference would call a label
array of `1..12` a `"Byte"` image and `ImageData` would then divide every label
by 255. Labels are indices, not brightnesses, and must not be scaled.

**Complexity / limits.** `O(width · height · α(width · height))` — effectively
linear in the pixel count with path-compressed union-find. Foreground is
`pixel > t`; a colour image's pixels are compared on their stored channel values.

**Attributes:** `Protected`.

## References

**See also:** [Image](../../image-processing/Image/), [ImageData](../../image-processing/ImageData/), [Max](../../data-structures/Max/), [List](../../other-advanced/List/)

- R. M. Haralick and L. G. Shapiro, *Computer and Robot Vision*, vol. 1 (Addison-Wesley, 1992), ch. 2 (connected components labelling).
- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`MorphologicalComponents[image]` labels the connected components of the
foreground, giving an **integer matrix** with background `0` and components
numbered `1..k` in raster order of first appearance.
`MorphologicalComponents[image, t]` takes pixels above `t` as foreground (default
`0`), and `CornerNeighbors -> False` uses 4-connectivity instead of the default
8.

The two-pass union-find assigns provisional labels in raster order, records
equivalences (a U-shape is the case that needs the second pass), then relabels to
`1..k` with no gaps — so `Max` of the result is the component count.
Connectivity is the only discriminating property: two pixels touching at a corner
are one component under 8 and two under 4. The result is a matrix, not an `Image`,
on purpose: `Image` would type-infer a `1..12` label array as `"Byte"` and
`ImageData` would divide every label by 255.
