# BoxMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BoxMatrix[r] gives a (2r+1) x (2r+1) matrix of 1s. It is NOT normalised, matching Mathematica, so ImageConvolve[image, BoxMatrix[1]] is nine times too bright; the normalised version is a mean filter. Kept faithful rather than helpfully rescaled, since a caller using BoxMatrix in arithmetic needs the ones.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (1)

```mathematica
In[1]:= BoxMatrix[1]
Out[1]= {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}}
```

### Applications (3)

A 3x3 block of ones

```mathematica
In[2]:= BoxMatrix[1]
Out[2]= {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}}
```

Radius 0 is the 1x1 matrix {{1}}

```mathematica
In[3]:= BoxMatrix[0]
Out[3]= {{1}}
```

Used as a morphology element

```mathematica
In[4]:= ImageData[Dilation[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}], BoxMatrix[1]]]
Out[4]= {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}}
```

## Implementation notes

**Algorithm.** `builtin_boxmatrix` builds a `(2r+1) × (2r+1)` matrix of `1`s for a
non-negative integer radius `r` (capped at 512). It is **not** normalised — that is
Mathematica's definition, and a trap worth naming: convolving with `BoxMatrix[1]` multiplies
brightness by the element count, so `ImageConvolve[img, BoxMatrix[1]]` is nine times too
bright. The normalised version is a mean filter (`MeanFilter`). It is kept faithful rather than
helpfully rescaled, because a caller reaching for `BoxMatrix` in an arithmetic expression, or
as the structuring element of a morphology op, needs the ones.

**Data structures.** Returns a nested `List` of `List`s of integer `1`s — an ordinary matrix,
**not** an image, so it is not a packed image buffer and the `make check-image-packing` audit
does not apply to it. It is consumed as a convolution kernel (`ker_load`) or a morphology
element (`morph_support`, which reads only its nonzero *support*, so the values beyond "present"
do not enter a flat max/min).

**Complexity / limits.** `O((2r+1)²)` construction. The radius must be a non-negative integer;
a fractional radius has no matrix size and declines.

**Attributes:** `Protected`.

## References

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`BoxMatrix[r]` is a `(2r+1) × (2r+1)` matrix of `1`s. It is **not** normalised, matching
Mathematica — so `ImageConvolve[image, BoxMatrix[1]]` is nine times too bright; the normalised
version is a mean filter, which is what `MeanFilter` provides.

The result is an ordinary integer matrix, not an image. Its common uses are as a convolution
kernel for `ImageConvolve` and as the structuring element for the morphology operators
(`Dilation`, `Erosion`, `Opening`, `Closing`), where only its nonzero support matters.
