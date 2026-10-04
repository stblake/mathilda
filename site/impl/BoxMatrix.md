---
source: src/imagefilter.c
---
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
