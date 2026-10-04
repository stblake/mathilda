---
references:
  - "J. Serra, *Image Analysis and Mathematical Morphology* (Academic Press, 1982)."
  - "M. van Herk, *A fast algorithm for local minimum and maximum filters on rectangular and octagonal kernels*, Pattern Recognition Letters **13** (1992) 517-521."
  - "J. Gil and M. Werman, *Computing 2-D min, median, and max filters*, IEEE TPAMI **15** (1993) 504-507."
source: src/imagefilter.c
---
**Algorithm.** `builtin_dilation` is `morph_builtin(res, MORPH_DILATE, two_pass =
false)` — the pointwise **maximum** over a structuring element.
`Dilation[image, r]` uses a `(2r+1) × (2r+1)` square; `Dilation[image, elem]`
uses the **support** of the matrix `elem` — its nonzero positions
(`morph_support`), so the element's values do not enter the maximum. This is
*flat* morphology, which is what keeps `Dilation[img, BoxMatrix[1]]` and
`Dilation[img, 1]` the same operation. A full rectangle is separable for the
maximum exactly as for a sum (max over a rectangle = max over rows of the maxima
over columns), so `morph_separable` does a row pass then a column pass, and each
1-D pass uses the **van Herk–Gil-Werman** trick: block the line into runs of `k`,
keep a prefix and a suffix running maximum per block, and any width-`k` window is
`max(suffix at its start, prefix at its end)` — three comparisons per pixel
*regardless of the radius*, which is what makes morphology usable at large `r`.
The fast path is exact (max is associative and idempotent) and a test asserts it
agrees bit-for-bit with the naive `morph_direct`. An arbitrary (non-rectangular)
element falls back to `morph_direct`. Padding replicates the border, the same
rule the convolutions use, so `Dilation >= image` holds at the edges too.

**Data structures.** Flat `double` buffers; the element is an `unsigned char`
support mask plus a `full` flag. The separable path uses padded scratch lines and
`pre`/`suf` arrays sized to the longer axis plus the kernel. A `"Bit"` image
gives a `"Bit"` image (`image_build_typed`, since a max of stored 0/1 values is
itself 0/1); every other type gives `"Real"`.

**Complexity / limits.** Separable full-rectangle: `O(width · height · channels)`
amortised, **independent of `r`**. Arbitrary element: `O(width · height ·
channels · |support|)`. A volume takes the rank-3 path, with an integer radius
only (an arbitrary 3-D element is not separable).
