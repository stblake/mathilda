---
source: src/imagefilter.c
---
**Algorithm.** `builtin_opening` is `morph_builtin(res, MORPH_ERODE, true)`: an **erosion
followed by a dilation** with the *same* structuring element. It removes bright features
smaller than the element while leaving larger ones close to their original size. Both passes go
through `morph_run` (the separable van Herk–Gil–Werman min/max, bit-exact with the direct
form); `two_pass` runs the second with the opposite operator (`MORPH_DILATE`).

Using the same element both times is what makes the pair **idempotent**:
`Opening[Opening[f]] = Opening[f]`, the defining property of an opening and the reason opening
twice is not a sharpening loop. A different second element would still smooth but would no
longer be an opening. Opening sits in the morphology ordering
`Erosion ≤ Opening ≤ f ≤ Closing ≤ Dilation` pointwise everywhere — the replicate padding is
what keeps that true at the border.

**Data structures.** A decoded unit buffer plus one scratch buffer for the intermediate
erosion; the result through `image_build_typed`/`_real`. A `"Bit"` image stays `"Bit"` (a
max/min of stored values is a stored value), other types become `"Real"`; either way it is a
packed buffer (`make check-image-packing`).

**Complexity / limits.** `O(pixels)` for a full rectangle (two van Herk passes),
`O(pixels · |support|)` for an arbitrary element. Radius ≤ 256 planar, ≤ 64 volumetric.
