---
source: src/imagefilter.c
---
**Algorithm.** `builtin_medianfilter` replaces each pixel with the median over a `(2r+1)`
square (cube for a volume). Per pixel it gathers the `(2r+1)²` neighbourhood with replicate
("Fixed") boundary clamping into a scratch window and selects the lower-middle order statistic
with `window_median`, a **quickselect** (median-of-three pivot, so a nearly-sorted window — the
common case in a smooth image — does not hit the quadratic case). For an even window the lower
middle is taken, never the average of two, because a rank filter's output must be one of its
inputs.

The median is the **one operator here that is not separable**, and that is the point: a sum, a
max and a min all decompose because they ignore grouping, but the median depends on a value's
*rank* within the whole window, and grouping destroys rank — the median of the row-medians of
`{{1,2,9},{3,4,5},{6,7,8}}` is 4 where the true median is 5. So the window really is gathered.
Its value is that it removes an isolated outlier *exactly* (a single bright pixel in a constant
field vanishes), where a Gaussian only attenuates and smears it — the filter for
salt-and-pepper noise.

**Data structures.** A decoded unit buffer plus a `(2r+1)²` scratch window; the result through
`image_build_real` as a packed `"Real"` image (every image head returns a packed buffer —
`make check-image-packing`).

**Complexity / limits.** `O(pixels · k² )` to gather plus `O(k²)` average per quickselect;
radius ≤ 64 planar, ≤ 16 volumetric (the cube `(2r+1)³` grows fast).
