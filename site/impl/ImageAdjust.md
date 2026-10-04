---
source: src/imagefilter.c
---
**Algorithm.** `builtin_imageadjust` has two modes. `ImageAdjust[image]` is a
**full-range stretch**: it finds the min `lo` and max `hi` over the whole buffer
and maps `v -> (v − lo)/(hi − lo)`, so the darkest pixel becomes exactly `0` and
the brightest exactly `1`. It is idempotent (a second stretch is the identity),
and a constant image has `span = 0` and is returned unchanged — dividing by zero
is not the answer and mapping the single value to either end would be arbitrary.
`ImageAdjust[image, {c, b}]` and `ImageAdjust[image, {c, b, g}]` apply contrast
`c`, brightness `b` and gamma `g` by a stated curve: `v' = (v − 1/2)(1 + c) + 1/2
+ b`, clipped to `[0, 1]`, then raised to the power `1/g`. Contrast pivots about
mid-grey so a contrast change does not also shift brightness; clipping precedes
gamma because a negative base has no real power; a non-positive gamma is rejected
(not a curve). This curve is Mathilda's documented choice, not a claim of
bit-compatibility with Mathematica. Both ranks are handled — a volume is loaded
with `image3d_load` and rebuilt with `image3d_build_real`.

**Data structures.** One flat unit-scale `double` buffer (`image_load` /
`image3d_load`) mutated in place and rebuilt as a `"Real"` image. No histogram and
no second buffer — the stretch is a single min/max pass then a linear map; the
parametric form is a single per-element pass.

**Complexity / limits.** `O(n)` where `n = width · height · channels` (× depth) —
one pass for min/max and one to apply the map. The parametric curve is applied
across all channels uniformly.
