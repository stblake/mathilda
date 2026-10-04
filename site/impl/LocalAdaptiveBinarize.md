---
source: src/imagefilter.c
---
**Algorithm.** `builtin_localadaptivebinarize` thresholds each pixel against a statistic of its
own `(2r+1) × (2r+1)` neighbourhood rather than a global level — the only way to binarize
unevenly lit content, since if one half of a page is darker than the other, no single number
separates ink from paper in both halves. The rule is `g > c1·mean + c2·stddev + c3`, with the
coefficients defaulting to `{1, 0, 0}` (Bradley's mean thresholding);
`LocalAdaptiveBinarize[image, r, {c1, c2, c3}]` sets them, and a negative `c2` gives Sauvola's
method, tightening the threshold where the neighbourhood is busy.

The window statistics are `O(1)` per pixel via **summed-area tables**: a padded prefix-sum
table `s1` for the mean, and a second table `s2` of squares (built only when `c2 ≠ 0`, since
the default never needs the standard deviation) for the variance `s2/area − mean²`, clamped at
0 before `sqrt` because cancellation can push a uniform window a hair negative. Without the
tables a radius-16 window would be 1089 taps per pixel.

**Data structures.** A grey plane, one or two `(PH+1)(PW+1)` summed-area tables, and a
`width · height` byte mask; the result is built by `bit_image_from_mask` as a `"Bit"` image —
the output is binary by construction, so it is typed `"Bit"` and handed back as a packed buffer
(every image head returns one — `make check-image-packing`).

**Complexity / limits.** `O(pixels)` regardless of radius. Radius 1–256; colour is reduced to
luminance first. The coefficient triple must be numeric.
