---
references:
  - "P. F. Felzenszwalb and D. P. Huttenlocher, *Distance Transforms of Sampled Functions*, Theory of Computing **8** (2012) 415-428."
source: src/imagefilter.c
---
**Algorithm.** `builtin_distancetransform` replaces each pixel by its **exact** Euclidean
distance to the nearest background pixel. `DistanceTransform[image, t]` takes pixels above `t`
as foreground (default 0). The image is reduced to grey and seeded: `0` at background,
`DT_INF` at foreground, so the transform reports the distance to the nearest *seed*. It then
runs one lower-envelope pass per axis (`dt_1d`): for a line, `D(x) = min_y ((x−y)² + f(y))` is
the lower envelope of identical parabolas translated to each `y` and raised by `f(y)`. Because
all parabolas share curvature, any two intersect once, so the envelope is built in a single
sweep maintaining a stack of visible parabolas — `O(n)` per row, no sorting.

Exactness comes from choosing the envelope method over the classic two-pass **chamfer**
transform, which cannot represent `√2` with integer steps and gets diagonal distances a few
percent wrong — invisible on a picture, fatal to a test. Separability is **exact** here (unlike
the median's) because squared Euclidean distance is a *sum* over axes, so minimising it
decomposes per axis; the square root is taken once at the end, not per pass (per-axis roots
would be wrong, not merely slow). A 3-4-5 triangle reads exactly 5. A volume adds a third
`dt_1d` pass.

**Data structures.** A grey seed buffer `f` plus scratch arrays `v` (parabola indices) and `z`
(intersections) per line; the result through `image_build_real` as a packed single-channel
`"Real"` image (every image head returns a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(pixels)` — linear in the number of pixels, independent of how far
the nearest seed is. An all-foreground image has no seed, so every pixel reads `DT_INF`.
