---
source: src/imagefilter.c
---
**Algorithm.** `builtin_meanfilter` averages each pixel over a `(2r+1) × (2r+1)`
neighbourhood. The radius must be a non-negative integer. The key design choice
is that this IS a convolution with a normalised box: the routine builds a `n × n`
kernel of `1/n²` (`n = 2r+1`) and runs it through the shared
`convolve_dispatch`, rather than through a separate averaging loop — two
implementations of one identity is how the identity quietly stops holding. A full
rectangle factorises, so `convolve_dispatch` takes the separable path and the
cost is `kw + kh` taps per pixel per axis rather than `kw · kh`. A volume takes
the rank-3 path (`mean3_run`). Border reads clamp to the nearest edge pixel, the
same `"Fixed"` padding every convolution here uses.

**Data structures.** A `n × n` flat `double` box kernel; `src`/`dst` flat
height · width · channels unit-scale buffers. The result is a `"Real"` image
(`image_build_real`). No bespoke accumulator — it is the convolution core.

**Complexity / limits.** `O(width · height · channels · (kw + kh))` via the
separable path, independent of the kernel area. Radius is capped at 256. The box
is normalised (unlike `BoxMatrix`, whose entries are `1`), so the output is a
true mean rather than a sum.
