---
source: src/imagecolor.c
---
**Algorithm.** `builtin_colorreplace` recolours pixels within a tolerance of a target colour.
The rule(s) — `old -> new`, or a list of them — are parsed by `read_rules`, which accepts
`RGBColor[r, g, b]`, `GrayLevel[v]`, a bare number, or `{r, g, b}` and expands grey specs to
three equal components so grey and colour pixels compare in one path. The image is decoded to a
unit buffer; the default tolerance is `0.02` (a stated constant, because at 0 only a
bit-identical colour matches, which after any filtering is nothing). For each pixel the
Euclidean RGB distance to every rule's source is measured, and the **nearest** rule within
tolerance wins — not the first, so overlapping rules do not make the answer depend on the order
they were written. The matched pixel takes the rule's target; an unmatched pixel is copied
through, and an alpha channel always passes through untouched (transparency is not a colour).

If any target is non-grey while the input is grey, the output is promoted to three channels —
flattening the new colour to its luminance would give grey where the caller asked for red.

**Data structures.** Fixed-capacity rule tables `from[64][3]`/`to[64][3]`; one decoded unit
buffer in, a fresh `width · height · out_channels` buffer out, wrapped by `image_build_real` as
a packed `"Real"` image (every image head returns a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(width · height · rules)`. At most 64 rules; the tolerance must be
non-negative. Distance is Euclidean in RGB only — there is no perceptual colour-space metric.
