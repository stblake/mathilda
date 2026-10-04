---
references:
  - "J. Serra, *Image Analysis and Mathematical Morphology* (Academic Press, 1982)."
source: src/imagefilter.c
---
**Algorithm.** `builtin_closing` is `morph_builtin(res, MORPH_DILATE, two_pass =
true)` — a **dilation followed by an erosion** with the *same* structuring
element. The shared `morph_builtin` runs `morph_run` once for the first pass and,
when `two_pass`, again with the opposite operation (`MORPH_DILATE` then
`MORPH_ERODE`), reusing the same support mask both times. Using the same element
both times is what makes the pair **idempotent** — `Closing[Closing[f]]` equals
`Closing[f]` — and what makes it a closing rather than merely two smoothings; a
different second element would still smooth but would not be a closing.
Geometrically it fills dark features smaller than the element while leaving larger
ones close to their original size. Each pass inherits `Dilation`'s machinery: the
separable van Herk–Gil-Werman max/min for a full rectangle, `morph_direct` for an
arbitrary element, and replicate padding. Because dilation and erosion are duals
and the padding is self-dual, the family brackets the image:
`Erosion <= Opening <= image <= Closing <= Dilation` pointwise everywhere.

**Data structures.** Two flat `double` scratch buffers `a` (first pass) and `b`
(second pass), plus the `unsigned char` support mask. A `"Bit"` image stays
`"Bit"` (`image_build_typed`), every other type gives `"Real"`.

**Complexity / limits.** Twice a single morphology pass: `O(width · height ·
channels)` amortised for a full rectangle (independent of `r`), `O(width ·
height · channels · |support|)` for an arbitrary element. A volume takes the
rank-3 two-pass path.
