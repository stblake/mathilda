---
references:
  - "N. Otsu, *A threshold selection method from gray-level histograms*, IEEE Trans. Systems, Man, and Cybernetics **9** (1979) 62-66."
source: src/imagefilter.c
---
**Algorithm.** `builtin_binarize` thresholds an image into a `"Bit"` image.
`Binarize[image]` chooses the threshold by Otsu's method (`img_otsu`, see
`FindThreshold`); `Binarize[image, t]` uses the explicit level `t`. A colour
image is reduced to a Rec. 601 luminance plane first (`img_grey_plane`). The
decision is **strictly above**: a pixel with value `> t` becomes `1`, so a pixel
exactly *at* `t` becomes `0`. That boundary matters because "above" and "at or
above" differ on exactly the pixels a threshold was chosen to sit between, and a
test pins it. The mask is built by `bit_image_from_mask` (shared with
`LocalAdaptiveBinarize`, so the two cannot disagree on the output type) into an
`Image[..., "Bit"]`, typed `"Bit"` rather than `"Real"` because the result is
binary *by construction* — a later `ImageData` scales nothing and the image
reports its binary nature. A volume takes the rank-3 path (`binarize3_run`). With
the default (Otsu) threshold, a degenerate single-value image declines, since no
threshold splits one cluster.

**Data structures.** A flat luminance `double` plane and an `unsigned char` 0/1
mask of `width · height`, packed into a `"Bit"` image buffer
(`image_build_bit`). No multi-channel buffer survives — binarisation collapses to
one plane.

**Complexity / limits.** `O(width · height)` for the threshold comparison, plus
`O(width · height + bins)` for an Otsu default. With an explicit `t` it never
declines; with the default it declines on a single-value image.
