---
source: src/imagefilter.c
---
**Algorithm.** `builtin_colorconvert` converts an image (plane or `Image3D`) to greyscale.
Only `"Grayscale"` (or `"Gray"`) is accepted: the other colour spaces Mathematica supports
(LAB, HSB, XYZ, …) each carry their own white point and transfer-function decisions, and
accepting the name while doing something approximate would be worse than declining it. The
reduction uses the **Rec. 601** luminance weights `0.299 R + 0.587 G + 0.114 B` (`img_to_grey`
for a plane, `img3_grey_volume` for a volume — one place, so both ranks agree), because the eye
is not equally sensitive across the spectrum: green carries most perceived brightness and blue
almost none, so a plain average would put pure red and pure blue on the same side of a
threshold when perceptually they are far apart.

What is exact is documented honestly. An **already-grey** image is copied through bit for bit.
An image whose three channels are merely **equal** is exact only to within an ulp, and *whether*
depends on the value: those weights sum to `0.9999999999999999` in the order they are applied
but exactly `1.0` in any order beginning with `0.114`, so the final rounding lands on the input
for some values and one ulp below for others. The weights are the standard's and are not
adjusted to compensate — a triple hand-tuned to sum to exactly `1.0` in double would no longer
be Rec. 601.

**Data structures.** One decoded unit buffer in; a fresh single-channel buffer out, wrapped by
`image_build_real` / `image3d_build_real` as a packed `"Real"` image (every image head returns
a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(pixels)`. Only the greyscale target is implemented.
