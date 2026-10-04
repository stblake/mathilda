---
source: src/imagegeom.c
---
**Algorithm.** `builtin_imageresize` resamples an image to a new size.
`ImageResize[image, {w, h}]` gives `w × h` pixels; `ImageResize[image, w]` gives
width `w` with the height following to preserve the aspect ratio
(`round(sh · w / sw)`). `Resampling -> "Nearest" | "Bilinear" | "Average"`
selects the method; the default `Automatic` uses **area averaging when either
axis shrinks** and bilinear otherwise. That default is about aliasing: Nyquist
requires every frequency above half the new sampling rate to be removed *before*
resampling, and point-sampling a shrinking image (nearest) destroys them with no
recovery — a fine checkerboard reduced by nearest returns a flat field. Area
averaging (`rs_average`) is a box prefilter and a resample in one pass, using
**true fractional coverage** so a 3→2 reduction is as correct as 4→2, and it is
exact for integer reduction factors. Enlarging has no frequencies to remove, so
bilinear (`rs_bilinear`) is used there. All three maps are **pixel-centred**:
a destination pixel `i` covers `[i·s, (i+1)·s)` with centre `(i + 0.5)·s`, so the
bilinear map is `sx = (i + 0.5)·s − 0.5` — avoiding the half-pixel shift and
edge asymmetry that the naive `sx = i·s` introduces at any scale other than 1:1.
A volume takes the rank-3 path (`resize3_run`). The result is a `"Real"` image.

**Data structures.** Flat `double` buffers: `src` of `sw · sh · channels`, `dst`
of `dw · dh · channels`, resampled by the chosen kernel (each channel handled
independently, edge indices clamped by `clamp_idx`). Built into a `"Real"` image
(`image_build_real`).

**Complexity / limits.** Nearest `O(dw · dh · c)`; bilinear `O(dw · dh · c)` with
four taps per pixel; area averaging `O(dst area × average source coverage)`.
Target sizes must be positive integers (a fractional pixel count has no meaning
and is declined rather than silently rounded), capped at `10⁶` per axis.
