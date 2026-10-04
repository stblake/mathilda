---
source: src/imagecompose.c
---
**Algorithm.** `builtin_removealphachannel` drops an image's alpha channel.
`RemoveAlphaChannel[image]` simply discards it, keeping each colour channel's
stored value. `RemoveAlphaChannel[image, b]` instead **composites over** a
background of brightness `b`: for each pixel it reads the alpha `al` (1.0 for an
image with no alpha) and writes `v · al + b · (1 − al)` per colour channel —
which is the difference between *forgetting* the transparency and *resolving*
it. A half-transparent white pixel over black (`b = 0`) becomes grey, where
plain dropping would leave it white. The output channel count is
`colour_count(im.c)` — the colour channels without the alpha — so a 4-channel
RGBA image becomes 3-channel RGB and a 2-channel grey-plus-alpha becomes
1-channel grey.

**Data structures.** The image is read into an `Img` struct via `img_read`;
`sample(&im, x, y, k, cc)` reads channel `k`, and `has_alpha`/`colour_count`
derive the alpha presence and colour-channel count from the channel count. The
output is a flat `double` buffer of `width · height · colour_count`, built into a
`"Real"` image (`image_build_real`).

**Complexity / limits.** `O(width · height · channels)` — a single per-pixel
pass. The composite form needs a scalar background brightness; an image with no
alpha channel is returned with its colour channels intact (the composite reads
`al = 1`, leaving values unchanged).
