---
source: src/imagecompose.c
---
**Algorithm.** `builtin_setalphachannel` attaches or replaces an image's opacity. The result
always carries one more channel than the source's colour-channel count (`colour_count` strips
any existing alpha first, so the output is colour+alpha, never doubled). The new opacity is
either:

- a constant — `SetAlphaChannel[image, a]` with `a` a number in `[0, 1]`, written into the
  alpha slot of every pixel; or
- a mask — `SetAlphaChannel[image, maskimage]`, which must match the source in size (a
  mismatched mask is declined, not silently resampled) and is read as **grey**: its own colour
  channels are averaged, so a colour mask is not quietly taken as its red channel alone; or
- fully opaque `1.0` when no second argument is given.

Colour samples are copied through `sample`, which replicates a grey source across colour
channels so the operation has one code path for grey and colour inputs.

**Data structures.** One or two decoded `Img` buffers in; a fresh
`width · height · (colour_count + 1)` buffer out, wrapped by `image_build_real` as a packed
`"Real"` image (every image head returns a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(width · height · channels)`, one pass. A mask of the wrong
dimensions is refused; the opacity constant must lie in `[0, 1]`.
