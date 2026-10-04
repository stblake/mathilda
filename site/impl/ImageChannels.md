---
source: src/image.c
---
**Algorithm.** `builtin_imagechannels` returns the number of colour channels as
an integer. It tries `image3d_info` first and then `image_info`, each of which
reports the channel count as a by-product of validating the shape. The count is
read off the pixel array's rank: a rank-2 (height × width) array is grey and
reports `1`; a rank-3 (height × width × channels) array reports its last
dimension — `3` for RGB, `4` when an alpha channel is present, `2` for
grey-plus-alpha. For a volume the same rule applies to the trailing axis of the
depth × height × width (× channels) buffer. The channel axis is the **innermost,
interleaved** dimension, which is the order every filter's flat `image_load`
buffer also uses.

**Data structures.** Reads only the shape metadata of the canonical image node;
no buffer is loaded. `ImageChannels` is on `pack.c`'s `AWARE` list, so a packed
image is inspected in place.

**Complexity / limits.** `O(1)` — the shape walk that `image_info` already does.
Returns unevaluated for a non-image.
