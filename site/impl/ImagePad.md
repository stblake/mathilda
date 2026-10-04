---
source: src/imagegeom.c
---
**Algorithm.** `builtin_imagepad` adds (or, for negative amounts, removes) a border.
`ImagePad[image, m]` pads `m` on every side; `ImagePad[image, {{left, right}, {bottom, top}}]`
names each side in Mathematica's **visual** order — so `top` adds rows at the *start* of the
array, since row 0 is the top of the image. The optional third argument chooses the fill:
a constant value (`PAD_VALUE`, default 0); `"Fixed"` replicates the edge pixel (the same
boundary rule the filters use, so padding then filtering composes with it); or `"Reflected"`,
which mirrors **without** repeating the edge (`{1,2,3}` padded by 1 → `{2,1,2,3,2}`), using a
period of `2n − 2` so arbitrarily deep padding works — doubling the edge sample would bias any
later average toward the border. `pad_src_index` maps a padded coordinate back into the source
or reports it outside.

The interior is a **row-block `memcpy`**, not a per-pixel map: every mode agrees on the
interior (it is the source, unshifted), so only the thin frame is computed per pixel. The first
version mapped every pixel and cost 0.57 ms on a 512 × 512 constant pad against NumPy's 0.055;
the span-per-row form also handles negative padding with no second code path. A volume takes
`pad3_run`, the same structure one rank up.

**Data structures.** One decoded unit buffer in; a fresh `nw · nh · c` buffer out, wrapped by
`image_build_real` as a packed `"Real"` image (every image head returns a packed buffer —
`make check-image-packing`).

**Complexity / limits.** `O(output pixels)`, dominated by the interior `memcpy`. Negative
padding may crop but not erase the image (`nw < 1` or `nh < 1` declines); amounts must be
integers.
