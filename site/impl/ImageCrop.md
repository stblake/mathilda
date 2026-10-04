---
source: src/imagegeom.c
---
**Algorithm.** `builtin_imagecrop` has two modes. `ImageCrop[image, {w, h}]` is a **centred**
crop: it takes the `w × h` window `x0 = (W − w)/2`, `y0 = (H − h)/2`, any odd remainder going
to the right and bottom. That is floor division, the same convention the kernel centres and
`ImagePad` use, which is what makes `ImageCrop[ImagePad[image, m], ImageDimensions[image]]`
exactly the original image. A crop may not enlarge (`w > W` or `h > H` declines).

`ImageCrop[image]` with no size instead **trims a uniform border**: `border_trim` shrinks each
edge inward one row/column at a time while that outer ring equals the reference colour, tested
by `span_uniform` per edge so a border uniform on three sides and not the fourth trims the
three. The reference colour is read from the top-left pixel rather than assumed black — a
scanned page's margin is white, and assuming black would trim nothing. An entirely uniform
image has no content to keep, so it comes back unchanged rather than as a zero-sized image. A
volume takes only the sized form (`crop3_run`); trimming a border in 3-D is ambiguous (which
faces? a shell or a box?) and is declined.

**Data structures.** One decoded unit buffer in; a fresh `cw · ch · c` buffer out, wrapped by
`image_build_real` as a packed `"Real"` image (every image head returns a packed buffer —
`make check-image-packing`). `border_trim` keeps four running edge indices.

**Complexity / limits.** The sized crop is `O(cw · ch · c)`; the border trim adds
`O(perimeter · trimmed-depth)`. Crop sizes must be positive integers no larger than the image.
