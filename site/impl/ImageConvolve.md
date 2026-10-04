---
source: src/imagefilter.c
---
**Algorithm.** `builtin_imageconvolve` convolves an image with a rank-2 numeric
kernel. It dispatches on the **image** rank — a volume takes the rank-3 path
(`ker3_load` + `convolve3_run`), a plane the rank-2 path — so a rank-3 kernel
handed to a plane is declined, not reinterpreted. The plane path loads the image
to a flat unit-scale buffer (`image_load`), loads the kernel (`ker_load`,
rejecting ragged matrices), and calls `convolve_dispatch`. True convolution
**reflects** the kernel on both axes before summing, which is the only difference
from `ImageCorrelate`; the two agree exactly on a symmetric kernel. Padding is
`"Fixed"`: out-of-range reads clamp to the nearest edge pixel (`clampi` on
`int64_t`, because the index goes negative near the top edge), so a constant
image convolved with a kernel summing to 1 comes back unchanged everywhere
including the border — zero padding would darken the edges. The core
`convolve_planes` splits each image into an **interior** (every tap in range, a
plain forward dot product over the once-reversed kernel, which vectorises) and a
thin clamped **border**; `convolve_dispatch` first tries a separable
factorisation (`ker_separable`, verified to ~1e-16 since treating a non-separable
kernel as separable computes a *different* filter), falling back to the direct
form or, when built with FFTW, a transform. Each colour channel is convolved
independently.

**Data structures.** Flat `double` buffers: `src`/`dst` are height · width ·
channels unit-scaled arrays, the kernel a `kh × kw` `double` matrix. The result
is always a `"Real"` image (`image_build_real`) — a filtered byte is not
generally a byte. All arithmetic is float64; vImage's float32 convolution is
deliberately not used, as it would drop six digits a CAS test asserts.

**Complexity / limits.** Direct form `O(width · height · channels · kw · kh)`;
the separable path drops it to `kw + kh` taps per pixel per axis when the kernel
factorises; the FFTW transform path is used only when its cost model wins on a
large non-separable kernel. The kernel must be a rectangular numeric matrix with
real entries.
