# ImageConvolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageConvolve[image, kernel] convolves image with the rank-2 numeric kernel. This is true convolution: the kernel is REFLECTED before summing, so it differs from correlation on an asymmetric kernel (the two agree exactly on a symmetric one such as a Gaussian or a box). Out-of-range reads clamp to the nearest edge pixel, replicating the border, so a constant image convolved with a kernel summing to 1 comes back unchanged everywhere including the edges -- zero padding would darken them. The result is always a "Real" image of the same dimensions, since a filtered byte is not generally a byte. Each colour channel is convolved independently.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= ImageData[ImageConvolve[Image[{{0., 1., 0.}}], {{1, 2, 3}}]]
Out[1]= {{1.0, 2.0, 3.0}}
```

### Applications (4)

```mathematica
In[2]:= ImageData[ImageConvolve[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}]]
Out[2]= {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}}

In[3]:= ImageData[ImageConvolve[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], {{1, 2, 3}}]]
Out[3]= {{0.0, 0.0, 0.0}, {1.0, 2.0, 3.0}, {0.0, 0.0, 0.0}}

In[4]:= ImageData[ImageCorrelate[Image[{{0., 0, 0}, {0, 1., 0}, {0, 0, 0}}], {{1, 2, 3}}]]
Out[4]= {{0.0, 0.0, 0.0}, {3.0, 2.0, 1.0}, {0.0, 0.0, 0.0}}

In[5]:= ImageData[ImageConvolve[Image[{{0.5, 0.5}, {0.5, 0.5}}], GaussianMatrix[1]]]
Out[5]= {{0.5, 0.5}, {0.5, 0.5}}
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [ImageCorrelate](../../image-processing/ImageCorrelate/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`ImageConvolve[image, kernel]` convolves with a rank-2 numeric kernel. It is
**true convolution**: the kernel is reflected on both axes before summing, so it
differs from `ImageCorrelate` on an asymmetric kernel (the two agree exactly on a
symmetric one such as a Gaussian or a box).

Out-of-range reads clamp to the nearest edge pixel (`"Fixed"` padding), so a
constant image convolved with a kernel summing to 1 comes back unchanged
everywhere, including the edges — zero padding would darken them. The result is
always a `"Real"` image of the same dimensions, since a filtered byte is not
generally a byte, and each colour channel is convolved independently. A separable
kernel is factorised automatically into two 1-D passes; a volume takes a rank-3
convolution, and a rank-3 kernel handed to a plane is declined rather than
guessed at.
