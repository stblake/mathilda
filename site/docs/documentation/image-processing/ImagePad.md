# ImagePad

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImagePad[image, m] pads m pixels on every side; ImagePad[image, {{left, right}, {bottom, top}}] pads each side separately, in Mathematica's VISUAL order -- so `top` adds rows at the start of the data, since row 1 is the top of the image. Negative amounts crop, but may not erase the image. ImagePad[image, m, v] fills with the value v (default 0); ImagePad[image, m, "Fixed"] replicates the edge pixel, the same boundary rule the filters use, so padding then filtering composes with it; ImagePad[image, m, "Reflected"] mirrors WITHOUT repeating the edge -- {1,2,3} padded by 1 gives {2,1,2,3,2}, not {1,1,2,3,3}, because doubling the edge sample biases any later average toward the border. Reflection uses a period of 2n-2, so padding deeper than the image still works.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

Pad one pixel on every side: 2x1 -> 4x3

```mathematica
In[1]:= ImageDimensions[ImagePad[Image[{{0., 1.}}], 1]]
Out[1]= {4, 3}
```

A 1x1 image padded with the default value 0

```mathematica
In[2]:= ImageData[ImagePad[Image[{{0.5}}], 1]]
Out[2]= {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.0}}
```

The Fixed mode replicates the edge pixel instead

```mathematica
In[3]:= ImageData[ImagePad[Image[{{0.5}}], 1, "Fixed"]]
Out[3]= {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}}
```

Reflection mirrors without repeating the edge

```mathematica
In[4]:= ImageData[ImagePad[Image[{{1., 2., 3.}}], {{1, 1}, {0, 0}}, "Reflected"]]
Out[4]= {{2.0, 1.0, 2.0, 3.0, 2.0}}
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [ImageCrop](../../image-processing/ImageCrop/), [Image3D](../../image-processing/Image3D/), [ImageConvolve](../../image-processing/ImageConvolve/)

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`ImagePad[image, m]` pads `m` pixels on every side; `ImagePad[image, {{left, right}, {bottom,
top}}]` names each side in Mathematica's visual order — so `top` adds rows at the *start* of
the data, since row 1 is the top of the image. Negative amounts crop but may not erase the
image.

The third argument sets the fill: a value (default 0); `"Fixed"`, which replicates the edge
pixel — the same boundary rule the filters use, so padding then filtering composes with it; or
`"Reflected"`, which mirrors *without* repeating the edge (`{1,2,3}` padded by 1 →
`{2,1,2,3,2}`, not `{1,1,2,3,3}`), because doubling the edge sample biases any later average
toward the border. `ImagePad` and `ImageCrop` are exact inverses.
