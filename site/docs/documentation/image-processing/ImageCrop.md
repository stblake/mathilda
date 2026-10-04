# ImageCrop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageCrop[image, {w, h}] crops to w x h about the centre, any odd remainder going to the right and bottom -- the same floor-division convention the kernel centres use, which is what makes ImageCrop[ImagePad[image, m], ImageDimensions[image]] exactly the original image. A crop may not enlarge. ImageCrop[image] instead TRIMS A UNIFORM BORDER, asking how much of the frame carries no information; the border colour is read from a corner rather than assumed black, since a scanned page's margin is white. An entirely uniform image comes back unchanged, there being no content to keep and a zero-sized image not being one.`**

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

Crop to 1x1 about the centre

```mathematica
In[1]:= ImageDimensions[ImageCrop[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], {1, 1}]]
Out[1]= {1, 1}
```

The centre pixel survives

```mathematica
In[2]:= ImageData[ImageCrop[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], {1, 1}]]
Out[2]= {{1.0}}
```

Trim the uniform black border

```mathematica
In[3]:= ImageDimensions[ImageCrop[Image[{{0., 0., 0., 0.}, {0., 1., 1., 0.}, {0., 0., 0., 0.}}]]]
Out[3]= {2, 1}
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [ImagePad](../../image-processing/ImagePad/), [Image3D](../../image-processing/Image3D/), [ImageConvolve](../../image-processing/ImageConvolve/)

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`ImageCrop[image, {w, h}]` crops to `w × h` about the centre, any odd remainder going to the
right and bottom — the same floor-division convention the kernel centres use, which is what
makes `ImageCrop[ImagePad[image, m], ImageDimensions[image]]` exactly the original image. A
crop may not enlarge.

`ImageCrop[image]` with no size trims a uniform border, asking how much of the frame carries no
information. The border colour is read from a corner rather than assumed black, so a scanned
page's white margin is trimmed too. An entirely uniform image comes back unchanged — there is
no content to keep and a zero-sized image is not an image.
