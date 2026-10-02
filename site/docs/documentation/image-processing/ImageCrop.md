# ImageCrop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageCrop[image, {w, h}] crops to w x h about the centre, any odd remainder going to the right and bottom -- the same floor-division convention the kernel centres use, which is what makes ImageCrop[ImagePad[image, m], ImageDimensions[image]] exactly the original image. A crop may not enlarge. ImageCrop[image] instead TRIMS A UNIFORM BORDER, asking how much of the frame carries no information; the border colour is read from a corner rather than assumed black, since a scanned page's margin is white. An entirely uniform image comes back unchanged, there being no content to keep and a zero-sized image not being one.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ImagePad](../../image-processing/ImagePad/), [Image3D](../../image-processing/Image3D/), [ImageConvolve](../../image-processing/ImageConvolve/)

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
