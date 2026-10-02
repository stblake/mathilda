# ImagePad

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImagePad[image, m] pads m pixels on every side; ImagePad[image, {{left, right}, {bottom, top}}] pads each side separately, in Mathematica's VISUAL order -- so `top` adds rows at the start of the data, since row 1 is the top of the image. Negative amounts crop, but may not erase the image. ImagePad[image, m, v] fills with the value v (default 0); ImagePad[image, m, "Fixed"] replicates the edge pixel, the same boundary rule the filters use, so padding then filtering composes with it; ImagePad[image, m, "Reflected"] mirrors WITHOUT repeating the edge -- {1,2,3} padded by 1 gives {2,1,2,3,2}, not {1,1,2,3,3}, because doubling the edge sample biases any later average toward the border. Reflection uses a period of 2n-2, so padding deeper than the image still works.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ImageCrop](../../image-processing/ImageCrop/), [Image3D](../../image-processing/Image3D/), [ImageConvolve](../../image-processing/ImageConvolve/)

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
