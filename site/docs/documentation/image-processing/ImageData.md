# ImageData

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageData[image] gives the pixel array as reals in [0, 1], scaling out the image's type -- a "Byte" 255 comes back as exactly 1.0. The array is height x width, or height x width x channels for a colour image, interleaved. ImageData[image, type] gives the stored values unscaled instead, where type must be the image's own type; converting between types is a separate operation with its own rounding, not something this does silently.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= ImageData[Image[{{0, 128, 255}}]]
Out[1]= {{0.0, 0.501961, 1.0}}

In[2]:= ImageData[Image[{{0, 128, 255}}], "Byte"]
Out[2]= {{0, 128, 255}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [FullForm](../../expression-information/FullForm/)

- Source: [`src/image.c`](https://github.com/stblake/mathilda/blob/main/src/image.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
