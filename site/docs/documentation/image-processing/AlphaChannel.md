# AlphaChannel

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlphaChannel[image] gives the image's opacity as a one-channel image. An image with no alpha channel answers with an all-opaque one rather than declining: "how transparent is this?" has an answer for every image, and it is "not at all". Two channels are read as grey+alpha and four as RGB+alpha.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= a = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= Union[Flatten[ImageData[AlphaChannel[a]]]]
Out[2]= {1.0}

In[3]:= Union[Flatten[ImageData[AlphaChannel[SetAlphaChannel[a, 0.25]]]]]
Out[3]= {0.25}

In[4]:= ImageChannels[AlphaChannel[Image[Table[{0.2, 0.4, 0.6}, {i, 1, 8}, {j, 1, 8}], "Real"]]]
Out[4]= 1
```

## Implementation notes

- `Protected`.
- An image with **no** alpha channel answers with an all-opaque one rather than declining: "how
  transparent is this?" has an answer for every image, and it is "not at all".
- Two channels are read as grey+alpha, four as RGB+alpha.

**Attributes:** `Protected`.

## References

- Source: [`src/imagecompose.c`](https://github.com/stblake/mathilda/blob/main/src/imagecompose.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
