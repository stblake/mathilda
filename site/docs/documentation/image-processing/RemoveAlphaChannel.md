# RemoveAlphaChannel

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RemoveAlphaChannel[image] drops the alpha channel. RemoveAlphaChannel[image, b] instead COMPOSITES over a background of brightness b, which is the difference between forgetting the transparency and resolving it: a half-transparent white pixel over black is grey, where dropping alpha would leave it white.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= t = SetAlphaChannel[Image[{{1.0, 1.0}, {1.0, 1.0}}, "Real"], 0.5];

In[2]:= ImageChannels[RemoveAlphaChannel[t]]
Out[2]= 1

In[3]:= Round[Max[Flatten[ImageData[RemoveAlphaChannel[t]]]], 0.001]
Out[3]= 1.0

In[4]:= Round[Max[Flatten[ImageData[RemoveAlphaChannel[t, 0.]]]], 0.001]
Out[4]= 0.5
```

## Implementation notes

- `Protected`.
- The two forms are genuinely different: a half-transparent white pixel over black is **grey**,
  where dropping alpha leaves it white. One forgets the transparency; the other resolves it.

**Attributes:** `Protected`.

## References

- Source: [`src/imagecompose.c`](https://github.com/stblake/mathilda/blob/main/src/imagecompose.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
