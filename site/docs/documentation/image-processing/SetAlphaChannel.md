# SetAlphaChannel

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SetAlphaChannel[image] attaches a fully opaque alpha channel. SetAlphaChannel[image, a] sets one opacity everywhere when a is a number in [0, 1], or per pixel when a is an image of the same dimensions (read as grey, so a colour mask is not taken as its red channel alone). A mask of the wrong size is declined rather than resampled.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (5)

```mathematica
In[1]:= a = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= ImageChannels[SetAlphaChannel[a, 0.5]]
Out[2]= 2

In[3]:= ImageChannels[SetAlphaChannel[Image[Table[{0.5, 0.2, 0.9}, {i, 1, 8}, {j, 1, 8}], "Real"], 0.5]]
Out[3]= 4

In[4]:= SetAlphaChannel[a, Image[Table[N[j/16], {i, 1, 16}, {j, 1, 16}], "Real"]]
Out[4]= -Image-

In[5]:= Head[SetAlphaChannel[a, Image[{{0.5}}, "Real"]]]
Out[5]= SetAlphaChannel
```

### Properties & Relations (3)

```mathematica
In[6]:= a = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];
```

A colour mask is averaged, not read as its red channel: {0.2, 0.4, 0.6} gives 0.4

```mathematica
In[7]:= Module[{m = Image[Table[{0.2, 0.4, 0.6}, {i, 1, 16}, {j, 1, 16}], "Real"], u}, u = Union[Flatten[ImageData[AlphaChannel[SetAlphaChannel[a, m]]]]]; Round[First[u], 0.0001]]
Out[7]= 0.4
```

Setting then removing gets back the channel count it started with

```mathematica
In[8]:= ImageChannels[RemoveAlphaChannel[SetAlphaChannel[a, 0.5]]] === ImageChannels[a]
Out[8]= True
```

## Implementation notes

- `Protected`.
- A mask is read as **grey** (its channels averaged), so a colour mask is not silently taken as
  its red channel alone.
- A mask of the wrong size is declined rather than resampled: a mismatch is a mistake, not a
  request to interpolate.

**Attributes:** `Protected`.

## References

- Source: [`src/imagecompose.c`](https://github.com/stblake/mathilda/blob/main/src/imagecompose.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
