# ImageAssemble

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageAssemble[{{a, b}, {c, d}}] tiles a grid of images into one; ImageAssemble[{a, b}] makes a single row. Each tile keeps its natural size -- a row is as tall as its tallest tile and a column as wide as its widest, and any gap is left blank rather than stretched, since stretching would resample an image the caller did not ask to resize. Alpha survives if any tile had it.`**

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= a = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[2]:= ImageAssemble[{a, a}]
Out[2]= -Image-
```

![32x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAAAwCAYAAADuFn/PAAAArklEQVR42u3UMQrEMBADQOcI2KW/7VcnTUh9lZWFUbWlYNEcc86rPRljtAp37/2911qtcv9fk2g8IJyz4myrdP6nvwUgCEHYCfa3AAQhCDvB/haAIARhJ9jfAhCEIOwE+1sAghCEnWB/C0AQgrAT7G8BCEIQdoL9LQBBCMJOsL8FIAhB2Al2sAAEIQg7CEKQfJkg7CAIQbKZIOwgCEGymSDsIAhBspkg7CAIQZLLDZFSMriIzhEcAAAAAElFTkSuQmCC)

```mathematica
In[3]:= ImageDimensions[ImageAssemble[{a, a}]]
Out[3]= {32, 16}

In[4]:= ImageDimensions[ImageAssemble[{{a, a}, {a, a}}]]
Out[4]= {32, 32}
```

### Applications (3)

```mathematica
In[5]:= a = Image[Table[N[Boole[(i - 16)^2 + (j - 16)^2 <= 100]], {i, 1, 32}, {j, 1, 32}], "Real"];
```

A contact sheet comparing one filter at four radii

```mathematica
In[6]:= ImageAssemble[{Table[GaussianFilter[a, r], {r, 1, 2}], Table[GaussianFilter[a, r], {r, 3, 4}]}]
Out[6]= -Image-
```

![64x64 result](data:image/png;base64,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)

The same image before and after, side by side

```mathematica
In[7]:= ImageAssemble[{a, EdgeDetect[a]}]
Out[7]= -Image-
```

![64x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEAAAAAgCAYAAACinX6EAAAAx0lEQVR42u2ZSxaAIAhFwdP+t0yjRnUSBNTiMS1+N1A0JiKhwtKouAAAABSXw6MsIsTM04MWua/bo3GwZRd4chwViMX3k4/ecxcATeJZIKyJWd/vAhhJPgrC5XvEjla3ZSUfoe/qbaXetrtABDyNnbY6gMwW0uhjENq1/KN2EWZ+rUJUQGb/Z9lLBxA9za0Yl3EY+jKA3sIVuaCiArL71mNnxii+bQXM+gBt9SiaVQVaPdwH4EbI8WOk3J0gBiEAAAAAAICfyQmr6XEzGTAoegAAAABJRU5ErkJggg==)

## Implementation notes

- `Protected`.
- Each tile keeps its **natural size**: a row is as tall as its tallest tile, a column as wide as
  its widest, and any gap is left blank rather than stretched — stretching would resample an image
  the caller did not ask to resize.
- Channels are promoted as in `ImageCompose`, and alpha survives if any tile had it.

**Attributes:** `Protected`.

## References

**See also:** [ImageCompose](../../image-processing/ImageCompose/)

- Source: [`src/imagecompose.c`](https://github.com/stblake/mathilda/blob/main/src/imagecompose.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
