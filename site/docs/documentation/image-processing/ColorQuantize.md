# ColorQuantize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ColorQuantize[image, n] reduces the image to at most n colours by MEDIAN CUT: the box with the widest single-channel spread is split at its median until n boxes remain, and each collapses to its mean colour. Widest spread rather than most pixels, since a large box of nearly identical colours does not need splitting and a small one spanning half the spectrum does. Median cut rather than k-means because it is DETERMINISTIC -- a palette that depended on the random stream could not be tested or documented. The channel count is preserved and alpha passes through.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= ramp = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[2]:= ColorQuantize[ramp, 4]
Out[2]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA+ElEQVR42u3Z0Q2EMBDEUEB0kg5SQ/pvCSrgDwjsPDdwAmtWFreOMY4FS+99yu9uXv1cCCCAABCQy6523qG1ZgFOEAggAASooMDasQAnCAQQAAJUUGbtWIATBAIIAAEqKLN2LMAJAgEEgAAVlFk7FuAEgQACQIAKyqydq+e1ACeIABBAAKpWUNXaueu5LMAJIgAEEIC/V9DTtfO1bzgW4ASBAAJAQNUK8o+VBThBIIAAEFC7gtSOBThBIIAAEOBbkNqxACcIBBAAAlSQ2rEAJwgEEAACVJDasQAnCAQQAAJyKkjtWAABIIAAEFC7gtSOBRAAAggAAZGcz88VHrVjc9UAAAAASUVORK5CYII=)

```mathematica
In[3]:= Table[Length[Union[Flatten[ImageData[ColorQuantize[ramp, n]]]]], {n, 1, 4}]
Out[3]= {1, 2, 3, 4}
```

```mathematica
In[4]:= ColorQuantize[Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"], 6]
Out[4]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAu0lEQVR42u3U0QnAIAxAwVgcxE/HcBNX6ajdxC7RQAP3BoiQI7Y57xOZrdTpMXaU7goBACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgD6vr5H8wlN6vAvwBQkAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAV69sOXAAAAQAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAQAD+3gt/JwXRJQKAGQAAAABJRU5ErkJggg==)

### Properties & Relations (4)

```mathematica
In[5]:= ramp = Image[Table[N[(i + j)/32], {i, 1, 16}, {j, 1, 16}], "Real"];
```

The same input gives the same palette, every time

```mathematica
In[6]:= ImageData[ColorQuantize[ramp, 5]] === ImageData[ColorQuantize[ramp, 5]]
Out[6]= True
```

```mathematica
In[7]:= ImageDimensions[ColorQuantize[ramp, 4]] === ImageDimensions[ramp]
Out[7]= True
```

More colours than the image holds cannot invent any

```mathematica
In[8]:= Length[Union[Flatten[ImageData[ColorQuantize[Image[{{0., 1.}, {0., 1.}}, "Real"], 8]]]]] <= 2
Out[8]= True
```

## Algorithm

imagecolor.c -- ColorReplace, ColorQuantize and HistogramTransform.

Three heads that act on an image's COLOURS rather than its geometry, and they share the one thing that makes such operations awkward: a decision made per pixel needs a global view first. Replacing a colour needs a distance rule, quantising needs a palette derived from every pixel, and equalising needs the whole distribution. So each of these makes a pass to gather, then a pass to write — which is why none of them fits the filter machinery in imagefilter.c.

## Implementation notes

- `Protected`. **Median cut**: the box with the widest single-channel spread is split at its
  median until `n` boxes remain, and each collapses to its mean colour. Widest spread rather
  than most pixels — a large box of nearly identical colours does not need splitting, and a
  small one spanning half the spectrum does.
- Median cut rather than k-means because it is **deterministic**: a palette that depended on the
  random stream could be neither tested nor documented.
- Channel count and dimensions are preserved; alpha passes through.
- Asking for more colours than the image holds cannot invent any.

**Attributes:** `Protected`.

## References

- Source: [`src/imagecolor.c`](https://github.com/stblake/mathilda/blob/main/src/imagecolor.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
