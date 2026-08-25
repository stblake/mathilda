# HistogramTransform

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HistogramTransform[image] equalises the histogram, spreading the brightness distribution toward uniform over 256 bins by mapping each value through the cumulative distribution. The mapping is computed from the LUMINANCE and applied to every channel as a ratio, so hue survives; equalising each channel independently would shift colour, since it removes exactly the imbalance that makes an image warm or cool. A black pixel has no ratio to scale and takes the new luminance in every channel. Alpha passes through.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= dark = Image[Table[N[(i + j)/64], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= HistogramTransform[dark]
Out[2]= -Image-
```

Nothing in the original is above 0.5; afterwards the range is covered

```mathematica
In[3]:= {Round[Max[Flatten[ImageData[dark]]], 0.01], Round[Max[Flatten[ImageData[HistogramTransform[dark]]]], 0.01], Round[Min[Flatten[ImageData[HistogramTransform[dark]]]], 0.01]}
Out[3]= {0.5, 1.0, 0.0}
```

```mathematica
In[4]:= HistogramTransform[Image[Table[{N[i/32], N[j/32], 0.25}, {i, 1, 16}, {j, 1, 16}], "Real"]]
Out[4]= -Image-
```

### Properties & Relations (4)

```mathematica
In[5]:= dark = Image[Table[N[(i + j)/64], {i, 1, 16}, {j, 1, 16}], "Real"];
```

Monotone: equalisation is a cumulative distribution, so it cannot reorder pixels

```mathematica
In[6]:= Module[{d = ImageData[HistogramTransform[dark]]}, d[[1, 1]] <= d[[8, 8]] && d[[8, 8]] <= d[[16, 16]]]
Out[6]= True
```

```mathematica
In[7]:= ImageDimensions[HistogramTransform[dark]] === ImageDimensions[dark]
Out[7]= True
```

An already-spread image is near a fixed point, which is what "toward uniform" means

```mathematica
In[8]:= Module[{a = HistogramTransform[dark]}, Max[Abs[Flatten[ImageData[HistogramTransform[a]] - ImageData[a]]]] < 0.2]
Out[8]= True
```

## Algorithm

imagecolor.c -- ColorReplace, ColorQuantize and HistogramTransform.

Three heads that act on an image's COLOURS rather than its geometry, and they share the one thing that makes such operations awkward: a decision made per pixel needs a global view first. Replacing a colour needs a distance rule, quantising needs a palette derived from every pixel, and equalising needs the whole distribution. So each of these makes a pass to gather, then a pass to write — which is why none of them fits the filter machinery in imagefilter.c.

## Implementation notes

- `Protected`. Each value is mapped through the cumulative distribution over 256 bins, spreading
  the histogram toward uniform.
- The mapping is computed from the **luminance** and applied to every channel as a ratio, so hue
  survives. Equalising each channel independently would shift colour — it removes exactly the
  imbalance that makes an image warm or cool.
- A black pixel has no ratio to scale and takes the new luminance in every channel.
- **Monotone**: a cumulative distribution can never reorder two pixels.
- Alpha passes through.

**Attributes:** `Protected`.

## References

- Source: [`src/imagecolor.c`](https://github.com/stblake/mathilda/blob/main/src/imagecolor.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
