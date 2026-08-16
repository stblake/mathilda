# ColorConvert

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ColorConvert[image, "Grayscale"] (or "Gray") reduces an image or an Image3D to a single channel using the Rec. 601 luminance weights 0.299 R + 0.587 G + 0.114 B, the same weights every filter here uses when it needs brightness. An image that is ALREADY GREY is returned unchanged, bit for bit, since no weighting happens. An image whose three channels are merely EQUAL is returned only to within an ulp, and whether it is exact depends on the value: those weights sum to 0.9999999999999999 when added in the order they are applied, though to exactly 1.0 in any order beginning with 0.114, so the final rounding lands on the input for some values and one ulp below it for others. The weights are the standard's and are not adjusted to compensate; a triple hand-tuned to sum to exactly 1.0 in double would no longer be Rec. 601.`**

## Examples (30)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (6)

```mathematica
In[1]:= ImageData[ColorConvert[Image[{{{1., 0., 0.}, {0., 1., 0.}}}], "Grayscale"]]
Out[1]= {{0.299, 0.587}}

In[2]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];
```

```mathematica
In[4]:= ColorConvert[rgb, "Grayscale"]
Out[4]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABHUlEQVR42u3dQQrDMAwFUaXYC/sOvv8xmwtkEyjIrt4sU0hLB32ErKbXWusbD4wx4s31OWecfJ+s9/0EUiGAAAJAQF1aa+3xhd77q+u/uk/W9azvQQWIIAJAAAHI6oJ0I7ndnQoQQQSAAAKwWxdUrRvJ6u5UgAgiAAQQgFNmQad3Nbt9fhUggggAAQTg9L0gJ2j2gkQQCCAABNTcC3KCFmZBIggEEAAC7AU5QbMXJIJAAAEgwF6QE7QwCxJBIIAAEOCX8mZNKkAEgQACQICnJpo1qQARBAIIIAD2gurOmlSACCIABBAAz46ueYKmAkQQASCAAHhqYs0TNBUggggAAQTAXlD4N1WIIAJAAAH4o70gsyZ7QSIIBBAAArbkBhyIET1hNe4tAAAAAElFTkSuQmCC)

```mathematica
In[5]:= ImageChannels[ColorConvert[rgb, "Grayscale"]]
Out[5]= 1
```

```mathematica
In[6]:= ColorConvert[sky, "Grayscale"]
Out[6]= -Image-
```

![24x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABACAYAAADlNHIOAAAAkklEQVR42u3cwQkAIQxFQRX7ryltJRctI4IzJfj4hL3sjIgzaLM8gQAC0GdnplewAAEQwA3AAgRAADcACxAAAT65AVXlFSxAAATwHYAFCIAAbgAWIAACuAFYgAAI4AZgAQIggBuABQiAAG4AFiAAArgBWIAACOAGYAECIIAbgAUIgABuABYgAAL4VwQWIAACvO0CBngvqKaFg9kAAAAASUVORK5CYII=)

### Scope (13)

```mathematica
In[7]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[8]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[9]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[10]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[11]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[12]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];
```

```mathematica
In[13]:= ColorConvert[chk, "Grayscale"]
Out[13]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzUlEQVR42u3aMQ6EMBAEwV7E/79sXkDgxARUxyeS0m1gzVSrjdba+nkzk++/d6VPAwAAgAAAEAAAAgBAAP7V7W3n2+/7BzhBAAQAgAAAEAAAAgBAB5vsgrILcoIEAIAAABAAAAIAQACyC8rbjl2QEyQAAAQAgAAAEAAAApBdkLcduyAnSAAACAAAAQAgAAAEILugvB3ZBTlBAgBAAAAIAAABACAA2QX5vn+AEyQAAAQAgAAAEAAAApBdUHZBcoIACAAAAQAgAAAEAIAO9QC292S1XLLeOgAAAABJRU5ErkJggg==)

```mathematica
In[14]:= ColorConvert[zone, "Grayscale"]
Out[14]= -Image-
```

![32x32 result](data:image/png;base64,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)

```mathematica
In[15]:= ColorConvert[byte, "Grayscale"]
Out[15]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA4klEQVR42u3YMQqAQAwF0ShCqnS5/1X1AjYWIm7elCkd9jO4dfcZN2RmuL9/3wOfQgABBICAuRxVpUZUkAkCAQSAABWkUlSQCQIBBIAAFeSugkwQCCAABKggdxVkgkAAASBABbmrIBMEAggAASrIXQWZIBBAAAhQQe4qyASBAAJAgApyV0EmCAQQAAJUkLsKMkEggAAQoILcvQATBAIIAAEqyN0LMEEE+AQEEAAVtMz96ff0AkwQASCAAKigNarGvyATBAIIAAEqaGbVqCATBAIIAAEqaGbVqCATBAIIAAE/4QJSnxhav9FdIQAAAABJRU5ErkJggg==)

```mathematica
In[16]:= ColorConvert[vol, "Grayscale"]
Out[16]= -Image-
```

![12x10x8 volume](data:image/png;base64,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)

```mathematica
In[17]:= ColorConvert[rgb, "Gray"]
Out[17]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABHUlEQVR42u3dQQrDMAwFUaXYC/sOvv8xmwtkEyjIrt4sU0hLB32ErKbXWusbD4wx4s31OWecfJ+s9/0EUiGAAAJAQF1aa+3xhd77q+u/uk/W9azvQQWIIAJAAAHI6oJ0I7ndnQoQQQSAAAKwWxdUrRvJ6u5UgAgiAAQQgFNmQad3Nbt9fhUggggAAQTg9L0gJ2j2gkQQCCAABNTcC3KCFmZBIggEEAAC7AU5QbMXJIJAAAEgwF6QE7QwCxJBIIAAEOCX8mZNKkAEgQACQICnJpo1qQARBAIIIAD2gurOmlSACCIABBAAz46ueYKmAkQQASCAAHhqYs0TNBUggggAAQTAXlD4N1WIIAJAAAH4o70gsyZ7QSIIBBAAArbkBhyIET1hNe4tAAAAAElFTkSuQmCC)

```mathematica
In[18]:= ImageChannels[ColorConvert[sky, "Grayscale"]]
Out[18]= 1

In[19]:= ImageDimensions[ColorConvert[rgb, "Grayscale"]]
Out[19]= {16, 16}
```

### Applications (4)

```mathematica
In[20]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[21]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];
```

```mathematica
In[22]:= EdgeDetect[ColorConvert[sky, "Grayscale"]]
Out[22]= -Image-
```

![24x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABACAYAAADlNHIOAAAAh0lEQVR42u3ZwQ0AIQwDQQfRf8uhC4LEbAk38iNHJelorOUTAAAgAAAEAIAAABAAALrX7vYryAIACMCnVbwHWAAAAXAHyAIACIA7QBYAQAAACAAAAQAgAAAEAIAAABCAeBOWBQAQgHgTlgUAEAB3gCwAgAC4A2QBAAQAgAAAEAAAAgBAAN7uADQkEXNZbOpXAAAAAElFTkSuQmCC)

```mathematica
In[23]:= Binarize[ColorConvert[rgb, "Grayscale"]]
Out[23]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzklEQVR42u3YwQ3CQBAEwQaRf8pHBDyNLU11CFday5pXddLPzrn2ed6e+N4AAAAgALt9PMF//nZcgE+QAAAQAAACAEAAAAhAtqBWNh8X4BMkAAAEAIAAABAAAAKQLaiVzccFABAAAAIAQAAACAAAAcgW1Mrm4wIACAAAAQAgAAAEAIAAZAtqZfNxAQAEAIAAABAAAAIAQACyBbWy+bgAAAIAQAAACAAAAQAgAA1vQWubjwsAIAAABACAAAAQAAAC0MAWZPNxAQAEAIAAPLIvVPYctP+7jrcAAAAASUVORK5CYII=)

### Properties & Relations (5)

```mathematica
In[24]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[25]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[26]:= ImageData[ColorConvert[chk, "Grayscale"]] === ImageData[chk]
Out[26]= True

In[27]:= ImageData[ColorConvert[rgb, "Gray"]] === ImageData[ColorConvert[rgb, "Grayscale"]]
Out[27]= True

In[28]:= ImageChannels[ColorConvert[rgb, "Grayscale"]] === 1
Out[28]= True
```

### Neat Examples (2)

```mathematica
In[29]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];
```

```mathematica
In[30]:= ColorConvert[zone, "Grayscale"]
Out[30]= -Image-
```

![32x32 result](data:image/png;base64,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)

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
