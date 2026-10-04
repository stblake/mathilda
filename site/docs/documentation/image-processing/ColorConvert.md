# ColorConvert

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ColorConvert[image, "Grayscale"] (or "Gray") reduces an image or an Image3D to a single channel using the Rec. 601 luminance weights 0.299 R + 0.587 G + 0.114 B, the same weights every filter here uses when it needs brightness. An image that is ALREADY GREY is returned unchanged, bit for bit, since no weighting happens. An image whose three channels are merely EQUAL is returned only to within an ulp, and whether it is exact depends on the value: those weights sum to 0.9999999999999999 when added in the order they are applied, though to exactly 1.0 in any order beginning with 0.114, so the final rounding lands on the input for some values and one ulp below it for others. The weights are the standard's and are not adjusted to compensate; a triple hand-tuned to sum to exactly 1.0 in double would no longer be Rec. 601.`**

## Examples (33)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (6)

```mathematica
In[1]:= ImageData[ColorConvert[Image[{{{1., 0., 0.}, {0., 1., 0.}}}], "Grayscale"]]
Out[1]= {{0.299, 0.587}}

In[2]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[4]:= ColorConvert[rgb, "Grayscale"]
Out[4]= -Image-

In[5]:= ImageChannels[ColorConvert[rgb, "Grayscale"]]
Out[5]= 1

In[6]:= ColorConvert[sky, "Grayscale"]
Out[6]= -Image-
```

### Scope (13)

```mathematica
In[7]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[8]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[9]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[10]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[11]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[12]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[13]:= ColorConvert[chk, "Grayscale"]
Out[13]= -Image-

In[14]:= ColorConvert[zone, "Grayscale"]
Out[14]= -Image-

In[15]:= ColorConvert[byte, "Grayscale"]
Out[15]= -Image-

In[16]:= ColorConvert[vol, "Grayscale"]
Out[16]= -Image-

In[17]:= ColorConvert[rgb, "Gray"]
Out[17]= -Image-

In[18]:= ImageChannels[ColorConvert[sky, "Grayscale"]]
Out[18]= 1

In[19]:= ImageDimensions[ColorConvert[rgb, "Grayscale"]]
Out[19]= {16, 16}
```

### Applications (4)

```mathematica
In[20]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[21]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[22]:= EdgeDetect[ColorConvert[sky, "Grayscale"]]
Out[22]= -Image-

In[23]:= Binarize[ColorConvert[rgb, "Grayscale"]]
Out[23]= -Image-
```

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

In[30]:= ColorConvert[zone, "Grayscale"]
Out[30]= -Image-
```

### Applications (3)

Pure red -> its Rec. 601 luminance 0.299

```mathematica
In[31]:= ImageData[ColorConvert[Image[{{{1., 0., 0.}}}], "Grayscale"]]
Out[31]= {{0.299}}
```

Green carries most of the perceived brightness

```mathematica
In[32]:= ImageData[ColorConvert[Image[{{{0., 1., 0.}}}], "Grayscale"]]
Out[32]= {{0.587}}
```

The result is one channel

```mathematica
In[33]:= ImageChannels[ColorConvert[Image[{{{1., 1., 0.}, {0., 0., 1.}}}], "Gray"]]
Out[33]= 1
```

## Implementation notes

**Algorithm.** `builtin_colorconvert` converts an image (plane or `Image3D`) to greyscale.
Only `"Grayscale"` (or `"Gray"`) is accepted: the other colour spaces Mathematica supports
(LAB, HSB, XYZ, …) each carry their own white point and transfer-function decisions, and
accepting the name while doing something approximate would be worse than declining it. The
reduction uses the **Rec. 601** luminance weights `0.299 R + 0.587 G + 0.114 B` (`img_to_grey`
for a plane, `img3_grey_volume` for a volume — one place, so both ranks agree), because the eye
is not equally sensitive across the spectrum: green carries most perceived brightness and blue
almost none, so a plain average would put pure red and pure blue on the same side of a
threshold when perceptually they are far apart.

What is exact is documented honestly. An **already-grey** image is copied through bit for bit.
An image whose three channels are merely **equal** is exact only to within an ulp, and *whether*
depends on the value: those weights sum to `0.9999999999999999` in the order they are applied
but exactly `1.0` in any order beginning with `0.114`, so the final rounding lands on the input
for some values and one ulp below for others. The weights are the standard's and are not
adjusted to compensate — a triple hand-tuned to sum to exactly `1.0` in double would no longer
be Rec. 601.

**Data structures.** One decoded unit buffer in; a fresh single-channel buffer out, wrapped by
`image_build_real` / `image3d_build_real` as a packed `"Real"` image (every image head returns
a packed buffer — `make check-image-packing`).

**Complexity / limits.** `O(pixels)`. Only the greyscale target is implemented.

**Attributes:** `Protected`.

## References

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`ColorConvert[image, "Grayscale"]` (or `"Gray"`) reduces an image or an `Image3D` to one
channel using the Rec. 601 weights `0.299 R + 0.587 G + 0.114 B` — the same weights every
filter here uses when it needs brightness. The weights are not a mean: green carries most of
the perceived brightness and blue almost none, so a plain average would make a saturated blue
and a saturated green look equally bright.

An already-grey image is returned unchanged, bit for bit. An image whose channels are merely
equal matches the input only to within an ulp, since those weights do not sum to exactly 1 in
binary; the weights are the standard's and are not adjusted to compensate. Only the greyscale
target is supported.
