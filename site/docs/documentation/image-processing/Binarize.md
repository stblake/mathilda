# Binarize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Binarize[image] thresholds image by Otsu's method (see FindThreshold), giving a "Bit" image. Binarize[image, t] thresholds at t. A pixel STRICTLY ABOVE the threshold becomes 1, so a pixel exactly at it becomes 0 -- which matters, because "above" and "at or above" differ on exactly the pixels a threshold was chosen to sit between. A colour image is reduced to luminance first.`**

## Examples (41)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (5)

```mathematica
In[1]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= Binarize[chk]
Out[3]= -Image-

In[4]:= Binarize[disk]
Out[4]= -Image-

In[5]:= ImageDimensions[Binarize[chk]]
Out[5]= {16, 16}
```

### Scope (22)

```mathematica
In[6]:= ImageData[Binarize[Image[{{0.4, 0.5, 0.6}}], 0.5], "Bit"]
Out[6]= {{0, 0, 1}}

In[7]:= ImageData[Binarize[Image[{{0.2,0.2,0.8,0.8},{0.2,0.2,0.8,0.8}}]], "Bit"]
Out[7]= {{0, 0, 1, 1}, {0, 0, 1, 1}}

In[8]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[9]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[10]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[11]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[13]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[14]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[15]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[16]:= volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[17]:= Binarize[ramp]
Out[17]= -Image-

In[18]:= Binarize[zone]
Out[18]= -Image-

In[19]:= Binarize[noise]
Out[19]= -Image-

In[20]:= Binarize[rgb]
Out[20]= -Image-

In[21]:= Binarize[sky]
Out[21]= -Image-

In[22]:= Binarize[bit]
Out[22]= -Image-

In[23]:= Binarize[byte]
Out[23]= -Image-

In[24]:= Binarize[vol]
Out[24]= -Image-

In[25]:= Binarize[volb]
Out[25]= -Image-

In[26]:= ImageChannels[Binarize[rgb]]
Out[26]= 1

In[27]:= ImageDimensions[Binarize[vol]]
Out[27]= {12, 10, 8}
```

### Applications (4)

```mathematica
In[28]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[29]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[30]:= Binarize[Binarize[zone]]
Out[30]= -Image-

In[31]:= Dilation[Binarize[disk], 1]
Out[31]= -Image-
```

### Properties & Relations (4)

```mathematica
In[32]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[33]:= ImageDimensions[Binarize[chk]] === ImageDimensions[chk]
Out[33]= True

In[34]:= Max[Flatten[ImageData[Binarize[chk]]]] <= 1.0
Out[34]= True

In[35]:= Min[Flatten[ImageData[Binarize[chk]]]] >= 0.0
Out[35]= True
```

### Neat Examples (2)

```mathematica
In[36]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[37]:= Binarize[zone]
Out[37]= -Image-
```

### Applications (4)

```mathematica
In[38]:= Binarize[Image[{{0.1, 0.9}, {0.8, 0.2}}]]
Out[38]= -Image-

In[39]:= ImageData[Binarize[Image[{{0.1, 0.9}, {0.8, 0.2}}]]]
Out[39]= {{0.0, 1.0}, {1.0, 0.0}}

In[40]:= ImageData[Binarize[Image[{{0.1, 0.4, 0.9}}], 0.5]]
Out[40]= {{0.0, 0.0, 1.0}}

In[41]:= ImageType[Binarize[Image[{{0.2, 0.9}, {0.8, 0.1}}]]]
Out[41]= "Bit"
```

## Implementation notes

**Algorithm.** `builtin_binarize` thresholds an image into a `"Bit"` image.
`Binarize[image]` chooses the threshold by Otsu's method (`img_otsu`, see
`FindThreshold`); `Binarize[image, t]` uses the explicit level `t`. A colour
image is reduced to a Rec. 601 luminance plane first (`img_grey_plane`). The
decision is **strictly above**: a pixel with value `> t` becomes `1`, so a pixel
exactly *at* `t` becomes `0`. That boundary matters because "above" and "at or
above" differ on exactly the pixels a threshold was chosen to sit between, and a
test pins it. The mask is built by `bit_image_from_mask` (shared with
`LocalAdaptiveBinarize`, so the two cannot disagree on the output type) into an
`Image[..., "Bit"]`, typed `"Bit"` rather than `"Real"` because the result is
binary *by construction* — a later `ImageData` scales nothing and the image
reports its binary nature. A volume takes the rank-3 path (`binarize3_run`). With
the default (Otsu) threshold, a degenerate single-value image declines, since no
threshold splits one cluster.

**Data structures.** A flat luminance `double` plane and an `unsigned char` 0/1
mask of `width · height`, packed into a `"Bit"` image buffer
(`image_build_bit`). No multi-channel buffer survives — binarisation collapses to
one plane.

**Complexity / limits.** `O(width · height)` for the threshold comparison, plus
`O(width · height + bins)` for an Otsu default. With an explicit `t` it never
declines; with the default it declines on a single-value image.

**Attributes:** `Protected`.

## References

- N. Otsu, *A threshold selection method from gray-level histograms*, IEEE Trans. Systems, Man, and Cybernetics **9** (1979) 62-66.
- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`Binarize[image]` thresholds by Otsu's method (see `FindThreshold`), giving a
`"Bit"` image; `Binarize[image, t]` thresholds at the explicit level `t`.

A pixel **strictly above** the threshold becomes `1`, so a pixel exactly at it
becomes `0` — which matters, because "above" and "at or above" differ on exactly
the pixels a threshold was chosen to sit between. A colour image is reduced to
luminance first. The result is typed `"Bit"` rather than `"Real"` because it is
binary by construction: `ImageData` then scales nothing and the image still
reports itself as binary. With the default (Otsu) threshold a single-value image
declines, since no threshold splits one cluster.
