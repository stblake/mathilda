# Dilation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Dilation[image, r] gives the maximum over a (2r+1) x (2r+1) square neighbourhood; Dilation[image, elem] uses the SUPPORT of the matrix elem -- its nonzero positions -- as the neighbourhood. This is flat morphology: the element's values do not enter the maximum, which is what keeps Dilation[img, BoxMatrix[1]] and Dilation[img, 1] the same operation. Padding replicates the border, the same rule the convolutions use, which is what makes Dilation >= image hold at the edges too. A full rectangle is separable for the maximum exactly as for a sum, so it costs kw + kh comparisons rather than kw * kh. A "Bit" image gives a "Bit" image; other types give "Real".`**

## Examples (43)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (7)

```mathematica
In[1]:= ImageData[Dilation[Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}], 1]]
Out[1]= {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}}

In[2]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[4]:= Dilation[disk, 1]
Out[4]= -Image-

In[5]:= Dilation[disk, 2]
Out[5]= -Image-

In[6]:= ImageDimensions[Dilation[disk, 2]]
Out[6]= {16, 16}

In[7]:= Dilation[bit, 1]
Out[7]= -Image-
```

### Scope (21)

```mathematica
In[8]:= ImageData[Dilation[Image[{{0, 1, 0, 0}, {0, 0, 0, 0}}], 1], "Bit"]
Out[8]= {{1, 1, 1, 0}, {1, 1, 1, 0}}

In[9]:= ImageType[Dilation[Image[{{0, 1, 0}, {0, 0, 0}}], 1]]
Out[9]= "Bit"

In[10]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[13]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[14]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[15]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[16]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[17]:= volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[18]:= Dilation[chk, 1]
Out[18]= -Image-

In[19]:= Dilation[ramp, 1]
Out[19]= -Image-

In[20]:= Dilation[noise, 2]
Out[20]= -Image-

In[21]:= Dilation[rgb, 1]
Out[21]= -Image-

In[22]:= Dilation[byte, 1]
Out[22]= -Image-

In[23]:= Dilation[vol, 1]
Out[23]= -Image-

In[24]:= Dilation[volb, 1]
Out[24]= -Image-

In[25]:= Dilation[disk, 3]
Out[25]= -Image-

In[26]:= Dilation[disk, 4]
Out[26]= -Image-

In[27]:= ImageChannels[Dilation[rgb, 1]]
Out[27]= 3

In[28]:= ImageDimensions[Dilation[vol, 2]]
Out[28]= {12, 10, 8}
```

### Applications (4)

```mathematica
In[29]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[30]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[31]:= Binarize[Dilation[noise, 1]]
Out[31]= -Image-

In[32]:= EdgeDetect[Dilation[disk, 1]]
Out[32]= -Image-
```

### Properties & Relations (6)

```mathematica
In[33]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[34]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[35]:= ImageData[Dilation[disk, 0]] === ImageData[disk]
Out[35]= True

In[36]:= ImageDimensions[Dilation[disk, 3]] === ImageDimensions[disk]
Out[36]= True

In[37]:= ImageData[Dilation[Dilation[disk, 1], 1]] === ImageData[Dilation[disk, 2]]
Out[37]= True

In[38]:= Max[Flatten[ImageData[Dilation[bit, 1]]]] <= 1.0
Out[38]= True
```

### Neat Examples (2)

```mathematica
In[39]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[40]:= Dilation[zone, 2]
Out[40]= -Image-
```

### Applications (3)

```mathematica
In[41]:= ImageData[Dilation[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}], 1]]
Out[41]= {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}}

In[42]:= ImageData[Dilation[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}, "Bit"], {{0, 1, 0}, {1, 1, 1}, {0, 1, 0}}]]
Out[42]= {{0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 0.0}}

In[43]:= ImageType[Dilation[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}, "Bit"], 1]]
Out[43]= "Bit"
```

## Implementation notes

**Algorithm.** `builtin_dilation` is `morph_builtin(res, MORPH_DILATE, two_pass =
false)` — the pointwise **maximum** over a structuring element.
`Dilation[image, r]` uses a `(2r+1) × (2r+1)` square; `Dilation[image, elem]`
uses the **support** of the matrix `elem` — its nonzero positions
(`morph_support`), so the element's values do not enter the maximum. This is
*flat* morphology, which is what keeps `Dilation[img, BoxMatrix[1]]` and
`Dilation[img, 1]` the same operation. A full rectangle is separable for the
maximum exactly as for a sum (max over a rectangle = max over rows of the maxima
over columns), so `morph_separable` does a row pass then a column pass, and each
1-D pass uses the **van Herk–Gil-Werman** trick: block the line into runs of `k`,
keep a prefix and a suffix running maximum per block, and any width-`k` window is
`max(suffix at its start, prefix at its end)` — three comparisons per pixel
*regardless of the radius*, which is what makes morphology usable at large `r`.
The fast path is exact (max is associative and idempotent) and a test asserts it
agrees bit-for-bit with the naive `morph_direct`. An arbitrary (non-rectangular)
element falls back to `morph_direct`. Padding replicates the border, the same
rule the convolutions use, so `Dilation >= image` holds at the edges too.

**Data structures.** Flat `double` buffers; the element is an `unsigned char`
support mask plus a `full` flag. The separable path uses padded scratch lines and
`pre`/`suf` arrays sized to the longer axis plus the kernel. A `"Bit"` image
gives a `"Bit"` image (`image_build_typed`, since a max of stored 0/1 values is
itself 0/1); every other type gives `"Real"`.

**Complexity / limits.** Separable full-rectangle: `O(width · height · channels)`
amortised, **independent of `r`**. Arbitrary element: `O(width · height ·
channels · |support|)`. A volume takes the rank-3 path, with an integer radius
only (an arbitrary 3-D element is not separable).

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- J. Serra, *Image Analysis and Mathematical Morphology* (Academic Press, 1982).
- M. van Herk, *A fast algorithm for local minimum and maximum filters on rectangular and octagonal kernels*, Pattern Recognition Letters **13** (1992) 517-521.
- J. Gil and M. Werman, *Computing 2-D min, median, and max filters*, IEEE TPAMI **15** (1993) 504-507.
- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`Dilation[image, r]` gives the maximum over a `(2r+1) × (2r+1)` square;
`Dilation[image, elem]` uses the support (nonzero positions) of the matrix
`elem`. This is **flat** morphology — the element's values do not enter the
maximum — which keeps `Dilation[img, BoxMatrix[1]]` and `Dilation[img, 1]` the
same operation.

A full rectangle is separable for the maximum, and each 1-D pass runs in three
comparisons per pixel by the van Herk–Gil-Werman algorithm, so the cost does not
grow with `r`. Padding replicates the border, so `Dilation >= image` holds at the
edges too. A `"Bit"` image gives a `"Bit"` image; other types give `"Real"`.
Together with `Erosion`, `Opening` and `Closing` it satisfies
`Erosion <= Opening <= image <= Closing <= Dilation` pointwise everywhere.
