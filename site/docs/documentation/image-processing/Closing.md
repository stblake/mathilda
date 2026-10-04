# Closing

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Closing[image, r] dilates then erodes with the same element, filling dark features smaller than it. Idempotent, like Opening, and the two bracket the image: Erosion <= Opening <= image <= Closing <= Dilation pointwise everywhere A "Bit" image stays "Bit"; other types give "Real".`**

## Examples (41)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (7)

```mathematica
In[1]:= ImageData[Closing[Image[{{1., 1., 1.}, {1., 0., 1.}, {1., 1., 1.}}], 1]]
Out[1]= {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}}

In[2]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[4]:= Closing[disk, 1]
Out[4]= -Image-

In[5]:= Closing[disk, 2]
Out[5]= -Image-

In[6]:= ImageDimensions[Closing[disk, 2]]
Out[6]= {16, 16}

In[7]:= Closing[bit, 1]
Out[7]= -Image-
```

### Scope (19)

```mathematica
In[8]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[9]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[10]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[12]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[13]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[14]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[15]:= volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[16]:= Closing[chk, 1]
Out[16]= -Image-

In[17]:= Closing[ramp, 1]
Out[17]= -Image-

In[18]:= Closing[noise, 2]
Out[18]= -Image-

In[19]:= Closing[rgb, 1]
Out[19]= -Image-

In[20]:= Closing[byte, 1]
Out[20]= -Image-

In[21]:= Closing[vol, 1]
Out[21]= -Image-

In[22]:= Closing[volb, 1]
Out[22]= -Image-

In[23]:= Closing[disk, 3]
Out[23]= -Image-

In[24]:= Closing[disk, 4]
Out[24]= -Image-

In[25]:= ImageChannels[Closing[rgb, 1]]
Out[25]= 3

In[26]:= ImageDimensions[Closing[vol, 2]]
Out[26]= {12, 10, 8}
```

### Applications (4)

```mathematica
In[27]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[28]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[29]:= Binarize[Closing[noise, 1]]
Out[29]= -Image-

In[30]:= EdgeDetect[Closing[disk, 1]]
Out[30]= -Image-
```

### Properties & Relations (6)

```mathematica
In[31]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[32]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[33]:= ImageData[Closing[disk, 0]] === ImageData[disk]
Out[33]= True

In[34]:= ImageDimensions[Closing[disk, 3]] === ImageDimensions[disk]
Out[34]= True

In[35]:= ImageData[Closing[Closing[disk, 1], 1]] === ImageData[Closing[disk, 2]]
Out[35]= True

In[36]:= Max[Flatten[ImageData[Closing[bit, 1]]]] <= 1.0
Out[36]= True
```

### Neat Examples (2)

```mathematica
In[37]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[38]:= Closing[zone, 2]
Out[38]= -Image-
```

### Applications (3)

```mathematica
In[39]:= ImageData[Closing[Image[{{1, 1, 1}, {1, 0, 1}, {1, 1, 1}}, "Bit"], 1]]
Out[39]= {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}}

In[40]:= ImageType[Closing[Image[{{1, 1, 1}, {1, 0, 1}, {1, 1, 1}}, "Bit"], 1]]
Out[40]= "Bit"

In[41]:= ImageData[Closing[Closing[Image[{{1, 1, 1}, {1, 0, 1}, {1, 1, 1}}, "Bit"], 1], 1]]
Out[41]= {{1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}}
```

## Implementation notes

**Algorithm.** `builtin_closing` is `morph_builtin(res, MORPH_DILATE, two_pass =
true)` — a **dilation followed by an erosion** with the *same* structuring
element. The shared `morph_builtin` runs `morph_run` once for the first pass and,
when `two_pass`, again with the opposite operation (`MORPH_DILATE` then
`MORPH_ERODE`), reusing the same support mask both times. Using the same element
both times is what makes the pair **idempotent** — `Closing[Closing[f]]` equals
`Closing[f]` — and what makes it a closing rather than merely two smoothings; a
different second element would still smooth but would not be a closing.
Geometrically it fills dark features smaller than the element while leaving larger
ones close to their original size. Each pass inherits `Dilation`'s machinery: the
separable van Herk–Gil-Werman max/min for a full rectangle, `morph_direct` for an
arbitrary element, and replicate padding. Because dilation and erosion are duals
and the padding is self-dual, the family brackets the image:
`Erosion <= Opening <= image <= Closing <= Dilation` pointwise everywhere.

**Data structures.** Two flat `double` scratch buffers `a` (first pass) and `b`
(second pass), plus the `unsigned char` support mask. A `"Bit"` image stays
`"Bit"` (`image_build_typed`), every other type gives `"Real"`.

**Complexity / limits.** Twice a single morphology pass: `O(width · height ·
channels)` amortised for a full rectangle (independent of `r`), `O(width ·
height · channels · |support|)` for an arbitrary element. A volume takes the
rank-3 two-pass path.

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- J. Serra, *Image Analysis and Mathematical Morphology* (Academic Press, 1982).
- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`Closing[image, r]` dilates then erodes with the same element, filling dark
features smaller than it. Using the **same** structuring element for both passes
is what makes it idempotent — `Closing[Closing[f]] = Closing[f]` — and what makes
it a closing rather than two unrelated smoothings.

Like `Opening`, it is idempotent, and the two bracket the image:
`Erosion <= Opening <= image <= Closing <= Dilation` pointwise everywhere. Each
pass reuses `Dilation`/`Erosion`'s separable van Herk–Gil-Werman machinery, so
the cost does not grow with the radius. A `"Bit"` image stays `"Bit"`; other
types give `"Real"`.
