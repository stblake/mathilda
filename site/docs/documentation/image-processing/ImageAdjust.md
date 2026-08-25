# ImageAdjust

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageAdjust[image] stretches to the full range: the darkest pixel becomes exactly 0 and the brightest exactly 1. It is IDEMPOTENT, a second stretch being the identity. A constant image has no range to stretch and comes back unchanged, since dividing by zero is not the answer and mapping the single value to either end would be arbitrary. ImageAdjust[image, {c, b}] and [image, {c, b, g}] apply contrast c, brightness b and gamma g by a curve stated here rather than inferred: v' = (v - 1/2)(1 + c) + 1/2 + b, clipped to [0, 1], then raised to the power 1/g. Contrast pivots about mid-grey so it does not also shift brightness; clipping precedes gamma because a negative base has no real power. This curve is Mathilda's documented choice, not a claim of bit-compatibility with Mathematica. Accepts volumes as well as planes.`**

## Examples (36)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (6)

```mathematica
In[1]:= ImageData[ImageAdjust[Image[{{0.25, 0.5}, {0.5, 0.75}}]]]
Out[1]= {{0.0, 0.5}, {0.5, 1.0}}

In[2]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= ImageAdjust[chk]
Out[4]= -Image-

In[5]:= ImageAdjust[disk]
Out[5]= -Image-

In[6]:= ImageDimensions[ImageAdjust[chk]]
Out[6]= {16, 16}
```

### Scope (20)

```mathematica
In[7]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[8]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[9]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[10]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[12]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[13]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[14]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[15]:= volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[16]:= ImageAdjust[ramp]
Out[16]= -Image-

In[17]:= ImageAdjust[zone]
Out[17]= -Image-

In[18]:= ImageAdjust[noise]
Out[18]= -Image-

In[19]:= ImageAdjust[rgb]
Out[19]= -Image-

In[20]:= ImageAdjust[sky]
Out[20]= -Image-

In[21]:= ImageAdjust[bit]
Out[21]= -Image-

In[22]:= ImageAdjust[byte]
Out[22]= -Image-

In[23]:= ImageAdjust[vol]
Out[23]= -Image-

In[24]:= ImageAdjust[volb]
Out[24]= -Image-

In[25]:= ImageChannels[ImageAdjust[rgb]]
Out[25]= 3

In[26]:= ImageDimensions[ImageAdjust[vol]]
Out[26]= {12, 10, 8}
```

### Applications (4)

```mathematica
In[27]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[28]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[29]:= Binarize[ImageAdjust[zone]]
Out[29]= -Image-

In[30]:= Dilation[ImageAdjust[disk], 1]
Out[30]= -Image-
```

### Properties & Relations (4)

```mathematica
In[31]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[32]:= ImageDimensions[ImageAdjust[chk]] === ImageDimensions[chk]
Out[32]= True

In[33]:= Max[Flatten[ImageData[ImageAdjust[chk]]]] <= 1.0
Out[33]= True

In[34]:= Min[Flatten[ImageData[ImageAdjust[chk]]]] >= 0.0
Out[34]= True
```

### Neat Examples (2)

```mathematica
In[35]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[36]:= ImageAdjust[zone]
Out[36]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ImageLevels](../../image-processing/ImageLevels/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
