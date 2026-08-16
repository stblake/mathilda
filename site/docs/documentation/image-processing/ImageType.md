# ImageType

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageType[image] gives the pixel type as "Bit", "Byte" or "Real". The type fixes the range of a stored value, which is what makes ImageData's scaling to the unit interval well defined.`**

## Examples (34)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (10)

```mathematica
In[1]:= ImageType[Image[{{0.5, 0.25}}]]
Out[1]= "Real"

In[2]:= ImageType[Binarize[Image[{{0.1, 0.9}}]]]
Out[2]= "Bit"

In[3]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[5]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[6]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[7]:= ImageType[chk]
Out[7]= "Real"

In[8]:= ImageType[rgb]
Out[8]= "Real"

In[9]:= ImageType[bit]
Out[9]= "Bit"

In[10]:= ImageType[byte]
Out[10]= "Byte"
```

### Scope (16)

```mathematica
In[11]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[13]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[14]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[15]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[16]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[17]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[18]:= volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[19]:= ImageType[disk]
Out[19]= "Real"

In[20]:= ImageType[ramp]
Out[20]= "Real"

In[21]:= ImageType[zone]
Out[21]= "Real"

In[22]:= ImageType[noise]
Out[22]= "Real"

In[23]:= ImageType[sky]
Out[23]= "Real"

In[24]:= ImageType[vol]
Out[24]= "Real"

In[25]:= ImageType[volb]
Out[25]= "Real"

In[26]:= ImageType[Import[Export["/tmp/mathilda_ex.png", rgb]]]
Out[26]= "Real"
```

### Applications (2)

```mathematica
In[27]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[28]:= Table[ImageType[GaussianFilter[chk, r]], {r, 1, 3}]
Out[28]= {"Real", "Real", "Real"}
```

### Properties & Relations (4)

```mathematica
In[29]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[30]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[31]:= ImageType[chk] === ImageType[GaussianFilter[chk, 1]]
Out[31]= True

In[32]:= ImageType[rgb] === ImageType[ImagePad[rgb, 2]]
Out[32]= True
```

### Neat Examples (2)

```mathematica
In[33]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[34]:= ImageType[zone]
Out[34]= "Real"
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ImageData](../../image-processing/ImageData/)

- Source: [`src/image.c`](https://github.com/stblake/mathilda/blob/main/src/image.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
