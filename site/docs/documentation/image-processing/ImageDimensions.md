# ImageDimensions

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageDimensions[image] gives {width, height}. This is TRANSPOSED relative to ImageData, which returns a height x width array -- the same convention Mathematica uses.`**

## Examples (31)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (9)

```mathematica
In[1]:= ImageDimensions[Image[{{1., 2., 3.}, {4., 5., 6.}}]]
Out[1]= {3, 2}

In[2]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[5]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[6]:= ImageDimensions[chk]
Out[6]= {16, 16}

In[7]:= ImageDimensions[rgb]
Out[7]= {16, 16}

In[8]:= ImageDimensions[bit]
Out[8]= {8, 8}

In[9]:= ImageDimensions[byte]
Out[9]= {16, 16}
```

### Scope (16)

```mathematica
In[10]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[13]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[14]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[15]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[16]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[17]:= volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[18]:= ImageDimensions[disk]
Out[18]= {16, 16}

In[19]:= ImageDimensions[ramp]
Out[19]= {16, 16}

In[20]:= ImageDimensions[zone]
Out[20]= {32, 32}

In[21]:= ImageDimensions[noise]
Out[21]= {32, 32}

In[22]:= ImageDimensions[sky]
Out[22]= {24, 16}

In[23]:= ImageDimensions[vol]
Out[23]= {12, 10, 8}

In[24]:= ImageDimensions[volb]
Out[24]= {12, 10, 8}

In[25]:= ImageDimensions[Import[Export["/tmp/mathilda_ex.png", rgb]]]
Out[25]= {16, 16}
```

### Applications (2)

```mathematica
In[26]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[27]:= Table[ImageDimensions[GaussianFilter[chk, r]], {r, 1, 3}]
Out[27]= {{16, 16}, {16, 16}, {16, 16}}
```

### Properties & Relations (2)

```mathematica
In[28]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[29]:= ImageDimensions[chk] === ImageDimensions[GaussianFilter[chk, 1]]
Out[29]= True
```

### Neat Examples (2)

```mathematica
In[30]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[31]:= ImageDimensions[zone]
Out[31]= {32, 32}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ImageData](../../image-processing/ImageData/)

- Source: [`src/image.c`](https://github.com/stblake/mathilda/blob/main/src/image.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
