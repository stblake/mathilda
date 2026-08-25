# ImageReflect

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageReflect[image] reflects top to bottom; ImageReflect[image, Left] or Right reflects left to right, and Top or Bottom is the vertical reflection again -- either name of a pair selects the same axis, since reflecting to the top and reflecting to the bottom are one operation. For an Image3D, Front or Back selects the DEPTH axis, the pair Mathematica uses for volumes; those two DECLINE on a plane, which has no depth axis, rather than being reinterpreted as some other axis and turning a mistake into a wrong picture. A reflection is a pure index permutation, so it interpolates nothing: reflecting twice about the same axis is the identity bit for bit, and reflections about different axes commute exactly.`**

## Examples (38)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (8)

```mathematica
In[1]:= ImageData[ImageReflect[Image[{{1., 2.}, {3., 4.}}]]]
Out[1]= {{3.0, 4.0}, {1.0, 2.0}}

In[2]:= ImageData[ImageReflect[Image[{{1., 2.}, {3., 4.}}], Left]]
Out[2]= {{2.0, 1.0}, {4.0, 3.0}}

In[3]:= Module[{img = Image[{{1., 2.}, {3., 4.}}]}, ImageReflect[ImageReflect[img]] === img]
Out[3]= True

In[4]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[5]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[6]:= ImageReflect[chk]
Out[6]= -Image-

In[7]:= ImageReflect[disk]
Out[7]= -Image-

In[8]:= ImageDimensions[ImageReflect[chk]]
Out[8]= {16, 16}
```

### Scope (20)

```mathematica
In[9]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[10]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[11]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[12]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[13]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[14]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[15]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[16]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[17]:= volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[18]:= ImageReflect[ramp]
Out[18]= -Image-

In[19]:= ImageReflect[zone]
Out[19]= -Image-

In[20]:= ImageReflect[noise]
Out[20]= -Image-

In[21]:= ImageReflect[rgb]
Out[21]= -Image-

In[22]:= ImageReflect[sky]
Out[22]= -Image-

In[23]:= ImageReflect[bit]
Out[23]= -Image-

In[24]:= ImageReflect[byte]
Out[24]= -Image-

In[25]:= ImageReflect[vol]
Out[25]= -Image-

In[26]:= ImageReflect[volb]
Out[26]= -Image-

In[27]:= ImageChannels[ImageReflect[rgb]]
Out[27]= 3

In[28]:= ImageDimensions[ImageReflect[vol]]
Out[28]= {12, 10, 8}
```

### Applications (4)

```mathematica
In[29]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[30]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[31]:= Binarize[ImageReflect[zone]]
Out[31]= -Image-

In[32]:= Dilation[ImageReflect[disk], 1]
Out[32]= -Image-
```

### Properties & Relations (4)

```mathematica
In[33]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[34]:= ImageDimensions[ImageReflect[chk]] === ImageDimensions[chk]
Out[34]= True

In[35]:= Max[Flatten[ImageData[ImageReflect[chk]]]] <= 1.0
Out[35]= True

In[36]:= Min[Flatten[ImageData[ImageReflect[chk]]]] >= 0.0
Out[36]= True
```

### Neat Examples (2)

```mathematica
In[37]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[38]:= ImageReflect[zone]
Out[38]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
