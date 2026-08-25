# ImageCrop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageCrop[image, {w, h}] crops to w x h about the centre, any odd remainder going to the right and bottom -- the same floor-division convention the kernel centres use, which is what makes ImageCrop[ImagePad[image, m], ImageDimensions[image]] exactly the original image. A crop may not enlarge. ImageCrop[image] instead TRIMS A UNIFORM BORDER, asking how much of the frame carries no information; the border colour is read from a corner rather than assumed black, since a scanned page's margin is white. An entirely uniform image comes back unchanged, there being no content to keep and a zero-sized image not being one.`**

## Examples (34)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (7)

```mathematica
In[1]:= ImageData[ImageCrop[Image[{{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}}], {1, 1}]]
Out[1]= {{5.0}}

In[2]:= ImageDimensions[ImageCrop[Image[{{0., 0., 0.}, {0., 0.5, 0.}, {0., 0., 0.}}]]]
Out[2]= {1, 1}

In[3]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[5]:= ImageCrop[chk, {8, 8}]
Out[5]= -Image-

In[6]:= ImageDimensions[ImageCrop[chk, {8, 8}]]
Out[6]= {8, 8}

In[7]:= ImageCrop[disk, {12, 12}]
Out[7]= -Image-
```

### Scope (16)

```mathematica
In[8]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[9]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[10]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[11]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[13]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[14]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[15]:= ImageCrop[rgb, {8, 8}]
Out[15]= -Image-

In[16]:= ImageCrop[sky, {12, 8}]
Out[16]= -Image-

In[17]:= ImageCrop[bit, {4, 4}]
Out[17]= -Image-

In[18]:= ImageCrop[byte, {8, 8}]
Out[18]= -Image-

In[19]:= ImageCrop[zone, {16, 16}]
Out[19]= -Image-

In[20]:= ImageCrop[noise, {16, 24}]
Out[20]= -Image-

In[21]:= ImageCrop[ramp, {8, 16}]
Out[21]= -Image-

In[22]:= ImageChannels[ImageCrop[rgb, {8, 8}]]
Out[22]= 3

In[23]:= ImageDimensions[ImageCrop[zone, {20, 10}]]
Out[23]= {20, 10}
```

### Applications (4)

```mathematica
In[24]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[25]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[26]:= Binarize[ImageCrop[zone, {16, 16}]]
Out[26]= -Image-

In[27]:= EdgeDetect[ImageCrop[disk, {12, 12}]]
Out[27]= -Image-
```

### Properties & Relations (5)

```mathematica
In[28]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[29]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[30]:= ImageDimensions[ImageCrop[chk, {8, 8}]] === {8, 8}
Out[30]= True

In[31]:= ImageChannels[ImageCrop[rgb, {4, 4}]] === ImageChannels[rgb]
Out[31]= True

In[32]:= ImageDimensions[ImageCrop[chk, {16, 16}]] === ImageDimensions[chk]
Out[32]= True
```

### Neat Examples (2)

```mathematica
In[33]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[34]:= ImageCrop[zone, {24, 8}]
Out[34]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
