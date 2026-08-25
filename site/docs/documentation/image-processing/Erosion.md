# Erosion

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Erosion[image, r] gives the minimum over a (2r+1) x (2r+1) square neighbourhood; Erosion[image, elem] uses the support of elem. Dual to Dilation: for a symmetric element, Erosion[f, k] equals 1 - Dilation[1 - f, k] exactly, which holds at the border only because the replicate padding is itself self-dual.`**

## Examples (38)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (7)

```mathematica
In[1]:= ImageData[Erosion[Image[{{1., 1., 1.}, {1., 0., 1.}, {1., 1., 1.}}], 1]]
Out[1]= {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}

In[2]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[4]:= Erosion[disk, 1]
Out[4]= -Image-

In[5]:= Erosion[disk, 2]
Out[5]= -Image-

In[6]:= ImageDimensions[Erosion[disk, 2]]
Out[6]= {16, 16}

In[7]:= Erosion[bit, 1]
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

In[16]:= Erosion[chk, 1]
Out[16]= -Image-

In[17]:= Erosion[ramp, 1]
Out[17]= -Image-

In[18]:= Erosion[noise, 2]
Out[18]= -Image-

In[19]:= Erosion[rgb, 1]
Out[19]= -Image-

In[20]:= Erosion[byte, 1]
Out[20]= -Image-

In[21]:= Erosion[vol, 1]
Out[21]= -Image-

In[22]:= Erosion[volb, 1]
Out[22]= -Image-

In[23]:= Erosion[disk, 3]
Out[23]= -Image-

In[24]:= Erosion[disk, 4]
Out[24]= -Image-

In[25]:= ImageChannels[Erosion[rgb, 1]]
Out[25]= 3

In[26]:= ImageDimensions[Erosion[vol, 2]]
Out[26]= {12, 10, 8}
```

### Applications (4)

```mathematica
In[27]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[28]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[29]:= Binarize[Erosion[noise, 1]]
Out[29]= -Image-

In[30]:= EdgeDetect[Erosion[disk, 1]]
Out[30]= -Image-
```

### Properties & Relations (6)

```mathematica
In[31]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[32]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[33]:= ImageData[Erosion[disk, 0]] === ImageData[disk]
Out[33]= True

In[34]:= ImageDimensions[Erosion[disk, 3]] === ImageDimensions[disk]
Out[34]= True

In[35]:= ImageData[Erosion[Erosion[disk, 1], 1]] === ImageData[Erosion[disk, 2]]
Out[35]= True

In[36]:= Max[Flatten[ImageData[Erosion[bit, 1]]]] <= 1.0
Out[36]= True
```

### Neat Examples (2)

```mathematica
In[37]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[38]:= Erosion[zone, 2]
Out[38]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Dilation](../../image-processing/Dilation/), [Image3D](../../image-processing/Image3D/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
