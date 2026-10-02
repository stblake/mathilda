# Dilation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Dilation[image, r] gives the maximum over a (2r+1) x (2r+1) square neighbourhood; Dilation[image, elem] uses the SUPPORT of the matrix elem -- its nonzero positions -- as the neighbourhood. This is flat morphology: the element's values do not enter the maximum, which is what keeps Dilation[img, BoxMatrix[1]] and Dilation[img, 1] the same operation. Padding replicates the border, the same rule the convolutions use, which is what makes Dilation >= image hold at the edges too. A full rectangle is separable for the maximum exactly as for a sum, so it costs kw + kh comparisons rather than kw * kh. A "Bit" image gives a "Bit" image; other types give "Real".`**

## Examples (40)

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

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
