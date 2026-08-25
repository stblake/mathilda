# Dilation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Dilation[image, r] gives the maximum over a (2r+1) x (2r+1) square neighbourhood; Dilation[image, elem] uses the SUPPORT of the matrix elem -- its nonzero positions -- as the neighbourhood. This is flat morphology: the element's values do not enter the maximum, which is what keeps Dilation[img, BoxMatrix[1]] and Dilation[img, 1] the same operation. Padding replicates the border, the same rule the convolutions use, which is what makes Dilation >= image hold at the edges too. A full rectangle is separable for the maximum exactly as for a sum, so it costs kw + kh comparisons rather than kw * kh.`**

## Examples (38)

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

In[16]:= Dilation[chk, 1]
Out[16]= -Image-

In[17]:= Dilation[ramp, 1]
Out[17]= -Image-

In[18]:= Dilation[noise, 2]
Out[18]= -Image-

In[19]:= Dilation[rgb, 1]
Out[19]= -Image-

In[20]:= Dilation[byte, 1]
Out[20]= -Image-

In[21]:= Dilation[vol, 1]
Out[21]= -Image-

In[22]:= Dilation[volb, 1]
Out[22]= -Image-

In[23]:= Dilation[disk, 3]
Out[23]= -Image-

In[24]:= Dilation[disk, 4]
Out[24]= -Image-

In[25]:= ImageChannels[Dilation[rgb, 1]]
Out[25]= 3

In[26]:= ImageDimensions[Dilation[vol, 2]]
Out[26]= {12, 10, 8}
```

### Applications (4)

```mathematica
In[27]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[28]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[29]:= Binarize[Dilation[noise, 1]]
Out[29]= -Image-

In[30]:= EdgeDetect[Dilation[disk, 1]]
Out[30]= -Image-
```

### Properties & Relations (6)

```mathematica
In[31]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[32]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[33]:= ImageData[Dilation[disk, 0]] === ImageData[disk]
Out[33]= True

In[34]:= ImageDimensions[Dilation[disk, 3]] === ImageDimensions[disk]
Out[34]= True

In[35]:= ImageData[Dilation[Dilation[disk, 1], 1]] === ImageData[Dilation[disk, 2]]
Out[35]= True

In[36]:= Max[Flatten[ImageData[Dilation[bit, 1]]]] <= 1.0
Out[36]= True
```

### Neat Examples (2)

```mathematica
In[37]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[38]:= Dilation[zone, 2]
Out[38]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
