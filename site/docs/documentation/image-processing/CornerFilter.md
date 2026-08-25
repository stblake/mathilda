# CornerFilter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CornerFilter[image] gives the corner strength at every pixel, from the eigenvalues of the Gaussian-weighted second-moment matrix of the gradient (the structure tensor). Both eigenvalues small is flat, one large is an edge, both large is a corner. CornerFilter[image, r] sets the window radius (default 2); CornerFilter[image, r, method] selects "MinimumEigenvalue" (the default -- Shi-Tomasi's lambda_min, which is directly "how much does the weaker direction vary" and is comparable across images) or "Harris" (det - 0.04 trace^2, cheaper since it needs no square root, and negative on edges). A STRAIGHT EDGE SCORES ZERO under both: every gradient in the window is parallel, so the matrix has rank 1 and its determinant and smaller eigenvalue vanish. Colour is reduced to luminance first, since a corner is a property of brightness.`**

## Examples (49)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (7)

```mathematica
In[1]:= Max[Flatten[ImageData[CornerFilter[Image[Table[If[j <= 4, 0., 1.], {i, 8}, {j, 8}]]]]]]
Out[1]= 0.0

In[2]:= Options[CornerFilter]
Out[2]= {Method -> "MinimumEigenvalue"}

In[3]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= CornerFilter[chk]
Out[4]= -Image-

In[5]:= CornerFilter[chk, 2]
Out[5]= -Image-

In[6]:= ImageDimensions[CornerFilter[chk, 2]]
Out[6]= {16, 16}

In[7]:= ImageType[CornerFilter[chk, 1]]
Out[7]= "Real"
```

### Scope (22)

```mathematica
In[8]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[9]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[10]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[12]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[13]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[14]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[15]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[16]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[17]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[18]:= CornerFilter[disk, 1]
Out[18]= -Image-

In[19]:= CornerFilter[ramp, 2]
Out[19]= -Image-

In[20]:= CornerFilter[zone, 2]
Out[20]= -Image-

In[21]:= CornerFilter[noise, 3]
Out[21]= -Image-

In[22]:= CornerFilter[rgb, 1]
Out[22]= -Image-

In[23]:= CornerFilter[sky, 2]
Out[23]= -Image-

In[24]:= CornerFilter[bit, 1]
Out[24]= -Image-

In[25]:= CornerFilter[byte, 2]
Out[25]= -Image-

In[26]:= CornerFilter[vol, 1]
Out[26]= -Image-

In[27]:= ImageChannels[CornerFilter[rgb, 2]]
Out[27]= 1

In[28]:= ImageDimensions[CornerFilter[vol, 1]]
Out[28]= {12, 10, 8}

In[29]:= CornerFilter[chk, 4]
Out[29]= -Image-
```

### Options (5)

```mathematica
In[30]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[31]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[32]:= CornerFilter[chk, 2, Method -> "Harris"]
Out[32]= -Image-

In[33]:= CornerFilter[chk, 2, Method -> "MinimumEigenvalue"]
Out[33]= -Image-

In[34]:= CornerFilter[disk, 1, Method -> "Harris"]
Out[34]= -Image-
```

### Applications (6)

```mathematica
In[35]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[36]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[37]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[38]:= Binarize[CornerFilter[noise, 2]]
Out[38]= -Image-

In[39]:= EdgeDetect[CornerFilter[zone, 2]]
Out[39]= -Image-

In[40]:= ImageDimensions[CornerFilter[Import[Export["/tmp/mathilda_ex.png", rgb]], 2]]
Out[40]= {16, 16}
```

### Properties & Relations (6)

```mathematica
In[41]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[42]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[43]:= ImageDimensions[CornerFilter[chk, 3]] === ImageDimensions[chk]
Out[43]= True

In[44]:= Max[Flatten[ImageData[CornerFilter[chk, 2]]]] <= 1.0
Out[44]= True

In[45]:= Min[Flatten[ImageData[CornerFilter[chk, 2]]]] >= 0.0
Out[45]= True

In[46]:= ImageDimensions[CornerFilter[vol, 2]] === ImageDimensions[vol]
Out[46]= True
```

### Neat Examples (3)

```mathematica
In[47]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[48]:= CornerFilter[zone, 4]
Out[48]= -Image-

In[49]:= CornerFilter[zone, 1]
Out[49]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
