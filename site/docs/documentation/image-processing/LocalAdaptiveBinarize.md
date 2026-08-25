# LocalAdaptiveBinarize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LocalAdaptiveBinarize[image, r] binarizes by comparing each pixel to the MEAN of its own (2r+1)x(2r+1) neighbourhood, and LocalAdaptiveBinarize[image, r, {c1, c2, c3}] to c1*mean + c2*stddev + c3. A global threshold cannot binarize unevenly lit content, and that is not a tuning problem: if one half of a page is darker than the other, no single number separates ink from paper in both halves at once. Mean alone (the default {1, 0, 0}) is Bradley's method; a negative c2 is Sauvola's, tightening the threshold where the neighbourhood is busy. Summed-area tables make the window statistics O(1) per pixel regardless of r -- without them a radius-16 window would be 1089 taps per pixel. The result is typed "Bit", since it is binary by construction. Colour is reduced to luminance first.`**

## Examples (42)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (5)

```mathematica
In[1]:= ImageData[LocalAdaptiveBinarize[Image[{{0.2, 0.3, 0.9}, {0.2, 0.8, 0.9}, {0.1, 0.2, 0.3}}], 1]]
Out[1]= {{0.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {0.0, 0.0, 0.0}}

In[2]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= LocalAdaptiveBinarize[chk, 2]
Out[3]= -Image-

In[4]:= ImageDimensions[LocalAdaptiveBinarize[chk, 2]]
Out[4]= {16, 16}

In[5]:= ImageType[LocalAdaptiveBinarize[chk, 1]]
Out[5]= "Bit"
```

### Scope (22)

```mathematica
In[6]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[7]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[8]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[9]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[10]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[11]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[13]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[14]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[15]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[16]:= LocalAdaptiveBinarize[disk, 1]
Out[16]= -Image-

In[17]:= LocalAdaptiveBinarize[ramp, 2]
Out[17]= -Image-

In[18]:= LocalAdaptiveBinarize[zone, 2]
Out[18]= -Image-

In[19]:= LocalAdaptiveBinarize[noise, 3]
Out[19]= -Image-

In[20]:= LocalAdaptiveBinarize[rgb, 1]
Out[20]= -Image-

In[21]:= LocalAdaptiveBinarize[sky, 2]
Out[21]= -Image-

In[22]:= LocalAdaptiveBinarize[bit, 1]
Out[22]= -Image-

In[23]:= LocalAdaptiveBinarize[byte, 2]
Out[23]= -Image-

In[24]:= LocalAdaptiveBinarize[vol, 1]
Out[24]= -Image-

In[25]:= ImageChannels[LocalAdaptiveBinarize[rgb, 2]]
Out[25]= 1

In[26]:= ImageDimensions[LocalAdaptiveBinarize[vol, 1]]
Out[26]= {12, 10, 8}

In[27]:= LocalAdaptiveBinarize[chk, 4]
Out[27]= -Image-
```

### Applications (6)

```mathematica
In[28]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[29]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[30]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[31]:= Binarize[LocalAdaptiveBinarize[noise, 2]]
Out[31]= -Image-

In[32]:= EdgeDetect[LocalAdaptiveBinarize[zone, 2]]
Out[32]= -Image-

In[33]:= ImageDimensions[LocalAdaptiveBinarize[Import[Export["/tmp/mathilda_ex.png", rgb]], 2]]
Out[33]= {16, 16}
```

### Properties & Relations (6)

```mathematica
In[34]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[35]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[36]:= ImageDimensions[LocalAdaptiveBinarize[chk, 3]] === ImageDimensions[chk]
Out[36]= True

In[37]:= Max[Flatten[ImageData[LocalAdaptiveBinarize[chk, 2]]]] <= 1.0
Out[37]= True

In[38]:= Min[Flatten[ImageData[LocalAdaptiveBinarize[chk, 2]]]] >= 0.0
Out[38]= True

In[39]:= ImageDimensions[LocalAdaptiveBinarize[vol, 2]] === ImageDimensions[vol]
Out[39]= True
```

### Neat Examples (3)

```mathematica
In[40]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[41]:= LocalAdaptiveBinarize[zone, 4]
Out[41]= -Image-

In[42]:= LocalAdaptiveBinarize[zone, 1]
Out[42]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [Binarize](../../image-processing/Binarize/), [ImagePad](../../image-processing/ImagePad/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
