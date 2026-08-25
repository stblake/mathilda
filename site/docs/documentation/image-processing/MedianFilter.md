# MedianFilter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MedianFilter[image, r] replaces each pixel with the median over a (2r+1) x (2r+1) neighbourhood. Unlike a Gaussian it removes an isolated outlier EXACTLY rather than attenuating and smearing it, which is what makes it the filter for salt-and-pepper noise. It is also the one filter here that is NOT separable: a sum, a maximum and a minimum all decompose because they ignore grouping, but a median depends on a value's rank within the whole window, and grouping destroys rank -- the median of row medians of {{1,2,9},{3,4,5},{6,7,8}} is 4 where the true median is 5. For an even window the lower middle is taken rather than the average of the two, so the output is always one of the inputs.`**

## Examples (46)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (5)

```mathematica
In[1]:= ImageData[MedianFilter[Image[{{0.25, 0.25, 0.25}, {0.25, 1., 0.25}, {0.25, 0.25, 0.25}}], 1]]
Out[1]= {{0.25, 0.25, 0.25}, {0.25, 0.25, 0.25}, {0.25, 0.25, 0.25}}

In[2]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= MedianFilter[chk, 2]
Out[3]= -Image-

In[4]:= ImageDimensions[MedianFilter[chk, 2]]
Out[4]= {16, 16}

In[5]:= ImageType[MedianFilter[chk, 1]]
Out[5]= "Real"
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

In[16]:= MedianFilter[disk, 1]
Out[16]= -Image-

In[17]:= MedianFilter[ramp, 2]
Out[17]= -Image-

In[18]:= MedianFilter[zone, 2]
Out[18]= -Image-

In[19]:= MedianFilter[noise, 3]
Out[19]= -Image-

In[20]:= MedianFilter[rgb, 1]
Out[20]= -Image-

In[21]:= MedianFilter[sky, 2]
Out[21]= -Image-

In[22]:= MedianFilter[bit, 1]
Out[22]= -Image-

In[23]:= MedianFilter[byte, 2]
Out[23]= -Image-

In[24]:= MedianFilter[vol, 1]
Out[24]= -Image-

In[25]:= ImageChannels[MedianFilter[rgb, 2]]
Out[25]= 3

In[26]:= ImageDimensions[MedianFilter[vol, 1]]
Out[26]= {12, 10, 8}

In[27]:= MedianFilter[chk, 4]
Out[27]= -Image-
```

### Applications (6)

```mathematica
In[28]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[29]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[30]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[31]:= Binarize[MedianFilter[noise, 2]]
Out[31]= -Image-

In[32]:= EdgeDetect[MedianFilter[zone, 2]]
Out[32]= -Image-

In[33]:= ImageDimensions[MedianFilter[Import[Export["/tmp/mathilda_ex.png", rgb]], 2]]
Out[33]= {16, 16}
```

### Properties & Relations (10)

```mathematica
In[34]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[35]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[36]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[37]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[38]:= ImageDimensions[MedianFilter[chk, 3]] === ImageDimensions[chk]
Out[38]= True

In[39]:= ImageChannels[MedianFilter[rgb, 2]] === ImageChannels[rgb]
Out[39]= True

In[40]:= ImageData[MedianFilter[ramp, 0]] === ImageData[ramp]
Out[40]= True

In[41]:= Max[Flatten[ImageData[MedianFilter[chk, 2]]]] <= 1.0
Out[41]= True

In[42]:= Min[Flatten[ImageData[MedianFilter[chk, 2]]]] >= 0.0
Out[42]= True

In[43]:= ImageDimensions[MedianFilter[vol, 2]] === ImageDimensions[vol]
Out[43]= True
```

### Neat Examples (3)

```mathematica
In[44]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[45]:= MedianFilter[zone, 4]
Out[45]= -Image-

In[46]:= MedianFilter[zone, 1]
Out[46]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
