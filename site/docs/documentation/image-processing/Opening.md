# Opening

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Opening[image, r] erodes then dilates with the same element, removing bright features smaller than it while leaving larger ones close to their original size. IDEMPOTENT: Opening[Opening[f]] equals Opening[f], which is the defining property and the reason opening twice is not a sharpening loop.`**

## Examples (37)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (6)

```mathematica
In[1]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[3]:= Opening[disk, 1]
Out[3]= -Image-

In[4]:= Opening[disk, 2]
Out[4]= -Image-

In[5]:= ImageDimensions[Opening[disk, 2]]
Out[5]= {16, 16}

In[6]:= Opening[bit, 1]
Out[6]= -Image-
```

### Scope (19)

```mathematica
In[7]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[8]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[9]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[10]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[11]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[13]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[14]:= volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[15]:= Opening[chk, 1]
Out[15]= -Image-

In[16]:= Opening[ramp, 1]
Out[16]= -Image-

In[17]:= Opening[noise, 2]
Out[17]= -Image-

In[18]:= Opening[rgb, 1]
Out[18]= -Image-

In[19]:= Opening[byte, 1]
Out[19]= -Image-

In[20]:= Opening[vol, 1]
Out[20]= -Image-

In[21]:= Opening[volb, 1]
Out[21]= -Image-

In[22]:= Opening[disk, 3]
Out[22]= -Image-

In[23]:= Opening[disk, 4]
Out[23]= -Image-

In[24]:= ImageChannels[Opening[rgb, 1]]
Out[24]= 3

In[25]:= ImageDimensions[Opening[vol, 2]]
Out[25]= {12, 10, 8}
```

### Applications (4)

```mathematica
In[26]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[27]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[28]:= Binarize[Opening[noise, 1]]
Out[28]= -Image-

In[29]:= EdgeDetect[Opening[disk, 1]]
Out[29]= -Image-
```

### Properties & Relations (6)

```mathematica
In[30]:= Module[{img = Image[{{0., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}}]}, ImageData[Opening[Opening[img, 1], 1]] === ImageData[Opening[img, 1]]]
Out[30]= True

In[31]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[32]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[33]:= ImageData[Opening[disk, 0]] === ImageData[disk]
Out[33]= True

In[34]:= ImageDimensions[Opening[disk, 3]] === ImageDimensions[disk]
Out[34]= True

In[35]:= Max[Flatten[ImageData[Opening[bit, 1]]]] <= 1.0
Out[35]= True
```

### Neat Examples (2)

```mathematica
In[36]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[37]:= Opening[zone, 2]
Out[37]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
