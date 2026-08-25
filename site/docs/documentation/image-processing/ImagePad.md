# ImagePad

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImagePad[image, m] pads m pixels on every side; ImagePad[image, {{left, right}, {bottom, top}}] pads each side separately, in Mathematica's VISUAL order -- so `top` adds rows at the start of the data, since row 1 is the top of the image. Negative amounts crop, but may not erase the image. ImagePad[image, m, v] fills with the value v (default 0); ImagePad[image, m, "Fixed"] replicates the edge pixel, the same boundary rule the filters use, so padding then filtering composes with it; ImagePad[image, m, "Reflected"] mirrors WITHOUT repeating the edge -- {1,2,3} padded by 1 gives {2,1,2,3,2}, not {1,1,2,3,3}, because doubling the edge sample biases any later average toward the border. Reflection uses a period of 2n-2, so padding deeper than the image still works.`**

## Examples (29)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (9)

```mathematica
In[1]:= ImageDimensions[ImagePad[Image[{{1., 2.}, {3., 4.}}], 1]]
Out[1]= {4, 4}

In[2]:= ImageData[ImagePad[Image[{{1., 2.}, {3., 4.}}], 1]]
Out[2]= {{0.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 2.0, 0.0}, {0.0, 3.0, 4.0, 0.0}, {0.0, 0.0, 0.0, 0.0}}

In[3]:= ImageData[ImagePad[Image[{{1., 2., 3.}}], {{1, 1}, {0, 0}}, "Reflected"]]
Out[3]= {{2.0, 1.0, 2.0, 3.0, 2.0}}

In[4]:= Module[{img = Image[{{1., 2.}, {3., 4.}}]}, ImageCrop[ImagePad[img, 2], ImageDimensions[img]] === img]
Out[4]= True

In[5]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[6]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[7]:= ImagePad[chk, 2]
Out[7]= -Image-

In[8]:= ImageDimensions[ImagePad[chk, 2]]
Out[8]= {20, 20}

In[9]:= ImagePad[disk, 1]
Out[9]= -Image-
```

### Scope (12)

```mathematica
In[10]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[13]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[14]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[15]:= ImagePad[rgb, 2]
Out[15]= -Image-

In[16]:= ImagePad[bit, 1]
Out[16]= -Image-

In[17]:= ImagePad[byte, 2]
Out[17]= -Image-

In[18]:= ImagePad[vol, 1]
Out[18]= -Image-

In[19]:= ImagePad[chk, {{1, 2}, {3, 4}}]
Out[19]= -Image-

In[20]:= ImageDimensions[ImagePad[vol, 1]]
Out[20]= {14, 12, 10}

In[21]:= ImageChannels[ImagePad[rgb, 2]]
Out[21]= 3
```

### Applications (2)

```mathematica
In[22]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[23]:= EdgeDetect[ImagePad[disk, 2]]
Out[23]= -Image-
```

### Properties & Relations (4)

```mathematica
In[24]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[25]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[26]:= ImageChannels[ImagePad[rgb, 2]] === ImageChannels[rgb]
Out[26]= True

In[27]:= ImageDimensions[ImagePad[chk, 0]] === ImageDimensions[chk]
Out[27]= True
```

### Neat Examples (2)

```mathematica
In[28]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[29]:= ImagePad[zone, 4]
Out[29]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
