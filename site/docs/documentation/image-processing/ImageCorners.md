# ImageCorners

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageCorners[image] gives the positions of corners. ImageCorners[image, r, t, d, n] sets the window radius (default 2), the threshold as a fraction of the largest response (0.05), the MINIMUM SEPARATION in pixels (0), and the maximum number of features (all). Three filters apply in that order because each removes what the others cannot: a threshold alone returns a blob of adjacent pixels per corner since the response is smooth; 3x3 non-maximum suppression alone returns a maximum in every flat region since a plateau of zeros has maxima; and separation is what makes the list usable, since the first two leave clusters a pixel apart -- 4104 of them on a noise-like 512x512 image. Separation is greedy in DESCENDING RESPONSE order, so the survivor of a cluster is its strongest member rather than whichever came first in raster order, and the feature limit is applied last: before separation it would return n positions from a single cluster. The result is sorted strongest first, ties broken by position so the same image always gives the same list. Positions are {row, column}, 1-based, so each indexes ImageData directly; that is NOT Mathematica's {x, y} from the bottom left, and Mathematica spells the feature limit as a MaxFeatures option where this takes it positionally -- both differences are stated rather than guessed.`**

## Examples (39)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (7)

```mathematica
In[1]:= ImageCorners[Image[Table[If[i >= 5 && j >= 5, 1., 0.], {i, 12}, {j, 12}]]]
Out[1]= {{5, 5}}

In[2]:= Options[ImageCorners]
Out[2]= {MaxFeatures -> Infinity}

In[3]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[5]:= ImageCorners[chk]
Out[5]= {{15, 2}, {15, 15}, {2, 2}, {2, 15}}

In[6]:= ImageCorners[disk]
Out[6]= {{5, 12}, {12, 5}, {12, 12}, {5, 5}, {4, 8}, {8, 4}, {8, 13}, {13, 8}}

In[7]:= ImageDimensions[ImageCorners[chk]]
Out[7]= ImageDimensions[{{15, 2}, {15, 15}, {2, 2}, {2, 15}}]
```

### Scope (20)

```mathematica
In[8]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[9]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[10]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[11]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[13]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[14]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[15]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[16]:= volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[17]:= ImageCorners[ramp]
Out[17]= {}

In[18]:= ImageCorners[zone]
Out[18]= {{31, 31}, {11, 10}, {11, 22}, {21, 10}, {21, 22}, {6, 7}, {6, 25}, {25, 6}, {26, 25}, {22, 28}, {28, 10}, {10, 28}, {28, 22}, {2, 16}, {16, 2}, {29, 12}, {29, 20}, {12, 29}, {20, 29}}

In[19]:= ImageCorners[noise]
Out[19]= {{13, 30}, {10, 3}, {18, 31}, {5, 2}, {4, 8}, {4, 14}, {8, 29}, {9, 9}, {9, 15}, {14, 10}, {14, 16}, {19, 11}, {19, 17}, {24, 12}, {24, 18}, {8, 23}, {13, 24}, {14, 4}, {18, 25}, {19, 5}, {23, 26}, {24, 6}, {28, 27}, {4, 18}, {9, 19}, {14, 20}, {19, 21}, {24, 22}, {29, 23}}

In[20]:= ImageCorners[rgb]
Out[20]= {}

In[21]:= ImageCorners[sky]
Out[21]= {}

In[22]:= ImageCorners[bit]
Out[22]= {}

In[23]:= ImageCorners[byte]
Out[23]= {{12, 15}, {14, 11}}

In[24]:= ImageCorners[vol]
Out[24]= {}

In[25]:= ImageCorners[volb]
Out[25]= {}

In[26]:= ImageChannels[ImageCorners[rgb]]
Out[26]= ImageChannels[{}]

In[27]:= ImageDimensions[ImageCorners[vol]]
Out[27]= ImageDimensions[{}]
```

### Options (6)

```mathematica
In[28]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[29]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[30]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[31]:= ImageCorners[chk, MaxFeatures -> 3]
Out[31]= {{15, 2}, {15, 15}, {2, 2}}

In[32]:= ImageCorners[disk, MaxFeatures -> 2]
Out[32]= {{5, 12}, {12, 5}}

In[33]:= ImageCorners[zone, MaxFeatures -> 5]
Out[33]= {{31, 31}, {11, 10}, {11, 22}, {21, 10}, {21, 22}}
```

### Applications (4)

```mathematica
In[34]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[35]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[36]:= Binarize[ImageCorners[zone]]
Out[36]= Binarize[{{31, 31}, {11, 10}, {11, 22}, {21, 10}, {21, 22}, {6, 7}, {6, 25}, {25, 6}, {26, 25}, {22, 28}, {28, 10}, {10, 28}, {28, 22}, {2, 16}, {16, 2}, {29, 12}, {29, 20}, {12, 29}, {20, 29}}]

In[37]:= Dilation[ImageCorners[disk], 1]
Out[37]= Dilation[{{5, 12}, {12, 5}, {12, 12}, {5, 5}, {4, 8}, {8, 4}, {8, 13}, {13, 8}}, 1]
```

### Neat Examples (2)

```mathematica
In[38]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[39]:= ImageCorners[zone]
Out[39]= {{31, 31}, {11, 10}, {11, 22}, {21, 10}, {21, 22}, {6, 7}, {6, 25}, {25, 6}, {26, 25}, {22, 28}, {28, 10}, {10, 28}, {28, 22}, {2, 16}, {16, 2}, {29, 12}, {29, 20}, {12, 29}, {20, 29}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [ImageData](../../image-processing/ImageData/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
