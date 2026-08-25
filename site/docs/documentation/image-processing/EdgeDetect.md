# EdgeDetect

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EdgeDetect[image] finds edges by the Canny algorithm, giving a "Bit" image. EdgeDetect[image, r] sets the Gaussian smoothing radius (default 2; 0 means no smoothing). EdgeDetect[image, r, t] sets the high threshold explicitly. Four stages: smooth, because a derivative amplifies noise; gradient by the normalised Sobel pair; non-maximum suppression along the gradient direction, which is what makes an edge ONE pixel wide rather than a thick band; and hysteresis, keeping any pixel above the high threshold plus any above 0.4 of it that is 8-connected to one, so a real edge survives its faint stretches while isolated weak responses do not. The high threshold defaults to Otsu's method applied to the SUPPRESSED magnitude, where the two classes really are edge against non-edge; on the raw magnitude it would be dominated by the ridge flanks.`**

## Examples (31)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (5)

```mathematica
In[1]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= EdgeDetect[chk]
Out[3]= -Image-

In[4]:= EdgeDetect[disk]
Out[4]= -Image-

In[5]:= ImageDimensions[EdgeDetect[chk]]
Out[5]= {16, 16}
```

### Scope (16)

```mathematica
In[6]:= Map[Total, ImageData[EdgeDetect[step, 0], "Bit"]]
Out[6]= ImageData[step, "Bit"]

In[7]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[8]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[9]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[10]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[12]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[13]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[14]:= EdgeDetect[ramp]
Out[14]= -Image-

In[15]:= EdgeDetect[zone]
Out[15]= -Image-

In[16]:= EdgeDetect[noise]
Out[16]= -Image-

In[17]:= EdgeDetect[rgb]
Out[17]= -Image-

In[18]:= EdgeDetect[sky]
Out[18]= -Image-

In[19]:= EdgeDetect[bit]
Out[19]= -Image-

In[20]:= EdgeDetect[byte]
Out[20]= -Image-

In[21]:= ImageChannels[EdgeDetect[rgb]]
Out[21]= 1
```

### Applications (4)

```mathematica
In[22]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[23]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[24]:= Binarize[EdgeDetect[zone]]
Out[24]= -Image-

In[25]:= Dilation[EdgeDetect[disk], 1]
Out[25]= -Image-
```

### Properties & Relations (4)

```mathematica
In[26]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[27]:= ImageDimensions[EdgeDetect[chk]] === ImageDimensions[chk]
Out[27]= True

In[28]:= Max[Flatten[ImageData[EdgeDetect[chk]]]] <= 1.0
Out[28]= True

In[29]:= Min[Flatten[ImageData[EdgeDetect[chk]]]] >= 0.0
Out[29]= True
```

### Neat Examples (2)

```mathematica
In[30]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[31]:= EdgeDetect[zone]
Out[31]= -Image-
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [DerivativeFilter](../../image-processing/DerivativeFilter/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
