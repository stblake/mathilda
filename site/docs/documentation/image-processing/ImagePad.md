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
```

```mathematica
In[7]:= ImagePad[chk, 2]
Out[7]= -Image-
```

![20x20 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFAAAABQCAYAAACOEfKtAAAArklEQVR42u3ZsQ0DIRQFwf0n998y7uACE4FnczsYgYTeTbXSzz0IAAIECFAAAZ7ZZ/cP1np/Rs7M1b93Al1hgAABCiDAW9+B//7OcwJdYYAABRAgwEub3e/C9kC5wgABAhRAgNkDswcKIECAAAUQYPZAe6BcYYAAAQogwOyB2QMFECBAgAIIMHtg9kABBAgQoAACzB6YPdAVFkCAAAEKYKftgU6gAAIECFAAAR7ZFxJEZaE9y/SYAAAAAElFTkSuQmCC)

```mathematica
In[8]:= ImageDimensions[ImagePad[chk, 2]]
Out[8]= {20, 20}
```

```mathematica
In[9]:= ImagePad[disk, 1]
Out[9]= -Image-
```

![18x18 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFoAAABaCAYAAAA4qEECAAAAxklEQVR42u3bUQrEMAhAQS29/5XNEfpjEyHzTrAMYgjNZkRU6PceBKBBCzRo0AINWqBBgxZo0AINGrRAgxbomb2TfkxV7+fLzDTRVodAg9b3ebHjAU33IdeOsOHQNNFWB2iBBg0aAWjQAn3PzXD6LfDUbdFEWx2gBRo0aASgQQs0aNACDVqgQYMWaNACHV6Tek0qqwO0QIMGLdCgBfram6E/3Zto0KAF2mEoEw0atECDFmjQoAUaNGgEoEELNGjQAg1aoKe2AEnRFatwryKVAAAAAElFTkSuQmCC)

### Scope (12)

```mathematica
In[10]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[13]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[14]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];
```

```mathematica
In[15]:= ImagePad[rgb, 2]
Out[15]= -Image-
```

![20x20 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFAAAABQCAYAAACOEfKtAAAAt0lEQVR42u3cMQqAMAADwFTc2yf3yf2B7i6Cg1J72RRBONLgZElyRB5nQwAQIECAAhDgnNnvHmitX27k3euaT98/StdARxjgyhu4ONDQQEcYoA0UDQQI0AaKBgIEaANFAwECtIGigQAB2kDRQIAAbaAGCkCAv93AykgDAQL0HSgaCBCgDRQNBAjQBooGAgRoA0UDAQK0gRqIACBAG6iBAhCgDdRAAQhwtpT4j7QGAgQIUAACBLhiTmjTB5x+7iMrAAAAAElFTkSuQmCC)

```mathematica
In[16]:= ImagePad[bit, 1]
Out[16]= -Image-
```

![10x10 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFoAAABaCAYAAAA4qEECAAAAwElEQVR42u3ZQQqEQBAEwWzx/18en6CCejHyvOwhGBqkplrp9TYEoEELNGjQAg1at9uv/Git82+amfnt/3jRTgdogQYt0KBB653mqYXFR40X7XSAFmjQAg0atLKwWFgEGjRogQYt0KCzsGRhycLidAg0aNACDVqgs7BkYZHTAVqgQYMWaNCysGRhcToEGrRAgwYt0FlYsrB40U4HaIEGLdCgQSsLSxYWgQYNWqBBCzRoHyzyokELNGjQAg1aoL/sAGLLZrERTEpoAAAAAElFTkSuQmCC)

```mathematica
In[17]:= ImagePad[byte, 2]
Out[17]= -Image-
```

![20x20 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFAAAABQCAYAAACOEfKtAAAA10lEQVR42u3cMQqFMBRE0acIqdJlx1lsNpBea6uAjZqc6R4I4mXmwm/+FhFnyOPsEAAIIIAACoAA/jPH6IFSyu1OKS1111o10IQBXNiBOeelHaiBJgwgB3KgBpowgBzIgRpowgByIAdqoAkDyIEcqIEmDCAHcqAGmjCAHMiBYsIAciAHigkDyIEcKCYMIAdyoAAIIAdyoAAIIAd+/R59X+9dA00YQA6c1nGju7WmgSYMIAdO6zi/hU0YQA78sgPffr8GmjCAU2cL/yOtgQACCKAACCCAK+YC8NAgVJZGgxkAAAAASUVORK5CYII=)

```mathematica
In[18]:= ImagePad[vol, 1]
Out[18]= -Image-
```

![14x12x10 volume](data:image/png;base64,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)

```mathematica
In[19]:= ImagePad[chk, {{1, 2}, {3, 4}}]
Out[19]= -Image-
```

![19x23 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEwAAABcCAYAAADEQVOwAAAAu0lEQVR42u3aQQqAMBAEwVnx/19efyDmkEBM9d1L4YIMVpKOPnchAAYMGDBgAgYMGDBgAgYMGLCdu0cf6H6fz6rq1897w5wkMGDAgAnYrO+w07+zvGFOEhgwYAIGbFE1+n+YPUxOEhgwYMAELPaw2MOcJDBgAgYMWOxh9jABAwYMGDBgij3MHuYkgQETMGDAYg+LPcxJChgwYMCAKfYwe5iTBAYMmIABy657mDdMwIABAwZMwIABAwbs+B4DOGW5eOB8AgAAAABJRU5ErkJggg==)

```mathematica
In[20]:= ImageDimensions[ImagePad[vol, 1]]
Out[20]= {14, 12, 10}

In[21]:= ImageChannels[ImagePad[rgb, 2]]
Out[21]= 3
```

### Applications (2)

```mathematica
In[22]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[23]:= EdgeDetect[ImagePad[disk, 2]]
Out[23]= -Image-
```

![20x20 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFAAAABQCAYAAACOEfKtAAAAwUlEQVR42u3aSwqAMAxAwUS8/5XjXhBFSWth3gEUhjTgJyOiQq/bEAAECBCgAAIECFAAAQIEKIAAAQIUQIArtI++YVXvJ5jMNIGOMEA9Xhnd34Xvdt7XndV9fRPoCAO0A1faSaPvbwIBAvQs/Kdn0/P9up+9TSBAgAABCiBAgAAFECBAgAIIMLwPjKv3cbO/iZhARxig/BtjAh1hgBqzA/0jLYAA7UATKIAAAQIUQIAAAQogQIAABRAgQIAAEQCc2gHetiebwaWH1gAAAABJRU5ErkJggg==)

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
```

```mathematica
In[29]:= ImagePad[zone, 4]
Out[29]= -Image-
```

![40x40 result](data:image/png;base64,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)

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
