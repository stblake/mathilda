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

![14x12x10 volume](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAQQAAAEECAYAAADOCEoKAAADCklEQVR42u3c0RXEIAhFQbD/nk0XinkzFUSFe/ZrqwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAqHYFXLTNpCAgAuZSEBAB8ykICIA5FQREwLwKAgJgZgUBETC7goAAmGFBQATMsiAgAuZZEBAAcy0IiID5dmEIgDl3UYiAWcclCYCZx+WIAGbfpYgAdsBlCAD2wEWIAPbBBQgA9sLBRQD74cACgB1xWBHArjikAGBnHE4EsDsOJQLYHwcSAOyRg4gA9skBBABh8OEigCj4aAFAGHysCCAMPlIE4KWdawEAe/dqEEQAYQgOggAgDOFBEAGEITgIAgBDdnKJAfiFPCUIwCCCAAgCIAiAIACCAAgCIAiAIACCAAgCIAiAIACCAAgCIAiAIACCAAgCIAiAIACCUP4OGwSh/B02CAIgCIAgAIIACAIgCIAgAAgCIAiAIACCAAgCIAiAIACCAAgCIAiAIACCAAgCIAiAIACCAAgCIAiAIACCAAgCIAiAIACCAAgCIAiAIACCAAgCgCAAggAIAiAIgCAAggAIAiAIgCAAggAIAiAIgCAAggAIAiAIgCAAggAIAiAIgCAAggAIAiAIgCAAggAIAiAIgCAACAIgCIAgAIIACAIgCIAgAIIACAIgCIAgAIIACAIgCIAgAIIACAIgCIAgAIIACAIgCIAgAIIACAIgCIAgAIIAIAiAIACCAAgCIAiAIACCAAgCIAjA/4LQngDm7MSUhdzmAAHwMcKACPgwYUAAfKgwIAI+WhQQAIcQBkTAgYQBEXA4YcCOOKwoYC8cXBiwBy5CGDD7LkUYMOsuSRQw3y5MGDDTLk8YMMMuM48wmFtcrCiYU1y0MJhNXLowmEU8gjCYPzyIKJg5PI4wmDE8ljCYKzycMJgjPKQomB08qjCYFTyyMJgPPLgwmAcMgCiYAQyDMHh3DIYweGcEQRi8LYLA4Sh4SwxReBi8HwYqOAzeC0EID4M3QhDCo+BdAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADglg/E2E9KPMbN6gAAAABJRU5ErkJggg==)

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
