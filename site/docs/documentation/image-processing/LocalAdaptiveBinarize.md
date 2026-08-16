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
```

```mathematica
In[3]:= LocalAdaptiveBinarize[chk, 2]
Out[3]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzUlEQVR42u3aMQ6EMBAEwV7E/79sXkDgxARUxyeS0m1gzVSrjdba+nkzk++/d6VPAwAAgAAAEAAAAgBAAP7V7W3n2+/7BzhBAAQAgAAAEAAAAgBAB5vsgrILcoIEAIAAABAAAAIAQACyC8rbjl2QEyQAAAQAgAAAEAAAApBdkLcduyAnSAAACAAAAQAgAAAEILugvB3ZBTlBAgBAAAAIAAABACAA2QX5vn+AEyQAAAQAgAAAEAAAApBdUHZBcoIACAAAAQAgAAAEAIAO9QC292S1XLLeOgAAAABJRU5ErkJggg==)

```mathematica
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
```

```mathematica
In[16]:= LocalAdaptiveBinarize[disk, 1]
Out[16]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABAklEQVR42u3cgQ3CMAxFwW/U/VcOM0Sq5aa9twCIkyOrASrJisb6+QgAABAAAAIAQAAACAAAAQAgAAAEAIB6u572htbqvaCrKhMgAAAEAICS6v5e0F1bze72MvW6JsARJAAABCCeBY1uF1NbkwlwBAkAAAH42hY0tUVk6Iburi3OBDiCAAgAAL31WdDTvoez+z67tzsT4AgCIAAABACAAAAQAADKi54Fdd8onX6jZwIcQQAEAIBO34KmbpRy+I2eCXAEARAAAPrab8R2t6NTfvNlAhxBAgBAAOL/gvxroglwBAkAAAGwBckEABAAAAIAQAAACAAAAQAAwEcAAIAAANBMf4lJK8dMhxxOAAAAAElFTkSuQmCC)

```mathematica
In[17]:= LocalAdaptiveBinarize[ramp, 2]
Out[17]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABBElEQVR42u3cQQrDIBRF0Wfo/rdsVpCBYNGv545LBz0o8k3akvRMqPcpXzNcay07Nfo7PNHSAAAAIAD39vMTrD31WQG2IAACAEBVTkH/nvnsNtuxAgAIAAABACAAAAQAgADk4lnQ16zma0Y0+vkcevNlBdiCBACAAAAQAAACAEAActDT0VWe59ltBmUFAAAgAADkHbFazbrpswJsQQAEAICcglJq5uNNeVuQAAAQAKeglJrV7HZTZgXYggAIAAABACAAAAQAgHLxLGjWTKbK09pWAAAAAgBATkFnvwtmBQAQAAACAEAAAAgAAAEAIAAAAAgAAMWNWFb+b48VYAsSAAACcFkvlmwuwZs+ZgUAAAAASUVORK5CYII=)

```mathematica
In[18]:= LocalAdaptiveBinarize[zone, 2]
Out[18]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABoElEQVR42u3dUY7DIAyEYWe1978yfc1Dq02zEGz45gAVYTS/DAVztNZavNFxHNFDH34+tXp9+5V5+AmaKgbMTtsZQd9GryJeMmDqPG8SAEGbIygi2mzUjKg6qnyLBEAQBLWRUX0SL9kwdeXbJQCC9tZv9apm9PhHV00SAEGqoNYjbtWxM7o6+jQ/EgBBEHQbQTtgZzSOJACCLMRSYWfEwmfEmHst1iQAgvwp3zKdk1n1/I8qCILose3oKueIzuOZtaiUAAiCoK2wkwFHjiZCEP2JoJ23mp/cI5IACIKgbSufDBWRBEAQA4gBDCAGMIAYwABiAAOIAQwgBjCAGMAAYgADiAHhT/lI2uyierMRCYAgCIqslxdigxv0EgBBELRM+68qlY9rqhBEU66pZsNRhupOAhgQWpbdbTSxc7+g/4xZFQRBdLsKGr1HVGXfqRcqJQCCIKh0vx29owmCPOKzcKOP0fMgARCkClqyG2GVZiMSAEGqIK+pek0Vgigzgrwp7015CKIHXtDoFSvYuT4/EgBBENRWOf9T8byQBEDQ3noBqlavtZX9o4IAAAAASUVORK5CYII=)

```mathematica
In[19]:= LocalAdaptiveBinarize[noise, 3]
Out[19]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABNElEQVR42u3csa7CQAxE0THi/3/ZNAglRaBBWYjPrV7DSk/IV57VLJWk86T79eeOqsonjj7rnPfn3IKl+AIWc9+OBkWcf44JoKDZ1HYLoojzzzEBFGQLGqWIX/tfTAAF2YLaJrPuHBNAQcO3IIpYe44JoCBBjCIEMQqCICaIgYIEMVfTekEUBL2g6AWBgvSCKCJ6QRQEvaB4oAEKmnAXRBGCGAVBEMukq2kTQEGCmE0mekEUBL2g6AWBgvSCKEIQoyAIYoIYKEgQu6Ai9ILgC/i7IGaTiV4QBcEDjegFgYL0gqjmm2HQBFCQLYgiBDEKggcaM6+mTQAFCWJtk9ELoiDoBUUvCBTkOpoiBDEKggcaghgoSBBzNa0XREHQC4ofbgUFeaBBNdELoiB4oBEPNEBBl+EBT0AuwpQgQq0AAAAASUVORK5CYII=)

```mathematica
In[20]:= LocalAdaptiveBinarize[rgb, 1]
Out[20]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABNElEQVR42u3c0Q3DIBAE0SNy/y2TBuIPJAhw96aBWB7tZQWYFhE9sJzef7/mj1ezFwIIIAAE1OWZ9W++mtbaFa1GAowgEEAACMjagkZbyltbWN1qdv2uBBhBIIAAEJC1BZ3WInat7Yy2PgkwgkAAASAg61rQ6raQ9TklwAgiAAQQgKwtKIqddpYAIwgEEAACLqGd9o1YtTUiCTCCCAABBCAOWwtafU5m9HzRLTtrEkAACCAABNgRi7+s+cw6rf32PM4FGUEggAAQoAXlZteX+xJgBBEAAgjAaS1o1g7UaLu4fYdLAggAAQSAAGtBkeKOaAkgAAQQAAK0oEh926EEGEEggAAQoAWd2Uay7pRJAAEEgAACcPuX8lHszmcJMIJAAAEgIJLemljtlLIEEAACCAABxfgCn8dUwC9HRZMAAAAASUVORK5CYII=)

```mathematica
In[21]:= LocalAdaptiveBinarize[sky, 2]
Out[21]= -Image-
```

![24x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABACAYAAADlNHIOAAABGUlEQVR42u3byw4CIQxG4WJ8/1fGPRtCuBXznaXRyHjSvxZmSq21Bq7x8RMQQAAIIAAEEIDzfNsXSimpFtiOKb319caa0eubHZN636cCRBAByNQDRjN3d6bfXt/u61UBIogAZOoBsxnZZlzv/bOZO/r51T2kd73mABEEAhJTnAnH0r0jFSCCQIDzgH17K3H5PGJ27lEBIggEvDwHzPaI1ecBr88JKkAEEYDMc8DqzF29/569R/TWowJEEAG4+Tc2IurJey1n54R/O75QAQQQgMw9IB7bn1cBIIAA7OsBoxmc/Zmz23thKkAEEYBMPeD2/+zTPWP1c8gqQASBAM8H5Okh2e8jUgEiiAA4D1ABIIAAEEAACCAAB/kB1Gl4ZHGDFOoAAAAASUVORK5CYII=)

```mathematica
In[22]:= LocalAdaptiveBinarize[bit, 1]
Out[22]= -Image-
```

![8x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzUlEQVR42u3asQ2AMAwAwTdi/5XNDGmg4L5GaU64sDy7ux00Myefd/j8796/0qcBAABAAAAIAAABACAA/2qqtdv57n1/gBEEQAAACAAAAQAgAAD0YrfdTu6CjCABACAAAAQAgAAAEIDcBWW3k7sgI0gAAAgAAAEAIAAABCB3QdntuAsyggQAgAAAEAAAAgBAAHIXZHfkLsgIEgAAAgBAAAAIAAAByF2Q9/0BRpAAABAAAAIAQAAACEDugnIXJCMIgAAAEAAAAgBAAADopR6fA2WxoEgxUAAAAABJRU5ErkJggg==)

```mathematica
In[23]:= LocalAdaptiveBinarize[byte, 2]
Out[23]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABIUlEQVR42u3dQQ7DIAxEURP1/lcmJ2CBFAo276+rqMrXTK0AaYuIHtjG4xYQQAAIIAB7+LkF/6H3LgEqCAQQAAJST0GjX/PW2tTnR4yuM/t9Zq8/ex0JUEEggAAQ4FnQt1PNadeXABUEAggAAQSAAAJAAAEggAAQQAAIIAAEhBWxiJjfLxSH7VKWABUEAggAAZfRspyUP216kQACQAABICA8Cyo97aw+syYBKogAEEAAqq6IZWHXCXoJUEEEgAACcNuK2Orp4rT3AkmACgIBBICAqlNQ1X07q1fQJEAFEQACCEDYF1RiqpEAFQQCCAABt/2DhqlGAlQQCCAABNSYgrJMO9nfRyQBKogAEEAAdk1B2U+mZ9/FLQEqiAAQQAA28QIiQkC7JZkvkQAAAABJRU5ErkJggg==)

```mathematica
In[24]:= LocalAdaptiveBinarize[vol, 1]
Out[24]= -Image-
```

![12x10x8 volume](data:image/png;base64,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)

```mathematica
In[25]:= ImageChannels[LocalAdaptiveBinarize[rgb, 2]]
Out[25]= 1

In[26]:= ImageDimensions[LocalAdaptiveBinarize[vol, 1]]
Out[26]= {12, 10, 8}
```

```mathematica
In[27]:= LocalAdaptiveBinarize[chk, 4]
Out[27]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzUlEQVR42u3aMQ6EMBAEwV7E/79sXkDgxARUxyeS0m1gzVSrjdba+nkzk++/d6VPAwAAgAAAEAAAAgBAAP7V7W3n2+/7BzhBAAQAgAAAEAAAAgBAB5vsgrILcoIEAIAAABAAAAIAQACyC8rbjl2QEyQAAAQAgAAAEAAAApBdkLcduyAnSAAACAAAAQAgAAAEILugvB3ZBTlBAgBAAAAIAAABACAA2QX5vn+AEyQAAAQAgAAAEAAAApBdUHZBcoIACAAAAQAgAAAEAIAO9QC292S1XLLeOgAAAABJRU5ErkJggg==)

### Applications (6)

```mathematica
In[28]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[29]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[30]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[31]:= Binarize[LocalAdaptiveBinarize[noise, 2]]
Out[31]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABdElEQVR42u3cwQrCMBCE4Yn4/q8cLyIV1Co1adJ+c/LiUij7szPstiSpuavWx88nlVKypnf/VedznUtoV3kBO+u6bA2I6F9HB0DQuVWWU9BREfHNf3vW0QEQRC+noBFasvXzjIZNHQBBpqCfODJCXjTa82ypowMg6ORTkHxm3zo6AIIYMYhgxCCIRjZioul2dXQABJmCpouU0zGa/lfuxIhBEK1mQWeLgmMviLwAe0HnrqMDIIgRg4g40IAgmsSIiaYZMQiiRkZstGh6tK3pFtG0DoAgRuww286ZMCrXARBkCqrymTjQgCByoBEHGgRBjBhEMGIQRFPsBc34QY/RcicdAEGmoENOILNE0zoAgmRBTFMcaEAQOdCIAw2CIEYMIrrmVzoAghix6b6xc6QPeugACGLEqm872wuCILIXFHtBBEHiaIiIAw0IIgcajBhBECMmmrYXBEE09V6QaFoHQBBtNGJMU+wFQRC1N2IQwYhBEPXVDalMIs5V94w5AAAAAElFTkSuQmCC)

```mathematica
In[32]:= EdgeDetect[LocalAdaptiveBinarize[zone, 2]]
Out[32]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAB30lEQVR42u3dwVLFIAyFYer4/q+MG3VY2LFcoAm331m5aaXNnP+GQMNRSqnlW7X+/lmO4yg9aq9t1XufDFrxLGfv9qNQqAQgWEdtvHHFYmf2XI2d3f/v2X04AIIejqA2C4rKCi4NdBIe7xz/lftzAARBUJ2NnaiMJUNm1TuZ5QAIgqC6EjtPqwX1XssBEARBL5ej3wk7UTjiAAhSjq49FnsCdlbjCIIgiH70mQE7vWXkIeYOjLm9th3zyGYGDoAgE7E6u8T6tAX9K+/HojwE0b9Z0C7b/GZlPrMymZEMigMgSBZU78p8MteLop6XAyBIFhQy8SmblJ1XY5MDIAiC0pZ/b00HTyZKHABBJAACQAIgACQAAkBpJ2K7lJ1HJmUrnosDIAiCtrNtUY4mARAAEgABIAEQABIAE7H32prIARBEWyMoao/NjnuZOACCIOi28nK2jGgWHkf6LHEABEFQyAL3nZ+vZmgk60t5CKI/ETTyC95r294v1jNPskbGqWUZBJGuibomkgDoHa13NAdAEAQ5QcMJGhBE79SyzDliWpZBEJV5K2K9v+xRJ59mO02VAyCIprSvzzDhyryw7kx5CKLbVsSiso5s2xSVoyGIwlfEdjlBI2psV94nB0DQs/UF0HlCwrRV3TMAAAAASUVORK5CYII=)

```mathematica
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
```

```mathematica
In[41]:= LocalAdaptiveBinarize[zone, 4]
Out[41]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABh0lEQVR42u3dwW7EIAyEYafq+7+ye+khl5W2zRoMfPMAEWE0vxwC5srMjF9d1xX/1e0x26pifr6CpooBs1MVEQk783AkARAEQTkbNU/i3A2bf30XCYAgCMpRUR2Jmg5oeud9JQCCztb36lVN9ZirqyYJgCBVUH4ibitiZ2Sl9Gp+JACCDkdQPsjVCdipxpEEQJAPsS2XjqvHeX/+k3FKAARB0JZLwa+e2a1ykwAG+BDLURVF531Es95RAiCIAcQABlC3KuiEyqdDRSQBEMQAYgADKDZajl797Nh9/NXL1xIAQQwgBjCAGMAAYgADiAEMIAYwgBgQlqNj8uGFWPynvARAEDGAAcQABzQc0JAACCIGMIAYsHq/oIp1kpHVUYfx6xcEQRSdmnVUNE1dpZ2aBDAAgtoeXjih3b0EQJDe0XpH6x0NQbTTJT474ah6HiQAgnyIHds0tcOGAQmAIFWQ21TdpgpB5E55d8oTBB2+L+hTsYKd9+dHAiAIgnLHLYirbFmUAAg6Wz8DPqycYJELmAAAAABJRU5ErkJggg==)

```mathematica
In[42]:= LocalAdaptiveBinarize[zone, 1]
Out[42]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABlElEQVR42u3dQQ7DIAxEUafq/a9Mt1k0Ummh2OHNBaIwmi8Tgn201lq80XEcMVoXj0qhVe/7CFoqBqxO3hlBn8QwM0aq4Ou8hhIAQZunJyLaCuzMqDoyVGW9GJcACIKgNjOeGVCzCk2fvLsEQBAENdhZhyMJgKC99bzTxmr2d54ZGzcJgCBV0NdV0F2xM7s6Oq+bBECQKqjE4XW2am1UdSQBEARBt8fO1XMz4EgCIGhzBP0TBdn+KZqBIwmAIBpSBY2KZJVfGUfhqLcikgAIgqCtLmJkq44kAIIgiBYiTgIgiAHEAAYQAxhADGAAMYABxAAGEAMYQAxgADGAAcSAcChPMeWemgRAEASVa3BRsY+0BEAQdSFo1YWF6thxTRWCqAtB/0RENhxl6HckAQxQBZVrcFFxWI8EQBAN+RY0G0exWU9pCYAgCCrdb6dis1YJgCAafiJ2VxzN3iRKAASF9vXmiJkjBkFkmqppqgRBZsqbKW+mPARR1BviswOafqncrtZHAiBIFdQybZqqb/r8HQ1B1KMXyGy4u/PLgMIAAAAASUVORK5CYII=)

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [Binarize](../../image-processing/Binarize/), [ImagePad](../../image-processing/ImagePad/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
