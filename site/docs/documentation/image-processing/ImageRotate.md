# ImageRotate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageRotate[image] rotates a quarter turn counterclockwise; ImageRotate[image, angle] rotates by angle in radians (use n Degree for degrees). A multiple of a right angle takes an EXACT index-permutation path -- every pixel lands on another pixel's position, nothing is interpolated, and four quarter turns are exactly the identity. An odd number of quarter turns swaps the dimensions. Any other angle interpolates bilinearly, sampling the source per destination pixel (inverse mapping, so every output is filled exactly once; forward mapping leaves holes wherever the rotation stretches). Area rotated in from outside reads as 0 rather than the replicated edge, because that area was never photographed and smearing the border across it would invent content.`**

## Examples (39)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (8)

```mathematica
In[1]:= ImageData[ImageRotate[Image[{{1., 2.}, {3., 4.}}]]]
Out[1]= {{3.0, 1.0}, {4.0, 2.0}}

In[2]:= Module[{img = Image[{{1., 2.}, {3., 4.}}]}, Nest[ImageRotate, img, 4] === img]
Out[2]= True

In[3]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[5]:= ImageRotate[chk]
Out[5]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzUlEQVR42u3asQ2AMAwAwTdi/5XNDGmg4L5GaU64sDy7ux00Myefd/j8796/0qcBAABAAAAIAAABACAA/2qqtdv57n1/gBEEQAAACAAAAQAgAAD0YrfdTu6CjCABACAAAAQAgAAAEIDcBWW3k7sgI0gAAAgAAAEAIAAABCB3QdntuAsyggQAgAAAEAAAAgBAAHIXZHfkLsgIEgAAAgBAAAAIAAAByF2Q9/0BRpAAABAAAAIAQAAACEDugnIXJCMIgAAAEAAAAgBAAADopR6fA2WxoEgxUAAAAABJRU5ErkJggg==)

```mathematica
In[6]:= ImageRotate[disk]
Out[6]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzklEQVR42u3ciw2AIBBEwT1j/y1jFSuQzKvAODlC8DNJVrStxy0AAEAAAAgAAAEAIAAABACAAAAQAADq9p52QWt1H9DNjAkQAAACAEDJtN8Lau9q6jeovGsyAZYgAAIAQAAACAAAAQCgXHgWdPuZz64zIhNgCQIgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIDiGzHfiAkAAAEAIAAABACAAABQzj0Lir8mmgBLkAAAEAC7IJkAAAIAQAAACAAAAQAAQAAACAAAAQCgv/sAoIQVtzrUII4AAAAASUVORK5CYII=)

```mathematica
In[7]:= ImageDimensions[ImageRotate[chk]]
Out[7]= {16, 16}
```

```mathematica
In[8]:= ImageRotate[chk, 0.4]
Out[8]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAEQElEQVR42u2dyyu1YRTF9+tSOimlUDoZMDCQCCmRIpQSyW3GSEKRiXLL5TBgIkJk4jIiBm4pxACFUChlYKIkjJRE4vgH1h7swZdzPmsNf7znlGXvVvu5vI6IeMUgj8cDeWdnJ+Rzc3OQn5+fQ56bmwt5SUkJ5C6XC/Ln52fIFxYWIG9qaoL89fVV/qUChPpV0QAaQAMoGvB3FWR94PT0FPKdnR3IZ2ZmIN/b24M8Pj4e8vX1dcgPDw9N6cjtdpvSTmRkJOSFhYWQz87OsgLYgigaQAMoGuAPcsQ4C9KUnp4O+fj4uCkdDQ0NmVJNUlIS5O/v75BXVlZCPjExIZYZUXd3N+STk5OQ19fXswLYgigaQAMoGiD+MAtaWlqC/OrqCvKIiAjI09LSIN/f34c8Ojoa8sbGRsg/Pz9NXFvR6+3tNaWgk5MTyJeXl1kBbEEUDaABFA3wixSkzVLKysogz8nJMaWgo6MjyD8+PiB/eXmBPDU1FfKBgQGxzGrW1tZM36vNso6Pj8Uyy9LSHSuALYgGUDSABlC/lYIuLi7gDwoKCiDPzs6G/Pv725QuVldXTSnl8fERcq8XL+hNTU1B3tfXB/ni4iLkLS0tkMfFxZlS4uXlJSuALYiiATSAogG+JMflcsEY8fb2Bh+4ubmBPDg4GPL5+XnI6+rqIB8eHoZ8cHAQ8traWlNK0XZTa6lMS3HabKempgby+/t7VgBbEEUDaABFA3xJQdqJ8pGREbHsOv76+hLLylRPT49YZjhZWVmQV1dXm35fS03aCqC2EqelxLu7O7GccWMFsAXRAIoG0ADql+R4lSWloqIi+MDGxoZY9uecnZ1BHh4eDvn09DTkBwcHkLe1tZnSkTaD0lbWrq+vxbL7Wvv7NDc3swLYgigaQAMoGuBLCtJuI4yJiRHLiXXtDFdYWBjkxcXFYln5ur29NaWjzMxMyB8eHsQygwoJCYE8MDAQ8vb2dsj7+/tZAWxBFA2gARQN8KkUpN3VrM0uEhISxHLSPCMjQyy7oLUT6Fp60e4pioqKMnHtNkhtRU87cZ+fnw/59vY2K4AtiKIBNICiAT6VgrQ3WWh3O2u3HWopIiUlRSy7kbXPcRxHLCtr2qypoqJCLDOrgAD8P9rQ0CCWN4kkJyezAtiCKBpAAygaIP5wd3ReXp5Y9v9o7wsbHR01pYvY2FixrEyVlpaK5Z4fbbajnYjXPkfbVZ6YmAh5eXk5K4AtiKIBNICiAeIPd0drZ6a0FajNzU2x3KaozVK0E/ca11JZR0cH5K2trWKZQWnpTlu5YwWwBVE0gAZQNED+0/eIaffhaCno6ekJ8q2tLbHsq9F2HWu7oLUT61r60m471L53bGwM8pWVFVYAWxBFA2gARQP+1NtUNYWGhoplpayqqkosMyXt3h4tjezu7opl/5J24l4789XV1cUKYAuiaAANoGiAP+gHKRdjfSzf4qEAAAAASUVORK5CYII=)

### Scope (19)

```mathematica
In[9]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[10]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[12]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[13]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[14]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[15]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[16]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];
```

```mathematica
In[17]:= ImageRotate[rgb]
Out[17]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAyUlEQVR42u3YMQ6AIBBFwcXQw5E9MpwALgC1EOeVVsbJTzamUd8Ri1qsO+15v+Q9d8+f0KcBAABAAP5bjuojWAAAAQAgAK4gWQAAAQAgAK4gWQAAAQAgAK4gWQAAAQAgAK4gWQAAAQAgAK4gWQAAAQAgAK4gWQAAAQAgAFddQcVHsAAAAgBAAPwLkgUAEAAAAuAKkgUAEAAAAuAKkgUAEAAAAuAKkgUAEAAAAuAKkgUAEAAAAuAKkgUAEAAAAuAKkgUAEAAAAnB4EyrOEmwBM6EnAAAAAElFTkSuQmCC)

```mathematica
In[18]:= ImageRotate[sky]
Out[18]= -Image-
```

![16x24 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEAAAABgCAYAAACtxXToAAAAr0lEQVR42u3QMQqAMAwF0OhiwZ7DO3l0r2IHXUsFIaP4/hZ+AuFN235c0WVdzn6MWlq89c+5pfpasvfJ/eH/OvRz/DwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAOALuQE/0RLXzgUb7QAAAABJRU5ErkJggg==)

```mathematica
In[19]:= ImageRotate[bit]
Out[19]= -Image-
```

![8x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzUlEQVR42u3aMQ6EMBAEwV7E/79sXkDgxARUxyeS0m1gzVSrjdba+nkzk++/d6VPAwAAgAAAEAAAAgBAAP7V7W3n2+/7BzhBAAQAgAAAEAAAAgBAB5vsgrILcoIEAIAAABAAAAIAQACyC8rbjl2QEyQAAAQAgAAAEAAAApBdkLcduyAnSAAACAAAAQAgAAAEILugvB3ZBTlBAgBAAAAIAAABACAA2QX5vn+AEyQAAAQAgAAAEAAAApBdUHZBcoIACAAAAQAgAAAEAIAO9QC292S1XLLeOgAAAABJRU5ErkJggg==)

```mathematica
In[20]:= ImageRotate[byte]
Out[20]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA+ElEQVR42u3cMQ6DMBAF0QUh+f6n80HcbQMXMKWxxL4pSTnia7RKcvTe75gwxgjP1z8/A1shgAACQEBdrtba9IPM9PyD594AE0QACCAAKkgFwQQRAAIIgApSQTBBBIAAAqCCVBBMEAEggACoIBUEE0QACCAAqyvo7TdN6kUFmSAQQAAI+HcFuRGpIBMEAggAASoo3IhUkAkCAQSAgPLfC3IjUkEmCAQQAAJUULgRqSATBAIIAAHlfyPmRqSCTBAIIAAEqKBwI1JBJggEEAACyv9fkBuRCjJBIIAAEKCCwo3IG2CCQAABIKD8f0e7EXkDTBAIIIAAbOUBpV2plvXO4+oAAAAASUVORK5CYII=)

```mathematica
In[21]:= ImageRotate[zone]
Out[21]= -Image-
```

![32x32 result](data:image/png;base64,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)

```mathematica
In[22]:= ImageRotate[ramp]
Out[22]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAuUlEQVR42u3WwQkAIAwEwSg2kP6LPauQPJwtIcNBVlWlNNZ2AgAABACAZjrd7QoWAEAAAAiAL0gWAEAAAAiAL0gWAEAAAAiAL0gWAEAAAAiAL0gWAEAAAAiAL0gWAEAAAAiAL0gWAEAAAAiAL0gWAEAAAAiAL0gWAEAAAAiAL0gWAEAAAAiAL0gWAACAEwAAIF+QBQgAAAEAIF+QBQgAAAEAIF+QBQgAAAEAoDetJHEGCwAgAAAE4LsuxOYIoKDv3ygAAAAASUVORK5CYII=)

```mathematica
In[23]:= ImageRotate[noise, 0.8]
Out[23]= -Image-
```

![32x32 result](data:image/png;base64,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)

```mathematica
In[24]:= ImageRotate[disk, 1.2]
Out[24]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAB9klEQVR42u3cMa6hYRSH8fdOJqJRyt0BorIBlsAGRHQqrU5Dd5U6pdgAS2ADKrnXFpQa0cxs4LyFiZk7+D3lv+J9ck5Ovu8931tK6VfCt/HDERBAAAggAAQQAAIIAAEE4N/x81F+aKPRCPPD4RDm1+s1zEulUpifz2cVoAWBAAJAgCnojhSLxTDv9Xph3u12w7zZbIb5ZrMJ81qtFubr9TrMR6ORCtCCQAABIOCVeEt/+V7QeDwO88lk8i1/+HQ6hfl8Pg/z2WwW5pfLRQVoQSCAABDwUlNQ7o1Sbtrp9/thXi6XH+KABoNBmC8WCxWgBYEAAkCAZ0Eppa+vrzCvVCoPfUC73S7MW62WCtCCQAABICC90r2gQqEQ5p+fn085BeXuI+Vua+/3exWgBYEAAkDAU05B9Xo9zNvt9lMeUO72dW43TQVoQSCAABCQXuyN2Ha7Tbc8S3kUjsdjmFerVRWgBYEAAkBAsimf0mq1eugpKLc7ltusv9d3h1SAFkQACCAAjz4FLZfLMH9/fw/z4XCY/qfdsdym/HQ6VQFaEAggAAQkb8Tux8fHR5h3Op10y+3r3H2k3G5X7plVboq713eBVIAWBAIIAAGmoPRHb5Ryu2m529q37mqpAC0IBBAAAkxBUAEEgAACQAABIIAAEEAACHhSfgNtuGvq8F5rfgAAAABJRU5ErkJggg==)

```mathematica
In[25]:= ImageRotate[zone, 0.3]
Out[25]= -Image-
```

![32x32 result](data:image/png;base64,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)

```mathematica
In[26]:= ImageChannels[ImageRotate[rgb]]
Out[26]= 3

In[27]:= ImageDimensions[ImageRotate[sky]]
Out[27]= {16, 24}
```

### Applications (4)

```mathematica
In[28]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[29]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];
```

```mathematica
In[30]:= EdgeDetect[ImageRotate[chk]]
Out[30]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA9ElEQVR42u3dSQ6DMBBFwd+I+1/ZuQC9AEHMUG+ZVUSpHeQYUUlGNhpj8+NUVdS397otLtncAAAAIADfbe1+tfWfuyMTYAkCIAAANKnq9oJkAgAIAAABACAAAAQAgAAAEAAAAgBAAF7R6hLEuSBLkAAAEIA4FyQTAEAAAAhAPv6MmCfiTQAAAQAgALEXJBMAQAAACECcC1KcC7IECQAAAYi9IJkAAAIAQADiXJBMAAABACAAsRckEwBAAAAIAAABACAAAAQgD3+bqn/KcugUtAmwBAkAAAG4211Qd1dz1rvmrz539JRzTd33MQGWIAACAECT+gHo1C234WvMQQAAAABJRU5ErkJggg==)

```mathematica
In[31]:= Binarize[ImageRotate[zone, 0.5]]
Out[31]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABmklEQVR42u3dW27DMAxEUano/rfM/vYnaFLrQYlnFuBYGMwF5VBUb61FG6CIIY9Zqt779nf4asSAyvquvPhX2FyJJgmAoNrqT6qgUZXPysg/eecZ7ykBEARByxCUYeOTbV0SAEEQFLCzb70SAEEQFCNieDpqdm3iJACCIChO+Gx7K44kAIIgKE6ofDJvAJ/gSAIgCIJid+WzsqdoNpo+XYsEQFDTF3Q7dl79rtZEYkCp1sRsHdQzcPT7Oe+sVwIgqDiCTjxYcVN1JAEQpAqauvmCOAlgADkjlrYikgAIYgAxgAHEAAYQAxhADGAAMYABxIBW+3N0tha+duk0RQmAIAhaFmd/0EsAA2jjn/KZcbSrD0oCIAiCjjjIcNP4esM6IIhStSbOxlG20/ESAEH0bwTN/kx966dvUxMhiP6cF5TtdolWYKa0BECQKujKOTynTG6UAAhSBcUtd4SdeO+YBEAQBMXt1xdmXpcEQBAERaZZQCcOiXWVIQTRFgRVHlM2CpUSAEEQFFVnR7tTnhhQCkGnoMkxVQiiVfoBIMKlvQPRNR0AAAAASUVORK5CYII=)

### Properties & Relations (6)

```mathematica
In[32]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[33]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[34]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[35]:= ImageData[ImageRotate[ImageRotate[ImageRotate[ImageRotate[chk]]]]] === ImageData[chk]
Out[35]= True

In[36]:= ImageChannels[ImageRotate[rgb]] === ImageChannels[rgb]
Out[36]= True

In[37]:= ImageData[ImageRotate[disk, 0.]] === ImageData[disk]
Out[37]= True
```

### Neat Examples (2)

```mathematica
In[38]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];
```

```mathematica
In[39]:= ImageRotate[zone, 0.7]
Out[39]= -Image-
```

![32x32 result](data:image/png;base64,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)

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
