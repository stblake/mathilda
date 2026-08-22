# ImageCrop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageCrop[image, {w, h}] crops to w x h about the centre, any odd remainder going to the right and bottom -- the same floor-division convention the kernel centres use, which is what makes ImageCrop[ImagePad[image, m], ImageDimensions[image]] exactly the original image. A crop may not enlarge. ImageCrop[image] instead TRIMS A UNIFORM BORDER, asking how much of the frame carries no information; the border colour is read from a corner rather than assumed black, since a scanned page's margin is white. An entirely uniform image comes back unchanged, there being no content to keep and a zero-sized image not being one.`**

## Examples (34)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (7)

```mathematica
In[1]:= ImageData[ImageCrop[Image[{{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}}], {1, 1}]]
Out[1]= {{5.0}}

In[2]:= ImageDimensions[ImageCrop[Image[{{0., 0., 0.}, {0., 0.5, 0.}, {0., 0., 0.}}]]]
Out[2]= {1, 1}

In[3]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[5]:= ImageCrop[chk, {8, 8}]
Out[5]= -Image-
```

![8x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAu0lEQVR42u3cwQnAMAwEwbuQ/ltWqhBOwmwBxjBYL6MmmSw2s3p82ubL97+iowEAAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAjAP7v92zl7fy/ACAIgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAWauwLsi/ICBIAAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAQABiX1DsC/ICjCABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAATg3T0wnxq7Fx7xeAAAAABJRU5ErkJggg==)

```mathematica
In[6]:= ImageDimensions[ImageCrop[chk, {8, 8}]]
Out[6]= {8, 8}
```

```mathematica
In[7]:= ImageCrop[disk, {12, 12}]
Out[7]= -Image-
```

![12x12 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA0UlEQVR42u3dgQnDMAxFwa+Q/VdWhjCODL43QckhUdeQVpKOxno8AgAABACAAADQ/73TH6B79hhSVSbAChIAAAIAQABuq1bvA6a/x48/wMVzhAmwggAIAAABACAAAAQAgAAAEAAAygH3Abf/3r/7vsAEWEEABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAMX7grwvSFYQAAEAIAAABACAAAAQAAACAEDZfx8Q/yNmAqwgAQAgAAAEwDlAJgCAAAAQAAACcHofh1wVtxnbmh4AAAAASUVORK5CYII=)

### Scope (16)

```mathematica
In[8]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[9]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[10]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[11]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[13]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[14]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];
```

```mathematica
In[15]:= ImageCrop[rgb, {8, 8}]
Out[15]= -Image-
```

![8x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAwElEQVR42u3aMQ6AIBREwcXQw9E4kjeXM1D5jfNqCsNkG2Jb635y0kxKnR/Fvufw/BW9GgAAAAQAgAAAEAAAAvCv+nQHFgBAAAAIAAABACAAAAQg3oJkAQAEAIAAABAAAAIAQADiLUgWAEAAAAgAAAEAIAAABOCb9eEOLACAAAAQAAACAEAAAAhA/BckCwAgAAAEAIAAABAAAAIQb0GyAAACAEAAAAgAAAEAIADxFiQLACAAAAQAgAAAEAAAAlC6DYjtA8t8eg1wAAAAAElFTkSuQmCC)

```mathematica
In[16]:= ImageCrop[sky, {12, 8}]
Out[16]= -Image-
```

![12x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABACAYAAADlNHIOAAAAk0lEQVR42u3bsQ2AIABFQTHC/j2VnQNKghVsIQHujcDlF4QQ7ie3Q8M6HQEAAAIAQAAA6P+u8iWnYAEABACAAAAQAPcAWQAAAQAgAAAEYPl7QHQKFgBAAAAIAAAB2Owe8FbvARYAQAAACAAAAfAeIAsAIAAABACAAPgfIAsAIAAABACAAHgPkAUAEAAAAgBAACavAy+VERhpQTLJAAAAAElFTkSuQmCC)

```mathematica
In[17]:= ImageCrop[bit, {4, 4}]
Out[17]= -Image-
```

![4x4 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAvUlEQVR42u3cwQnAMAwEwVNI/y1fqhA2YbYAYxisl9G0bRabmc3js3z99fs/0dEAAAAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAE4J9Nkvq3c+7+XoARBEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAC30+rdjX5ARJAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAMS+oNgX5AUYQQIAQAAACAAAAQAgAAAEAIAAABAAAAIAQAAACMDlfRirG7fp1IJFAAAAAElFTkSuQmCC)

```mathematica
In[18]:= ImageCrop[byte, {8, 8}]
Out[18]= -Image-
```

![8x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAyElEQVR42u3YMQqAMBQFwVWEVLll7n8FPYMQsHC2TuXgL96x1rp70Rgj7/e9P9OnAQAAQAAACAAAAQAgAP/qmnNmq8kW5AQJAAABACAAAAQAgABkC8q2YwtyggQAgAAAEAAAAgBAALIFeZ8tyAkSAAACAEAAAAgAAAHIFuS9P8AJEgAAAgBAAAAIAAAByBbkvT/ACRIAAAIAQAAACAAAAD5BtqBsQXKCAAgAAAEAIAAABACAsgVlC5ITBEAAAAgAAAEAIAAAtK8HXqgIounc3X4AAAAASUVORK5CYII=)

```mathematica
In[19]:= ImageCrop[zone, {16, 16}]
Out[19]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAADJElEQVR42u2dyaoqQRBEy6teR1ScFoIouvD/f0gUwVlRcR7fD2QsEnx0wz2xjIVWZ3QmQVZXVqLdbn+CgU6nY9Gh2+2afLvdNvlms2nylUrF5HO5nMknk0mTf71eJn+5XEx+t9uZ/HK5NPnJZGLyo9HI5MfjsckvFguT/wkgUiAAAiAAQIC/i5RyO4PBwOR7vZ7LBdXrdZMvlUomn81mv+KCrteryR8OB5Ov1WomXywWTT6TybgC/Xw+yQBKEEAABAAIECsXpHo7yu30+32Tb7VaJl+tVk2+UCiY/O/vr/2m/Njvyvv9Nvn7/R48PSjldtR6FG63m8mfz2cygBIEEAABAALEygWpHo7ildtRO1/lcjl4ej7pdNrkE4mEyX8+5oZeeDwerv9NpVKuwCmXdTweTX673ZIBlCCAAAgAECBWLki5F7WTpXo7yu14ez5q58vrgpSrUb8fnDtZyu2oeDYaDTKAEgQQAAEAAsTKBakdIvXdjnI1qsei3I5yKWrnywvlmhTU90XqeVV8VDwVTwZQghAAIAACgKhckDqT5XU1aidL9V6U2/G6F/lmid9X61HrV8+r4qPimc/nyQBKEEAABAAIECsXpFzBt9zLt1zNt+Bdp9dNuePJO0gJQgCAAAgAInJBaidI8epMlvo+R/FRwbtO9bzeuCmeDKAEIQBAAAQAUbkgNV1QzdtRZ6PUmSz1/Y+39+KF172o9avnVfFR8eSkPCUIIAACAASImwtSs5TVdEH1la932qHCt86IKbfjdTWn08kVHxVPxZMBlCAEAAiAACAqF6RujvDOUvbO21Eu5X/PC1JuZ7/fB8+cn/V6HTw3caxWKzKAEgQQAAEAAsTKBal7sr41S1nN24lqdrTq7Si3M51Og+d+McXP53MygBIEEAABAALEygWpW0G992R5ZynH7R4x1dtRrmY4HAbPLauz2YwMoAQBBEAAgACxckHqDvTgvCfLO0v5r90pv9lsyABKEEAABAAIECsXtFgsgmcnS511UjtK6uYI5YLUdEGvC1LrVC5IfbejdrJUb0e5HdUrIwMoQQgAEAABQET4B16ETk14Z7siAAAAAElFTkSuQmCC)

```mathematica
In[20]:= ImageCrop[noise, {16, 24}]
Out[20]= -Image-
```

![16x24 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEAAAABgCAYAAACtxXToAAABsUlEQVR42u2cLa6EMBRGmYdDI1kBP3LWgEMjkWgca8CxBiQaNxKNnC2MZAd37HvPEQptcw8OMU3mS84h/XLbR13XEvx6siwLjrzneR6c+f3R9dM0Nbr+T6D8IQDtATyWZZEzjJ9l8m6H4AAQIIC/DojjWGwyeNYpZx0CAgSg3QFd18kZJm1/x4++//8/IEAA2h2wbZvYZND0+kfXAwEC0O6AoijE5e/21R0hCBCAdgcMwyAmGfetIwQBAtDugM/nIy7v76/uCEGAALQ7oCxLMcmkbx0hCBCAdgdM0yQu9fR3d4QgQADaHRCGoVw5h+d6RwgCBKDdAU3TXDor7HpHCAIEoN0Br9dLfNrfm+4IQYAAtDsgSRIxuZ/3rSMEAQLQ7oC+78XlWd6rO0IQIADtDni/33LnWV3XOkIQIADtDng+n7fOCrvWEYIAAWh3wDiO4tN5f9MdIQgQgHYH7PtudFbYt44QBAhAuwOqqhKfzvuzFwABAjDrgHmexeZ9frY7QhAgAO0OiKLI6l1itjtCECAA7Q5o21Zs3udnuyMEAQLQ7oB1XS+dFXa9IwQBAlD+fAFryvB8xqIf/AAAAABJRU5ErkJggg==)

```mathematica
In[21]:= ImageCrop[ramp, {8, 16}]
Out[21]= -Image-
```

![8x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADAAAABgCAYAAABbjPFwAAAAcElEQVR42u3PMREAIAwEwYCD9xT/lsAAJU1m9srrdnX3qUdJasLfNTwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACA/12fMQPwbQB2NQAAAABJRU5ErkJggg==)

```mathematica
In[22]:= ImageChannels[ImageCrop[rgb, {8, 8}]]
Out[22]= 3

In[23]:= ImageDimensions[ImageCrop[zone, {20, 10}]]
Out[23]= {20, 10}
```

### Applications (4)

```mathematica
In[24]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[25]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];
```

```mathematica
In[26]:= Binarize[ImageCrop[zone, {16, 16}]]
Out[26]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA3ElEQVR42u3d0Q2AIAxF0Vf337ksoRTMuQv4cVLSmCiVpDNQ98hjU1U5qScCAEAAAGhmKXhrC5raam7fmkyAIwiAAADQLVvQX7edqe3IBDiCAAgAAJ22Bdl29mxHJsARBEAAAAgAAAEAIAAAtLFqL31G3xGZAEcQAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAACKL+V9KS9HEAABACAAAAQAgAAAUPw72r+jBQCAAACQe8TcIyYAAAQAgL7dgtwpbwIcQQIAQABsQTIBAAQAgAAc3wI6CCWteXz19gAAAABJRU5ErkJggg==)

```mathematica
In[27]:= EdgeDetect[ImageCrop[disk, {12, 12}]]
Out[27]= -Image-
```

![12x12 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA6klEQVR42u3d2Q2EMAxAQRvRf8uhBD6sgA3zGthj5CjRApsRsUKvdfgKAAAQAAACAEDPd+5+gbVmHzMy0wRYggQAgAAAEICvldXfA6r7/N377O7v3wRYggAIAAABAKCO54Dp+/zu5wQTYAkCIAAABACAYuB1QdP3+dXPVz0nmABLEAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBA0eD+gLvr4/9+n7AJsAQJAAABACAA4bmhnhtqAixBAgBAAAAIwAfOAf5HzARYggQAgAAAEADnAJkAAAIAQAAACED7LnpfJ7uLft2WAAAAAElFTkSuQmCC)

### Properties & Relations (5)

```mathematica
In[28]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[29]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[30]:= ImageDimensions[ImageCrop[chk, {8, 8}]] === {8, 8}
Out[30]= True

In[31]:= ImageChannels[ImageCrop[rgb, {4, 4}]] === ImageChannels[rgb]
Out[31]= True

In[32]:= ImageDimensions[ImageCrop[chk, {16, 16}]] === ImageDimensions[chk]
Out[32]= True
```

### Neat Examples (2)

```mathematica
In[33]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];
```

```mathematica
In[34]:= ImageCrop[zone, {24, 8}]
Out[34]= -Image-
```

![24x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAAAgCAYAAADtwH1UAAAB50lEQVR42u2aSasCMRCEEx33bcCreBK8+v9/hVfBk0c3dMR9e1erhDQBxYFXdfsQM93pSfVMMt57/3Qvqtfrr+gGgwHwaDQK8nA4BO73+8Ddbhe40WgAl8tl4GKxCOy9B34+IXx3v9+BL5cL8H6/B16tVsCz2Qx4MpkAj8fjIE+nU+DD4QDc6/WAC076qVQAFeB/K2EPrVQqwK1WCzhN0yC32+2gx1er1aDnJ0mCd0gh7h7hHsHiHsHxcfxWvjw/PH/cc7gnaAXIglQA6Zc9gD2zVCoB12o1F3pP4N8tj+fx+TmfPd/y9Lc7iv7P4/P1OT6OPzZ/Hp/jv91uWgGyIEkFyE0PiPXQWLY8PdbjY2Vd79v5sh6Ph1aALEhSAXLbA9ijeO8klnk83nti/rSs6307X6vHagXIglQAKU/nAdfrFfh4PLrQfjb/fjqdXOhMlsfn/X/rOd2S5el8fY6P44/Nn8fn+X0779A9KAtSAaQcnQecz2fg3W4HvNlsgpxlmQudofJ+O++lsD79XRB7PJ/ZcvxWvjw/PH8cL58naAXIglQA6Zc9gJ9L2cPY8xaLBfB8Pnehbz+bzWbwOdgZ3+1YZ6zWewx7/na7BV6v18DL5TKYH+fP88Pzx/l2Oh2tAFmQpALkRX+SfvONqIRlcwAAAABJRU5ErkJggg==)

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
