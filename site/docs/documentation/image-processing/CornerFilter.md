# CornerFilter

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`CornerFilter[image] gives the corner strength at every pixel, from the eigenvalues of the Gaussian-weighted second-moment matrix of the gradient (the structure tensor). Both eigenvalues small is flat, one large is an edge, both large is a corner. CornerFilter[image, r] sets the window radius (default 2); CornerFilter[image, r, method] selects "MinimumEigenvalue" (the default -- Shi-Tomasi's lambda_min, which is directly "how much does the weaker direction vary" and is comparable across images) or "Harris" (det - 0.04 trace^2, cheaper since it needs no square root, and negative on edges). A STRAIGHT EDGE SCORES ZERO under both: every gradient in the window is parallel, so the matrix has rank 1 and its determinant and smaller eigenvalue vanish. Colour is reduced to luminance first, since a corner is a property of brightness.`**

## Examples (49)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (7)

```mathematica
In[1]:= Max[Flatten[ImageData[CornerFilter[Image[Table[If[j <= 4, 0., 1.], {i, 8}, {j, 8}]]]]]]
Out[1]= 0.0

In[2]:= Options[CornerFilter]
Out[2]= {Method -> "MinimumEigenvalue"}

In[3]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[4]:= CornerFilter[chk]
Out[4]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABT0lEQVR42u3dvQ6CMBiGURAJAwn3f6MMhD9Xh+9lxeB5xg4m9khtamLbaZrOpmhd12q42fe9HD/P8mUeW9u25XjXdeV43/fl+KvRrQEAAEAA/rd32u2M41iOL8tSjh/H8V+f3Ff92R2GoRyf59kTYAkSAAAC8FO7oHS2k3Y76Uwjvc5TS2c+ad7S/HgCLEEABACA7toFpV+y0tlO+jZPvxA9tTQPad7SPHsCLEEABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEBfvdN/Pqd7stLNEW7QaC7/OzrOs8+gJQiAAADQTbug9G2ebgV1m+r1+03ztm2bJ8ASJAAABOCndkHpdtR0B3o680n3ZD21dLaTdjtpnj0BliAAAgBAN/UB+zhVlXCPC3cAAAAASUVORK5CYII=)

```mathematica
In[5]:= CornerFilter[chk, 2]
Out[5]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABT0lEQVR42u3dvQ6CMBiGURAJAwn3f6MMhD9Xh+9lxeB5xg4m9khtamLbaZrOpmhd12q42fe9HD/P8mUeW9u25XjXdeV43/fl+KvRrQEAAEAA/rd32u2M41iOL8tSjh/H8V+f3Ff92R2GoRyf59kTYAkSAAAC8FO7oHS2k3Y76Uwjvc5TS2c+ad7S/HgCLEEABACA7toFpV+y0tlO+jZPvxA9tTQPad7SPHsCLEEABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEBfvdN/Pqd7stLNEW7QaC7/OzrOs8+gJQiAAADQTbug9G2ebgV1m+r1+03ztm2bJ8ASJAAABOCndkHpdtR0B3o680n3ZD21dLaTdjtpnj0BliAAAgBAN/UB+zhVlXCPC3cAAAAASUVORK5CYII=)

```mathematica
In[6]:= ImageDimensions[CornerFilter[chk, 2]]
Out[6]= {16, 16}

In[7]:= ImageType[CornerFilter[chk, 1]]
Out[7]= "Real"
```

### Scope (22)

```mathematica
In[8]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[9]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[10]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[12]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[13]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[14]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[15]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[16]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[17]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];
```

```mathematica
In[18]:= CornerFilter[disk, 1]
Out[18]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABJ0lEQVR42u3dWQrDMAxFUdl0/0u2u4HmwyBq2Tn3M4QMveghPKQtImZgG91PQAABIIAAEEAACCAABBCA//HZdePW2tL5c85S11cBIggEEAACDqdlzYhldR2r11kl675ZXZMKEEEEgAACEJeOBfX+2/EYY8tY0NPxrOdUASIIBBAAAuLSsaCnLqLaGMuuMavV51cBIogAEEAATu+CImlm6pR1SlnvpQJEEAEggADcui6oWrdT7b1UgAgiAAQQAAIIAAEEgAACEGbEzIhBBBEAAghAnLFHbHUV8a4ZtGo7+lWACCIABBCAuHRdULWd8lk7+q0LEkEggAAQoAvy1UQVIIJAAAEg4M1dkH/QUAEiCAQQAAJ0QVABBIAAAkAAASCAABBAAAgoyRdD42S5FRy4rwAAAABJRU5ErkJggg==)

```mathematica
In[19]:= CornerFilter[ramp, 2]
Out[19]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAm0lEQVR42u3RAQ0AAAjDsIN/zyCDhHQS1koy0VltAQAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIwIcWi54BvzXkqDQAAAAASUVORK5CYII=)

```mathematica
In[20]:= CornerFilter[zone, 2]
Out[20]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA7klEQVR42u3dsQqFMAxA0UT8/1+uawdLVQoReu78ph4S6luamdliUmvTn+imzJz+5nBMtQEo7rR2TAAA/WQFWTtr689zdCMyAVbQ5t8KEWHvFH6gmQAryAqygkwAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAACAIwAAQAAACEB4TVXhHTErSN4RMwECsMstaPTQsMed3YKsIBV+iPVjYh2tXTsmwArSq/+CrKNv5/PkhmkCrKC9uwBtwR/O6Zud/wAAAABJRU5ErkJggg==)

```mathematica
In[21]:= CornerFilter[noise, 3]
Out[21]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABq0lEQVR42u3dy46EIBRFUTDt//+x9qQHDpooIg9l7XHlxgo5O3iK0rgsyx7+iDGGM/Z9DzmkZo42J3dm7pzUzCWgKxagMz8lcT5+vkQRNebk6qI2qWuQAAqioNvaKaHXnJTuUnOufP7KHAmgINxW0MzaKZkjARSEx3dBb2cEnR71JQEUREGfifYbNSgBFERB06ogDFBNSwAFUdBpTL6qo5Y/0EsABeH2jViujnqdvZEAWAB1dFBN5yhUAihobuLxdLReyNFECoI6eqquSQIoiIKmraNHqKYlgIIoqOpfRFXTEmABoI7uUk27EaMgZNXR+h91NAWh8S6IdtTRFISOCpqtdm5ZTUsABaFqHa2algAKwoMKsiNqX01LAAVR0O1HeHmghwRYAPhR/nU7Kw/roCBQUBijmpYACqKgV7wE56vVtARYAApCY315cCsFYdqX+EgALIDT0UG/JAEUBDdiA1TTEkBB/qDxmUp55Gra0UQKQpaCck/55r7seLQ5Nd6OceXaJICCJr+3WNf13xxu2/bI49lH3gWVfC+7IAqCLqhBt1NyoyoBFIQzfgGBsB7fd5wO5AAAAABJRU5ErkJggg==)

```mathematica
In[22]:= CornerFilter[rgb, 1]
Out[22]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAm0lEQVR42u3RAQ0AAAjDsIN/zyCDhHQS1koy0VltAQAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIwIcWi54BvzXkqDQAAAAASUVORK5CYII=)

```mathematica
In[23]:= CornerFilter[sky, 2]
Out[23]= -Image-
```

![24x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABACAYAAADlNHIOAAAAcElEQVR42u3RMQEAMAyAMDr/nlsZe4IEMtWmbz0LAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAAVAdTDAF/JaOFigAAAABJRU5ErkJggg==)

```mathematica
In[24]:= CornerFilter[bit, 1]
Out[24]= -Image-
```

![8x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAm0lEQVR42u3RAQ0AAAjDsIN/zyCDhHQS1koy0VltAQAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIwIcWi54BvzXkqDQAAAAASUVORK5CYII=)

```mathematica
In[25]:= CornerFilter[byte, 2]
Out[25]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA5ElEQVR42u3cSQ6DMBBFwW/E/Y9scwGzQGphhnrLLCm10zJRWpIRLWvzCAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAmrd7BLW11qafjzFMgCNIAAAIwKO+tP13dO22YwIcQQIAQADiLsh2ZAIcQQIAQABsQZ/r7A2XCXAECQAAAYg3YqV3JlVbhwkQAAACAEBJbrgLOtt2qn5X8/btyAQAACAAAPS3N2JXt6Pe+ye3IxMAAIAAANCiDu8PFcijAIuAAAAAAElFTkSuQmCC)

```mathematica
In[26]:= CornerFilter[vol, 1]
Out[26]= -Image-
```

![12x10x8 volume](data:image/png;base64,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)

```mathematica
In[27]:= ImageChannels[CornerFilter[rgb, 2]]
Out[27]= 1

In[28]:= ImageDimensions[CornerFilter[vol, 1]]
Out[28]= {12, 10, 8}
```

```mathematica
In[29]:= CornerFilter[chk, 4]
Out[29]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABlElEQVR42u3dQaqEMBBFURXFufvfqCj6N2ANHtjE8M8d1qxzO5WHJmbctu0eHjjPc0jq13U91u/7juqtGMcxqk/T9Fif5zmqTwOaQgABBICA/8tcpZp1XaPVP01H3fxDw7SzLMtjfd93M0ALAgEEgIAuUlC1+lerfEWVgr6Wjqrfm47DcRxRSjQDtCACQAABaJWCqjRSrdppiuj9mU9VT9NONc5mgBZEAAggAK1SULU/561nOOl+m+Fj+4LScUj3R5kBWhABIIAAfC0FDeEbrt5TUJpe0nGQgrQgEEAACCAABBAAAggAAQSAAAIIAAEEgAACQAABIIAAEEAACCAABBAAAggAAQSAAALwA+b0rFb6XZ30BPrXvhf01hm68iYO/0EtiAAQQAA+loLStFPV0xPorShPsodfg0zHwQzQgggAAQSgVQp6656s6lvKvdygkfLWTSJmgBZEAAggAK1S0Ftpp/fbVMuUEj4LqsZNCtKCQAABIKCXFFTdgZ6mnepN09feiP36Bo1qnM0ALYgAEEAAGvEHouKpxz5zHb8AAAAASUVORK5CYII=)

### Options (5)

```mathematica
In[30]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[31]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[32]:= CornerFilter[chk, 2, Method -> "Harris"]
Out[32]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAz0lEQVR42u3dMQ6DMBBE0RnE/a+8qSlSxYgA7x/ABU9rISPhtp0saGbJMrep7ZJ1tujSAAAAIADvbf+3twITIAAABACAAAAQAAACAECn1STjbOf8vn0xNAG2IAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAoBx+Ed3xGEwAAAEAIAAABACAAAAQgLhTPr/fk5XnHqKZAFuQAAAQgDv3AdNWDMhzYT5jAAAAAElFTkSuQmCC)

```mathematica
In[33]:= CornerFilter[chk, 2, Method -> "MinimumEigenvalue"]
Out[33]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABT0lEQVR42u3dvQ6CMBiGURAJAwn3f6MMhD9Xh+9lxeB5xg4m9khtamLbaZrOpmhd12q42fe9HD/P8mUeW9u25XjXdeV43/fl+KvRrQEAAEAA/rd32u2M41iOL8tSjh/H8V+f3Ff92R2GoRyf59kTYAkSAAAC8FO7oHS2k3Y76Uwjvc5TS2c+ad7S/HgCLEEABACA7toFpV+y0tlO+jZPvxA9tTQPad7SPHsCLEEABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEBfvdN/Pqd7stLNEW7QaC7/OzrOs8+gJQiAAADQTbug9G2ebgV1m+r1+03ztm2bJ8ASJAAABOCndkHpdtR0B3o680n3ZD21dLaTdjtpnj0BliAAAgBAN/UB+zhVlXCPC3cAAAAASUVORK5CYII=)

```mathematica
In[34]:= CornerFilter[disk, 1, Method -> "Harris"]
Out[34]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAm0lEQVR42u3RAQ0AAAjDsIN/zyCDhHQS1koy0VltAQAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIwIcWi54BvzXkqDQAAAAASUVORK5CYII=)

### Applications (6)

```mathematica
In[35]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[36]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[37]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[38]:= Binarize[CornerFilter[noise, 2]]
Out[38]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABpUlEQVR42u3dwY6DMAyEYbPq+79y9rbisCg0hGCSb449uK3Q/MQTB7ZSSolGbdv27+ffllytzl4/QY/KBXhYnzuK7q13gXDd6mTT/r9wAARBUPMd/IzFMtTJjEQOgCAI6mL/I0s+VSez9r+fAyAIgpbS0UrmW6xdqaMRgyD6c8xRHD1TI8YB5AK8AkGz2jwzNjkAghZHUESUWg4z64roTKR8Jb/iAAiimj5n7COa5gAIInH0NNG0OBqCqBpHi6bDaCIEkTh6lmiaAyCImhuxbCuZkdF0r/yKAyCIqggSTXMABJE4espo2gENCKLqXNDKjViIoyGIsiIoRNO3RdMcAEF2xEprJmNq2o4YBNHgRuxuS46sA0HkAkSmAxrmgsZH0xwAQRoxjViIoyGIXrIp38vm2eo8hTsOgCAImmYlk3lq2lwQBFEXBGWbdn7j1LRGDILIdPTgaNoqCIKo+0l50TQHQBCZjn7d1LRGDIKoeyMmmm77Lg6AIKOJXe7momkOgCDq+Lwgk9L3vb6QAyCIur/QeaY6zohBENmUjzUe4sEBELS2fgFKZ42pfTizcAAAAABJRU5ErkJggg==)

```mathematica
In[39]:= EdgeDetect[CornerFilter[zone, 2]]
Out[39]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAm0lEQVR42u3RAQ0AAAjDsIN/zyCDhHQS1koy0VltAQAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIwIcWi54BvzXkqDQAAAAASUVORK5CYII=)

```mathematica
In[40]:= ImageDimensions[CornerFilter[Import[Export["/tmp/mathilda_ex.png", rgb]], 2]]
Out[40]= {16, 16}
```

### Properties & Relations (6)

```mathematica
In[41]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[42]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[43]:= ImageDimensions[CornerFilter[chk, 3]] === ImageDimensions[chk]
Out[43]= True

In[44]:= Max[Flatten[ImageData[CornerFilter[chk, 2]]]] <= 1.0
Out[44]= True

In[45]:= Min[Flatten[ImageData[CornerFilter[chk, 2]]]] >= 0.0
Out[45]= True

In[46]:= ImageDimensions[CornerFilter[vol, 2]] === ImageDimensions[vol]
Out[46]= True
```

### Neat Examples (3)

```mathematica
In[47]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];
```

```mathematica
In[48]:= CornerFilter[zone, 4]
Out[48]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABDUlEQVR42u3dwQqDMBQF0fuK///L6arQTdFiYgyeWZcgDG+IWklVVcsOre3+5CdV1WWd0fS6zu91jvAKpkLAZLYR2RmxzmhmXacJkCAJWi4Xq6TsyI7IBEjQs6kkunPxjZ4JkCBI0E1yZAIkSIIkyAQQAAIIAAEEgAACcCHbKv/bMQEgIJ4FwQQQAAIIAAEEgAACQAABIIAAEEAACCAABBAAAggAAQSAAAJAAAEggAAQQAAIIAAEEAACCAABBCDxpbwJIADxpbwJAAEEgAACQAABICCOMoQTNCQIEuQcMUjQg96IxcmqJoAA3GMX5GX9uJ2PCZAg/HUjJkfnU2MCJAhdngXJ0fnsmAAJwoc3kigw1VKYci0AAAAASUVORK5CYII=)

```mathematica
In[49]:= CornerFilter[zone, 1]
Out[49]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAwUlEQVR42u3dsQ3AIAxFQf+I/Vd2ikgoC4CbexVpOdmiS6qqS1dLss+P65gNwHDLFdxfO91tAqwgfZPhFWQCAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAACAAAAQAAACAEAAAAgAAAEAIAAABACAAADQsdb/I8k+d/vVvAkAoOMryNoxAQA01wstowvG6g8UxAAAAABJRU5ErkJggg==)

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
