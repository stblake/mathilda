# DistanceTransform

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DistanceTransform[image] replaces each pixel by its EXACT Euclidean distance to the nearest background pixel; background pixels are 0, so the value rises toward the interior of a blob. DistanceTransform[image, t] takes pixels above t as foreground (default 0). Exact rather than the classic two-pass chamfer approximation, which cannot represent sqrt(2) with integer steps and so gets diagonal distances a few percent wrong -- invisible on a picture and fatal to a test. Uses Felzenszwalb and Huttenlocher's lower-envelope-of-parabolas method, O(n) per row with no sorting. Separability is EXACT here because squared Euclidean distance is a sum over the axes, so minimising it decomposes per axis; the square root is taken once at the end rather than per pass.`**

## Examples (33)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (6)

```mathematica
In[1]:= ImageData[DistanceTransform[Image[{{1., 1., 1., 1., 1.}, {1., 1., 1., 1., 1.}, {1., 1., 0., 1., 1.}, {1., 1., 1., 1., 1.}, {1., 1., 1., 1., 1.}}]]]
Out[1]= {{2.82843, 2.23607, 2.0, 2.23607, 2.82843}, {2.23607, 1.41421, 1.0, 1.41421, 2.23607}, {2.0, 1.0, 0.0, 1.0, 2.0}, {2.23607, 1.41421, 1.0, 1.41421, 2.23607}, {2.82843, 2.23607, 2.0, 2.23607, 2.82843}}

In[2]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[3]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[4]:= DistanceTransform[chk]
Out[4]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzUlEQVR42u3aMQ6EMBAEwV7E/79sXkDgxARUxyeS0m1gzVSrjdba+nkzk++/d6VPAwAAgAAAEAAAAgBAAP7V7W3n2+/7BzhBAAQAgAAAEAAAAgBAB5vsgrILcoIEAIAAABAAAAIAQACyC8rbjl2QEyQAAAQAgAAAEAAAApBdkLcduyAnSAAACAAAAQAgAAAEILugvB3ZBTlBAgBAAAAIAAABACAA2QX5vn+AEyQAAAQAgAAAEAAAApBdUHZBcoIACAAAAQAgAAAEAIAO9QC292S1XLLeOgAAAABJRU5ErkJggg==)

```mathematica
In[5]:= DistanceTransform[disk]
Out[5]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzklEQVR42u3ciw2AIBBEwT1j/y1jFSuQzKvAODlC8DNJVrStxy0AAEAAAAgAAAEAIAAABACAAAAQAADq9p52QWt1H9DNjAkQAAACAEDJtN8Lau9q6jeovGsyAZYgAAIAQAAACAAAAQCgXHgWdPuZz64zIhNgCQIgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIDiGzHfiAkAAAEAIAAABACAAABQzj0Lir8mmgBLkAAAEAC7IJkAAAIAQAAACAAAAQAAQAAACAAAAQCgv/sAoIQVtzrUII4AAAAASUVORK5CYII=)

```mathematica
In[6]:= ImageDimensions[DistanceTransform[chk]]
Out[6]= {16, 16}
```

### Scope (20)

```mathematica
In[7]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[8]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[9]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[10]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[11]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[12]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[13]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[14]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[15]:= volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];
```

```mathematica
In[16]:= DistanceTransform[ramp]
Out[16]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAoklEQVR42u3RsREAIAwDMcP+O5sFUlLl9KVL6yRphtpx1ueuCwAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAADb0AMYNBLyt8NpwAAAAAElFTkSuQmCC)

```mathematica
In[17]:= DistanceTransform[zone]
Out[17]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAm0lEQVR42u3RAQ0AAAjDMMC/5yODhHQS1k6S0lljAQAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIwIcWDdQEvLmeTUQAAAAASUVORK5CYII=)

```mathematica
In[18]:= DistanceTransform[noise]
Out[18]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABBklEQVR42u3cSw7CQBADUTfi/ldudtwgmU+/WrGFkq3EI6a6u4NlfPwEBBAAAggAAQSAgHF8p33hqvp/3uEdVAIImE0l6Z0iKQEgYNRTkNqRAAJAAAEggAAQQADM0Zk0TUuAClJB254WPY0TMRCwVQWZpiWAABBAAAggAAQQAHN0bp+mJUAFqaAj/siQS6dpCVBBKmjb0yIJAAEEgAACQAABIMAcHdO0BKggvF5BpmkJUEFYWEGmaQkgAAQQAAIIAAEEwBydW6ZpCVBBKujo699z+DQtASpIBR13x44EgAACQAABIIAAEGCOjmlaAlQQXq8g07QEqCAsrCDTtAQQAAIIAAHH8gNdukq4IViudgAAAABJRU5ErkJggg==)

```mathematica
In[19]:= DistanceTransform[rgb]
Out[19]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAm0lEQVR42u3RAQ0AAAjDMMC/5yODhHQS1k6S0lljAQAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIwIcWDdQEvLmeTUQAAAAASUVORK5CYII=)

```mathematica
In[20]:= DistanceTransform[sky]
Out[20]= -Image-
```

![24x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABACAYAAADlNHIOAAAAcElEQVR42u3RMQEAAAjDsIF/z0MGTyqhmbaN3loLAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAAlBz9IwR8WwehAwAAAABJRU5ErkJggg==)

```mathematica
In[21]:= DistanceTransform[bit]
Out[21]= -Image-
```

![8x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzUlEQVR42u3asQ2AMAwAwTdi/5XNDGmg4L5GaU64sDy7ux00Myefd/j8796/0qcBAABAAAAIAAABACAA/2qqtdv57n1/gBEEQAAACAAAAQAgAAD0YrfdTu6CjCABACAAAAQAgAAAEIDcBWW3k7sgI0gAAAgAAAEAIAAABCB3QdntuAsyggQAgAAAEAAAAgBAAHIXZHfkLsgIEgAAAgBAAAAIAAAByF2Q9/0BRpAAABAAAAIAQAAACEDugnIXJCMIgAAAEAAAAgBAAADopR6fA2WxoEgxUAAAAABJRU5ErkJggg==)

```mathematica
In[22]:= DistanceTransform[byte]
Out[22]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAm0lEQVR42u3RAQ0AAAjDMMC/5yODhHQS1k6S0lljAQAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIwIcWDdQEvLmeTUQAAAAASUVORK5CYII=)

```mathematica
In[23]:= DistanceTransform[vol]
Out[23]= -Image-
```

![12x10x8 volume](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAQQAAAEECAYAAADOCEoKAAAFmUlEQVR42u3dzVHcMBiAYYnhSCXuZ6ukAHrYArjTA+coF5IBBjIB1uvv53lmOHDE1r7IluwdAwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAf5kOARE8Pz+v17/f3d0Zm4JA1wB8RBQEgeYREAZBQACEQRAQgS9bL2G4caQFgcYR+GC2IAqCQNcA/OMyQhgEgc4REAZBQACEQRAQga+PY6sRgkC8AFw1ApYpBYGYEdgrCtP+BUEgZwQuFYZpY5MgUCMA34nC7uNRFARBBGKH4ZAxKAyCIABiIAyCUCoCv5KevxV9HAqDIGQMQKbzuLKNRVEQhMwRiHo+V/Yx2TkMghDvXsAqcE5XhbHZMQyCEPeGYLcwhB2LncIgCLFXBVaRc7yyj78uURCEHMuClcIwMy+XVn+aUhB6PSw0A8wOSiyVVg2DIOTcHBQ9DKvJUmm5MAhC3h2Ckdb5V9IxuLyURRCqbRFeAc5vpnsc3tYkCC2eE4hww24lvoS52Gco84qEINR6WCjT/oOsEZiVlyqnAKQPQKa9B9PbmgRBBOpOxasev9l1Y9MUgHIREIbAN2Cjh2E2eVpwNgtAlfO812XHOvq4RA3D9M4AEUh6TNKP3YhRmN4ZIADJwygMHQfRu/sBNsLUuuRb3Y9DlCjMpu8MmAJgNeVoDw8Pb34/nU6H/x23yVcFpo0wo9qzNavyMXkfgWhCBOF8Pq9t2y7xoTnyFd+z2LsLjp61roqzgOjCzBAeHx/HGGNcKAwZ/ssJgFmAIFwpDFH/y4lAQZkDED4Ir8NwQBQuvRFmFnh1GMUjkCIIV5wtrLH/rjYxEABBCBqGNZJud0UEBGGfG48/WQkQBhEQhGL3F6bHihGAIkHYYZlyCYAIkDgIV7iMEAEBEITGYRAAERCEan/QQfsXEABBaL5/AREQBGFABARBGBAAQXB/AREQBLMFBEAQhAEREITWhEEAEAT3F0QAQTBbGDu/Svz+/n45EoIgDGYBQwwEQRhcCiAI7i8IAIJgttAmDCIgCDQOgwAgCM3DIAIIwuXDsLZtmwKAIPA3Ci+zhSkCCAKhwiACCML40fcszMxhEAAEYZ+vc5tZ7i+IAIKQMAyXWo0QAASh1mXEl8MgAghC/cuIT8MgAAhC4zBs2yYChHLjEHwpDGuP2QIIQv4ZAwgC+8wWQBCsRoAgFDSFAUHAV8ojCD7EZgsM+xAaRGGZLWCGgA8xCIIpPwiCMIAgCAMIgvsLIAhmCyAIwgCC4DICBAEQBEAQAEEABAEQBEAQAEEABAEQBEAQAEEABAEQBEAQAEEABIHh6+MRBHx9PIKA2QKCgNkCgoAwIAgAggAIAiAIgCAAggAIAiAIgCAAggAIAiAIgCAAggAIAiAIgCAAggAIAiAIgCDw3nz5AUHgTRhAEDBbQBAQBgQBlxEIAiAIgCAAggAIAiAIgCAAggAIAoAgAIIACAIgCIAgAIIACAIgCIAgAIIACAIgCIAgAEe6dQjgWKfTaQoCiIAZAojA8M1B/+N8Pq9ug+Tp6cknRQAEQRgEQQQEQRQEQQAEQRgEQQQEQRgEQQAEQRgEQQQEQRQEQQAEQRgEQQQEQRgEQQQEQRgEQQAEQRQEQQQEQRgEQQAEQRgEQQQEQRgEQQAEYYiCIIgA7YIQMQyCIAKCIAyCIACCIAyCIAKCIArNgyAAgiAMzYMgAoIgDI2DIACCIAzNgyACgiAKjYMgAILAzmGIHgQREASuGIaIQRABQeCgMEQIggAgCEGicFQQRABBCBiGawVBABCEBGHYMwgigCAkC8MlgyAACELyKPw0CCKAIBQKw3eCIAIIQtEw/E8QBABBaBKGz4IgAghCwyj8CYIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAfL8BzYvCXivSMBEAAAAASUVORK5CYII=)

```mathematica
In[24]:= DistanceTransform[volb]
Out[24]= -Image-
```

![12x10x8 volume](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAQQAAAEECAYAAADOCEoKAAAEMklEQVR42u3cMW7iUBCA4ecopbtdlDtwAN+EU/oAvoMPQO9yRUrqZZuwIogkMtjY7833SZGSJkpI+DUzSKQEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAISoPwbKOx+Ph8uu6rjceFQQhcARuEQYEIXgERAFBKDcAp4sv3x/5XsKAIOQfgWvvj35/YUAQ8g2AMCAIIvCcMIgCgpB3AIQBQRABawSCIALCgCAIgDUCQRABYUAQBMAagSCIgDAgCAJgjUAQREAYEAQRsEYgCAIgDAiCCAiDKJDSyxp+iL7vT2Iwid8fH3e/ics9b+SCCWHyIJw/3263/irWCAThM2EQBoKuDLfs93t/nRWsEfe+HyQmhEknBNOCwyOCkITBGoEgjHqFQRiEgWA3hOS+kNwXMCGMZFpwX0AQhMEagSAIgzDghpDcF5L7AiYE04L7AoIgDNYIBGFuwiAMBLshJPeF5L6ACcG04L7Aw16j/uLnaUEYJpsYRkWh67pfHjZBEIayo/DltCAAgpBdGERh2jCIgCCYFgKr67pq2/Z8DxIDQRCGaLqu+//5RQwQBGGIGAEEwX1BABAEUYg0LYiAIBA4DAKAIAQPgwggCIHvCwKAIASfFkQAQQgeBhFAEAKHQQAQhOD3BRFAEAJPCwKAIAQPgwggCIHDIAAIQvD7ggiQvMkq3vQVQfja6eMDEIRPYQAEwbQAgiAMIAjCAILgvgCCYFoAQRAGEARhAEFwXwBBMC2AIAgDCMJMDsIAgnAdhYP7AgjC3GEwLUDmNwRrBAjC09YIYYBMX2VwXwBBcF8AQXBfAEFwXwBBcF8AQXBfAEFwXwBB+E7TNG/uCyAIc0bBfQFyXhmapnnLKAymBQShgDBYIyDHo+JMYXBfgJxfZXBfAEFwXwBBcF+Apbzm+oOfo9D3/Z8ZorCZaY2oAvxPCaAJwX0h8JPldPWBCaHIaWGuiaGEacET34TgvhD4vmAKMCGYGJ50XzAFYEJwX1j9JIAJwbQQbGLwxCfOhJDxfcEUgAnBfcEUwLKqyL/8DGtEGhuGYRiqqBHY7XYvnoKCIAzTByGbKUAEBCF6GDYzBcEUgCCUFoaRQRABBKHkNeKHIAgAghApDDeCIAIIQtQ1YhiGJAAIgjCcV4aNCCAIorC6IAgAgrBwGJYOggggCCsKwxJBEAEEYaVheEYQBABByCQKcwVBBBCEDMMwVRAEAEEoIAyPBEEEEITCwjAmCAKAIBQehZ+CIAIIQqAw3AqCCEDQMLRt+9cjAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAFCIf79uUq5nWIhqAAAAAElFTkSuQmCC)

```mathematica
In[25]:= ImageChannels[DistanceTransform[rgb]]
Out[25]= 1

In[26]:= ImageDimensions[DistanceTransform[vol]]
Out[26]= {12, 10, 8}
```

### Applications (2)

```mathematica
In[27]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[28]:= Dilation[DistanceTransform[disk], 1]
Out[28]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzklEQVR42u3dAQqEMAxFwR/x/leOl9A20nkXEBxSShbcStLRti6vAAAAAQAgAAAEAIAAANC67l0P7p61gqoqE+AIEgAAAnBS9dYvYtNuNX+5NZkARxAAAQAgAAAEAIAAAFAG74JO2/l8vSMyAY4gAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAIrvBflekAAAEAAAAgBAAAAIAABl3S4o/kHDBDiCBACAALgFyQQAEAAAAgBAAAAIAAABmNgDNvgVt7gWgnMAAAAASUVORK5CYII=)

### Properties & Relations (3)

```mathematica
In[29]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[30]:= ImageDimensions[DistanceTransform[chk]] === ImageDimensions[chk]
Out[30]= True

In[31]:= Min[Flatten[ImageData[DistanceTransform[chk]]]] >= 0.0
Out[31]= True
```

### Neat Examples (2)

```mathematica
In[32]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];
```

```mathematica
In[33]:= DistanceTransform[zone]
Out[33]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAm0lEQVR42u3RAQ0AAAjDMMC/5yODhHQS1k6S0lljAQAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIwIcWDdQEvLmeTUQAAAAASUVORK5CYII=)

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Image3D](../../image-processing/Image3D/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
