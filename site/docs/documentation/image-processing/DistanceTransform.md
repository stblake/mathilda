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

![12x10x8 volume](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAQQAAAEECAYAAADOCEoKAAAELElEQVR42u3cwW3CQBBA0SHiSCUUQCdbpQtwD1sAd3rgDocoUpASCyVgjz3v3XxDK/M1HlmOAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABgys4RkMH1er19vz4cDu5NQaBqAH4iCoJA8QgIgyAgAMIgCIiAMAgCIiAKgoAACIMgIALCIAgIgDAIAiIgCoKAAAiDICACwiAIiIAwCAICIAqCgAgIgyAgAMIgCIiAMAgCAoAoCIIIIAyCIAAIgyCIAMIgCCKAKAiCACAMgiACCIMgCADCIAgigCgIggAgDIIgAgiDIIgAwiAIAoAoCIIIIAyCIAAIgyCIAMIgCAKAKCSMwk4EQBhWFwQRQBQKB0EA2LpxHB+uW2uL/x/3IgDLRSCbFEHovd+Ox6O7BQEQhE/n8zkiIoQBERAEYUAABGE6DKKACAiCaQEBEARhQAQEQRgQAUGwX0AABMG0gAgIgjAgAIIgDIiAINgvIACCYFpABARBGBABQRAGBEAQ7BcQAUEwLSAAgiAMiIAgCAMCIAjVwyAKIiAImBbidZ8SH4bBh3IFQRhMARFiIAjC4FEAQbBfEAAEwbRQJgwiIAgUDoMAIAjFwyACCELh/YIAIAjFpwURQBCKh0EEEITCYRAABKH4fkEEEITC04IAIAjFwyACCELhMAgAglB8vyACZPLhCJafFkAQAEEABAEQBEAQAEEABAEQBEAQAEEABAEQBEAQAEEAEARAEABBAAQBEARAEABBAAQBEIQZnU4nh4Ag8BgFYUAQEAYEAY8RCAKmBQQBQBAAQQAEARAEQBAAQQAE4Y167w4BQUAUEAQmoiAMCALCQBl7R+AxAkwIgCAAggAIAiAIgCAAggAIAiAILKK15hDCm4qIAIKAACAIiABhhwAgCIAgAIIACAIgCIAgAIIACAIQ3lRkJsMwOARBQAAQBESAsENADEAQAEEABAEQBEAQAEEABAEQBEAQAEEABAEQBEAQAEEABAEQBEAQAEEABAEIH1kFntda2wkCiIAJAUQgr1Q/sPd+q3aTXC4X/xQBEARhEAQREARREAQBEARhEAQREARhEAQBEARhEAQREARREAQBEARhEAQREARhEAQREARhEAQBEARREAQREARhEAQBEARhEAQREARhEAQBEIQQBUEQAcoFIWMYBEEEBEEYBEEABEEYBEEEBEEUigdBAARBGIoHQQQEQRgKB0EABEEYigdBBARBFAoHQQAEgTeHIXsQREAQmDEMGYMgAoLAQmHIEAQBQBCSRGGpIIgAgpAwDHMFQQAQhBWE4Z1BEAEEYWVheGUQBABBWHkU/hsEEUAQNhSGvwRBBBCEjYbhmSAIAIJQJAy/BUEEEISCUfgKggAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEB+d7M5ePBywqciAAAAAElFTkSuQmCC)

```mathematica
In[24]:= DistanceTransform[volb]
Out[24]= -Image-
```

![12x10x8 volume](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAQQAAAEECAYAAADOCEoKAAAEoklEQVR42u3dwXHqMBRAUZlhSSUUQCeukgLoQQV47x68J6ssMsPAJAHrSe+cBr4tPd3R94KUAgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADAM5MlIIi72RQEBMB8CgIiYE4FAQEwr4KACJhbQUAEzK4gIABmWBAQAbMsCAiAmRYERMBcCwICYL4tGCJgzi0UImDeLRACYOYtDiJg9i0KAuAMYDFEAGfBIggAzoMFEAGcCy8e3rZtPwJwOp0sivPhhTNH4BFhcE68aPIIiIKz4iUFQBicGS8nAsLg7HgpARAGZ8jLiIAoiIIXEQBhEAYvIALCIAweXASEgc7P1yQA4xEFZ0wQkkdAGJw1QRAAYRAGQRABYRAGQRAAURAFQRABYRAGQRABYRAFQRAAYRAFQRABUSDAeTxEWIVa610MQgS5bNtmIRI7RHmQZVnKsix2RBgQBGEQBgThRRiIEwYEwW0BtwVBEAaEQRCEAWEQhB4fWhR8X0AQ3BbcFhAEYRAGBEEYhAFB8H3B9wUEwW3BbQFBEAZhQBCEQRgQBN8X8H1BENwWcFvoxDHri39H4Xw+m4IGbrebRRAEYRAABKGTMIiCCAgCbgvlbT/SOl2vV7+NKQjC4BZQihgIgjD4rwCC4PuCACAIbgtpwiACgkDiMAgAgpA8DCKAICT+viAACELy24IIIAjJwyACCELiMAgAgpD8+4IIIAiJbwsCgCAkD4MIIAiJwyAACELy7wsiQPGbivjRVwQBEARAEABBAAQBEARAEABBAAQBEARAEABBAAQBQBAAQQAEARAEQBAAQQAEARAEQBAAQQAEAUgWhMvlYidAEEQBBOFJFIQBBEEYQBCEAQTB9wUQBLcFEARhAEEQBhAE3xdAENwWQBCEAQRBGEAQfF8AQXBbAEEQBhAEYYCspiDPcW/xj9Zam7/4uq5ph2+eZycw2HlMHYQIYcgWBBEQhPBBaBmG0YMgAILQbRBahGHEIIiAIAwVhD2jMEIQBEAQhg/CXmHoNQgiIAgpg/DpMPQSBAEQBEHYIQyRgyACgiAIO0chUhAEQBAEoXEYWgdBBARBEAKFoUUQREAQBCFoGPYIggAIgiB0EoVPBUEEBEEQOgzDu4IgAIIgCAOE4T9BEAFBEITBwvCbIAiAIAjC4FF4FQQREARBSBSGR0EQAUEQhIRqrWVdVwFAEAQB4p3Hgz0ABAEQBEAQAEEABAEQBEAQAEEABAEQBEAQAEEABAEQBEAQAEEABAEQBEAQAEEABAEYwtESQPHnEAQBREAQQAS6eUB/wQkB8MDCgAh4eFFAALyMMODceDFhwDnxosKAs+GlRQFnwSIIA+bfgggDZt7iCANm3GKJAubawgkD5thCCgNm16IKA2bVIosC5tOCCwNmkuyLLwxmEJshDOYOGyMK5gwbJQxmC5smDGYJmygM5gcbKgrmBRssDGYEmy0M5gIbLwzmAIMgCvYeQyEM9hpDIgz2FwMjDPYTQRAFe4ggCIN9QxCEwT4hCPwhDPYGQUgcBXsBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADBfQEX4ydM+vJRcwAAAABJRU5ErkJggg==)

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
