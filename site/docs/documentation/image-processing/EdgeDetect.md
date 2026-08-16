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
```

```mathematica
In[3]:= EdgeDetect[chk]
Out[3]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA9UlEQVR42u3dQQ6DMAxFwW/E/a+cXgAvqNQahXnLLrrIyCmiQVSSlYvWuvw4VRX13V23w5LNBgAAAAF4b2f3q63/XB2ZAFsQAAEAoKGquxckEwBAAAAIAAABACAAAAQAgAAAEAAAApCtzwU5BW0CAAgAAAGIc0EyAQAEAIAA7NZpCeIZMVuQAAAQgLgXJBMAQAAACECcC7I6JgCAAAAQgLgXJBMAQAAACECcC1KcC7IFCQAAAYh7QTIBAAQAgAAAEAAAAgBAALLn21Sdms5X/3yZAFuQAAAQgKddBXVXNVPvmv/1M2tTz8R1328CbEEABACAhvoANRwtt204DmQAAAAASUVORK5CYII=)

```mathematica
In[4]:= EdgeDetect[disk]
Out[4]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA6ElEQVR42u3cwQ2DMAxAURux/8phAg45WLHh/QWK+uTIEmkzIlboWJevAAAAAQAgAAAEAIAAABAAAAIAQAAAqLb71Aev1etFXGaaAEeQAAAQgD+V1feCdred6m2k2/OYAEcQAAEAoOlbULftYsrzmwBHEAABAKD46BuxbtvO7nNWv7kzAY4gAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAQAAAKAbeC3q7VzP9drQJcAQJAAABCL8R80t5E+AIEgAAAmAL8q+JJsARJAAABMAWJBMAQAAACAAAAQAgAAAEAIAAABAAAAJwvgcg2Ce7XwXFrwAAAABJRU5ErkJggg==)

```mathematica
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
```

```mathematica
In[14]:= EdgeDetect[ramp]
Out[14]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAArklEQVR42u3RsQkAMAhFwW/239kskCZgSHOvFgWvknQG6j6vqaqR+am7r/fctqKvAQAAQAAACAAAAQAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAAC8AAAAAQAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAJppA7+RDb+K9RSmAAAAAElFTkSuQmCC)

```mathematica
In[15]:= EdgeDetect[zone]
Out[15]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAB00lEQVR42u3dwY7DIAxFUTLq//8ys82iUZuCAYfz1g0itd6Vgw0ctdZa3ug4jvJJ50e/+X12tbzv1bN/haZKACbr1YKdWfa/0kgM9sIRB0DQ3jpKKfVX+7dY/i7KIvAYMf+7Y3IABEFQjbRYBL6iM6uR78gBEARBtYeVVvtQmvWhd/c/5AAI2hxBvZajs6BmZKb0zbMcAEGyoPqUD6ssOIIgCKLyriJ2m18B2ImouPWa53mclnmqiEEQpaiIZRz/bkbEARAkC0pXyG6Z53kOLb09vTIiDoAgWVCN7HzeoSjf8r9xAATJgsbxbrEl616ZDAdAEC2FoFl2XmHzBQdAEC2VBWUp1vfKiK7G0ZoIQfQRQbvtfOcAEgABEAASAAEgARAAWuq8IB9lHCAA9IDl6Cwo61XF0xcEQTQdQSv022Sp4nEABEHQI/ttsvQvcYAAFBs0RqFjhc0aK2zKsE0VgujnLChia+dTd+JzAARR6JFlu50X1Cvz4QAIoi5rQRE42qFq5tRECCJnRzs7mgTADRpu0OAACIIg94i5RwyCyNnRZZuzozkAgijFJT69zu1ZeV2IAyDIcnTddf1nhXUhDoAgWVAd1QWdZS1oZNc0B0DQ3voHkCw9whE0ZbMAAAAASUVORK5CYII=)

```mathematica
In[16]:= EdgeDetect[noise]
Out[16]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABz0lEQVR42u3dwW7rMAxEUfqh///L6ua18KIGrIqsLPnMKggSwkEwN+Qoko+IaPFfrX0/7NZxHKHOvTpn/QuaKl/AZH1c2eqO9a5sdfX8rDq9NavrnF/DARD0bh3nLqhCvdbOqtOL1hHcjXxGDoAgXVBKBzLSLVTXGUFWRR0OgCD6sQu6+jVfHVNZQ1xFHQ6AIF3QltZepQ4HQNDLEdS74gM1uXU4AIIMYm13RDw5v+IACHp5F7RKpBzF0fT5vVl14saKHgdAkC6o7YKXFfMrDoAgXZB8ZuIwyAEQZEXMKtjEOhwAQQaxphuZNxhyAAQZxLo2Goimc6+BAyDo5V1QO/lkJ7ysMlRyAAQZxJp8pr6ORXkIovQ4GmpieNc/B0CQDRoWyosGwzt5GgdAkCyoZZ8EuGs0XXENHABBsqCW0XWIpn93zBoHQJAsqPSgVHU4AILIBo14ajTNARAEQRAxcTDkAAgKWdCbMp8nRNMcAEHUdXa0aNrBrRBEE+NoOU/YoAFBZINGbJg7cQAE6YIg4o8HOv+OhiBKOS9INC2OhiAqOr4+6/aFq9epuDk1B0AQpdzKsNeSK9ZxcCsE0aMW5Qd2tUbvbRN3rQNBEERf+gSF3UPUT7bfCwAAAABJRU5ErkJggg==)

```mathematica
In[17]:= EdgeDetect[rgb]
Out[17]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABKElEQVR42u3d2w3DIBBE0SVK/y2TBvKDtXh5nFuAFeVqRyPAuEVEjwR6//+Y1lrJc6p+/yifQCkEEEAACLiXltWCbiOrNZkAEUQACCAARXz9Bc/ajgkQQSCAABCgBb3TIrJ21rLI2ikzASKIABBAAG5bC8pqEaNtZ/Y5HxMggkAAASDgthY0e+2lau3I6WgRBAIIAAFa0JptZLW1HRNAAAggAARoQXbKTIAIAgEEgAAt6P37hXbZKTMBIogAEEAArAXFkmtKJkAEgQACQIAWdCdOR4sgEEAACNCCatvCqV/WMAEiiAAQQAC8KR9H3AVtAkQQCCAABMSh3xHb/ZZCEwACCAABBCAiNtoRMwEggAAQQADiiG/Kj563mf1mvfuCRBAIIAAEXMYPYeNJ4B6nOlEAAAAASUVORK5CYII=)

```mathematica
In[18]:= EdgeDetect[sky]
Out[18]= -Image-
```

![24x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABACAYAAADlNHIOAAAAh0lEQVR42u3ZwQ0AIQwDQQfRf8uhC4LEbAk38iNHJelorOUTAAAgAAAEAIAAABAAALrX7vYryAIACMCnVbwHWAAAAXAHyAIACIA7QBYAQAAACAAAAQAgAAAEAIAAABCAeBOWBQAQgHgTlgUAEAB3gCwAgAC4A2QBAAQAgAAAEAAAAgBAAN7uADQkEXNZbOpXAAAAAElFTkSuQmCC)

```mathematica
In[19]:= EdgeDetect[bit]
Out[19]= -Image-
```

![8x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAtklEQVR42u3dwQ2AQAwDQRvRf8uhAh75gBCzJXh0901nZrKobXTfcs4cJns3AAAACAAAAQAgAAAE4F81yZjBCwAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAQAAACAAAAQAgAAAEAIAAABAAAAIAQAAAABAAAAIAQAAACAAAPdq5vYHupnzclPcFCQAAAQAgAAAEAIAAfLUL8ngOuYkOyg0AAAAASUVORK5CYII=)

```mathematica
In[20]:= EdgeDetect[byte]
Out[20]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAvElEQVR42u3WwQ0AIAgEQTD23zL+LMAPUWdLcIK5jIgKtTU8AQAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAAEAAAAgBAAAAIAAABACAAAAQAgAAA+L3pCc6qKhfgCxIAAAJgBVk1u8x0Ab4gAQAgAFbQh6vGBfiCBACAAFhBveula9W4AF+QAAAQACvo7fXiAgAIAAABuLwFLokNya3JGqAAAAAASUVORK5CYII=)

```mathematica
In[21]:= ImageChannels[EdgeDetect[rgb]]
Out[21]= 1
```

### Applications (4)

```mathematica
In[22]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[23]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];
```

```mathematica
In[24]:= Binarize[EdgeDetect[zone]]
Out[24]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAB00lEQVR42u3dwY7DIAxFUTLq//8ys82iUZuCAYfz1g0itd6Vgw0ctdZa3ug4jvJJ50e/+X12tbzv1bN/haZKACbr1YKdWfa/0kgM9sIRB0DQ3jpKKfVX+7dY/i7KIvAYMf+7Y3IABEFQjbRYBL6iM6uR78gBEARBtYeVVvtQmvWhd/c/5AAI2hxBvZajs6BmZKb0zbMcAEGyoPqUD6ssOIIgCKLyriJ2m18B2ImouPWa53mclnmqiEEQpaiIZRz/bkbEARAkC0pXyG6Z53kOLb09vTIiDoAgWVCN7HzeoSjf8r9xAATJgsbxbrEl616ZDAdAEC2FoFl2XmHzBQdAEC2VBWUp1vfKiK7G0ZoIQfQRQbvtfOcAEgABEAASAAEgARAAWuq8IB9lHCAA9IDl6Cwo61XF0xcEQTQdQSv022Sp4nEABEHQI/ttsvQvcYAAFBs0RqFjhc0aK2zKsE0VgujnLChia+dTd+JzAARR6JFlu50X1Cvz4QAIoi5rQRE42qFq5tRECCJnRzs7mgTADRpu0OAACIIg94i5RwyCyNnRZZuzozkAgijFJT69zu1ZeV2IAyDIcnTddf1nhXUhDoAgWVAd1QWdZS1oZNc0B0DQ3voHkCw9whE0ZbMAAAAASUVORK5CYII=)

```mathematica
In[25]:= Dilation[EdgeDetect[disk], 1]
Out[25]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABAklEQVR42u3d0Q3DIAxAQTvK/iuTCfiDgOHeAK3Uk5GVpG1GRAst6/ERAAAgAAAEAIAAABAAAPqvd/YbtFb7UlNmmgBHkAAAEIATy1F3xKpvO6u2IxPgCAIgAABUZQu6bduZvR2ZAEcQAAEAoDj0jtiq7aLKtmYCAAAQAACyBY2909R7nd22IxMAAIAAANBtW9Dsp46rbEcmwBEEQAAACAAAAQAgAAAUF1wL6l17mX2NyB0xAQAgAAC053NBo7YjT0cLAAABAKAzviN26jf0TQAAAAIAQOFXE/1qohxBAAQAgApuQf5BwwQ4ggQAgADYgmQCAAgAAAEAIAAABACAAGzZB+nULcBa3WRrAAAAAElFTkSuQmCC)

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
```

```mathematica
In[31]:= EdgeDetect[zone]
Out[31]= -Image-
```

![32x32 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAB00lEQVR42u3dwY7DIAxFUTLq//8ys82iUZuCAYfz1g0itd6Vgw0ctdZa3ug4jvJJ50e/+X12tbzv1bN/haZKACbr1YKdWfa/0kgM9sIRB0DQ3jpKKfVX+7dY/i7KIvAYMf+7Y3IABEFQjbRYBL6iM6uR78gBEARBtYeVVvtQmvWhd/c/5AAI2hxBvZajs6BmZKb0zbMcAEGyoPqUD6ssOIIgCKLyriJ2m18B2ImouPWa53mclnmqiEEQpaiIZRz/bkbEARAkC0pXyG6Z53kOLb09vTIiDoAgWVCN7HzeoSjf8r9xAATJgsbxbrEl616ZDAdAEC2FoFl2XmHzBQdAEC2VBWUp1vfKiK7G0ZoIQfQRQbvtfOcAEgABEAASAAEgARAAWuq8IB9lHCAA9IDl6Cwo61XF0xcEQTQdQSv022Sp4nEABEHQI/ttsvQvcYAAFBs0RqFjhc0aK2zKsE0VgujnLChia+dTd+JzAARR6JFlu50X1Cvz4QAIoi5rQRE42qFq5tRECCJnRzs7mgTADRpu0OAACIIg94i5RwyCyNnRZZuzozkAgijFJT69zu1ZeV2IAyDIcnTddf1nhXUhDoAgWVAd1QWdZS1oZNc0B0DQ3voHkCw9whE0ZbMAAAAASUVORK5CYII=)

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [DerivativeFilter](../../image-processing/DerivativeFilter/)

- Source: [`src/imagefilter.c`](https://github.com/stblake/mathilda/blob/main/src/imagefilter.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
