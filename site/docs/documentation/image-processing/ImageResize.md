# ImageResize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageResize[image, {w, h}] resizes to w x h pixels; ImageResize[image, w] gives width w with the height following to preserve the aspect ratio. Resampling -> "Nearest" | "Bilinear" | "Average" selects the method; the default Automatic uses AREA AVERAGING when either axis shrinks and bilinear otherwise. That default is about aliasing: point-sampling a shrinking image destroys every frequency above half the new sampling rate -- a fine checkerboard reduced by nearest-neighbour comes back a flat field -- and no interpolation afterwards can restore what point-sampling discarded. Area averaging is a box prefilter and a resample in one pass, exact for integer reduction factors, using true fractional coverage so a 3 -> 2 reduction is as correct as 4 -> 2. Enlarging has no frequencies to remove, so bilinear is used there; area averaging on an enlargement would degenerate to nearest. Coordinates are centre-aligned, avoiding the half-pixel shift that sx = i * scale introduces at any scale other than 1:1. The result is a "Real" image; sizes must be positive integers.`**

## Examples (35)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (5)

```mathematica
In[1]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[2]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];
```

```mathematica
In[3]:= ImageResize[chk, {8, 8}]
Out[3]= -Image-
```

![8x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAzUlEQVR42u3aMQ6EMBAEwV7E/79sXkDgxARUxyeS0m1gzVSrjdba+nkzk++/d6VPAwAAgAAAEAAAAgBAAP7V7W3n2+/7BzhBAAQAgAAAEAAAAgBAB5vsgrILcoIEAIAAABAAAAIAQACyC8rbjl2QEyQAAAQAgAAAEAAAApBdkLcduyAnSAAACAAAAQAgAAAEILugvB3ZBTlBAgBAAAAIAAABACAA2QX5vn+AEyQAAAQAgAAAEAAAApBdUHZBcoIACAAAAQAgAAAEAIAO9QC292S1XLLeOgAAAABJRU5ErkJggg==)

```mathematica
In[4]:= ImageDimensions[ImageResize[chk, {8, 8}]]
Out[4]= {8, 8}
```

```mathematica
In[5]:= ImageResize[disk, {12, 12}]
Out[5]= -Image-
```

![12x12 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABCElEQVR42u3c0QmDMBRA0ZfSwXQSM1oySdws3SAiYiXtub+BUj28gFWbIqKHHuvlFAAAIAAABACAAAAQAAACAEAAAOj+3k9/gZzzcH3btkufX2sdrpdSTIAtSAAACAAAAfi3Utz8XFDvcz92lFIyAbYgAQAgAAAEINwPiFO/58/e0fFdvZ9gAmxBAAQAgAAAUEx4P6C1NlxflmXqE7Tv+3B9XVcTYAsSAAACAEAAAAgAAAEAIAAABCB+/bmgo/dwZ78fcHR8JsAWJAAABACAAIT3hMN7wibAFiQAAAQAgAC4Dgj/G2oCbEECAEAAAAiA6wCZAAACAEAAAAgAAAACAEAAAAgAAH27D6WqJ2Fooz5uAAAAAElFTkSuQmCC)

### Scope (19)

```mathematica
In[6]:= chk = Image[{{0.,1.,0.,1.},{1.,0.,1.,0.},{0.,1.,0.,1.},{1.,0.,1.,0.}}];
```

Area averaging

```mathematica
In[7]:= ImageData[ImageResize[chk, {2, 2}]]
Out[7]= {{0.5, 0.5}, {0.5, 0.5}}
```

```mathematica
In[8]:= ImageData[ImageResize[chk, {2, 2}, Resampling -> "Nearest"]]
Out[8]= {{0.0, 0.0}, {0.0, 0.0}}

In[9]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[10]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[11]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[12]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[13]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[14]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[15]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];
```

```mathematica
In[16]:= ImageResize[rgb, {8, 8}]
Out[16]= -Image-
```

![8x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAwElEQVR42u3aMQ6AIBREwcVwCDo4qjeXM1D5jfNqCsNkG2Ib435y0kpKnZ/Fvufw/BW9GgAAAAQAgAAAEAAAAvCv+nIHFgBAAAAIAAABACAAAAQg3oJkAQAEAIAAABAAAAIAQADiLUgWAEAAAAgAAAEAIAAABOCb9ekOLACAAAAQAAACAEAAAAhA/BckCwAgAAAEAIAAABAAAAIQb0GyAAACAEAAAAgAAAEAIADxFiQLACAAAAQAgAAAEAAAAlC6DcsWBDuYXpm0AAAAAElFTkSuQmCC)

```mathematica
In[17]:= ImageResize[sky, {12, 8}]
Out[17]= -Image-
```

![12x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABACAYAAADlNHIOAAAAl0lEQVR42u3awQ2AIAxAUTHecHDW4uYqcMMFZAtMyvsj8NJDQ9NTy3fot05PAACAAAAQAABa39Xf7BVMAAABACAAAARgsz2g2QNMAAABACAAAARgv/+AYQ8wAQAEAIAAABAAd0EyAQAEAIAAABCA+HdBt1cwAQAEAIAAABAAd0EyAQAEAIAAABAAd0EyAQAEAIAAABCAWE0vshDsUdqb1QAAAABJRU5ErkJggg==)

```mathematica
In[18]:= ImageResize[bit, {4, 4}]
Out[18]= -Image-
```

![4x4 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAnElEQVR42u3RAQ0AMAgEsTHl7xxkkJCehGsl6ae1vgUAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAAAgAAAEAIAAXGiVVAz+1HLpUAAAAAElFTkSuQmCC)

```mathematica
In[19]:= ImageResize[byte, {8, 8}]
Out[19]= -Image-
```

![8x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA2UlEQVR42u3auw3EIBRE0bHlmAIIKJmWcQ1I/ghxbky0R37BaI/W2shEpZR4/9z7M/o1AAAACAAAAQAgAAAEYK+uWmtsNbEFOUECAEAAAAgAAAEAIACxBcW2YwtyggQAgAAAEAAAAgBAAGIL8j62ICdIAAAIAAABACAAAAQgtiDvfQFOkAAAEAAAAgBAAAAIQGxBW7zvvfsCnCABACAAAAQAgAAAEIDYgrLkVjP73v+CnCABACAAAAQAgAAAEIDYgrLqVvP27+MLcIIACAAAAQAgAAAEAIA+7AaMmxUdJ8bm2QAAAABJRU5ErkJggg==)

```mathematica
In[20]:= ImageResize[zone, {16, 16}]
Out[20]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAEWklEQVR42u2dzSv8URjF70xDhIRJY+OlWJFSrCwmFCEv2bGwULPQzMbSQqxsLOzUbJSFBUs7KSnJS3nJRqwQkoSYmYzB+AfOWTyr+f76nbP8ZMa3ztzb6bnPc7++gYGBrAPa29tD2I2NjUEeiUQgr6+vh/z9/R3yh4cH09/n5eVBXllZCXkoFIL86+sL8s3NTcgXFxchv7u7g3x+fh5yv5NyKhkgA2SAJAP+X/nKy8thCpqamoIfmJychPzz8xPyra0tyLe3tyG/urqC/OPjA/JAIAB5VVUV5K2trZD39vZC3tLSAvnp6Snk09PTkF9fX2sFaAuSZIAMkGSAp1JQLBaDKWhubg5+4OnpCfJ4PA75xsYG5Pf39/iBfD7I8/PzIf/9/XWW2k5hYaEpHU1MTEA+NDQE+f7+PuTRaFQrQFuQJANkgCQDvKQAO8litR2WdlZXVyFPJpOQNzQ0OMsJWllZmSntsJOpy8tLyHd3d03PX1RUBHl3dzfk4+PjWgHagiQZIAMkGeCpFMRSx/r6urPUdlhaaGtrg7yjowPyxsZGUwrKZDKQ397eQn5wcOAsJ3RnZ2eQr62tQd7c3Az54OCgVoC2IEkGyABJBngqBbGuY5YK2EkWq+2wtNPZ2Ql5XV2ds9Refn5+IK+pqYG8uLgY8kQi4Szd0YeHh5AfHx9D3tPToxWgLUiSATJAkgGeSkFsJot1KbO+HVZTYrUdlnbYbBebBWN9QayPKJVKmZ7z5OTEWfqjLi4uIO/v79cK0BYkyQAZIMmAf6IWxGayWLpgJ1aMs9oOSzt+v9/ECwoKIC8tLYW8oqIC8pKSEmepib28vDhLzUorQFuQDJBkgAyQcpWCWOpgE+jWmSzWt8NSAft+lnay2azpe9jzpNNpyL+/v03Pw1Iiq6FpBWgLkgGSDJABUq5SEDuBYvftsBMfNpPFupRZ3w5LEay2w9LO29ubs9zKeHNzA/nz87Oz1LKqq6tttSz9BrUFyQBJBsgAKUcpiN2lzO7PYffhsAl0NpPFupRZ3w47yWK1HZZ2jo6OnKWr+fX11VluU2ScnTxqBWgLkgGSDJABUq5SEDvJYncpsxTE7tths2ZsJot1KbO+HXaSxWo7LO2cn587S19TX18f5E1NTZDv7OxoBWgLkmSADJBkgKdSEJsEHx4edpa7lNl9Qey+HfZ/2UwW61JmfTvsJIvVdljaGRkZgXx0dNSU7lZWVrQCtAVJMkAGSDLAUymIvRWU9e2wN0ewPhl2uyC7b4dNoLOZLNZvw56HnVix2g5LO8Fg0FnesspqYloB2oJkgCQDZICUI/lCoVDWcgvizMwM5O3t7ZA/Pj46y8kU675mE+ism5p1KbMUxE6yWG2HvUlkaWkJ8nA4rBWgLUiSATJAkgGeSkHLy8swBc3OzjrLrBZ7TxZ7c0Rtba2zTOizyXo6gU5qRKxLmdWm2EkWq+10dXVBvrCwoBWgLUiSATJAkgFe0h/R0YzATk+5GAAAAABJRU5ErkJggg==)

```mathematica
In[21]:= ImageResize[noise, {16, 24}]
Out[21]= -Image-
```

![16x24 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEAAAABgCAYAAACtxXToAAACKElEQVR42u2cMY4CMQxFh9VeYeAUlHAKKAFxBkTHVFBQINEiJKDnDoiTAKeghdLbMi4SRRlmJ/JLZ7C+0Jf+d+LEtA6Hg2Qf6/F4fIbZ6/Uqxf1+vxTf7/eg/Kbh/2TGFwRYJ+BXa6zdbjvj2+0WlN90fCQAAdY9oNfrOeum1oyuozpfa9KX78M/nU7O73V8PB6D8pEABBhfrW63WzoLjMdjZx3VGvPlaw03DR8JQIB1Dzifz+Kq4zrWGtPf6/O4L9+Hr/FCf48vHwlAgHUPmEwm4qqbWjO67vrydd1tGj4SgADrHrDZbKTK83XV5/dv3yMgAQiw7gHT6VRC+upak1X37evGRwIQkHEvUGufv2n4SAACrO8DlsulxOy13+93KdaeEruX1/i6rms8HfvykQAEWPeA3W4nKe3dq76nQAIQYP0scLlcSh90Oh1nXff17XW+7vnF4hdFEZS/WCy4F0ACEODYB+z3e/nPNzs6fzQa1YqPBCDAugdst1tp0psdH37VZxMkAAHWPWC1WomrBxfaY/Plx+L7+vwan3sBJAAB7n7A8/lMavZ3OBwG4Q8GA94HIAEIcHhAnudJz/7G4iMBCMAD8qTfBMW+TUYCEGC9H7BeryXl2d9YfCQAAdY9oCgK+ebsb9Nni5EABFj3gNlsJinP/sbiIwEIsO4B1+tV6vzvrqbNFiMBCLDuAfP5XFKe/eVeAAlAQFbpvUBqs7+x+EgAAoyvP3tdQn29t4yzAAAAAElFTkSuQmCC)

```mathematica
In[22]:= ImageResize[ramp, {8, 16}]
Out[22]= -Image-
```

![8x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADAAAABgCAYAAABbjPFwAAAAcElEQVR42u3PoREAIAwEwYD7NtJ/j9AAEpOZPXluV5JTj7q7JvxdwwMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAD43wX+VwSkuSehUgAAAABJRU5ErkJggg==)

```mathematica
In[23]:= ImageChannels[ImageResize[rgb, {8, 8}]]
Out[23]= 3

In[24]:= ImageDimensions[ImageResize[zone, {20, 10}]]
Out[24]= {20, 10}
```

### Applications (4)

```mathematica
In[25]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[26]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];
```

```mathematica
In[27]:= Binarize[ImageResize[zone, {16, 16}]]
Out[27]= -Image-
```

![16x16 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAABRklEQVR42u3d0Q3DIAxFUVN1/5XpBP0gMjIm5w7QNr2y9WQIjIiYkcCcKR/ThjFGyv/zCZRCAAEEgID3MlZTUFbaWU0Ru6l6LhWgBREAAgjAaSloNRWclmq6pCMVoAURAAIIQBFfaefZ82alIxWgBREAAghAVQqKS2cvu1PcajqyL0gLAgEEgIAuKWj3zKdqN/W/762acakALYgAEEAAbp0FeXfMLEgLAgEEgAArYtFwF/TuGZEK0IIIAAEEgAACQAABIIAAhFlQnH4KYiStAKoALYgAEEAAbk1Bu980735aowoggAAQQABOS0G798lk3UDR/VwjFUAAASCAABQxZtJQ5m2nKWalMhWgBREAAghA9xs0bk1Hu2dQKkALIgAEEIC33aZ62i7lqt+pArQgAkAAAThtRaxLGun+7pgK0IIIAAEEoIgfMkFevbzKQAQAAAAASUVORK5CYII=)

```mathematica
In[28]:= EdgeDetect[ImageResize[disk, {12, 12}]]
Out[28]= -Image-
```

![12x12 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAA2ElEQVR42u3ciQnDMBBFwb/G/bcsFyEWHZ7XQEyGFTJKVElGtKzHVwAAgAAAEAAAAgBAAAAIAAABAKD+3u4PGOPs44aqMgGWIAEAIAAABOC2avZ3QbP7/O599u7PbwIsQQAEAIAAAFAuPA9Yvc+ffb7u8wwTYAkCIAAABACAAAAQAAACAEAAAAgAAAEAIAAABACAAAAQAAACAEAAACgX3Bd0+v+ETYAlSAAACAAAAcgP7wtyb6gJsAQJAAABACAA3gNkAgAIAAABACAAAAQAgAAAEAAAArBdH2b9Gb8UFg+WAAAAAElFTkSuQmCC)

### Properties & Relations (5)

```mathematica
In[29]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[30]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[31]:= ImageDimensions[ImageResize[chk, {8, 8}]] === {8, 8}
Out[31]= True

In[32]:= ImageChannels[ImageResize[rgb, {4, 4}]] === ImageChannels[rgb]
Out[32]= True

In[33]:= ImageDimensions[ImageResize[chk, {16, 16}]] === ImageDimensions[chk]
Out[33]= True
```

### Neat Examples (2)

```mathematica
In[34]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];
```

```mathematica
In[35]:= ImageResize[zone, {24, 8}]
Out[35]= -Image-
```

![24x8 result](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGAAAAAgCAYAAADtwH1UAAAC6klEQVR42u2azUsyURTGr1+llBYqYkRF0IeKCeJntIiiVuK2hctW9b/0V0TQRmgT2koEF0Fa4EI0ClqJYg1EU5mU+m57zuK4THjP2f3mzsCdOdzn3PPMNcTj8aH6Fe/v779RBYNB4FQqxY7PzMwAv729AbfbbeBWq8WO0/l8f38Dm0wm4KmpKWCn0wm8sLAAPDc3Bzw7Owvc6/WA7+/vgUulEnCxWFRcbG9vAxuVxJ+GJEAS8H+H4eTkZMhpPNX0m5sb4EKhwI4/PT0B//z8AFutVmCHwwFss9mALRYLcL/fB+52u2wN+fj4wA9gMAC73W7gZDIJvLOzA7y5uQk8MTHBfp/T01NZASJBEpKAsakBuq5DDahWq3DD+fk58MXFBfDr6yuwx+MBjsViwBsbG8ArKyvsPt1ut7M1gNYUqvEvLy/AjUYDuFarKe79Hx4e2JqUSCSAj46OgLe2ttg+R1aASJAkQOIPw9xsNuHC5eUlcD6fB35+fmY1nGre/v4+cCAQUJwXQzXfbDaz+/bhcMj2BV9fX8Dr6+vAfr8feHFxETiXywHX63XFeUFerxfY5XIBh8NhWQEiQRKSgLGpAZOTk3Bhenqa3ffSoBpLvZdRXgx9nnpDVPNH1QDaF1A///Pzk52fruuK85YGg4HivJ9R34/OV1aASJAkQOIPw6BpGoiSpmlsX3B2dgb8+PiouH+0a2trwD6fj+0jlpaWFPd/gGou1Xyq4Z1OR3FeEGX6PtTrol7V3t4e8OHhIfDy8jJwpVKRFSASJCEJGJsakMlkoAYcHBzADZFIRHHnfKimXV9fA9/e3rKaTPfZdJ9Pa4rRaGT35dQLouO076FezerqKvDu7i5wPB4Hnp+fV9y5oWw2C3x1dSUrQCRIQhIwNjXg+Ph4yJ3roZoajUZZTaT+Oj1nQzWZejO0xtAaQff9tEZQL4meFaV9BX2ezof6/+VyGfju7o7tG0KhEHA6nZYVIBIkIQkYl/gHtSQpq9iDjaIAAAAASUVORK5CYII=)

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [NDArray](../../linear-algebra/NDArray/), [ImageData](../../image-processing/ImageData/), [List](../../other-advanced/List/), [ImageConvolve](../../image-processing/ImageConvolve/), [GaussianFilter](../../image-processing/GaussianFilter/), [Table](../../lists-and-iteration/Table/)

- Source: [`src/imagegeom.c`](https://github.com/stblake/mathilda/blob/main/src/imagegeom.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)
