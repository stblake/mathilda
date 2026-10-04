### Worked examples

```mathematica
In[1]:= ImageChannels[SetAlphaChannel[Image[{{0., 1.}}]]]  (* a grey image gains an alpha channel: grey+alpha = 2 *)
```

```mathematica
In[1]:= ImageChannels[SetAlphaChannel[Image[{{{1., 0., 0.}}}]]]  (* an RGB image gains alpha: 4 channels *)
```

```mathematica
In[1]:= ImageData[SetAlphaChannel[Image[{{0., 1.}}], 0.5]]  (* one opacity everywhere, written into the new channel *)
```

```mathematica
In[1]:= ImageData[SetAlphaChannel[Image[{{0., 1.}}], Image[{{0.25, 0.75}}]]]  (* a same-size mask sets per-pixel opacity *)
```

### Notes

`SetAlphaChannel[image]` attaches a fully opaque alpha; `SetAlphaChannel[image, a]` sets one
opacity everywhere (a number in `[0, 1]`); `SetAlphaChannel[image, mask]` sets it per pixel
from a same-size image read as grey — a colour mask is averaged rather than taken as its red
channel alone. A mask of the wrong size is declined, not resampled silently.

The result always carries exactly one alpha channel: an existing alpha is replaced, not
duplicated, so grey → grey+alpha (2 channels) and RGB → RGBA (4).
