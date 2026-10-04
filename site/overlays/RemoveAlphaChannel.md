### Worked examples

```mathematica
(* an RGBA pixel loses its alpha, leaving 3 colour channels *)
In[1]:= ImageChannels[RemoveAlphaChannel[Image[{{{1., 0, 0, 0.5}}}]]]
```

```mathematica
(* compositing a half-transparent white pixel over black (b=0) resolves it to grey *)
In[1]:= ImageData[RemoveAlphaChannel[Image[{{{1., 1., 1., 0.5}}}], 0]]
```

```mathematica
(* without a background, the colour values pass through unchanged *)
In[1]:= ImageData[RemoveAlphaChannel[Image[{{{1., 0., 0., 0.5}}}]]]
```

### Notes

`RemoveAlphaChannel[image]` drops the alpha channel.
`RemoveAlphaChannel[image, b]` instead **composites** over a background of
brightness `b`, which is the difference between forgetting the transparency and
resolving it: a half-transparent white pixel over black is grey, where dropping
alpha would leave it white.

The composite writes `v · a + b · (1 − a)` per colour channel. The output keeps
only the colour channels — a 4-channel RGBA image becomes 3-channel RGB, a
2-channel grey-plus-alpha becomes 1-channel grey.
