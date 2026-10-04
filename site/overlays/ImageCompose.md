### Worked examples

```mathematica
In[1]:= base = Image[{{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}];  (* a 3x3 black field *)
```

```mathematica
In[1]:= ImageDimensions[ImageCompose[base, Image[{{1.}}]]]  (* compositing keeps the base's size *)
```

```mathematica
In[1]:= ImageData[ImageCompose[base, Image[{{1.}}]]]  (* a single white pixel lands at the centre *)
```

```mathematica
In[1]:= ImageChannels[ImageCompose[Image[{{{0., 0., 0.}}}], Image[{{1.}}]]]  (* grey overlay on colour base stays colour *)
```

### Notes

`ImageCompose[base, over]` draws `over` onto `base` with the standard over operator, centred
and keeping the base's dimensions; whatever falls outside is clipped.
`ImageCompose[base, over, {x, y}]` places the overlay's centre at image coordinate `{x, y}` —
x from the left, **y from the bottom**, Mathematica's convention. `ImageCompose[base, {over,
a}]` fades the overlay by a constant opacity `a`.

A grey overlay on a colour base produces colour: grey means the same value in every channel,
so it is replicated rather than zero-padded (which would turn a grey pixel red). The result
carries alpha only if the base did.
