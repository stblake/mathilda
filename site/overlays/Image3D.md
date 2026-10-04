### Worked examples

```mathematica
(* two 2x2 slices; ImageDimensions reports {width, height, depth}, reversed from storage *)
In[1]:= ImageDimensions[Image3D[{{{0., 1.}, {1., 0.}}, {{1., 0.}, {0., 1.}}}]]
```

```mathematica
(* a volume is Image3DQ-valid, not ImageQ-valid *)
In[1]:= Image3DQ[Image3D[{{{0., 1.}}, {{1., 0.}}}]]
```

```mathematica
(* ImageQ is False for a volume *)
In[1]:= ImageQ[Image3D[{{{0., 1.}}, {{1., 0.}}}]]
```

```mathematica
(* integer voxels in {0, 1} infer "Bit"; a stated type is honoured too *)
In[1]:= ImageType[Image3D[{{{0, 1}}, {{1, 0}}}, "Bit"]]
```

### Notes

`Image3D[data]` is a volumetric image, normalising to `Image3D[data, type]`. The
data is a depth × height × width array of voxels (or depth × height × width ×
channels for colour), indexed `data[[z, y, x]]` with **slices outermost**.

`ImageDimensions` reports `{width, height, depth}` — fully reversed from that
storage order, which is Mathematica's convention and the 3-D version of the same
transpose the 2-D accessors carry. Type inference and coercion match `Image`:
all-integer voxels in `{0, 1}` give `"Bit"`, in `0..255` give `"Byte"`, else
`"Real"`, and a stated integer type rounds and clips. `ImageQ` is `False` for a
volume — use `Image3DQ` to test validity — while `ImageDimensions`,
`ImageChannels`, `ImageType` and `ImageData` all accept either rank.
