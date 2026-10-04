### Worked examples

```mathematica
In[1]:= Image3DQ[Image3D[{{{0., 1.}, {1., 0.}}, {{1., 0.}, {0., 1.}}}]]  (* a 2-slice volume is a valid Image3D *)
```

```mathematica
In[1]:= Image3DQ[Image[{{0., 1.}, {1., 0.}}]]  (* a plane is not a volume *)
```

```mathematica
In[1]:= Image3DQ[{{{0., 1.}, {1., 0.}}}]  (* a bare array is not yet an Image3D *)
```

### Notes

`Image3DQ` is the volumetric counterpart of `ImageQ`, and the two are disjoint: `ImageQ`
answers `False` for a volume and `Image3DQ` `False` for a plane.

It is the way to test validity because malformed input to `Image3D` stays **unevaluated**
rather than failing — a bare nested array has not yet been through the `Image3D` constructor,
so it is not in canonical form and `Image3DQ` reports `False`.

A volume is stored `depth × height × width` (slices outermost, indexed `data[[z, y, x]]`),
while `ImageDimensions` reports it fully reversed as `{width, height, depth}`.
