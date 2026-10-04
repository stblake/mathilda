# Image3DQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Image3DQ[expr] gives True if expr is a valid volumetric image in canonical form. Malformed input to Image3D stays unevaluated, so this is how validity is tested.`**

## Examples (37)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (10)

```mathematica
In[1]:= Image3DQ[Image3D[{{{0., 1.}, {1., 0.}}}]]
Out[1]= True

In[2]:= Image3DQ[Image[{{0., 1.}}]]
Out[2]= False

In[3]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[5]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[6]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[7]:= Image3DQ[chk]
Out[7]= False

In[8]:= Image3DQ[rgb]
Out[8]= False

In[9]:= Image3DQ[bit]
Out[9]= False

In[10]:= Image3DQ[byte]
Out[10]= False
```

### Scope (16)

```mathematica
In[11]:= disk = Image[Table[N[Boole[(i - 8.5)^2 + (j - 8.5)^2 <= 25]], {i, 1, 16}, {j, 1, 16}], "Real"];

In[12]:= ramp = Image[Table[N[(j - 1)/15], {i, 1, 16}, {j, 1, 16}], "Real"];

In[13]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[14]:= noise = Image[Table[N[Mod[i*37 + j*17, 101]]/101, {i, 1, 32}, {j, 1, 32}], "Real"];

In[15]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[16]:= sky = Image[Table[{N[0.15 + 0.7 (16 - i)/16], N[0.35 + 0.45 (16 - i)/16], N[0.85 - 0.35 (16 - i)/16]}, {i, 1, 16}, {j, 1, 24}], "Real"];

In[17]:= vol = Image3D[Table[N[Mod[z*7 + y*13 + x*3, 97]]/97, {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[18]:= volb = Image3D[Table[N[Boole[x <= 6 && y <= 5]], {z, 1, 8}, {y, 1, 10}, {x, 1, 12}], "Real"];

In[19]:= Image3DQ[disk]
Out[19]= False

In[20]:= Image3DQ[ramp]
Out[20]= False

In[21]:= Image3DQ[zone]
Out[21]= False

In[22]:= Image3DQ[noise]
Out[22]= False

In[23]:= Image3DQ[sky]
Out[23]= False

In[24]:= Image3DQ[vol]
Out[24]= True

In[25]:= Image3DQ[volb]
Out[25]= True

In[26]:= Image3DQ[Import[Export["/tmp/mathilda_ex.png", rgb]]]
Out[26]= False
```

### Applications (2)

```mathematica
In[27]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[28]:= Table[Image3DQ[GaussianFilter[chk, r]], {r, 1, 3}]
Out[28]= {False, False, False}
```

### Properties & Relations (4)

```mathematica
In[29]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[30]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[31]:= Image3DQ[chk] === Image3DQ[GaussianFilter[chk, 1]]
Out[31]= True

In[32]:= Image3DQ[rgb] === Image3DQ[ImagePad[rgb, 2]]
Out[32]= True
```

### Neat Examples (2)

```mathematica
In[33]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[34]:= Image3DQ[zone]
Out[34]= False
```

### Applications (3)

A 2-slice volume is a valid Image3D

```mathematica
In[35]:= Image3DQ[Image3D[{{{0., 1.}, {1., 0.}}, {{1., 0.}, {0., 1.}}}]]
Out[35]= True
```

A plane is not a volume

```mathematica
In[36]:= Image3DQ[Image[{{0., 1.}, {1., 0.}}]]
Out[36]= False
```

A bare array is not yet an Image3D

```mathematica
In[37]:= Image3DQ[{{{0., 1.}, {1., 0.}}}]
Out[37]= False
```

## Implementation notes

**Algorithm.** `builtin_image3dq` is a predicate: it calls `image3d_info` on its argument and
returns `True` or `False`. `image3d_info` accepts only the canonical volumetric form
`Image3D[data, "type"]` — the head must be `Image3D`, the second argument one of the canonical
type names, and the data a rectangular `depth × height × width` array (or
`depth × height × width × channels` for colour). `img3_shape_fast` checks rectangularity by
walking only slice and row lengths, and `img_data_storable` rejects complex storage; a packed
buffer answers its shape in `O(1)` from its dims. A malformed `Image3D[...]` returns `NULL`
from its constructor and so stays unevaluated, which is exactly why a separate predicate is
needed: an unevaluated `Image3D["hello"]` and a valid one would otherwise be
indistinguishable.

`Image3DQ` is the volumetric twin of `ImageQ`; `ImageQ` is `False` for a volume and
`Image3DQ` `False` for a plane, so the two partition the image heads cleanly.

**Data structures.** Reads the `Expr` tree without copying it. For a packed `Image3D` the
voxels live in a visible rank-3 (grey) or rank-4 (colour) `NDArray`, slices outermost, indexed
`data[[z, y, x]]`.

**Complexity / limits.** `O(1)` for a packed volume's shape plus the storable scan; a nested
volume pays `O(depth · height)` for the rectangularity check and `O(voxels)` for the storable
check. Returns a boolean, so it is registered packed-aware (the probe never materialises the
buffer it inspects).

**Attributes:** `Protected`.

## References

- Source: [`src/image.c`](https://github.com/stblake/mathilda/blob/main/src/image.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`Image3DQ` is the volumetric counterpart of `ImageQ`, and the two are disjoint: `ImageQ`
answers `False` for a volume and `Image3DQ` `False` for a plane.

It is the way to test validity because malformed input to `Image3D` stays **unevaluated**
rather than failing — a bare nested array has not yet been through the `Image3D` constructor,
so it is not in canonical form and `Image3DQ` reports `False`.

A volume is stored `depth × height × width` (slices outermost, indexed `data[[z, y, x]]`),
while `ImageDimensions` reports it fully reversed as `{width, height, depth}`.
