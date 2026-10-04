# ImageQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ImageQ[expr] gives True if expr is a valid image in canonical form, and False otherwise. Malformed input to Image stays unevaluated, so ImageQ is how validity is tested.`**

## Examples (38)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (10)

```mathematica
In[1]:= ImageQ[Image[{{0., 1.}, {1., 0.}}]]
Out[1]= True

In[2]:= ImageQ[{{0., 1.}, {1., 0.}}]
Out[2]= False

In[3]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[4]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[5]:= bit = Image[Table[Boole[Mod[i + j, 2] == 0], {i, 1, 8}, {j, 1, 8}]];

In[6]:= byte = Image[Table[Mod[i*13 + j*7, 256], {i, 1, 16}, {j, 1, 16}]];

In[7]:= ImageQ[chk]
Out[7]= True

In[8]:= ImageQ[rgb]
Out[8]= True

In[9]:= ImageQ[bit]
Out[9]= True

In[10]:= ImageQ[byte]
Out[10]= True
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

In[19]:= ImageQ[disk]
Out[19]= True

In[20]:= ImageQ[ramp]
Out[20]= True

In[21]:= ImageQ[zone]
Out[21]= True

In[22]:= ImageQ[noise]
Out[22]= True

In[23]:= ImageQ[sky]
Out[23]= True

In[24]:= ImageQ[vol]
Out[24]= False

In[25]:= ImageQ[volb]
Out[25]= False

In[26]:= ImageQ[Import[Export["/tmp/mathilda_ex.png", rgb]]]
Out[26]= True
```

### Applications (2)

```mathematica
In[27]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[28]:= Table[ImageQ[GaussianFilter[chk, r]], {r, 1, 3}]
Out[28]= {True, True, True}
```

### Properties & Relations (4)

```mathematica
In[29]:= chk = Image[Table[If[Mod[Quotient[i - 1, 2] + Quotient[j - 1, 2], 2] == 0, 0., 1.], {i, 1, 16}, {j, 1, 16}], "Real"];

In[30]:= rgb = Image[Table[{N[i/16], N[j/16], 0.5}, {i, 1, 16}, {j, 1, 16}], "Real"];

In[31]:= ImageQ[chk] === ImageQ[GaussianFilter[chk, 1]]
Out[31]= True

In[32]:= ImageQ[rgb] === ImageQ[ImagePad[rgb, 2]]
Out[32]= True
```

### Neat Examples (2)

```mathematica
In[33]:= zone = Image[Table[N[(1 + Cos[((i - 16)^2 + (j - 16)^2)/40.])/2], {i, 1, 32}, {j, 1, 32}], "Real"];

In[34]:= ImageQ[zone]
Out[34]= True
```

### Applications (4)

```mathematica
In[35]:= ImageQ[Image[{{0., 1.}, {1., 0.}}]]
Out[35]= True

In[36]:= ImageQ[5]
Out[36]= False

In[37]:= ImageQ[Image[{{0., 1.}, {1.}}]]
Out[37]= False

In[38]:= ImageQ[Image3D[{{{0., 1.}}, {{1., 0.}}}]]
Out[38]= False
```

## Implementation notes

**Algorithm.** `builtin_imageq` is a one-argument predicate that returns
`True`/`False` and nothing symbolic: it calls `image_info(arg, NULL, NULL, NULL,
NULL)` and wraps the boolean. `image_info` is the single validity gate the whole
subsystem shares — it checks that the argument is a two-argument `Image[data,
type]` whose second argument is a type string (`"Bit"`, `"Byte"`, `"Bit16"`,
`"Real"`), that `img_shape_fast` finds a rectangular height × width (× channels)
array, and that `img_data_storable` confirms every leaf is in the range the type
fixes. `ImageQ` is the companion to the fact that malformed input to `Image[...]`
is left **unevaluated** rather than erroring — the constructor returns `NULL` for
ragged, non-numeric or complex data, so the head stays `Image[...]` and `ImageQ`
is how a caller tests whether that happened. A volume is deliberately `False`
here (it is `Image3DQ` that accepts one), since the two ranks are distinct
objects.

**Data structures.** Reads the canonical `Image` node — an `EXPR_FUNCTION` with
head `Image`, argument 0 the pixel array (normally a packed NDArray buffer),
argument 1 the type string. No buffer is loaded or copied; validation walks only
the shape and leaf scalars. `ImageQ` is on `pack.c`'s `AWARE` list, so a packed
pixel buffer is inspected in place rather than being unpacked into boxed `Expr`
nodes.

**Complexity / limits.** `O(1)` for the structural checks plus one pass over the
leaves for `img_data_storable`; no allocation. Returns `False`, never
unevaluated, for every non-image — it is a total predicate.

**Attributes:** `Protected`.

## References

**See also:** [Head](../../structural-manipulation/Head/), [Image](../../image-processing/Image/)

- Source: [`src/image.c`](https://github.com/stblake/mathilda/blob/main/src/image.c)
- Specification: [`docs/spec/builtins/image-processing.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/image-processing.md)
- Tests: [`tests/test_image.c`](https://github.com/stblake/mathilda/blob/main/tests/test_image.c)

## Notes & additional examples

### Notes

`ImageQ[expr]` gives `True` exactly when `expr` is a valid image in canonical
`Image[data, type]` form, and `False` otherwise — it never returns unevaluated,
so it is a total predicate usable in a pattern test or condition.

Validity is decided by the subsystem's shared `image_info` gate: a two-argument
`Image`, a recognised type string, a rectangular pixel array, and every stored
value inside the range the type fixes. This matters because the `Image[...]`
constructor leaves **malformed** input (ragged, non-numeric, complex) unevaluated
rather than erroring — the head stays `Image[...]`, and `ImageQ` is the way to
test whether the construction actually succeeded. A volumetric `Image3D` is
`False`; its own predicate is `Image3DQ`.
