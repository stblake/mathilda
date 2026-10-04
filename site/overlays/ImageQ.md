### Worked examples

```mathematica
In[1]:= ImageQ[Image[{{0., 1.}, {1., 0.}}]]
```

```mathematica
(* a plain number is not an image *)
In[1]:= ImageQ[5]
```

```mathematica
(* malformed Image[...] stays unevaluated, so ImageQ reports it False *)
In[1]:= ImageQ[Image[{{0., 1.}, {1.}}]]
```

```mathematica
(* a volume is False here -- use Image3DQ for Image3D *)
In[1]:= ImageQ[Image3D[{{{0., 1.}}, {{1., 0.}}}]]
```

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
