### Worked examples

```mathematica
In[1]:= ImageData[ImageReflect[Image[{{0., 1.}, {2., 3.}}]]]  (* default: flip top to bottom *)
```

```mathematica
In[1]:= ImageData[ImageReflect[Image[{{0., 1.}, {2., 3.}}], Left]]  (* the Left side flips left to right *)
```

```mathematica
In[1]:= ImageData[ImageReflect[ImageReflect[Image[{{0., 1.}, {2., 3.}}]]]]  (* reflecting twice is the identity, bit for bit *)
```

### Notes

`ImageReflect[image]` reflects top-to-bottom; `ImageReflect[image, Left]` (or `Right`) reflects
left-to-right, and `Top`/`Bottom` is the vertical reflection again — either name of a pair
selects the same axis, since reflecting to the top and to the bottom are one operation. For an
`Image3D`, `Front`/`Back` select the depth axis; those two decline on a plane rather than being
reinterpreted.

A reflection is a pure index permutation, so it interpolates nothing: reflecting twice about
the same axis is the identity bit for bit, and reflections about different axes commute
exactly.
