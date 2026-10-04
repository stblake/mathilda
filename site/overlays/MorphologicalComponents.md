### Worked examples

```mathematica
(* four corner pixels, each isolated under 8-connectivity: labels 1..4 in raster order *)
In[1]:= MorphologicalComponents[Image[{{1, 0, 1}, {0, 0, 0}, {1, 0, 1}}]]
```

```mathematica
(* a diagonal touch is ONE component under the default 8-connectivity *)
In[1]:= MorphologicalComponents[Image[{{1, 0}, {0, 1}}]]
```

```mathematica
(* ... and TWO components under 4-connectivity *)
In[1]:= MorphologicalComponents[Image[{{1, 0}, {0, 1}}], 0, CornerNeighbors -> False]
```

```mathematica
(* contiguous labels mean Max is the component count *)
In[1]:= Max[MorphologicalComponents[Image[{{1, 0, 1}, {1, 0, 1}}]]]
```

### Notes

`MorphologicalComponents[image]` labels the connected components of the
foreground, giving an **integer matrix** with background `0` and components
numbered `1..k` in raster order of first appearance.
`MorphologicalComponents[image, t]` takes pixels above `t` as foreground (default
`0`), and `CornerNeighbors -> False` uses 4-connectivity instead of the default
8.

The two-pass union-find assigns provisional labels in raster order, records
equivalences (a U-shape is the case that needs the second pass), then relabels to
`1..k` with no gaps — so `Max` of the result is the component count.
Connectivity is the only discriminating property: two pixels touching at a corner
are one component under 8 and two under 4. The result is a matrix, not an `Image`,
on purpose: `Image` would type-infer a `1..12` label array as `"Byte"` and
`ImageData` would divide every label by 255.
