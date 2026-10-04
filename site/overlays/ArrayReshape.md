### Worked examples

```mathematica
In[1]:= ArrayReshape[Range[12], {3, 4}]  (* fill a 3x4 matrix row by row *)
```

```mathematica
In[1]:= ArrayReshape[Range[6], {2, 2}]  (* too many elements: the extras are dropped *)
```

```mathematica
In[1]:= ArrayReshape[{1, 2}, {2, 3}, 0]  (* too few: the tail is padded with 0 *)
```

```mathematica
In[1]:= Flatten[ArrayReshape[Range[24], {2, 3, 4}]] === Range[24]  (* reshaping preserves the flattened order *)
```

### Notes

`ArrayReshape[list, dims]` lays the fully flattened elements of `list` into a
rectangular `dims` array in row-major order. If `list` has more elements than the
shape needs, the extras are dropped; if it has fewer, the tail is filled with the
padding (default `0`, or a named scheme as a third argument). Up to the shared
length, `Flatten[ArrayReshape[list, dims]] == Flatten[list]`. Dimensions must be
non-negative integers.
