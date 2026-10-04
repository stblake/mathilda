### Worked examples

```mathematica
In[1]:= ArrayPad[{{1, 2}, {3, 4}}, 1]  (* one layer of zeros around a matrix *)
```

```mathematica
In[1]:= ArrayPad[{1, 2, 3, 4}, {2, 0}]  (* two before, none after *)
```

```mathematica
In[1]:= ArrayPad[{1, 2, 3}, 2, "Periodic"]  (* wrap the array around cyclically *)
```

```mathematica
In[1]:= ArrayPad[{1, 2, 3}, 2, "Reflected"]  (* mirror the edge values outward *)
```

### Notes

`ArrayPad[array, m]` pads `m` elements on every side of every level;
`ArrayPad[array, {m, n}]` puts `m` before and `n` after; and
`ArrayPad[array, {{m1,n1}, ...}]` gives per-level amounts. A **negative** amount
removes elements from that side instead.

The optional third argument is the padding: a constant (default `0`), a cyclic
list of constants, or a named scheme — `"Fixed"`, `"Periodic"`, `"Reflected"`,
`"Reversed"`, `"ReversedNegation"`, `"ReflectedDifferences"`,
`"ReversedDifferences"`, or `"Extrapolated"` (which honours
`InterpolationOrder`). The difference schemes need an axis of length at least 2
and otherwise emit `ArrayPad::mindimsize`.
