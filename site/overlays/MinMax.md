### Worked examples

```mathematica
In[1]:= MinMax[{3, 1, 4, 1, 5, 9, 2, 6}]  (* {smallest, largest} in one shot *)
In[2]:= MinMax[<|a -> 3, b -> 1, c -> 4|>]  (* over an association, uses the values *)
```

### Notes

`MinMax[list]` is `{Min[list], Max[list]}` computed without writing the two calls by
hand. It threads onto the packed-array fast path, so a large numeric vector stays on
the buffer and both extrema come from one data pass each rather than from 10⁶ boxed
elements.
