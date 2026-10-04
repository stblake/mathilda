### Worked examples

```mathematica
In[1]:= RankedMax[{12, 13, 11}, 1]  (* the largest element is Max *)
```

```mathematica
In[2]:= RankedMax[{2.5, E, 12, 15, 485}, -2]  (* the 2nd smallest, via a negative index *)
```

```mathematica
In[3]:= RankedMax[{Infinity, 5, Infinity, -Infinity}, 2]  (* ties and infinities rank by value *)
```

### Notes

`RankedMax[list, n]` is the `n`-th largest element, the mirror of `RankedMin`:
`RankedMax[list, n]` equals `RankedMin[list, -n]`, so `RankedMax[list, 1]` is
`Max[list]` and `RankedMax[list, -1]` is `Min[list]`; a negative index
`RankedMax[list, -n]` gives the `n`-th smallest. It returns a definite result
whenever every element is a real number (symbolic real constants order by value,
`±Infinity` as `±∞`) and returns the element in its exact form; a non-real
element, an empty list, or an out-of-range index leaves the call unevaluated.
Selection is an `O(n)` quickselect sharing `RankedMin`'s core, with a
packed-array fast path and a `Compile[]` lowering.
