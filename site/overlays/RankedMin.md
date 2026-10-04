### Worked examples

```mathematica
In[1]:= RankedMin[{12, 13, 11}, 2]  (* the 2nd smallest element *)
```

```mathematica
In[2]:= RankedMin[{Pi, Sqrt[2], E, 3}, 3]  (* symbolic reals order by value *)
```

```mathematica
In[3]:= RankedMin[{12, 13, 11}, -1]  (* the 1st largest is Max *)
```

### Notes

`RankedMin[list, n]` is the `n`-th smallest element — the order statistic between
`Min` and `Max`. A negative index counts from the top, so `RankedMin[list, -n]`
is the `n`-th largest, `RankedMin[list, 1]` is `Min[list]`, and
`RankedMin[list, -1]` is `Max[list]`. It returns a definite result whenever every
element is a real number, including symbolic real constants (`Pi`, `E`,
`Sqrt[2]`), which order by value, with `±Infinity` ranking as `±∞`; the element
comes back in its exact form. A non-real element, an empty list, or `|n|` out of
range leaves the call unevaluated. Selection is an `O(n)` quickselect (int64
exact), with a packed-array fast path and a `Compile[]` lowering.
