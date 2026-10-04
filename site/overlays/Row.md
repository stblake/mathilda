### Worked examples

```mathematica
In[1]:= Row[{a, b, c}]  (* elements concatenated with no separator *)
```

```mathematica
In[1]:= Row[{1, 2, 3}, ", "]  (* a string inserted between successive elements *)
```

```mathematica
In[1]:= Row[{"x", "=", 5}]  (* strings print without their quotes *)
```

### Notes

`Row[{e1, e2, ...}]` displays the elements concatenated left to right, with
strings shown without quotes; `Row[{...}, s]` inserts the string `s` between
successive elements. It is a pure presentation head with no computational value,
and its main internal use is assembling a custom display for `NumberForm`'s
`NumberFormat` option — the pieces `(mantissa, "10", exponent)` are handed to the
`NumberFormat` function, which typically wraps them in a `Row`.
