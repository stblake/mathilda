### Worked examples

```mathematica
In[1]:= UnitVector[4, 2]  (* a 1 in position 2 of a length-4 vector *)
```

```mathematica
In[1]:= UnitVector[2]  (* one argument means dimension 2 *)
```

```mathematica
In[1]:= UnitVector[5, 5]  (* the last basis vector *)
```

```mathematica
In[1]:= Total[UnitVector[100, 7]]  (* exactly one nonzero component *)
```

### Notes

`UnitVector[n, k]` is the length-`n` vector with a `1` in position `k` and `0`
everywhere else; `UnitVector[k]` is the two-dimensional case (`UnitVector[2, k]`).
Both `n` and `k` must be positive integers with `1 <= k <= n`.

Components are exact integers by default. `WorkingPrecision -> MachinePrecision`
gives machine reals, and a digit count gives arbitrary-precision reals. At scale
the integer and machine-real forms are returned as a packed array, so
constructing a very long unit vector costs a buffer fill rather than one boxed
element per zero.
