### Worked examples

```mathematica
In[1]:= NoneTrue[{1, 3, 5}, EvenQ]  (* no element is even, so True *)
```

```mathematica
In[1]:= NoneTrue[{1, 2, 3}, EvenQ]  (* one even element makes it False *)
```

```mathematica
In[1]:= NoneTrue[Range[5], # > 10 &]  (* nothing exceeds 10 *)
```

### Notes

`NoneTrue[list, test]` is `True` exactly when `test[e]` is `True` for no element
`e` — the negation of `AnyTrue` — and short-circuits to `False` on the first
match. Over an association the values are tested. The empty list gives `True`.

A `test[e]` that is neither `True` nor `False` leaves the whole call unevaluated.
Over a boolean packed array with `TrueQ` or `Identity` as the test, the check is a
single early-exit scan of the raw buffer.
