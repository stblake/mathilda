### Worked examples

```mathematica
In[1]:= AllTrue[{2, 4, 6}, EvenQ]  (* every element passes the test *)
```

```mathematica
In[1]:= AllTrue[{1, 2, 3, 4}, # > 0 &]  (* a pure function as the test *)
```

```mathematica
In[1]:= AllTrue[{2, 3, 4}, EvenQ]  (* one odd element is enough for False *)
```

```mathematica
In[1]:= AllTrue[<|a -> 2, b -> 4|>, EvenQ]  (* over an association the VALUES are tested *)
```

```mathematica
In[1]:= AllTrue[{}, EvenQ]  (* vacuously True on the empty list *)
```

### Notes

`AllTrue[list, test]` is `True` exactly when `test[e]` is `True` for every element
`e`, and short-circuits on the first `False`. It is the universal quantifier to
`AnyTrue`'s existential and `NoneTrue`'s negation. Over an association the values
are tested, not the keys. The empty list is vacuously `True`.

If some `test[e]` returns a value that is neither `True` nor `False`, the whole
call is left unevaluated rather than guessed — so a partially-symbolic list flows
through unchanged. Over a boolean packed array with `TrueQ` or `Identity` as the
test, the answer comes from an early-exit scan of the raw buffer.
