### Worked examples

```mathematica
In[1]:= AnyTrue[{1, 3, 5}, EvenQ]  (* no element passes, so False *)
```

```mathematica
In[1]:= AnyTrue[{1, 2, 3}, EvenQ]  (* one even element is enough for True *)
```

```mathematica
In[1]:= AnyTrue[{1, 3, 5}, # > 4 &]  (* a pure function as the test *)
```

### Notes

`AnyTrue[list, test]` is `True` exactly when `test[e]` is `True` for at least one
element `e`, and short-circuits on the first match. It is the existential
quantifier dual to `AllTrue`; `NoneTrue` is its negation. Over an association the
values are tested. The empty list gives `False`.

A `test[e]` whose value is neither `True` nor `False` leaves the whole call
unevaluated. Over a boolean packed array with `TrueQ` or `Identity` as the test,
the answer is an early-exit scan of the raw byte buffer.
