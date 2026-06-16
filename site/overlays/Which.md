### Worked examples

```mathematica
In[1]:= Which[False, 1, True, 2, True, 3]
Out[1]= 2
```

`Which` has attribute `HoldAll`, so it is the natural tool for piecewise function definitions — the tests are examined in order and only the matching branch's value is evaluated:

```mathematica
In[1]:= sign[x_] := Which[x < 0, -1, x == 0, 0, x > 0, 1];
        {sign[-7], sign[0], sign[42]}
Out[1]= {-1, 0, 1}
```

When a test cannot be decided (neither `True` nor `False`), `Which` returns itself unevaluated from that test onward, preserving the symbolic conditional rather than guessing:

```mathematica
In[1]:= Which[x > 0, pos, x < 0, neg, True, zero]
Out[1]= Which[x > 0, pos, x < 0, neg, True, zero]
```

### Notes

`Which[test1, value1, test2, value2, ...]` evaluates each test in turn and returns the value for the first test yielding `True`. With `HoldAll`, tests and values stay unevaluated until examined. A leftover undecidable test returns an unevaluated `Which`; if all tests are `False` (or none are given) the result is `Null`. Use `True` as the final test for a default.
