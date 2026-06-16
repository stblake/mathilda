# Which

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Which[test1, value1, test2, value2, ...]
    evaluates each test_i in turn, returning the corresponding value_i
    for the first test that yields True.
Which has attribute HoldAll, so the tests and values are not
evaluated until Which examines them.
If a test evaluates to neither True nor False, a Which object
containing that test (in evaluated form) and the remaining
elements is returned unevaluated.
If all tests evaluate to False (or no tests are supplied), Which
returns Null.
Use True as the final test to supply a default value.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Which[False, a, True, b]
Out[1]= b

In[2]:= Which[1 < 0, a, x == 0, b, 0 < 1, c]
Out[2]= Which[x == 0, b, 0 < 1, c]

In[3]:= Which[a == 1, x, a == 2, b] /. a -> 2
Out[3]= b

In[4]:= sign[x_] := Which[x < 0, -1, x > 0, 1, True, Indeterminate]; sign /@ {-2, 0, 3}
Out[4]= {-1, Indeterminate, 1}
```

## Implementation notes

**Algorithm.** `builtin_which` is `ATTR_HOLDALL`, so every argument arrives unevaluated. It requires an even argument count (test/value pairs); an odd count returns `NULL` (unevaluated usage error) and `Which[]` returns `Null`. It walks the pairs in order, calling `evaluate` on each test in turn. A test that reduces to the interned symbol `True` makes the handler return a copy of the corresponding held value (which the outer evaluator then reduces). A test that reduces to `False` is dropped and iteration continues. An inconclusive test (anything else) short-circuits the scan: the builtin rebuilds `Which[t_i_eval, v_i, ...remaining...]` with the inconclusive test in its already-evaluated form and the remaining arguments copied unevaluated, so re-evaluation does not redo earlier `False` tests. If every test is `False`, it returns `Null`.

- Has attribute `HoldAll`; tests and values are held until `Which` examines them.
- If every `test_i` evaluates to `False`, `Which` returns `Null`. `Which[]` (no arguments) likewise yields `Null`.
- If a `test_i` evaluates to something other than `True` or `False`, a `Which` containing that test (in evaluated form) plus the remaining elements is returned unevaluated.
- A trailing test of `True` acts as a default clause.
- An odd number of arguments is a usage error; the expression is returned unevaluated.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/cond.c`](https://github.com/stblake/mathilda/blob/main/src/cond.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)

## Notes & additional examples

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
