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

- Has attribute `HoldAll`; tests and values are held until `Which` examines them.
- If every `test_i` evaluates to `False`, `Which` returns `Null`. `Which[]` (no arguments) likewise yields `Null`.
- If a `test_i` evaluates to something other than `True` or `False`, a `Which` containing that test (in evaluated form) plus the remaining elements is returned unevaluated.
- A trailing test of `True` acts as a default clause.
- An odd number of arguments is a usage error; the expression is returned unevaluated.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
