# Which

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Which[test1, value1, test2, value2, ...]`**

evaluates each test\_i in turn, returning the corresponding value\_i for the first test that yields True.

<details>
<summary>Notes</summary>

Which has attribute HoldAll, so the tests and values are not evaluated until Which examines them. If a test evaluates to neither True nor False, a Which object containing that test (in evaluated form) and the remaining elements is returned unevaluated. If all tests evaluate to False (or no tests are supplied), Which returns Null. Use True as the final test to supply a default value.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

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

### Applications (2)

```mathematica
In[5]:= Which[False, 1, True, 2, True, 3]
Out[5]= 2

In[6]:= Which[x > 0, pos, x < 0, neg, True, zero]
Out[6]= Which[x > 0, pos, x < 0, neg, True, zero]
```

## Implementation notes

**Algorithm.** `builtin_which` is `ATTR_HOLDALL`, so every argument arrives unevaluated. It requires an even argument count (test/value pairs); an odd count returns `NULL` (unevaluated usage error) and `Which[]` returns `Null`. It walks the pairs in order, calling `evaluate` on each test in turn. A test that reduces to the interned symbol `True` makes the handler return a copy of the corresponding held value (which the outer evaluator then reduces). A test that reduces to `False` is dropped and iteration continues. An inconclusive test (anything else) short-circuits the scan: the builtin rebuilds `Which[t_i_eval, v_i, ...remaining...]` with the inconclusive test in its already-evaluated form and the remaining arguments copied unevaluated, so re-evaluation does not redo earlier `False` tests. If every test is `False`, it returns `Null`.

- Has attribute `HoldAll`; tests and values are held until `Which` examines them.
- If every `test_i` evaluates to `False`, `Which` returns `Null`. `Which[]` (no arguments) likewise yields `Null`.
- If a `test_i` evaluates to something other than `True` or `False`, a `Which` containing that test (in evaluated form) plus the remaining elements is returned unevaluated.
- A trailing test of `True` acts as a default clause.
- An odd number of arguments is a usage error; the expression is returned unevaluated.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/cond.c`](https://github.com/stblake/mathilda/blob/main/src/cond.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
- Tests: [`tests/test_catch_throw.c`](https://github.com/stblake/mathilda/blob/main/tests/test_catch_throw.c)
- Tests: [`tests/test_cond.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cond.c)

## Notes & additional examples

### Notes

`Which[test1, value1, test2, value2, ...]` evaluates each test in turn and returns the value for the first test yielding `True`. With `HoldAll`, tests and values stay unevaluated until examined. A leftover undecidable test returns an unevaluated `Which`; if all tests are `False` (or none are given) the result is `Null`. Use `True` as the final test for a default.
