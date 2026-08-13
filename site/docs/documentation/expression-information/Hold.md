# Hold

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Hold[expr]`**

maintains expr in an unevaluated form.

**`Evaluate[expr] inside Hold overrides the hold and evaluates expr once.`**

<details>
<summary>Notes</summary>

Hold has attribute HoldAll: its arguments are not evaluated. Sequence expressions inside Hold are flattened; use HoldComplete to prevent this.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= Hold[1+1]
Out[1]= Hold[1 + 1]

In[2]:= ReleaseHold[Hold[1+1]]
Out[2]= 2

In[3]:= Hold[1+1, 2+2]
Out[3]= Hold[1 + 1, 2 + 2]
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| 100000 ops on 10 elements | 6.24 s | 18.8 s | 35.4 s |
| 10000 elementwise ops on 10^2 | 1.91 s | 1.54 s | 1.84 s |
| 100000 scalar additions | 0.411 s | 10.9 s | 1.21 s |
| 1000 ops on 10^3 elements | 0.176 s | 0.234 s | 0.474 s |
| 1 elementwise op on 10^6 | 0.166 s | 0.059 s | 0.14 s |
| 1 op on 10^6 elements | 0.114 s | 0.037 s | 0.105 s |

## Implementation notes

`Hold` is not a C builtin. It is registered in `attr.c`'s attribute table as `{"Hold", ATTR_HOLDALL | ATTR_PROTECTED}`, i.e. the symbol `Hold` simply carries the `ATTR_HOLDALL` attribute. When the evaluator (`eval.c`) processes `Hold[args...]` it reads that attribute before evaluating arguments and therefore leaves every argument unevaluated, returning the `Hold[...]` wrapper as-is. There is no per-head logic — holding falls entirely out of the generic attribute-driven argument-evaluation step. `HoldComplete`, `HoldPattern`, and `Unevaluated` are registered the same way with `ATTR_HOLDALLCOMPLETE`/`ATTR_HOLDALL`.

**Attributes:** `HoldAll`, `Protected`.

## References

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_catch_throw.c`](https://github.com/stblake/mathilda/blob/main/tests/test_catch_throw.c)
- Tests: [`tests/test_core_algebra.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core_algebra.c)
- Tests: [`tests/test_eval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval.c)
- Tests: [`tests/test_eval_timestamps.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval_timestamps.c)

## Notes & additional examples

### Notes

`Hold` has attribute `HoldAll`, so it keeps every argument unevaluated; wrap an argument in `Evaluate[...]` to force a single evaluation, and use `ReleaseHold` to strip the `Hold` and evaluate the contents.
