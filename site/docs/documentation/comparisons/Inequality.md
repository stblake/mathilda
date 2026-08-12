# Inequality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Inequality[v0, op0, v1, op1, v2, ...] is the canonical form for a chained comparison such as a < b <= c. It returns True if every adjacent pair holds, False if any pair fails, and otherwise the residual chain with decidable pairs dropped.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= 1 < 2 < 3
Out[1]= True

In[2]:= 1 < 5 < 3
Out[2]= False

In[3]:= 2 <= 2 < 5
Out[3]= True
```

```mathematica
In[1]:= 1 < x < 3 < 5
Out[1]= 1 < x < 3
```

```mathematica
In[1]:= 2 < 3 < x < 1
Out[1]= 3 < x < 1
```

## Implementation notes

**Algorithm.** `builtin_inequality` handles the variadic chained-comparison head the parser emits for mixed chains like `a < b <= c == d`. The argument layout is strict: `2k+1` args, even slots values, odd slots bare operator symbols (`Less`/`LessEqual`/`Greater`/`GreaterEqual`/`Equal`); a malformed layout returns NULL. Each adjacent pair `(v_i, op_i, v_{i+1})` is resolved by `decide_pair`, which routes through `compare_numeric`/`expr_eq` and returns 1 (True), 0 (False), or -1 (undecidable). Any False pair short-circuits the whole chain to `False`. If all pairs are True → `True`. Otherwise the proven-True pairs are dropped and a residual `Inequality[...]` is rebuilt from the undecided pairs (collapsing to a single binary `op[a,b]` head when exactly one survives, and de-duplicating a shared boundary value so the residual stays a well-formed `2k+1` chain).

**Data structures.** A `bool* undecided` scratch array (one slot per pair) flags which pairs to keep; survivors are deep-copied into a fresh argument array that becomes the residual chain.

- Adjacent comparisons that evaluate to `True` are dropped; the result collapses
  to `True`, `False`, or a reduced `Inequality` over the undecided operands.
- An `Indeterminate` operand makes its adjacent pairs `False`, so the whole
  chain is `False`.

**Attributes:** `Protected`.

## References

- Source: [`src/comparisons.c`](https://github.com/stblake/mathilda/blob/main/src/comparisons.c)
- Specification: [`docs/spec/builtins/comparisons.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/comparisons.md)
- Tests: [`tests/test_deriv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_deriv.c)
- Tests: [`tests/test_expand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expand.c)
- Tests: [`tests/test_facpoly.c`](https://github.com/stblake/mathilda/blob/main/tests/test_facpoly.c)
- Tests: [`tests/test_factor_baseline.c`](https://github.com/stblake/mathilda/blob/main/tests/test_factor_baseline.c)

## Notes & additional examples

### Notes

A chained comparison such as `a < b <= c` parses to the canonical `Inequality` form and holds only when every adjacent pair holds. Mixed operators (e.g. `<=` and `<`) may be combined in a single chain. When some pairs are undecidable (because a value is symbolic), the decidable-and-true pairs are dropped and only the residual chain is returned — `1 < x < 3 < 5` collapses the true `3 < 5` and keeps `1 < x < 3`.
