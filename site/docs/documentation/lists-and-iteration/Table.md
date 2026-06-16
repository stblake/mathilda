# Table

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Table[expr, n]
    generates a list of n copies of expr.
Table[expr, {i, imax}]
    generates a list of the values of expr with i running from 1 to imax.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Table[i^2, {i, 4}]
Out[1]= {1, 4, 9, 16}
```

## Implementation notes

**Algorithm.** `builtin_table` (in `src/list.c`) implements iterator-driven list construction. `Table` carries `HoldAll | Protected` (set in `src/attr.c`), so the body and iterator specs arrive unevaluated. Multi-iterator forms `Table[expr, spec1, ..., specN]` are handled by *rewriting*: the handler peels off the outermost iterator, wraps the remaining iterators in a fresh `Table[Table[expr, ..., specN], spec1]`, and recurses through the evaluator — so nesting depth equals iterator count and the rightmost spec varies fastest.

For a single iterator the spec is parsed by the shared `iter_spec_parse` (`src/iter.c`) into an `IterSpec` discriminated by `kind`: `ITER_KIND_COUNT` (`{n}` or bare `n` — repeat the body), `ITER_KIND_LIST` (`{i, {v1, v2, ...}}` — iterate over explicit values), or `ITER_KIND_RANGE` (`{i, imax}`, `{i, imin, imax}`, `{i, imin, imax, di}`). Range bounds are resolved to doubles via `iter_spec_resolve_numeric` (with `allow_inf=false`, since `Table` never iterates to `Infinity`).

**Localization.** The iteration variable is dynamically scoped, not renamed. `iter_spec_shadow(var)` saves and clears the symbol's `own_values` chain; each step calls `symtab_add_own_value` to bind the current index, `evaluate`s the body, then the next index is computed symbolically (`evaluate(Plus[curr, di])`) so exact integers/rationals are preserved while a parallel `double` accumulator drives loop termination (with a `1e-14` tolerance and a 1,000,000-step safety cap). `iter_spec_restore` reinstalls the saved `own_values` afterward, so the variable's outer value is untouched.

**Data structures.** Results accumulate in a geometrically grown `Expr**` buffer, finally wrapped as `List[...]`.

- `HoldAll`: `expr` is evaluated once for each step.
- Supports nested iterators to create matrices.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Table[i^2, {i, 1, 5}]
Out[1]= {1, 4, 9, 16, 25}
```

```mathematica
In[1]:= Table[i + j, {i, 1, 2}, {j, 1, 3}]
Out[1]= {{2, 3, 4}, {3, 4, 5}}
```

```mathematica
In[1]:= Table[x, 4]
Out[1]= {x, x, x, x}
```

```mathematica
In[1]:= Table[i, {i, 0, 1, 1/2}]
Out[1]= {0, 1/2, 1}
```

```mathematica
In[1]:= Table[Fibonacci[n], {n, 1, 12}]
Out[1]= {1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144}
```

```mathematica
In[1]:= Table[Sum[1/k, {k, 1, n}], {n, 1, 5}]
Out[1]= {1, 3/2, 11/6, 25/12, 137/60}
```

```mathematica
In[1]:= Table[If[i == j, 1, 0], {i, 1, 3}, {j, 1, 3}]
Out[1]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
```

### Notes

The single-argument iterator `{i, imax}` runs `i` from 1 to `imax`; `{i, imin,
imax}` and `{i, imin, imax, step}` give explicit bounds and a step. The step may
be an exact rational, so the values stay exact (`{0, 1/2, 1}`). Multiple iterator
specifications nest: the leftmost varies slowest, producing a list of lists.
`Table[expr, n]` with a plain count simply repeats `expr` `n` times. `Table` holds
its arguments, so the body is only evaluated as each iterator value is assigned —
the body may itself be a `Sum`, `D`, or any other computation (giving exact
harmonic numbers, identity matrices, and the like).
