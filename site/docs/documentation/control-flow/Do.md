# Do

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Do[expr, n] evaluates expr n times.
Do[expr, {i, imax}] evaluates expr with i successively taking on values 1 through imax.
Do[expr, {i, imin, imax}] starts with i = imin.
Do[expr, {i, imin, imax, di}] uses steps di.
Do[expr, {i, {i1, i2, ...}}] uses the successive values i1, i2, ....
Do[expr, {n}] evaluates expr n times with no iteration variable.
Do[expr, iter1, iter2, ...] iterates over multiple variables, with the rightmost varying fastest.
Do has attribute HoldAll: expr and the iterator specifications are held unevaluated until each iteration.
Break[] inside expr exits the innermost Do loop.
Continue[] inside expr skips the rest of expr and proceeds to the next iteration.
Return[v] inside expr causes the enclosing function to yield v; Do itself returns Null.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_do` is `ATTR_HOLDALL`, so the body and iterator spec are re-evaluated on every pass rather than once at call time. A multi-spec `Do[expr, s1,...,sk]` is rewritten as nested two-spec `Do[Do[expr, sk], s1,...,s_{k-1}]` and handed back to the evaluator, reducing every case to a single iterator. The spec is parsed by the shared `iter_spec_parse` into an `IterSpec` (kinds `COUNT` for `n`/`{n}`, `RANGE` for `{i,imin,imax,di}` with defaults, `LIST` for `{i,{...}}`); numeric bounds are resolved to doubles by `iter_spec_resolve_numeric`, with `Infinity` allowed for unbounded loops.

**Variable localization.** Before looping, `iter_spec_shadow` saves and clears the iterator symbol's `own_values`; each iteration binds the current value via `symtab_add_own_value`, and `iter_spec_restore` frees the per-iteration binding chain and restores the original OwnValue afterward — a manual mimic of Mathematica's iterator scoping. For exact ranges the bound value `curr_e` is advanced with `Plus[curr_e, di_e]` evaluated each step (keeping integer/rational exactness), while a parallel `double val` drives the termination comparison (with a `1e-14` tolerance).

**Control flow.** After each `evaluate(body)`, `iter_flow_classify` (boundary head `SYM_Do`) maps the result to break / continue / return-value / propagate (Throw/Abort/Quit/foreign Return). `Continue` in a range loop still advances the counter before re-testing. Returns the Return payload if any, else `Null`.

- `HoldAll`, evaluating its body only after arguments are substituted.
- Employs exact dynamic iteration identical to `Table` but discards the evaluated results, returning `Null`.
- Supports explicit break states (`Return`, `Break`, `Continue`, `Throw`, `Abort`, `Quit`).
- Can execute an infinite loop using `Do[expr, Infinity]`.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/iter.c`](https://github.com/stblake/mathilda/blob/main/src/iter.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= (s = 0; Do[s = s + i, {i, 5}]; s)
Out[1]= 15

In[2]:= (p = 1; Do[p = p*2, {3}]; p)
Out[2]= 8
```

### Notes

`Do[expr, {i, imax}]` runs `expr` with `i` taking values `1` through `imax`; `Do[expr, n]` simply repeats `expr` `n` times. `Do` returns `Null`, so it is used for side effects — read the accumulated value out of a variable afterwards.
