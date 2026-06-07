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

- `HoldAll`, evaluating its body only after arguments are substituted.
- Employs exact dynamic iteration identical to `Table` but discards the evaluated results, returning `Null`.
- Supports explicit break states (`Return`, `Break`, `Continue`, `Throw`, `Abort`, `Quit`).
- Can execute an infinite loop using `Do[expr, Infinity]`.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
