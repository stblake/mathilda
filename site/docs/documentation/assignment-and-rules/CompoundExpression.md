# CompoundExpression

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
expr1; expr2; ... evaluates its arguments in sequence, returning the last result.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_compoundexpression` (in `src/core.c`) walks its argument list in order, calling `evaluate()` on each `expr;` clause and freeing the previous result before storing the next, so only the value of the last clause survives. An empty `CompoundExpression[]` returns the symbol `Null`. The head carries `ATTR_HOLDALL` (registered in `src/attr.c`), so the evaluator does not pre-evaluate the arguments — the builtin evaluates them itself, one at a time, which is what makes `a; b; c` sequence side effects (e.g. assignments installed by `Set`) in source order rather than all at once. After each evaluation the result head is inspected for the control-flow markers `Return`, `Break`, `Continue`, `Throw`, `Abort`, and `Quit`; if any is seen the loop breaks early and returns that wrapper unevaluated so an enclosing construct can act on it.

**Data structures.** Plain iteration over `res->data.function.args`; a single `last_val` accumulator holds the running result.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
