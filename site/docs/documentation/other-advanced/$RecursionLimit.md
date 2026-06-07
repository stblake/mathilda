# $RecursionLimit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$RecursionLimit
    gives the maximum length of the evaluation stack -- the maximum
    number of nested invocations of the evaluator that can occur.

Assigning a positive integer N (>= 20) updates the limit; smaller
values are rejected with a $RecursionLimit::limset message.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

A user-visible system variable seeded in `recursion_limit_init` (`src/eval.c`) with an integer OwnValue (the default C-level limit). The evaluator enforces it with a static depth counter `g_eval_depth` incremented per `evaluate_step`; when depth exceeds the limit it aborts the call and emits `$RecursionLimit::reclim`. Assignment is special-cased in the `Set` handler: `$RecursionLimit = expr` routes through `recursion_limit_set` (called from `eval.c`'s assignment path), which rejects values below a hard floor with `$RecursionLimit::limset` and otherwise reinstalls the new integer as the symbol's OwnValue so the C counter and the symbol stay in sync.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
