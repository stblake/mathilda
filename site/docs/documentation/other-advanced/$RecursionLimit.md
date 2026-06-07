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

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
