# $Pre

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$Pre
    is a global variable whose value, if set, is applied to every
    input expression after parsing and before standard evaluation.

Unless $Pre is assigned to a head with HoldAll, the wrapped
expression is evaluated before $Pre sees it -- which makes the
effect indistinguishable from $Post.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/repl_hooks.c`](https://github.com/stblake/mathilda/blob/main/src/repl_hooks.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
