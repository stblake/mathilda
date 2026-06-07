# BeginPackage

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
BeginPackage["ctx`"] sets the current context to "ctx`" and restricts
$ContextPath to {"ctx`", "System`"}, matching Mathematica's package
prologue.
BeginPackage["ctx`", {"need1`", ...}] additionally prepends the
listed contexts to $ContextPath.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/context.c`](https://github.com/stblake/mathilda/blob/main/src/context.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
