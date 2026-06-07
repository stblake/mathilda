# HoldAllComplete

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
HoldAllComplete
    is an attribute which specifies that all arguments to a function are not to be modified or looked at in any way in the process of evaluation.
HoldAllComplete prevents argument evaluation, Sequence flattening inside arguments, Unevaluated wrapper stripping, and application of Evaluate.
Evaluate cannot override HoldAllComplete.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
