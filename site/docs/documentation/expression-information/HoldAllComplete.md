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

`HoldAllComplete` is an attribute name, not a function. It maps to the bitflag `ATTR_HOLDALLCOMPLETE` (`attr_name_to_flag` / `get_attributes` in `src/attr.c`); when set on a symbol the evaluator holds all arguments and additionally bypasses Sequence flattening, `Unevaluated` stripping, and upvalue lookup for that head.

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/attr.c`](https://github.com/stblake/mathilda/blob/main/src/attr.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
