# InverseFunctions

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
InverseFunctions is an option for Solve that enables the
    inverse-function specialist for elementary invertible heads
    (Log, Exp, Sin, Cos, Tan, ArcSin, ArcCos, Sinh, ..., and
    integer Power).  Default: Automatic (enabled).  Setting it
    to False disables the specialist; equations that can only
    be solved through inversion then return unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/solve.c`](https://github.com/stblake/mathilda/blob/main/src/solve.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
