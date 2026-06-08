# Degree

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
Degree
    gives the number of radians in one degree, with numerical value
    Pi/180 (~= 0.0174533).
Multiply by Degree to convert degrees to radians, so 30 Degree is 30
degrees. It is a mathematical constant: it has attributes Constant and
Protected, NumericQ[Degree] is True, and D[Degree, x] is 0. N[Degree,
prec] evaluates it to any precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[Degree] = {Constant,

**Attributes:** `Constant`, `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)
