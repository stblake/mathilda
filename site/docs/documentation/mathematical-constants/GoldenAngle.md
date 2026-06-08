# GoldenAngle

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
GoldenAngle
    is the golden angle (3 - Sqrt[5]) Pi = 2 Pi / GoldenRatio^2, with
    numerical value ~= 2.39996 radians (~= 137.5 degrees).
GoldenAngle is a mathematical constant: it has attributes Constant and
Protected, NumericQ[GoldenAngle] is True, and D[GoldenAngle, x] is 0.
N[GoldenAngle, prec] evaluates it to any precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[GoldenAngle] = {Constant,

**Attributes:** `Constant`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)
