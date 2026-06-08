# GoldenRatio

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
GoldenRatio
    is the golden ratio phi = (1 + Sqrt[5])/2, with numerical value
    ~= 1.61803.
GoldenRatio is the positive root of x^2 == x + 1. It is a mathematical
constant: it has attributes Constant and Protected, NumericQ[GoldenRatio]
is True, and D[GoldenRatio, x] is 0. N[GoldenRatio, prec] evaluates it to
any precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[GoldenRatio] = {Constant,

**Attributes:** `Constant`, `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)
