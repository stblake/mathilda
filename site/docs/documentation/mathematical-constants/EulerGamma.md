# EulerGamma

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
EulerGamma
    is Euler's constant gamma, with numerical value ~= 0.5772156649.
EulerGamma is the Euler-Mascheroni constant, the limit of
HarmonicNumber[n] - Log[n] as n -> Infinity. It is a mathematical
constant: it has attributes Constant and Protected, NumericQ[EulerGamma]
is True, and D[EulerGamma, x] is 0. N[EulerGamma, prec] evaluates it to
any precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[EulerGamma] = {Constant,

**Attributes:** `Constant`, `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)
