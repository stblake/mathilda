# Glaisher

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Glaisher
    is the Glaisher-Kinkelin constant A, with numerical value ~= 1.28243.
Glaisher's constant satisfies Log[A] == 1/12 - Zeta'[-1], where Zeta is
the Riemann zeta function. It is a mathematical constant: it has
attributes Constant and Protected, NumericQ[Glaisher] is True, and
D[Glaisher, x] is 0. N[Glaisher, prec] evaluates it to any precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[Glaisher] = {Constant,

**Attributes:** `Constant`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)
