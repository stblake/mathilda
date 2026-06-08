# Khinchin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Khinchin
    is Khinchin's constant K (also Khintchine's constant), with numerical
    value ~= 2.68545.
Khinchin's constant is the limiting geometric mean of the partial
quotients in the continued-fraction expansion of almost every real
number, given by the product over s >= 1 of (1 + 1/(s (s + 2)))^Log2[s].
It is a mathematical constant: it has attributes Constant and Protected,
NumericQ[Khinchin] is True, and D[Khinchin, x] is 0. N[Khinchin, prec]
evaluates it to any precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[Khinchin] = {Constant,

**Attributes:** `Constant`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)
