# Pi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Pi
    is pi, with numerical value ~= 3.14159.
Pi is a mathematical constant: it has attributes Constant and Protected,
NumericQ[Pi] is True, and D[Pi, x] is 0. N[Pi, prec] evaluates it to any
precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[Pi] = {Constant, Protected}`;

**Attributes:** `Constant`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)
