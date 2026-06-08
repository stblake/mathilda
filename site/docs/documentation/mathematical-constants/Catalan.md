# Catalan

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
Catalan
    is Catalan's constant, with numerical value ~= 0.915966.
Catalan is the sum over k >= 0 of (-1)^k (2 k + 1)^-2. It is a
mathematical constant: it has attributes Constant and Protected,
NumericQ[Catalan] is True, and D[Catalan, x] is 0. N[Catalan, prec]
evaluates it to any precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `Constant`, `Protected`. `Attributes[Catalan] = {Constant,

**Attributes:** `Constant`, `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/mathematical-constants.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/mathematical-constants.md)
