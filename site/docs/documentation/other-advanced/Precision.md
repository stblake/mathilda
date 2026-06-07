# Precision

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Precision[x]
    Returns the number of decimal digits of precision in x.
    Exact numbers return Infinity; machine-precision reals return
    the symbol MachinePrecision; MPFR values return their declared
    precision in decimal digits.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/precision.c`](https://github.com/stblake/mathilda/blob/main/src/precision.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
