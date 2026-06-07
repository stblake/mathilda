# Accuracy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Accuracy[x]
    Returns the number of digits of accuracy in x — equal to
    Precision[x] − Log10[Abs[x]]. Exact numbers (including exact 0)
    return Infinity. Inexact zeros are finite: machine 0. returns
    ≈ 323.607; MPFR 0 of precision p returns p digits.
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
