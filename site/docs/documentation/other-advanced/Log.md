# Log

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Log[z]
    gives the principal natural logarithm of z, with branch cut along
    the negative real axis.
Log[b, z]
    gives the logarithm to base b, i.e. Log[z] / Log[b].
Log is Listable. Log[1] = 0, Log[E] = 1, Log[E^n] = n for symbolic n.
Numeric inputs route to libm / MPFR; negative reals yield I Pi + Log[|z|].
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
