# FourierParameters

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
FourierParameters is an option for Fourier and InverseFourier that
specifies the {a, b} convention for the transform. The default {0, 1}
uses the symmetric 1/Sqrt[n] normalisation; {-1, 1} and {1, -1} give the
data-analysis and signal-processing conventions.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
