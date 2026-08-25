# InexactNumberQ

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`InexactNumberQ[expr]`**

gives True if expr is an inexact number, False otherwise.

<details>
<summary>Notes</summary>

Inexact numbers are machine reals, arbitrary-precision (MPFR) reals, and Complex numbers with an inexact part. The complement of ExactNumberQ among numbers.

</details>

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
