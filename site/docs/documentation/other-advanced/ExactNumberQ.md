# ExactNumberQ

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`ExactNumberQ[expr]`**

gives True if expr is an exact number, False otherwise.

<details>
<summary>Notes</summary>

Exact numbers are integers, rationals, and Complex numbers whose parts are exact. Reals and MPFR numbers are inexact, so ExactNumberQ is False.

</details>

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
