# MachineIntegerQ

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`MachineIntegerQ[expr]`**

gives True if expr is a machine-word (64-bit) integer, False otherwise.

<details>
<summary>Notes</summary>

Unlike IntegerQ, returns False for a BigInt: MachineIntegerQ\[2^100\] is False because that value has been promoted out of a 64-bit word.

</details>

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
