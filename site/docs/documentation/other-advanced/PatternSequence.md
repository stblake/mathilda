# PatternSequence

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`PatternSequence[p1, p2, ...] is a pattern object that matches a sequence of arguments, each in turn matching p1, p2, ....`**

**`PatternSequence[] matches an empty (zero-length) sequence of arguments.`**

<details>
<summary>Notes</summary>

x:PatternSequence\[...\] binds x to the matched sequence.

</details>

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
