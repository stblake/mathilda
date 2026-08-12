# Verbatim

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Verbatim[expr] is a pattern object that matches expr taken literally: the pattern constructs inside expr (Blank, Pattern, ...) are not interpreted, so Verbatim[x_] matches only the literal expression x_.`**

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
