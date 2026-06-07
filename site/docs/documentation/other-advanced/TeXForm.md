# TeXForm

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
TeXForm[expr]
    prints as a TeX version of expr.
TeXForm produces AMS-LaTeX-compatible TeX output.
When an input evaluates to TeXForm[expr], TeXForm does not appear in the output.
TeXForm translates standard mathematical functions and operations.
Following standard mathematical conventions, single-character symbol names are given in italic font, while multiple character names are given in roman font.
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
