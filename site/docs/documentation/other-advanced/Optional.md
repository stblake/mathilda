# Optional

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
patt:def or Optional[patt, def]
    is a pattern object that matches patt if it is present; if patt is
    omitted from the argument sequence, def is used in its place.
patt_. (sugar for Optional[patt_, Default[f]]) draws the default value
from Default[f] at the call site.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
