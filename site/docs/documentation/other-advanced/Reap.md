# Reap

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Reap[expr]
    Evaluates expr and returns {value, {sown...}}, collecting
    every expression sown by Sow during the evaluation. Reap[expr, patt]
    reaps only tags matching patt; Reap[expr, {p1, ...}] makes one sublist
    per pattern; Reap[expr, patt, f] returns f[tag, {e...}] per tag.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
