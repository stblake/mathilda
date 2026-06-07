# RootSum

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RootSum[Function[t, p[t]], Function[t, body[t]]]
    The formal sum of body[\[Alpha]] over the roots \[Alpha] of
    p[\[Alpha]] == 0.  Held symbolic form, used by the rational
    integrator's NaiveLogPart fallback when the logarithmic part
    cannot be expressed in closed-form real elementary functions.
    Differentiation threads through the body Function:
      D[RootSum[f1, Function[t, body]], x]
        == RootSum[f1, Function[t, D[body, x]]].
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/root.c`](https://github.com/stblake/mathilda/blob/main/src/root.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
