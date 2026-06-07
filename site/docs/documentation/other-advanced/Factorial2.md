# Factorial2

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Factorial2[n] (also typeset n!!) gives the double factorial of n.
For non-negative integer n: n!! = n * (n-2) * (n-4) * ... down to 2 (n even) or 1 (n odd).
Special values: 0!! = 1, (-1)!! = 1.
Negative even integers and negative odd integers below -1 give ComplexInfinity.
Factorial2 stays unevaluated on symbolic arguments.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
