# Limit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Limit[f, x -> a]
    finds the limit of f as x approaches a.
Limit[f, {x1 -> a1, ..., xn -> an}]
    iterated limit, applied rightmost-first.
Limit[f, {x1, ..., xn} -> {a1, ..., an}]
    multivariate (joint) limit.
Limit[f, x -> a, Direction -> d]
    specifies the direction of approach:
      Reals or "TwoSided" -- default two-sided limit
      "FromAbove" or -1   -- approach from above (x -> a^+)
      "FromBelow" or +1   -- approach from below (x -> a^-)
      Complexes           -- limit over all complex directions

May return a finite value, Infinity, -Infinity, ComplexInfinity,
Indeterminate, Interval[{lo, hi}], or the original unevaluated
expression when the limit cannot be determined.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldAll`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 3.
- Source: [`src/calculus/limit.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/limit.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Limit[(Cos[x] - 1)/x^2, x -> 0]
Out[1]= -1/2
```

```mathematica
In[1]:= Limit[(1 + 1/x)^x, x -> Infinity]
Out[1]= E
```

```mathematica
In[1]:= Limit[(x^2 - 1)/(x - 1), x -> 1]
Out[1]= 2
```

```mathematica
In[1]:= Limit[Tan[x]/x, x -> 0]
Out[1]= 1
```

### Notes

`Limit[f, x -> a]` resolves the standard removable-singularity and indeterminate forms, including the classic `(1 + 1/x)^x -> E` and `0/0` cancellations such as `(x^2 - 1)/(x - 1)`. The `Direction` option selects one-sided (`"FromAbove"`/`"FromBelow"`) or complex approaches; the default is two-sided. Results may be a finite value, `Infinity`, `ComplexInfinity`, `Indeterminate`, an `Interval`, or the original expression unevaluated when the limit cannot be determined. Iterated and joint multivariate limits are supported through the list forms of the second argument.
