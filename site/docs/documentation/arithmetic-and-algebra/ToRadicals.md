# ToRadicals

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ToRadicals[expr]
    attempts to express all Root objects in expr in terms of radicals.

ToRadicals can always give expressions in terms of radicals when the
    highest degree of the polynomial that appears in any Root object is
    four.  Binomial Root objects of the form Root[Function[a #^n + b], k]
    are also reduced to radicals for any degree n.  Other Root objects
    of degree five or higher are returned unchanged.
If Root objects in expr contain parameters, ToRadicals[expr] may yield
    a result that is not equal to expr for all values of the parameters.
ToRadicals automatically threads over lists, equations, inequalities,
    and logic functions.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`.
- Closed-form radicals are always returned when the polynomial has degree

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
