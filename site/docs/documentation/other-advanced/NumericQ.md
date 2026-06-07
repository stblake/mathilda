# NumericQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NumericQ[expr] gives True if expr is a numeric quantity, and False otherwise.
An expression is considered a numeric quantity if it is either an explicit number or a mathematical constant such as Pi, or is a function that has attribute NumericFunction and all of whose arguments are numeric quantities.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_numericq` (1-arg) returns `True`/`False` by calling the recursive predicate `is_numeric_quantity`. That predicate returns true for `EXPR_INTEGER`/`EXPR_REAL`/`EXPR_BIGINT`/`EXPR_MPFR`; for the named numeric constants `Pi`, `E`, `I`, `Infinity`, `ComplexInfinity`, `EulerGamma`, `GoldenRatio`, `Catalan`, `Degree`; for `Complex[...]` and `Rational[...]` heads; and for any function whose head carries `ATTR_NUMERICFUNCTION` *provided every argument is itself numeric* (recursive check). Everything else — bare symbols, non-numeric heads — yields `False`. Unlike `NumberQ`, this resolves the "would evaluate to a number" question structurally via the attribute system rather than by numericalizing.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
