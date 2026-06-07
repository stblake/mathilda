# Divide

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
x / y or Divide[x, y] represents x / y; rewritten by the evaluator to
Times[x, Power[y, -1]] so it inherits Times's flattening and ordering.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_divide` (`src/arithmetic.c`) evaluates `Divide[num, den]`. It special-cases numerics: if either operand is a `Real`, it computes a `double` quotient (coercing integer/bigint operands via `mpz_get_d`), emitting `Power::infy` and returning `ComplexInfinity` on a zero denominator. If both are exact integers/rationals (`is_rational`), it forms the reduced rational `n1*d2 / d1*n2` via `make_rational`, handling `x/0 -> ComplexInfinity` and `0/0 -> Indeterminate` with the matching diagnostics. Otherwise it falls back to the symbolic canonical form `Times[num, Power[den, -1]]`, which the evaluator then simplifies. (Divide carries no Hold attributes, so arguments arrive pre-evaluated.)

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/arithmetic.c`](https://github.com/stblake/mathilda/blob/main/src/arithmetic.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
