# Divide

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Times[x, Power[y, -1]] so it inherits Times's flattening and ordering.`**

<details>
<summary>Notes</summary>

x / y or Divide\[x, y\] represents x / y; rewritten by the evaluator to

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= 10/4
Out[1]= 5/2

In[2]:= Divide[10, 4]
Out[2]= 5/2

In[3]:= a/b
Out[3]= a/b
```

## Implementation notes

`builtin_divide` (`src/arithmetic.c`) evaluates `Divide[num, den]`. It special-cases numerics: if either operand is a `Real`, it computes a `double` quotient (coercing integer/bigint operands via `mpz_get_d`), emitting `Power::infy` and returning `ComplexInfinity` on a zero denominator. If both are exact integers/rationals (`is_rational`), it forms the reduced rational `n1*d2 / d1*n2` via `make_rational`, handling `x/0 -> ComplexInfinity` and `0/0 -> Indeterminate` with the matching diagnostics. Otherwise it falls back to the symbolic canonical form `Times[num, Power[den, -1]]`, which the evaluator then simplifies. (Divide carries no Hold attributes, so arguments arrive pre-evaluated.)

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/arithmetic.c`](https://github.com/stblake/mathilda/blob/main/src/arithmetic.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_eval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval.c)
- Tests: [`tests/test_numeric_domain.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_domain.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`x / y` is rewritten by the evaluator to `Times[x, Power[y, -1]]`, so it inherits Times's flattening and ordering; integer quotients auto-reduce to an exact `Rational`.
