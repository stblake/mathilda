# NumericQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NumericQ[expr] gives True if expr is a numeric quantity, and False otherwise.`**

<details>
<summary>Notes</summary>

An expression is considered a numeric quantity if it is either an explicit number or a mathematical constant such as Pi, or is a function that has attribute NumericFunction and all of whose arguments are numeric quantities.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= NumericQ[Pi]
Out[1]= True

In[2]:= NumericQ[x]
Out[2]= False

In[3]:= NumericQ[Gamma[1/2] + Zeta[3]]
Out[3]= True

In[4]:= NumericQ[x + 1]
Out[4]= False

In[5]:= NumericQ[Sin[2] + Log[3]]
Out[5]= True
```

## Implementation notes

`builtin_numericq` (1-arg) returns `True`/`False` by calling the recursive predicate `is_numeric_quantity`. That predicate returns true for `EXPR_INTEGER`/`EXPR_REAL`/`EXPR_BIGINT`/`EXPR_MPFR`; for the named numeric constants `Pi`, `E`, `I`, `Infinity`, `ComplexInfinity`, `EulerGamma`, `GoldenRatio`, `Catalan`, `Degree`; for `Complex[...]` and `Rational[...]` heads; and for any function whose head carries `ATTR_NUMERICFUNCTION` *provided every argument is itself numeric* (recursive check). Everything else — bare symbols, non-numeric heads — yields `False`. Unlike `NumberQ`, this resolves the "would evaluate to a number" question structurally via the attribute system rather than by numericalizing.

**Attributes:** `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric.c)

## Notes & additional examples

### Notes

An expression is numeric if it is an explicit number, a constant such as `Pi`, or
a `NumericFunction` whose arguments are all numeric. The test is structural and
recursive — it never evaluates the expression to a number — so `Gamma[1/2] +
Zeta[3]` is reported numeric while a single symbolic leaf such as `x` makes the
whole expression non-numeric.
