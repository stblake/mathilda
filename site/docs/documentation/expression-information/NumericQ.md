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
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= NumericQ[Pi]
Out[1]= True

In[2]:= NumericQ[x]
Out[2]= False
```

A deep tree of constants and transcendental functions is recognised as numeric
without any value being computed:

```mathematica
In[1]:= NumericQ[Gamma[1/2] + Zeta[3]]
Out[1]= True
```

One non-numeric leaf is enough to spoil the whole expression:

```mathematica
In[1]:= NumericQ[x + 1]
Out[1]= False
```

The classification looks through `NumericFunction` heads recursively, so mixed
elementary and special functions of numeric arguments still qualify:

```mathematica
In[1]:= NumericQ[Sin[2] + Log[3]]
Out[1]= True
```

### Notes

An expression is numeric if it is an explicit number, a constant such as `Pi`, or
a `NumericFunction` whose arguments are all numeric. The test is structural and
recursive — it never evaluates the expression to a number — so `Gamma[1/2] +
Zeta[3]` is reported numeric while a single symbolic leaf such as `x` makes the
whole expression non-numeric.
