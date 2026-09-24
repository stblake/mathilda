# AlgebraicNumberDenominator

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlgebraicNumberDenominator[a]`**

gives the smallest positive integer n such that n a is an algebraic integer.  a may be a rational, a radical, a Root object, or an AlgebraicNumber object; for an algebraic integer the denominator is 1. Threads over lists.  Requires FLINT.  Attributes: Listable, Protected.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= AlgebraicNumberDenominator[1/Sqrt[3]]
Out[1]= 3

In[2]:= AlgebraicNumberDenominator[(1 + Sqrt[5])/2]
Out[2]= 1

In[3]:= AlgebraicNumberDenominator[AlgebraicNumber[Sqrt[2], {1/5, 1}]]
Out[3]= 5

In[4]:= AlgebraicNumberDenominator[{Sqrt[2], 1/Sqrt[2], 1/3}]
Out[4]= {1, 2, 3}

In[5]:= AlgebraicIntegerQ[AlgebraicNumberDenominator[1/Sqrt[1 + I]] / Sqrt[1 + I]]
Out[5]= True
```

## Algorithm

algebraicnumberdenominator.c — AlgebraicNumberDenominator[a].

See algebraicnumberdenominator.h. The value — the smallest positive integer n with n a an algebraic integer — is computed exactly by flint_qqbar_algebraic_number_denominator (a per-prime valuation over the minimal polynomial). This file only checks the argument shape, maps the tri-state engine result to Integer / message / unevaluated, and (being Listable) lets the evaluator thread over a list of algebraic numbers.

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [Root](../../solutions-of-equations/Root/), [AlgebraicNumber](../../algebra/AlgebraicNumber/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumberdenominator.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberdenominator.c)
