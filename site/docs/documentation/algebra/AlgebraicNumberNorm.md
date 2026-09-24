# AlgebraicNumberNorm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlgebraicNumberNorm[a]`**

gives the field norm of the algebraic number a: the product of the conjugates of a over the rationals, equivalently the product of the roots of a's minimal polynomial.  a may be an integer, a rational, a radical, GoldenRatio, a Root object, or an AlgebraicNumber object. AlgebraicNumberNorm\[a, Extension -\> theta\] gives the norm relative to the field Q(theta), for a an element of Q(theta).  Threads over lists. Requires FLINT.  Attributes: Listable, Protected.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= AlgebraicNumberNorm[Sqrt[2]]
Out[1]= -2

In[2]:= AlgebraicNumberNorm[GoldenRatio]
Out[2]= -1

In[3]:= AlgebraicNumberNorm[1/Sqrt[Sqrt[2] + 3]]
Out[3]= 1/7

In[4]:= AlgebraicNumberNorm[{2 Sqrt[2], E^(Pi I/8), 1 + I}]
Out[4]= {-8, 1, 2}
```

### Options (2)

```mathematica
In[5]:= AlgebraicNumberNorm[Sqrt[2], Extension -> E^(Pi I/4)]
Out[5]= 4

In[6]:= AlgebraicNumberNorm[{2, Sqrt[5]}, Extension -> Sqrt[5]]
Out[6]= {4, -5}
```

## Algorithm

algebraicnumbernorm.c — AlgebraicNumberNorm[a].

```text
See algebraicnumbernorm.h.  The value — the field norm of an algebraic number,
```

optionally relative to a field Q(theta) given by Extension -> theta — is computed exactly by flint_qqbar_algebraic_number_norm (both cases reduce to a

```text
minimal-polynomial read; see flint_qqbar.c).  This file only separates the
```

positional argument from the Extension option, checks arity, maps the tri/quad-state engine result to Integer/Rational / message / unevaluated, and (being Listable) lets the evaluator thread over a list of algebraic numbers.

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [GoldenRatio](../../mathematical-constants/GoldenRatio/), [Root](../../solutions-of-equations/Root/), [AlgebraicNumber](../../algebra/AlgebraicNumber/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumbernorm.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumbernorm.c)
