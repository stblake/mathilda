# AlgebraicNumberTrace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlgebraicNumberTrace[a]`**

gives the field trace of the algebraic number a: the sum of the conjugates of a over the rationals, equivalently the sum of the roots of a's minimal polynomial.  a may be an integer, a rational, a radical, GoldenRatio, a Root object, or an AlgebraicNumber object. AlgebraicNumberTrace\[a, Extension -\> theta\] gives the trace relative to the field Q(theta), for a an element of Q(theta).  Threads over lists. Requires FLINT.  Attributes: Listable, Protected.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= AlgebraicNumberTrace[5 + Sqrt[2]]
Out[1]= 10

In[2]:= AlgebraicNumberTrace[GoldenRatio]
Out[2]= 1

In[3]:= AlgebraicNumberTrace[Root[#1^4 + 11 #1^3 + #1^2 + #1 + 1 &, 1]]
Out[3]= -11

In[4]:= AlgebraicNumberTrace[{5 + Sqrt[2], E^(Pi I/8)}]
Out[4]= {10, 0}
```

### Options (2)

```mathematica
In[5]:= AlgebraicNumberTrace[5, Extension -> Sqrt[2]]
Out[5]= 10

In[6]:= AlgebraicNumberTrace[E^(2 Pi I/3), Extension -> ToNumberField[{2^(1/3), E^(2 Pi I/3)}, All][[1, 1]]]
Out[6]= -3
```

## Algorithm

algebraicnumbertrace.c — AlgebraicNumberTrace[a].

```text
See algebraicnumbertrace.h.  The value — the field trace of an algebraic number,
```

optionally relative to a field Q(theta) given by Extension -> theta — is computed exactly by flint_qqbar_algebraic_number_trace (both cases reduce to a

```text
minimal-polynomial read; see flint_qqbar.c).  This file only separates the
```

positional argument from the Extension option, checks arity, maps the tri/quad-state engine result to Integer/Rational / message / unevaluated, and (being Listable) lets the evaluator thread over a list of algebraic numbers.

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [GoldenRatio](../../mathematical-constants/GoldenRatio/), [Root](../../solutions-of-equations/Root/), [AlgebraicNumber](../../algebra/AlgebraicNumber/), [AlgebraicNumberNorm](../../algebra/AlgebraicNumberNorm/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumbertrace.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumbertrace.c)
