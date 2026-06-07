# Denominator

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Denominator[expr]
    gives the denominator of expr regarded as a rational expression.
    Collects factors of expr that carry a superficially negative
    exponent, inverted; returns 1 when no such factors exist.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Denominator[(x-1)(x-2)/(x-3)^2]
Out[1]= (-3 + x)^2

In[2]:= Denominator[3/7 + I/11]
Out[2]= 77
```

## Implementation notes

- `Protected`, `Listable`.
- Picks out terms which have superficially negative exponents.
- Can be used on rational and complex numbers.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on rational normal forms.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Denominator[6/8]
Out[1]= 4
```

```mathematica
In[1]:= Denominator[(x+1)/(x-1)]
Out[1]= -1 + x
```

```mathematica
In[1]:= Denominator[3]
Out[1]= 1
```

```mathematica
In[1]:= Denominator[a/b + c/d]
Out[1]= 1
```

### Notes

Denominator returns the bottom of the structural rational form. Rational constants
are reduced first, so `Denominator[6/8] = 4`. A non-fractional expression such as
an integer has denominator `1`. For symbolic quotients it returns the literal
denominator, e.g. `-1 + x` for `(x+1)/(x-1)`. Because Mathilda does not implicitly
`Together` a sum, `Denominator[a/b + c/d]` returns `1` (the expression is a `Plus`
with no overall denominator); call `Together` first to obtain `b d`.
