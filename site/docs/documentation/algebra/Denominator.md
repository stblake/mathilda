# Denominator

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Denominator[expr]`**

gives the denominator of expr regarded as a rational expression. Collects factors of expr that carry a superficially negative exponent, inverted; returns 1 when no such factors exist.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Denominator[(x-1)(x-2)/(x-3)^2]
Out[1]= (-3 + x)^2

In[2]:= Denominator[3/7 + I/11]
Out[2]= 77
```

### Applications (6)

```mathematica
In[3]:= Denominator[6/8]
Out[3]= 4

In[4]:= Denominator[(x+1)/(x-1)]
Out[4]= -1 + x

In[5]:= Denominator[a/b + c/d]
Out[5]= 1

In[6]:= Denominator[Together[a/b + c/d]]
Out[6]= b d

In[7]:= Denominator[(x^2-1)/((x-2)^3 (x+5))]
Out[7]= (5 + x) (-2 + x)^3

In[8]:= Denominator[x^(-2) y^3 z^(-1)]
Out[8]= x^2 z
```

## Implementation notes

`builtin_denominator` calls `extract_num_den` and returns the denominator part (freeing the numerator). `extract_num_den` recognises `Rational[n, d]` (returns `d`); `Complex` (clears to a common integer denominator); `Power[b, e]`/`Exp[e]` with a superficially-negative exponent or a `Plus` exponent split into positive/negative pieces (the negative-exponent base becomes the denominator); and `Times`, which recurses into each factor and multiplies the collected denominators. A factor with no denominator contributes `1`. `Numerator` in the same file is the symmetric accessor.

- `Protected`, `Listable`.
- Picks out terms which have superficially negative exponents.
- Can be used on rational and complex numbers.

**Attributes:** `Listable`, `Protected`.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on rational normal forms.
- Source: [`src/rat.c`](https://github.com/stblake/mathilda/blob/main/src/rat.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_bignum_rational_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bignum_rational_numeric.c)
- Tests: [`tests/test_integrate_newton_leibniz.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_newton_leibniz.c)
- Tests: [`tests/test_integrate_risch_transcendental.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integrate_risch_transcendental.c)
- Tests: [`tests/test_rat.c`](https://github.com/stblake/mathilda/blob/main/tests/test_rat.c)

## Notes & additional examples

### Notes

Denominator returns the bottom of the structural rational form. Rational constants
are reduced first, so `Denominator[6/8] = 4`. A non-fractional expression such as
an integer has denominator `1`. For symbolic quotients it returns the literal
denominator, e.g. `-1 + x` for `(x+1)/(x-1)`. Because Mathilda does not implicitly
`Together` a sum, `Denominator[a/b + c/d]` returns `1` (the expression is a `Plus`
with no overall denominator); call `Together` first to obtain `b d`.
