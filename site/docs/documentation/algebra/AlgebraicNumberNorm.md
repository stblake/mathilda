# AlgebraicNumberNorm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlgebraicNumberNorm[a]`**

gives the field norm of the algebraic number a: the product of the conjugates of a over the rationals, equivalently the product of the roots of a's minimal polynomial.  a may be an integer, a rational, a radical, GoldenRatio, a Root object, or an AlgebraicNumber object. AlgebraicNumberNorm\[a, Extension -\> theta\] gives the norm relative to the field Q(theta), for a an element of Q(theta).  Threads over lists. Requires FLINT.  Attributes: Listable, Protected.

## Examples (11)

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

### Applications (5)

Product of the roots of x^2 - 2

```mathematica
In[7]:= AlgebraicNumberNorm[Sqrt[2]]
Out[7]= -2
```

(1 + Sqrt[2])(1 - Sqrt[2])

```mathematica
In[8]:= AlgebraicNumberNorm[1 + Sqrt[2]]
Out[8]= -1
```

The norm of the golden ratio

```mathematica
In[9]:= AlgebraicNumberNorm[GoldenRatio]
Out[9]= -1
```

A degree-4 field: product of four conjugates

```mathematica
In[10]:= AlgebraicNumberNorm[Sqrt[2] + Sqrt[3]]
Out[10]= 1
```

The relative norm over Q(Sqrt[2])

```mathematica
In[11]:= AlgebraicNumberNorm[1 + Sqrt[2], Extension -> Sqrt[2]]
Out[11]= -1
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

**Algorithm.** `builtin_algebraicnumbernorm` separates a trailing
`Extension -> theta` option from the positional argument
(`extract_extension_option`), checks arity, and delegates to
`flint_qqbar_algebraic_number_norm`. The **absolute norm** `N_{Q(a)/Q}(a)` is
read straight off `a`'s primitive integer minimal polynomial
`P(x) = c_n x^n + … + c_0` (content 1, `c_n > 0`): the product of its `n` roots
is `(−1)^n · c_0/c_n` (`qqbar_abs_norm`). With `Extension -> theta` the
**relative norm** `N_{Q(theta)/Q}(a)` is computed by transitivity of the norm in
the tower `Q ⊆ Q(a) ⊆ K = Q(theta)`:

```
N_{K/Q}(a) = N_{Q(a)/Q}(a)^{[K:Q(a)]} = (absolute norm)^{n/d}
```

with `n = deg minpoly(theta)` (via the algebraic-integer generator `phi`) and
`d = deg minpoly(a)`; membership `a ∈ Q(theta)` is decided by
`qqbar_express_in_field` (escalating working precision). The engine's return
codes map to: success (`Integer`/`Rational`); `AlgebraicNumberNorm::ext` (`a`
not an element of `Q(theta)`); `AlgebraicNumberNorm::nalg` (not a constant
algebraic number); silent decline (FLINT off). Both messages route through
`mth_message`.

**Data structures.** FLINT `qqbar_t` for `a` and `theta`, an `fmpq` accumulator
for the norm, and an `fmpq_poly` for the field-membership expression in the
relative case. No `AlgebraicNumber` object is constructed.

**Complexity / limits.** The absolute norm is `O(1)` once `a` is converted (one
minimal-polynomial read); the relative case adds a primitive-element build and an
express-in-field solve. `Listable` (a trailing `Extension` option is repeated
across the list), `Protected`; bounded by the degree cap
`QQBAR_DEGREE_CAP = 120`.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [GoldenRatio](../../mathematical-constants/GoldenRatio/), [Root](../../solutions-of-equations/Root/), [AlgebraicNumber](../../algebra/AlgebraicNumber/)

- H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.3 (norm and trace of algebraic numbers; transitivity in a tower).
- The FLINT library (https://flintlib.org), `qqbar` module — minimal polynomial and field membership.
- Source: [`src/poly/algebraicnumbernorm.c`](https://github.com/stblake/mathilda/blob/main/src/poly/algebraicnumbernorm.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumbernorm.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumbernorm.c)

## Notes & additional examples

### Notes

The absolute norm `N_{Q(a)/Q}(a)` is the product of the conjugates of `a` —
equivalently the product of the roots of `a`'s minimal polynomial, read off as
`(-1)^deg` times the monic constant term. `AlgebraicNumberNorm[a, Extension ->
theta]` gives the relative norm `N_{Q(theta)/Q}(a)` for `a` an element of
`Q(theta)`, equal to the absolute norm raised to the tower index `[Q(theta):Q(a)]`.

`Listable` (a trailing `Extension` option is repeated across the list) and
`Protected`; requires FLINT. If `a` is not an element of the stated extension,
the call reports and stays unevaluated.
