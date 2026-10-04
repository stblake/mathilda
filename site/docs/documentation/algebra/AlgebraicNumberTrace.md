# AlgebraicNumberTrace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlgebraicNumberTrace[a]`**

gives the field trace of the algebraic number a: the sum of the conjugates of a over the rationals, equivalently the sum of the roots of a's minimal polynomial.  a may be an integer, a rational, a radical, GoldenRatio, a Root object, or an AlgebraicNumber object. AlgebraicNumberTrace\[a, Extension -\> theta\] gives the trace relative to the field Q(theta), for a an element of Q(theta).  Threads over lists. Requires FLINT.  Attributes: Listable, Protected.

## Examples (10)

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

### Applications (4)

Sqrt[2] + (-Sqrt[2]) = 0

```mathematica
In[7]:= AlgebraicNumberTrace[Sqrt[2]]
Out[7]= 0
```

(1 + Sqrt[2]) + (1 - Sqrt[2])

```mathematica
In[8]:= AlgebraicNumberTrace[1 + Sqrt[2]]
Out[8]= 2
```

The trace of the golden ratio

```mathematica
In[9]:= AlgebraicNumberTrace[GoldenRatio]
Out[9]= 1
```

Trace of 1 + 3 Sqrt[2]

```mathematica
In[10]:= AlgebraicNumberTrace[AlgebraicNumber[Sqrt[2], {1, 3}]]
Out[10]= 2
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

**Algorithm.** `builtin_algebraicnumbertrace` mirrors `AlgebraicNumberNorm`
structurally: it separates a trailing `Extension -> theta` option, checks arity,
and delegates to `flint_qqbar_algebraic_number_trace`. The **absolute trace**
`Tr_{Q(a)/Q}(a)` is the sum of the roots of `a`'s primitive integer minimal
polynomial `P(x) = c_n x^n + … + c_0`, equal to `−c_{n−1}/c_n`
(`qqbar_abs_trace`). With `Extension -> theta` the **relative trace**
`Tr_{Q(theta)/Q}(a)` follows from transitivity of the trace in the tower
`Q ⊆ Q(a) ⊆ K = Q(theta)`:

```
Tr_{K/Q}(a) = [K:Q(a)] · Tr_{Q(a)/Q}(a) = (n/d) · (absolute trace)
```

with `n = deg minpoly(theta)`, `d = deg minpoly(a)`, and `a ∈ Q(theta)` decided
by `qqbar_express_in_field`. The contrast with the norm is exact: the trace is
**additive**, so it *scales* by the tower index `n/d`, whereas the multiplicative
norm is *raised to the power* `n/d`. Return codes map to the
`AlgebraicNumberTrace::ext` (`a` not in `Q(theta)`) and
`AlgebraicNumberTrace::nalg` (not a constant algebraic number) messages through
`mth_message`, a success `Integer`/`Rational`, or a silent decline (FLINT off).

**Data structures.** FLINT `qqbar_t` for `a` and `theta`, an `fmpq` accumulator
for the trace, and an `fmpq_poly` for membership in the relative case.

**Complexity / limits.** `O(1)` minimal-polynomial read for the absolute trace;
the relative case adds a primitive-element build and an express-in-field solve.
`Listable`, `Protected`; bounded by the degree cap `QQBAR_DEGREE_CAP = 120`.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [GoldenRatio](../../mathematical-constants/GoldenRatio/), [Root](../../solutions-of-equations/Root/), [AlgebraicNumber](../../algebra/AlgebraicNumber/), [AlgebraicNumberNorm](../../algebra/AlgebraicNumberNorm/)

- H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.3 (norm and trace of algebraic numbers; transitivity in a tower).
- The FLINT library (https://flintlib.org), `qqbar` module — minimal polynomial and field membership.
- Source: [`src/poly/algebraicnumbertrace.c`](https://github.com/stblake/mathilda/blob/main/src/poly/algebraicnumbertrace.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumbertrace.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumbertrace.c)

## Notes & additional examples

### Notes

The absolute trace `Tr_{Q(a)/Q}(a)` is the sum of the conjugates of `a` —
equivalently the sum of the roots of `a`'s minimal polynomial, read off as
`-(coeff of x^{deg-1}) / (leading coeff)`. `AlgebraicNumberTrace[a, Extension ->
theta]` gives the relative trace `Tr_{Q(theta)/Q}(a)` for `a` in `Q(theta)`.

Where the norm is multiplicative and *raises* the absolute value to the tower
index, the trace is additive and *scales* by it. `Listable` and `Protected`;
requires FLINT.
