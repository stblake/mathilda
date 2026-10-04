# AlgebraicIntegerQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlgebraicIntegerQ[x]`**

gives True if x is an algebraic integer (a root of a monic polynomial with integer coefficients) and False otherwise.  Rational integers are algebraic integers; non-integer rationals and non-algebraic quantities are not.  Requires FLINT.  Attribute: Protected.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= AlgebraicIntegerQ[Sqrt[2]]
Out[1]= True

In[2]:= AlgebraicIntegerQ[1/2]
Out[2]= False

In[3]:= AlgebraicIntegerQ[(1 + Sqrt[5])/2]
Out[3]= True

In[4]:= AlgebraicIntegerQ[2^(1/3)/2]
Out[4]= False
```

### Applications (6)

A root of the monic x^2 - 2

```mathematica
In[5]:= AlgebraicIntegerQ[Sqrt[2]]
Out[5]= True
```

A non-integer rational is not an algebraic integer

```mathematica
In[6]:= AlgebraicIntegerQ[1/2]
Out[6]= False
```

A root of the monic x^2 - x - 1

```mathematica
In[7]:= AlgebraicIntegerQ[GoldenRatio]
Out[7]= True
```

The same number, spelled as a radical

```mathematica
In[8]:= AlgebraicIntegerQ[(1 + Sqrt[5])/2]
Out[8]= True
```

Gaussian integers are algebraic integers

```mathematica
In[9]:= AlgebraicIntegerQ[2 + 3 I]
Out[9]= True
```

Scaling by 1/3 leaves the monogenic ring

```mathematica
In[10]:= AlgebraicIntegerQ[Sqrt[2]/3]
Out[10]= False
```

## Implementation notes

**Algorithm.** `builtin_algebraicintegerq` is a thin wrapper: it checks arity 1
and hands the argument to `flint_qqbar_algebraic_integer_q`, mapping the
tri-state result to `True` (1), `False` (0), or unevaluated (−1, FLINT compiled
out). The decision is exact. The engine converts `x` to a FLINT `qqbar_t` via
`to_qqbar` — the shared converter that accepts integers, rationals, radicals
`Power[base, p/q]`, roots of unity, the imaginary unit, `Root[]` objects and
`GoldenRatio`, combined by `+ - * / ^`. `x` is an algebraic integer **iff the
leading coefficient of its primitive integer minimal polynomial is 1** (monic):
`fmpz_poly_get_coeff_fmpz(lead, QQBAR_POLY(v), deg)` and test `fmpz_is_one`. A
rational `p/q` has minimal polynomial `q x − p`, so only true integers (`q = 1`)
qualify. Anything that is not a constant algebraic number — a free symbol, `Pi`,
`Log[2]` — fails the `to_qqbar` conversion and returns `0` (`False`), matching
WL.

**Data structures.** The FLINT `qqbar_t` (an exact algebraic number: its
primitive integer minimal polynomial as an `fmpz_poly` plus an isolating complex
enclosure). No `AlgebraicNumber[...]` object is built — this is a predicate that
only reads the minimal polynomial's leading coefficient.

**Complexity / limits.** Cost is dominated by the `to_qqbar` conversion of the
argument (field arithmetic over the compositum, bounded by the degree cap
`QQBAR_DEGREE_CAP = 120`). `Protected`, and deliberately **not** `Listable` —
`AlgebraicIntegerQ[list]` asks whether the list itself is an algebraic integer
(`False`), not whether its elements are. Declines (stays unevaluated) only when
FLINT is unavailable.

**Attributes:** `Protected`.

## References

**See also:** [Pi](../../mathematical-constants/Pi/)

- H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.1 (algebraic integers and minimal polynomials).
- The FLINT library (https://flintlib.org), `qqbar` module — exact real and complex algebraic numbers via minimal polynomial plus isolating enclosure.
- Source: [`src/poly/algebraicintegerq.c`](https://github.com/stblake/mathilda/blob/main/src/poly/algebraicintegerq.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumberdenominator.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberdenominator.c)
- Tests: [`tests/test_numberfieldintegralbasis.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numberfieldintegralbasis.c)

## Notes & additional examples

### Notes

The test is exact: `x` is an algebraic integer iff the leading coefficient of its
primitive integer minimal polynomial is `1`. A rational `p/q` has minimal
polynomial `q x - p`, so only ordinary integers qualify among the rationals.

Anything that is not a constant algebraic number — a free symbol, `Pi`,
`Log[2]` — is simply not an algebraic integer, and returns `False`. The predicate
is **not** `Listable`: `AlgebraicIntegerQ[list]` asks about the list itself.
Requires FLINT; with FLINT compiled out the call stays unevaluated.
