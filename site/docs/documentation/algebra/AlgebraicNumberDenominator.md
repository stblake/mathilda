# AlgebraicNumberDenominator

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlgebraicNumberDenominator[a]`**

gives the smallest positive integer n such that n a is an algebraic integer.  a may be a rational, a radical, a Root object, or an AlgebraicNumber object; for an algebraic integer the denominator is 1. Threads over lists.  Requires FLINT.  Attributes: Listable, Protected.

## Examples (9)

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

### Applications (4)

An algebraic integer: denominator 1

```mathematica
In[6]:= AlgebraicNumberDenominator[Sqrt[2]]
Out[6]= 1
```

The smallest n making n a an algebraic integer

```mathematica
In[7]:= AlgebraicNumberDenominator[1/5 + Sqrt[2]]
Out[7]= 5
```

A plain rational denominator

```mathematica
In[8]:= AlgebraicNumberDenominator[1/2]
Out[8]= 2
```

Threads over a list

```mathematica
In[9]:= AlgebraicNumberDenominator[{1/5 + Sqrt[2], Sqrt[2], 1/2}]
Out[9]= {5, 1, 2}
```

## Algorithm

algebraicnumberdenominator.c — AlgebraicNumberDenominator[a].

See algebraicnumberdenominator.h. The value — the smallest positive integer n with n a an algebraic integer — is computed exactly by flint_qqbar_algebraic_number_denominator (a per-prime valuation over the minimal polynomial). This file only checks the argument shape, maps the tri-state engine result to Integer / message / unevaluated, and (being Listable) lets the evaluator thread over a list of algebraic numbers.

## Implementation notes

**Algorithm.** `builtin_algebraicnumberdenominator` checks arity 1 and delegates
to `flint_qqbar_algebraic_number_denominator`, which computes exactly the
smallest positive integer `d` such that `d·x` is an algebraic integer. It is
**not** `qqbar_denominator` (the leading coefficient `a_n` of the primitive
integer minimal polynomial), which only over-estimates. Writing the monic
minimal polynomial of `d·x`, its `x^i` coefficient is `(a_i/a_n)·d^{n−i}`, so
`d·x` is an algebraic integer iff for every `i < n` the denominator `q_i` of
`a_i/a_n` divides `d^{n−i}`. Since `q_i | a_n`, the minimal `d` divides `a_n`.
The engine factors `a_n` once (`fmpz_factor`) and, per prime `P | a_n`, takes

```
v_P(d) = max_{i<n} ceil( v_P(q_i) / (n−i) ),   v_P(q_i) = max(0, v_P(a_n) − v_P(a_i))
```

(a zero coefficient `a_i` contributes no constraint). Worked example from the
source: `1/5 + Sqrt[2]` has `p = 25x² − 10x − 49`, `a_n = 25`, yet `d = 5`.

A non-constant-algebraic argument gives the engine result `0`, which routes a
`AlgebraicNumberDenominator::nalg` message through `mth_message` (Quiet/Check
funnel) and leaves the expression unevaluated; `−1` (FLINT off) is a silent
decline.

**Data structures.** FLINT `qqbar_t` (minimal polynomial as `fmpz_poly`), an
`fmpz_factor_t` of the leading coefficient, and an `fmpz` valuation scan over the
lower coefficients. The result `d` is returned as an `Integer`/`BigInt` `Expr`.

**Complexity / limits.** One factorisation of the leading coefficient `a_n`
(small in practice) plus an `O(n)` valuation scan. `Listable`, `Protected`;
bounded by the degree cap `QQBAR_DEGREE_CAP = 120`. For any algebraic integer the
answer is `1`.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [Root](../../solutions-of-equations/Root/), [AlgebraicNumber](../../algebra/AlgebraicNumber/)

- H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.1–4.2 (algebraic integers; denominators of algebraic numbers).
- The FLINT library (https://flintlib.org), `qqbar` module — minimal polynomial and `fmpz_factor`.
- Source: [`src/poly/algebraicnumberdenominator.c`](https://github.com/stblake/mathilda/blob/main/src/poly/algebraicnumberdenominator.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumberdenominator.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberdenominator.c)

## Notes & additional examples

### Notes

`AlgebraicNumberDenominator[a]` is the smallest positive integer `n` such that
`n a` is an algebraic integer. It is computed exactly, and is **not** simply the
leading coefficient of the minimal polynomial, which only bounds it from above:
`1/5 + Sqrt[2]` has minimal polynomial `25 x^2 - 10 x - 49` (leading coefficient
`25`), yet its denominator is `5`. The engine factors the leading coefficient
once and takes a per-prime valuation over the lower coefficients.

`Listable` and `Protected`; requires FLINT. For any algebraic integer the answer
is `1`.
